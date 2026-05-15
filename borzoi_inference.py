import numpy as np
import pandas as pd
import pysam
from typing import Dict, Tuple

# baskerville imports
from baskerville import seqnn #seqnn initialises model
from baskerville import dataset # contains functions to deal with targets df and untransformation of predicted values
import json #required for reading the paramters file

import argparse

class BorzoiInputError(Exception):
    """Base exception for Borzoi input creation errors."""


class GeneNotFoundError(BorzoiInputError):
    """Raised when the requested gene cannot be found in the GTF."""


class ReferenceAlleleMismatchError(BorzoiInputError):
    """Raised when the provided reference allele does not match the FASTA."""


class InvalidIntervalError(BorzoiInputError):
    """Raised when the requested sequence window is invalid."""


_DNA_TO_INDEX: Dict[str, int] = {"A": 0, "C": 1, "G": 2, "T": 3}


def load_borzoi(model_file: str, 
                model_parameters_file: str, 
                model_targets_file: str,
                rc=True,
                shifts = [0,1,2]):
    
    # open the parameters file read the model, train paramters and get the sequence length
    with open(model_parameters_file) as model_parameters:
        parameters = json.load(model_parameters)
        parameters_model = parameters["model"]
        # parameters_train = parameters["train"]
        parameter_seqlen = parameters_model["seq_length"]

    # read the model targets file
    borzoi_model_targets = pd.read_csv(model_targets_file, sep="\t", index_col=0)
    # can be put into a separate file
    borzoi_model_targets_stranded_collapsed = dataset.targets_prep_strand(borzoi_model_targets) 
    targets_new_index = dict(zip(borzoi_model_targets.index, np.arange(borzoi_model_targets.shape[0])))
    targets_strand_pair = np.array([targets_new_index[ti] for ti in borzoi_model_targets.strand_pair])
    parameters_model["strand_pair"] = [targets_strand_pair]

    borzoi_seqnn_model = seqnn.SeqNN(parameters_model)
    borzoi_seqnn_model.restore(model_file)
    borzoi_seqnn_model.build_slice(borzoi_model_targets.index)
    borzoi_seqnn_model.build_ensemble(rc, shifts)

    return borzoi_seqnn_model, parameter_seqlen

        
# function to take the one hot encoded sequence, model parameters as input and return the ref and alt predictions (untransformed)
def inference(sequence1hot_ref: np.array,
              sequence1hot_alt: np.array,
              borzoi_seqnn_model):
    
    ref_preds = borzoi_seqnn_model(sequence1hot_ref)
    alt_preds = borzoi_seqnn_model(sequence1hot_alt)
    
    return ref_preds, alt_preds


def untransform_predictions(ref_preds, alt_preds, model_targets_file):
    
    borzoi_model_targets = pd.read_csv(model_targets_file, sep="\t", index_col=0)
    ref_preds = dataset.untransform_preds(ref_preds, borzoi_model_targets)
    alt_preds = dataset.untransform_preds(alt_preds, borzoi_model_targets)
    return ref_preds, alt_preds


def l2Norm(ref_array, alt_array):
    return np.round(np.sqrt(np.sum(np.power(alt_array - ref_array, 2), 0)), 4)

def calculate_l2_norm_for_gene(ref_values, alt_values, gene_start, gene_end):
    pass

def normalise_chromosome(chromosome: str, chr: bool= True) -> str:
    
    chromosome = chromosome.strip()
    if chr:
        chromosome = chromosome if chromosome.startswith("chr") else f"chr{chromosome}"
    else:
        chromosome = chromosome.replace("chr", "")
    return chromosome


def compute_centered_window(position: int, sequence_length: int) -> Tuple[int, int, int]:
    variant_index = sequence_length // 2
    window_start_1based = position - variant_index
    window_end_1based = window_start_1based + sequence_length - 1

    if window_start_1based < 1:
        raise InvalidIntervalError(
            f"Computed window start is < 1: {window_start_1based}. "
            "The requested variant is too close to the start of the chromosome."
        )

    return window_start_1based, window_end_1based, variant_index
  
def fetch_sequence_from_fasta(
    fasta_path: str,
    chromosome: str,
    start_1based: int,
    end_1based: int,
    ) -> str:
    chromosome = normalise_chromosome(chromosome, False)
    with pysam.FastaFile(fasta_path) as fasta:
        references = set(fasta.references)
        if chromosome not in references:
            raise BorzoiInputError(
                f"Chromosome '{chromosome}' not found in FASTA references."
            )

        chrom_length = fasta.get_reference_length(chromosome)
        if end_1based > chrom_length:
            raise InvalidIntervalError(
                f"Computed window end {end_1based} exceeds chromosome length "
                f"{chrom_length} for {chromosome}."
            )

        sequence = fasta.fetch(
            reference=chromosome,
            start=start_1based - 1,
            end=end_1based,
        ).upper()

    expected_length = end_1based - start_1based + 1
    if len(sequence) != expected_length:
        raise BorzoiInputError(
            f"Fetched sequence length {len(sequence)} does not match expected length "
            f"{expected_length}."
        )

    return sequence


def one_hot_encode_sequence(sequence: str) -> np.ndarray:
    sequence = sequence.upper()
    encoded = np.zeros((len(sequence), 4), dtype=np.float32)

    for i, base in enumerate(sequence):
        base_index = _DNA_TO_INDEX.get(base)
        if base_index is not None:
            encoded[i, base_index] = 1.0

    return encoded


def create_borzoi_input(
    chromosome: str,
    position: int,
    ref_allele: str,
    alt_allele: str,
    hg38_fasta: str,
    model_sequence_length: int = 524288,
) -> Tuple[np.ndarray, np.ndarray]:
    
    # handle chromosome input
    chromosome = normalise_chromosome(chromosome)

    # make ref and alt_alleles upper case
    ref_allele, alt_allele = ref_allele.upper(), alt_allele.upper()

    # get the coordinate of the window of sequence length, centered around the SNP
    window_start_1based, window_end_1based, variant_index = compute_centered_window(
        position=position,
        sequence_length=model_sequence_length,
    )

    # get the reference sequence for the window
    ref_sequence = fetch_sequence_from_fasta(
        fasta_path=hg38_fasta,
        chromosome=chromosome,
        start_1based=window_start_1based,
        end_1based=window_end_1based,
    )

    # check that we got the correct sequence
    observed_ref = ref_sequence[variant_index].upper()
    if observed_ref != ref_allele:
        raise ReferenceAlleleMismatchError(
            f"Reference allele mismatch at {chromosome}:{position}. "
            f"Provided ref_allele='{ref_allele}', FASTA has '{observed_ref}'."
        )

    # create alt sequence
    alt_sequence = (
        ref_sequence[:variant_index] + alt_allele + ref_sequence[variant_index + 1 :]
    )

    # create one_hot_encoded_arrays
    ref_array = one_hot_encode_sequence(ref_sequence)[None, :, :]
    alt_array = one_hot_encode_sequence(alt_sequence)[None, :, :]


    return ref_array, alt_array


def main(snp_vcf_file, output_file, model_file:str, model_parameters:str, targets_file:str, hg38_fasta:str, column_suffix:str = ""):

    # open the snp vcf file
    snps = pd.read_csv(snp_vcf_file)
    snps["CHROM"] = snps["CHROM"].astype(str)
    Target_identifiers = pd.read_csv(targets_file, sep="\t", index_col=0)["identifier"]
    #load the borzoi model
    model, model_seq_len = load_borzoi(model_file=model_file,
                                   model_parameters_file=model_parameters,
                                   model_targets_file=targets_file)    
    # for each snp, run the inference function, calculate two different L2 norm values for both strands
    Outputs = []
    for snp_row in snps.itertuples():
        print(f"> Procesing {snp_row.CHROM}:{snp_row.POS} {snp_row.REF}>{snp_row.ALT}")
        ref_array, alt_array = create_borzoi_input(chromosome=snp_row.CHROM,
                                                    position=snp_row.POS,
                                                    ref_allele=snp_row.REF,
                                                    alt_allele=snp_row.ALT,
                                                    hg38_fasta=hg38_fasta,
                                                    model_sequence_length=model_seq_len)
        # what is the shape of the ref and alt arrays? (1, sequence_length, num_targets) (i think)
        ref_preds, alt_preds = inference(sequence1hot_ref=ref_array,
                                            sequence1hot_alt=alt_array,
                                            borzoi_seqnn_model=model)
    
        ref_preds_untransformed, alt_preds_untransformed = untransform_predictions(ref_preds=ref_preds[0],
                                                                                    alt_preds=alt_preds[0],
                                                                                    model_targets_file=targets_file)
        # for each SNP, there will be two vectors of size (1, sequence_length, num_targets) for the ref and alt predictions
        Outputs.append([l2Norm(ref_preds_untransformed[:, i], alt_preds_untransformed[:, i]) for i in range(len(Target_identifiers))])
        print(f"Finished processing {snp_row.CHROM}:{snp_row.POS} {snp_row.REF}>{snp_row.ALT}")
    # outputs will be a list opf lists, of shape (num_snps, num_targets), convert this to dataframe
    Outputs_df = pd.DataFrame(Outputs, columns=[f"L2NORM_{ti}_{column_suffix}" for ti in Target_identifiers])
    # add the outputs to the snps dataframe
    snps = pd.concat([snps, Outputs_df], axis=1)
    
    snps.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--snp-vcf-file", help="Path to the SNP VCF file")
    parser.add_argument("--output-file", help="Path to the output file")
    parser.add_argument("--model-file", help="Path to the Borzoi model file")
    parser.add_argument("--model-parameters", help="Path to the model parameters file")
    parser.add_argument("--targets-file", help="Path to the targets file")
    parser.add_argument("--hg38-fasta", help="Path to the HG38 fasta file")
    parser.add_argument("--column-suffix", help="Suffix for the column names")
    args = parser.parse_args()

    main(snp_vcf_file=args.snp_vcf_file,
         output_file=args.output_file,
         model_file=args.model_file,
         model_parameters=args.model_parameters,
         targets_file=args.targets_file,
         hg38_fasta=args.hg38_fasta,
         column_suffix=args.column_suffix)