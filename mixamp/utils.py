from regex import regex
import math
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess
import argparse
import pandas as pd
import os
import itertools
import numpy as np
from tqdm import tqdm


def generate_random_values(N):
    # Generate N random values
    random_values = np.random.rand(N)
    
    # Normalize the values so that they sum to 1
    normalized_values = random_values / np.sum(random_values)
    
    # Convert the numpy array to a list
    return normalized_values.tolist()

def merge_fastq_files(fastq_file, output_file):
    """
    Merges a FASTQ file into an output FASTQ file using subprocess.call.

    Parameters:
    fastq_file (str): Path to the input FASTQ file.
    output_file (str): Path to the output FASTQ file.
    """
    command = f'cat "{fastq_file}" >> "{output_file}"'
    subprocess.call(command, shell=True)

def create_valid_primer_combinations(df):
    # create df with all the valid amplicon coordinates
    df["valid_combinations"] = ""
    d = pd.DataFrame()
    for i in range(len(df)):
        df["valid_combinations"][i] = evaluate_matches(df["left_primer_loc"][i], df["right_primer_loc"][i])
        for j in range(len(df["valid_combinations"][i])):
            primer_start = df["valid_combinations"][i][j][0]
            primer_end = df["valid_combinations"][i][j][1]
            amplicon_number = df["amplicon_number"][i]
            amplicon_length = primer_end - primer_start
            temp = pd.DataFrame(
                {
                    "amplicon_number": [amplicon_number],
                    'primer_start': [primer_start],
                    'primer_end': [primer_end],
                    'amplicon_length': [amplicon_length]
                }
            )

            d = pd.concat([d, temp])
    all_amplicons = pd.merge(d, df[["amplicon_number","primer_seq_x","primer_seq_y"]], how='outer', sort=False, on='amplicon_number')
    return all_amplicons

def preprocess_primers(primer_file):
    # define column names to read primers bed file
    col_names = ["ref","start", "end", "left_right", "primer_pool","strand", "primer_seq"]
    # read the primer bed file
    primer_bed = pd.read_csv(primer_file, sep= "\t", names=col_names)
        # split the amplicon name into number and left/right
    primer_bed["amplicon_number"] = primer_bed["left_right"].str.split('_').str[1]
    # merge the df with itself to have right and left primer on one row
    df = pd.merge(
        primer_bed.loc[primer_bed["left_right"].str.contains("LEFT")],
        primer_bed.loc[primer_bed["left_right"].str.contains("RIGHT")],
        on=["amplicon_number"]
    )
    # select needed columns
    df = df[["amplicon_number","primer_seq_x","primer_seq_y"]]
    # get complementary reverse sequence of the right primer
    df["comp_rev"] = df.apply(lambda row: Seq(row['primer_seq_y']).reverse_complement(), axis=1)
    return df
    
    


def run_wgsim_on_fasta(fasta_file, output_dir, read_length, error_rate, mutation_rate, outer_distance, read_cnt, indel_fraction, indel_extend_probability):
    """Runs wgsim on a single FASTA file with the given parameters."""
    output_prefix = os.path.splitext(os.path.basename(fasta_file))[0]
    output1 = os.path.join(output_dir, f"{output_prefix}_1.fastq")
    output2 = os.path.join(output_dir, f"{output_prefix}_2.fastq")
    merged_output = os.path.join(output_dir, "merged_reads.fastq")

    command = [
        "wgsim",
        "-e", str(error_rate),
        "-r", str(mutation_rate),
        "-d", str(outer_distance),
        "-N", str(read_cnt),
        "-R", str(indel_fraction),
        "-X", str(indel_extend_probability),
        "-1", str(read_length),
        "-2", str(read_length),
        fasta_file,
        output1,
        output2
    ]

    subprocess.run(command, check=True,capture_output = True,
    text = True)
    command = f'cat "{output1}" >> "{merged_output}"'
    subprocess.call(command, shell=True)
    



def find_closest_primer_match(pattern,reference_seq):
    """function to find a string allowing up to 1 mismatches"""
    matches = [m.start() for m in regex.finditer(r"\L<primer_string>{s<1}",
                                str(reference_seq), primer_string=[pattern])]
    # if the primer not found, try finding it in the complimentary reverse strand
    if len(matches) == 0:
        reference_seq = str(reference_seq.reverse_complement())
        matches = [m.start() for m in regex.finditer(r"\L<primer_string>{s<1}",
                                                     reference_seq, primer_string=[pattern])]
        return matches
    else:
        return matches
    
    
def make_amplicon(left_primer_loc,right_primer_loc,
                  primer_seq_y, reference):
    """function to create amplicons based on the string match location"""
    if math.isnan(left_primer_loc) or math.isnan(right_primer_loc):
        # if there is either no left and right primer match
        amplicon = ""
        # length of the right primer is added to include the 
        # right primer in the amplicon
    else:
        amplicon = str(reference[int(left_primer_loc - 1): int(right_primer_loc + len(primer_seq_y))])
    return amplicon


def evaluate_matches(left_primer_coordinates, right_primer_coordinates):
    """function to evaluate which coordinates found for each primer makes a valid amplicon"""
    # if both left and right primer string matches exist,
    # find left and right primer pairs that can create an amplicon
    if len(left_primer_coordinates) !=0 and len(right_primer_coordinates)!=0:
        valid_combinations = []
        combinations = list(itertools.product(left_primer_coordinates, right_primer_coordinates))
        for combination in combinations:
            amplicon_length = combination[1] - combination[0]
            if 0 < amplicon_length <= 2000:
                valid_combinations.append(combination)
            else:
                pass
        return valid_combinations
    else:
        return []
    

def write_fasta_group(group, amplicon_number, output_dir):

    fasta_filename = os.path.join(output_dir, f'amplicon_{amplicon_number}.fasta')

    filtered_records = [
         SeqRecord(Seq(seq), id=f"{amplicon_number}_{i}", description="")
         for i, seq in enumerate(group['amplicon_sequence'])
         if 1 < len(seq) < 10000
    ]

    if filtered_records:
        SeqIO.write(filtered_records, fasta_filename, 'fasta')
    else:
        print(f"No valid sequences for amplicon {amplicon_number}, no file written.\
             please checkout the amplicon_stats.csv file for more information.")