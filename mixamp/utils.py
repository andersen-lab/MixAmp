import math
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess
import pandas as pd
import os
import itertools
import numpy as np
from tqdm import tqdm
import regex as re


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
import pandas as pd

def create_valid_primer_combinations(df):
    # create df with all the valid amplicon coordinates
    df["valid_combinations"] = ""
    d = pd.DataFrame()  # Initialize empty dataframe
    for i in range(len(df)):
        df["valid_combinations"][i] = evaluate_matches(df["left_primer_loc"][i], df["right_primer_loc"][i])
        for j in range(len(df["valid_combinations"][i])):
            primer_start = df["valid_combinations"][i][j][0]
            primer_end = df["valid_combinations"][i][j][1]
            amplicon_number = df["amplicon_number"][i]
            temp = pd.DataFrame(
                {
                    "amplicon_number": [amplicon_number],
                    'primer_start': [primer_start],
                    'primer_end': [primer_end]
                }
            )
            d = pd.concat([d, temp], ignore_index=True)    # Check if the dataframe `d` is empty and exit with an error message if true
    if d.empty:
        raise ValueError("No primer matches found, please check your primer file")
    else:
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
        on=["amplicon_number","primer_pool"]
    )
    # select needed columns
    df = df[["amplicon_number","primer_seq_x","primer_seq_y"]]
    # get complementary reverse sequence of the right primer
    df["comp_rev"] = df.apply(lambda row: Seq(row['primer_seq_y']).reverse_complement(), axis=1)
    mask = df.duplicated(subset=['amplicon_number'], keep=False)
    # Apply a function to append an index for duplicated amplicon_number values
    df.loc[mask, 'amplicon_number'] = df.loc[mask, 'amplicon_number'].astype(str) + "_" +  df.loc[mask].groupby('amplicon_number').cumcount().astype(str) 
    return df
    
    

def count_contigs(fasta_file):
    """Counts the number of contigs in the FASTA file."""
    contig_count = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                contig_count += 1
    return contig_count

def run_wgsim_on_fasta(fasta_file, output_dir, read_length, error_rate, mutation_rate, outer_distance, read_cnt, indel_fraction, indel_extend_probability, haplotype):
    """Runs wgsim on a single FASTA file with the given parameters, simulating read_cnt / number_of_contigs reads per contig."""
    
    # Count the number of contigs in the FASTA file
    num_contigs = count_contigs(fasta_file)
    
    if num_contigs == 0:
        raise ValueError("No contigs found in the FASTA file.")
    
    # Calculate the number of reads per contig
    reads_per_contig = read_cnt // num_contigs  # Integer division
    
    # Prepare output filenames
    output_prefix = os.path.splitext(os.path.basename(fasta_file))[0]
    merged_output1 = os.path.join(output_dir, "merged_reads_1.fastq")
    merged_output2 = os.path.join(output_dir, "merged_reads_2.fastq")
    
    # Initialize merged output files
    with open(merged_output1, 'a') as f1, open(merged_output2, 'a') as f2:
        pass  # This creates empty files

    # Loop through the contigs and run wgsim for each
    for contig_idx in range(num_contigs):
        # Temporary output files for each contig
        output1 = os.path.join(output_dir, f"{output_prefix}_contig{contig_idx + 1}_1.fastq")
        output2 = os.path.join(output_dir, f"{output_prefix}_contig{contig_idx + 1}_2.fastq")
        
        # Build wgsim command for this contig
        command = [
            "wgsim",
            "-e", str(error_rate),
            "-r", str(mutation_rate),
            "-d", str(outer_distance),
            "-N", str(reads_per_contig),
            "-R", str(indel_fraction),
            "-X", str(indel_extend_probability),
            "-1", str(read_length),
            "-2", str(read_length),
            fasta_file,
            output1,
            output2
        ]

        # Add the "-h" flag if haplotype is True
        if haplotype:
            command.append("-h")

        # Run the wgsim command
        subprocess.run(command, check=True, capture_output=True, text=True)
        
        # Merge the contig-specific output into the final merged output files
        command_merge = f'cat "{output1}" >> "{merged_output1}" && cat "{output2}" >> "{merged_output2}"'
        subprocess.call(command_merge, shell=True)
    
    

def find_closest_primer_match(pattern,reference_seq,maxmismatch):
    """function to find a string allowing up to 1 mismatches"""
    # Define the fuzzy regex pattern with a maximum number of mismatches (substitutions)
    primer_pattern = f"({pattern}){{s<={maxmismatch}}}"
    
    matches = [match.start() for match in re.finditer(primer_pattern, reference_seq, re.IGNORECASE)]

    # if the primer not found, try finding it in the complimentary reverse strand
    if len(matches) == 0:
        matches = [match.start() for match in re.finditer(primer_pattern, str(Seq(reference_seq).reverse_complement()), re.IGNORECASE)]
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