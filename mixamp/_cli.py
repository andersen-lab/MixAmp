import os

from Bio import SeqIO
import click
from tqdm import tqdm

@click.group(context_settings={'show_default': True})
@click.version_option('1.0.0')
def cli():
    pass

@cli.command()
@click.argument('genomes',type=str)
@click.argument('primers', type=click.Path(exists=True))
@click.option('--proportions', default="NA",
              help='Read proportions for each sample, e.g.(0.8,0.2)\
                  must sum to 1.0, if not provided, program will\
                      randomly assign proportions', type = str)
@click.option('--outdir', default='results',
              help='Output directory',
              type=click.Path(exists=False), show_default=True)
@click.option('--outerdistance', default=100,
              help='outerdistance for simulation')
@click.option('--readcnt', default=100000,
              help='Total number of reads to be simulated')
@click.option('--read_length', default=250,
              help='Read length for simulation')
@click.option('--error_rate', default=0,
              help='Base error rate (e.g., 0.02) for simulation',
              show_default=True)
@click.option('--mutation_rate', default=0,
              help='Mutation rate (e.g., 0.001) for simulation',
              show_default=True)
@click.option('--indel_fraction', default=0,
              help='Fraction of indels (e.g., 0.15) for simulation', 
              show_default=True)
@click.option('--indel_extend_probability', default=0,
              help='Probability an indel is extended (e.g., 0.3) for simulation',
              show_default=True)
  
def simulate_proportions(genomes,proportions,primers,outdir,read_length,
                               error_rate, mutation_rate,outerdistance,
                               readcnt,indel_fraction, indel_extend_probability):
    from mixamp.utils import (preprocess_primers,
                       create_valid_primer_combinations,
                       make_amplicon,write_fasta_group,
                       run_wgsim_on_fasta,merge_fastq_files,
                       find_closest_primer_match, generate_random_values)
    # the sequence used to create amplicons
    os.makedirs(outdir)
    sample_names = [fp.split("/")[-1].split(".")[0] for fp in str(genomes).split(",")]
    sample_paths = str(genomes).split(",")
    if proportions == "NA":
        print("Read simulation proportions not provided."
            "Will generate proportions randomly..")
        if len(sample_names) == 1:
            print("Only one sample provided...")
            proportions = [1]
        else:
            proportions = generate_random_values(len(sample_names))
            with open(os.path.join(outdir,"sample_proportions.txt"), 'w') as file:
                for name, proportion in zip(sample_names, proportions):
                    file.write(f"{name}: {proportion}\n")

            
    else:       
        proportions = list(map(float, str(proportions).split(",")))
    print("Reading and preprocessing the primer file...")
    df = preprocess_primers(primers)
    # error if the proportion counts do not match the sample counts
    if len(sample_names) != len(proportions) != len(sample_paths) != len(sample_paths):
        raise Exception("Number of samples and proportions should match!")
    # error if proportions do not add up to 1.0
    if sum(proportions) != 1.0:
        raise Exception("Sum of all proportions should equal to 1.0!")
    # Number of distinct amplicons
    amplicon_cnt = int(df['amplicon_number'].nunique())
    # Number of reads per sample
    read_cnts = [i * int(readcnt) for i in proportions]
    # Number of reads per amplicon
    read_cnts_per_amp = [int(i /amplicon_cnt) for i in read_cnts]
    # loop over each isolate and create reads
    for name,path,count in zip(sample_names,sample_paths, read_cnts_per_amp):
            # find the left and right primer
        genome_seq = next(SeqIO.parse(path, "fasta"))
        print(f"extracting amplicons for sample {name}...")
        df["left_primer_loc"] = df.apply(lambda row: find_closest_primer_match(str(row["primer_seq_x"]),genome_seq.seq), axis=1)
        df["right_primer_loc"] = df.apply(lambda row: find_closest_primer_match(str(row["comp_rev"]),genome_seq.seq), axis=1)
        all_amplicons = create_valid_primer_combinations(df)
        # If there is no primer match, or the length of amplicon
        # substitute 0 instead
        all_amplicons = all_amplicons.fillna(0)
        # Create a folder for amplicons
        os.makedirs(os.path.join(outdir,name,"amplicons"))
        # write df containing information about amplicons
        all_amplicons.to_csv(os.path.join(outdir,name,"amplicons/amplicon_stats.csv"))
        # create amplicon sequence
        all_amplicons["amplicon_sequence"] = all_amplicons.apply(lambda row: make_amplicon(row["primer_start"],
                                                                          row["primer_end"],
                                                                          row["primer_seq_y"],genome_seq.seq), axis=1)
        # write amplicons
        for amplicon_number, group in all_amplicons.groupby('amplicon_number'):
            fasta_file = write_fasta_group(group, amplicon_number, os.path.join(outdir,name,"amplicons"))
                        
        print("Starting read simulation...")
        # Iterate over all files in the directory
        if not os.path.exists(os.path.join(outdir,name,"reads")):
            os.makedirs(os.path.join(outdir,name,"reads"))
        fasta_files = [os.path.join(outdir, name, "amplicons", f) for f in os.listdir(os.path.join(outdir, name, "amplicons")) if f.endswith(".fasta") or f.endswith(".fa")]
        # Use tqdm to track progress
        with tqdm(total=len(fasta_files), desc="Processing FASTA files") as pbar:
            for fasta_file in fasta_files:
                run_wgsim_on_fasta(fasta_file, os.path.join(outdir, name, "reads"),read_length,
                        error_rate, mutation_rate,outerdistance,
                        readcnt,indel_fraction, indel_extend_probability)
                pbar.update(1)
        # create paths to merge all reads after simualtion
        read_path = os.path.join(os.path.abspath(outdir),name,"reads/merged_reads.fastq")
        output_path = os.path.join(os.path.abspath(os.path.abspath(outdir)),"reads.fastq")
        # merge fastq files 
        print("Merging all reads...")
        merge_fastq_files(read_path,output_path)
        print("Finished!")
        
if __name__ == '__main__':
    cli()