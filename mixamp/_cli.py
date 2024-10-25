import os
from Bio import SeqIO
import click
from tqdm import tqdm


@click.group(context_settings={'show_default': True})
@click.version_option('1.0.0')
def cli():
    pass


@cli.command()
@click.argument('genomes', type=str)
@click.argument('primers', type=click.Path(exists=True))
@click.option(
    '--proportions', default="NA", type=str,
    help=(
        'Read proportions for each sample, e.g.(0.8,0.2) must sum to 1.0. '
        'If not provided, the program will randomly assign proportions'
    )
)
@click.option(
    '--outdir', default='results', type=click.Path(exists=False),
    help='Output directory', show_default=True
)
@click.option(
    '--simulator', default="wgsim", type=click.Choice(['wgsim', 'mason_simulator'], case_sensitive=False),
    help='Select the simulator to use (wgsim or mason_simulator)'
)
@click.option(
    '--outerdistance', default=150,
    help='Outer distance for simulation'
)
@click.option(
    '--readcnt', default=500,
    help='Number of reads per amplicon.'
)
@click.option(
    '--read_length', default=150,
    help='Read length for simulation'
)
@click.option(
    '--error_rate', default=0, type=float,
    show_default=True,
    help='Base error rate (e.g., 0.02) for simulation'
)
@click.option(
    '--mutation_rate', default=0, type=float,
    show_default=True,
    help='Mutation rate (e.g., 0.001) for simulation'
)
@click.option(
    '--indel_fraction', default=0, type=float,
    show_default=True,
    help='Fraction of indels (e.g., 0.15) for simulation'
)
@click.option(
    '--indel_extend_probability', default=0, type=float,
    show_default=True,
    help='Probability an indel is extended (e.g., 0.3) for simulation'
)
@click.option(
    '--maxmismatch', default=1, show_default=True,
    help='Maximum number of mismatches allowed in primer region'
)
@click.option(
    '--haplotype', is_flag=True, default=True,
    help='use this to simulate reads for a haploid organism.'
)

def simulate_proportions(
    genomes, proportions, primers, outdir, read_length,
    error_rate, mutation_rate, outerdistance, readcnt,
    indel_fraction, indel_extend_probability, maxmismatch,
    haplotype,simulator
):
    from mixamp.utils import (
        preprocess_primers, create_valid_primer_combinations, make_amplicon,
        write_fasta_group, run_simulation_on_fasta, merge_fastq_files,
        find_closest_primer_match, generate_random_values,
    )

    os.makedirs(outdir)
    sample_names = [fp.split("/")[-1].split(".")[0]
                    for fp in str(genomes).split(",")]
    sample_paths = str(genomes).split(",")

    if proportions == "NA":
        print("Read simulation proportions not provided.\
            Will generate proportions randomly...")
        if len(sample_names) == 1:
            print("Only one sample provided...")
            proportions = [1]
        else:
            proportions = generate_random_values(len(sample_names))
            with open(os.path.join(outdir, "sample_proportions.txt"),
                      'w') as file:
                for name, proportion in zip(sample_names, proportions):
                    file.write(f"{name}: {proportion}\n")
    else:
        proportions = list(map(float, str(proportions).split(",")))

    print("Reading and preprocessing the primer file...")
    df = preprocess_primers(primers)

    if len(sample_names) != len(proportions) != len(sample_paths):
        raise Exception("Number of samples and proportions should match!")
    if sum(proportions) != 1.0:
        raise Exception("Sum of all proportions should equal to 1.0!")
    read_cnts = [i * int(readcnt) for i in proportions]

    for name, path, cnt in zip(sample_names,
                                 sample_paths,read_cnts):
        genome_seq = next(SeqIO.parse(path, "fasta"))
        print(f"Extracting amplicons for sample {name}...")

        df["left_primer_loc"] = df.apply(
            lambda row: find_closest_primer_match(str(row["primer_seq_x"]),
                                                  str(genome_seq.seq),
                                                  maxmismatch),
            axis=1
        )
        df["right_primer_loc"] = df.apply(
            lambda row: find_closest_primer_match(str(row["comp_rev"]),
                                                  str(genome_seq.seq),
                                                  maxmismatch),
            axis=1
        )
        all_amplicons = create_valid_primer_combinations(df)
        all_amplicons = all_amplicons.fillna(0)

        os.makedirs(os.path.join(outdir, name, "amplicons"))
        all_amplicons.to_csv(os.path.join(outdir,
                                          name,
                                          "amplicons/amplicon_stats.csv"))

        all_amplicons["amplicon_sequence"] = all_amplicons.apply(
            lambda row: make_amplicon(row["primer_start"],
                                      row["primer_end"],
                                      row["primer_seq_y"],
                                      genome_seq.seq),
            axis=1
        )

        for amplicon_number, group in all_amplicons.groupby('amplicon_number'):
            fasta_file = write_fasta_group(group,
                                           amplicon_number,
                                           os.path.join(outdir,
                                                        name,
                                                        "amplicons"))

        print("Starting read simulation...")
        if not os.path.exists(os.path.join(outdir, name, "reads")):
            os.makedirs(os.path.join(outdir, name, "reads"))

        fasta_files = [
            os.path.join(outdir, name, "amplicons", f)
            for f in os.listdir(os.path.join(outdir, name, "amplicons"))
            if f.endswith(".fasta") or f.endswith(".fa")
        ]

        with tqdm(total=len(fasta_files),
                  desc="Processing FASTA files") as pbar:
            for fasta_file in fasta_files:
                run_simulation_on_fasta(
                    fasta_file, os.path.join(outdir, name, "reads"),
                    read_length, error_rate, mutation_rate, outerdistance,
                    cnt, indel_fraction, indel_extend_probability, haplotype,simulator
                )
                pbar.update(1)

        read_path1 = os.path.join(os.path.abspath(outdir),
                                 name, "reads/merged_reads_1.fastq")
        read_path2 = os.path.join(os.path.abspath(outdir),
                                 name, "reads/merged_reads_2.fastq")
        
        output_path1 = os.path.join(os.path.abspath(outdir), "reads_1.fastq")
        output_path2 = os.path.join(os.path.abspath(outdir), "reads_2.fastq")
        print("Merging all reads...")
        merge_fastq_files(read_path1, output_path1)
        merge_fastq_files(read_path2, output_path2)
        print("Finished!")


if __name__ == '__main__':
    cli()
