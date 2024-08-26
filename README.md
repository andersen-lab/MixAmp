# MixAmp: Amplicon Read Simulator


A tool for Amplicon read simulation for waste water sequencing or other aplications. Users can easily simulate reads from mutiple samples with different proportions using the tool.

## Usage
If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <https://github.com/andersen-lab/MixAmp> repository.

## Installation

Please follow these steps to install required dependencies:

1. `pip install git+https://github.com/andersen-lab/MixAmp`
2. Install dependencies in an activated conda environment:

`conda env update --file environment.yml`


## Example commands

* Run the tool using the following command.
 ```
    mixamp simulate-proportions [SAMPLE1.fasta,SAMPLE2.fasta,..] [primer.bed] --proportions [0.8,0.2,..] --outdir [output_directory]
 ```
Imagine we have two different samples with their whole genome fasta files and we are trying to simulate reads from these two.

* Simulate reads without defining proportions (will be assigned randomly, proportions can be found in `results/sample_proportions.txt`)
 ```
    mixamp simulate-proportions sample.fasta,sample2.fasta primer.bed --outdir results/
 ```
* Simulate reads with user-defined proportions.
 ```
    mixamp simulate-proportions sample.fasta,sample2.fasta primer.bed --proportions 0.2,0.8
 ```
* Simulate reads with user-defined proportions and total number of reads.
 ```
    mixamp simulate-proportions sample.fasta,sample2.fasta primer.bed --proportions 0.2,0.8 --readcnt 20000
 ```

 * Simulate reads with additional parameters such as base error rate, read length and indels fraction
 ```
    mixamp simulate-proportions sample.fasta,sample2.fasta primer.bed --proportions 0.2,0.8 --readcnt 20000 --error_rate 0.001 --read_length 400 --indel_fraction 0.001
 ```
### Notes
* Please remember that the primer file must contain a column containing primer sequence. The maximum number of mismatches allowed for each primer sequence is 1 SNP. To change this number, you may use the `--maxmismatches` flag.

* To learn more about how to adjust other parameters use `mixamp simulate-proportions --help`
* Simulated reads from all samples are located in `provided_output_path/reads.fastq`
* The pipeline will automatically generate random proportions if not provided by the user.
* In order to find more about amplicon dropouts, please refer to `provided_output_path/sample_name/amplicon_stats.csv` file. 
* If there are no right or left primer match, the start/end and length of the amplicon will equal to zero.