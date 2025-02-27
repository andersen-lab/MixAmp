# MixAmp: Amplicon Read Simulator


A tool for Amplicon read simulation for waste water sequencing or other aplications. Users can easily simulate reads from mutiple samples with different proportions using the tool.

## Usage
If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <https://github.com/andersen-lab/MixAmp> repository.

## Installation
Mixamp is written in python 3 but it requires <a href="https://github.com/lh3/wgsim">wgsim</a> and 
<a href="https://github.com/seqan/seqan/blob/main/apps/mason2/README.mason_simulator">mason simulator</a>
to simulate reads.

### Local build from source
```
git clone https://github.com/andersen-lab/MixAmp
cd mixamp
pip install -e .
```
Please note that pip does not install all the requirements,
some packages need to be installed via Conda or be built from source.

### Installing via Conda
1. `pip install git+https://github.com/andersen-lab/MixAmp`
2. Create a conda environment as mixamp and install the dependencies:


```
conda create -n mixamp
conda activate mixamp
conda env update --file environment.yml
```


## Example commands

Run the tool using the following command.
 ```
mixamp simulate-proportions [SAMPLE1.fasta,SAMPLE2.fasta,..] [primer.bed] --proportions [0.8,0.2,..] --outdir [output_directory]
 ```

Simulate reads from different samples without defining proportions (will be assigned randomly, proportions can be found in `results/sample_proportions.txt`) and allowing upto 2 SNPs mistmatches in the primer regions.
 ```
mixamp simulate-proportions sample.fasta,sample2.fasta primer.bed --outdir results/ --maxmismatch 2
 ```
Simulate reads with user-defined proportions and specifing read simulator.
Mixamp uses wgsim as a simulator but you can change it to mason.
 ```
mixamp simulate-proportions sample.fasta,sample2.fasta primer.bed --proportions 0.2,0.8 --simulator mason
 ```
Simulate reads with user-defined proportions and number of reads per amplicon.
 ```
mixamp simulate-proportions sample.fasta,sample2.fasta primer.bed --proportions 0.2,0.8 --readcnt 1000
 ```

Simulate reads with additional parameters such as base error rate, read length and indels fraction
 ```
mixamp simulate-proportions sample.fasta,sample2.fasta primer.bed --proportions 0.2,0.8 --readcnt 1000 --error_rate 0.001 --read_length 400 --indel_fraction 0.001
 ```
## Notes
#### Number of reads per amplicon
It is recommended to define the number of reads per amplicon to be greater than the number of contigs in your amplicon file. This is particularly important when your primers are designed for whole genome sequencing, where each amplicon may contain a substantial number of contigs. Setting too few reads per amplicon may result in empty read files for certain amplicons, leading to incomplete simulated reads.
#### Primer bed file
Please remember that the primer file must contain a column containing primer sequence. The maximum number of mismatches allowed for each primer sequence is 1 SNP. To change this number, you may use the `--maxmismatches` flag.
#### Complete set of available parameters
To learn more about how to adjust other parameters use `mixamp simulate-proportions --help`
#### Simulated reads output
Simulated reads from all samples are located in `provided_output_path/reads.fastq`
#### Information about amplicon dropouts
In order to find more about amplicon dropouts, please refer to `provided_output_path/sample_name/amplicon_stats.csv` file. This file will have right/left primer matching coordinates as zero if no matches found.