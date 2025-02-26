import unittest
import os

def file_exists(directory, filename):
    file_path = os.path.join(directory, filename)
    return os.path.exists(file_path)


class CommandLineTests(unittest.TestCase):
    def test_version(self):
        os.system("mixamp --version")

    def test_simulation(self):
        os.system(
            "mixamp simulate-proportions mixamp/tests/data/ATM-2FFMD73N3.fasta \
            mixamp/tests/data/ARTIC_V4-1.bed"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_simulator(self):
        os.system(
            "mixamp simulate-proportions mixamp/tests/data/ATM-2FFMD73N3.fasta \
            mixamp/tests/data/ARTIC_V4-1.bed --simulator wgsim"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_max_mismatch(self):
        os.system(
            "mixamp simulate-proportions mixamp/tests/data/ATM-2FFMD73N3.fasta \
            mixamp/tests/data/ARTIC_V4-1.bed --maxmismatch 2"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_readcnt(self):
        os.system(
            "mixamp simulate-proportions mixamp/tests/data/ATM-2FFMD73N3.fasta \
            mixamp/tests/data/ARTIC_V4-1.bed --readcnt 200"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_readlength(self):
        os.system(
            "mixamp simulate-proportions mixamp/tests/data/ATM-2FFMD73N3.fasta \
            mixamp/tests/data/ARTIC_V4-1.bed --readlength 130"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))
    def test_simulation_with_2genomes(self):
        os.system(
            "mixamp simulate-proportions mixamp/tests/data/ATM-2FFMD73N3.fasta,\
                mixamp/tests/data/KR-SEARCH-120354.fasta \
            mixamp/tests/data/ARTIC_V4-1.bed"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))
    def test_simulation_with_2genomes_proportions(self):
        os.system(
            "mixamp simulate-proportions mixamp/tests/data/ATM-2FFMD73N3.fasta,\
                mixamp/tests/data/KR-SEARCH-120354.fasta \
            mixamp/tests/data/ARTIC_V4-1.bed --proportions 0.8,0.2"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))


if __name__ == "__main__":
    unittest.main()

