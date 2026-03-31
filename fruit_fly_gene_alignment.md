# Testing align.py with real genes

Thomas Hunt Morgan ran experiments on fruit fly to study genetics. Morgan identified the mutant white-eyed gene by comparing it to the wildtype red-eyed gene population.

The `test_files` folder contains the nucleotide sequences of these 2 alleles.

Code output:

```
py align.py test_files\wildtype.fasta test_files\mutant.fasta
alignment 1: GTTTCGTGACGAAGCTCCAAGCGGTTTACGCCATCAATTAAACACAAAGTGCTGTGCC--AA--------
alignment 2: GTTTCGTGACGAAGCTCCAAGCGGTTTACGCCATCAATTAAACACAAAGTGCTGTGCCAAAACTCCTCTC
optimal score: 50.0
```

This indicates that the white-eyed mutant allele is an insertion mutation, adding on to the end of the wildtype allele. 