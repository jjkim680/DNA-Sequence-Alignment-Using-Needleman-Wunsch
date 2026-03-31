# Custom Genomic Sequence Aligner

## Overview
A Python-based bioinformatics tool engineered from scratch to parse raw FASTA files and process high-throughput genetic sequence data. This project computes the optimal global alignment between complex DNA sequences using the Needleman-Wunsch Algorithm, optimized with NumPy matrix operations.

## Features
* **Custom FASTA Parser:** Extracts sequence data and metadata directly from raw `.fasta` files without relying on external bioinformatics libraries (like Biopython).
* **Global Alignment:** Implements the classic Needleman-Wunsch algorithm to find the mathematically optimal alignment across the entire length of two sequences.
* **Optimized Traceback:** Leverages NumPy vectorization and matrix operations to efficiently handle gap penalty calculations and traceback paths.

## Algorithm & Mathematical Formulation

### The Needleman-Wunsch Algorithm
The core of the aligner relies on dynamic programming to build a scoring matrix $F$. The algorithm calculates the optimal alignment score for sequences of length $m$ and $n$ by evaluating matches, mismatches, and gap penalties recursively. 

The scoring matrix is populated using the following recurrence relation:

$$F(i,j) = \max \begin{cases} F(i-1, j-1) + S(x_i, y_j) & \text{(Match/Mismatch)} \\ F(i, j-1) + d & \text{(Gap in Sequence 1)} \\ F(i-1, j) + d & \text{(Gap in Sequence 2)} \end{cases}$$

Where $S$ is the substitution matrix score and $d$ is the linear gap penalty.

### Complexity
By utilizing NumPy to manage the 2D matrix allocations, the algorithm strictly adheres to constraints of **$O(mn)$** time complexity and **$O(mn)$** space complexity, ensuring scalability for larger genetic sequences.

## Installation & Usage

Clone the repository and navigate to the directory:

```bash
git clone [https://github.com/jjkim680/DNA-Sequence-Alignment-Using-Needleman-Wunsch.git](https://github.com/jjkim680/DNA-Sequence-Alignment-Using-Needleman-Wunsch.git)
cd DNA-Sequence-Alignment-Using-Needleman-Wunsch
```

Run the aligner from the command line, passing in the file paths of your two FASTA files as arguments:

```python
py align.py sequence1.fasta sequence2.fasta
```

Example output code for `test_wildtype.fasta` and `test_mutated.fasta`, which are provided in the `test_files` folder:

```python
alignment 1: ATGCGTACGTAGC
alignment 2: ATGCAT-CGTAGC
optimal score: 9.0
```

There are 11 matches, 1 mismatch, and 1 indel, so the optimal score is $11 - 2 = 9$

