from pathlib import Path
import numpy as np

def parse_FASTA(file_path):
    """
    Returns sequence of nucleotides after reading FASTA file
    """

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except FileNotFoundError:
        print("Error: The file was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    print(content)

    lines = content.split('\n')

    dna_sequences = []
    for line in lines:
        if line and line[0] != '>':
            dna_sequences.append(line)

    return dna_sequences


def check(a, b):
    match = 1
    mismatch = -1

    if a == b:
        return match
    else:
        return mismatch


def needleman_wunsch(a, b):
    '''
    Returns optimal alignment of 2 DNA sequences using Needleman-Wunsch algorithm
    '''

    d = -1 # gap penalty score

    num_rows = len(a) + 1
    num_cols = len(b) + 1

    dp_table = np.zeros((num_rows, num_cols))

    # initialize first column 
    dp_table[:, 0] = -1 * np.arange(0, num_rows)

    # initialize first row
    dp_table[0, :] = -1 * np.arange(0, num_cols)

    # fill in dp table
    for i in range(1, num_rows):
        for j in range(1, num_cols):
            dp_table[i, j] = max(
                dp_table[i-1, j-1] + check(a[i-1], b[j-1]), # diagonal upper left box + maybe_match
                dp_table[i, j-1] + d, # left box + indel; insert
                dp_table[i-1, j] + d # upper box + indel; delete
            )

    # search dp table
    alignment_A = ''
    alignment_B = ''

    # start from bottom right. remember that dp table has padding
    i = len(a) 
    j = len(b) 

    while (i > 0 or j > 0):
        # remember that strings don't have padding
        if (i > 0 and j > 0 and dp_table[i, j] == dp_table[i-1, j-1] + check(a[i-1], b[j-1])):
            # i should go up left
            alignment_A = a[i-1] + alignment_A
            alignment_B = b[j-1] + alignment_B

            # update position
            i -= 1
            j -= 1
        elif (i > 0 and dp_table[i, j] == dp_table[i-1, j] + d):
            # go up
            alignment_A = a[i-1] + alignment_A
            alignment_B = '-' + alignment_B

            # update position
            i -= 1
        else:
            # go left
            alignment_A = '-' + alignment_A
            alignment_B = b[j-1] + alignment_B

            # update position
            j -= 1

    optimal_score = dp_table[-1, -1]
    return alignment_A, alignment_B, optimal_score
    

def main():
    BASE_DIR = Path(__file__).resolve().parent # gets directory file is running in
    file_path = BASE_DIR / "test.fasta"
    dna_sequences = parse_FASTA(file_path)

    dna_1, dna_2 = dna_sequences
    
    alignment_A, alignment_B, optimal_score = needleman_wunsch(dna_1, dna_2)

    print('alignment A: ' + alignment_A)
    print('alignment B: ' + alignment_B)
    print('optimal score: ' + str(optimal_score))

main()