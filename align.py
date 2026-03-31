import numpy as np
import argparse

def parse_fasta(file_path):
    """
    Returns sequence of nucleotides after reading FASTA file
    """

    # open file if possible
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except FileNotFoundError:
        print("Error: The file was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    # split file content into lines
    lines = content.split('\n')

    # return dna sequence, not header
    for line in lines:
        if line and line[0] != '>':
            return line

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
    dp_table[:, 0] = d * np.arange(0, num_rows)

    # initialize first row
    dp_table[0, :] = d * np.arange(0, num_cols)

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

    # retrace optimal path
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
        elif (j > 0 and dp_table[i, j] == dp_table[i, j-1] + d):
            # go left
            alignment_A = '-' + alignment_A
            alignment_B = b[j-1] + alignment_B
            j -= 1

    optimal_score = dp_table[-1, -1]
    return alignment_A, alignment_B, optimal_score
    

def main():
    # initialize parser
    parser = argparse.ArgumentParser(description="A custom Needleman-Wunsch global sequence aligner.")
    # positional arguments
    parser.add_argument("sequence1", help="First FASTA file path")
    parser.add_argument("sequence2", help="Second FASTA file path")

    args = parser.parse_args()
    file1_path = args.sequence1
    file2_path = args.sequence2

    # parse the fasta files
    sequence1 = parse_fasta(file1_path)
    sequence2 = parse_fasta(file2_path)
    
    alignment_A, alignment_B, optimal_score = needleman_wunsch(sequence1, sequence2)

    print('alignment 1: ' + alignment_A)
    print('alignment 2: ' + alignment_B)
    print('optimal score: ' + str(optimal_score))

if __name__ == "__main__":
    main()