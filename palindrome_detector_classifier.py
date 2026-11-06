# =============================================================================
# Minimal script to take input of 1) a nucleotide sequence, and 2) a dot-
# bracket string describing its secondary structure. The code returns a list
# of detected palindromes, as well as their per-residue structural classif-
# ication based on the work from Ahn et al.
# =============================================================================

import numpy as np
import re

def is_valid_pair(base1, base2):
    pairs = {'A': ['T', 'U'], 'T': ['A', 'G'], 'U': ['A', 'G'], 'C': ['G'], 'G': ['C', 'U', 'T']}
    return base2 in pairs.get(base1, '')

def find_rna_palindromes(rna_sequence):
    global palindromes
    palindromes = []
    length_list = []
    sequence_length = len(rna_sequence)

    for i in range(sequence_length - 1):
        if is_valid_pair(rna_sequence[i], rna_sequence[i + 1]):
            length_list.append(2)
            palindrome_found = True
            start = i
            end = i + 1

            while palindrome_found:
                start -= 1
                end += 1
                if start < 0 or end > sequence_length - 1:
                    break
                if is_valid_pair(rna_sequence[start], rna_sequence[end]):
                    palindromes.append((start, end + 1, rna_sequence[start:end + 1]))
                    length_list.append(end - start + 1)
                else:
                    palindrome_found = False

    length_hist = np.histogram(length_list, bins=np.arange(1, np.array(length_list).max() + 1, 2))
    print((np.array(length_list) == 4).sum())
    print(palindromes)
    return palindromes, length_hist

def classify_sequence(dot_bracket, palindromes):
    classifications = []
    total_class_1 = 0
    total_class_2 = 0
    total_class_3 = 0

    for start, stop, seq_fragment in palindromes:
        subseq = dot_bracket[start:stop]
        # Initialize per_res_class with the length of subseq, not the palindrome
        per_res_class = np.zeros(len(subseq), dtype=int)

        # Set everything to class 2 if a bracket, and class 3 if a dot
        per_res_class[:] = 2
        per_res_class[np.array(list(subseq)) == '.'] = 3

        # Check if it's a perfect stem-loop (class 1)
        matches = re.finditer(r'\(\.+\)', subseq)
        for match in matches:
            startstem, endstem = match.span()

            # Extend startstem and endstem until you no longer have base pairs
            while startstem > 0 and endstem < len(subseq) and subseq[startstem - 1] == '(' and subseq[endstem] == ')':
                startstem -= 1
                endstem += 1
            per_res_class[startstem:endstem] = 1

        # Calculate the fraction of residues in each class
        f1 = (per_res_class == 1).sum() / per_res_class.shape[0]
        f2 = (per_res_class == 2).sum() / per_res_class.shape[0]
        f3 = (per_res_class == 3).sum() / per_res_class.shape[0]

        # Set a class for the palindrome to describe the majority of the residues
        if (f1 > f2) and (f1 > f3):
            classification = "Class 1"
        elif f2 > f3:
            classification = "Class 2"
        else:
            classification = "Class 3"

        # Add all data to the master dictionary
        classifications.append({
            'location': (start, stop),
            'sequence': seq_fragment,
            'dot_bracket': subseq,
            'f1': f1,
            'f2': f2,
            'f3': f3,
            'class': classification,
            'per_res_class': per_res_class
        })

        # Update total counts for each class
        total_class_1 += (per_res_class == 1).sum()
        total_class_2 += (per_res_class == 2).sum()
        total_class_3 += (per_res_class == 3).sum()

    return classifications, total_class_1, total_class_2, total_class_3

# Define palindromic residues using a statistical equation for expected values
def palindromeRes(i):
    return 2 * i * (6 / 16) ** i

def expected_palindromic_residues(sequence_length):
    expected = 0
    for i in range(2, sequence_length + 1):
        expected += palindromeRes(i)
    return expected

def process_large_sequence(rna_seq, dot_bracket):
    sequence_length = len(rna_seq)
    all_classifications = []
    dot_bracket_str = "".join(dot_bracket)

    # Find palindromes
    palindromes, length_histogram = find_rna_palindromes(rna_seq)
    print(length_histogram)
    print(len(palindromes))

    # Classify the identified palindromic sequences
    classifications, total_class_1, total_class_2, total_class_3 = classify_sequence(dot_bracket_str, palindromes)

    # Calculate the total number of palindromic residues
    total_palindromic_residues = total_class_1 + total_class_2 + total_class_3

    # Calculate the expected number of palindromic residues
    expected_palindromes = expected_palindromic_residues(sequence_length)

    return classifications, total_class_1, total_class_2, total_class_3, total_palindromic_residues, expected_palindromes

# RNA sequence and dot-bracket notation for example RNA
rna_seq = "GAAAUUGUGCUGAAGCCGGCACUGGAGUCUGGGCGAGGCUGCUUACUUCCAAACGCUCGAUCUAUACCGCUCUUCCUGUAAGGUCUGCAUUAUUAAGGAGUAGUCCCUCUGAUAAUGCUCCUUGAGGGUCUCGAAGG"
dot_bracket = "......((((((....))))))(((((..((((((....)))))).)))))................((....((((.(((((...(((((((((.(((....)))...))))))))).)))))))))...))...."

classifications, total_class_1, total_class_2, total_class_3, total_palindromic_residues, expected_palindromes = process_large_sequence(rna_seq, dot_bracket)

# Output the classifications
# Uncomment the following lines if you want to see all palindromes and their per-residue classes.
"""
for classification in classifications:
    loc = classification['location']
    seq = classification['sequence']
    dot_bracket = classification['dot_bracket']
    f1 = classification['f1']
    f2 = classification['f2']
    f3 = classification['f3']
    cls = classification['class']
    per_res_class = classification['per_res_class']
    print('%s\n%s\n%s' % (seq, dot_bracket, ''.join(per_res_class.astype(str))))
    print('f1 = %.2f f2 = %.2f f3 = %.2f' % (f1, f2, f3))
    print('\n')
"""

# Output the summary of palindrome classes
print("Summary of Palindrome Classes:")
print(f"Sequence length: {len(rna_seq)}")
print(f"Total residues in Class 1: {total_class_1}")
print(f"Total residues in Class 2: {total_class_2}")
print(f"Total residues in Class 3: {total_class_3}")
print(f"Total palindromic residues: {total_palindromic_residues}")
print(f"Expected palindromic residues: {1.17 * len(rna_seq):.2f}")


