
####Goal: Implement Needleman-Wunsch algorithm for global alignment of DNA sequences


### Description 
This program finds the best pairwise alignment of sequences in a FASTA file using the Needleman-Wunsch algorithm for global alignment of DNA sequences. It reports the best alignment (reports only one alignment if there are multiple "bests"), the score, and the FASTA description line for the two sequences.

### How to use program
```$ perl align.plx [-o 1] <input_file>```

##### Options #####
[-o 1] is option to output the dynamic programming table and traceback matrix. If [-o] is not included or not followed by a 1, then default is to NOT print the dynamic programming table and traceback matrix

##### Input & Output #####
Input file should be FASTA format and output is to standard out. It's best to redirect the standard out to a file if you want to print the dynamic programing tables, especially with long sequences, because they don't fit well in the terminal window.

##### Example #####

```$ perl align.plx hw1_sequence.fasta```

```
Best Alignment Score:	66
Sequence 1:	GCTTTTTATGAACAA-AATTATAGACATTTTAGTTCTTA-TAATAAATAATAGATATTAAAGAAAATAAAAAAATAGAAATA
Sequence 2:	AATATC-ATAACCCTTGATAACCCAGAAATTAATACTTAATCAAAAATGA-AAATATTAATTAATAAAAGTGAATTGAATAA
Sequence 1 ID: 	>seq1
Sequence 2 ID:	>seq2
```
