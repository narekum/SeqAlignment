# Pairwise sequence alignment tools

Implementation of dynamic programming algorithms in python. 

## NeedlemanWunsch.py

This script uses [Needleman and Wunsch](https://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm) algorithm for performing global alignment of two sequences.

Usage:
```sh
$ ./Needleman-Wunsch.py --help
```
```
usage: NeedlemanWunsch.py [-h] [-f] [-g GAP] [-m {blosum62,blosum45,blosum80}]
                          seq1 seq2

Needleman and Wunsch Algorithmn for Global Sequence Alignment

positional arguments:
  seq1                  First sequeunce
  seq2                  Second sequence

optional arguments:
  -h, --help            show this help message and exit
  -f, --fasta           seq1 and seq2 are in fasta formatted file
  -g GAP                GAP is gap penalty (default: -4)
  -m {blosum62,blosum45,blosum80}
                        Substitution matrix: default is blosum62
```
## SmithWaterman.py

This script uses [Smith and Waterman](https://en.wikipedia.org/wiki/Smith-Waterman_algorithm) algorithm for performing local alignment of two sequences.

usage:
```sh
$ ./SmithWaterman.py
```
```
usage: SmithWaterman.py [-h] [-f] [-g GAP] [-m {blosum62,blosum45,blosum80}]
                        seq1 seq2

Needleman and Wunsch Algorithmn for Global Sequence Alignment

positional arguments:
  seq1                  First sequeunce
  seq2                  Second sequence

optional arguments:
  -h, --help            show this help message and exit
  -f, --fasta           seq1 and seq2 are in fasta formatted file
  -g GAP                GAP is gap penalty (default: -4)
  -m {blosum62,blosum45,blosum80}
                        Substitution matrix: default is blosum62
```

