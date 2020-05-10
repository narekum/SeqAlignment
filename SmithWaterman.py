#!/usr/bin/python3

################################################################################
# Program name      : SmithWaterman.py
# Original Author   : Narendra Kumar, narekum@gmail.com
# Date created      : 03 May 2020
# Purpose           : Smith-Waterman Algorithm for local sequence alignment
#                     https://en.wikipedia.org/wiki/Smith-Waterman_algorithm
#                     Uses linear gap penalty
################################################################################

import argparse
from lib.substitution_matrixes import *
from lib.align import N_W, S_W, PAIR_ALIGNED
from lib.readfasta import readSeqFromFastaFile

def parse_args():
    parser = argparse.ArgumentParser(
                        description="Needleman and Wunsch Algorithmn for Global Sequence Alignment",
                        #epilog="Cite: narekum@gmail.com if you use"
    )

    parser.add_argument("seq1", help="First sequeunce")
    parser.add_argument("seq2", help="Second sequence")
    parser.add_argument("-f","--fasta", help="seq1 and seq2 are in fasta formatted file", action="store_true")
    parser.add_argument('-g', dest="GAP",type=int, default=-4, help="GAP is gap penalty (default: -4)")
    parser.add_argument('-m', dest="substitution_matrix", default='blosum62', choices=['blosum62','blosum45','blosum80'], help="Substitution matrix: default is blosum62")
    return parser.parse_args()


def do_stuff(args):
    seq1=args.seq1
    seq2=args.seq2
    gap=args.GAP
    fasta=args.fasta
    substitution_matrix=args.substitution_matrix

    if fasta is True :
        (header1, seq1 ) = readSeqFromFastaFile(seq1)
        (header2, seq2 ) = readSeqFromFastaFile(seq2)
    else:
        header1 = ""
        header2 = ""

    # import substitution matrix. Could be read from a file 
    if substitution_matrix == "blosum62":
        sm=blosum62
    elif substitution_matrix == "blosum45":
        sm=blosum45
    elif substitution_matrix == "blosum80":
        sm=blosum80
    else :
        print ( "Substitution matrix not defined")

    (aligned1,aligned2,indicators,identicals,similars,positives,gaps,score,seq1_chunk_spent,seq2_chunk_spent)=S_W(seq1,seq2,sm,gap)
    ali = PAIR_ALIGNED(aligned1,aligned2,indicators,header1,header2,identicals,similars,positives,gaps,gap,score,substitution_matrix,len(seq1),len(seq2),seq1_chunk_spent,seq2_chunk_spent,"Smith-Waterman")
    ali.print_alignment()

def main():
    args = parse_args()
    do_stuff(args)

if __name__ == '__main__':
    main()

