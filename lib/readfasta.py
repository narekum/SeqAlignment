#!/usr/bin/python3

################################################################################
# Program name      : Library files for sequence alignment
# Original Author   : Narendra Kumar, narekum@gmail.com
# Date created      : 03 May 2020
################################################################################

def readSeqFromFastaFile(fastafile):
        """Reads fasta file and returns a tuple of header and sequence"""
        seq=''
        with open(fastafile, 'r') as readfasta:
                for line in readfasta.readlines():
                        line = line.strip()
                        if line[0] == '>':
                                header=line
                        else:
                                seq += line
        return (header,seq)


