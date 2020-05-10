#!/usr/bin/python3

################################################################################
# Program name      : readmatrix.py
# Original Author   : Narendra Kumar, narekum@gmail.com
# Date created      : 03 May 2020
# Purpose           : Smith-Waterman Algorithm for local sequence alignment
#                     Reads substitution matrix
################################################################################

# TODO: use this code to include ability to read a user defined substitution matrix
# Currently if just reformats the matrix to the python dictionary code which must be pasted in 
# substitution_matrices.py in lib directory

import sys
import re
substitution_matrix = sys.argv[1]

sm={}
residues=[]
with open(substitution_matrix) as fh:
	for line in fh:
		line=line.strip()
		if line.startswith("#"):
			continue
		if (not residues) is True:
			residues=re.split(r"\s+",line)
		else:
			scores=re.split(r"\s+",line)
			res2=scores.pop(0)
			len(scores)
			for i in range(len(scores)):
				print( "('%s','%s'): %2d" %  ( res2,residues[i],int(scores[i])), end=", ")
				sm[(res2,residues[i])]=scores[i]
			print()


