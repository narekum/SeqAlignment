#!/usr/bin/python3

################################################################################
# Program name      : Library files for sequence alignment
# Original Author   : Narendra Kumar, narekum@gmail.com
# Date created      : 03 May 2020
################################################################################

def N_W (seq1,seq2,sm,gap):
	"""Needleman and Wunsch Algorith for Global Seqeunce Alignment"""

	# declare dictionaries
	scores = {(0,0):0}
	pointers = {(0,0):"n"}

	# Initialize matrix. Fill first row and first column i.e. gap penalties
	for j in range(1,len(seq1)+1):
		scores[(0,j)] = j * gap
		pointers[(0,j)] = "l"

	for i in range(1,len(seq2)+1):
		scores[i,0] = i * gap
		pointers[(i,0)] = "u"

	# Fill the matrix
	for i in range(1, len(seq2)+1):
		for j in range(1,len(seq1)+1):
			character1=seq1[j-1]
			character2=seq2[i-1]
			character1_diagonal = seq1[j-2]
			character2_diagonal = seq2[i-2]
		
			#Get diagonal and gap scores
			diagonal = scores[(i-1,j-1)] + sm[(character1_diagonal,character2_diagonal)]    
			up=scores[(i-1,j)] + gap
			left = scores[(i,j-1)] + gap
	    
			# select best from diagonal, up and left
			if diagonal >= up :
				if diagonal >= left :
					scores[(i,j)] = diagonal
					pointers[(i,j)] = "d"
				else :
					scores[(i,j)] = left
					pointers[(i,j)] = "l"
			else :
				if up >= left :
					scores[(i,j)] = up
					pointers[(i,j)] = "u"
				else :
					scores[(i,j)] = left
					pointers[(i,j)] = "l"
		    
	# Begin traceback
	align_seq1=[]
	align_seq2=[]
	align_indicators=[]

	j=len(seq1)
	i=len(seq2)

	score=scores[i,j]

	while pointers[(i,j)] != "n" :
		#print (i,j, pointers[(i,j)])
		if pointers[(i,j)] == "d":
			align_seq1.append(seq1[j-1])
			align_seq2.append(seq2[i-1])

			if seq1[j-1] == seq2[i-1] :
				align_indicators.append("|")
			else:
				align_indicators.append(".")

			i=i-1
			j=j-1
		
		elif pointers[(i,j)] == "l":
			align_seq1.append(seq1[j-1])
			align_seq2.append("-")
			align_indicators.append(" ")
			j=j-1
		elif pointers[(i,j)] == "u":
			align_seq1.append("-")
			align_seq2.append(seq2[i-1])
			align_indicators.append(" ")
			i=i-1
		
	# Get aligned sequences
	aligned1= "".join(list(reversed(align_seq1)))
	aligned2= "".join(list(reversed(align_seq2)))
	indicators= "".join(list(reversed(align_indicators)))

	identicals = indicators.count("|") * 100 / len(indicators) 
	similars = indicators.count(".") * 100 / len(indicators)

	# Calculate statistics
	identicals = indicators.count("|") * 100 / len(indicators) 
	similars = indicators.count(".") * 100 / len(indicators)
	positives = identicals + similars
	gaps = indicators.count(" ") * 100 / len(indicators)

	seq1_chunk_spent=0
	seq2_chunk_spent=0

	return (aligned1,aligned2,indicators,identicals,similars,positives,gaps,score,seq1_chunk_spent,seq2_chunk_spent)

def S_W (seq1, seq2, sm, gap):
	# declare dictionaries
	scores = {(0,0):0}
	pointers = {(0,0):"n"}

	# Initialize matrix. Fill first row and first column is all set to zero
	for j in range(1,len(seq1)+1):
		scores[(0,j)] = 0
		pointers[(0,j)] = "n"

	for i in range(1,len(seq2)+1):
		scores[i,0] = 0
		pointers[(i,0)] = "n"

	# Fill the matrix and record highest score and its i and j value
	highest_score = 0
	highest_score_i = 0
	highest_score_j = 0

	for i in range(1, len(seq2)+1):
		for j in range(1,len(seq1)+1):
			character1=seq1[j-1]
			character2=seq2[i-1]
			character1_diagonal = seq1[j-2]
			character2_diagonal = seq2[i-2]
		
			#Get diagonal and gap scores
			diagonal = scores[(i-1,j-1)] + sm[(character1_diagonal,character2_diagonal)]    
			up=scores[(i-1,j)] + gap
			left = scores[(i,j-1)] + gap

			# Set scores to zero if its less than zero
			if ( diagonal <= 0 ) and ( up <= 0) and (left <= 0):
				scores[(i,j)] = 0
				pointers[(i,j)] = 'n'
				continue 

			# select best from diagonal, up and left
			if diagonal >= up :
				if diagonal >= left :
					scores[(i,j)] = diagonal
					pointers[(i,j)] = "d"
				else :
					scores[(i,j)] = left
					pointers[(i,j)] = "l"
			else :
				if up >= left :
					scores[(i,j)] = up
					pointers[(i,j)] = "u"
				else :
					scores[(i,j)] = left
					pointers[(i,j)] = "l"

			# Reset highest score if its more than previous highest
			if scores[(i,j)] > highest_score :
				highest_score = scores[(i,j)] 
				highest_score_i = i
				highest_score_j = j
		    
	# Begin traceback
	align_seq1=[]
	align_seq2=[]
	align_indicators=[]

	# Start from the highest score's i and j 
	j=highest_score_j
	i=highest_score_i
	score=scores[i,j]

	while pointers[(i,j)] != "n" :
		#print (i,j, pointers[(i,j)])
		if pointers[(i,j)] == "d":
			align_seq1.append(seq1[j-1])
			align_seq2.append(seq2[i-1])

			if seq1[j-1] == seq2[i-1] :
				align_indicators.append("|")
			else:
				align_indicators.append(".")

			i=i-1
			j=j-1
		
		elif pointers[(i,j)] == "l":
			align_seq1.append(seq1[j-1])
			align_seq2.append("-")
			align_indicators.append(" ")
			j=j-1

		elif pointers[(i,j)] == "u":
			align_seq1.append("-")
			align_seq2.append(seq2[i-1])
			align_indicators.append(" ")
			i=i-1
		
	i_start = i
	j_start = j

	# Get aligned sequences
	aligned1="".join(list(reversed(align_seq1)))
	aligned2="".join(list(reversed(align_seq2)))
	indicators="".join(list(reversed(align_indicators)))

	identicals = indicators.count("|") * 100 / len(indicators) 
	similars = indicators.count(".") * 100 / len(indicators)

	# Calculate statistics
	identicals = indicators.count("|") * 100 / len(indicators) 
	similars = indicators.count(".") * 100 / len(indicators)
	positives = identicals + similars
	gaps = indicators.count(" ") * 100 / len(indicators)

	seq1_chunk_spent = j_start 
	seq2_chunk_spent = i_start

	return (aligned1,aligned2,indicators,identicals,similars,positives,gaps,score,seq1_chunk_spent,seq2_chunk_spent)

class PAIR_ALIGNED:
	"""A class containing aligned pairwise sequences. 
	and method to print them"""
	def __init__ (self,aligned1,aligned2,indicators,header1,header2,identicals,similars,positives,gaps,gap_penalty,score,substitution_matrix,seq1_len,seq2_len,seq1_chunk_spent,seq2_chunk_spent,algorithm):
		self.aligned1=aligned1  # first sequence aligned
		self.aligned2=aligned2  # second sequence aligned
		self.indicators=indicators # indicators of alignment
		self.header1=header1	
		self.header2=header2
		self.identicals=identicals	
		self.similars=similars
		self.positives=positives
		self.gaps=gaps
		self.gap_penalty=gap_penalty
		self.score=score
		self.substitution_matrix=substitution_matrix
		self.seq1_len=seq1_len	# length of original sequence1
		self.seq2_len=seq2_len  # length of original sequence2
		self.seq1_chunk_spent=seq1_chunk_spent  # sequence1 not included in alignemnt
		self.seq2_chunk_spent=seq2_chunk_spent  # sequence2 not included in alignemnt
		self.algorithm=algorithm

	seq_chunks=60 # Default length in each line in output alignment

	def set_printsize(self,n):
		self.seq_chunks=n

	def print_alignment(self):
		
		print ()
		print ( self.algorithm , " algorithm for sequence alignment")
		print ("written by: narekum@gmail.com")
		print ()
		print ("Substn matrix: %s" % (self.substitution_matrix))
		print ("Gap penalty  : %s" % (self.gap_penalty))
		print ()
		print ("Sequence1    : %s length %d" % (self.header1, self.seq1_len))
		print ("Sequence2    : %s length %d" % (self.header2, self.seq2_len))

		print ()
		print ("Identicals: %3.2f%%    Positives: %3.2f%%    Gaps: %3.2f%%    Score:%d" % (self.identicals, self.positives, self.gaps, self.score))
		print ()

		for i in range(0,len(self.indicators),self.seq_chunks):

		    seq1_chunk = len(self.aligned1[i:i+self.seq_chunks]) - self.aligned1[i:i+self.seq_chunks].count("-")
		    seq2_chunk = len(self.aligned2[i:i+self.seq_chunks]) - self.aligned2[i:i+self.seq_chunks].count("-")

		    seq1_chunk_start = self.seq1_chunk_spent
		    seq2_chunk_start = self.seq2_chunk_spent

		    seq1_chunk_end = self.seq1_chunk_spent + seq1_chunk
		    seq2_chunk_end = self.seq2_chunk_spent + seq2_chunk


		    print ( "seq1: %4d  %s  %-4d" % ( seq1_chunk_start + 1, self.aligned1[i:i+self.seq_chunks], seq1_chunk_end ) )
		    print ( "            %s" % self.indicators[i:i+self.seq_chunks])
		    print ( "seq2: %4d  %s  %-4d" % (seq2_chunk_start + 1,self.aligned2[i:i+self.seq_chunks], seq2_chunk_end))
		    print ()

		    self.seq1_chunk_spent = seq1_chunk_start + seq1_chunk
		    self.seq2_chunk_spent = seq2_chunk_start + seq2_chunk

