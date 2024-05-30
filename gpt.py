from Bio import SeqIO
from pathlib import Path
import numpy as np
import time 



def calculate_C(X):
    """ Calculate the C array using numpy for better performance. """
    unique, counts = np.unique(np.array(list(X)), return_counts=True)
    total = 0
    C = {}
    i = 0
    for char in sorted(unique):
        C[char] = total
        total += counts[i]
        i += 1
    return C

def calculate_O(B):
    """ Calculate the O array using numpy. """
    B = np.array(list(B))
    unique = np.unique(B)
    O = {char: np.zeros(len(B), dtype=int) for char in unique}
    for char in unique:
        O[char] = np.cumsum(B == char)  #- (B == char)
        if O[char][0] == 1:
            O[char][0] = 0
    return O

def backward_search(W,B,C,O,n):
    """ Perform backward search using numpy arrays. """
    m = len(W)
    i = m - 1
    char = W[i]
    k = 1
    l = len(B) -1
    while k <= l and i >= (n):
        char = W[i]
        k = C[char] + O[char][k - 1] + 1
        l = C[char] + O[char][l]
        i -= 1
    return k, l

def inexact_search(W, z, B, C, O):
    """ Implement iterative inexact search to avoid recursion limits. """
    stack = [(len(W) - 1, z, 0, len(B) - 1)]
    results = []
    while stack:
        i, z, l, r = stack.pop()
        if z < 0:
            continue
        if i < 50:
            results.append((l, r))
            continue
        for b in set(W):
            l_new = C[b] + (O[b][l - 1] if l > 0 else 0)
            r_new = C[b] + O[b][r] - 1
            if l_new <= r_new:
                stack.append((i - 1, z - (W[i] != b), l_new, r_new))
    if len(results) != 0:
        lr = results[0]
        return lr
    if len(results) == 0:
        return None
        


def bwt_transform(X):
    """ Simplified BWT transform using sorted rotations and numpy. """
    X = X + '$'
    rotations = sorted(X[i:] + X[:i] for i in range(len(X)))
    r = ''.join(rot[-1] for rot in rotations)
    rotations.clear()
    return r,X



FASTA_FILE = Path("/Users/sasakitakato/Desktop/itohlab/研修/mapping/mapping_program/ref_SE11.fasta")
FASTQ_FILE = Path("/Users/sasakitakato/Desktop/itohlab/研修/mapping/mapping_program/Illumina_SE11.fastq")
BWT_FILE = Path("/Users/sasakitakato/Desktop/itohlab/研修/BWT.txt")
NN_FILE = Path("/Users/sasakitakato/Desktop/itohlab/研修/nn.txt")

def read_fasta_file(fastafile):
    fastaseq = []
    for record in SeqIO.parse(fastafile, "fasta"):
        fastaseq.append([record.id,record.description,record.seq])
    return fastaseq
def read_fastq_file(fastqfile):
    sequences = []
    # SeqIO.parseはFASTAやFASTQファイルのパースに使用される
    for records in SeqIO.parse(fastqfile,"fastq"):
        sequences.append([records.id,records.seq,records.letter_annotations['phred_quality']])
    return sequences

fasta = read_fasta_file(FASTA_FILE)
fastq = read_fastq_file(FASTQ_FILE)


fasta_C = []
fasta_O = []
nn = []



# Python3 program for building suffix 
# array of a given text

# Class to store information of a suffix
class suffix:
	def __init__(self):
		self.index = 0
		self.rank = [0, 0]

# This is the main function that takes a string 'txt' of size n as an argument, builds and return the suffix array for the given string
def buildSuffixArray(txt, n):
    # A structure to store suffixes and their indexes
	suffixes = [suffix() for _ in range(n)]

	# Store suffixes and their indexes in an array of structures. The structure is needed to sort the suffixes alphabetically and maintain their old indexes while sorting
	for i in range(n):
		suffixes[i].index = i
		suffixes[i].rank[0] = (ord(txt[i]) -ord("a"))
		suffixes[i].rank[1] = (ord(txt[i + 1]) -ord("a")) if ((i + 1) < n) else -1

	# Sort the suffixes according to the rank and next rank
	suffixes = sorted(suffixes, key = lambda x: (x.rank[0], x.rank[1]))

	# At this point, all suffixes are sorted according to first 2 characters. Let us sort suffixes according to first 4 characters, then first 8 and so on
	ind = [0] * n # This array is needed to get the index in suffixes[] from original index.This mapping is needed to get next suffix.
	k = 4
	while (k < 2 * n):
		# Assigning rank and index values to first suffix
		rank = 0
		prev_rank = suffixes[0].rank[0]
		suffixes[0].rank[0] = rank
		ind[suffixes[0].index] = 0

		# Assigning rank to suffixes
		for i in range(1, n):
			# If first rank and next ranks are same as that of previous suffix inarray, assign the same new rank to this suffix
			if (suffixes[i].rank[0] == prev_rank and suffixes[i].rank[1] == suffixes[i - 1].rank[1]):
				prev_rank = suffixes[i].rank[0]
				suffixes[i].rank[0] = rank
			# Otherwise increment rank and assign 
			else: 
				prev_rank = suffixes[i].rank[0]
				rank += 1
				suffixes[i].rank[0] = rank
			ind[suffixes[i].index] = i

		# Assign next rank to every suffix
		for i in range(n):
			nextindex = suffixes[i].index + k // 2
			suffixes[i].rank[1] = suffixes[ind[nextindex]].rank[0] \
				if (nextindex < n) else - 1

		# Sort the suffixes according to first k characters
		suffixes = sorted(suffixes, key = lambda x: (x.rank[0], x.rank[1]))
		k *= 2

	# Store indexes of all sorted suffixes in the suffix array
	# suffixArr = [0] * n
	suffixArr = ""
	number = ""
	for i in range(n):
		# suffixArr[i] = txt[suffixes[i].index-1]
		suffixArr += txt[suffixes[i].index-1]
		number += str(suffixes[i].index)
		# suffixArr[i] = suffixes[i].rank[0]
	# Return the suffix array
	return suffixArr,number

# A utility function to print an array of given size
def printArr(arr, n):
	for i in range(n):
		print(arr[i], end = " ")
	print()

start = time.time()



bwtfile = open(BWT_FILE,"r")
BWT  = (bwtfile.read().splitlines())
bwtfile.close()

# file = open(NN_FILE,"r")
# a = (file.readline().split("	"))
# b = (file.readline().split("	"))
# c = (file.readline().split("	"))
# d = (file.readline().split("	"))
# e = (file.readline().split("	"))
# f = (file.readline().split("	"))
# g = (file.readline().split("	"))
# nn = [a,b,c,d,e,f,g]
# file.close




for i in range(len(fasta)):
# Driver code
	if __name__ == "__main__":
		txt = fasta[i][2] + "$"
		n = len(txt)
		B,number = buildSuffixArray(txt, n)
		BWT.append(B)
		nn.append(number)
		fasta_C.append(calculate_C(fasta[i][2]))
		fasta_O.append(calculate_O(BWT[i]))
en = time.time()
print(en - start)

# with open("BWT.txt","w") as file0:
#     for i in BWT:
#         file0.write(i)
#         file0.write("\n")
# file0.close

# with open("nn.txt","w") as file0:
#     for i in nn
#         file0.write(i)
#         file0.write("\n")
# file0.close

print(fasta_C)
print(fasta_O)

with open("untitled.txt", "w") as file1:
    for seq in fastq:
        file1.write(seq[0])
        rseq = seq[1].reverse_complement()
        rrseq = rseq[:70]
        for i in range(len(BWT)):
            k,l = backward_search(seq[1],BWT[i],fasta_C[i],fasta_O[i],75) 
            kr,lr = backward_search(rseq,BWT[i],fasta_C[i],fasta_O[i],75) 
            kk,ll = backward_search(seq[1][:70],BWT[i],fasta_C[i],fasta_O[i],50) 
            kkr,llr = backward_search(rrseq,BWT[i],fasta_C[i],fasta_O[i],50) 
            if k <= l:
                file1.write("\t")
                a = int(nn[i][k])
                file1.write(fasta[i][0])
                file1.write("\t")
                file1.write(str(a-74))
                file1.write("\t")
                file1.write("+")
                break
                
            elif kr <= lr:
                file1.write("\t")
                a = int(nn[i][kr])
                file1.write(fasta[i][0])
                file1.write("\t")
                file1.write(str(a-74))
                file1.write("\t")
                file1.write("-")
                break

            elif kk <= ll:
                file1.write("\t")
                a = int(nn[i][kk])
                file1.write(fasta[i][0])
                file1.write("\t")
                file1.write(str(a-49))
                file1.write("\t")
                file1.write("+")
                break
            
            elif kkr <= llr:
                file1.write("\t")
                a = int(nn[i][kkr])
                file1.write(fasta[i][0])
                file1.write("\t")
                file1.write(str(a-49))
                file1.write("\t")
                file1.write("-")
                break
            else:
                continue
        file1.write("\n")
file1.close

end = time.time()
time_diff = end - start
print(time_diff) 

