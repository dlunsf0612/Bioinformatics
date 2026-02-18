Lab6 MSA in R


# in command line combine files into one file before begining 
# cat file1.fasta file2.fasta > combined_files.fasta
#Brassmtdna.fasta is the combined file

#load all libraries needed

library(Biostrings)
library(seqinr)
library(msa)

Brmtdna = readDNAStringSet("Brassmtdna.fasta")
names(Brmtdna) = c("juncea", "napus", "nigra", "oleracea", "rapa") # change the names of Brmtdna using Biostrings function names()
Brmtdna

clustalaln = msa(Brmtdna) # msa using Clustal algorithm 
clustalaln

msclaln = msaMuscle(Brmtdna) # msa using MUSCLE algorithm 
msclaln

print(msclaln, show = "complete") # shows the entire alignment using the MUSCLE algorithm 

# converts msclaln into a Biostrings compatible object using as() function 

align = as(msclaln, "DNAStringSet") 

gaps = letterFrequency(align, "-")  # shows total number of gaps 
gaps 

# find the GC content as a percentage 
# width() is a Biostrings function that returns the width of a sequences read into Biostrings

GC = letterFrequency(align, "GC")
GCpercent = GC/width(align) * 100 
GCpercent

# seqinr function must be converted to an alignment using msaConvert()

aln = msaConvert(msclaln,"seqinr::alignment") # necessary to use any functions in seqinr 

# distance alignment 

d = dist.alignment(aln, matrix = c("identity","similarity")) 
as.matrix(d)

# translate one of the sequences

Br1rna = Biostrings::translate(Brmtdna[1]) # translate is a masked function in both seqinr and Biostring this code uses the Biostrings function
Br1rna


Alignment_phyDat = msaConvert(Alignment, type="phangorn::phyDat")  # convverts msa to phangorn package

write.phyDat(Alignment_phyDat, "alignment.fasta", format = "fasta") # writes the msa to a fasta file within the working directory

list.files() # check that the new file is within your directory  


