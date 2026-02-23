# Midterm

library(Biostrings)
library(msa)
library(pwalign)
library(annotate)
library(seqinr)

# 1. import and align sequences 

seqs = readDNAStringSet("midterm.fasta")
seqs

names(seqs) = c("Homo1", "Homo2", "Homo3", "Homo4", "Homo5", "Homo6", "Homo7", "Homo8", "Homo9", "Homo10", "Homo11", "Homo12", "Homo13", "Homo14", "Homo15", "Homo16", "Homo17", "Homo18", "Homo19", "Homo20")
seqs

aln = msaMuscle(seqs)
aln

# 2. conduct distance alignment to determine how well the alignment is 
# and see how different each sample is from one another

daln = msaConvert(aln,"seqinr::alignment") 
d = dist.alignment(daln, matrix = c("identity","similarity")) 
as.matrix(d)


# The alignment is good because most of the values in the matrix are 0
# so these sequences are identical 
# Homo6 (Sample 6) is the most different from the other samples
# Homo4 and Homo10 also have higher distance alignment scores 

# 3. Calculate the consensus sequence

conseq = msaConsensusSequence(aln)
conseq

# 2. Counting the number of gaps in the consensus sequence 
# also shows how good the alignment is 

consensus = as(conseq, "DNAStringSet")
letterFrequency(consensus, "-")

# no gaps are present in the consensus sequence
# the alignment is good 


# 4. Calculate the GC content for all sequences in the alignment

# creates a data.frame with the GC content of all sequences in the alignment 

alnBio = as(aln, "DNAStringSet")

consensus = as(conseq, "DNAStringSet")

for(i in 1:length(alnBio)){

	if(i == 1){

		GC = letterFrequency(alnBio[1], "GC")/width(alnBio[1]) * 100
	}

	if(i > 1){

		GC = c(letterFrequency(alnBio[i], "GC")/width(alnBio[i]) * 100)
	}

}

GCcontent = data.frame(Name = names(alnBio), GCperecent = GC)
GCcontent


# 5. Determine how different the samples are from one another 
# Determine the type of mutations they have 

# the function msapid converts an msa object
# into a Biostrings DNAStringSet
# then it does a pairwise alignment of each sequence
# comparing them to the argument "sequence" 
# the sequence arguemnt must be a DNAString or similar object  
# in this case it is returned from 
# the msaConsensusSequence function
# it returns at data.frame that displays 
# the names of each sequence compared
# and their given percent identity score  

msapid = function(alignment, sequence){

consensus = as(sequence, "DNAStringSet")	

alnBio = as(alignment, "DNAStringSet")

for (i in 1:length(alnBio)){

	if (i == 1){

		paln = (pairwiseAlignment(consensus, alnBio[1]))

	}

	if (i > 1){

		paln = c(paln, pairwiseAlignment(consensus, alnBio[i]))
	
	}

} 

piddata = data.frame(Name = names(alnBio), pid = pid(paln))

return(piddata)

}

msapid(aln, conseq)

# Homo6, Homo4, and Homo10 all have mutations


# findsnp is a function that finds the positions where nucleotides 
# are dissimilar between aligned sequences
# uses a nested loop so it takes awhile to run 
# i is the nucleotide position
# n is the sequence in the DNAStringSet
# the argument alignment is an msa or Biostrings object
# an alignment must be done before using the function 


# converts the argument to a Biostrings DNAStringSet
# creates an empty vector
# for all sequences in the DNAStringSet and for all nucleotides in a given sequence 
# if a nucleotide at the position i in the sequence n is different from
# the nucleotide at position in in the sequence n+1  
# then i is put into a vector
# findsnp returns a vector of integers


findsnp = function(alignment){
	
	aligninBio = as(alignment, "DNAStringSet") 

	positions = c()  

	for(n in 1:length(aligninBio)){
		
		for(i in 1:length(aligninBio[1][[1]])){

		 	
		 	if (i <= length(aligninBio[1][[1]]) && n != 20){

		 		
		 		if (aligninBio[n][[1]][i] != aligninBio[n + 1][[1]][i]){
		 		
		 			positions = c(positions, i)  

		 		}
		 	} 
	 	}
	}

	return(positions)
}

# check the values at the positions found in the vector positions from the previous code
# deterime if there are mismatching nucleotides 

# shownts displays the nucleotides of a Biostrings DNAStringSet 
# in this project it functions to show the SNPs at a given position 
# this position is determined by the findsnp function 

# The argument alignment is an aligned msa or Biostrings object 
# the argument x is the position that will be examined 
# The code traverses the alignment and creates a character vector
# that includes all nucleotides at a given position for each sequence in the alignment 


shownts = function(alignment, x){

	aligninBio = as(alignment, "DNAStringSet")
	
	for (i in 1:length(aligninBio)){
		
		if(i == 1){
			
			difference = c()
			difference = as(aligninBio[[i]][x], "character") 
		}

		if(i > 1){

			difference = c(difference, as(aligninBio[[i]][x], "character"))
		}

		if(i == length(aligninBio)){
			
			diffDataFrame = data.frame(Name = names(aligninBio), SNP = difference)
			print(diffDataFrame)
		}
	}
}

# save the returned values the variable snps

snps = findsnp(alnBio)

# displays a data.frame for each position in snps
# the dat frame includes the name of the sequence and the nucleotide 

for(i in 1:length(snps)){

	shownts(alnBio, snps[i])
}


# Individual 4 has one silent mutation C --> A at position 39 of the alignment
# Individual 10 has two silent mutations C --> G at position 39 A -> T at position 45 
# the silent mutations were determined by translating the aligned sequences
# Individual 1, 4, and 10 have the same amino acid sequence
# 1 doesnt have a mutation so the mutations in 4 and 10 must be silent 

# Individual 6 has a deletion at position 3 causing a frame-shift mutation
# Moreover when aligned sequence of individual 6 has a number of substitution mutations 

# A --> G at base 47
# T --> C at base 134
# A --> T at base 145
# C --> G at base 152
# G --> C at base 586
# A --> G at base 623



# 6. BLAST a sequence to compare it to a database 
# determine the identity of the gene 
# blastSequences uses the annotate package 

blast = blastSequences(seqs[1], hitListSize = 3, as = "data.frame")
blast

# the sequence comes from the beta subunit of the HBB gene
# the accession number of the best match is LC121775
# this sequence has 100% pid with query 
# LC121775 has a silent mutation (HbLimassol Cd8(AAG>AAC))

# 7. translate the most mutant sequence and write it to a fasta file 

aaseq = Biostrings::translate(seqs[6])
aaseq

# write into fasta file

write.fasta(aaseq, names = "Homo_sapien_sample_6", file.out = "mutant_aa.fasta")
list.files() # ensure the file was properly written and put in the directory 

# 8. The protein matches to a beta subunit of hemoglobin in humans 
# the accession number is KAI2558340


# 9. The HBB gene is associated with sickle cell anemia and beta-thalassemia
# individual 6 has both sickle cell and beta-thalassemia  
# This is called sickle-beta thalassemia. AY356351 is the acession number 
# Both individual 4 and 6 express the wild type because thei mutations are silent 



