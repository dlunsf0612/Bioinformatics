# convert to protein data

library(Biostrings)

orchids = readDNAStringSet("orchids.fasta")

# which samples have ambiguous nucleotides

letterFrequency(orchids, "N")

# samples 60 and 88 have N

# sample 60 has other ambiguous characters
# Letters S, M, and K were found in the sample and must be removed

letterFrequency(orchids, "RYDWSKMBHV")

letterFrequency(orchids[60], "S")

letterFrequency(orchids[60], "M")

letterFrequency(orchids[60], "K")

clean = function(seq, x){

	new = DNAStringSet(c())

	positions = c()
	
	for(i in 1:length(seq[x][[1]])){
	
		if (as(seq[x][[1]][i], "character") == "N"){

			positions = c(positions, i)
		 		
		 }
		
	 }

	new = c(new, replaceAt(seq[x], IRanges(positions, positions), "")) 

 return(new)

}

clean88 = clean(orchids, 88)



# written for sample 60

newclean = function(seq, x){

	new = DNAStringSet(c())

	positions = c()
	
	for(i in 1:length(seq[x][[1]])){
	
		if (as(seq[x][[1]][i], "character") == "K" || as(seq[x][[1]][i], "character") == "M" || 
			as(seq[x][[1]][i], "character") == "S" || as(seq[x][[1]][i], "character") == "N"){

			positions = c(positions, i)
		 		
		 }
		
	 }

	new = c(new, replaceAt(seq[x], IRanges(positions, positions), "")) 

 return(new)

}

clean60 = newclean(orchids, 60)

# add all samples without ambiguous characters to a DNAStringSet

cleanDNA = DNAStringSet(c())

for (i in 1:length(orchids)){

	if (i != 60 && i != 88){

	cleanDNA = c(cleanDNA, orchids[i])

	}
	
	if (i == 60){

	cleanDNA = c(cleanDNA, clean60)

	}

	if (i == 88){

	cleanDNA = c(cleanDNA, clean88)

	}

}

# translate to protein data

orchidsaa = translate(cleanDNA)

# write new fasta file

writeXStringSet(orchidsaa, filepath = "orchid_aa.fasta")

