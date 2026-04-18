# convert to protein data

library(Biostrings)

orchids = readDNAStringSet("orchids.fasta")

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

set = DNAStringSet(c())

for (i in 1:length(orchids)){

	cleaned = clean(orchids, i)
	cleanDNA = DNAStringSet(c(set, cleaned))
}

cleanDNA