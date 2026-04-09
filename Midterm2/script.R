library(Biostrings)
library(protti)
library(Uniprot)

aln = readDNAStringSet("metazoa_alignment.gene.fasta")

# view entire alignment 

as.character(aln)

# view gaps and masked nucleotides

letterFrequency(aln, "-")

letterFrequency(aln, "N")

# view length of aligned sequences

width(aln)

# translate the Homo sapiens sample 

homo = aln$"Homo_sapiens"

# remove gaps from sequence and returns a DNAStringSet

nogap_homo =  as(gsub("-", "", homo), "DNAStringSet")

homoaa = translate(nogap_homo)

# writes the translated sequence to a fasta file in the working directory

writeXStringSet(homoaa, filepath = "Homo_sapiens_aa.fasta")

# after BLAST retrieves an accession number, write to a variable 

acc = "P54098"

# gets gene ontology data and formats into .csv and .tiff files in the working directory

info = GetProteinGOInfo(acc, getwd())

data = PlotGoInfo(info, getwd())
