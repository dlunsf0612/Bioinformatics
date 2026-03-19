
library(Biostrings)
library(protti)
library(UniprotR)
library(seqinr)

# read in DNA data then translate to amino acid sequence
# write the translated sequence into a fasta file
# confirm the fasta file is in the working drectory 

Brapa = readDNAStringSet("Brapmt.fasta")

aaBrapa = Biostrings::translate(Brapa)

write.fasta(aaBrapa, names = "Rapa", file.out = "Brassica_rapa_aa.fasta")

list.files()

# after BLASTing the created fasta file
# create a text file with the accession numbers
# read the text file into R
# change the table to a character vector 

blast = read.csv("Accession.txt", header = FALSE)

blast = as(blast, "vector")

# Get the Gene Ontology (GO) terms and plot them into a data frame
# using functions from UniprotR
# My accession numbers don't work so the provided accession numbers P0A799 and P08839 were used
# The GETSeqFastaUniprot function was used to write the accession to a fasta file

blast = c("P0A799", "P08839")

GETSeqFastaUniprot(blast, getwd(), "protein.fasta")

info = GetProteinGOInfo(blast, getwd())

data = PlotGoInfo(info, getwd())



# Find information on the pathology and related diseases 

pathology = GetPathology_Biotech(blast)

Get.diseases(pathology, getwd())

 # No allergenic properties were found
 # Mutagens 122, 189, and 502 were found in P08839
 # No other data was found using the GetPathology_Biotech
 # Neither accession is associated with any diseases

 # fetch_uniprot function uses the UniProt accession numbers to retrieve more data on the proteins
 # One data column is the xref_pdb (Protein Database ID)
 # pdbs can be used in the fetch_pdb function to pull available structural information about the protein
 # both these function return dataframes 


fetch = fetch_uniprot(blast)

fetch$xref_pdb

pdb = c("1EZA", "1EZB", "1EZC", "1EZD", "1ZYM", "2EZA", "2EZB", "2EZC", "2HWG", "2KX9", "2L5H", "2MP0", "2N5T", "2XDF", "3EZA", "3EZB", "3EZE", "6V9K", "6VU0", "1ZMR")


table = fetch_pdb(pdb, batchsize = 100)

# the fetch_alphafold_prediction uses UniProt accession numbers to get alphafold prediction data 
# this can be used to develop a 3D model of the protein and returns a dataframe 


alphafold = fetch_alphafold_prediction(blast, return_data_frame = TRUE)


