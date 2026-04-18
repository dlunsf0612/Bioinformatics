# MAFFT script

mafft --retree 1 --memsave orchids.fasta > msa_100.fasta &&  echo -e "\a" \
 && osascript -e 'display notification "Task Finished" with title "Terminal"'