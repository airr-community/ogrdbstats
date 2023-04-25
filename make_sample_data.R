# Make sample data files in the current directory, using data shipped with Tigger

library(tigger)

# Save the sample germline sequences to v_germline_gapped.fasta
writeFasta(SampleGermlineIGHV, 'v_germline_gapped.fasta')

# Save the sample repertoire to repertoire.tsv
write.table(SampleGermlineIGHV, 'rep_genotyped.tsv', sep='\t', row.names=F)
