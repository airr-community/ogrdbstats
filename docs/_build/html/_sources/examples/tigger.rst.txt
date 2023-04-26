.. _tigger:

Using OGRDBstats with TIgGER
============================

To conduct a V-gene analysis with TIgGER:

- Use findNovelAlleles to identify novel alleles in a Change-O-formatted data set. Write these to a FASTA file.
- Use inferGenotype or inferGenotypeBayesian to infer the genotype.
- Use reassignAlleles to correct allele calls in the data set, based on the inferred genotype

The following R code, based on the `TIgGER vignette <https://tigger.readthedocs.io/en/stable/vignettes/Tigger-Vignette/>`_, will perform these steps and save the output. It uses the sample repertoire provided with Tigger:
   
.. code-block:: R

	# Load packages required for this example
	library(tigger)
	library(dplyr)

	# Save the sample germline sequences to v_germline_gapped.fasta
	writeFasta(SampleGermlineIGHV, 'v_germline_gapped.fasta')

	# Detect novel alleles in the sample repertoire
	novel <- findNovelAlleles(AIRRDb, SampleGermlineIGHV, nproc=1)

	# Extract and rows that contain successful novel allele calls
	novel_rows <- selectNovel(novel)

	# Infer the individual's genotype, using only unmutated sequences and checking
	# for the use of the novel alleles inferred in the earlier step.
	geno <- inferGenotype(AIRRDb, germline_db=SampleGermlineIGHV, novel=novel, find_unmutated=TRUE)
	
	# Save the genotype sequences to a vector
	genotype_db <- genotypeFasta(geno, SampleGermlineIGHV, novel)	
						  
	# Use the personlized genotype to determine corrected allele assignments
	# Updated genotype will be placed in the v_call_genotyped column
	sample_db <- reassignAlleles(AIRRDb, genotype_db)

	# Save the repertoire with corrected allele calls to rep_genotype.tsv
	write.table(sample_db, 'rep_genotyped.tsv', sep='\t', row.names=F)

	# Save the sequences used in the corrected repertoire, including novel allele sequences
	writeFasta(genotype_db, 'v_genotyped_seqs.fasta')					  

This code creates the following files:

- `v_germline_gapped.fasta` - the germline sequences used in the annotation
- `rep_genotyped.tsv` - the annotated reads, with a column V_GERMLINE_GAPPED containing the corrected (genotyped) V-call
- `v_genotyped_seqs.fasta` - the set of germline sequences referenced in the annotated reads, including novel sequences

These files can be provided to OGRDBstats in the following command:

.. code-block:: bash

   $ Rscript ogrdbstats.R --inf_file v_genotyped_seqs.fasta v_germline_gapped.fasta human rep_genotyped.tsv IGHV

OGRDBstats will produce two files: `rep_genotyped_ogrdb_plots.csv` and `rep_genotyped_ogrdb_report.csv`.
