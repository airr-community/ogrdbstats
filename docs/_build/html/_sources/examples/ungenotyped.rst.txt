.. _ungenotyped:

Using OGRDBstats with any repertoire
====================================

This section demonstrates how to conduct a V-gene analysis on a repertoire that has not been processed by a genotyping tool, and therefore for which no novel alleles have been infered. 

We will use the sample repertoire provided with the TIgGER package for the example. An IMGT-gapped V germline set is required. 
This must match the germline set used to annotate the repertoire. Here we use the example germline set included with TIgGER. This R
code will save the sample files to the current directory:

.. code-block:: R

	library(tigger)

	# Save the sample germline sequences to v_germline_gapped.fasta
	writeFasta(SampleGermlineIGHV, 'v_germline_gapped.fasta')
	# Save the sample repertoire to repertoire.tsv
	write.table(SampleGermlineIGHV, 'rep_genotyped.tsv', sep='\t', row.names=F)
	
If you are using the OGRDBstats Docker image, you can create these files in the current directory with the command

.. code-block:: bash

    $ docker run -v $(pwd):/scratch ogrdbstats:stable make_sample_data
	
With these commands, the following sample files should be in the current directory and can be used to demonstrate the use of ogrdbstats on a repertoire:

- `v_germline_gapped.fasta` - The germline set used for annotation
- `rep_genotyped.tsv` - The annotated repertoire, in AIRR format.

The following command will run OGRDBstats on these files:

.. code-block:: bash

    $ Rscript ogrdbstats.R v_germline_gapped.fasta Human rep_genotyped.tsv IGHV
	
The species is defined as Human.

OGRDBstats will provide status as it runs. When processing is complete, there will be two output files: a CSV containing statistics and
a pdf containing plots. The filenames follow the filename of the repertoire: in this case they will be called `rep_genotyped_ogrdb_report.csv`
and `rep_genotyped_ogrdb_plots.pdf`. 

You can produce a report for D or J- genes by changing the final argument. V, (D) and J gene analysis in all 7 loci is supported.

Please note that the germline set must contain all germline sequences used for annotation. For V-gene analysis, the 
germline set must be gapped. The repertoire must be in AIRR or CHANGE-O format. 


A full description of the parameters and options available is given in :ref:`ogrdbstats_r`.


