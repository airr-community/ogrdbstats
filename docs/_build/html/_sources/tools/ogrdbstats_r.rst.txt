.. _ogrdbstats_r:

ogrdbstats.R - command-line description
=======================================

The description applies equally to the Dockerised command `run_ogrdbstats`.

.. code-block:: bat

	usage: ogrdbstats.R [--] [--help] [--plot_unmutated] [--all_novel]
		   [--opts OPTS] [--inf_file INF_FILE] [--hap_gene HAP_GENE]
		   [--format FORMAT] REF_FILE SPECIES READ_FILE CHAIN

	Create genotype statistics

	positional arguments:
	  REF_FILE		germline set filename
	  SPECIES		species name used in field 3 of the IMGT
	    			germline set header, with spaces removed, e.g.
	    			Homosapiens for Human
	  READ_FILE		name of file containing annotated reads in
	    			AIRR, CHANGEO, IMPRE or IgDiscover format
	  CHAIN			one of IGHV, IGKV, IGLV, IGHD, IGHJ, IGKJ,
				IGLJ, TRAV, TRAJ, TRBV, TRBD, TRBJ, TRGV, TRGJ,
				TRDV, TRDD, TRDJ
					
	flags:
	  -h, --help            show this help message and exit
	  -p, --plot_unmutated  Plot base composition using only unmutated
	                        sequences (V-genes only)
	  -a, --all_novel       Treat all alleles in the germline set as if novel
							
	optional arguments:
	  -x, --opts            RDS file containing argument values
	  -i, --inf_file        sequences of inferred novel alleles (FASTA
	                        format)
	  --hap_gene            haplotyping gene, e.g. IGHJ6
	  -f, --format          Output report format: pdf, html or none
							[default: pdf]

Positional Arguments:

`REF_FILE` - pathname of a FASTA file containing IMGT gap-aligned germline sequences. Usually this would be downloaded from IMGT.

`SPECIES` - should contain the species name used in field 3 of the IMGT REF_FILE FASTA header, with spaces removed, e.g. Homosapiens for Human. If you are not using an IMGT REF_FILE, you can use any single word for the species here, and REF_FILE should only contain genes for that species.

`READ_FILE` - pathname of a tab-separated file containing the annotated reads used to infer the genotype, in AIRR, CHANGEO or IgDiscover format

`CHAIN` - specifies the sequence type to be analysed. 

Optional Arguments:

`INF_FILE` - pathname of a FASTA file containing sequences of inferred novel alleles. This file must be provided if the read file contains assignments to alleles that are not listed in REF_FILE.

`HAP_GENE` - the gene to be used for haplotyping analysis (see haplotyping section below)

