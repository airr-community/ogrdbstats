.. _partis:

Using OGRDBstats with partis
============================

Although partis can provide annotated repertoires in AIRR format, some information required by OGRDBstats is only contained
in the default yaml format. A python script, `convert_partis.py <https://github.com/airr-community/ogrdbstats/blob/master/convert_partis.py>`_, is provided for download.
This will combine output from partis’s yaml and presto annotations, producing CHANGE-O format annotations and a FASTA file of genotype V-sequences. 

Usage of convert_partis.py:

.. code-block::

	python convert_partis.py [-h] partis_yaml partis_tsv ogrdb_recs ogrdb_vs

	positional arguments:
	  partis_yaml  .yaml file created by partis
	  partis_tsv   .tsv file created by partis presto-output mode
	  ogrdb_recs   annotation output file (.tsv)
	  ogrdb_vs     v_gene sequences (.fasta)

	optional arguments:
	  -h, --help   show this help message and exit

The following shell commands will annotate with partis in both formats, merge the output, and process with OGRDBstats:

.. code-block:: bash

	# Run partis to produce annotations in YAML format
	partis annotate --extra-annotation-columns cdr3_seqs:invalid:in_frames:stops --infname TW01A.fasta --outfname TW01A.yaml --n-procs 5
	# Run partis again with additional --presto-output option. This will produce TSV-formatted output from cached data
	partis annotate --extra-annotation-columns cdr3_seqs:invalid:in_frames:stops --infname TW01A.fasta --outfname TW01A.tsv --presto-output \
	 --aligned-germline-fname IMGT_REF_GAPPED_DEDUPED.fasta --n-procs 5
	# Extract and merge required information from YAML and TSV files
	python convert_partis.py TW01A.yaml TW01A.tsv TW01A_OGRDB.tsv TW01A_V_OGRDB.fasta
	# Process the resulting output to produce the genotye file and plots
	Rscript ogrdbstats.R --inf_file TW01A_V_OGRDB.fasta IMGT_REF_GAPPED.fasta Homosapiens TW01A_OGRDB.tsv VH
	
Although partis must be run twice - once without the presto-output option, and once with it - it will use cached information provided 
other parameters remain the same, so that the overall impact on run time is low. 

Partis by default uses germline sequences supplied
with the package. `–-presto-output` requires an IMGT-aligned V-gene germline file, which you may be able to obtain from the same source
that partis uses.  Otherwise, please see :ref:`gapped_v`. This describes a tool that can create a gapped file from the ungapped sequences.
partis will report as an error any duplicated identical sequences in the file: duplicates must be 
removed before processing will complete successfully. The `identical_sequences <https://williamdlees.github.io/receptor_utils/_build/html/gap_sequences.html>`_ 
command from the `receptor_utils <https://williamdlees.github.io/receptor_utils/_build/html/fix_macaque_gaps.html>`_ package will identify duplicates.


