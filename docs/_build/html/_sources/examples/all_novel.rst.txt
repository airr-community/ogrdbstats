.. _all_novel:

Producing a more complete report
================================

By default, the first section of the OGRDBstats pdf, and columns C to I in the CSV, contain data for novel alleles only. These can be extended to cover all 
alleles by using the --all_novel flag, which instructs OGRDBstats to treat all alleles as though they are novel, for example:

.. code-block:: bash

   $ Rscript --all_novel ogrdbstats.R v_germline_gapped.fasta Human rep_genotyped.tsv IGHV

Processing time (and the size of the report) may be considerably extended. The report will contain warnings for any alleles present in the germline set that 
are not observed in the genotyped repertoire. These warnings can be ignored, as one would not expect to observe every allele in the germline set in a single repertoire. 





