.. _igdiscover:

Using OGRDBstats with IgDiscover
================================

When an IgDiscover run has completed, the results will be found in the `final/ directory <http://docs.igdiscover.se/en/stable/guide.html#final-results>`_. 

The read file file for OGRDBstats analysis is `final.tsv.gz`. This contains the repertoire reads, annotated with the alleles IgDiscover has determined are
present in the personalised germline database. It must be decompressed before use by OGRDBstats:

.. code-block:: bash

    $ unzip final.tsv.gz

OGRDBstats also requires the personalised germline database, containing the alleles IgDiscover has inferred to be in the repertoire. This is located in `final/database`, 
for example `final/database/V.fasta` for V-genes.

For V-gene analysis, OGRDBstats requires the germline set to be aligned according to the `IMGT unique numbering <https://www.imgt.org/IMGTScientificChart/Numbering/IMGTnumbering.html>`_.
All alleles presented to IgDiscover in the `initial V.fasta <http://docs.igdiscover.se/en/stable/guide.html#obtaining-a-v-d-j-database>`_ must be present in this 'gapped' file.

If the `V.fasta` file you used was downloaded from IMGT, you can `download a gapped version <https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP>`_. 
This command downloads the current gapped set from IMGT and names it `IMGT_REF_GAPPED.fasta`:

.. code-block:: bash

    $ wget -O IMGT_REF_GAPPED.fasta https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP

If you do not have ready access to an IMGT-aligned version,
please see :ref:`gapped_v`. This describes a tool that can create a gapped file from your sequences, that would be specified in place of IMGT_REF_GAPPED.fasta in the example below.

The entire IMGT file, containing sequences for multiple species, can be presented to OGRDBstats without further processing. The species must be specified to enable OGRDBstats to find the appripriate records. 
The two words in the species name, for example `Homo sapiens` are run together without spaces. 

The following command, run in the `final/` directory, would process V genes in the IgDiscover output, using the downloaded IMGT gapped sequences:

.. code-block:: bash

    $ Rscript ogrdbstats.R --inf_file database/V.fasta IMGT_REF_GAPPED.fasta Homosapiens filtered.tsv IGHV
	
Output files will be `V_ogrdb_report.csv` and `V_ogrdb_plots.pdf`.
	


Creating reports for D- or J-genes
----------------------------------

If you have `discovered D- or J- genes <http://docs.igdiscover.se/en/stable/changes.html?highlight=discoverj#v0-12-2020-01-20>`_, the process is simpler
as the germline set does not need to be gapped. Again, from the final/ directory:

.. code-block:: bash

    $ Rscript ogrdbstats.R --inf_file database/J.fasta ../database/J.fasta Homosapiens filtered.tab IGHJ

