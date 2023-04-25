.. _gapped_v:

Creating and using a gapped V reference file
============================================

If the starting V.fasta file you used with IgDiscover was not based on an IMGT set, and you do not have ready access to an IMGT-aligned version, you can create a
gapped version of V.fasta using the `gap_sequences` command from the `receptor_utils <https://williamdlees.github.io/receptor_utils/_build/html/fix_macaque_gaps.html>`_ package. 
Installation requires a recent version of Python. `gap_sequences` requires a gapped set of sequences to use as a template. The receptor_utils package 
can download a gapped set from IMGT. The following commands will install the package, download a gapped human IGHV set, and use it to create 
a gapped set from the starting V.fasta database (again the commands are written to be run in the final/ directory):

.. code-block:: bash

	$ # install receptor_utils latest version (requires Python)
	$ pip install -U receptor_utils
	
	$ # download human IGHV reference set from IMGT
	$ extract_refs -L IGH "Homo sapiens" 
	
	$ # create a gapped version of V.fasta, using the downloaded gapped reference set 
	$ # (Homo_sapiens_IGHV_gapped.fasta) as a template
	$ # Note that the V.fasta file used here should be the starting database passed to 
	$ # IgDiscover, i.e. database/V.fasta not final/database/V.fasta
	$ gap_sequences ../database/V.fasta V_starting_gapped.fasta
	
You may see some warnings relating to sequences present in Homo_sapiens_IGHV_gapped.fasta. These can be ignored: gap_sequences will use the
functionally correct sequences in the database and ignore the ones it warns about.


Creating and using a gapped V reference file for Rhesus macaque
---------------------------------------------------------------

IMGT's alignment of rhesus macaque IGH,K and L reference sequences is non-canonical: all three alignments contain inserted codons relative to those defined in the `IMGT unique 
numbering <https://www.imgt.org/IMGTScientificChart/Numbering/IMGTnumbering.html>`_. These cause problems for the functionality checking built in to ogrdbstats, and will cause 
problems with other tools that expect gapped sequences to follow the canonical alignment.

`receptor_utils <https://williamdlees.github.io/receptor_utils/_build/html/fix_macaque_gaps.html>`_ contains a command `fix_macaque_gaps <https://williamdlees.github.io/receptor_utils/_build/html/fix_macaque_gaps.html>`_ which will correct the alignment.
This set of commands will download a Rhesus IGH gapped set from IMGT, and fix the alignment. Receptor_utils requires a recent version of Python.
database, such as one based on `kimdb <http://kimdb.gkhlab.se/>`_:

.. code-block:: bash

	$ # install receptor_utils latest version (requires Python)
	$ pip install -U receptor_utils
	
	$ # download rhesus IGHV reference set from IMGT
	$ extract_refs -L IGH "Rhesus macaque"
	
	$ # fix the non-canonical alignment
	$ fix_macaque_gaps Rhesus_macaque_IGHV_gapped.fasta Rhesus_macaque_IGHV_gapped_fixed.fasta IGH
	
The resulting file `Rhesus_macaque_IGHV_gapped_fixed.fasta` will follow the canonical IMGT alignment. `IGH` can be replaced with `IGK` or `IGL`.

`fix_macaque_gaps` will notify you that a small number of sequences have been removed from the set because they contain the inserted codon. These sequences have not, to the author's knowledge at the time of writing, been
reported as being present in significant levels in AIRR-seq repertoires, or demonstrated to form functional antibodies. 