.. _gapped_v:

Creating and using a gapped V germline set
===========================================

If the V-sequences you annotated with do not contain gaps, and you do not have ready access to an IMGT-aligned version, you can create a
one using the `gap_sequences <https://williamdlees.github.io/receptor_utils/_build/html/gap_sequences.html>`_ command from the `receptor_utils <https://williamdlees.github.io/receptor_utils/_build/html/fix_macaque_gaps.html>`_ package. 
Installation requires a recent version of Python. `gap_sequences` requires a gapped set of sequences to use as a template. The receptor_utils package 
can download a gapped set from IMGT. The following shell commands will install the package, download a gapped human IGHV set to use as a template, and create 
a gapped set of sequences:

.. code-block:: bash

	# install receptor_utils latest version (requires Python)
	pip install -U receptor_utils
	
	# download human IGHV germline set from IMGT
	extract_refs -L IGH "Homo sapiens" 
	
	# create a gapped version of V.fasta, using the downloaded gapped germline set 
	# (Homo_sapiens_IGHV_gapped.fasta) as a template
	gap_sequences V.fasta Homo_sapiens_IGHV_gapped.fasta V_gapped.fasta
	
The downloaded germline set is used only as a template to determine the gapping. The resulting gapped reference set will contain only the sequences that you provide to be gapped.
You will be warned if any sequences in your set appear to be non-functional. They will still be gapped.

You may see warnings relating to a few sequences present in the template file. These can be ignored: `gap_sequences` will use complete and 
functionally correct sequences in the database as a template and ignore the ones it warns about.


Creating and using a gapped V germline set for Rhesus macaque
---------------------------------------------------------------

IMGT's alignment of rhesus macaque IGH, K and L germline sequences is non-canonical: all three contain inserted codons relative to those defined in the `IMGT unique 
numbering <https://www.imgt.org/IMGTScientificChart/Numbering/IMGTnumbering.html>`_. These cause problems for the functionality checking built in to OGRDBstats, and will cause 
problems with other tools that expect gapped sequences to follow the canonical alignment.

`receptor_utils <https://williamdlees.github.io/receptor_utils/_build/html/fix_macaque_gaps.html>`_ contains a command `fix_macaque_gaps <https://williamdlees.github.io/receptor_utils/_build/html/fix_macaque_gaps.html>`_.
This set of shell commands will download a Rhesus IGH gapped set from IMGT, and convert to the canonical alignment. Receptor_utils requires a recent version of Python.

.. code-block:: bash

	# download rhesus IGHV germline set from IMGT
	extract_refs -L IGH "Rhesus macaque"
	
	# fix the non-canonical alignment
	fix_macaque_gaps Rhesus_macaque_IGHV_gapped.fasta Rhesus_macaque_IGHV_gapped_fixed.fasta IGH
	
The resulting file `Rhesus_macaque_IGHV_gapped_fixed.fasta` will follow the canonical IMGT alignment. `IGH` can be replaced with `IGK` or `IGL`.

`fix_macaque_gaps` will notify you that a small number of sequences have been removed from the set because they contain an inserted codon. These sequences have not, to the author's knowledge at the time of writing, been
reported as being present in significant levels in AIRR-seq repertoires, or demonstrated to form functional antibodies. 