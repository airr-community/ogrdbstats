.. _formats:


Detailed Description of Input Files
===================================

`REF_FILE` - FASTA file containing the IMGT gap-aligned germline sequences.

    - All germlines that are called in the repertoire file (apart from those of novel alleles) should be included.
    - The V-sequences must be aligned according to the `IMGT unique numbering <https://www.imgt.org/IMGTScientificChart/Numbering/IMGTnumbering.html>`_.
    - The must simply consist of the allele name, unless the file is downloaded from IMGT as in the next point
    - The IMGT set can be `downloaded <https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP>`_ and used as-is: the script will filter out the records for the nominated species. As the IMGT set changes from time to time, please make sure that the same version is used by the inference tool and by this script.
    - A warning will be given if any calls in the repertoire file do not have a corresponding sequence in the germline set or the inferred novel alleles file. Unmutated counts will not be provided for those sequences.

`READ_FILE` - A tab-separated file containing annotated reads, in AIRR, CHANGEO or IgDiscover format.

    - The format will be determined automatically by OGRDBstats.
    - AIRR format files must contain at least the following columns: sequence_id, v_call_genotyped, d_call, j_call, sequence_alignment, cdr3. For J or D inferences they must also contain J_sequence_start, J_sequence_end, J_germline_start, J_germline_end, or the equivalent fields for D genes. IgBLAST’s --format airr creates compatible AIRR format files.
    - CHANGEO files must contain at least the following columns: SEQUENCE_ID, V_CALL_GENOTYPED, D_CALL, J_CALL, SEQUENCE_IMGT, CDR3_IMGT, V_MUT_NC, D_MUT_NC, J_MUT_NC, SEQUENCE, JUNCTION_START, V_SEQ, D_SEQ, J_SEQ. 

    - D- related fields are only required for heavy chain records.

    - In AIRR and CHANGEO formats, the v_call_genotyped/V_CALL_GENOTYPED column should contain the V calls made after the subject’s V-gene genotype has been inferred (including calls of the novel alleles). Sequences may be either unagpped or IMGT gap-aligned. If there is no v_call_genotyped/V_CALL_GENOTYPED column, ogrdbstats will fall back to using V-CALL, with a warning.

    - For IgDiscover, `final/filtered.tsv` should be used - see :ref:`igdiscover`

`INF_FILE` - FASTA file containing the inferred novel alleles

    - Sequences in INF_FILE should all be of the same type (V, D or J) and from the same locus. 
    - The header should consist solely of the allele name as assigned by the inference tool.
    - V-gene sequences may either be IMGT gap-aligned or ungapped. If they are ungapped, the script will determine the nearest germline gene and use it as a template. If you are not satisfied with the resulting alignment, just align the sequence in the file as you prefer.
    - If a gene with the same name is present in both the germline file and the inferred file, its presence in the inferred file will be ignored. This makes it easier to use the script with inference tools that do not write the inferred sequences to a separate file.
