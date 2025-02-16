.. _overview_label:

Overview
========

Scope and Features
******************

OGRDBstats can provide allele usage statistics and plots for AIRR-seq repertoires.

The statistics provide overall sequence count for each allele, mutation levels, and recombination statistics, for example the number of J-alleles found
to combine with each allele, and the number of associated CDR3 sequences. An example can be found `here <https://github.com/airr-community/ogrdbstats/blob/master/example_ogrdbstats_genotype.csv>`_.

The available plots provide read pile-up information for each allele, and plots showing the mutation levels and counts associated with each assigned allele. Only
reads with an assignment to a single allele are considered in the analysis. If an allele sequence is equally distant to two or more germline alleles, the aligner
will call these as equally likely in the assignment. Such sequences are not considered in the OGRDBstats analysis, as the analysis is focussed on the level of support
for each allele in the annotated repertoire, and such calls are inconclusive.

The statistics and plots can help to assess the overall level of support for the allele, and the likely valididty of the allele call in the repertoire under study.
For example, if a B-receptor allele is only observed at a low level in the annotated repertoire compared to other alleles of the same gene, and the majority of sequences
assigned to it are mutated compared to the germline, one might suspect that the allele is mis-called, and that the assigned sequences represent mutated sequences
that would be more accurately assigned to another germline allele, that is present in the repertoire at a higher level.


.. figure:: read_pile_up.jpg
   :width: 600

   Read pile-up for a specific allele. As well as showing the number of supporting alleles and the degree of diversity at specific positions, this can indicate
   aretefacts that might otherwise be missed. In this case we can see that some supporting reads are incomplete at the 5-prime end, and some sequences contain
   gaps, likely as a result of non-overlapping paired-end reads.


.. figure:: mutation_barplot.jpg
   :width: 600

   Mutation bar plot for a specific allele. The format of this plot is modelled after one used by `IgDiscover <https://igdiscover.se>`_.


Usage with Genotype Inference tools
***********************************

A number of tools have been created to infer a subject's genotype from a repertoire. Because aligners will typically assign each read to the most similar
germline allele, mis-assignments are common, particularly in highly mutated repertoires. These inference tools use a variety of methods to identify
and correct such mis-assignments. This can result in a substantially lower number of alleles being called in the assignments. The set of remaining
alleles represents the subject's inferred genotype. The tools can also infer 'novel' alleles: that is, germline alleles that are not listed in the
germline set, but can be inferred from the repertoire.

The information that OGRDBstats provides can be particularly useful to assess the degree of support for alleles listed in the inferred genotype -
both those listed in the germline set, and those that are inferred. Because the tools provide information in a variety of formats, we have included
specific usage information for tested tools. The tools we have currently tested are:

- `IgDiscover`_
- `TIgGER <https://tigger.readthedocs.io/en/stable/>`_
- `partis <https://github.com/psathyrella/partis>`_
- `IMPre <https://github.com/zhangwei2015/IMPre>`_


References
**********

Please cite this paper if you use OGRDBstats in your work:

Lees et al. 2019. OGRDB: A Reference Database of Inferred Immune Receptor Genes. *Nucleic Acids Research*. `doi: 10.1093/nar/gkz822 <https://doi.org/10.1093/nar/gkz822>`_.

Guidance on the use of OGRDBstats statistics and plots to assess the confidence of allele calls can be found in this publication:

Ohlin et al. 2019. Inferred Allelic Variants of Immunoglobulin Receptor Genes: A System for Their Evaluation, Documentation, and Naming. *Frontiers in Immunology* `doi: 10.3389/fimmu.2019.00435 <https://doi.org/10.3389/fimmu.2019.00435>`_.

Acknowledgements
****************

Example data is taken from:

Rubelt et al. 2016. Individual Heritable Differences Result in Unique Cell Lymphocyte Receptor Repertoires of Naïve and Antigen-Experienced Cells. *Nature Communications*. `doi: 10.1038/ncomms11112 <https://doi.org/10.1038/ncomms11112>`_.

Some functions in OGRDBstats are adapted from `TIgGER	 <https://tigger.readthedocs.io/en/stable/>`_, with thanks to the authors.
