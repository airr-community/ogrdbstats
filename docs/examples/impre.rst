.. _impre:

Using OGRDBstats with IMPre
===========================

IMPre does not provide a set of sequences annotated with the novel allele calls. The sequences must be annotated by a 
separate tool in order to provide the information needed for the OGRDB genotype. One possible approach is as follows:

Annotate with IgBLAST using a custom germline set that includes the novel alleles inferred by IMPre. Details for creating 
IgBLAST’s germline database are given in the `setup notes <https://ncbi.github.io/igblast/cook/How-to-set-up.html>`_. 
The novel alleles inferred by IMPre can be added to the germline sequences downloaded from IMGT, before running makeblastdb.

Select the AIRR output format, and provide the resulting annotion file to OGRDBstats, along with the novel 
inferences provided by IMPre.


