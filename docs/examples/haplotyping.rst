.. _haplotyping:


Haplotyping
===========

The V, (D) and J genes that rearrange are all taken from the same chromosome. If, say, an individual has two different alleles of the same J-gene,
the alleles of V- and D- genes can be separated into two haplotypes, by considering which J-allele each is seen to rearrange with. Here, the J-gene 
concerned is called the `anchor`.

OGRDBstats will report on haplotyping opportunities: when analysing V or D genes, it looks for potential J anchors. When analysing J genes, it looks
for potential V anchors. The identification is based on allele usage as summarised in the Allele Usage plot in section 3, 
`Allele usage in potential haplotype anchor genes`. Where an anchor gene is discovered, a corresponding haplotype report will be generated
in section 4, `Haplotype plots`.

When calling ogrdbstats, you can also request haplotyping ratios to be added to the csv report, by specifying the gene in the option `--hap_gene`,
for example `--hap_gene IGHJ6`.