Using ogrdbstats
================
William Lees
2025-02-16

# Conducting a Genotype Analysis with ogrdbstats

[ogrdbstats](https://github.com/airr-community/ogrdbstats) is an R
package that can be used to create an analysis of gene usage in a
adaptive immune receptor sequencing repertoire. The analysis consists of
[usage
statistics](https://github.com/airr-community/ogrdbstats/blob/master/example_ogrdbstats_genotype.csv)
and
[plots](https://github.com/airr-community/ogrdbstats/blob/master/example_ogrdbstats_plots.pdf).
The package is intended to be used in conjunction with a tool that
infers a ‘personalised genotype’. Currently the following tools are
supported:

- [TIgGER](https://tigger.readthedocs.io/en/stable/)
- [IgDiscover](https://igdiscover.se)
- [partis](https://github.com/psathyrella/partis)
- [IMPre](https://github.com/zhangwei2015/IMPre)

### Package Installation

Ogrdbstats requires a recent installation of
[Pandoc](https://pandoc.org/). If Rstudio is installed on your machine,
pandoc will already be installed. Otherwise please follow the
installation instructions on the Pandoc website.

Once Pandoc is installed, please install ogrdbstats from CRAN:

``` r
install.packages('ogrdbstats')
```

Alternatively, you can install the latest development version from
Github:

``` r
devtools::install_github("https://github.com/airr-community/ogrdbstats")
```

The package requires R version 4.2.2 or above.

### Using The Package from the Command Line

It’s easiest to use the package from the command line. To do this,
download
[ogrdbstats.R](https://raw.githubusercontent.com/airr-community/ogrdbstats/master/ogrdbstats.R)
and copy to the directory in which you wish to conduct the analysis.

*Command Syntax (see following section for detailed description of file
formats)*

``` bash
Rscript ogrdbstats.R [--inf_file INF_FILE] [--hap_gene HAP_GENE] REF_FILE SPECIES READ_FILE CHAIN
```

Positional Arguments:

`REF_FILE` - pathname of a FASTA file containing IMGT gap-aligned
reference germline sequences. Usually this would be [downloaded from
IMGT](https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP).

`SPECIES` - should contain the species name used in field 3 of the IMGT
REF_FILE FASTA header, with spaces removed, e.g. Homosapiens for Human.
If you are not using an IMGT REF_FILE, you can use any single word for
the species here, and the reference file should only contain genes for
that species.

`READ_FILE` - pathname of a tab-separated file containing the annotated
reads used to infer the genotype, in MiAIRR, CHANGEO or IgDiscover
format

`CHAIN` specifies the sequence type to be analysed. It must be one of
`VH, VK, VL, D, JH, JK, JL`.

Optional Arguments:

`INF_FILE` - pathname of a FASTA file containing sequences of inferred
novel alleles. This file must be provided if the read file contains
assignments to alleles that are not listed in the reference file.

`HAP_GENE` - the gene to be used for haplotyping analysis (see
haplotyping section below)

Detailed descriptions of the required input files are given in the next
section, but for quick usage with a supported tool, please skip to the
Usage Notes for that tool towards the end of this document.

The package provides functions to allow you to create the reports
programmatically, or use the information it reads from input files for
your own purposes. Please refer to the ‘ogrdbstats API’ section for a
brief overview.

### Detailed Description of Input Files

#### REF_FILE - FASTA file containing the IMGT gap-aligned reference germline sequences.

- All germlines that are called in the read file (apart from those of
  novel alleles) should be included.
- The sequences must be IMGT gap-aligned
- The FASTA header can either be in IMGT’s germline library format, or
  simply consist of the allele name
- The IMGT set can be
  [downloaded](https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP)
  and used as-is: the script will filter out the records for the
  nominated species. As the IMGT set changes from time to time, please
  make sure that the same version is used by the inference tool and by
  this script.
- A warning will be given if any calls in the read file do not have a
  corresponding sequence in the reference file or the inferred novel
  alleles file. Unmutated counts will not be provided for those
  sequences.

#### READ_FILE - A tab-separated file containing the annotated reads used to infer the genotype, in MiAIRR, CHANGEO or IgDiscover format.

- The format will be determined automatically by the script.

- MiAIRR format files must contain at least the following columns:
  `sequence_id, v_call_genotyped, d_call, j_call, sequence_alignment, cdr3`.
  For J or D inferences they must also contain `J_sequence_start`,
  `J_sequence_end`, `J_germline_start`, `J_germline_end`, or the
  equivalent fields for D genes.
  [IgBLAST](https://www.ncbi.nlm.nih.gov/igblast/)’s `--format airr`
  creates compatible MiAIRR format files.

- CHANGEO files must contain at least the following columns:
  `SEQUENCE_ID, V_CALL_GENOTYPED, D_CALL, J_CALL, SEQUENCE_IMGT, CDR3_IMGT`,
  `V_MUT_NC`, `D_MUT_NC`, `J_MUT_NC`, `SEQUENCE`, `JUNCTION_START`,
  `V_SEQ`, `D_SEQ`, `J_SEQ`. If you would like to process files from
  IMGT V-Quest, please [parse them with
  CHANGEO](https://changeo.readthedocs.io/en/stable/examples/imgt.html)
  to convert them to CHANGEO format.

- D- related fields are only required for heavy chain records.

- In both the above file formats, `v_call_genotyped/V_CALL_GENOTYPED`
  should contain the V calls made after the subject’s V-gene genotype
  has been inferred (including calls of the novel alleles). Sequences
  may be either unagpped or IMGT gap-aligned. Determining the
  personalised V-gene genotype is recommended when processing D or J
  gene inferences, so that V-gene usage counts are accurate. However,
  this step can be omitted for D or J gene processing by providing a
  V_CALL field instead of V_CALL_GENOTYPED.

- For IgDiscover, the file ‘final/filtered.tab’ should be used - see
  section on IgDiscover below.

#### INF_FILE - FASTA file containing the inferred novel alleles

- Sequences in the inferred file should all be of the same type: VH, VK,
  VL, D, JH, JK, or JL
- The header should simply consist of the allele name as assigned by the
  tool.
- V-gene sequences may either be IMGT gap-aligned or not aligned. If
  they are not aligned, the script will determine the nearest reference
  gene and use it as a template. If you are not satisfied with the
  resulting alignment, just align the sequence in the inferred file as
  you prefer.
- If a gene with the same name is present in both the germline file and
  the inferred file, its presence in the inferred file will be ignored.
  This makes it easier to use the script with inference tools that do
  not write the inferred sequences to a separate file.

### Output

- `<READ_FILE>_ogrdb_report.csv` - the Genotype File ready to be
  uploaded to OGRDB.
- `<READ_FILE>_ogrdb_plots.csv` - plots (see next section for details).

`READ_FILE` is used as a prefix to the output file names. They will be
written to the directory containing the read file. If you dould like to
specify a different prefix, use the `--file_prefix` option.

If you are submitting inferences to OGRDB, you will be prompted to
upload the genotype file. Please also upload the plots file as an
attachment to the Notes section of your submission.

### Plots

The script produces the following plots:

- For each allele used in the the read file, a histogram showing the
  number of mutated and unmutated sequences
- Barcharts showing nucleotide usage at locations in the IMGT-aligned
  sequence: both across the sequence as a whole, and in more detail at
  the 3’ end
- A barchart showing usage of the alleles of potential haplotyping
  genes, across the whole genotype. This can be used to identify a
  suitable gene for haplotyping analysis.
- For each potential haplotyping candidate (selected by the tool from
  the usage chart above), a plot comparing the usage of the two most
  frequently used alleles of that gene.

The nucleotide usage plots are not produced from IgDiscover output, as
aligned V-sequences are not available.

### Haplotyping

The script should first be run without the optional `HAP_GENE`
parameter. If, having consulted the plots, you identify a suitable gene
for haplotyping, please run the script again, with this gene specified
as `HAP_GENE`. The haplotyping_gene and haplotyping_ratio columns of the
genotype file will be appropriately populated. A J-gene should be used
with V- and D- gene inferences, and a V-gene with J-gene inferences.

### Example

To download the IMGT reference file and complete an analysis using
[example
files](https://github.com/airr-community/ogrdbstats/tree/master/testdata/VH_tigger),
run the following commands:

``` bash
wget -O IMGT_REF_GAPPED.fasta https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP
Rscript ogrdbstats.R --inf_file TWO01A_naive_novel_ungapped.fasta --hap_gene IGHJ6 IMGT_REF_GAPPED.fasta Homosapiens TWO01A_naive_genotyped.tsv VH
```

### Usage Notes

Usage notes are indicative only and are not intended to discount other
approaches. Notes for other tools will follow.

### Usage Notes - TIgGER

To conduct a V-gene analysis with TIgGER:

- Use `findNovelAlleles` to identify novel alleles in a
  Change-O-formatted data set. Write these to a FASTA file.
- Use `inferGenotype` or `inferGenotypeBayesian` to infer the genotype.
- Use `reassignAlleles` to correct allele calls in the data set, based
  on the inferred genotype
- Provide the resulting Change-O file, together with the FASTA file
  containing the novel alleles, to `ogrdbstats`.

Note that `inferGenotype` will not necessarily include every inferred
allele produced by `findNovelAlleles` in the genotype that it produces.
Only those alleles included in the genotype will be considered by
`genotype_statistics.R` because, leaving other considerations aside, no
sequences are assigned to other alleles.

TIgGER provides additonal information, including its own plots and
statistics We encourage you to take these into consideration, and to
upload them as attachments to your submission if they are informative.

### Usage Notes - IgDiscover

Following an IgDiscover run, please copy ogrdbstats.R to IgDiscover’s
`final` directory. The commands below can then be used to download the
IMGT reference file and run a VH gene analysis. All commands should be
run in the `final` directory.

``` bash
$ wget -O IMGT_REF_GAPPED.fasta https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP
$ unzip final.tab.gz
$ Rscript ogrdbstats.R --inf_file database/V.fasta IMGT_REF_GAPPED.fasta Homosapiens filtered.tab VH
```

alternatively, to produce a JH gene analysis:

``` bash
$ Rscript ogrdbstats.R --inf_file database/J.fasta IMGT_REF_GAPPED.fasta Homosapiens filtered.tab JH
```

### Usage Notes - partis

The information required by generate_statstics.R is split between
partis’s normal yaml output and that provided by the ‘presto-output’
mode. A python script,
[convert_partis.py](https://github.com/airr-community/ogrdbstats/blob/master/convert_partis.py),
is provided. This will combine output from partis’s yaml and presto
annotations, producing CHANGEO format annotations and a FASTA file of
genotype V-sequences. These files can then passed to
generate_statistics.R. convert_partis.py is written in python 2.7 for
compatibility with partis, and can be run from the command line in the
same virtual environment.

Usage of convert_partis.py:

``` bash
python convert_partis.py [-h] partis_yaml partis_tsv ogrdb_recs ogrdb_vs

positional arguments:
  partis_yaml  .yaml file created by partis
  partis_tsv   .tsv file created by partis presto-output mode
  ogrdb_recs   annotation output file (.tsv)
  ogrdb_vs     v_gene sequences (.fasta)

optional arguments:
  -h, --help   show this help message and exit
```

Although partis must be run twice - once without the presto-output
option, and once with it - it will use cached information provided other
parameters remain the same, so that the overall impact on run time is
low. Typical processing steps are shown below. Note that –presto-output
requires an IMGT-gapped V-gene germline file. This can be extracted from
the full germline library downloaded from IMGT (see ‘Prerequisites’
above), but partis will report as an error any duplicated identical
sequences in the library: duplicates must be removed from the file
before processing will complete successfully. Note in the examples below
the –extra-annotations option. The CDR3 is required by
generate_statistics.R: the other fields are included for reference.

``` bash
# Run partis to produce annotations in YAML format
partis annotate --extra-annotation-columns cdr3_seqs:invalid:in_frames:stops --infname TW01A.fasta --outfname TW01A.yaml --n-procs 5
# Run partis again with additional --presto-output option. This will produce TSV-formatted output from cached data
partis annotate --extra-annotation-columns cdr3_seqs:invalid:in_frames:stops --infname TW01A.fasta --outfname TW01A.tsv --presto-output \
 --aligned-germline-fname IMGT_REF_GAPPED_DEDUPED.fasta --n-procs 5
# Extract and merge required information from YAML and TSV files
python convert_partis.py TW01A.yaml TW01A.tsv TW01A_OGRDB.tsv TW01A_V_OGRDB.fasta
# Process the resulting output to produce the genotye file and plots
Rscript ogrdbstats.R --inf_file TW01A_V_OGRDB.fasta IMGT_REF_GAPPED.fasta Homosapiens TW01A_OGRDB.tsv VH
```

### Usage Notes - IMPre

IMPre does not provide a set of sequences annotated with the novel
allele calls. The sequences must be annotated by a separate tool in
order to provide the information needed for the OGRDB genotype. One
possible approach is as follows:

- Annotate with IgBLAST using a custom germline set that includes the
  novel alleles inferred by IMPre. Details for creating IgBLAST’s
  germline database are given in the [setup
  notes](https://ncbi.github.io/igblast/cook/How-to-set-up.html). The
  novel alleles inferred by IMPre can be added to the germline sequences
  downloaded from IMGT, before running makeblastdb.
- Select the AIRR output format, and provide the resulting annotion file
  to genotype_statistics.R, along with the novel inferences provided by
  IMPre

### ogrdbstats API

#### Creating reports programmatically

`generate_ogrdb_report()` takes equivalent arguments to ogrdbstats.R and
generates the genotype file and plots.

#### Customising and Extending the Reports

`read_input_files()` reads and parses the input files. It returns a
representation of the genotype file in memory, as well as a number of
other structures containing related information.

`make_barplot_grobs()`, `make_novel_base_grobs()` and
`make_haplo_grobs()` create the plots: each function returns a list of
grobs that can be used as you wish, or passed to `write_plot_file()` to
create a the standard file of plots.

The use of the grob functions is demonstrated below with example data
included in the package.

``` r
library(ogrdbstats)


reference_set = system.file("extdata/ref_gapped.fasta", package = "ogrdbstats")
inferred_set = system.file("extdata/novel_gapped.fasta", package = "ogrdbstats")
repertoire = system.file("extdata/ogrdbstats_example_repertoire.tsv", package = "ogrdbstats")

rd = suppressMessages(
  read_input_files(reference_set, inferred_set, 'Homosapiens', repertoire, 'IGHV', NA, 'V', 'H', FALSE)
)

barplot_grobs = make_barplot_grobs(rd$input_sequences, rd$genotype_db, rd$inferred_seqs, 
                                   rd$genotype, 'V', rd$calculated_NC)
base_grobs = make_novel_base_grobs(rd$inferred_seqs, rd$input_sequences, 'V', FALSE)
gridExtra::grid.arrange(grobs=list(barplot_grobs[2][[1]], base_grobs$end[1][[1]], 
                                   base_grobs$conc[1][[1]]),ncol=1)
```

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />

Please use help(package=“ogrdbstats”) or ? to find function-level
documentation within R.

### Acknowledgements

Some functions are adapted from [TIgGER](https://tigger.readthedocs.io)
with thanks to the authors.

The example annotated reads and inferences linked from this description
are taken from the data of [Rubelt et
al](https://pubmed.ncbi.nlm.nih.gov/?term=27005435) and were downloaded
from [VDJServer](https://vdjserver.org/). The genotype was inferred by
[TIgGER](https://tigger.readthedocs.io). A small number of light-chain
records were removed from the data set.
