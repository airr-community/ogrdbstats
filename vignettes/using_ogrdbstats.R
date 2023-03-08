## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(eval = FALSE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)

## ----eval=TRUE, fig.width=5, fig.asp=3----------------------------------------
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
gridExtra::grid.arrange(grobs=list(barplot_grobs[3][[1]], base_grobs$end[1][[1]], 
                                   base_grobs$conc[1][[1]]),ncol=1)



