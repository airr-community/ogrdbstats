# This software is licensed under the CC BY-SA 4.0 licence: https://creativecommons.org/licenses/by-sa/4.0/

# Command and test drivers for genotype_statistics

#' Collect parameters from the command line and use them to create a report and CSV file for OGRDB
#' @param args A string vector containing the command line arguments. If NULL, will take them from the command line
#' @return Nothing
#' @export
#' @examples
#' # Prepare files for example
#' reference_set = system.file("extdata/ref_gapped.fasta", package = "ogrdbstats")
#' inferred_set = system.file("extdata/novel_gapped.fasta", package = "ogrdbstats")
#' repertoire = system.file("extdata/repertoire.tsv", package = "ogrdbstats")
#' file.copy(repertoire, tempdir())
#' repfile = file.path(tempdir(), 'repertoire.tsv')
#'
#' genotype_statistics_cmd(c(
#'               reference_set,
#'               'Homosapiens',
#'               repfile,
#'               'IGHV',
#'               '--inf_file', inferred_set))
#'
#'# clean up
#' outfile = file.path(tempdir(), 'repertoire_ogrdb_report.csv')
#' plotfile = file.path(tempdir(), 'repertoire_ogrdb_plots.pdf')
#' file.remove(repfile)
#' file.remove(outfile)
#' file.remove(plotfile)
genotype_statistics_cmd = function(args=NULL) {

  p = argparser::arg_parser('Create genotype statistics')
  p = argparser::add_argument(p, 'REF_FILE', help='reference set filename')
  p = argparser::add_argument(p, 'SPECIES', help=' species name used in field 3 of the IMGT reference set header, with spaces removed, e.g. Homosapiens for Human')
  p = argparser::add_argument(p, 'READ_FILE', help='name of file containing annotated reads in AIRR, CHANGEO, IMPRE or IgDiscover format')
  p = argparser::add_argument(p, 'CHAIN', help='one of IGHV, IGKV, IGLV, IGHD, IGHJ, IGKJ, IGLJ, TRAV, TRAJ, TRBV, TRBD, TRBJ, TRGV, TRGj, TRDV, TRDD, TRDJ')
  p = argparser::add_argument(p, '--inf_file', help='sequences of inferred novel alleles (FASTA format)')
  p = argparser::add_argument(p, '--hap_gene', help='haplotyping gene, e.g. IGHJ6')
  p = argparser::add_argument(p, '--plot_unmutated', flag=TRUE, help='Plot base composition using only unmutated sequences (V-chains only)')
  p = argparser::add_argument(p, '--all_novel', flag=TRUE, help='Treat all alleles in reference set as if novel')


  if (is.null(args)) {
    argv = argparser::parse_args(p, commandArgs(trailingOnly=TRUE))
  } else {
    argv = argparser::parse_args(p, args)
  }

  ref_filename = argv$REF_FILE
  species = argv$SPECIES
  inferred_filename = argv$inf_file
  filename = argv$READ_FILE
  chain = argv$CHAIN
  hap_gene = argv$hap_gene
  plot_unmutated = argv$plot_unmutated
  all_inferred = argv$all_novel

  # convert legacy chain names

  if(chain =='VH') {
    chain = 'IGHV'
  } else if(chain == 'VK') {
    chain = 'IGKV'
  } else if(chain == 'VL') {
    chain = 'IGLV'
  } else if(chain == 'JH') {
    chain = 'IGHJ'
  } else if(chain == 'JK') {
    chain = 'IGKJ'
  } else if(chain == 'JL') {
    chain = 'IGLJ'
  } else if(chain == 'DH') {
    chain = 'IGHD'
  }

  if(!(chain %in% c('IGHV', 'IGKV', 'IGLV', 'IGHD', 'IGHJ', 'IGKJ', 'IGLJ', 'TRAV', 'TRAJ', 'TRBV', 'TRBD', 'TRBJ', 'TRGV', 'TRGj', 'TRDV', 'TRDD', 'TRDJ'))) {
    stop('Unrecognised chain name.')
  }

  if(!(chain %in% c('IGHV', 'IGKV', 'IGLV', 'TRBV', 'TRDV')) && plot_unmutated) {
    stop('The plot_unmutated option can only be used with V chains.')
  }

  segment = substr(chain, 4, 4)

  if(substr(chain, 3, 3) %in% c('H', 'B', 'D')) {
    chain_type = 'H'
  } else {
    chain_type = 'L'
  }

  if(!is.na(inferred_filename)) {
    inferred_filename = argv$inf_file
  } else {
    inferred_filename = '-'
  }

  generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, plot_unmutated, all_inferred)
}


