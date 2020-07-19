# This software is licensed under the CC BY-SA 4.0 licence: https://creativecommons.org/licenses/by-sa/4.0/

# Command and test drivers for genotype_statistics

#' Collect parameters from the command line and use them to create a report and CSV file for OGRDB
#' @param test Run a simple preconfigured test (files must be in the current directory)
#' @return Nothing
#' @export
genotype_statistics_cmd = function(test = F) {

  p = arg_parser('Create genotype statistics')
  p = add_argument(p, 'REF_FILE', help='reference set filename')
  p = add_argument(p, 'SPECIES', help=' species name used in field 3 of the IMGT reference set header, with spaces removed, e.g. Homosapiens for Human')
  p = add_argument(p, 'READ_FILE', help='name of file containing annotated reads in AIRR, CHANGEO, IMPRE or IgDiscover format')
  p = add_argument(p, 'CHAIN', help='one of VH, VK, VL, D, JH, JK, JL')
  p = add_argument(p, '--inf_file', help='sequences of inferred novel alleles (FASTA format)')
  p = add_argument(p, '--hap_gene', help='haplotyping gene, e.g. IGHJ6')
  p = add_argument(p, '--plot_unmutated', flag=T, help='Plot base composition using only unmutated sequences (V-chains only)')


  if(!test) {
    argv = parse_args(p, commandArgs(trailingOnly=TRUE))
  } else {
    #argv = parse_args(p, c('IMGT_REF_GAPPED.fasta', 'Homosapiens', 'TWO01A_naive_genotyped.tsv', 'VH', '--inf_file', 'TWO01A_naive_novel_ungapped.fasta', '--hap_gene', 'IGHJ6', '--plot_unmutated'))
    #setwd('D:\\Research\\ogrdbstats\\testdata\\VH_tigger')

    #argv = parse_args(p, c('..\\..\\IMGT_REF_GAPPED.fasta', 'Homosapiens', 'filtered.tab', 'VL', '--inf_file', 'database/V.fasta'))
    #setwd('D:\\Research\\from_martin_corcoran\\VL_second_analysis_2-14_combined\\final')

    argv = parse_args(p, c('IMGT_REF_GAPPED.fasta', 'Homosapiens', 'P8_I1_S1_airr.tsv', 'VH'))
    setwd('D:\\Research\\ogrdbstats\\testdata\\VH_no_novel')

  }

  ref_filename = argv$REF_FILE
  species = argv$SPECIES
  inferred_filename = argv$inf_file
  filename = argv$READ_FILE
  chain = argv$CHAIN
  hap_gene = argv$hap_gene
  plot_unmutated = argv$plot_unmutated

  if(!(chain %in% c('VH', 'VK', 'VL', 'DH', 'JH', 'JK', 'JL'))) {
    stop('Unrecognised chain name.')
  }

  if(!(chain %in% c('VH', 'VK', 'VL')) && plot_unmutated) {
    stop('The plot_unmutated option can only be used with V chains.')
  }

  segment = substr(chain, 1, 1)

  if(segment == 'D') {
    chain_type = 'H'
  } else {
    if(substr(chain, 2, 2) == 'H') {
      chain_type = 'H'
    } else {
      chain_type = 'L'
    }
  }

  if(!is.na(inferred_filename)) {
    inferred_filename = argv$inf_file
  } else {
    inferred_filename = '-'
  }

  generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, plot_unmutated)
}

#' Generate OGRDB reports using test data
#' @param testdir Directory containing test data sets
#' @param full If false, run a single set. Otherwise run all sets we have.
#' @return Nothing
#' @export
genotype_statistics_test = function(testdir, full=F) {
  Sys.setenv('R_MAX_VSIZE'=32000000000)  # this setting may not be effective in Rstudio
  setwd(testdir)

  # VH - tigger (Example in Readme)
  report('VH - tigger')
  setwd('VH_tigger')
  ref_filename = 'IMGT_REF_GAPPED.fasta'
  species = 'Homosapiens'
  inferred_filename = 'TWO01A_naive_novel_ungapped.fasta'
  filename = 'TWO01A_naive_genotyped.tsv'
  chain = 'VH'
  hap_gene = 'IGHJ6'
  segment = 'V'
  chain_type = 'H'
  generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
  setwd('..')

  if(full) {

    # VH - tigger with truncated sequences in IMGT alignment
    report('VH - tigger with truncated sequences in IMGT alignment')
    setwd('V_tigger_truncated')
    ref_filename = 'IMGT_REF_GAPPED.fasta'
    species = 'Homosapiens'
    inferred_filename = 'TWO01A_naive_novel_ungapped.fasta'
    filename = 'TWO01A_naive_genotyped.tsv'
    chain = 'VH'
    hap_gene = 'IGHJ6'
    segment = 'V'
    chain_type = 'H'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
    setwd('..')
    gc()

    # JH - tigger
    report('JH - tigger')
    setwd('JH_tigger')
    ref_filename = 'IMGT_REF_GAPPED_J_CHANGES.fasta'
    species = 'Homosapiens'
    inferred_filename = 'TWO01A_naive_novel.fasta'
    filename = 'TWO01A_naive.airr.tab'
    chain = 'JH'
    hap_gene = 'IGHV2-5'
    segment = 'J'
    chain_type = 'H'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
    setwd('..')
    gc()

    # D - tigger
    report('D - tigger')
    setwd('D_tigger')
    ref_filename = 'IMGT_REF_GAPPED_D_1_26_01_removed.fasta'
    species = 'Homosapiens'
    inferred_filename = 'TWO01A_naive_novel.fasta'
    filename = 'TWO01A_naive.airr.tab'
    chain = 'DH'
    hap_gene = 'IGHJ6'
    segment = 'D'
    chain_type = 'H'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
    setwd('..')
    gc()

    # JH - igdiscover
    report('JH - igdiscover')
    setwd('JH_igdiscover')
    ref_filename = 'IMGT_REF_GAPPED_fake_j.fasta'
    species = 'Homosapiens'
    inferred_filename = 'J.fasta'
    filename = 'filtered.tab'
    chain = 'JH'
    hap_gene = 'IGHV2-5'
    segment = 'J'
    chain_type = 'H'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
    setwd('..')
    gc()

    # JL - igdiscover
    # just a fake based on JH
    report('JL - igdiscover')
    setwd('JL_igdiscover')
    ref_filename = 'IMGT_REF_GAPPED_fake_j_fake_JK.fasta'
    species = 'Homosapiens'
    inferred_filename = 'J.fasta'
    filename = 'filtered.tab'
    chain = 'JK'
    hap_gene = 'IGHV2-5'
    segment = 'J'
    chain_type = 'L'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
    setwd('..')
    gc()

    # VH - partis
    # This takes a lot of memory and does not run comfortably in RStudio
    # report('VH - partis')
    # setwd('VH_partis')
    # ref_filename = 'IMGT_REF_GAPPED.fasta'
    # species = 'Homosapiens'
    # inferred_filename = 'TW02A_V_OGRDB.fasta'
    # filename = 'TW02A_OGRDB.tsv'
    # chain = 'VH'
    # hap_gene = 'IGHJ6'
    # segment = 'V'
    # chain_type = 'H'
    # generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
    # setwd('..')
    gc()

    # JK - igDiscover
    report('JK - igdiscover')
    setwd('private/PRJEB30386 - Kappa')
    ref_filename = 'IMGT_REF_GAPPED.fasta'
    species = 'Homosapiens'
    inferred_filename = 'Inferred_file.fasta'
    filename = 'Read_file.tab'
    chain = 'VK'
    hap_gene = 'IGHJ6'
    segment = 'V'
    chain_type = 'K'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
    setwd('..')
  }
}
