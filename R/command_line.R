# This software is licensed under the CC BY-SA 4.0 licence: https://creativecommons.org/licenses/by-sa/4.0/

# Command and test drivers for genotype_statistics

#' Collect parameters from the command line and use them to create a report and CSV file for OGRDB
#' @param test Run a simple preconfigured test (files must be in the current directory)
#' @return Nothing
genotype_statistics_cmd = function(test = F) {

  p = arg_parser('Create genotype statistics')
  p = add_argument(p, 'ref_filename', help='reference set filename')
  p = add_argument(p, 'species', help=' species name used in field 3 of the IMGT reference set header, with spaces removed, e.g. Homosapiens for Human')
  p = add_argument(p, 'filename', help='name of file containing annotated reads in AIRR, CHANGEO, IMPRE or IgDiscover format')
  p = add_argument(p, 'chain', help='one of VH, VK, VL, D, JH, JK, JL')
  p = add_argument(p, '--inf', help='sequences of inferred novel alleles (FASTA format)')
  p = add_argument(p, '--hap', help='haplotyping gene, e.g. IGHJ6')


  if(!test) {
    argv = parse_args(p, commandArgs(trailingOnly=TRUE))
  } else {
    argv = parse_args(p, c('IMGT_REF_GAPPED.fasta', 'Homosapiens', 'TWO01A_naive_genotyped.tsv', 'VH', '--inf', 'TWO01A_naive_novel_ungapped.fasta', '--hap', 'IGHJ6'))
    # setwd('D:\\Research\\ogrdbstats\\testdata\\VH_tigger')
  }

  ref_filename = argv$ref_filename
  species = argv$species
  inferred_filename = argv$inf
  filename = argv$filename
  chain = argv$chain
  hap_gene = argv$hap

  if(!(chain %in% c('VH', 'VK', 'VL', 'DH', 'JH', 'JK', 'JL'))) {
    stop('Unrecognised chain name.')
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

  if(!is.na(argv$inf)) {
    inferred_file = argv$inf
  } else {
    inferred_file = '-'
  }

  generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
}

#' Generate OGRDB reports using test data
#' @param full If false, run a single set. Otherwise run all sets we have.
#' @return Nothing
genotype_statistics_test = function(full=F) {

  setwd('D:\\Research\\ogrdbstats')

  # VH - tigger (Example in Readme)
  report('VH - tigger')
  setwd('testdata/VH_tigger')
  ref_filename = 'IMGT_REF_GAPPED.fasta'
  species = 'Homosapiens'
  inferred_filename = 'TWO01A_naive_novel_ungapped.fasta'
  filename = 'TWO01A_naive_genotyped.tsv'
  chain = 'VH'
  hap_gene = 'IGHJ6'
  segment = 'V'
  chain_type = 'H'
  generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
  setwd('D:\\Research\\ogrdbstats')

  if(full) {

    # VH - tigger with truncated sequences in IMGT alignment
    report('VH - tigger with truncated sequences in IMGT alignment')
    setwd('testdata/V_tigger_truncated')
    ref_filename = 'IMGT_REF_GAPPED.fasta'
    species = 'Homosapiens'
    inferred_filename = 'TWO01A_naive_novel_ungapped.fasta'
    filename = 'TWO01A_naive_genotyped.tsv'
    chain = 'VH'
    hap_gene = 'IGHJ6'
    segment = 'V'
    chain_type = 'H'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
    setwd('D:\\Research\\ogrdbstats')

    # JH - tigger
    report('JH - tigger')
    setwd('testdata/JH_tigger')
    ref_filename = 'IMGT_REF_GAPPED_J_CHANGES.fasta'
    species = 'Homosapiens'
    inferred_filename = 'TWO01A_naive_novel.fasta'
    filename = 'TWO01A_naive.airr.tab'
    chain = 'JH'
    hap_gene = 'IGHV2-5'
    segment = 'J'
    chain_type = 'H'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
    setwd('D:\\Research\\ogrdbstats')

    # D - tigger
    report('D - tigger')
    setwd('testdata/D_tigger')
    ref_filename = 'IMGT_REF_GAPPED_D_1_26_01_removed.fasta'
    species = 'Homosapiens'
    inferred_filename = 'TWO01A_naive_novel.fasta'
    filename = 'TWO01A_naive.airr.tab'
    chain = 'DH'
    hap_gene = 'IGHJ6'
    segment = 'D'
    chain_type = 'H'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
    setwd('D:\\Research\\ogrdbstats')

    # JH - igdiscover
    report('JH - igdiscover')
    setwd('testdata/JH_igdiscover')
    ref_filename = 'IMGT_REF_GAPPED_fake_j.fasta'
    species = 'Homosapiens'
    inferred_filename = 'J.fasta'
    filename = 'filtered.tab'
    chain = 'JH'
    hap_gene = 'IGHV2-5'
    segment = 'J'
    chain_type = 'H'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
    setwd('D:\\Research\\ogrdbstats')

    # JL - igdiscover
    # just a fake based on JH
    report('JL - igdiscover')
    setwd('testdata/JL_igdiscover')
    ref_filename = 'IMGT_REF_GAPPED_fake_j_fake_JK.fasta'
    species = 'Homosapiens'
    inferred_filename = 'J.fasta'
    filename = 'filtered.tab'
    chain = 'JK'
    hap_gene = 'IGHV2-5'
    segment = 'J'
    chain_type = 'L'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
    setwd('D:\\Research\\ogrdbstats')

    # VH - partis
    # This takes a lot of memory and does not run comfortably in RStudio
    # report('VH - partis')
    # setwd('testdata/VH_partis')
    # ref_filename = 'IMGT_REF_GAPPED.fasta'
    # species = 'Homosapiens'
    # inferred_filename = 'TW02A_V_OGRDB.fasta'
    # filename = 'TW02A_OGRDB.tsv'
    # chain = 'VH'
    # hap_gene = 'IGHJ6'
    # segment = 'V'
    # chain_type = 'H'
    # generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
    # setwd('D:\\Research\\ogrdbstats')

    # JK - igDiscover
    report('JK - igdiscover')
    setwd('testdata/private/PRJEB30386 - Kappa')
    ref_filename = 'IMGT_REF_GAPPED.fasta'
    species = 'Homosapiens'
    inferred_filename = 'Inferred_file.fasta'
    filename = 'Read_file.tab'
    chain = 'VK'
    hap_gene = 'IGHJ6'
    segment = 'V'
    chain_type = 'L'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type)
    setwd('D:\\Research\\ogrdbstats')
  }

  setwd('../..')
}
