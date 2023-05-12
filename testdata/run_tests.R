# Standalone script to run ogrdb against test data


library(ogrdbstats)

#argv = parse_args(p, c('IMGT_REF_GAPPED.fasta', 'Homosapiens', 'TWO01A_naive_genotyped.tsv', 'IGHV', '--inf_file', 'TWO01A_naive_novel_ungapped.fasta', '--hap_gene', 'IGHJ6', '--plot_unmutated'))
#setwd('D:\\Research\\ogrdbstats\\testdata\\VH_tigger')

#argv = parse_args(p, c('..\\..\\IMGT_REF_GAPPED.fasta', 'Homosapiens', 'filtered.tab', 'IGLV', '--inf_file', 'database/V.fasta'))
#setwd('D:\\Research\\from_martin_corcoran\\VL_second_analysis_2-14_combined\\final')

#argv = parse_args(p, c('IMGT_REF_GAPPED.fasta', 'Homosapiens', 'P8_I1_S1_airr.tsv', 'IGHV'))
#setwd('D:\\Research\\ogrdbstats\\testdata\\VH_airr_no_novel')

#argv = parse_args(p, c('IMGT_REF_GAPPED.fasta', 'Homosapiens', 'TW02A_OGRDB.tsv', 'IGHV', '--inf_file', 'TW02A_V_OGRDB.fasta', '--hap_gene', 'IGHJ6'))
#setwd('D:\\Research\\ogrdbstats\\testdata\\VH_partis')

#argv = parse_args(p, c('IMGT_REF_GAPPED.fasta', 'Homosapiens', 'Read_file.tab', 'IGKV', '--inf_file', 'Inferred_file.fasta', '--hap_gene', 'IGHJ6'))
#setwd('D:\\Research\\ogrdbstats\\testdata\\private\\PRJEB30386 - Kappa')

#argv = parse_args(p, c('human_gl_IGHmakedb_F+ORF+in-frame_P.fasta', 'Homosapiens', 'P1_I1_S1.tsv', 'IGHV', '--inf_file', 'P1_I1_S1_novel_gapped.fasta', '--hap_gene', 'IGHJ6'))
#setwd('D:\\Research\\ogrdbstats\\testdata\\private\\ogrdbstats_in_vdjbase')

#setwd('D:\\Research\\ogrdbstats\\testdata\\private\\hamster')
#argv = parse_args(p, c('hamster_IGH_VDJ.fasta', 'Hamster', 'igblast_001_db-pass.tsv', 'IGHV', '--all_novel'))

# setwd('D:\\Research\\ogrdbstats\\testdata\\trb')
# argv = parse_args(p, c('imgt_gapped.fasta', 'Homosapiens', 'SC9.tab', 'TRBV', '--inf_file', 'SC9_novel.fasta'))

# setwd('D:\\Research\\new_ham\\16. Leaders')
# argv = parse_args(p, c('ref/hamster_IGH_VDJ.fasta', 'hamster', 'leader_allele_matches.tsv', 'VH', '--inf_file', 'leader_aliases.fasta'))

# setwd('D:\\Research\\new_ham\\2. j_allele_search\\igblast')
# argv = parse_args(p, c('../ref/hamster_IGH_VDJ.fasta', 'hamster', 'igblast_009_x-clones.tsv', 'DH', '--all_novel'))

# setwd('D:\\Research\\ogrdbstats\\testdata\\trb2')
# argv = parse_args(p, c('personal_repo.fasta', 'Homosapiens', 'HIP00110_genotyped.tab', 'TRBV', '--inf_file', 'only_novel_repo.fasta'))

# setwd('D:\\Research\\ogrdbstats\\testdata\\private\\JK_changeo')
# argv = parse_args(p, c('IMGT-IGKJ-inc-novels.fasta', 'mouse', '129_igblast_db-pass.tsv', 'JK', '--all_novel'))

setwd('D:\\Research\\ogrdbstats\\testdata\\private\\ogrdbstats_in_vdjbase')
genotype_statistics_cmd(c('IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP', 'Homosapiens', 'P1_I28_S1.tsv', 'IGHV', '--hap_gene', 'IGHJ6', '--plot_unmutated', '--inf_file', 'P1_I28_novels.fasta'))


#' Generate OGRDB reports using test data
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
  generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, FALSE)
  setwd('..')

  if(full) {

    # VH - AIRR with no novel sequences
    report('VH - AIRR no novel')
    setwd('VH_airr_no_novel')
    ref_filename = 'IMGT_REF_GAPPED.fasta'
    species = 'Homosapiens'
    inferred_filename = '-'
    filename = 'P8_I1_S1_airr.tsv'
    chain = 'IGHV'
    hap_gene = 'IGHJ6'
    segment = 'V'
    chain_type = 'H'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, FALSE)
    setwd('..')
    gc()

    # VH - tigger with truncated sequences in IMGT alignment
    report('VH - tigger with truncated sequences in IMGT alignment')
    setwd('V_tigger_truncated')
    ref_filename = 'IMGT_REF_GAPPED.fasta'
    species = 'Homosapiens'
    inferred_filename = 'TWO01A_naive_novel_ungapped.fasta'
    filename = 'TWO01A_naive_genotyped.tsv'
    chain = 'IGHV'
    hap_gene = 'IGHJ6'
    segment = 'V'
    chain_type = 'H'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, FALSE)
    setwd('..')
    gc()

    # JH - tigger
    report('JH - tigger')
    setwd('JH_tigger')
    ref_filename = 'IMGT_REF_GAPPED_J_CHANGES.fasta'
    species = 'Homosapiens'
    inferred_filename = 'TWO01A_naive_novel.fasta'
    filename = 'TWO01A_naive.airr.tab'
    chain = 'IGHJ'
    hap_gene = 'IGHV2-5'
    segment = 'J'
    chain_type = 'H'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, FALSE)
    setwd('..')
    gc()

    # D - tigger
    report('D - tigger')
    setwd('D_tigger')
    ref_filename = 'IMGT_REF_GAPPED_D_1_26_01_removed.fasta'
    species = 'Homosapiens'
    inferred_filename = 'TWO01A_naive_novel.fasta'
    filename = 'TWO01A_naive.airr.tab'
    chain = 'IGHD'
    hap_gene = 'IGHJ6'
    segment = 'D'
    chain_type = 'H'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, FALSE)
    setwd('..')
    gc()

    # JH - igdiscover
    report('JH - igdiscover')
    setwd('JH_igdiscover')
    ref_filename = 'IMGT_REF_GAPPED_fake_j.fasta'
    species = 'Homosapiens'
    inferred_filename = 'J.fasta'
    filename = 'filtered.tab'
    chain = 'IGHJ'
    hap_gene = 'IGHV2-5'
    segment = 'J'
    chain_type = 'H'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, FALSE)
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
    chain = 'IGKJ'
    hap_gene = 'IGHV2-5'
    segment = 'J'
    chain_type = 'L'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, FALSE)
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
    # generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, FALSE)
    # setwd('..')
    gc()

    # JK - igDiscover
    report('JK - igdiscover')
    setwd('private/PRJEB30386 - Kappa')
    ref_filename = 'IMGT_REF_GAPPED.fasta'
    species = 'Homosapiens'
    inferred_filename = 'Inferred_file.fasta'
    filename = 'Read_file.tab'
    chain = 'IGKJ'
    hap_gene = 'IGHJ6'
    segment = 'V'
    chain_type = 'L'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, FALSE)
    setwd('..')

    # VH - IgDiscover 1.51
    # Recent versions of IgDiscover use IgBLAST output in filtered.tsv

    # VH - IgDiscover 1.51
    # Recent versions of IgDiscover use IgBLAST output in filtered.tsv
    report('VK_igdiscover_151')
    setwd('VK_igdiscover_151')
    ref_filename = 'V_ref_gapped.fasta'
    species = 'Homosapiens'
    inferred_filename = 'V_head.fasta'
    filename = 'filtered_head.tsv'
    inferred_filename = 'V.fasta'
    filename = 'filtered.tsv'
    chain = 'IGKV'
    hap_gene = 'IGKJ3'
    segment = 'V'
    chain_type = 'L'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, FALSE)
    setwd('..')


    report('VL_igdiscover_151')
    setwd('VL_igdiscover_151')
    ref_filename = 'V_gapped.fasta'
    species = 'Homosapiens'
    #inferred_filename = 'V_head.fasta'
    #filename = 'filtered_1_01.tsv'
    inferred_filename = 'V.fasta'
    filename = 'filtered.tsv'
    chain = 'IGLV'
    hap_gene = ''
    segment = 'V'
    chain_type = 'L'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, FALSE)
    setwd('..')

    report('VL_mouse_igblast')
    setwd('VL_mouse_igblast')
    ref_filename = 'V_gapped.fasta'
    species = 'Mouse'
    filename = 'annotated_head.tsv'
    #filename = 'annotated_full.tsv'
    inferred_filename = '-'
    chain = 'IGLV'
    hap_gene = NA
    segment = 'V'
    chain_type = 'L'
    generate_ogrdb_report(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, TRUE)
    setwd('..')

  }
}
