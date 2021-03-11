# This software is licensed under the CC BY-SA 4.0 licence: https://creativecommons.org/licenses/by-sa/4.0/
#
# Some functions are adapted from TIgGER (https://tigger.readthedocs.io) with thanks to the authors.


#' @import tigger
#' @import alakazam
#' @import stringr
#' @import grid
#' @import tidyr
#' @import dplyr
#' @import stringdist
#' @import RColorBrewer
#' @import argparser
#' @import ggplot2
#' @import ComplexHeatmap

#' @importFrom Biostrings pairwiseAlignment
#' @importFrom gridExtra grid.arrange
#' @importFrom gridExtra marrangeGrob
#' @importFrom utils write.csv

globalVariables(c('A_CALL', 'CDR3_nt', 'D_CALL', 'D_MUT_NC', 'D_SEQ', 'D_errors', 'D_gene', 'D_region', 'GENE',
                  'J_CALL', 'J_MUT_NC', 'J_SEQ', 'J_errors', 'J_gene', 'J_nt', 'SEG_CALL', 'SEG_MUT_NC', 'VDJ_nt',
                  'V_CALL', 'V_CALL_GENOTYPED', 'V_CDR3_start', 'V_MUT_NC', 'V_SEQ', 'V_errors', 'V_gene',
                  'V_nt', 'a_allele', 'a_gene', 'aa_diff', 'aa_substitutions', 'aggregate',
                  'allelic_percentage', 'assigned_unmutated_frequency', 'closest_host',
                  'closest_reference', 'dev.off', 'haplotyping_gene', 'haplotyping_ratio',
                  'host_aa_difference', 'host_aa_subs', 'host_closest', 'host_difference',
                  'host_nt_diffs', 'name', 'nt_diff', 'nt_diff_host', 'nt_sequence', 'nt_sequence_gapped', 'nt_substitutions',
                  'nuc', 'pdf', 'percent', 'read.delim', 'reference_aa_difference', 'reference_aa_subs',
                  'reference_closest', 'reference_difference', 'reference_nt_diffs', 'sequence_id',
                  'sequences', 'setNames', 'triplet', 'triplet_grobs', 'unique_cdr3s', 'unique_cdr3s_unmutated',
                  'unique_ds', 'unique_ds_unmutated', 'unique_js', 'unique_js_unmutated', 'unique_vs',
                  'unique_vs_unmutated', 'unmutated_frequency', 'unmutated_sequences',
                  'unmutated_umis', 'warnings'))


# Timestamped progress report

report = function(x) {
  cat(paste0(Sys.time(), ':', x, '\n'))
}

# Warning messages (echoed to pdf)

pkg.globals = new.env()
pkg.globals$report_warnings = ""

report_warn = function(x) {
  cat(x)
  pkg.globals$report_warnings = paste(pkg.globals$report_warnings, x, sep='\n', collapse='\n')
}

#' Generate OGRDB reports from specified files
#' @param ref_filename Name of file containing IMGT-aligned reference genes in FASTA format
#' @param inferred_filename Name of file containing sequences of inferred novel alleles, or '-' if none
#' @param species Species name used in field 3 of the IMGT germline header with spaces omitted, if the reference file is from IMGT. Otherwise ''
#' @param filename Name of file containing annotated reads in AIRR, CHANGEO or IgDiscover format. The format is detected automatically
#' @param chain one of VH, VK, VL, D, JH, JK, JL
#' @param hap_gene The haplotyping columns will be completed based on the usage of the two most frequent alleles of this gene. If NA, the column will be blank
#' @param segment one of V, D, J
#' @param chain_type one of H, L
#' @param plot_unmutated Plot base composition using only unmutated sequences (V-chains only)
#' @param all_inferred Treat all alleles as novel
#' @return Nothing
#' @export
generate_ogrdb_report = function(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, plot_unmutated, all_inferred=F) {
  report('Processing started')
  pdf(NULL) # this seems to stop an empty Rplots.pdf from being created. I don't know why.

  rd = read_input_files(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, all_inferred)

  file_prefix = basename(strsplit(filename, '.', fixed=T)[[1]][1])

  report('writing genotype file')
  write_genotype_file(paste0(file_prefix, '_ogrdb_report.csv'), segment, chain_type, rd$genotype)

  report('plotting bar charts')
  barplot_grobs = make_barplot_grobs(rd$input_sequences, rd$genotype_db, rd$inferred_seqs, rd$genotype, segment, rd$calculated_NC)

  report('plotting novel base composition')
  if(segment == 'V' && plot_unmutated) {
    nbgrobs = make_novel_base_grobs(rd$inferred_seqs, rd$input_sequences[rd$input_sequences$SEG_MUT_NC==0,], segment, all_inferred)
  } else {
    nbgrobs = make_novel_base_grobs(rd$inferred_seqs, rd$input_sequences, segment, all_inferred)
  }

  report('plotting haplotyping charts')
  haplo_grobs = make_haplo_grobs(segment, rd$haplo_details)

  report('writing plot file')

  if(plot_unmutated) {
    report_warn("Base plots are calculated from unmutated sequences only.")
  }

  write_plot_file(paste0(file_prefix, '_ogrdb_plots.pdf'), rd$input_sequences, nbgrobs$cdr3_dist, nbgrobs$end, nbgrobs$whole, nbgrobs$triplet, barplot_grobs, haplo_grobs$aplot, haplo_grobs$haplo, pkg.globals$report_warnings)
}

#' Read input files into memory
#' @param ref_filename Name of file containing IMGT-aligned reference genes in FASTA format
#' @param inferred_filename Name of file containing sequences of inferred novel alleles, or '-' if none
#' @param species Species name used in field 3 of the IMGT germline header with spaces omitted, if the reference file is from IMGT. Otherwise ''
#' @param filename Name of file containing annotated reads in AIRR, CHANGEO or IgDiscover format. The format is detected automatically
#' @param chain one of VH, VK, VL, D, JH, JK, JL
#' @param hap_gene The haplotyping columns will be completed based on the usage of the two most frequent alleles of this gene. If NA, the column will be blank
#' @param segment one of V, D, J
#' @param chain_type one of H, L
#' @param all_inferred Treat all alleles as novel
#' @return A named list containing the following elements:
#' \tabular{ll}{
#' ref_genes  \tab named list of IMGT-gapped reference genes \cr
#' inferred_seqs \tab named list of IMGT-gapped inferred (novel) sequences. \cr
#' input_sequences \tab data frame with one row per annotated read, with CHANGEO-style column names
#'                           One key point: the column SEQ_CALL is the gene call for the segment under analysis. Hence if segment is 'V',
#'                           'V_CALL' will be renamed 'SEG_CALL' whereas is segment is 'J', 'J_CALL' is renamed 'SEG_CALL'. This simplifies
#'                           downstream processing.
#'                           Rows in the input file with ambiguous SEG_CALLs, or no call, are removed. \cr
#' genotype_db \tab named list of gene sequences referenced in the annotated reads (both reference and novel sequences) \cr
#' haplo_details \tab data used for haplotype analysis, showing allelic ratios calculated with various potential haplotyping genes \cr
#' genotype \tab data frame containing information provided in the OGRDB genotype csv file \cr
#' calculated_NC \tab a boolean that is TRUE if mutation counts were calculated by this library, FALSE if they were read from the annotated read file \cr
#' }
#' @export
read_input_files = function(ref_filename, inferred_filename, species, filename, chain, hap_gene, segment, chain_type, all_inferred) {
  report('reading reference genes')
  ref_genes = read_reference_genes(ref_filename, species, chain, segment)

  report('reading inferred sequences')
  inferred_seqs = read_inferred_sequences(inferred_filename, segment, ref_genes)

  if (all_inferred) {
    report('Treating all sequences as inferred for the purpose of reporting')
    inferred_seqs = c(inferred_seqs, ref_genes[!(names(ref_genes) %in% names(inferred_seqs))])
  }

  report('reading input sequences')
  input_sequences = read_input_sequences(filename, segment, chain_type)
  genotype_db = make_genotype_db(input_sequences, inferred_seqs, ref_genes)

  report('checking and fixing read alignment')
  input_sequences = gap_input_sequences(input_sequences, inferred_seqs, ref_genes)

  report('calculating haplotype details')
  haplo_details = calc_haplo_details(segment, input_sequences)

  report('calculating genotype')
  c = calc_genotype(segment, chain_type, input_sequences, ref_genes, inferred_seqs, genotype_db, hap_gene, haplo_details)

  return(list('ref_genes'=ref_genes, 'inferred_seqs'=inferred_seqs, 'input_sequences'=c$input_sequences, 'genotype_db'=genotype_db,
           'haplo_details'=haplo_details, 'genotype'=c$genotype, 'calculated_NC'=c$calculated_NC))
}


#' Read the reference set from the specified file
#' @param ref_filename Name of file containing IMGT-aligned reference genes in FASTA format
#' @param species Species name used in field 3 of the IMGT germline header with spaces omitted, if the reference file is from IMGT. Otherwise ''
#' @param chain one of VH, VK, VL, D, JH, JK, JL
#' @param segment one of V, D, J
#' @return Named list of gene sequences
read_reference_genes = function(ref_filename, species, chain, segment) {
  # get the reference set

  ref_genes = readIgFasta(ref_filename, strip_down_name =F)
  set = paste0('IG', substr(chain, 2, 2), segment)
  region = paste0(segment, '-REGION')

  # process IMGT library, if header is in corresponding format
  if(grepl('|', names(ref_genes)[1], fixed=T)) {
    ref_genes = ref_genes[grepl(species, names(ref_genes),fixed=T)]
    ref_genes = ref_genes[grepl(region, names(ref_genes),fixed=T)]  # restrict to V-D-J regions
    ref_genes = ref_genes[grepl(set, names(ref_genes),fixed=T)]

    gene_name = function(full_name) {
      return(strsplit(full_name, '|', fixed=T)[[1]][2])
    }

    names(ref_genes) = sapply(names(ref_genes), gene_name)
  } else {
    ref_genes = ref_genes[grepl(set, names(ref_genes),fixed=T)]
  }

  # Crude fix for two misaligned human IGHJ reference genes

  if(set == 'IGHJ') {
    for(g in c('IGHJ6*02', 'IGHJ6*03')) {
      if(g %in% names(ref_genes) && str_sub(ref_genes[g], start= -1) == 'A') {
        ref_genes[g] = paste0(ref_genes[g], 'G')
        print(paste0('Modified truncated reference gene ', g, ' to ', ref_genes[g]))
      }
    }
  }

  # Check that the reference set is IMGT-aligned by looking for dots in the V-gene sequences

  misaligned_v = function(ref, name) {
    if(grepl('V', name, fixed=T)) {
      return(!grepl('.', ref, fixed=T))
    } else {
      return(NA)
    }
  }

  misaligned = mapply(misaligned_v, ref_genes, names(ref_genes))

  if(any(misaligned, na.rm=T)) {
    cat(paste0('Error: the following reference V-genes are not IMGT aligned:\n'))
    print(names(ref_genes)[misaligned])
    stop('The reference genes must be IMGT gap-aligned. Please consult the Prerequisites section in the README file for details.\n')
  }

  if(length(ref_genes) < 1) {
    report_warn('Warning: no reference genes were found.')
  }

  return(ref_genes)
}

#' Read inferred sequences from a file. If ungapped, they will be gapped using the closest available reference gene
#' @param inferred_filename Name of file containing sequences of inferred novel alleles, or '-' if none
#' @param segment one of V, D, J
#' @param ref_genes Named list of reference sequences
#' @return Named list of gapped inferred sequences
read_inferred_sequences = function(inferred_filename, segment, ref_genes) {
  # get the genotype and novel alleles in this set

  if(inferred_filename != '-') {
    inferred_seqs = readIgFasta(inferred_filename, strip_down_name=F)
  } else {
    inferred_seqs = c()
  }

  # ignore inferred genes that are listed in the reference. gap any that aren't already gapped.

  inferred_seqs = inferred_seqs[!(names(inferred_seqs) %in% names(ref_genes))]

  if(segment == 'V') {
    inferred_seqs = sapply(names(inferred_seqs), imgt_gap_inferred, seqs=inferred_seqs, ref_genes=ref_genes)
  }

  return(inferred_seqs)
}



#' Read annotated input sequences and convert format to CHANGEO-style, filling in as many gaps in what we get from the inference tool as we can
#' @param filename Name of file containing annotated reads in AIRR, CHANGEO or IgDiscover format. The format is detected automatically
#' @param segment one of V, D, J
#' @param chain_type one of H, L
#' @return Data Frame containing sequence annotations, with CHANGEO format headers.
read_input_sequences = function(filename, segment, chain_type) {
  # Read the sequences. Changeo format is assumed unless airr or IgDiscover format is identified
  # TODO - check and give nice error message if any columns are missing

  s = read.delim(filename, stringsAsFactors=F)

  # Fallback in airr,changeo if no v_call_genotyped fied is present

  if(!('V_CALL_GENOTYPED' %in% names(s)) && !('v_call_genotyped' %in% names(s)) && ( ('V_CALL' %in% names(s)) || ('v_call' %in% names(s)))) {
    report_warn('No v_call_genotyped field found. Falling back to v_call (please see notes on genotyping).')
    if('V_CALL' %in% names(s))
    {
      s = rename(s, V_CALL_GENOTYPED=V_CALL)
    } else {
      s = rename(s, v_call_genotyped=v_call)
    }
  }

  if('sequence_id' %in% names(s))
  {
    # airr format - convert wanted fields to changeo
    # TODO - would be better to get the field ranges from the CIGAR string, given that it's a required field
    # https://www.bioconductor.org/packages/devel/bioc/manuals/GenomicAlignments/man/GenomicAlignments.pdf

    # add a dummy D_CALL to light chain repertoires, for ease of processing

    if(chain_type == 'L' && !('d_call' %in% names(s))) {
      s$d_call = ''
    }

    req_names = c('sequence', 'sequence_id', 'v_call_genotyped', 'd_call', 'j_call', 'sequence_alignment', 'cdr3')
    col_names = c('SEQUENCE_INPUT', 'SEQUENCE_ID', 'V_CALL_GENOTYPED', 'D_CALL', 'J_CALL', 'SEQUENCE_IMGT', 'CDR3_IMGT')

    if(segment != 'V') {
      a_seg = tolower(segment)
      req_names = c(req_names, paste0(a_seg, '_sequence_start'), paste0(a_seg, '_sequence_end'), paste0(a_seg, '_germline_start'), paste0(a_seg, '_germline_end'))
      col_names = c(col_names, 'SEG_SEQ_START', 'SEG_SEQ_END', 'SEG_GERM_START', 'SEG_GERM_END')
    }

    s = select(s, req_names)
    names(s) = col_names

    if(segment != 'V') {
      s$SEG_SEQ_LENGTH = s$SEG_SEQ_END - s$SEG_SEQ_START + 1
      s$SEG_GERM_LENGTH = s$SEG_GERM_END - s$SEG_GERM_START + 1
    }

  } else if('V_gene' %in% names(s)) {
    # IgDiscover format
    #  s = uncount(s, count)  for consistency with IgDiscover results, count each unique record in the file only once, regardless of 'count'

    # add a dummy D_CALL to light chane repertoires, for ease of processing

    if(chain_type == 'L' && !('D_gene' %in% names(s))) {
      s$D_gene = ''
      s$D_errors = ''
      s$D_region = ''
    }

    if(chain_type == 'L' && !('VDJ_nt' %in% names(s))) {
      s$VDJ_nt = s$VJ_nt
    }

    col_names = c('SEQUENCE_ID', 'V_CALL_GENOTYPED', 'D_CALL', 'J_CALL', 'CDR3_IMGT', 'V_MUT_NC', 'D_MUT_NC', 'J_MUT_NC', 'SEQUENCE', 'JUNCTION_START', 'V_SEQ', 'D_SEQ', 'J_SEQ')
    s = select(s, name, V_gene, D_gene, J_gene, CDR3_nt, V_errors, D_errors, J_errors, VDJ_nt, V_CDR3_start, V_nt, D_region, J_nt)
    names(s) = col_names
    # adjust IgDiscover's V_CDR3_START to the 1- based location of the conserved cysteine
    s$JUNCTION_START = s$JUNCTION_START - 2
    s$SEQUENCE = toupper(s$SEQUENCE)

    if(segment == 'V') {
      s = rename(s, SEG_MUT_NC=V_MUT_NC, SEG_SEQ=V_SEQ)
    } else if(segment == 'D') {
      s = rename(s, SEG_MUT_NC=D_MUT_NC, SEG_SEQ=D_SEQ)
    } else {
      s = rename(s, SEG_MUT_NC=J_MUT_NC, SEG_SEQ=J_SEQ)
    }
  }

  if(segment == 'V') {
    s = rename(s, SEG_CALL=V_CALL_GENOTYPED)
  } else if(segment == 'D') {
    s = rename(s, SEG_CALL=D_CALL)
  } else {
    s = rename(s, SEG_CALL=J_CALL)
  }

  if('V_CALL_GENOTYPED' %in% names(s)) {
    if('V_CALL' %in% names(s)) {
      s = subset(s, select = -c('V_CALL'))
    }
    s = rename(s, V_CALL=V_CALL_GENOTYPED)
  }

  s$CDR3_IMGT = toupper(s$CDR3_IMGT)

  #remove sequences with ambiguous calls, or no call, in the target segment

  s = s[!grepl(',', s$SEG_CALL),]
  s = s[!(s$SEG_CALL == ''),]

  # At this point, s contains Changeo-named columns with all data required for calculations

  return(s)
}

#' Build a list of gene sequences that matches the personalised genotype
#' @param input_sequences the input_sequences data frame
#' @param inferred_seqs named list of novel gene sequences
#' @param ref_genes named list of reference genes
#' @return named list of gene sequences
make_genotype_db = function(input_sequences, inferred_seqs, ref_genes) {
  # Warn if we don't have genotype statistics for any of the inferred alleles
  # this can happen, for example, with Tigger, if novel alleles are detected but do not pass subsequent criteria for being included in the genotype.

  genotype_alleles = unique(input_sequences$SEG_CALL)
  missing = inferred_seqs[!(names(inferred_seqs) %in% genotype_alleles)]

  if(length(missing) >= 1) {
    report_warn(paste('Novel sequence(s)', paste0(names(missing), collapse=' '), 'are not listed in the genotype and will be ignored.', sep=' ', '\n'))
    inferred_seqs = inferred_seqs[(names(inferred_seqs) %in% genotype_alleles)]
  }

  genotype_alleles = genotype_alleles[!(genotype_alleles %in% names(inferred_seqs))]
  genotype_seqs = lapply(genotype_alleles, function(x) {ref_genes[x]})
  genotype_db = setNames(c(genotype_seqs, inferred_seqs), c(genotype_alleles, names(inferred_seqs)))


  # Check we have sequences for all alleles named in the reads - either from the reference set or from the inferred sequences
  # otherwise - one of these two is incomplete!

  if(any(is.na(genotype_db))) {
    report_warn(paste0("Warning: sequence(s) for allele(s) ", names(genotype_db[is.na(genotype_db)]), " can't be found in the reference set or the novel alleles file.\n"))
  }

  # CHeck that a reasonable number of genes in the reference set are not called in the sequence set

  unused_ref = length(setdiff(names(ref_genes), names(genotype_db)))

  if(unused_ref < (length(names(ref_genes))/10)) {
    report_warn("Warning: Over 90% of reference genes are called in the repertoire. This suggests that either the reference set is incomplete, or the genotype has not been personalised.")
  }

  return(genotype_db)
}

#' Check whether input sequences are gapped, and gap them if necessary
#' @param input_sequences the input_sequences data frame
#' @param inferred_seqs named list of novel gene sequences
#' @param ref_genes named list of reference genes
#' @return a revised input_sequences structure, guaranteed to contain gapped sequences in SEQUENCE_IMGT
gap_input_sequences = function(input_sequences, inferred_seqs, ref_genes) {
  find_template = function(call) {
    tem = ref_genes[call]
    if(is.na(tem))
      tem = inferred_seqs[call]

    return(tem)
  }

  input_sequences$SEG_REF_IMGT = sapply(input_sequences$SEG_CALL, find_template)

  if(!('SEQUENCE_IMGT' %in% names(input_sequences))) {
    input_sequences$SEQUENCE_IMGT = mapply(imgt_gap, input_sequences$SEQUENCE,input_sequences$CDR3_IMGT, input_sequences$JUNCTION_START, input_sequences$SEG_REF_IMGT)
  } else {
    # remove any sequences that do not have an aligned sequence

    input_sequences$SEQ_LEN=str_length(input_sequences$SEQUENCE_IMGT)
    count_zero = length(input_sequences$SEQ_LEN[input_sequences$SEQ_LEN==0])

    if(count_zero > 0) {
      print(paste0('Warning: removing ', count_zero, ' sequences with no SEQUENCE_IMGT'))
      input_sequences = input_sequences[input_sequences$SEQ_LEN>0,]
    }
  }

  input_sequences$SEQUENCE_IMGT = toupper(input_sequences$SEQUENCE_IMGT)

  return(input_sequences)
}

#' Calculate the data used for haplotype analysis, showing allelic ratios calculated with various potential haplotyping genes
#' @param segment one of V, D, J
#' @param input_sequences the input_sequences data frame
#' @return A named list containing the following elements:
#' \tabular{ll}{
#'     sa \tab           a variant of the input_sequences data frame, with fields for haplotyping analysis \cr
#'     a_genes \tab       list of potential haplotyping genes (if the segment under analysis is J, these are V genes, otherwise J genes) \cr
#'     a_props \tab       proportion of each allele of each potential haplotyping gene in the input sequences\cr
#'}
calc_haplo_details = function(segment, input_sequences) {
  if(segment == 'V' || segment == 'D') {
    sa = input_sequences[!grepl(',', input_sequences$J_CALL),]            # unique J-calls only
    sa = rename(sa, A_CALL=J_CALL)
  } else {
    sa = input_sequences[!grepl(',', input_sequences$V_CALL),]
    sa = rename(sa, A_CALL=V_CALL)
  }

  sa$a_gene = sapply(sa$A_CALL, function(x) {strsplit(x, '*', fixed=T)[[1]][[1]]})
  sa$a_allele = sapply(sa$A_CALL, function(x) {
      if (grepl('*', x, fixed = T)) {
        strsplit(x, '*', fixed=T)[[1]][[2]]
      } else {
        'X'
      }
    })
  sa$a_gene = factor(sa$a_gene, sort_alleles(unique(sa$a_gene)))

  su = select(sa, A_CALL, a_gene, a_allele)
  a_genes = sort(unlist(unique(su$a_gene)))

  sa = sa[!grepl(',', sa$SEG_CALL),]        # remove ambiguous V-calls
  sa$SEG_CALL = factor(sa$SEG_CALL, sort_alleles(unique(sa$SEG_CALL)))

  # calc percentage of each allele in a gene
  allele_props = function(gene, su) {
    alleles = su %>% filter(a_gene==gene) %>% group_by(a_allele) %>% summarise(count=n())
    alleles$a_gene = gene
    total = sum(alleles$count)
    alleles$percent = 100*alleles$count/total
    return(alleles)
  }

  a_props = do.call('rbind', lapply(a_genes, allele_props, su=su))

  return(list('sa'=sa, 'a_props'=a_props, 'a_genes'=a_genes))
}


#' Build the genotype data required in the OGRDB genotype file
#' @param segment one of V, D, J
#' @param chain_type one of H, L
#' @param s the input_sequences data frame
#' @param ref_genes named list of reference genes
#' @param inferred_seqs named list of novel gene sequences
#' @param genotype_db named list of gene sequences in the personalised genotype
#' @param hap_gene The haplotyping columns will be completed based on the usage of the two most frequent alleles of this gene. If NA, the column will be blank
#' @param haplo_details Data structure created by create_haplo_details
#' @return A named list containing the following elements:
#' \tabular{ll}{
#'     input_sequences \tab     updated data frame, guaranteed to include mutation counts, which have been calculated if necessary \cr
#'     calculated_NC  \tab      a boolean, TRUE if mutation counts had to be calculated, FALSE otherwise \cr
#'     genotype \tab            data frame containing the information required for the OGRDB genotype file \cr
#'}
calc_genotype = function(segment, chain_type, s, ref_genes, inferred_seqs, genotype_db, hap_gene, haplo_details) {

  # unmutated count for each allele

  calculated_NC = F

  if(!('SEG_MUT_NC' %in% names(s))) {
    if(segment == 'V') {
      # We take the count up to the 2nd CYS at 310
      # This matches IgDiscover practice and facilitates Tigger's reassignAlles approach which does not re-analyse the junction with the novel V-allele
      s$SEQUENCE_IMGT_TRUNC = sapply(s$SEQUENCE_IMGT, substring, first=1, last=309)
      s$SEG_MUT_NC = unlist(getMutCount(s$SEQUENCE_IMGT_TRUNC, s$SEG_CALL, genotype_db))
    } else {
      s$SEG_SEQ = mapply(substr, s$SEQUENCE_INPUT, s$SEG_SEQ_START, s$SEG_SEQ_START+s$SEG_SEQ_LENGTH-1)
      s$SEG_REF_SEQ = mapply(substr, s$SEG_REF_IMGT, s$SEG_GERM_START, s$SEG_GERM_START+s$SEG_GERM_LENGTH-1)
      s$SEG_MUT_NC = stringdist(s$SEG_SEQ, s$SEG_REF_SEQ, method="hamming")
    }

    calculated_NC = T
    s[,is.na(s$SEG_MUT_NC)]$SEG_MUT_NC = 0
  }

  report('calculated mutation count')

  # make a space-alignment of D sequences for the usage histogram (V and J use the IMGT alignment)

  if(segment == 'D') {
    s$SEG_SEQ_ALIGNED = mapply(paste0, sapply(s$SEG_GERM_START, function(x) {paste(rep(' ', x), collapse='')}), s$SEG_SEQ)
    width = max(str_length(s$SEG_SEQ_ALIGNED))
    s$SEG_SEQ_ALIGNED = sapply(s$SEG_SEQ_ALIGNED, function(x) {str_pad(x, width, side='right')})
  }

  genotype = s %>% group_by(SEG_CALL) %>% summarize(sequences = n())
  s0 = s[s$SEG_MUT_NC == 0,] %>% group_by(SEG_CALL) %>% summarize(unmutated_sequences = n())
  genotype = merge(genotype, s0, all=T)

  if(any(is.na(genotype$unmutated_sequences))) {
    genotype[is.na(genotype$unmutated_sequences),]$unmutated_sequences = 0
  }

  total_unmutated = sum(genotype$unmutated_sequences)
  genotype$unmutated_frequency = round(100*genotype$unmutated_sequences/total_unmutated, digits=2)

  s_totals = s %>% group_by(SEG_CALL) %>% summarize(sequences = n())
  genotype = merge(genotype, s_totals, all=T)

  genotype$GENE = sapply(genotype$SEG_CALL, function(x) {if(grepl('*', x, fixed=T)) {strsplit(x, '*', fixed=T)[[1]][1]} else {x}})
  allelic_totals = genotype %>% group_by(GENE) %>% summarise(allelic_total=sum(sequences))
  genotype = merge(genotype, allelic_totals, all=T)
  genotype$allelic_percentage = round(100*genotype$sequences/genotype$allelic_total)

  su=s[s$SEG_MUT_NC==0,]

  if(segment != 'V') {
    genotype$unique_vs = sapply(genotype$SEG_CALL, unique_calls, segment='V_CALL', seqs=s)
    genotype$unique_vs_unmutated = sapply(genotype$SEG_CALL, unique_calls, segment='V_CALL', seqs=su)
  }

  if(segment != 'D' && chain_type == 'H') {
    genotype$unique_ds = sapply(genotype$SEG_CALL, unique_calls, segment='D_CALL', seqs=s)
    genotype$unique_ds_unmutated = sapply(genotype$SEG_CALL, unique_calls, segment='D_CALL', seqs=su)
  }

  if(segment != 'J') {
    genotype$unique_js = sapply(genotype$SEG_CALL, unique_calls, segment='J_CALL', seqs=s)
    genotype$unique_js_unmutated = sapply(genotype$SEG_CALL, unique_calls, segment='J_CALL', seqs=su)
  }

  genotype$unique_cdr3s = sapply(genotype$SEG_CALL, unique_cdrs, segment='CDR3_IMGT', seqs=s)
  genotype$unique_cdr3s_unmutated = sapply(genotype$SEG_CALL, unique_cdrs, segment='CDR3_IMGT', seqs=su)

  genotype$assigned_unmutated_frequency = round(100*genotype$unmutated_sequences/genotype$sequences, digits=2)

  # closest in genotype and in reference (inferred alleles only)
  # Inferred D alleles should be aligned for best match (if this is an allele of an existing D-gene, align against a knon allele of that gene)

  if (length(inferred_seqs) == 0) {
    report_warn('Warning - no inferred sequences found.\n')

    genotype$reference_closest = NA
    genotype$host_closest = NA
    genotype$reference_difference = NA
    genotype$reference_nt_diffs = NA
    genotype$reference_aa_difference = NA
    genotype$reference_aa_subs = NA
    genotype$host_difference = NA
    genotype$host_nt_diffs = NA
    genotype$host_aa_difference = NA
    genotype$host_aa_subs = NA
  } else {
    nearest_ref = data.frame(t(sapply(seq_along(inferred_seqs), find_nearest, ref_genes=ref_genes, prefix='reference', inferred_seqs=inferred_seqs, segment=segment)))
    nearest_ref$SEG_CALL = names(inferred_seqs)
    genotype = merge(genotype, nearest_ref, by='SEG_CALL', all.x=T)

    nearest_ref = data.frame(t(sapply(seq_along(inferred_seqs), find_nearest, ref_genes=ref_genes[names(ref_genes) %in% genotype$SEG_CALL], prefix='host', inferred_seqs=inferred_seqs, segment=segment)))
    nearest_ref$SEG_CALL = names(inferred_seqs)
    genotype = merge(genotype, nearest_ref, by='SEG_CALL', all.x=T)
  }


  genotype$nt_sequence = sapply(genotype$SEG_CALL, function(x) genotype_db[x][[1]])

  genotype = rename(genotype, sequence_id=SEG_CALL, closest_reference=reference_closest, closest_host=host_closest,
                           nt_diff=reference_difference, nt_diff_host=host_difference, nt_substitutions=reference_nt_diffs, aa_diff=reference_aa_difference,
                           aa_substitutions=reference_aa_subs)

  genotype$unmutated_umis = ''
  genotype$nt_sequence_gapped = genotype$nt_sequence
  genotype$nt_sequence = gsub('-', '', genotype$nt_sequence, fixed=T)
  genotype$nt_sequence = gsub('.', '', genotype$nt_sequence, fixed=T)
  genotype = unnest(genotype, cols = c(closest_reference, nt_diff, nt_substitutions, aa_diff, aa_substitutions,
                                       closest_host, nt_diff_host, host_nt_diffs, host_aa_difference,
                                       host_aa_subs))

  # If we had to calculate MUT_NC, set the mutated counts to NA for any allele for which we don't have a sequence

  if(calculated_NC && any(is.na(genotype$nt_sequence))) {
    genotype[is.na(genotype$nt_sequence),]$unmutated_frequency = NA
    genotype[is.na(genotype$nt_sequence),]$assigned_unmutated_frequency = NA
    genotype[is.na(genotype$nt_sequence),]$unmutated_sequences = NA
  }

  # Check for duplicate germline sequences

  concat_names = function(x) {
    paste(x, collapse=', ')
  }

  warn_dupes = function(x) {
    report_warn(paste0('Warning: ', x, ' have identical germline sequences.\n'))
  }

  dupes = aggregate(genotype["sequence_id"], by=genotype["nt_sequence"], FUN=concat_names)
  dupes = dupes$sequence_id[grepl(',', dupes$sequence_id, fixed=T)]

  if(length(dupes) > 0) {
    x = lapply(dupes, warn_dupes)
  }

  # Add haplo details
  genotype = add_hap_stats(genotype, hap_gene, haplo_details)

  return(list("input_sequences"=s, "genotype"=genotype, "calculated_NC"=calculated_NC))
}


#' If the haplotyping gene has been specified, add haplotyping ratios to the genotype data frame, provided the ratios for that gene are suitable
#' @param genotype genotype data frame
#' @param hap_gene haplotyping gene for which ratios should be calculated
#' @param haplo_details Data structure created by create_haplo_details
#' @return modified genotype with haplotyping ratios added
add_hap_stats = function(genotype, hap_gene, haplo_details) {
  sa = haplo_details$sa
  a_props = haplo_details$a_props

  genotype$haplotyping_gene = ''
  genotype$haplotyping_ratio = ''

  if(!is.na(hap_gene)) {
    ap = a_props[a_props$a_gene==hap_gene,]
    ap = ap[order(ap$percent, decreasing=T),]

    if(nrow(ap) < 2 || ap[1,]$percent > 75 || ap[2,]$percent < 20 || (ap[1,]$percent+ap[2,]$percent < 90)) {
      cat(paste0('Alelleic ratio is unsuitable for haplotyping analysis based on ', hap_gene, '\n'))
    } else
    {
      genotype$haplotyping_gene = hap_gene
      genotype = select(genotype, -c(haplotyping_ratio))

      a1 = ap[1,]$a_allele
      a2 = ap[2,]$a_allele
      cat(paste0('Haplotyping analysis is based on gene ', hap_gene, ' alleles ', a1, ':', a2, '\n'))
      recs = sa %>% filter(a_gene==hap_gene) %>% filter(a_allele==a1 | a_allele==a2)
      recs = recs %>% select(sequence_id=SEG_CALL, a_allele) %>% group_by(sequence_id, a_allele) %>% summarise(count=n()) %>% spread(a_allele, count)
      recs[is.na(recs)] = 0
      names(recs) = c('sequence_id', 'a1', 'a2')
      recs$totals = recs$a1 + recs$a2
      recs$a1 = round(100 * recs$a1/ recs$totals, 0)
      recs$a2 = round(100 * recs$a2/ recs$totals, 0)
      recs$haplotyping_ratio = paste0(' ', recs$a1, ':', recs$a2, ' ')
      recs = recs %>% select(sequence_id, haplotyping_ratio)
      genotype = merge(genotype, recs, all=T)
    }
  }

  return(genotype)
}

#' Write the genotype file required by OGRDB
#' @export
#' @param filename name of file to create (csv)
#' @param segment one of V, D, J
#' @param chain_type one of H, L
#' @param genotype genotype data frame
#' @return nothing
write_genotype_file = function(filename, segment, chain_type, genotype) {
  genotype = genotype[order_alleles(data.frame(genes=as.character(genotype$sequence_id), stringsAsFactors = F)),]

  if(chain_type == 'H') {
    if(segment == 'V') {
      g = select(genotype, sequence_id, sequences, closest_reference, closest_host, nt_diff, nt_diff_host, nt_substitutions, aa_diff,
               aa_substitutions, assigned_unmutated_frequency, unmutated_frequency, unmutated_sequences, unmutated_umis, allelic_percentage, unique_ds,
               unique_js,unique_cdr3s, unique_ds_unmutated, unique_js_unmutated, unique_cdr3s_unmutated, haplotyping_gene, haplotyping_ratio, nt_sequence, nt_sequence_gapped)
    } else if(segment == 'D') {
      g = select(genotype, sequence_id, sequences, closest_reference, closest_host, nt_diff, nt_diff_host, nt_substitutions, aa_diff,
                 aa_substitutions, assigned_unmutated_frequency, unmutated_frequency, unmutated_sequences, unmutated_umis, allelic_percentage, unique_vs,
                 unique_js,unique_cdr3s, unique_vs_unmutated, unique_js_unmutated, unique_cdr3s_unmutated, haplotyping_gene, haplotyping_ratio, nt_sequence)
    } else if(segment == 'J') {
      g = select(genotype, sequence_id, sequences, closest_reference, closest_host, nt_diff, nt_diff_host, nt_substitutions, aa_diff,
                 aa_substitutions, assigned_unmutated_frequency, unmutated_frequency, unmutated_sequences, unmutated_umis, allelic_percentage, unique_vs,
                 unique_ds,unique_cdr3s, unique_vs_unmutated, unique_ds_unmutated, unique_cdr3s_unmutated, haplotyping_gene, haplotyping_ratio, nt_sequence)
    }
  } else {
    if(segment == 'V') {
      g = select(genotype, sequence_id, sequences, closest_reference, closest_host, nt_diff, nt_diff_host, nt_substitutions, aa_diff,
                 aa_substitutions, assigned_unmutated_frequency, unmutated_frequency, unmutated_sequences, unmutated_umis, allelic_percentage,
                 unique_js,unique_cdr3s, unique_js_unmutated, unique_cdr3s_unmutated, haplotyping_gene, haplotyping_ratio, nt_sequence, nt_sequence_gapped)
    } else if(segment == 'J') {
      g = select(genotype, sequence_id, sequences, closest_reference, closest_host, nt_diff, nt_diff_host, nt_substitutions, aa_diff,
                 aa_substitutions, assigned_unmutated_frequency, unmutated_frequency, unmutated_sequences, unmutated_umis, allelic_percentage, unique_vs,
                 unique_cdr3s, unique_vs_unmutated, unique_cdr3s_unmutated, haplotyping_gene, haplotyping_ratio, nt_sequence)
    }
  }

  g[] = lapply(g, as.character)
  g[is.na(g)] = ''

  write.csv(g, filename, row.names=F)
}


