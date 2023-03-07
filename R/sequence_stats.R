# This software is licensed under the CC BY-SA 4.0 licence: https://creativecommons.org/licenses/by-sa/4.0/
#
# Some functions are adapted from TIgGER (https://tigger.readthedocs.io) with thanks to the authors.



# count unique calls
unique_calls = function(gene, segment, seqs) {
  calls = unique(seqs[seqs$SEG_CALL==gene,][segment])
  calls = calls[!grepl(',', calls[,segment]),]       # count unambiguous calls only
  calls = calls[grepl(substr(segment,1,1), calls)]                 #don't count blank calls
  return(length(calls))
}

# count unique CDRs
unique_cdrs = function(gene, segment, seqs) {
  calls = unique(seqs[seqs$SEG_CALL==gene,][segment])
  calls = calls[nchar(calls[,segment]) > 0,]
  return(length(calls))
}

# create list of nt differences between two strings
build_nt_diff_string = function(seq1, seq2, bias) {
  nt_diff_string = ""
  max_length = max(nchar(seq1), nchar(seq2))
  seq1 = stringr::str_pad(seq1, max_length, side="right", pad="-")
  seq2 = stringr::str_pad(seq2, max_length, side="right", pad="-")
  nt_diff = unlist(tigger::getMutatedPositions(seq1, seq2,ignored_regex="\\."))
  ref_nt <- strsplit(seq1,"")[[1]][nt_diff]
  novel_nt <- strsplit(seq2,"")[[1]][nt_diff]

  if(length(nt_diff) > 0) {
    nt_diff = sapply((nt_diff-bias), function(x) if(x <= 0) {x-1} else {x} )
    nt_diff_string <- paste(paste(
      ref_nt,
      nt_diff,
      replace(novel_nt, is.na(novel_nt), "-"),
      sep=""), collapse=",")
  }
  return(nt_diff_string)
}

# Drop any stray nucleotides not on a codon boundary

fixNonImgtGaps <- function (seq) {
  len <- ceiling(nchar(seq)/3)*3
  codons <- substring(seq, seq(1, len-2, 3), seq(3, len, 3))
  gaps_lengths <- nchar(gsub("[^\\.\\-]", "", codons))
  if (any(gaps_lengths %% 3 != 0)) {
    codons[gaps_lengths %% 3 != 0] = "..."
    seq = paste(codons, collapse='')
  }

  seq
}

# Find non triplet gaps in a nucleotide sequence

hasNonImgtGaps <- function (seq) {
  len <- ceiling(nchar(seq)/3)*3
  codons <- substring(seq, seq(1, len-2, 3), seq(3, len, 3))
  gaps_lengths <- nchar(gsub("[^\\.\\-]", "", codons))
  if (any(gaps_lengths %% 3 != 0)) {
    TRUE
  } else {
    FALSE
  }
}


# Use a gapped IMGT 'template' to apply gaps to a similar sequence

insert_at = function(seq, ins, loc) {
  if(nchar(seq) >= loc-1) {
    return(paste0(substr(seq, 1, loc-1), ins, substr(seq, loc, nchar(seq))))
  }


  return(seq)
}

apply_gaps = function(seq, tem) {
  tem = strsplit(tem, '')[[1]]
  res = seq

  for(i in 1:length(tem)) {
    if(tem[i] == '.') {
      res = insert_at(res, '.', i)
    }
  }
  return(res)
}

# Gap an ungapped sequence, using the nominated reference v_gene as template
# This is used where the inference tool does not provide imgt-aligned sequences.
# It will not handle indels, but it will try to spot them and set the aligned sequence to NA
# junction_start is the location of the first nucleotide of the cysteine preceding the CRD3
# this is location 310 in the IMGT numbering scheme
imgt_gap = function(sequence, cdr3_sequence, junction_start, ref_gene) {
  if(is.na(ref_gene) || is.na(cdr3_sequence) || is.na(sequence)) {
    return(NA)
  }

  # Find the cdr3_start in the un-aligned reference gene
  ref_junction_start = 310 - (nchar(ref_gene) - nchar(gsub('.', '', ref_gene, fixed=TRUE)))

  # Trim or pad this sequence to match the unaligned ref gene
  if(junction_start < ref_junction_start) {
    pad = strrep('-', ref_junction_start - junction_start)
    sequence = paste0(pad, sequence)
  } else if(junction_start > ref_junction_start) {
    chop = junction_start - ref_junction_start
    sequence = substring(sequence, chop+1)
  }

  gapped = apply_gaps(sequence, ref_gene)

  # if the alignment is correct, the cdr3 should start at location 313
  # if it doesn't, the sequence probably has indels - set the alignment to NA
  # Ideally we would check CDR1 as well, but as partis doesn't give us CDR1...

  aligned_cdr3 = substring(gapped, 313, 313+nchar(cdr3_sequence)-1)

  if(cdr3_sequence != aligned_cdr3) {
    gapped = NA
  }

  return(gapped)
}

# Gap an ungapped inferred germline sequence, using the closest reference gene as a template
# This is used to IMGT-gap inferred alleles, where the gapped sequence is not available
# It assumes full-length sequences and no indels.

imgt_gap_inferred = function(seqname, seqs, ref_genes) {

  # Do we need to gap?
  if(grepl('.', seqs[seqname], fixed=TRUE))
    return(unname(seqs[seqname]))

  # Find the closest reference gene
  r = data.frame(GENE=names(ref_genes),SEQ=ref_genes, stringsAsFactors = FALSE)
  r$SEQ = sapply(r$SEQ,stringr::str_replace_all,pattern='\\.',replacement='')
  r$dist=sapply(r$SEQ, Biostrings::pairwiseAlignment, subject=seqs[seqname], scoreOnly=TRUE)
  r = r[order(r$dist, decreasing=TRUE),]

  # Gap the sequence
  tem = ref_genes[r[1,]$GENE]
  gapped = apply_gaps(seqs[seqname], tem)
  report_note(paste0('Inferred gene ', seqname, ' gapped using ', r[1,]$GENE, ': ', gapped,'\n'))
  return(gapped)
}


# Compare two IMGT gapped sequences and find AA mutations
getMutatedAA <- function(ref_imgt, novel_imgt, ref_name, seq_name, segment, bias) {
  if (grepl("N", ref_imgt)) {
    report_note(paste0("Unexpected N in reference sequence ", ref_name, ": replacing with gap\n"))
    stringr::str_replace(ref_imgt, "N", "-")
  }
  if (grepl("N", novel_imgt)) {
    report_warn(paste0("Unexpected N in novel sequence ", seq_name, ": replacing with gap\n"))
    stringr::str_replace(novel_imgt, "N", "-")
  }

  if (segment == 'V' && hasNonImgtGaps(ref_imgt)) {
    report_warn(paste0("warning: non codon-aligned gaps in reference sequence ", ref_name, "\n"))
    ref_imgt = fixNonImgtGaps(ref_imgt)
  }

  if (segment == 'V' && hasNonImgtGaps(novel_imgt)) {
    report_warn(paste0("warning: non codon-aligned gaps were found in novel sequence ", seq_name, "\n"))
    novel_imgt = fixNonImgtGaps(novel_imgt)
  }

  ref_imgt <- alakazam::translateDNA(ref_imgt)
  novel_imgt <- alakazam::translateDNA(novel_imgt)
  max_length = max(nchar(ref_imgt), nchar(novel_imgt))
  ref_imgt = stringr::str_pad(ref_imgt, max_length, side="right", pad="-")
  novel_imgt = stringr::str_pad(novel_imgt, max_length, side="right", pad="-")
  ref_imgt = strsplit(ref_imgt, "")[[1]]
  novel_imgt = strsplit(novel_imgt, "")[[1]]

  mutations <- c()
  diff_idx <- which(ref_imgt != novel_imgt)
  if (length(diff_idx)>0) {
    index = sapply((diff_idx-bias), function(x) if(x <= 0) {x-1} else {x} )
    mutations <- paste0(ref_imgt[diff_idx], index,
                        replace(novel_imgt[diff_idx], is.na(novel_imgt[diff_idx]),"-"))
  }
  mutations
}

# Find nearest reference sequences and enumerate differences
find_nearest = function(sequence_ind, ref_genes, prefix, inferred_seqs, segment) {

  if (length(ref_genes) == 0) {
    l = list(closest=' ', difference=0, nt_diffs=' ', aa_difference=' ', aa_subs=' ')
    names(l) = (paste(prefix, names(l), sep='_'))
    return(l)
  }

  sequence = inferred_seqs[[sequence_ind]]
  seq_name = names(inferred_seqs)[[sequence_ind]]
  r = data.frame(GENE=names(ref_genes),SEQ=ref_genes, stringsAsFactors = FALSE)
  r = r[r$GENE != seq_name, ]

  # pad all Js so that they align on the right

  if(segment == 'J') {
    w = as.integer(max(max(nchar(r$SEQ)), nchar(sequence))/3)*3
    r$SEQ = stringr::str_pad(r$SEQ, w, 'left', '-')
    sequence = stringr::str_pad(sequence, w, 'left', '-')
  }

  r$diff = tigger::getMutatedPositions(r$SEQ, sequence)
  r$num_diff = sapply(r$diff, length)
  r = r[order(r$num_diff),]
  r = r[1,]

  if(segment == 'J') {
    # 'bias' the index position so that the first nucleotide of the reference is 1
    bias = nchar(r$SEQ) - nchar(ref_genes[r$GENE])
  } else {
    bias = 0
  }

  nt_diff_string = build_nt_diff_string(r$SEQ, sequence, bias)

  diff_aa = getMutatedAA(r$SEQ, sequence, r$GENE, seq_name, segment, (bias/3))

  if(length(diff_aa) > 0) {
    aa_diffs = length(diff_aa)
    aa_subs = paste(diff_aa,collapse=",")
  } else {
    aa_diffs = 0
    aa_subs = ""
  }

  l = list(closest=r[1,]$GENE, difference=r[1,]$num_diff, nt_diffs=nt_diff_string, aa_difference=aa_diffs, aa_subs=aa_subs)
  names(l) = (paste(prefix, names(l), sep='_'))
  return(l)
}
