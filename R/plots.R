# This software is licensed under the CC BY-SA 4.0 licence: https://creativecommons.org/licenses/by-sa/4.0/

#' Create a barplot for each allele, showing number of reads distributed by mutation count
#' @export
#' @param input_sequences the input_sequences data frame
#' @param genotype_db named list of gene sequences in the personalised genotype
#' @param inferred_seqs named list of novel gene sequences
#' @param genotype data frame created by calc_genotype
#' @param segment one of V, D, J
#' @param calculated_NC a boolean, TRUE if mutation counts had to be calculated, FALSE otherwise
#' @return list of grobs
make_barplot_grobs =function(input_sequences, genotype_db, inferred_seqs, genotype, segment, calculated_NC) {
  # bar charts for all alleles

  # If we had to calculate MUT_NC, omit any plots for alleles for which we don't have a sequence

  if(calculated_NC) {
    barplot_grobs = lapply(sort_alleles(names(genotype_db[!is.na(genotype_db)])), plot_allele_seqs, s=input_sequences, inferred_seqs=inferred_seqs, genotype=genotype, segment=segment)
  } else {
    barplot_grobs = lapply(sort_alleles(names(genotype_db)), plot_allele_seqs, s=input_sequences, inferred_seqs=inferred_seqs, genotype=genotype, segment=segment)
  }

  barplot_grobs=barplot_grobs[!is.na(barplot_grobs)]

  return(barplot_grobs)
}

#' Create plots showing base usage at selected locations in sequences based on novel alleles
#' @export
#' @param inferred_seqs named list of novel gene sequences
#' @param input_sequences the input_sequences data frame
#' @param segment one of V, D, J
#' @param all_inferred true if user has requested all alleles in reference set plotted - will suppress some warnings
#' @return named list containing the following elements:
#' \tabular{ll}{
#'     cdr3_dist \tab cdr3 length distribution plots
#'     whole \tab     whole-length usage plots \cr
#'     end  \tab     3' end usage plots \cr
#'     triplet \tab   3' end triplet usage plots \cr
#'}
make_novel_base_grobs = function(inferred_seqs, input_sequences, segment, all_inferred) {

  if(length(inferred_seqs) < 1) {
    return(list('whole'=list(), 'end'=list(), 'triplet'=list()))
  }

  whole_composition_grobs = c()
  end_composition_grobs = c()
  triplet_grobs = c()
  composition_heatmaps = c()

  cdr3_distribution_grobs = sapply(names(inferred_seqs), plot_cdr3_lengths, seqs=input_sequences)
  cdr3_distribution_grobs = cdr3_distribution_grobs[!is.na(cdr3_distribution_grobs)]

  if('SEQUENCE_IMGT' %in% names(input_sequences)) {
    if(segment == 'V') {
      recs = lapply(names(inferred_seqs), function(x) {input_sequences[input_sequences$SEG_CALL==x,]$SEQUENCE_IMGT})
      refs = lapply(names(inferred_seqs), function(x) {inferred_seqs[x]})

      end_composition_grobs = mapply(plot_base_composition, names(inferred_seqs), recs, refs, pos=313, filter=T, all_inferred=all_inferred)
      end_composition_grobs = end_composition_grobs[!is.na(end_composition_grobs)]

      whole_composition_grobs = mapply(plot_base_composition, names(inferred_seqs), recs, refs, pos=1, filter=F, all_inferred=all_inferred)
      whole_composition_grobs = whole_composition_grobs[!is.na(whole_composition_grobs)]

      triplet_grobs = mapply(plot_trailing_triplet, names(inferred_seqs), recs, refs)
      triplet_grobs = triplet_grobs[!is.na(triplet_grobs)]

      # mapply(plot_base_heatmap, names(inferred_seqs), recs, refs, pos=1, end_pos=318)
    } else if(segment == 'J') {
      recs = lapply(names(inferred_seqs), function(x) {input_sequences[input_sequences$SEG_CALL==x,]$SEG_SEQ})
      refs = lapply(names(inferred_seqs), function(x) {inferred_seqs[x]})

      whole_composition_grobs = mapply(plot_segment_composition, names(inferred_seqs), recs, refs, pos=1, filter=F, r_justify=T)
      whole_composition_grobs = whole_composition_grobs[!is.na(whole_composition_grobs)]
    } else if(segment == 'D') {
      recs = lapply(names(inferred_seqs), function(x) {input_sequences[input_sequences$SEG_CALL==x,]$SEG_SEQ_ALIGNED})
      refs = lapply(names(inferred_seqs), function(x) {inferred_seqs[x]})

      whole_composition_grobs = mapply(plot_segment_composition, names(inferred_seqs), recs, refs, pos=1, filter=F)
      whole_composition_grobs = whole_composition_grobs[!is.na(whole_composition_grobs)]
    }
  }

  return(list('cdr3_dist'=cdr3_distribution_grobs, 'whole'=whole_composition_grobs, 'end'=end_composition_grobs, 'triplet'=triplet_grobs))
}


#' Create haplotyping plots
#' @export
#' @param segment one of V, D, J
#' @param haplo_details Data structure created by create_haplo_details
#' @return named list containing the following elements:
#' \tabular{ll}{
#'     a_allele_plot \tab   plot showing allele usage for each potential haplotyping gene \cr
#'     haplo_grobs \tab     differential plot of allele usage for each usable haplotyping gene \cr
#'}
make_haplo_grobs = function(segment, haplo_details) {
  a_props = haplo_details$a_props
  a_genes = haplo_details$a_genes
  sa = haplo_details$sa

  if(segment == 'V') {
    theme = theme(legend.title = element_blank(), plot.margin=margin(1,4,19,4, 'cm'), axis.text.x = element_text(angle = 270, hjust = 0,vjust=0.5))
  } else {
    theme = theme(legend.title = element_blank(), plot.margin=margin(4,1,15,1, 'cm'), axis.text.x = element_text(angle = 270, hjust = 0,vjust=0.5, size=7))
  }

  # color_brewer doesn't work with large numbers of categories

  if(length(unique(a_props$a_allele)) > 6) {
    a_allele_plot = ggplot() + geom_bar(aes(y=percent, x=a_gene, fill=a_allele), data=a_props, stat='identity') +
      labs(x='Gene',
           y='Allele %',
           title='Allele Usage') +
      theme_classic(base_size=15) +
      theme
  } else {
    a_allele_plot = ggplot() + geom_bar(aes(y=percent, x=a_gene, fill=a_allele), data=a_props, stat='identity') +
      scale_fill_brewer(palette='Dark2') +
      labs(x='Gene',
           y='Allele %',
           title='Allele Usage') +
      theme_classic(base_size=15) +
      theme
  }

  haplo_grobs = lapply(a_genes, plot_differential, a_props=a_props, sa=sa, segment=segment)
  haplo_grobs = haplo_grobs[!is.na(haplo_grobs)]

  return(list('haplo'=haplo_grobs, 'aplot'=a_allele_plot))
}


#' Create the OGRDB style plot file
#' @export
#' @param filename name of file to create (pdf)
#' @param input_sequences the input_sequences data frame
#' @param cdr3_dist_grobs cdr3 length distribution grobs created by make_novel_base_grob
#' @param end_composition_grobs end composition grobs created by make_novel_base_grobs
#' @param whole_composition_grobs whole composition grobs created by make_novel_base_grobs
#' @param triplet_composition_grobs triplet composition grobs created by make_novel_base_grobs
#' @param barplot_grobs barplot grobs created by make_barplot_grons
#' @param a_allele_plot a_allele_plot grob created by make_haplo_grobs
#' @param haplo_grobs haplo_grobs created by make_haplo_grobs
#' @param message text message to display at end of report
#' @return nothing
write_plot_file = function(filename, input_sequences, cdr3_dist_grobs, end_composition_grobs, whole_composition_grobs, triplet_composition_grobs, barplot_grobs, a_allele_plot, haplo_grobs, message) {
  # Save all graphics to plot file

  x=pdf(filename, width=210/25,height=297/25)

  if('SEQUENCE_IMGT' %in% names(input_sequences)) {
    if(length(end_composition_grobs) > 0) {
      x=print(marrangeGrob(end_composition_grobs, nrow=3, ncol=2,top=NULL))
    }
    if(length(whole_composition_grobs) > 0) {
      x=print(marrangeGrob(whole_composition_grobs, nrow=3, ncol=1,top=NULL))
    }
    if(length(triplet_composition_grobs) > 0) {
      x=print(marrangeGrob(triplet_composition_grobs, nrow=3, ncol=1,top=NULL))
    }
  }

  if(length(barplot_grobs) > 0) {
    x = print(marrangeGrob(barplot_grobs, nrow=3, ncol=3,top=NULL))
  }

  if(length(cdr3_dist_grobs) > 0) {
    x=print(marrangeGrob(cdr3_dist_grobs, nrow=3, ncol=3,top=NULL))
  }

  #x=print(marrangeGrob(snap_composition_grobs, nrow=3, ncol=1,top=NULL))
  grid.arrange(a_allele_plot)
  x=print(marrangeGrob(haplo_grobs, nrow=1, ncol=1,top=NULL))

  if(str_length(message) > 0) {
    grid.arrange(textGrob(message, just="left", x = unit(0, "npc"), y = unit(0.9, "npc"), gp=gpar(fontsize=5)))
  }

  x=dev.off()

}


# -- Individual allele bar charts --

plot_allele_seqs = function(allele, s, inferred_seqs, genotype, segment) {
  g = genotype[genotype$sequence_id==allele,]
  recs = s[s$SEG_CALL==allele,]
  recs = recs[recs$SEG_MUT_NC < 21,]

  if(nrow(recs) == 0) {
    return(NA)
  }

  if(is.na(g$unmutated_sequences)) {
    g$unmutated_sequences = 0
  }

  if(g$unmutated_sequences != 0) {
    if(segment == 'V') {
      label_text = paste0(g$unmutated_sequences, ' (', round(g$unmutated_sequences*100/g$sequences, digits=1), '%) exact matches, in which:\n',
                          g$unique_cdr3s_unmutated, ' unique CDR3\n',
                          g$unique_js_unmutated, ' unique J')
    } else {
      label_text = paste0(g$unmutated_sequences, ' (', round(g$unmutated_sequences*100/g$sequences, digits=1), '%) exact matches, in which:\n',
                          g$unique_cdr3s_unmutated, ' unique CDR3\n',
                          g$unique_vs_unmutated, ' unique V')
    }
  } else {
    label_text = 'No exact matches.\n'
  }

  g = ggplot(data=recs, aes(x=SEG_MUT_NC)) +
    geom_bar(width=1.0) +
    labs(x='Nucleotide Difference',
         y='Count',
         title=allele,
         subtitle=paste0(g$sequences, ' sequences assigned\n', label_text)) +
    theme_classic(base_size=12) +
    theme(aspect.ratio = 1/1, plot.subtitle=element_text(size=8))


  return(ggplotGrob(g))
}

# -- CDR3 length distribution --

plot_cdr3_lengths = function(allele, seqs) {
  r = seqs[seqs$SEG_CALL==allele,]
  r = r[r$SEG_MUT_NC==0,]

  if(nrow(r) == 0) {
    return(NA)
  }

  r$CDR3_LEN = sapply(r$CDR3_IMGT, nchar)
  r = r[!is.na(r$CDR3_LEN),]
  r$CDR3_LEN = round(r$CDR3_LEN/3)
  r$CDR3_LEN = as.factor(r$CDR3_LEN)

  g = ggplot(data=r, aes(x=CDR3_LEN)) +
    geom_bar(width=1.0) +
    labs(x='CDR3 AA Length (unmutated)',
         y='Count',
         title=allele,
         theme_classic(base_size=12)) +
           theme(aspect.ratio = 1/1, plot.subtitle=element_text(size=8))

  return(ggplotGrob(g))
}

# -- Base Composition --

nuc_at = function(seq, pos, filter) {
  if(length(seq) >= pos) {
    if(filter) {
      if(seq[pos] %in% c('N', 'X', '.', '-')) {
        return(NA)
      }
    }
    return(seq[pos])
  } else {
    return(NA)
  }
}

nucs_at = function(seqs, pos, filter) {
  if(filter) {
    ret = data.frame(pos=as.character(c(pos)), nuc=(factor(sapply(seqs, nuc_at, pos=pos, filter=filter), levels=c('A', 'C', 'G', 'T'))))
  } else {
    ret = data.frame(pos=as.character(c(pos)), nuc=(factor(sapply(seqs, nuc_at, pos=pos, filter=filter), levels=c('A', 'C', 'G', 'T', 'N', 'X', '.', '-'))))
  }
  ret = ret[!is.na(ret$nuc),]
  return(ret)
}

label_nuc = function(pos, ref) {
  return(paste0(pos, "\n", ref[[1]][pos]))
}


label_5_nuc = function(pos, ref) {
  if(pos %% 5 == 0) {
    n = pos
  } else {
    n = ''
  }

  return(paste0(n, "\n", ref[[1]][pos]))
}


# Plot base composition from nominated nucleotide position to the end or to optional endpos.
# Only include gaps, n nucleotides if filter=F
# if pos is negative, SEQUENCE_IMGT contains a certain number of trailing nucleotides. Plot them all.

plot_base_composition = function(gene_name, recs, ref, pos=1, filter=T, end_pos=999, r_justify=F, all_inferred=F) {
  max_pos = nchar(ref)


  if(max_pos < pos) {
    return(NA)
  }

  if(length(recs) < 1) {
    if (!all_inferred) {
      cat(paste0("Warning: no sequences found for ", gene_name, "(check SEQUENCE_IMGT in input file)"))
    }
    return(NA)
  }

  max_pos = min(max_pos, end_pos)
  min_pos = max(pos, 1)

  if(r_justify) {
    recs = str_pad(recs, max_pos-min_pos+1, 'left')
  }

  recs = strsplit(recs, "")
  ref = strsplit(as.character(ref), "")

  x = do.call('rbind', lapply(seq(min_pos,max_pos), nucs_at, seqs=recs, filter=filter))

  g = ggplot(data=x, aes(x=pos, fill=nuc)) +
    scale_fill_brewer(palette='Dark2') +
    geom_bar(stat="count") +
    labs(x='Position', y='Count', fill='', title=paste0('Gene ', gene_name)) +
    theme_classic(base_size=12) +
    scale_y_continuous(expand=c(0,0)
    )

  if(filter) {
    b =sapply(seq(pos, max_pos), label_5_nuc, ref=ref)
    g = g + scale_x_discrete(labels=b)
  } else {
    g = g +   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }

  return(ggplotGrob(g))
}

# As above, but plot a heatmap
plot_base_heatmap = function(gene_name, recs, ref, pos=1, end_pos=999, r_justify=F) {
  max_pos = nchar(ref)

  if(max_pos < pos || length(recs) < 1) {
    return(NA)
  }

  max_pos = min(max_pos, end_pos)
  min_pos = max(pos, 1)

  if(r_justify) {
    recs = str_pad(recs, max_pos-min_pos+1, 'left')
  }

  recs = strsplit(recs, "")
  ref = strsplit(as.character(ref), "")

  m=sapply(recs,function(x) {head(x,end_pos)})
  m=t(m)
  n=m[sample(nrow(m),size=200,replace=TRUE),]
  h=Heatmap(n, name='', clustering_distance_rows=function(x,y){sum(stringdist(x,y,method='hamming'))}, row_dend_width = unit(50, "mm"), column_title=gene_name)

  pdf(file=paste0('heatmap_', gsub('*', '_', gene_name, fixed=T), '.pdf'))
  draw(h)
  dev.off()
}


# Plot composition of a segment rather than the whole IMGT-aligned sequence
plot_segment_composition = function(gene_name, recs, ref, pos=1,  filter=T, end_pos=999, r_justify=F) {
  max_pos = nchar(ref)

  if(max_pos < pos || length(recs) < 1) {
    return(NA)
  }

  max_pos = min(max_pos, end_pos)
  min_pos = max(pos, 1)

  if(r_justify) {
    recs = str_pad(recs, max_pos-min_pos+1, 'left')
  }

  recs = strsplit(recs, "")
  ref = strsplit(as.character(ref), "")

  x = do.call('rbind', lapply(seq(min_pos,max_pos), nucs_at, seqs=recs, filter=filter))

  g = ggplot(data=x, aes(x=pos, fill=nuc)) +
    scale_fill_brewer(palette='Dark2') +
    geom_bar(stat="count") +
    labs(x='Position', y='Count', fill='', title=paste0('Gene ', gene_name)) +
    theme_classic(base_size=12) +
    scale_y_continuous(expand=c(0,0)
    )

  b =sapply(seq(pos, max_pos), label_5_nuc, ref=ref)
  g = g + scale_x_discrete(labels=b)

  return(ggplotGrob(g))
}

# -- Trailing triplet composition --

# convert specified section of a nucleotide sequence to a number. Presence of non-nucleotides results in 0

nuc_value = list()
nuc_value['T'] = 0
nuc_value['C'] = 1
nuc_value['A'] = 2
nuc_value['G'] = 3

seq_to_num = function(rec) {
  rec = strsplit(rec, "")[[1]]

  if(!all(rec=="T" || rec=="C" || rec=="A" || rec=="G")) {
    return(0)
  }

  rec = rev(rec)

  return(sum(unlist(lapply(seq(length(rec)), function(x) nuc_value[rec[x]][[1]]*(4**(x-1))))) + 1)
}


# Extract a sub-sequence of specified position and length. Return the sequence or some text if non-nucleotides are present.

extract_frag = function(rec, start_pos, length) {
  if(length(rec) < start_pos + length - 1) {
    return('Truncated')
  }

  rec = rec[start_pos : (start_pos+length - 1)]

  if(!all(rec=="T" | rec=="C" | rec=="A" | rec=="G")) {
    return('Non-nucleotide')
  }

  return(paste(rec, collapse=""))
}

# Plot occurrence of each possible nucleotide sequence in the trailing triplet

plot_trailing_triplet = function(gene_name, recs, ref) {
  start_pos = nchar(ref) - 2
  length = 3

  if(start_pos < 1 || length(recs) < 1) {
    return(NA)
  }

  recs = strsplit(recs, "")
  ref = strsplit(as.character(ref), "")

  x = do.call('rbind', lapply(recs, extract_frag, start_pos=start_pos, length=3))
  x = data.frame(triplet = x[,1])
  x = aggregate(x, by=list(x$triplet), FUN=length)
  names(x)=c("triplet", "count")
  x$triplet=as.character(x$triplet)
  x$order = sapply(x$triplet, seq_to_num)
  x = x[order(x$order),]
  x$triplet=factor(x$triplet, levels=x$triplet)


  g = ggplot(data=x, aes(x=triplet, y=count)) +
    geom_bar(stat="identity") +
    labs(x='Triplet', y='Count', fill='', title=paste0(gene_name, '- Final 3 nucleotides as a triplet')) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.ticks = element_blank(), legend.position=c(0.9, 0.9),
          axis.text=element_text(size=6), axis.title =element_text(size=15), axis.text.x = element_text(angle = 270, hjust = 0,vjust=0.5))

  return(ggplotGrob(g))
}


# -- Differential plot by allele usage - if we have good alleles for this gene --

plot_differential = function(gene, a_props, sa, segment) {
  ap = a_props[a_props$a_gene==gene,]
  ap = ap[order(ap$percent, decreasing=T),]

  if(nrow(ap) < 2 || ap[1,]$percent > 75 || ap[2,]$percent < 20) {
    return(NA)
  }

  if(segment == 'V') {
    margins = 1
  } else {
    margins = 4
  }


  a1 = ap[1,]$a_allele
  a2 = ap[2,]$a_allele
  recs = sa %>% filter(a_gene==gene) %>% filter(a_allele==a1 | a_allele==a2)
  recs = recs %>% select(SEG_CALL, a_allele) %>% group_by(SEG_CALL, a_allele) %>% summarise(count=n())
  recs$pos = sapply(recs$a_allele, function(x) {if(x==a1) {1} else {-1}})
  recs$count = recs$count * recs$pos

  # don't let columns get too broad
  ncols = nrow(recs)
  if(ncols < 15) {
    width = 21 - 2*margins
    width = width * ncols/15
    margins = (21 - width)/2
  }

  g = ggplot(recs, aes(x=SEG_CALL, y=count, fill=a_allele)) +
    geom_bar(stat='identity', position='identity') +
    scale_fill_brewer(palette='Dark2') +
    labs(x='',
         y='Count',
         title=paste0('Sequence Count by ', gene, ' allele usage'),
         fill = 'Allele') +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.ticks = element_blank(), legend.position=c(0.9, 0.9),
          axis.text=element_text(size=4), axis.title =element_text(size=15), axis.text.x = element_text(angle = 270, hjust = 0,vjust=0.5),
          plot.margin=margin(1,margins,15,margins, 'cm'))
  return(ggplotGrob(g))
}



