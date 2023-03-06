# This software is licensed under the CC BY-SA 4.0 licence: https://creativecommons.org/licenses/by-sa/4.0/


# Functions to provide reasonable sorting of gene names

gene_family = function(gene_name) {
  if(!grepl('-', gene_name, fixed=TRUE)) {
    return( '000')
  }
  fam = strsplit(gene_name, '-')[[1]][[1]]
  return(substr(fam, nchar(fam), nchar(fam)))
}

gene_number = function(gene_name) {
  if(!grepl('-', gene_name, fixed=TRUE) && grepl('_S', gene_name, fixed=TRUE)) {
    num = strsplit(gene_name, '_S')[[1]][[2]]
  } else {
    spl = strsplit(gene_name, '-')

    if(length(spl[[1]]) > 1) {
      num = spl[[1]][[2]]
    } else {
      num = spl[[1]][[1]]
    }
  }

  if(grepl('*', num, fixed=TRUE)) {
    num = strsplit(num, '*', fixed=TRUE)[[1]][[1]]
  }

  return(stringr::str_pad(num, 3, side='left', pad='0'))
}

allele_number = function(gene_name) {
  if(!grepl('*', gene_name, fixed=TRUE)) {
    return('000')
  }

  return(stringr::str_pad(strsplit(gene_name, '*', fixed=TRUE)[[1]][[2]], 3, side='left', pad='0'))
}

order_alleles = function(allele_names) {

  allele_names$family = sapply(allele_names$genes, gene_family)
  allele_names$number = sapply(allele_names$genes, gene_number)
  allele_names$allele = sapply(allele_names$genes, allele_number)

  alleles = unique(allele_names$allele)[order(unique(allele_names$allele))]
  allele_names$allele_ind = sapply(allele_names$allele, function(x){which(alleles==x)})

  families = unique(allele_names$family)[order(unique(allele_names$family))]
  allele_names$family_ind = sapply(allele_names$family, function(x){which(families==x)})

  numbers = unique(allele_names$number)[order(unique(allele_names$number))]
  allele_names$number_ind = sapply(allele_names$number, function(x){which(numbers==x)})

  allele_names$index = allele_names$allele_ind + 1000*allele_names$number_ind + 1000000*allele_names$family_ind
  return(order(allele_names$index))
}

sort_alleles = function(allele_vec) {
  allele_names=data.frame(genes=as.character(unique(allele_vec)), stringsAsFactors = FALSE)
  return(allele_names$genes[order_alleles(allele_names)])
}
