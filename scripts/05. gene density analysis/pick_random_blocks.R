#! /usr/bin/env Rscript

library(GenomicRanges)
library(tidyverse)
library(argparse)
library(fitdistrplus)
library(cowplot)

# Sample input
# args <- list(n=1, only_genome='HUMAN', GENOME=Sys.glob('orthofinder/chrom_orthofinder/*.chrom'), OUTPUT_DIR='euprymna', CLUSTERS='orthofinder/5_allspecies.clus')
# args <- list(n=1000, only_genome='aq', GENOME=c('inst/extdata/aq.chrom'), OUTPUT_DIR='large_samples', CLUSTERS='../processing/04_blocks/nmax5.clust')
handle_args <- function() {
  parser <- ArgumentParser(description="Pick random synteny blocks from a given genome.")
  parser$add_argument("-n", type="integer", default=100,
                      help="sample this many times [default: 100]")
  parser$add_argument("-o", "--plot-intervening", type="character",
                      help="plot the intervening genes diagnostic plot")
  parser$add_argument("-g", "--only-genome", type="character",
                      help="Sample blocks for only the given genome")
  parser$add_argument("-d", "--use-data-set", action="store_true",
                      help="If set, use the clusters in the stored data set from the package. CLUSTERS is interpreted as a genome name.")
  parser$add_argument("-s", "--strip-prefix", action="store_true",
                      help="If set, strip the gene name prefixes from the genome.")
  parser$add_argument("OUTPUT_DIR", type="character", help="a directory to output the randomized blocks")
  parser$add_argument("CLUSTERS", type="character", help="a clusters file containing the blocks of interest")
  parser$add_argument("GENOME", type="character", nargs='+',
                      help="a .chrom-formatted file containing the locations of all genes in the genome. Note: the file is assumed to be named where the prefix before the first '.' is the genome name")

  args <- parser$parse_args()

  args
}

main <- function() {
  args <- handle_args()

  # parse the blocks file
  if(args$use_data_set) {
    library(devtools)
    load_all('.')
    all_blocks <- load_data_set('all_blocks')
    block_objs <- all_blocks[[args$CLUSTERS]]$observed
    blocks <- map(block_objs, block_to_tibble) %>% bind_rows()
  } else {
    blocks <- parse_blocks(args$CLUSTERS)
  }

  # parse the genome files
  genomes <- map(args$GENOME, parse_chrom_file, remove_prefix=args$strip_prefix) %>%
    set_names(basename(args$GENOME) %>% str_split('[.]') %>% map(~.[1]))

  # pre-compute the gene neighbor to each transcript in the genomes
  # (this is faster than) calling it individually
  following_genes <- map(genomes, index_following_genes) %>%
    set_names(names(genomes))

  # if we are only doing this for one genome, subset the blocks
  if(!is.null(args$only_genome)) {
    blocks <- blocks %>% filter(species == args$only_genome)
  }

  # now compute the number of intervening genes for each block and
  # add this back to the blocks tibble
  pb <- progress_estimated(nrow(blocks))
  intervening <- blocks %>%
    dplyr::select(gene_names, species) %>%
    pmap(~intervening_sizes(..1, genomes[[..2]], pb=pb))
  blocks <- blocks %>% mutate(sizes=intervening)

  # assume the nmax parameter is the highest number of intervening genes
  # between any two genes in any block
  nmax <- blocks %>% dplyr::select(sizes) %>% unnest() %>% max()

  # get a fit of the models to the negative binomial distribution
  fit_models <- blocks %>% group_split(species) %>%
    set_names(group_keys(blocks, species) %>% mutate_if(is.factor, as.character) %>% flatten_chr()) %>%
    map(~get.fit(.$sizes))

  # plot the fits if we are doing that
  if(!is.null(args$plot_intervening)) {
    plot_fits(fit_models, blocks, nmax, output=args$plot_intervening)
  }

  # generate the samples
  map(unique(blocks$species), function(genome_name) {
    intervening.params <- if(is.null(fit_models[[genome_name]])) NULL else as.list(fit_models[[genome_name]]$estimate)
    blocks %>%
      filter(species == genome_name) %>%
      sample_blocks(genomes[[genome_name]], following_genes[[genome_name]], nmax,
                    verbose=T, n=args$n, intervening.params = intervening.params) %>%
      mutate_if(is.list, ~map_chr(., ~paste(., collapse=','))) %>%
      write_delim(file.path(args$OUTPUT_DIR, paste(genome_name, "tsv", sep=".")), delim="\t")
  })
}

plot_fits <- function(fit_models, blocks, nmax, output=NULL, width=9, height=7) {
  fits <- crossing(species=names(fit_models), intervening=seq(0,nmax)) %>%
    add_column(estimate=pmap_dbl(.,  ~dnbinom(..2, fit_models[[..1]]$estimate[1], fit_models[[..1]]$estimate[2])))

  data <- blocks %>% dplyr::select(species,sizes) %>% unnest() %>% group_by(species) %>%
     mutate(frequency=1/n()) %>% group_by(species,sizes) %>% summarize(frequency=sum(frequency), count=n()) %>%
     add_column(estimate=pmap_dbl(
       dplyr::select(., species, sizes),
       ~dnbinom(..2, fit_models[[..1]]$estimate[1], fit_models[[..1]]$estimate[2])))

  mses <- data %>% group_by(species) %>% summarize(mse=(sum(frequency-estimate)**2)/n())

  fit_models_tib <- tibble(species=names(fit_models),
                           r=map_dbl(fit_models, ~.$estimate[1]),
                           prob=map_dbl(fit_models, ~.$estimate[2])) %>%
    inner_join(mses) %>%
    add_column(
      label=pmap_chr(., function(species, r, prob, mse)
        paste("r =", round(r,2),"\np =", round(prob,2), "\nMSE =", formatC(mse, format="e", digits=2))))

  p <- data %>%
      ggplot(aes(x=sizes, y=frequency)) +
         geom_bar(stat='identity') +
         geom_point(aes(y=estimate), color='red') +
         geom_text(data=fit_models_tib, aes(x=1,y=1,vjust=1, hjust=0, label=label), size=3) +
         scale_x_continuous(breaks=seq(0,nmax)) +
         lims(y=c(0,1)) +
         facet_wrap(.~species)

  if(!is.null(output)) {
    ggsave(output, width=width, height=height)
  }
}

#' Create an index of following genes to speed up sampling
#'
#' @param genome a GRanges object containing the genes in the genome
#'
#' @return an environment mapping all gene names to their following gene name or NA if none exists
index_following_genes <- function(genome) {
  #strand(genome) <- Rle(c('+'), c(length(genome)))
  following_idxs <- precede(genome, ignore.strand=T)
  following_genes <- new.env(hash=T)
  map(seq(length(following_idxs)),
      ~assign(x=names(genome)[.],
              value=if(is.na(following_idxs[.])) NA else names(genome)[following_idxs[.]],
              envir=following_genes)
      ) %>% invisible()
  following_genes
}

#' Sample blocks from a genome.
#'
#' @param nonrandom.blocks a set of blocks to avoid overlapping in generating the random sample sapce
#' @param genome a GRanges object containing all genes in a genome
#' @param following_genes an environment containing names of the following genes in the genome (speedup)
#' @param nmax maximum number of intervening seqeunces
#' @param genes.subset if not NULL, only sample blocks which have the given genes (useful if you would, for example, want only genes which have an ortholog to another gene)
#' @param prefilter remove some genes from the genome before sampling (e.g. genes with no expression data)
#' @param intervening.rfunc random variable generator function for intervening genes (default is negative binomial)
#' @param intervening.params parameters for the random function. Pass NULL to avoid generating intervening genes
#' @param max.tries maximum tries to generate a block. if overstepped, the script is stopped.
#' @param seed random seed
#' @param n number of samples to generate
#' @param verbose echo the ETA while computing
#'
#' @return a list of synBlock objects. If n > 1, a list of lists of synBlock objects
sample_blocks <- function(nonrandom.blocks, genome, following_genes, nmax, prefilter=NULL, genes.subset=NULL, max.tries=1000000000,
                          intervening.rfunc=rnbinom, intervening.params=list(size=1, prob=1), seed=030918, n=1,
                          check.visited=F, verbose=T) {
  set.seed(seed)

  filtered <- if(is.null(prefilter)) genome else genome[names(genome) %in% as.character(prefilter)]

  sorted <- sortSeqlevels(filtered)
  sorted <- sort(sorted, ignore.strand=T)
  sorted.by.chr <- as.list(split(sorted, seqnames(sorted)))

  shortest.block <-  min(nonrandom.blocks$num_genes)
  unusable.seqnames <- names(which(sapply(sorted.by.chr, length) < shortest.block))

  # remove chromosomes which do not have enough genes on them to contain a block
  filtered <- dropSeqlevels(sorted, unusable.seqnames, pruning.mode = 'coarse')
  filtered <- sortSeqlevels(filtered)
  filtered <- sort(filtered, ignore.strand=T)
  filtered.by.chr <- as.list(split(filtered, seqnames(filtered)))

  # get all ranges
  if(check.visited) {
    visited <- suppressWarnings(do.call(c, sapply(unlist(nonrandom.blocks, use.names = F), grange)))
    seqlevels(visited) <- seqlevels(genome)
  }

  genome.size <- length(filtered)
  chrsizes <- sapply(filtered.by.chr, length)
  # index of the chrom implied by a random number from 1:genome.size
  chr.idxs <- rep(1:length(chrsizes), chrsizes)
  # the global start index of each scaffold
  chr.starts <- as.integer(cumsum(chrsizes))

  sample_block <- function(block, sample_id=1) {
    size <- block$num_genes
    for(i in seq(max.tries)) {
      idx <- sample.int(n=genome.size, size=1)
      chr.idx <- chr.idxs[idx]
      # the position chromosome starts at is the cumulative sum of all the chromosomes prior to the current one
      # find the index of this current chromosome to start on
      start.on.chr <- 1+idx- max(0,chr.starts[chr.idx-1])
      cur.chr <- filtered.by.chr[[chr.idx]]

      # get all indices of the block by sampling the poisson distribution as many times as the size of the block
      .next.size <- function() { n <- do.call(intervening.rfunc, args=c(n=1,intervening.params)) ; if(n <= nmax)  n  else .next.size() }
      intervening.sizes <- replicate(size-1, if(is.null(intervening.params)) 1 else .next.size())

      # get the next genes
      genes <- reduce(intervening.sizes,
                   function(genes, size)  {
                     next_gene <- reduce(
                       1:(size+1),
                       function(gene, n) {
                         .next = if(is.na(gene)) NA else following_genes[[gene]]
                       },
                      .init=genes[length(genes)])
                     if(is.na(next_gene)) NA else c(genes, next_gene)
                   },
                   .init=names(genome)[idx])

      # if this goes off the boundaries of the chromosome then give up
      if(anyNA(genes)) { next }
      gene_objs <- genome[genes]
      rchrom <- as.character(seqnames(gene_objs))[1]
      rstart <- min(start(gene_objs))
      rend <- max(end(gene_objs))
      gr <- GRanges(rchrom, ranges=IRanges(start=rstart, end=rend))

      if(check.visited) {
        if(length(findOverlaps(gr, visited))) {  next }

        # add the interval to the overlaps
        visited <<- suppressWarnings(c(visited, gr))
      }

      # check if all the genes are present in the genes.subset, if given
      if(!is.null(genes.subset) && !all(sapply(gene.ids, function(g) g %fin% genes.subset))) { next }

      # return the randomized block
      return(mutate(block,
                    block_id=as.character(str_glue('{block$block_id}.{sample_id}')),
                    linked_blocks="None(Random)",
                    linked_species="None(Random)",
                    num_linked_species=0,
                    sizes=list(intervening.sizes),
                    location=paste0(as.character(rchrom), ':', rstart, '..', rend),
                    length=size,
                    gene_names=list(genes)))
    }
    warning(paste("Could not find a block for", name(block)), "after", max.tries, "tries.")
  }

  block_list <- nonrandom.blocks %>% group_split(block_id)
  if(n > 1) {
    if(verbose) pb <- progress_estimated(n)
    res <- map(1:n, function(x) {
      df <- map(block_list, sample_block, x) %>% bind_rows()
      if(verbose) { pb$tick()$print(); flush.console() }
      df
    }) %>% bind_rows()
  } else {
    res <- map(block_list, sample_block) %>% bind_rows()
  }

  res
}


#' Parse a chrom file into a GRanges object.
#'
#' @param filename name of the file
#' @param remove_prefix remove the genome prefix from the gene name ()
#'
#' @return GRanges object
#'
#' @export
parse_chrom_file <- function(filename, remove_prefix=F) {
  tab <- read.table(filename, sep='',
                    col.names = c('genome', 'geneName', 'chrom', 'strand', 'start', 'end'))
  if(remove_prefix) {

    tab$geneName <- remove_gene_prefixes(tab$geneName)
  }
  ranges <- IRanges(start=tab$start, end=tab$end, names=tab$geneName)

  # prepare for synteny analysis: sort sequence levels and sort by gene position, ignoring strand
  gs <- GRanges(seqnames=tab$chrom, ranges=ranges, strand=tab$strand, genome=tab$genome)
  gs <- sortSeqlevels(gs)
  gs <- sort(gs, ignore.strand=T)
  gs
}

#' Parse a cluster file into a tibble
#'
#' @param filename a cluster-formatted file containing all blocks of interest
#'
#' @return a tibble containing the contents of the file. The tibble's gene names are nested into a list.
parse_blocks <- function(filename) {
  read_delim(filename,
                    col_names=c("block_id","species",
                                "num_links","linked_blocks",
                                "num_linked_species","linked_species",
                                "strand","location","length","gene_names"),
                    comment="#",
                    delim = "\t") %>%
    mutate(gene_names=str_split(gene_names, ",")) %>%
    mutate(num_genes=map_int(gene_names,~length(.x)))
}

#' Convert a ablock object to a tibble
#'
#' @param block a synBlock object from the SCS library
#'
#' @return a tibble with one row containing a the converted block object. (several tibbles can in turn be bound into a single tibble)
block_to_tibble <- function(block) {
  tibble(block_id=str_split(block@.Name, '\\(', simplify=T)[[1]],
         species=block@.Species,
         num_links=length(block@.Conns),
         linked_blocks=if(length(block@.Conns)) str_flatten(block@.Conns, ',') else "",
         num_linked_species=length(block@.ConnSpecies),
         linked_species=str_flatten(Filter(Negate(is.null), block@.ConnSpecies), ','),
         strand='+',
         location=block@.Location,
         length=block@.Length,
         gene_names=list(genes(block)),
         num_genes=length(genes(block))
         )
}

#' Fit data to the negative binomial distribution
#'
#' This is a convenience function to fit the negative binomial distribution as done for intervening sizes.
#'
#' @param data a list of data points to fit.
#'
#' @return a fitdist object a la fitdistrplus or NULL if the fitting failed
get.fit <- function(data) {
  tryCatch(
    {
      fitdist(unlist(data), distr='nbinom', method='mle', lower=c(0,0),
              start=list(size=.Machine$double.eps,prob=.Machine$double.eps))
    },
    error = function(e) NULL)
  
}

plot.fit <- function(data, fit, name, end=4) {
  fits <- sapply(0:end, function(x) dnbinom(x, fit$estimate[1], fit$estimate[2]))

  ggplot(data.frame(InterveningSize=unlist(data)[which(data <= 4)]), aes(InterveningSize, stat(density))) +
    geom_histogram(bins=end+1) + lims(y=c(0,1)) +
    geom_point(data=data.frame(x=0:end,y=fits), aes(x=x,y=y),color='red') + ggtitle(name) +
    annotate('text', x=0, y=1, vjust=1, hjust=0, size=4.5,
             label=paste('r =', round(fit$estimate[1],2), '\np =', round(fit$estimate[2],2)))
}

#' Get the number of intervening genes between each gene in a synteny block
#'
#' @param gene_names names of genes which appear in the genome
#' @param genome the underlying GRanges object
#'
intervening_sizes <- function(gene_names, genome, pb=NULL) {
  ranges <- genome[unlist(gene_names)]
  ranges <- sort(ranges, ignore.strand=TRUE)
  ranges_seqnames <- unique(seqnames(ranges))
  if(length(ranges_seqnames) > 1) {
    stop(paste("Block with genes",
               paste(gene_names, collapse=", "),
               "has multiple seqnames",
               paste(ranges_seqnames, collapse=", ")))
  }
  overlap_free_genome <- append(
    subsetByOverlaps(genome[seqnames(genome) == ranges_seqnames[1]], ranges, invert=T),
    ranges)
  overlap_free_genome <- sort(overlap_free_genome, ignore.strand=TRUE)
  res <- sapply(1:(length(ranges)-1),
         function(i) {
           which(names(overlap_free_genome) == names(ranges)[i+1]) -
             which(names(overlap_free_genome) == names(ranges)[i]) - 1
         }
  )
  if(!is.null(pb)) { pb$tick()$print() }
  res
}

if(sys.nframe() == 0) {
  main()
}
