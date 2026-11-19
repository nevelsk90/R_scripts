extendGenes <- function(genes_gr,
                        prm_up        = 2000,
                        upstream_max  = 250000,
                        downstream_max = 5000,
                        seqlens        = seqlengths(si)) {
  
  stopifnot(is(genes_gr, "GRanges"))
  if (is.null(seqlens))
    stop("Please provide chromosome lengths via 'seqlens' or seqlengths(genes_gr)")
  
  ## 1. sort once per chromosome --------------------------------------------
  ord <- order(seqnames(genes_gr), start(genes_gr))
  g <- genes_gr[ord]
  
  prev_end   <- c(NA, end(g)[-length(g)])
  next_start <- c(start(g)[-1], NA)
  
  same_chr      <- as.character(seqnames(g))
  same_chr_prev <- c(FALSE, same_chr[-1]  == same_chr[-length(g)])
  same_chr_next <- c(same_chr[-1] == same_chr[-length(g)], FALSE)
  
  prev_end[ !same_chr_prev ]   <- NA
  next_start[ !same_chr_next ] <- NA
  
  ## 2. build extended regions ----------------------------------------------
  out <- vector("list", length(g))
  
  for (i in seq_along(g)) {
    chr    <- as.character(seqnames(g[i]))
    s_len  <- seqlens[chr]
    g_str  <- as.character(strand(g[i]))
    
    ## promoter --------------------------------------------------------------
    if (g_str == "-") {
      tss        <- end(g[i])
      prom_start <- tss + 1
      prom_end   <- min(tss + prm_up, s_len)
    } else {                           # "+" or "*"
      tss        <- start(g[i])
      prom_start <- max(1, tss - prm_up)
      prom_end   <- tss - 1
    }
    promoter <- GRanges(chr, IRanges(prom_start, prom_end), strand = g_str)
    
    ## upstream (before promoter) -------------------------------------------
    if (g_str == "-") {                # grows to the right
      boundary <- if (is.na(next_start[i])) s_len else next_start[i] - 1
      up_end   <- min(prom_end + upstream_max, boundary)
      up_start <- prom_end + 1
      upstream <- if (up_end >= up_start)
        GRanges(chr, IRanges(up_start, up_end), strand = g_str) else GRanges()
    } else {                           # "+" strand -> grows to the left
      boundary <- if (is.na(prev_end[i])) 1 else prev_end[i] + 1
      up_start <- max(boundary, prom_start - upstream_max)
      up_end   <- prom_start - 1
      upstream <- if (up_end >= up_start)
        GRanges(chr, IRanges(up_start, up_end), strand = g_str) else GRanges()
    }
    
    ## downstream (after gene end) ------------------------------------------
    if (g_str == "-") {                # grows to the left
      boundary   <- if (is.na(prev_end[i])) 1 else prev_end[i] + 1
      down_start <- max(boundary, start(g[i]) - downstream_max)
      down_end   <- start(g[i]) - 1
      downstream <- if (down_end >= down_start)
        GRanges(chr, IRanges(down_start, down_end), strand = g_str) else GRanges()
    } else {                           # "+" strand -> grows to the right
      boundary   <- if (is.na(next_start[i])) s_len else next_start[i] - 1
      down_start <- end(g[i]) + 1
      down_end   <- min(end(g[i]) + downstream_max, boundary)
      downstream <- if (down_end >= down_start)
        GRanges(chr, IRanges(down_start, down_end), strand = g_str) else GRanges()
    }
    
    ## union everything for this gene ---------------------------------------
    pieces <- c(g[i], promoter, upstream, downstream)
    ext    <- GenomicRanges::reduce(pieces, min.gapwidth = 1)   # merge adjacent parts
    # ext    <- ext[1]                             # should be a single range
    names(ext) <- names(g[i])
    out[[i]]  <- ext
  }
  
  ## restore original gene order --------------------------------------------
  extended <- do.call( base::c, out)
  extended <- extended[ match(names(genes_gr), names(extended)) ]
  mcols(extended) <- mcols(genes_gr)             # keep metadata
  
  extended
}