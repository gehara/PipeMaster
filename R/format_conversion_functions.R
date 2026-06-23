# internal function
# @description transforms fasta alignments in ms-like.
# @export
fasta2ms<-function(path.to.fasta,fasta.files,write.file=T){

  setwd(path.to.fasta)
  ms.out<-list()
  for(u in 1:length(fasta.files)){
    fas<-read.dna(file=fasta.files[u], format="fasta")
    fas<-as.character(fas)
    bin<-NULL
    pos<-NULL

    if(ncol(fas)==0){
      stop(paste(fasta.files[u],"has 0 base pairs! Check the alignment"))
    }

    for(i in 1:ncol(fas)){
      a<-length(grep("a",fas[,i]))
      c<-length(grep("c",fas[,i]))
      g<-length(grep("g",fas[,i]))
      t<-length(grep("t",fas[,i]))
      n<-length(grep("n",fas[,i]))
      gap<-length(grep("-",fas[,i]))
      if(nrow(fas) %in% c(a,c,g,t)){
      } else if (gap>0){
      } else if(n>0){
      } else {bin<-cbind(bin,fas[,i])
      pos<-c(pos,i)
      }
    }

    if(is.null(pos)){
      if(write.file==T){
        write(file=paste(strsplit(fasta.files[u],".fas")[[1]],".ms",sep=""),paste("ms",nrow(fas),1))
        write(file=paste(strsplit(fasta.files[u],".fas")[[1]],".ms",sep=""),"//",append=T)
        write(file=paste(strsplit(fasta.files[u],".fas")[[1]],".ms",sep=""),paste("segsites:",0),append=T)
        write(file=paste(strsplit(fasta.files[u],".fas")[[1]],".ms",sep=""),paste("positions:"),append=T)
      }
      next
    }

    pos<-pos/ncol(fas)
    for(j in 1:ncol(bin)){
      a<-length(grep("a",bin[,j]))/nrow(bin)
      c<-length(grep("c",bin[,j]))/nrow(bin)
      g<-length(grep("g",bin[,j]))/nrow(bin)
      t<-length(grep("t",bin[,j]))/nrow(bin)
      bases<-c(a,c,g,t)
      names(bases)<-c("a","c","g","t")
      base<-which.max(bases)
      bin[,j]<-gsub(names(base),0,bin[,j])
      bases<-bases[-base]
      for(i in 1:3){
        bin[,j]<-gsub(names(bases[i]),1,bin[,j])
      }
    }

    seqs<-NULL
    for(i in 1:nrow(bin)){
      seqs<-c(seqs,paste(bin[i,],collapse=""))
    }

    if(write.file==T){
      write(file=paste(strsplit(fasta.files[u],".fas")[[1]],".ms",sep=""),paste("ms",nrow(fas),1))
      write(file=paste(strsplit(fasta.files[u],".fas")[[1]],".ms",sep=""),"//",append=T)
      write(file=paste(strsplit(fasta.files[u],".fas")[[1]],".ms",sep=""),paste("segsites:",ncol(bin)),append=T)
      write(file=paste(strsplit(fasta.files[u],".fas")[[1]],".ms",sep=""),paste("positions:   ",paste(pos,collapse="    ")),append=T)
      write(file=paste(strsplit(fasta.files[u],".fas")[[1]],".ms",sep=""),seqs,sep="\n",append=T)
    }

    if(write.file==F){
      ms.out[[u]]<-paste("ms",nrow(fas),1)
      ms.out[[u]]<-c(ms.out[[u]],paste("//"))
      ms.out[[u]]<-c(ms.out[[u]],paste("segsites:",ncol(bin)))
      ms.out[[u]]<-c(ms.out[[u]],paste("positions:   ",paste(pos,collapse="    ")))
      ms.out[[u]]<-c(ms.out[[u]],paste(seqs,sep="\n"))
    }
  }
  if(write.file==F){
    return(ms.out)
  }
}

# internal function of the obs.snp.sumstat function
# @description transforms fasta alignments in ms-like.
# @export
fasta.snp.2ms<-function(path.to.fasta,fasta.files,write.file=T,pop.assign){

  ms.out<-list()
  for(u in 1:length(fasta.files)){
    fas<-read.dna(file=fasta.files[u], format="fasta")
    fas<-as.character(fas)
    if(is.matrix(fas)==F){
      stop(cat(paste("Something is wrong with alignment",fasta.files[u]),
               paste("Some potential problems:"),
               paste("1) sequences are not aligned"),
               paste("2) uknown character in the alignemt (only IUPAC nucleotide codes allowed, no question marks)"),
               paste("3) something else..."),sep="\n"))
    }

    pops<-pop.assign
    pops<-pops[with(pops, order(pops[,2])), ]


    if(length(grep(-9,match(rownames(fas),pops[, 1],nomatch=-9)))>0){
      rownames(fas)[grep(-9,match(rownames(fas),pops[, 1],nomatch=-9))]
      stop(cat(paste("There is at least one sequence in the alignment",fasta.files[u],"without assignment."),
               paste("The sequence name", rownames(fas)[grep(-9,match(rownames(fas),pops[, 1],nomatch=-9))]),
               paste("has no match in the assignment file."),sep="\n"))
    }

    if(length(unique(pops[pops[, 1] %in% rownames(fas),2])) != length(unique(pops[,2]))){
      stop(paste("locus",fasta.files[u],"does not have samples for all",length(unique(pops[,2])),"populations"))
    }



    fasta<-NULL
    p<-NULL
    for (j in 1:nrow(pops)) {
      x <- match(rownames(fas),pops[j, 1])
      x <- which(x==1)
      if (length(x)!=0) {
        p <- rbind(p, pops[rep(j,length(x)), ])
        fasta <- rbind(fasta, fas[x, ])
      }
    }
    fas<-fasta
    #ape:::write.dna(as.DNAbin(fas),fasta.files[u],format = "fasta", colw = 10000)
    pops<-p

    npops<-length(unique(pops[,2]))
    pop.list<-list()
    for(i in 1:npops){
      pop.list[[i]]<-length(which(pops[,2]==i))
    }

    string<-paste("-I",npops,paste(unlist(pop.list),collapse=" "))

    # keep only variable sites and get their position in the alignment
    bin<-NULL
    pos<-NULL
    for(i in 1:ncol(fas)){
      a<-length(grep("a",fas[,i]))
      c<-length(grep("c",fas[,i]))
      g<-length(grep("g",fas[,i]))
      t<-length(grep("t",fas[,i]))

      r<-length(grep("r",fas[,i]))
      y<-length(grep("y",fas[,i]))
      m<-length(grep("m",fas[,i]))
      k<-length(grep("k",fas[,i]))
      s<-length(grep("s",fas[,i]))
      w<-length(grep("w",fas[,i]))

      h<-length(grep("h",fas[,i]))
      b<-length(grep("b",fas[,i]))
      v<-length(grep("v",fas[,i]))
      d<-length(grep("d",fas[,i]))
      n<-length(grep("n",fas[,i]))
      gap<-length(grep("-",fas[,i]))
      if(gap>0){next
      } else if (n>0){next}
      if(!nrow(fas) %in% c(a,c,g,t,r,y,m,k,s,w,h,b,v,d)){
        bin<-cbind(bin,fas[,i])
        pos<-c(pos,i)
      }
    }
    pos<-pos/ncol(fas)

    if(!(is.null(bin))){
      for(j in 1:ncol(bin)){
        g<-length(grep("g",bin[,j]))/nrow(bin)
        a<-length(grep("a",bin[,j]))/nrow(bin)
        t<-length(grep("t",bin[,j]))/nrow(bin)
        c<-length(grep("c",bin[,j]))/nrow(bin)

        r<-length(grep("r",bin[,j]))/nrow(bin)
        y<-length(grep("y",bin[,j]))/nrow(bin)
        m<-length(grep("m",bin[,j]))/nrow(bin)
        k<-length(grep("k",bin[,j]))/nrow(bin)
        s<-length(grep("s",bin[,j]))/nrow(bin)
        w<-length(grep("w",bin[,j]))/nrow(bin)

        h<-length(grep("h",bin[,j]))/nrow(bin)
        b<-length(grep("b",bin[,j]))/nrow(bin)
        v<-length(grep("v",bin[,j]))/nrow(bin)
        d<-length(grep("d",bin[,j]))/nrow(bin)
        bases<-c(a,c,g,t,r,y,m,s,k,w,h,b,v,d)
        names(bases)<-c("a","c","g","t","r","y","m","s","k","w","h","b","v","d")
        base<-which.max(bases[1:4])
        bin[,j]<-gsub(names(base),0,bin[,j])
        bases<-bases[-base]
        for(i in 1:length(bases)){
          bin[,j]<-gsub(names(bases[i]),1,bin[,j])
        }
      }
      seqs<-NULL
      for(i in 1:nrow(bin)){
        seqs<-c(seqs,paste(bin[i,],collapse=""))
      }
      ss<-ncol(bin)

      x<-as.vector(bin)
      if(length(c(grep("0",x),grep("1",x)))!=length(x)){
        stop(paste("Something is wrong with alignment",fasta.files[u]))
      }

      ### if there is no variation
    }else{seqs<-NULL
    ss<-0}


    if(write.file==T){
      write(file=paste(strsplit(fasta.files[u],".",fixed=T)[[1]][1],".ms",sep=""),paste("ms",nrow(fas),1,string))
      write(file=paste(strsplit(fasta.files[u],".",fixed=T)[[1]][1],".ms",sep=""),"//",append=T)
      write(file=paste(strsplit(fasta.files[u],".",fixed=T)[[1]][1],".ms",sep=""),paste("segsites:",ss),append=T)
      write(file=paste(strsplit(fasta.files[u],".",fixed=T)[[1]][1],".ms",sep=""),paste("positions:   ",paste(pos,collapse="    ")),append=T)
      write(file=paste(strsplit(fasta.files[u],".",fixed=T)[[1]][1],".ms",sep=""),seqs,sep="\n",append=T)
    }
    if(npops>1){
      ms.out[[u]]<-paste("ms",nrow(fas),1,string)
    } else { ms.out[[u]]<-paste("ms",nrow(fas),1)}
    ms.out[[u]]<-c(ms.out[[u]],paste("//"))
    ms.out[[u]]<-c(ms.out[[u]],paste("segsites:",ss))
    ms.out[[u]]<-c(ms.out[[u]],paste("positions:   ",paste(pos,collapse="    ")))
    ms.out[[u]]<-c(ms.out[[u]],paste(seqs,sep="\n"))

 # print(u)
  }
  return(ms.out)
}

#' Read a multi-locus sequential PHYLIP file into per-locus character matrices
#'
#' @description Reads a multi-locus sequential PHYLIP file and returns the loci
#'   as a list of character matrices (samples x sites), equivalent to
#'   \code{ape::read.dna(format = "fasta")} followed by \code{as.character()}.
#'   Each locus block starts with a \code{ntax nchar} header line. Tab-delimited
#'   sample names (needed for sample IDs longer than the standard 10 characters)
#'   are supported.
#'
#' @param filepath Path to the sequential PHYLIP file.
#' @return A list of character matrices, one per locus. Row names are sample
#'   names; columns are individual nucleotide positions (lowercase).
#' @author Marcelo Gehara
#' @seealso \code{\link{alleles2phylip}}, \code{\link{observed.sumstats}}
#' @export
read.phylip.loci <- function(filepath) {
  lines <- readLines(filepath)
  n_lines <- length(lines)
  loci <- list()
  i <- 1
  locus_idx <- 0

  while(i <= n_lines) {
    # Skip blank lines
    if(nchar(trimws(lines[i])) == 0) {
      i <- i + 1
      next
    }

    # Parse header: " ntax nchar"
    header <- as.integer(strsplit(trimws(lines[i]), "\\s+")[[1]])
    ntax <- header[1]
    nchar_seq <- header[2]
    i <- i + 1
    locus_idx <- locus_idx + 1

    # Read ntax sequence lines.
    # Detect tab-delimited (long names) vs fixed-width (10-char names),
    # matching the dual-mode parser in phylip_to_ms_file()
    # (R/observed.sumstat.ngs.R:122-131).
    names_vec <- character(ntax)
    mat <- matrix("", nrow = ntax, ncol = nchar_seq)
    for(j in 1:ntax) {
      line <- lines[i]
      if (grepl("\t", line, fixed = TRUE)) {
        parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
        name <- trimws(parts[1])
        seq_str <- gsub("\\s", "", parts[2])
      } else {
        # Try relaxed PHYLIP first (name + whitespace run + sequence).
        # Fall back to classic fixed-width (10-char name) only if relaxed
        # split doesn't yield the expected sequence length.
        m <- regexpr("\\s+", line)
        relaxed_ok <- FALSE
        if (m > 1L) {
          cand_name    <- trimws(substr(line, 1, m - 1L))
          cand_seq_str <- gsub("\\s", "", substr(line, m, nchar(line)))
          if (nchar(cand_seq_str) == nchar_seq) {
            name       <- cand_name
            seq_str    <- cand_seq_str
            relaxed_ok <- TRUE
          }
        }
        if (!relaxed_ok) {
          name <- trimws(substr(line, 1, 10))
          seq_str <- gsub("\\s", "", substr(line, 11, nchar(line)))
        }
      }
      names_vec[j] <- name
      mat[j, ] <- strsplit(tolower(seq_str), "")[[1]]
      i <- i + 1
    }
    rownames(mat) <- names_vec
    loci[[locus_idx]] <- mat
  }
  loci
}

# internal function
# @description Converts a character matrix (samples x sites) to ms-like format.
#   This is the PHYLIP equivalent of fasta.snp.2ms(): instead of reading a file,
#   it takes an already-parsed alignment matrix (e.g., from read.phylip.loci()).
# @param fas Character matrix with nucleotide sequences (samples as rows, sites as columns).
#   Row names must be sample names.
# @param pop.assign Two-column data frame: col1 = sample names, col2 = population numbers.
# @return A list with one element: a character vector of ms-like output lines.
mat.snp.2ms <- function(fas, pop.assign) {

  pops <- data.frame(pop.assign)
  pops <- pops[with(pops, order(pops[,2])), ]

  # Reorder rows by population assignment
  fasta <- NULL
  p <- NULL
  for(j in 1:nrow(pops)) {
    x <- match(rownames(fas), pops[j, 1])
    x <- which(x == 1)
    if(length(x) != 0) {
      p <- rbind(p, pops[rep(j, length(x)), ])
      fasta <- rbind(fasta, fas[x, ])
    }
  }
  fas <- fasta
  pops <- p

  npops <- length(unique(pops[,2]))
  pop.list <- list()
  for(i in 1:npops) {
    pop.list[[i]] <- length(which(pops[,2] == i))
  }
  string <- paste("-I", npops, paste(unlist(pop.list), collapse = " "))

  # Identify variable sites (exclude gaps, N, monomorphic)
  bin <- NULL
  pos <- NULL
  for(i in 1:ncol(fas)) {
    a <- length(grep("a", fas[,i]))
    c <- length(grep("c", fas[,i]))
    g <- length(grep("g", fas[,i]))
    t <- length(grep("t", fas[,i]))

    n <- length(grep("n", fas[,i]))
    gap <- length(grep("-", fas[,i]))
    if(gap > 0) next
    if(n > 0) next
    if(!nrow(fas) %in% c(a,c,g,t)) {
      bin <- cbind(bin, fas[,i])
      pos <- c(pos, i)
    }
  }
  pos <- pos / ncol(fas)

  if(!is.null(bin)) {
    # Encode: major allele -> 0, everything else -> 1
    for(j in 1:ncol(bin)) {
      a <- length(grep("a", bin[,j])) / nrow(bin)
      c <- length(grep("c", bin[,j])) / nrow(bin)
      g <- length(grep("g", bin[,j])) / nrow(bin)
      t <- length(grep("t", bin[,j])) / nrow(bin)
      bases <- c(a, c, g, t)
      names(bases) <- c("a", "c", "g", "t")
      base <- which.max(bases[1:4])
      bin[,j] <- gsub(names(base), 0, bin[,j])
      bases <- bases[-base]
      for(i in 1:length(bases)) {
        bin[,j] <- gsub(names(bases[i]), 1, bin[,j])
      }
    }
    seqs <- NULL
    for(i in 1:nrow(bin)) {
      seqs <- c(seqs, paste(bin[i,], collapse = ""))
    }
    ss <- ncol(bin)
  } else {
    seqs <- NULL
    ss <- 0
  }

  ms.out <- list()
  if(npops > 1) {
    ms.out[[1]] <- paste("ms", nrow(fas), 1, string)
  } else {
    ms.out[[1]] <- paste("ms", nrow(fas), 1)
  }
  ms.out[[1]] <- c(ms.out[[1]], "//")
  ms.out[[1]] <- c(ms.out[[1]], paste("segsites:", ss))
  ms.out[[1]] <- c(ms.out[[1]], paste("positions:   ", paste(pos, collapse = "    ")))
  ms.out[[1]] <- c(ms.out[[1]], paste(seqs, sep = "\n"))

  return(ms.out)
}

#' Converts ms output to DNAbin file
#' @description This function will take a ms-like output and convert it to DNAbin format.
#' @param ms.output A list of strings representing a ms output.
#' @param bp.length The number of base pairs used to calculate theta for the ms simulation.
#' @return a DNAbin object
#' @note This function is used internally by all the main functions used to simulate coexpantion models.
#' One can read in an ms output with readLines().
#' @author Marcelo Gehara
#' @examples # theta = 4Ne x mi x bp
#' Ne = 100000 # effective pop size
#' mi = 1e-8  # per base pair per generation mutation rate
#' bp = 1000  # number of base pairs
#' theta<-4*Ne*mi*bp
#' x<-ms(nsam=10, nrep=1, opts = paste("-t",theta))
#' y<-ms.to.DNAbin(x,bp=1000)
#' nuc.div(y)
#' @noRd
ms.to.DNAbin<-function(ms.output, bp.length){
  ss<-as.numeric(strsplit(ms.output[3]," ")[[1]][2]) # get seg sites

  if(ss>0){
    x<-ms.output[5:length(ms.output)]
    x<-gsub("0","A",x)
    x<-gsub("1","C",x)
  } else {
    x<-vector(mode="character",length=as.numeric(strsplit(ms.output[1]," ")[[1]][2]))
  }

  if(bp.length>0){
    for(i in 1:length(x)){
      x[i]<-paste(x[i],paste(rep("G",(bp.length-ss)),collapse=""),sep="")
    }
  }
  se<-list(NULL,NULL,NULL,NULL)
  names(se)<-c("nb","seq","nam","com")
  se$nb<-length(x) # number of samples
  se$seq<-x # sequencies
  se$nam<-c(1:length(x)) # sequence names, just numbers
  se$com<-NA
  class(se)<-"alignment" # this is alignment
  x<-ape::as.DNAbin(se) # convert to DNAbin
  return(x)
}

#' Convert a pyRAD/ipyrad .alleles file to sequential PHYLIP
#' @description Converts a pyRAD/ipyrad \code{.alleles} file (per-locus
#'              sequence blocks separated by \code{//} lines) into a
#'              multi-locus sequential PHYLIP file that can be loaded
#'              through the PHYLIP slot of \code{main.menu.gui()} or
#'              passed to \code{observed.sumstats(path.to.phylip=...)}.
#'
#'              Each locus is written as a PHYLIP block: a header
#'              \code{<ntax> <nchar>} followed by one tab-delimited
#'              \code{name<TAB>sequence} line per allele. Per-locus
#'              allele counts may vary across loci (samples absent from a
#'              locus are simply skipped — the PipeMaster PHYLIP parser
#'              re-orders by \code{pop.assign} and drops missing samples).
#'              Gaps (\code{-}) and Ns are preserved; \code{observed.sumstats}
#'              masks them from SNP calling.
#'
#'              The \code{//} separator (and any SNP-annotation suffix it
#'              carries) is discarded. The \code{.alleles} format encodes
#'              phased diploid data with sample names ending in \code{_0}
#'              and \code{_1}; population assignment is per allele, so the
#'              pop-assign file passed alongside the resulting PHYLIP must
#'              list both alleles of each individual.
#'
#' @param path.to.alleles Path to the input \code{.alleles} file.
#' @param output Path to the output PHYLIP file.
#' @param verbose Logical. If TRUE prints a summary (default TRUE).
#' @return Invisibly returns the number of loci written.
#' @author Marcelo Gehara
#' @seealso \code{\link{observed.sumstats}}, \code{\link{one.snp.per.locus}}
#' @examples
#' \dontrun{
#' alleles2phylip(
#'   path.to.alleles = "tests/empirical_data/subset4_alleles_filtrado.alleles",
#'   output          = "tests/empirical_data/subset4_alleles_filtrado.phy")
#' }
#' @export
alleles2phylip <- function(path.to.alleles, output, verbose = TRUE) {

  if (!file.exists(path.to.alleles))
    stop("Missing input: ", path.to.alleles)

  t0 <- proc.time()
  if (verbose) cat(sprintf("PipeMaster:: alleles2phylip — reading %s\n",
                            path.to.alleles))
  lines <- readLines(path.to.alleles, warn = FALSE)

  sep_idx <- grep("^//", lines, perl = TRUE)
  n_loci  <- length(sep_idx)
  if (n_loci == 0L)
    stop("No '//' separators found - is this a pyRAD/ipyrad .alleles file?")
  if (verbose) cat(sprintf("PipeMaster::   %d lines read, %d loci detected\n",
                            length(lines), n_loci))

  # Pre-allocate output buffer: one header line per locus + all data lines
  out_buf <- character((length(lines) - n_loci) + n_loci)
  out_pos <- 1L
  prev    <- 0L
  skipped <- 0L

  for (i in seq_len(n_loci)) {
    s <- sep_idx[i]
    block <- lines[(prev + 1L):(s - 1L)]
    block <- block[nzchar(trimws(block))]
    if (length(block) == 0L) { prev <- s; skipped <- skipped + 1L; next }

    # Split on the FIRST whitespace run only (name vs sequence)
    m <- regexpr("[ \t]+", block, perl = TRUE)
    if (any(m < 1L))
      stop("Locus ", i, ": malformed line(s) without name/sequence separator.")
    nm   <- substring(block, 1L, m - 1L)
    rest <- substring(block, m + attr(m, "match.length"))
    seqs <- gsub("[ \t]", "", rest, perl = TRUE)

    L <- nchar(seqs[1L])
    if (any(nchar(seqs) != L))
      stop(sprintf("Locus %d: unequal sequence lengths (%s)",
                   i, paste(unique(nchar(seqs)), collapse = ", ")))

    n_seq <- length(seqs)
    out_buf[out_pos] <- sprintf("%d %d", n_seq, L)
    out_buf[(out_pos + 1L):(out_pos + n_seq)] <- paste(nm, seqs, sep = "\t")
    out_pos <- out_pos + n_seq + 1L
    prev <- s
  }

  out_buf <- out_buf[seq_len(out_pos - 1L)]
  if (verbose) cat(sprintf("PipeMaster::   writing %s (%d lines)\n",
                            output, length(out_buf)))
  writeLines(out_buf, output)

  written <- n_loci - skipped
  if (verbose) cat(sprintf("PipeMaster::   %d loci written (%d empty skipped, %.1f sec)\n",
                            written, skipped, (proc.time() - t0)[3]))
  invisible(written)
}

