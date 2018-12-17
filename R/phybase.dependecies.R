
#### phybase dependesies
read.tree.nodes <- function (str, name = "")
{
  str <- gsub("\\[.*\\]", "", str)
  nobrlens <- 0
  if (length(grep(":", str)) == 0) {
    nobrlens <- 1
    str <- gsub(",", ":1.0,", str)
    str <- gsub(")", ":1.0)", str)
    str <- gsub(";", ":1.0;", str)
  }
  string <- unlist(strsplit(str, NULL))
  leftpar <- which(string == "(")
  rightpar <- which(string == ")")
  if (length(leftpar) != length(leftpar))
    stop("The number of left parenthesis is NOT equal to the number of right  parenthesis")
  speciesname <- sort(species.name(str))
  nspecies <- length(speciesname)
  {
    if (length(leftpar) == (nspecies - 1))
      rooted <- TRUE
    else if (length(leftpar) == (nspecies - 2))
      rooted <- FALSE
    else stop("The number of comma in the tree string is wrong!")
  }
  if (length(name) > 1 & (nspecies != length(name)))
    stop("Wrong number of species names!")
  if (length(name) > 1)
    speciesname <- name
  {
    if (rooted)
      nNodes <- 2 * nspecies - 1
    else nNodes <- 2 * nspecies - 2
    nodes <- matrix(-9, nrow = nNodes, ncol = 7)
  }
  str1 <- str
  if (length(grep("[a-z]", speciesname, ignore.case = TRUE)))
    str1 <- name2node(str1, speciesname)
  father <- nspecies + 1
  while (father < (nNodes + 1)) {
    string <- unlist(strsplit(str1, NULL))
    leftpar <- which(string == "(")
    rightpar <- which(string == ")")
    colon <- which(string == ":")
    {
      if (length(leftpar) == 1)
        substr <- paste(string[leftpar[sum(leftpar <
                                             rightpar[1])]:rightpar[1]], sep = "", collapse = "")
      else substr <- paste(string[leftpar[sum(leftpar <
                                                rightpar[1])]:(colon[which(colon > rightpar[1])[1]] -
                                                                 1)], sep = "", collapse = "")
    }
    substring <- unlist(strsplit(substr, NULL))
    colon <- which(substring == ":")
    comma <- which(substring == ",")
    pound <- which(substring == "#")
    percent <- which(substring == "%")
    combine <- which(substring == "," | substring == ")" |
                       substring == "#" | substring == "%")
    node1 <- as.integer(paste(substring[2:(colon[1] - 1)],
                              sep = "", collapse = ""))
    node2 <- as.integer(paste(substring[(comma[1] + 1):(colon[2] -
                                                          1)], sep = "", collapse = ""))
    if (length(comma) > 1)
      node3 <- as.integer(paste(substring[(comma[2] + 1):(colon[3] -
                                                            1)], sep = "", collapse = ""))
    if (length(colon) == 0) {
      node1Branch <- -9
      node2Branch <- -9
    }
    if (length(colon) > 0) {
      x1 <- combine[sum(combine < colon[1]) + 1] - 1
      x2 <- combine[sum(combine < colon[2]) + 1] - 1
      if (length(colon) == 3) {
        x3 <- combine[sum(combine < colon[3]) + 1] -
          1
        nodes[node3, 4] <- as.double(paste(substring[(colon[3] +
                                                        1):x3], sep = "", collapse = ""))
      }
      node1Branch <- as.double(paste(substring[(colon[1] +
                                                  1):x1], sep = "", collapse = ""))
      node2Branch <- as.double(paste(substring[(colon[2] +
                                                  1):x2], sep = "", collapse = ""))
    }
    if (length(percent) == 0) {
      node1mu <- -9
      node2mu <- -9
    }
    if (length(percent) == 1) {
      if (percent[1] < comma[1]) {
        node1mu <- as.double(paste(substring[(percent[1] +
                                                1):(comma[1] - 1)], sep = "", collapse = ""))
        node2mu <- -9
      }
      else {
        node2mu <- as.double(paste(substring[(percent[1] +
                                                1):(length(substring) - 1)], sep = "", collapse = ""))
        node1mu <- -9
      }
    }
    if (length(percent) == 2) {
      node1mu <- as.double(paste(substring[(percent[1] +
                                              1):(comma[1] - 1)], sep = "", collapse = ""))
      node2mu <- as.double(paste(substring[(percent[2] +
                                              1):(length(substring) - 1)], sep = "", collapse = ""))
    }
    if (length(percent) == 3) {
      node1mu <- as.double(paste(substring[(percent[1] +
                                              1):(comma[1] - 1)], sep = "", collapse = ""))
      node2mu <- as.double(paste(substring[(percent[2] +
                                              1):(comma[2] - 1)], sep = "", collapse = ""))
      node3mu <- as.double(paste(substring[(percent[3] +
                                              1):(length(substring) - 1)], sep = "", collapse = ""))
      nodes[node3, 5] <- node3mu
    }
    if (length(percent) == 0) {
      if (length(pound) == 0) {
        node1theta <- -9
        node2theta <- -9
      }
      if (length(pound) == 1) {
        if (pound[1] < comma[1]) {
          node1theta <- as.double(paste(substring[(pound[1] +
                                                     1):(comma[1] - 1)], sep = "", collapse = ""))
          node2theta <- -9
        }
        else {
          node2theta <- as.double(paste(substring[(pound[1] +
                                                     1):(length(substring) - 1)], sep = "", collapse = ""))
          node1theta <- -9
        }
      }
      if (length(pound) == 2) {
        node1theta <- as.double(paste(substring[(pound[1] +
                                                   1):(comma[1] - 1)], sep = "", collapse = ""))
        node2theta <- as.double(paste(substring[(pound[2] +
                                                   1):(length(substring) - 1)], sep = "", collapse = ""))
      }
      if (length(pound) == 3) {
        node1theta <- as.double(paste(substring[(pound[1] +
                                                   1):(comma[1] - 1)], sep = "", collapse = ""))
        node2theta <- as.double(paste(substring[(pound[2] +
                                                   1):(comma[2] - 1)], sep = "", collapse = ""))
        node3theta <- as.double(paste(substring[(pound[3] +
                                                   1):(length(substring) - 1)], sep = "", collapse = ""))
        nodes[node3, 5] <- node3theta
      }
    }
    if (length(percent) > 0) {
      if (length(pound) == 0) {
        node1theta <- -9
        node2theta <- -9
      }
      if (length(pound) == 1) {
        if (pound[1] < comma[1]) {
          node1theta <- as.double(paste(substring[(pound[1] +
                                                     1):(percent[1] - 1)], sep = "", collapse = ""))
          node2theta <- -9
        }
        else {
          node2theta <- as.double(paste(substring[(pound[1] +
                                                     1):(percent[2] - 1)], sep = "", collapse = ""))
          node1theta <- -9
        }
      }
      if (length(pound) == 2) {
        node1theta <- as.double(paste(substring[(pound[1] +
                                                   1):(percent[1] - 1)], sep = "", collapse = ""))
        node2theta <- as.double(paste(substring[(pound[2] +
                                                   1):(percent[2] - 1)], sep = "", collapse = ""))
      }
      if (length(pound) == 3) {
        node1theta <- as.double(paste(substring[(pound[1] +
                                                   1):(percent[1] - 1)], sep = "", collapse = ""))
        node2theta <- as.double(paste(substring[(pound[2] +
                                                   1):(percent[2] - 1)], sep = "", collapse = ""))
        node3theta <- as.double(paste(substring[(pound[3] +
                                                   1):(percent[3] - 1)], sep = "", collapse = ""))
        nodes[node3, 5] <- node3theta
      }
    }
    nodes[node1, 1] <- father
    nodes[node1, 4] <- node1Branch
    nodes[node1, 5] <- node1theta
    nodes[node1, 6] <- node1mu
    nodes[node2, 1] <- father
    nodes[node2, 4] <- node2Branch
    nodes[node2, 5] <- node2theta
    nodes[node2, 6] <- node2mu
    if (length(comma) > 1) {
      nodes[node3, 1] <- father
      nodes[father, 4] <- node3
    }
    nodes[father, 2] <- node1
    nodes[father, 3] <- node2
    rightpar1 <- which(substring == ")")
    if (rightpar1 < length(substring)) {
      postprob <- paste(substring[(rightpar1 + 1):length(substring)],
                        sep = "", collapse = "")
      nodes[father, 7] <- as.numeric(postprob)
    }
    substr <- gsub("[(]", "[(]", substr)
    substr <- gsub("[)]", "[)]", substr)
    substr <- gsub("\\+", "", substr)
    str1 <- gsub("\\+", "", str1)
    str1 <- gsub(substr, father, str1)
    father <- father + 1
  }
  if (length(grep("%", str1)))
    nodes[nNodes, 6] <- as.double(gsub(";", "", gsub(".*\\%",
                                                     "", str1)))
  if (length(grep("#", str1))) {
    if (length(grep("%", str1)))
      nodes[nNodes, 5] <- as.double(gsub(".*\\#", "", gsub("\\%.*",
                                                           "", str1)))
    else nodes[nNodes, 5] <- as.double(gsub(";", "", gsub(".*\\#",
                                                          "", str1)))
  }
  if (!rooted)
    nodes[nNodes, 1] <- -8
  if (nobrlens == 1)
    nodes[, 4] <- -9
  z <- list(nodes = matrix(0, nNodes, 5), names = "", root = TRUE)
  z$nodes <- nodes
  z$names <- speciesname
  z$root <- rooted
  z
}

#### phybase dependecies
species.name<-function (str)
{
  str <- gsub("[ ]", "", str)
  str <- gsub("[;.]", "", str)
  str <- gsub(":[0-9e-]*", "", str)
  str <- gsub("#[0-9e-]*", "", str)
  str <- gsub(")[0-9e-]*", "", str)
  str <- gsub("[()]", "", str)
  str <- gsub("[ ]", "", str)
  name <- sort(unlist(strsplit(str, split = ",")))
  return(name)
}
