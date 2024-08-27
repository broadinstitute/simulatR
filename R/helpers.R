# Inverse CDF of the time of the first denovo mutation
F_denovo_inv <- function(x, mu, N_eff){
  log(1 - log(1-x)/mu) / N_eff
}

# Write VCF
write_vcf <- function(pos, af, ref, alt, id, outdir, sample_dp, sample_sb){

  if(length(pos) == 0){

    # Write empty VCF, if no within-host data
    vcf <- data.frame(
      CHROM = character(0),
      POS = integer(0),
      ID = character(0),
      REF = character(0),
      ALT = character(0),
      QUAL = character(0),
      FILTER = character(0),
      INFO = character(0)
    )

    colnames(vcf)[1] <- "#CHROM"

    write.table(vcf, file = paste0("./", outdir, "/vcf/person_", id, ".vcf"), quote = F, col.names = T, row.names = F)

  }else{

    # Read depths
    dp <- sample_dp(length(ref))

    # Strand biases
    sb <- sample_sb(length(ref))

    alt_reads <- round(dp * af)

    info <- paste0(
      "DP=",
      dp,
      ";AF=",
      format(round(alt_reads / dp, 6), scientific = F),
      ";SB=",
      sb
    )
    vcf <- data.frame(
      CHROM = "ref_genome",
      POS = pos,
      ID = ".",
      REF = ref,
      ALT = alt,
      QUAL = ".",
      FILTER = "PASS",
      INFO = info
    )

    colnames(vcf)[1] <- "#CHROM"

    write.table(vcf, file = paste0("./", outdir, "/vcf/person_", id, ".vcf"), quote = F, col.names = T, row.names = F)

  }
}

# DFS traversal of tree, lowest number first
dfs <- function(h){
  stack <- 1
  explored <- c()
  while (length(stack) > 0) {
    who <- sort(which(h == stack[1]))
    explored <- c(explored, stack[1])
    stack <- stack[-1]
    stack <- c(who, stack)
  }
  return(explored)
}

# BFS traversal of included hosts in a tree tree
bfs <- function(i, h, included){
  out <- i
  frontier <- intersect(which(h == i), included)
  while (length(frontier) > 0) {
    out <- c(out, frontier)
    frontier <- intersect(which(h %in% frontier), included)
  }
  return(out)
}

# Ancestry of a node, leading back to 1
ancestry <- function(i, h){
  if(i == 1){
    return(i)
  }else{
    return(c(ancestry(h[i], h), i))
  }
}

