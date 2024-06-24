# Evolve genome through exponential growth phase, starting from bottleneck
evolve_site <- function(site, mu, delta_t, b, p_growth_mut, p){

  # Evolve per JC model
  site <- 1/4 + (site - 1/4)*exp((-4/3)*mu*delta_t)

  # Sample bottleneck
  site <- rmultinom(1, b, site)[,1]

  # Does a mutation occur in the expo growth phase?
  denovo <- runif(1) < p_growth_mut
  if(denovo){
    # Proportion of mutated particles
    prop_mut <- rbeta(1, 1, sample(1:round(1/sqrt(p)), 1))

    # From which bottleneck particle is the mutation inherited?
    from <- sample(1:4, 1, prob = site)

    # And, to which nucleotide have we mutated?
    to <- sample((1:4)[-from], 1)
  }

  # Proportions of each allele, not yet accounting for de novo mutations
  if(sum(site == 0) == 3){
    site <- site / b
  }else{
    site <- LaplacesDemon::rdirichlet(1, site)[1, ]
  }



  if(denovo){
    # Change in composition of viral population
    change <- site[from] * prop_mut
    site[from] <- site[from] - change
    site[to] <- site[to] + change
  }

  site
}


# Evolve genome via Jukes-Cantor, then whatever happens in the exponential growth phase
evolve <- function(comp, mu, delta_t, lambda_b, p_growth_mut, p){

  # Bottleneck size
  b <- rpois(1, lambda_b - 1) + 1

  lapply(comp, evolve_site, mu=mu, delta_t=delta_t, b=b, p_growth_mut=p_growth_mut, p=p)

}

# Convert genome list to nucleotides
to_dna <- function(comp){
  c("A","C","G","T")[sapply(comp, which.max)]
}

# Convert vector of nucleotides to nucleotide composition at each site
to_comp <- function(genome){
  lapply(1:length(genome), function(i){as.numeric(c("A","C","G","T") == genome[i])})
}

# Write VCF
write_vcf <- function(comp, id, outdir, init_genome, sample_dp, sample_sb){
  pos <- c()
  af <- c()
  ref <- c()
  alt <- c()

  # Assuming that if no iSNV at a site at end of expo growth phase, no iSNV detected
  for (i in 1:length(comp)) {
    if(sum(comp[[i]] > 0) > 1){
      which_ref <- which(c("A","C","G","T") == init_genome[i])
      which_alt <- setdiff(which(comp[[i]] > 0), which_ref)
      pos <- c(pos, rep(i, length(which_alt)))
      af <- c(af, comp[[i]][which_alt])
      ref <- c(ref, rep(which_ref, length(which_alt)))
      alt <- c(alt, which_alt)
    }

  }

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

    write.table(vcf, file = paste0("./", outdir, "/vcf/person_", id, ".vcf"), quote = F, col.names = T, row.names = F)

  }else{
    ref <- c("A","C","G","T")[ref]
    alt <- c("A","C","G","T")[alt]

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

    write.table(vcf, file = paste0("./", outdir, "/vcf/person_", id, ".vcf"), quote = F, col.names = T, row.names = F)

  }
}


