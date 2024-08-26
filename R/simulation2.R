## New attempt at simulation using TransPhylo

epi_sim_2 <- function(
    a_g = 5,
    lambda_g = 1,
    a_s = 5,
    lambda_s = 1,
    R = 1.5,
    rho = Inf, # Overdispersion parameter. Inf means Poisson distribution.
    mu = 1e-5,
    p = 5e-6,
    v = 1000, # virions produced per replication cycle
    lambda_b = 1.5, # Mean bottleneck size, minus 1. Shifted Poisson distribution assumed.
    init_genome = sample(c("A","C","G","T"), 10000, replace = T),
    sample_dp = function(n){rep(10000, n)},
    sample_sb = function(n){rep(0, n)},
    N = 1e6, # Population size
    p_samp = 0.5, # Probability of sampling
    n_obs = 100, # Number of sampled individuals to simulate
    include_root = TRUE,
    outdir = "my_epidemic",
    seed = NA,
    seed_degrees = FALSE
){

  a_g = 5
  lambda_g = 1
  a_s = 5
  lambda_s = 1
  R = 2
  rho = 2 # Overdispersion parameter. Inf means Poisson distribution.
  mu = 1e-5
  p = 1e-6
  v = exp(1) # virions produced per replication cycle
  lambda_b = 1.01 # Mean bottleneck size, minus 1. Shifted Poisson distribution assumed.
  init_genome = rep("A", 10000)
  sample_dp = function(n){rep(10000, n)}
  sample_sb = function(n){rep(0, n)}
  N = 1e6 # Population size
  p_samp = 0.5 # Probability of sampling
  n_obs = 100 # Number of sampled individuals to simulate
  include_root = FALSE
  outdir = "run-juniper/input-data"
  seed = 20
  seed_degrees = T


  # "p" parameter in the Negative-Binomially distributed offspring distribution
  if(!is.infinite(rho)){
    psi = rho / (R + rho)
  }

  t_max <- 50

  trans <- TransPhylo::simulateOutbreak(
    off.r = rho,
    off.p = 1 - psi,
    pi = p_samp,
    w.shape = a_g,
    w.scale = 1 / lambda_g,
    ws.shape = a_s,
    ws.scale = 1 / lambda_s,
    dateStartOutbreak = 0,
    dateT = t_max,
    nSampled = n_obs
  )

  ttree<-extractTTree(trans)$ttree

  offset <- min(ttree[,1])
  ttree[,1] <- ttree[,1] - offset
  ttree[,2] <- ttree[,2] - offset
  t_max <- t_max - offset


  if(!is.na(seed)){
    set.seed(seed)
  }


  # Probability of a mutation in exponential growth phase, at a given site
  p_growth_mut = 1 - (1-p)^(1/sqrt(p))

  # Length of viral genome
  N_bases = length(init_genome)

  # Set up output directories
  if(dir.exists(outdir)){
    unlink(outdir, recursive = T)
  }
  dir.create(paste0("./", outdir))
  dir.create(paste0("./", outdir, "/vcf"))

  ### Simulation


  root <- which(ttree[,1] == 0)
  # Proportions of particles of each nucleotide at each site IN BOTTLENECK
  bot <- list()
  bot[[root]] <- to_comp(init_genome)

  props <- list()
  props[[root]] <- evolve_exp_growth(bot[[root]], p_growth_mut, p)

  checklist <- which(ttree[,3] == which(ttree[,1] == 0))

  t <- ttree[,1]
  s <- ttree[,2]
  h <- ttree[,3]
  id <- 1:length(t)

  # Write vcf for 1
  if(include_root){
    write_vcf(props[[root]], id[root], outdir, init_genome, sample_dp, sample_sb)
  }



  while (length(checklist) > 0) {

    print(checklist)

    i <- checklist[1]
    checklist <- checklist[-1]

    ## Inherited genotype

    # Evolve this per JC
    delta_t <- t[i] - (t[h[i]] + log(1/sqrt(p)) / (mu / p) / log(v))

    # Duration of exponential growth phase
    g <- log(1/sqrt(p)) / (mu / p) / log(v)

    # Evolution
    bot[[i]] <- evolve_quiescent(props[[h[i]]], bot[[h[i]]], mu, delta_t, lambda_b, g)
    props[[i]] <- evolve_exp_growth(bot[[i]], p_growth_mut, p)

    # If nobody else in the queue is infected by h[i], convert to nucleotides, to save memory
    # if(length(intersect(ttree[checklist, 3], ttree[i, 3])) == 0){
    #   props[[h[i]]] <- to_dna(props[[h[i]]])
    #   bot[[h[i]]] <- character(0)
    # }

    checklist <- c(checklist, which(h == i))

    if(!is.na(ttree[i, 2])){
      # Write vcf for i
      write_vcf(props[[i]], id[i], outdir, init_genome, sample_dp, sample_sb)
    }
  }

  # Cases that we report: sampled cases before (or equal) to t_max
  complete <- which(!is.na(ttree[,2]))

  # Remove root, if necessary
  if(!include_root){
    complete <- setdiff(complete, 1)
  }

  # Convert "props" to FASTA
  for (i in complete) {
    if(typeof(props[[i]]) != "character"){
      props[[i]] <- to_dna(props[[i]])
    }
  }

  names <- paste0("person_", id[complete])

  # Write test date table
  dates <- cbind(names, round(s[complete]))
  write.csv(dates, file = paste0("./", outdir, "/date.csv"), row.names = F, quote = F)

  # Write the true transmission network and the times at which the transmissions occurred
  write.csv(ttree, file = paste0("./", outdir, "/transmission.csv"), row.names = F, quote = F)

  props <- props[complete]
  names(props) <- names

  # Write aligned fasta
  if(length(complete) > 0){
    ape::write.dna(props, file = paste0("./", outdir, "/aligned.fasta"), format = "fasta")
  }else{
    file.create(paste0("./", outdir, "/aligned.fasta"))
  }

  # Write ref genome
  ref <- list(init_genome)
  names(ref) <- "reference_genome"
  ape::write.dna(ref, file = paste0("./", outdir, "/ref.fasta"), format = "fasta")


  plot(trans)






}
