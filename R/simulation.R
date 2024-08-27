# MIT License
#
# Copyright (c) 2024 Ivan Specht
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#' Simulate an epidemic.
#'
#' This function simulates pathogen spread through a population, given user-specified parameters. It outputs genomic and epidemiological data for each case, including read-level sequencing data.
#'
#' @param a_g Shape parameter of the Gamma-distributed generation interval. Defaults to 5.
#' @param lambda_g Rate parameter of the Gamma-distributed generation interval. Defaults to 1.
#' @param a_s Shape parameter of the Gamma-distributed sojourn interval (time from inoculation to receiving a diagnostic test. Defaults to 5.
#' @param lambda_s Rate parameter of the Gamma-distributed sojourn interval. Defaults to 1.
#' @param R Basic reproductive number. Defaults to 1.5.
#' @param psi Second parameter of the Negative Binomial distribution.
#' @param mu Evolution rate, in substitutions per site per day. Defaults to 1e-6.
#' @param N_eff The within-host effective population size, expressed as the number of virions at time t ie exp(N_eff * t). Defaults to log(1000)
#' @param init_genome The initial genome present in the seed of the epidemic. Defaults to a genome consisting of 1e4 random draws from A, C, G, T.
#' @param sample_dp A function with one argument, n, that randomly samples n read depths. Defaults to the function that always returns a constant depth of 10,000 reads.
#' @param sample_sb A function with one argument, n, that randomly samples n strand biases. Defaults to the function that always returns a constant strand bias of 0.
#' @param min_af Limit of detection for within-host variant frequencies.
#' @param N The size of the population. The population is assumed to be well-mixed, and to start with a single infectious individual. Defaults to 1e6.
#' @param p_samp The probability of being sampled. Defaults to 0.5.
#' @param n_obs The number of sampled cases to generate. Defaults to 100.
#' @param include_root Should the root (i.e. the first case to be seeded in the population) be included in the output FASTA and VCF files? Defaults to TRUE.
#' @param outdir Name of the output directory. Defaults to "my_epidemic."
#' @param seed If integer, the seed is set to that integer. If NA, no seed is set. Defaults to NA.
#' @param seed_degrees If TRUE, the seed is set before each draw of the number of outgoing transmissions per person. Defaults to FALSE. If TRUE, seed must have an integer value.
#' @return A directory containing a single FASTA for all simulated cases, a FASTA consisting of the reference genome, VCF files for each simulated case, and the times of sample collection relative to the inoculation time of the first case (in days).
#' @export

### Simulate outbreak

#library(ape)
#library(LaplacesDemon)

#source("helpers.R")

### Generate transmission network

epi_sim <- function(
  a_g = 5,
  lambda_g = 1,
  a_s = 5,
  lambda_s = 1,
  R = 2,
  psi = 0.5,
  mu = 2e-5,
  N_eff = log(1000), # Within-host N_eff
  init_genome = sample(c("A","C","G","T"), 10000, replace = T),
  sample_dp = function(n){rep(10000, n)},
  sample_sb = function(n){rep(0, n)},
  min_af = 0.03, # Limit of detection for within host variants
  N = 1e6, # Population size
  p_samp = 0.5, # Probability of sampling
  n_obs = 100, # Number of sampled individuals to simulate
  include_root = TRUE,
  outdir = "my_epidemic",
  seed = NA,
  seed_degrees = FALSE
){

  # a_g = 5
  # lambda_g = 1
  # a_s = 5
  # lambda_s = 1
  # R = 2
  # psi = 0.5
  # mu = 2e-5
  # N_eff = log(100)
  # init_genome = sample(c("A","C","G","T"), 10000, replace = T)
  # sample_dp = function(n){rep(10000, n)}
  # sample_sb = function(n){rep(0, n)}
  # N = 1e6
  # p_samp = 0.5
  # n_obs = 100
  # include_root = FALSE
  # outdir = "my_epidemic"
  # seed = 6
  # seed_degrees = TRUE

  if(is.na(seed) & seed_degrees){
    stop("A seed must be specified when seed_degrees = TRUE.")
  }

  if(!is.na(seed)){
    set.seed(seed)
  }

  rho <- R * psi / (1 - psi)

  # Length of viral genome
  N_bases = length(init_genome)

  # Set up output directories
  if(dir.exists(outdir)){
    unlink(outdir, recursive = T)
  }
  dir.create(paste0("./", outdir))
  dir.create(paste0("./", outdir, "/vcf"))

  ### Simulation

  # Hosts whose data needs to be updated
  checklist <- 1

  # Number of hosts
  n <- 1

  # Vector of ancestors
  h <- NA

  # Vector of IDs (element of population)
  # Root case is (WLOG) case 1
  id <- 1

  # Times of infection
  t <- 0

  # Time of test
  s <- rgamma(1, a_s, lambda_s)

  sampled <- runif(1) < p_samp

  # If not including root, we need to generate an extra case, since this one's getting eliminated
  if(!include_root & sampled){
    n_obs <- n_obs + 1
  }

  # True reproductive number
  true_Rs <- c()

  # Generate kids
  if(seed_degrees){
    set.seed(seed + id)
  }
  if(is.infinite(rho)){
    n_kids <- rpois(1, R)
  }else{
    n_kids <- rnbinom(1, rho, psi)
  }
  true_Rs <- c(true_Rs, n_kids)

  # Initialize maximum time of epidemic to Inf. We will revise this as we come close to the target number of sampled cases
  t_max <- Inf


  if(n_kids > 0){

    # What are their indices?
    who <- (n + 1):(n + n_kids)

    id[who] <- sample(1:N, n_kids, replace = T)
    h[who] <- 1
    t[who] <- t[1] + rgamma(length(who), a_g, lambda_g)
    s[who] <- t[who] + rgamma(length(who), a_s, lambda_s)
    sampled[who] <- runif(length(who)) < p_samp

    # Add kids to people who need to be accounted for
    checklist <- c(checklist, who)

    # Number of people
    n <- n + length(who)

  }

  # We filled in the data for host 1, so remove it from the checklist
  checklist <- setdiff(checklist, 1)

  # Number of cases generated thus far
  n_gen <- 0

  while (min(t[checklist]) < ifelse(length(s[sampled]) >= n_obs, sort(s[sampled])[n_obs], Inf)) {

    # Next host to update is the one with the earliest infection time
    i <- checklist[which.min(t[checklist])]
    checklist <- setdiff(checklist, i)

    # If already infected previously, nothing to do here
    if(id[i] %in% id[which(t < t[i])]){
      #print("browut")
      checklist <- setdiff(checklist, i)
    }else{
      ## Inherited genotype

      # Generate kids
      if(seed_degrees){
        set.seed(seed + id[i]) # To ensure seed changes after each iteration
      }
      if(is.infinite(rho)){
        n_kids <- rpois(1, R)
      }else{
        n_kids <- rnbinom(1, rho, psi)
      }

      true_Rs <- c(true_Rs, n_kids)

      if(n_kids > 0){

        # What are their indices?
        who <- (n + 1):(n + n_kids)

        id[who] <- sample(1:N, n_kids, replace = T)
        h[who] <- i
        t[who] <- t[i] + rgamma(length(who), a_g, lambda_g)
        s[who] <- t[who] + rgamma(length(who), a_s, lambda_s)
        sampled[who] <- runif(length(who)) < p_samp

        # Add kids to people who need to be accounted for
        checklist <- c(checklist, who)

        # Number of people in network
        n <- n + length(who)

      }

      n_gen <- n_gen + 1

      # Report progress
      if(n_gen %% 10 == 0){
        message(n_gen, " cases generated")
      }
    }
  }

  t_max <- sort(s[sampled])[n_obs]

  # Alert the user if the epidemic ends due to no children left
  if(length(checklist) == 0){
    message("The epidemic ended before reaching ", n_obs, " people. Consider re-running with a higher basic reproductive number (R) to avoid this behavior.")
  }

  # Cases that we report: sampled cases before (or equal) to t_max
  complete <- which(sampled & s <= t_max)

  print(paste("With perfect sampling,", sum(s <= t_max), "cases would have been sampled by time t_max."))
  print(paste("Realized reproductive number:", mean(true_Rs)))

  # Remove root, if necessary
  if(!include_root){
    complete <- setdiff(complete, 1)
  }

  names <- paste0("person_", id[complete])

  # Write test date table
  dates <- cbind(names, round(s[complete]))
  write.csv(dates, file = paste0("./", outdir, "/date.csv"), row.names = F, quote = F)

  # Write the true transmission network and the times at which the transmissions occurred
  # Include intermediates who are not in "complete" because they were sampled after t_max
  trans <- data.frame(from = paste0("person_", id[h]), to = paste0("person_", id), time = round(t, digits = 3))

  trans <- trans[!is.na(trans$from), ]
  write.csv(trans, file = paste0("./", outdir, "/transmission.csv"), row.names = F, quote = F)

  ## Plotting

  # Figure out who's included
  included <- 1
  for (i in complete) {
    anc <- i
    while(!(anc %in% included)){
      included <- c(included, anc)
      anc <- h[anc]
    }
  }

  # BFS order of included hosts on tree
  ord <- bfs(1, h, included)

  # Mutations
  mut_times <- list()
  mut_sites <- list()
  mut_from <- list()
  mut_to <- list()

  mut_times[[1]] <- numeric(0)
  mut_sites[[1]] <- integer(0)
  mut_from[[1]] <- character(0)
  mut_to[[1]] <- character(0)

  # Record total evolutionary time and number of mutations
  tot_evo_time <- 0
  tot_n_mut <- 0

  for (i in ord[2:length(ord)]) {

    # Time difference, bottleneck to bottleneck
    delta_t <- t[i] - t[h[i]]
    tot_evo_time <- tot_evo_time + delta_t

    # Number of mutations
    n_mut <- rpois(1, mu * delta_t * length(init_genome))
    tot_n_mut <- tot_n_mut + n_mut

    if(n_mut == 0){
      mut_times[[i]] <- numeric(0)
      mut_sites[[i]] <- integer(0)
      mut_from[[i]] <- character(0)
      mut_to[[i]] <- character(0)
    }else{
      # Times at which they occur
      mut_times[[i]] <- sort(runif(n_mut, t[h[i]], t[i]))

      # Sites at which they occur
      mut_sites[[i]] <- sample(1:length(init_genome), n_mut, replace = T)

      # For from and to, need to get current nucleotide at that position
      anc <- ancestry(h[i], h)

      mut_from[[i]] <- character(0)
      mut_to[[i]] <- character(0)

      for (j in 1:n_mut) {
        prev_sites <- unlist(mut_sites[anc])
        if(j > 1){
          prev_sites <- c(prev_sites, mut_sites[[i]][1:(j-1)])
        }
        if(mut_sites[[i]][j] %in% prev_sites){
          # Index of latest occurence of the site mutating
          ind_latest <- max(which(prev_sites == mut_sites[[i]][j]))
          prev_to <- unlist(mut_to[anc])
          if(j > 1){
            prev_to <- c(prev_to, mut_to[[i]][1:(j-1)])
          }
          mut_from[[i]][j] <- prev_to[ind_latest]
          mut_to[[i]][j] <- sample(setdiff(c("A", "C", "G", "T"), prev_to[ind_latest]), 1)
        }else{
          mut_from[[i]][j] <- init_genome[mut_sites[[i]][j]]
          mut_to[[i]][j] <- sample(setdiff(c("A", "C", "G", "T"), init_genome[mut_sites[[i]][j]]), 1)
        }
      }
    }
  }

  print(paste("Total evolutionary time:", round(tot_evo_time), "days"))
  print(paste("Total number of mutations:", tot_n_mut))

  cons <- list()

  ## Now, generate within-host variants. Only needed now for cases we output
  for (i in complete) {
    # Time of first mutation at each site, after the bottleneck
    t_1st_mut <- F_denovo_inv(runif(length(init_genome)), mu, N_eff)

    # But some of these positions will have already mutated based on what we drew in the previous step
    # Children of i
    js <- intersect(which(h == i), included)

    # Sites at which mutations occur in these js
    sites <- unlist(mut_sites[js])
    # Times at which they occur past the bottleneck
    ts <- unlist(mut_times[js]) - t[i]
    # What nucleotide we mutate to
    tos <- unlist(mut_to[js])

    # Which sites have detected iSNVs that get passed on?
    trans_isnv_sites <- c()

    # And what did said sites mutate into?
    trans_isnv_tos <- c()

    if(length(sites) > 0){
      for (s in unique(sites)) {
        # First time of a mutation at site s that gets transmitted onwards

        # This is a vector of all times we get a mutation at site s in host i that gets passed on
        t_mut <- ts[which(sites == s)]

        # The one recorded is the minimum of t_mut, which occurs at which.min(t_mut)
        # The "to" allele is tos[which(sites== s)], the "to" allele for sites matching s
        to <- (tos[which(sites == s)])[which.min(t_mut)]

        # We only care about the first time of the mutation past the bottleneck
        t_mut <- min(t_mut)

        # If it's before what we already sampled, update t_1st_mut, and record this fact
        if(t_mut < t_1st_mut[s]){
          # print(t_mut)
          # print(t_1st_mut[s])
          t_1st_mut[s] <- t_mut
          trans_isnv_sites <- c(trans_isnv_sites, s)
          trans_isnv_tos <- c(trans_isnv_tos, to)
        }
      }
    }

    # Proportion of the viral population that's mutated
    prop_mut <- 1 / exp(N_eff * t_1st_mut)

    # Which of these sites have proportions that are above LOD?
    # This step just saves time: we would never record proporitons under LOD anyway
    above_lod <- which(prop_mut > min_af)

    ## Get genome at bottleneck

    # Ancestry
    anc <- ancestry(i, h)

    # Sites that mutate (can replace old vectors called sites and tos; now unneccary)
    sites <- unlist(mut_sites[anc])
    tos <- unlist(mut_to[anc])
    bot <- init_genome

    if(length(sites) > 0){
      for (s in 1:length(sites)) {
        bot[sites[s]] <- tos[s]
      }
    }

    # Loop over sites with detected iSNVs
    isnv_from <- c()
    isnv_pos <- c()
    isnv_to <- c()
    isnv_af <- c()
    for (s in above_lod) {
      # If s is a site in transmitted iSNVs, we already know the "to" allele
      if(s %in% trans_isnv_sites){
        print("Hooray! There's iSNV evidence for a transmission!")
        to <- trans_isnv_tos[which(trans_isnv_sites == s)]
      }else{
        to <- sample(setdiff(c("A", "C", "G", "T"), bot[s]), 1)
      }


      # If the "to" is the same as the root sequence, from <- to, to <- bot
      if(to == init_genome[s]){
        to <- bot[s]
        af <- 1 - prop_mut[s]
      }else{
        af <- prop_mut[s]
      }
      from <- init_genome[s]

      # If iSNV frequency within LOD, report it
      if(af > min_af & af < 1 - min_af){
        isnv_to <- c(isnv_to, to)
        isnv_from <- c(isnv_from, from)
        isnv_pos <- c(isnv_pos, s)
        isnv_af <- c(isnv_af, af)
      }

      # If necessary, update bottleneck to consensus sequence
      if(af > 0.5){
        bot[s] <- to
      }
    }

    # Write VCF
    write_vcf(isnv_pos, isnv_af, isnv_from, isnv_to, id[i], outdir, sample_dp, sample_sb)

    # Append consensus sequence
    cons[[i]] <- bot
  }

  cons_complete <- cons[complete]
  names(cons_complete) <- paste0("person_", id[complete])

  # Write aligned fasta
  if(length(complete) > 0){
    ape::write.dna(cons_complete, file = paste0("./", outdir, "/aligned.fasta"), format = "fasta")
  }else{
    file.create(paste0("./", outdir, "/aligned.fasta"))
  }

  # Write ref genome
  ref <- list(init_genome)
  names(ref) <- "reference_genome"
  ape::write.dna(ref, file = paste0("./", outdir, "/ref.fasta"), format = "fasta")

  print(paste("Number of included hosts:", length(included)))

  h_comp <- match(h, included)[included]
  t_comp <- t[included]
  n <- length(h_comp)
  ord <- rev(dfs(h_comp))

  leaves <- which(!(1:n %in% h_comp))
  n_leaves <- length(leaves)

  # Angle of each node
  thetas <- c()
  leaf_count <- 0

  for (i in ord) {
    if(!(i %in% leaves)){
      kids <- which(h_comp == i)
      thetas[i] <- mean(thetas[kids])
    }else{
      thetas[i] <- leaf_count * 2 * pi / n_leaves
      leaf_count <- leaf_count + 1
    }
  }

  df_standard <- data.frame(x =t_comp, y = thetas)


  # vertical segments
  xs <- c()
  ystart <- c()
  yend <- c()
  if(include_root){
    who <- 1:n
  }else{
    who <- 2:n
    df_standard <- df_standard[-1, ]
  }
  for (i in who) {
    kids <- which(h_comp == i)
    if(length(kids) > 0){
      xs <- c(xs, t_comp[i])
      ystart <- c(ystart, min(thetas[kids]))
      yend <- c(yend, max(thetas[kids]))
    }
  }

  colors <- rep('black', length(who))
  colors[!(included[who] %in% complete)] <- 'gray'

  big <- ggplot2::ggplot() +
    ggplot2::geom_segment(mapping = ggplot2::aes(x = (t_comp[h_comp])[-1], xend = t_comp[-1], y = thetas[-1], yend = thetas[-1]), linewidth = 0.5) +
    ggplot2::geom_segment(mapping = ggplot2::aes(x = xs, xend = xs, y = ystart, yend = yend), linewidth = 0.5) +
    ggplot2::geom_point(mapping = ggplot2::aes(x = df_standard$x, y = df_standard$y, color = colors), size = 1) +
    ggplot2::xlab("Evolutionary Time (days)") +
    ggplot2::scale_y_continuous(breaks = NULL) +
    ggplot2::scale_color_manual(values = c("black", "gray")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      legend.position='none'
    )
  print(big)

}

