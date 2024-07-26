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
#' @param rho Overdispersion parameter of the offspring distribution. Defaults to Inf, indicating the offspring distribution is Poisson. If finite, the offspring distribution is taken to be Negative Binomial with a mean of R.
#' @param mu Evolution rate, in substitutions per site per day. Defaults to 1e-6.
#' @param p Within-host mutation rate, in mutations per site per replication cycle. Note that mu/p gives equals number of replication cycles per day. Defaults to 5e-6.
#' @param v The number of offspring produced per replication cycle.
#' @param lambda_b The mean number of particles per transmission bottleneck. A shifted Poisson distribution is assumed, i.e. the bottleneck size is 1 + Poisson(lambda_b - 1). Defaults to 1.5.
#' @param init_genome The initial genome present in the seed of the epidemic. Defaults to a genome consisting of 1e4 random draws from A, C, G, T.
#' @param sample_dp A function with one argument, n, that randomly samples n read depths. Defaults to the function that always returns a constant depth of 10,000 reads.
#' @param sample_sb A function with one argument, n, that randomly samples n strand biases. Defaults to the function that always returns a constant strand bias of 0.
#' @param N The size of the population. The population is assumed to be well-mixed, and to start with a single infectious individual. Defaults to 1e6.
#' @param t_max The number of days of epidemic to simulate. Defaults to 50.
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
  t_max = 50, # Number of days to simulate
  include_root = TRUE,
  outdir = "my_epidemic",
  seed = NA,
  seed_degrees = FALSE
){

  # a_g = 5
  # lambda_g = 1
  # a_s = 5
  # lambda_s = 1
  # R = 1.5
  # rho = 1 # Overdispersion parameter. Inf means Poisson distribution.
  # mu = 1e-5
  # p = 1e-6
  # v = exp(1) # virions produced per replication cycle
  # lambda_b = 1.5 # Mean bottleneck size, minus 1. Shifted Poisson distribution assumed.
  # init_genome = rep("A", 10000)
  # sample_dp = function(n){rep(10000, n)}
  # sample_sb = function(n){rep(0, n)}
  # N = 1e6 # Population size
  # t_max = 50 # Number of days to simulate
  # include_root = F
  # outdir = "my_epidemic"
  # seed = 15
  # seed_degrees = T

  if(is.na(seed) & seed_degrees){
    stop("A seed must be specified when seed_degrees = TRUE.")
  }

  if(!is.na(seed)){
    set.seed(seed)
  }

  # "p" parameter in the Negative-Binomially distributed offspring distribution
  if(!is.infinite(rho)){
    psi = rho / (R + rho)
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

  # Times of test
  s <- rgamma(1, a_s, lambda_s)

  # Proportions of particles of each nucleotide at each site IN BOTTLENECK
  bot <- list(to_comp(init_genome))
  props <- list(evolve_exp_growth(bot[[1]], p_growth_mut, p))

  # Generate kids
  if(seed_degrees){
    set.seed(seed + id)
  }
  if(is.infinite(rho)){
    n_kids <- rpois(1, R)
  }else{
    n_kids <- rnbinom(1, rho, psi)
  }

  diagnose <- c()
  diagnose2 <- c()

  if(n_kids > 0){

    # What are their indices?
    who <- (n + 1):(n + n_kids)

    id[who] <- sample(1:N, n_kids, replace = T)
    h[who] <- 1
    t[who] <- t[1] + rgamma(length(who), a_g, lambda_g)
    s[who] <- t[who] + rgamma(length(who), a_s, lambda_s)


    # Add kids to people who need to be accounted for
    checklist <- c(checklist, who)

    # Number of people
    n <- n + length(who)

  }

  # We filled in the data for host 1, so remove it from the checklist
  checklist <- setdiff(checklist, 1)

  # Write vcf for 1
  if(include_root){
    write_vcf(props[[1]], id[1], outdir, init_genome, sample_dp, sample_sb)
  }

  # People to be included in final outbreak output
  complete <- 1

  while (ifelse(length(checklist) > 0, min(t[checklist]), Inf) < t_max) {

    # Next host to update is the one with the earliest infection time
    i <- checklist[which.min(t[checklist])]
    checklist <- setdiff(checklist, i)

    # If already infected previously, nothing to do here
    if(id[i] %in% id[complete]){
      print("browut")
      checklist <- setdiff(checklist, i)
    }else{
      ## Inherited genotype

      # Evolve this per JC
      delta_t <- t[i] - (t[h[i]] + log(1/sqrt(p)) / (mu / p) / log(v))

      # Expected number of mutations X->Y per site per day
      # (1 - (1-p)^(1/sqrt(p))) * mean(rbeta(10000, 1, sample(1:round(1/sqrt(p)), 10000, replace = T))) / (log(1/sqrt(p)) / (mu / p))

      # Duration of exponential growth phase

      # Expected number of mutations X->Y in exponential growth phase under JC:
      # (1/4 - 1/4*exp((-4/3)*mu*log(1/sqrt(p)) / (mu / p) / log(v))) * 10000

      # Duration of exponential growth phase
      g <- log(1/sqrt(p)) / (mu / p) / log(v)

      # Evolution
      bot[[i]] <- evolve_quiescent(props[[h[i]]], bot[[h[i]]], mu, delta_t, lambda_b, g)
      props[[i]] <- evolve_exp_growth(bot[[i]], p_growth_mut, p)

      diagnose[i] <- sum(to_dna(props[[i]]) != to_dna(props[[h[i]]])) / N_bases / (t[i] - t[h[i]])
      diagnose2[i] <- sum(to_dna(props[[i]]) != "A") / N_bases / t[i]

      # Generate kids
      if(seed_degrees){
        set.seed(seed + id[i]) # To ensure seed changes after each iteration
      }
      if(is.infinite(rho)){
        n_kids <- rpois(1, R)
      }else{
        n_kids <- rnbinom(1, rho, psi)
      }

      if(n_kids > 0){

        # What are their indices?
        who <- (n + 1):(n + n_kids)

        id[who] <- sample(1:N, n_kids, replace = T)
        h[who] <- i
        t[who] <- t[i] + rgamma(length(who), a_g, lambda_g)
        s[who] <- t[who] + rgamma(length(who), a_s, lambda_s)


        # Add kids to people who need to be accounted for
        checklist <- c(checklist, who)

        # Number of people in network
        n <- n + length(who)

      }

      # If nobody else in the queue is infected by h[i], convert to nucleotides, to save memory
      if(length(intersect(h[checklist], h[i])) == 0){
        props[[h[i]]] <- to_dna(props[[h[i]]])
        bot[[h[i]]] <- character(0)
      }

      # People in output
      complete <- c(complete, i)

      if(s[i] < t_max){
        # Write vcf for i
        write_vcf(props[[i]], id[i], outdir, init_genome, sample_dp, sample_sb)


      }

      # Report progress
      if(length(complete) %% 10 == 0){
        message(length(complete), " cases generated")
      }
    }
  }



  # Alert the user if the epidemic ends due to no children left
  if(length(checklist) == 0){
    message("The epidemic ended before reaching ", t_max, " days. Consider re-running with a higher basic reproductive number (R) to avoid this behavior.")
  }

  # Remove from "complete" the people who were sampled after t_max
  complete <- setdiff(complete, which(s > t_max))

  # hist(diagnose[complete])
  # mean(diagnose[complete], na.rm = T)
  #
  # hist(diagnose2[complete])
  # mean(diagnose2[complete], na.rm = T)

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
  # Include intermediates who are not in "complete" because they were sampled after t_max
  trans <- data.frame(from = paste0("person_", id[h]), to = paste0("person_", id), time = round(t, digits = 3))
  #trans <- trans[complete, ]
  # if(include_root){
  #   trans <- trans[-1, ]
  # }
  trans <- trans[!is.na(trans$from), ]
  write.csv(trans, file = paste0("./", outdir, "/transmission.csv"), row.names = F, quote = F)

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
}

