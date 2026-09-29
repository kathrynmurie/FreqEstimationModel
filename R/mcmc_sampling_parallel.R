#'@title MCMC
#'@name mcmc_sampling_parallel
#'@description Function to implement mcmc sampler chains in parallel, followed retrospective posterior calculations
#'@param processed_data_list List of processed data
#'@param moi_list List of MOI settings
#'@param frequency_list List of frequency settings
#'@param mcmc_variable_list List of MCMC settings
#'@param cores_max Maximum number of cores to register with doMC
#'@param seed Integer seed for the per-chain RNG streams built by RNGseq().
#'@export
mcmc_sampling_parallel <- function(processed_data_list,    # List of processed data
                                   moi_list,
                                   frequency_list,
                                   mcmc_variable_list,
                                   cores_max,
                                   seed = 1234352)

{
  #==============================================================
  # If not already, install packages for parallising code
  installed_packages <- rownames(installed.packages())
  desired_packages <- c('foreach','doMC','rngtools')
  uninstalled_desired_packages <- desired_packages[!desired_packages %in% installed_packages]
  if(length(uninstalled_desired_packages) > 0){
    install.packages(uninstalled_desired_packages)
  }

  # Load packages if not already
  loaded_packages <- search()
  unloaded_desired_packages <- desired_packages[!paste('package:', desired_packages, sep = '') %in% loaded_packages]
  if(length(unloaded_desired_packages) > 0){
    for(i in 1:length(unloaded_desired_packages)){
      require(unloaded_desired_packages[i], character.only = TRUE)
    }
  }
  #==============================================================

  #===================================================================================================================================
  # Set likelihood function
  if(mcmc_variable_list$NGS & !mcmc_variable_list$log_like_zero){# Assignment for NGS

    log_likelihood_function <- log_likelihood_NGS

  } else { # Assignment for !NGS and QC check for NGS

    log_likelihood_function <- log_likelihood_prevalence
  }
  #===================================================================================================================================

  #===================================================================================================================================
  # Choose log_prior_moi_ratio_function
  if(moi_list$moi_prior == 'Uniform'){
    log_prior_moi_ratio_function <- log_prior_moi_ratio_unif
  }
  if(moi_list$moi_prior == 'Poisson'){
    log_prior_moi_ratio_function <- log_prior_moi_ratio_pois
  }
  if(moi_list$moi_prior == 'Geometric'){
    log_prior_moi_ratio_function <- log_prior_moi_ratio_geom
  }
  if(moi_list$moi_prior == 'nBinomial'){
    log_prior_moi_ratio_function <- log_prior_moi_ratio_nbinom
  }
  #===================================================================================================================================

  #==================================================================================================================================
  # Preallocate stores for MCMC chains
  # Store include burn in (mixed and non_mixed have to be separated first)
  genotype_count_store_chains <- array(0, dim = c(mcmc_variable_list$no_traces_preburnin,
                                                  processed_data_list$no_total,
                                                  processed_data_list$no_haplotypes,
                                                  mcmc_variable_list$no_mcmc_chains),
                                       dimnames = list(NULL,
                                                       processed_data_list$datasampleID,
                                                       processed_data_list$comp_genotypes_character,
                                                       paste('chain', 1:mcmc_variable_list$no_mcmc_chains)))

  genotype_freq_store_chains <- array(NA, dim = c(mcmc_variable_list$no_traces_preburnin,
                                                  processed_data_list$no_haplotypes,
                                                  mcmc_variable_list$no_mcmc_chains),
                                      dimnames = list(NULL,
                                                      processed_data_list$comp_genotypes_character,
                                                      paste('chain', 1:mcmc_variable_list$no_mcmc_chains)))

  log_likelihood_store_chains <- array(NA, dim = c(mcmc_variable_list$no_traces_preburnin,
                                                   processed_data_list$no_total,
                                                   mcmc_variable_list$no_mcmc_chains))

  log_genotype_count_prior_store_chains <- array(NA, dim = c(mcmc_variable_list$no_traces_preburnin,
                                                             processed_data_list$no_total,
                                                             mcmc_variable_list$no_mcmc_chains))

  log_moi_prior_store_chains <- array(NA, dim = c(mcmc_variable_list$no_traces_preburnin,
                                                  processed_data_list$no_total,
                                                  mcmc_variable_list$no_mcmc_chains)) # Not populated till after the mcmc is run

  acceptance_store_chains <- array(NA, dim = c(mcmc_variable_list$no_traces_preburnin,
                                               processed_data_list$no_total,
                                               mcmc_variable_list$no_mcmc_chains))

  run_time_store_chains <- array(NA, dim = c(mcmc_variable_list$no_mcmc_chains,
                                             length(system.time(NULL))),
                                 dimnames = list(paste('chain', 1:mcmc_variable_list$no_mcmc_chains),
                                                 names(system.time(NULL))))
  #==================================================================================================================================

  #==================================================================================================================================
  registerDoMC(cores = min(mcmc_variable_list$no_mcmc_chains, cores_max)) # Register number of cores
  rng <- RNGseq(mcmc_variable_list$no_mcmc_chains, seed) # Pre-specify seed per chain

  # Initialise chains
  parallel_return <- foreach(chain = 1:mcmc_variable_list$no_mcmc_chains) %dopar% {

    rngtools::setRNG(rng[[chain]]) # Launch seed per chain
    #--------------------------------------------------------------------------------------------
    # Initialise the frequency vector
    if(is.null(frequency_list$frequency_initial)){

      frequency_current <- rdirichlet(1, frequency_list$frequency_hyperparameter)

    }else{

      frequency_current <- frequency_list$frequency_initial

    }

    # Initialise MOIs
    if(is.null(moi_list$moi_initial)){

      moi_current <- draw_moi_initial(moi_list,
                                      processed_data_list$datasampleID,
                                      y_no_mxed = processed_data_list$y_no_mxed,
                                      y_mxed = processed_data_list$y_mxed)

    } else {

      moi_current <- moi_list$moi_initial

    }

    # Check point (not inside MCMC chain, hence negligable effect on the time)
    if(length(moi_current) != processed_data_list$no_total | # One moi per blood sample
       any(moi_current%%1 != 0) |                            # Integers
       any(moi_current > moi_list$moi_max) |                          # Less than or equal to maximum
       any(moi_current < 1) |
       any(moi_current[c(processed_data_list$y_mxed, moi_list$moi_prior_min2)] < 2)){
      stop('Illegal initial MOIs')
    }

    # Initialise genotype counts
    genotype_counts_current <- draw_genotype_counts_initial(moi_current,
                                                            processed_data_list,
                                                            frequency_current = if(any(frequency_list$frequency_hyperparameter < 1)){
                                                              rep(1,processed_data_list$no_haplotypes)
                                                            } else {
                                                              frequency_current
                                                            })

    # Check point (not inside MCMC chain, hence negligable effect on the time)
    if(any(rowSums(genotype_counts_current) != moi_current)){
      stop('Illegal initial genotype counts')
    }


    # Initialise the acceptance rate
    acceptance_count <- rep(0, processed_data_list$no_total) # To keep track of M-Hastings acceptance rate, set count to zero for each datasample at the start of each chain

    # Initialise the loglikelihood
    log_likelihood_current <- log_likelihood_function(processed_data_list, genotype_counts = genotype_counts_current)

    # Initialise log_genotype_count_prior
    log_genotype_count_prior_current <- my_log_dmultinom(x = genotype_counts_current, prob = frequency_current, pop_fequency = TRUE)

    # Work out current_genotype_counts_truncated to guide forward removal of a clone and for mh_genotype_counts_truncated
    current_genotype_counts_truncated <- calculate_genotype_counts_truncated(genotype_counts = genotype_counts_current,
                                                                             processed_data_list$comp_genotypes,
                                                                             processed_data_list$raw_data,
                                                                             moi_list$moi_prior_min2)
    #-------------------------------------------------------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------------------------------------------------------
    # Extract stores
    genotype_count_store <- genotype_count_store_chains[,,,chain]
    genotype_freq_store <- genotype_freq_store_chains[,,chain]
    log_likelihood_store <- log_likelihood_store_chains[,,chain]
    log_genotype_count_prior_store <- log_genotype_count_prior_store_chains[,,chain]
    acceptance_store <- acceptance_store_chains[,,chain]
    #-------------------------------------------------------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------------------------------------------------------
    # Store initial values (don't store moi)
    genotype_count_store[1,,] <- genotype_counts_current # Store initial genotype counts (we don't store m as it can be derived from genotype counts)
    genotype_freq_store[1,] <- frequency_current # Store initial frequency
    acceptance_store[1,] <- acceptance_count # (used to access metropolis hastings)
    log_likelihood_store[1,] <- log_likelihood_current # (used in metroposlis hastings update)
    log_genotype_count_prior_store[1,] <- log_genotype_count_prior_current # (used in metroposlis hastings update)
    #-------------------------------------------------------------------------------------------------------------------------------

    #----------------------------------------------------------------------------------------------
    run_time <- system.time(

      for(trace_preburnin in 2:mcmc_variable_list$no_traces_preburnin){

        for(j in 1:mcmc_variable_list$thinning_interval){

          #----------------------------------------------------------------------------------------------------------------
          # Now update frequency and, since the frequency changes, also log_genotype_counts_prior
          frequency_current <- rdirichlet(n = 1, alpha = (colSums(genotype_counts_current) + frequency_list$frequency_hyperparameter))
          #----------------------------------------------------------------------------------------------------------------

          #----------------------------------------------------------------------------------------------------------------
          # Work out which blood samples require clonal addition (because moi_min = {moi_i^t}_min, 2 or 1)
          moi_current_equal_moi_min <- names(which(rowSums(current_genotype_counts_truncated) == 0))

          # Propose new moi
          moi_proposed <- draw_moi_proposed(moi_current, moi_current_equal_moi_min, moi_list$moi_max)

          # Calculate which data samples have an additional clone (+1) and which have one less (-1)
          moi_diff <- moi_proposed - moi_current

          # Check point
          if(any(moi_current > moi_list$moi_max) |
             any(moi_proposed < 1) |
             any(moi_proposed[c(processed_data_list$y_mxed, moi_list$moi_prior_min2)] < 2) |
             any(abs(moi_diff)!=1)){
            stop('Illegal MOIs proposed')
          }

          # Propose new genotype counts
          genotype_counts_proposed <- update_genotype_counts(moi_diff,
                                                             genotype_counts_current,
                                                             current_genotype_counts_truncated,
                                                             processed_data_list)

          # Calculate genotype_counts diff for check point and genotype_count proposal ratio (see below)
          genotype_counts_diff <- genotype_counts_proposed - genotype_counts_current

          # Check point
          if(any(abs(rowSums(genotype_counts_diff))!=1) |
             any(rowSums(genotype_counts_proposed)!= moi_proposed)){
            stop('Illegal genotype update proposed')
          }


          # Calculate the would be genotype counts truncated (changes be moi_diff +1 or -1)
          proposed_genotype_counts_truncated <- calculate_genotype_counts_truncated(genotype_counts = genotype_counts_proposed,
                                                                                    processed_data_list$comp_genotypes,
                                                                                    processed_data_list$raw_data,
                                                                                    moi_list$moi_prior_min2)

          # Work out which proposals are at the lower boundary (either moi_min = {moi_i^t}_min, 2 or 1))
          moi_proposed_equal_moi_min <- names(which(rowSums(proposed_genotype_counts_truncated) == 0))

          # For metropolis hastings, need current boundary for forward step of clone removal and proposed boundary for backward step clone addition
          mh_genotype_counts_truncated <- proposed_genotype_counts_truncated
          mh_genotype_counts_truncated[moi_diff == -1, ] <- current_genotype_counts_truncated[moi_diff == -1, ]

          #----------------------------------------------------------------------------------------------------------------
          # Calculate the target ratio for the MH ratio
          # Note that the prior over the frequency doesn't feature in the target ratio when updating genotype counts and mois
          # Note that log_likelihood_current is an online store
          log_likelihood_proposed <- log_likelihood_function(processed_data_list, genotype_counts = genotype_counts_proposed)
          log_genotype_count_prior_current <- my_log_dmultinom(x = genotype_counts_current, prob = frequency_current, pop_fequency = TRUE)
          log_genotype_count_prior_proposed <- my_log_dmultinom(x = genotype_counts_proposed, prob = frequency_current, pop_fequency = TRUE)
          log_prior_moi_ratio <- log_prior_moi_ratio_function(moi_proposed, moi_current, moi_list$moi_hyperparameter, moi_list$moi_size_hyperparameter)
          log_target_ratio <- (log_likelihood_proposed + log_genotype_count_prior_proposed) - (log_likelihood_current + log_genotype_count_prior_current) + log_prior_moi_ratio

          #----------------------------------------------------------------------------------------------------------------
          # Calculate the moi proposal ratio
          log_moi_proposal_ratio <- calculate_log_moi_proposal_ratio(moi_list$moi_max,
                                                                     moi_current,
                                                                     moi_proposed,
                                                                     moi_current_equal_moi_min,
                                                                     moi_proposed_equal_moi_min,
                                                                     processed_data_list$datasampleID,
                                                                     processed_data_list$no_total)

          #----------------------------------------------------------------------------------------------------------------
          # Calculate the genotype counts proposal ratio
          log_genotype_counts_proposal_ratio <- calculate_log_genotype_counts_proposal_ratio(moi_diff,
                                                                                             genotype_counts_diff,
                                                                                             processed_data_list$neg_log_sum_frequency_truncated,
                                                                                             mh_genotype_counts_truncated,
                                                                                             y_character = processed_data_list$masked_data_character,
                                                                                             log_target_ratio,
                                                                                             processed_data_list$datasampleID,
                                                                                             processed_data_list$no_total)



          #----------------------------------------------------------------------------------------------------------------
          # Calculate acceptance probability
          MHratio <- log_target_ratio + log_genotype_counts_proposal_ratio + log_moi_proposal_ratio
          #----------------------------------------------------------------------------------------------------------------

          #----------------------------------------------------------------------------------------------------------------
          # Evaluate acceptance
          accept <- (log(runif(processed_data_list$no_total)) < MHratio)
          #----------------------------------------------------------------------------------------------------------------

          #----------------------------------------------------------------------------------------------------------------
          # Update latent variables
          genotype_counts_current[accept,] <- genotype_counts_proposed[accept,]
          moi_current[accept] <- moi_proposed[accept] # Accept moi_proposed elements where accept = TRUE

          # Update online store
          log_likelihood_current[accept] <- log_likelihood_proposed[accept] # Update elements where accept = TRUE

          # Update acceptance
          acceptance_count[accept] <- acceptance_count[accept]+1 # To calculate acceptance rate, update elements of count where accept = TRUE

          # Update matrix used to update a_i and m_i
          current_genotype_counts_truncated[accept,] <- proposed_genotype_counts_truncated[accept,] # CHECK boodsampleID order
          #-----------------------------------------------------------------------------------------------------------------
        }

        # Save updated values regardless of change or not
        genotype_count_store[trace_preburnin,,] <- genotype_counts_current
        genotype_freq_store[trace_preburnin,] <- frequency_current
        log_likelihood_store[trace_preburnin,] <- log_likelihood_current
        log_genotype_count_prior_store[trace_preburnin,] <- log_genotype_count_prior_current
        acceptance_store[trace_preburnin,] <- acceptance_count/(trace_preburnin*mcmc_variable_list$thinning_interval)
      }
    )

    chain_return <- list(genotype_count_store = genotype_count_store,
                         genotype_freq_store = genotype_freq_store,
                         log_likelihood_store = log_likelihood_store,
                         log_genotype_count_prior_store = log_genotype_count_prior_store,
                         acceptance_store = acceptance_store,
                         run_time = run_time)

    chain_return
  }

  #-----------------------------------------------------------------------------------------------------------------------------------
  # Update MCMC chain stores
  for(chain in 1:mcmc_variable_list$no_mcmc_chains){
    genotype_count_store_chains[,,,chain] <- parallel_return[[chain]]$genotype_count_store
    genotype_freq_store_chains[,,chain] <- parallel_return[[chain]]$genotype_freq_store
    log_likelihood_store_chains[,,chain] <- parallel_return[[chain]]$log_likelihood_store
    log_genotype_count_prior_store_chains[,,chain] <- parallel_return[[chain]]$log_genotype_count_prior_store
    acceptance_store_chains[,,chain] <- parallel_return[[chain]]$acceptance_store
    run_time_store_chains[chain,] <- parallel_return[[chain]]$run_time
  }
  #-----------------------------------------------------------------------------------------------------------------------------------

  #========================================================================================================================================
  # Calculate the log moi prior
  # First choose function according to moi_list$moi_prior
  if(moi_list$moi_prior=='Uniform'){log_moi_prior_calculation <- log_moi_prior_calculation_unif}
  if(moi_list$moi_prior=='Poisson'){log_moi_prior_calculation <- log_moi_prior_calculation_pois}
  if(moi_list$moi_prior=='Geometric'){log_moi_prior_calculation <- log_moi_prior_calculation_geom}
  if(moi_list$moi_prior=='nBinomial'){log_moi_prior_calculation <- log_moi_prior_calculation_nbinom}

  # Calculate per log moi prior over chains
  for(chain in 1:mcmc_variable_list$no_mcmc_chains){

    log_moi_prior_store <- log_moi_prior_calculation(moi_store = apply(genotype_count_store_chains[,,,chain],c(1,2),sum),
                                                     moi_list,
                                                     y_no_mxed = processed_data_list$y_no_mxed,
                                                     y_mxed = processed_data_list$y_mxed)

    log_moi_prior_store_chains[,,chain] <- log_moi_prior_store

  }
  #========================================================================================================================================


  #========================================================================================================================================
  # Calculate the log posterior
  log_posterior_store_chains <- matrix(NA, nrow = mcmc_variable_list$no_traces_preburnin, ncol = mcmc_variable_list$no_mcmc_chains) # Allocate space

  for(chain in 1:mcmc_variable_list$no_mcmc_chains){

    log_likelihood_store <- log_likelihood_store_chains[,,chain]
    log_genotype_count_prior_store <- log_genotype_count_prior_store_chains[,,chain]
    log_moi_prior_store <- log_moi_prior_store_chains[,,chain]
    log_frequency_prior_store <- log(ddirichlet(genotype_freq_store_chains[,,chain], frequency_list$frequency_hyperparameter)) # prior density on pi

    # Log posterior
    log_posterior_store_chains[,chain] <- apply(log_likelihood_store + log_genotype_count_prior_store + log_moi_prior_store, 1, sum) + log_frequency_prior_store
  }
  #========================================================================================================================================


  #================e========================================================================================================================
  # Return MCMC chain stores and run time
  return(list(genotype_count_store_chains = genotype_count_store_chains,
              genotype_freq_store_chains = genotype_freq_store_chains,
              log_likelihood_store_chains = log_likelihood_store_chains,
              log_posterior_store_chains = log_posterior_store_chains,
              acceptance_store_chains = acceptance_store_chains,
              run_time_store_chains = run_time_store_chains))
  #========================================================================================================================================
}

