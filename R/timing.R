create_sparsebayes_timing_log <- function(){

  timing <- tibble(NULL)

  overall_mcmc_time <- 0.0
  new_row <- tibble(Name = "Overall MCMC Time",
                    Type = "Overall",
                    Summary = "Yes",
                    Value = 0)
  timing <- rbind(timing, new_row)

  # Mixing Matrix
  new_row <- tibble(Name = "Overall Mixing Matrix Time",
                    Type = "Mixing Matrix",
                    Summary = "Yes",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "YiSi Update",
                    Type = "Mixing Matrix",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "YSi Reduction",
                    Type = "Mixing Matrix",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "mixing matrix partition export",
                    Type = "Mixing Matrix",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "mixing matrix sampler",
                    Type = "Mixing Matrix",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "mixing matrix collection",
                    Type = "Mixing Matrix",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "mixing matrix export",
                    Type = "Mixing Matrix",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)


  # Cluster Membership
  new_row <- tibble(Name = "Overall Cluster Membership Time",
                    Type = "Cluster Membership",
                    Summary = "Yes",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "Stick Breaking Weights Time",
                    Type = "Cluster Membership",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "U Time",
                    Type = "Cluster Membership",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "Cluster Membership Update Time",
                    Type = "Cluster Membership",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "Cluster Cleanup Time",
                    Type = "Cluster Membership",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  # Spatial Maps
  new_row <- tibble(Name = "Overall Spatial Map Time",
                    Type = "Spatial Maps",
                    Summary = "Yes",
                    Value = 0)
  timing <- rbind(timing, new_row)

  # DPM Hyperparameters
  new_row <- tibble(Name = "Overall DPM Hyperparameter Time",
                    Type = "DPM Hyperparameter",
                    Summary = "Yes",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "SuffStats Calculate Time",
                    Type = "DPM Hyperparameter",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "SuffStats Node Summary Time",
                    Type = "DPM Hyperparameter",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "Sample cluster parameters Time",
                    Type = "DPM Hyperparameter",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "Sample concentration time",
                    Type = "DPM Hyperparameter",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  # Variance Terms
  new_row <- tibble(Name = "Overall Variance Term Time",
                    Type = "Variance",
                    Summary = "Yes",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "SSE Time",
                    Type = "Variance",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "Scale Time",
                    Type = "Variance",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "Sigma2 Node Reduction Time",
                    Type = "Variance",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  new_row <- tibble(Name = "Sigma Sq Q",
                    Type = "Variance",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)
  new_row <- tibble(Name = "Local Shrinkage",
                    Type = "Variance",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)
  new_row <- tibble(Name = "Global Shrinkage",
                    Type = "Variance",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)

  # Intermediate Variable Management Terms
  new_row <- tibble(Name = "Overall Variable Management Time",
                    Type = "Variable Management",
                    Summary = "Yes",
                    Value = 0)
  timing <- rbind(timing, new_row)
  new_row <- tibble(Name = "At x Y Time",
                    Type = "Variable Management",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)
  new_row <- tibble(Name = "Eij and Sij Time",
                    Type = "Variable Management",
                    Summary = "No",
                    Value = 0)
  timing <- rbind(timing, new_row)


  # Storage
  new_row <- tibble(Name = "Overall Sample Storage Time",
                    Type = "Storage",
                    Summary = "Yes",
                    Value = 0)
  timing <- rbind(timing, new_row)

  return(timing)

}
