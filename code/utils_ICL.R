##  ICL computation for a model m adjusted thanks to moveHMM

log_emission_d <- function( m, hid){
  NA_pos <- ( is.na(m$data$step) & is.na(m$data$angle) )
  hid <- hid[ !NA_pos]
  dta <- m$data[!NA_pos,]
  if( sum(dta$step == 0 ) > 0){
    error('ICL not implemented with the zero inflation.\n')
  }
  ## 
  mean_step <- m$mle$stepPar[1, hid]
  sd_step <- m$mle$stepPar[2, hid]
  
  shape = mean_step^2 /sd_step^2
  scale = sd_step^2 / mean_step
  
  
  ##  handling zero  lengths steps
  step_contrib <- rep(0, length(sd_step))
  
  step_contrib <-   sum ( dgamma(dta$step, shape = shape, scale = scale, log = TRUE) * (!is.na(dta$step )))
  
  ## dvonmises is not vectorized over parameter
  mu <- m$mle$anglePar[1, hid]
  kappa <- m$mle$anglePar[2,hid]
  angle <- as.circular(dta$angle, )
  angle_contrib <-  as.numeric(
    sapply(1:length(mu), function(i){
      dvonmises(angle[i], mu = mu[i], kappa[i], log = TRUE)
    })
  )
  angle_contrib <-  sum(angle_contrib, na.rm = TRUE)
  return(step_contrib + angle_contrib)
}  


completeloglike <- function(m,  hid){
  K <- ncol(m$mle$stepPar)
  ID <- unique(m$data$ID)
  beta <- m$mle$beta
  hid_contrib <- lapply(ID, function(id_){
    dta_ind <- m$data %>% filter(ID == id_)
    n_id <- nrow(dta_ind)
    
    var_names <- row.names(m$mle$beta)
    dta_ind$hid_loglike <- rep(NA, n_id)
    dta_ind$hid_loglike[1]  <- 1/K
    for( time in 2:n_id ){
      ## compute transition 
      i <- hid[time-1]
      j <- hid[time]
      col <- (K-1) * (i-1) + (1: (K-1))
     ## avec covariable eta <- exp(beta[1, col] + beta[2, col] * dta_ind$dist.scaled[time-1] + beta[2, col] * dta_ind$dist.scaled2[time-1])
      eta <- eta <- exp(beta[1, col]) 
      if(i == 1){
        eta <- c(1, eta)
      } else if( i == K){
        eta <- c(eta, 1)
      } else{
        eta <- c(eta[1:(i-1)], 1, eta[i:(K-1)])
      }
      eta <- eta/sum(eta)
      ## compute the contribution of the transition to the log likelihood 
      dta_ind[time, 'hid_loglike'] = log(eta[j])
    }
    return(sum(dta_ind$hid_loglike))  
  }
  ) 
  hid_contrib <- Reduce('+', hid_contrib)
  emission_contrib <- log_emission_d(m, hid) 
  
  cat('\n Hidden contribution : ', hid_contrib)
  cat('\n Emission contribution : ', emission_contrib, '\n')
  return( hid_contrib + emission_contrib )
}



completeloglike_depmixS4 <- function(m,  dta, hid){
  K <- depmix_fit@nstates
  ID <- unique(dta$ID)
  n <- nrow(dta)
  hid_contrib <- lapply(ID, function(id_){
    dta_ind <- dta %>% filter(ID == id_)
    n_id <- nrow(dta_ind)
    
    dta_ind$hid_loglike <- rep(NA, n_id)
    dta_ind$hid_loglike[1]  <-depmix_fit@prior@parameters$coefficients[hid[1]]
    trans <- t(depmix_fit@trDens[1,,])
    for( time in 2:n_id ){
      ## compute transition 
      dta_ind[time, 'hid_loglike'] = log(trans[hid[time-1], hid[time]])
    }
    return(sum(dta_ind$hid_loglike))  
  }
  ) 
  hid_contrib <- Reduce('+', hid_contrib)
  emission_contrib <- sum( sapply(1:n, function(i){
    sum(log(m@dens[i,,hid[i]]))
  }), na.rm = TRUE)
  
  
  cat('\n Hidden contribution : ', hid_contrib)
  cat('\n Emission contribution : ', emission_contrib, '\n')
  return( hid_contrib + emission_contrib )
}


ICL_moveHMM <- function(m)
{
  - 2* completeloglike(m, moveHMM::viterbi(m)) + 0.5 * ( 
    length(m$mle$stepPar)+length(m$mle$stepPar)+length(m$mle$beta) + length(m$mle$delta) - 1 ) * log( nrow(m$data)) 
}

ICL_depmixS4 <- function(m, dta){
  - 2* completeloglike_depmixS4(m =m ,dta = dta, hid =  depmixS4::viterbi(m)[,1]) + 0.5 * (m@npars) * log( nrow(dta)) 
}

ICL <- function(m, dta){
  if( 'depmix.fitted' %in% class(m))
    return(ICL_depmixS4(m, dta))
  
  if('moveHMM' %in% class(m) )
    return(ICL_moveHMM(m)) 
  
}