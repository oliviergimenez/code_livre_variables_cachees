# function that applies the Viterbi algorithm and calculates prevalence
viterbi_phispdelta <- function(data, fc, k, pi, phi1, phi2, p, delta, n.states){

# build HMM with estimated parameters
hmm_estim <- initHMM(States = c("W", "H", "D"), 
                     Symbols = c("1", "2", "3", "0"), 
                     startProbs = c(pi, 1 - pi, 0), 
                     transProbs = matrix(c(phi1,    0,  1 - phi1,
                                              0, phi2,  1 - phi2,
                                              0,    0,  1), 
                                         nrow = n.states, 
                                         byrow = TRUE),
                     emissionProbs = matrix(c(p * delta,         0, p * (1 - delta), 1 - p,
                                                      0, p * delta, p * (1 - delta), 1 - p,
                                                      0,         0,               0, 1), 
                                            nrow = n.states, 
                                            byrow = TRUE))

# run Viterbi on all encounter histories
viterbi_res <- matrix(NA, nrow(data), ncol(data))
for (i in 1:nrow(data)){
	current_encounter_history <- data[i, fc[i]:k]
	if (length(current_encounter_history) == 1){
		if (current_encounter_history == 1) viterbi_res[i, fc[i]:k] <- 'W'
		if (current_encounter_history == 2) viterbi_res[i, fc[i]:k] <- 'H'
		if (current_encounter_history == 3){
			random_draw = rbinom(1, 1, pi)
			if (random_draw == 1) viterbi_res[i, fc[i]:k] <- 'W'
			if (random_draw == 0) viterbi_res[i, fc[i]:k] <- 'H'
		}	
	}
	if (length(current_encounter_history) > 1){
		# Calculate Viterbi path
		current_obs <- as.character(current_encounter_history) 
		viterbi_res[i, fc[i]:k] <- viterbi(hmm_estim, current_obs)
	}
}

# replace in data the 3s by the reconstituted state
data_filled <- data
for (i in 1:nrow(data_filled)){
	for (j in 1:ncol(data_filled)){
	  if (data_filled[i,j] == 3 & viterbi_res[i,j] == 'W') data_filled[i,j] <- 1
	  if (data_filled[i,j] == 3 & viterbi_res[i,j] == 'H') data_filled[i,j] <- 2
	}
}

estimated_prevalence <- rep(NA, k)
abu_w <- rep(NA, k)
abu_h <- rep(NA, k)
abu_tot <- rep(NA, k)
for (j in 1: k){
  nb_w <- sum(data_filled[,j] == 1, na.rm = T) / p
  abu_w[j] <- nb_w
  nb_h <- sum(data_filled[,j] == 2, na.rm = T) / p
  abu_h[j] <- nb_h
  nb_tot <- sum(data_filled[,j] == 2, na.rm = T) + sum(data_filled[,j] == 1, na.rm = T)
  abu_tot[j] <- nb_tot / p
  estimated_prevalence[j] <- nb_h / (nb_w + nb_h)
}
estimated_prevalence
}
