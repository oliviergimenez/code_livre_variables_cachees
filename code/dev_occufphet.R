# function that fits a static occupancy model with false positives and heterogeneity
dev_occufphet <- function(b, data, eff, nh, k) {
  
  # b = vector of parameters on real scale
  # data = site histories
  # eff = nb of sites with that particular history
  # nh = nb sites
  # k = nb surveys
  
  # logit link
  pi <- 1 / (1 + exp(-b[1]))
  psi1 <- 1 / (1 + exp(-b[2]))
  pA10 <- 1 / (1 + exp(-b[3]))
  pA11 <- 1 / (1 + exp(-b[4]))
  pB10 <- 1 / (1 + exp(-b[5]))
  pB11 <- 1 / (1 + exp(-b[6]))
  delta <- 1 / (1 + exp(-b[7]))

  #---------------------------------------------
  # initial states matrices
  PI1 <- c(pi, 1 - pi) # allocation of class A and B
  PI2 <- matrix(c(1 - psi1,        0, psi1, 0,
                         0, 1 - psi1,    0, psi1),
                nrow = 2,
                ncol = 4,
                byrow = T) # occupancy probability in sites of class A and sites of class B
  PI <- PI1 %*% PI2 # initial states matrix
  
  # transition matrix, static model so epsilon = gamma = 0
  A <- diag(4)
  
  # observation matrix
  # first, detections probabilities
  B1 <- matrix(c(
    1 - pA10,    0, pA10, # false positives for sites A
    1 - pB10,    0, pB10, # false positives for sites B
    1 - pA11, pA11,    0, # true positives for sites A
    1 - pB11, pB11,    0),
    nrow = 4,
    ncol = 3,
    byrow = T) # true positives for sites B
  
  # second, classification into certain or uncertain
  B2 <- matrix(c(
    1,     0, 0,
    0, delta, 1 - delta, # probability of classifying a true positive detection into certain
    0,     0, 1),
    nrow = 3,
    ncol = 3,
    byrow = T)
  B <- t(B1 %*% B2) # detection matrix
  
  # calculate deviance 
  garb <- data[1,] # initial states
  l <- 0
  for (i in 1:nh) # loop on sites
  {
    oe <- garb[i] + 1 # first obs
    evennt <- data[,i] + 1 # non-det/det -> 1/2
    ALPHA <- PI * B[oe,]
    for (j in 2:k){ # loop on time
      ALPHA <- (ALPHA %*% A) * B[evennt[j],]
    }
    l <- l + log(sum(ALPHA)) * eff[i]
  }
  l <- - 2 * l
  l
}
