# function that calculates the deviance of model (pi, phi_state, p, delta)
dev_phispdelta_age <- function(b, data, eff, e, garb, nh, km1){

pi <- 1 / (1 + exp(-b[1]))
phi1 <- 1 / (1 + exp(-b[2]))
phi2 <- 1 / (1 + exp(-b[3]))
p  <- 1 / (1 + exp(-b[4]))
delta1 <- 1 / (1 + exp(-b[5]))

# prob of obs (rows) cond on states (col)
B1 <- matrix(c(1-p, p, 0,
               1-p, 0, p,
                 1, 0, 0),
             nrow = 3,
             ncol = 3,
             byrow = T)
delta2 <- 0
B2 <- matrix(c(1,      0,      0,          0,
               0, delta2,      0, 1 - delta2,
               0,      0, delta2, 1 - delta2),
             nrow = 3,
             ncol = 4,
             byrow = T)
B <- t(B1 %*% B2)

# first encounter
BE1 <- matrix(c(0, 1, 0,
                0, 0, 1,
                1, 0, 0),
              nrow = 3,
              ncol = 3,
              byrow = T)
BE2 <- matrix(c(1,      0,      0,          0,
                0, delta1,      0, 1 - delta1,
                0,      0, delta1, 1 - delta1),
              nrow = 3,
              ncol = 4,
              byrow = T)
BE <- t(BE1 %*% BE2) 

# prob of states at t + 1 given states at t
A <- matrix(c(phi1,    0, 1 - phi1,
                 0, phi2, 1 - phi2,
                 0,    0, 1),
            nrow = 3,
            ncol = 3,
            byrow = T)

# init states
PI <- c(pi, 1 - pi, 0)

# likelihood
   l <- 0
   for (i in 1:nh) # loop on ind
   {
      ei <- e[i] # date of first det
      oe <- garb[i] + 1 # init obs
      evennt <- data[,i] + 1 # add 1 to obs to avoid 0s in indexing
      ALPHA <- PI * BE[oe,]
      for (j in (ei + 1):(km1 + 1)) # cond on first capture
      {
      	if ((ei + 1) > (km1 + 1)) {break} # sous MATLAB la commande >> 8:7 rend >> null, alors que sous R, ça rend le vecteur c(8,7)!
        ALPHA <- (ALPHA %*% A) * B[evennt[j],]
      }
      l <- l + log(sum(ALPHA)) * eff[i]
   }
    l <- - 2 * l
    l
}

