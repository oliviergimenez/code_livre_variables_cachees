
pZXprespast = function(X,param)      #forward algorithm for one time series
{
  T = length(X)
  result = matrix(0,3,T)
  for (t in 1:T)
  {
    if (t == 1) pZXpast = param$pZ0 else pZXpast = t(param$QZ) %*% result[,t-1] 
    result[,t] = pZXpast
    if (X[t] == 1) result[1:2,t] = 0 else result[3,t] = 0
  }
  return(result)
}

pXfuturegivenZ = function(X,param)   #backward algorithm for one time series
{
  T = length(X)
  result = matrix(0,3,T)
  result[,T] = 1
  for (t in T:2)
  {
    PZXgivenOldZ = param$QZ
    if (X[t] == 1) PZXgivenOldZ[,1:2] = 0 else PZXgivenOldZ[,3] = 0
    result[,t-1] = PZXgivenOldZ %*% result[,t]
  }
  return(result)
}

vectForwardBackward = function(X,param) #log-likelihood and lambda coefficients for one time series
{
  T = length(X)
  pZXpp = pZXprespast(X,param)
  pXfgZ = pXfuturegivenZ(X,param)
  pZX = pZXprespast(X,param)*pXfuturegivenZ(X,param)   #P(Zt,whole vector X)
  pZZX = array(rep(param$QZ,T-1),c(3,3,T-1))           #P(Zt,Zt+1,whole vector X)
  for (j in 1:3) pZZX[,j,] = pZZX[,j,]*pZXpp[,1:(T-1)]
  for (i in 1:3) pZZX[i,,] = pZZX[i,,]*pXfgZ[,2:T]
  pZZX[,1:2,X[2:T] == 1] = 0
  pZZX[,3,X[2:T] == 0] = 0
  likelihood = mean(apply(pZX,2,sum))
  result = list()
  result$ll = log(likelihood)
  result$coeffZ = pZX[,1]/likelihood
  result$coeffQ = apply(pZZX,1:2,sum)/likelihood
  return(result)
}

MatrixForwardBackward = function(X,param) #same function extended for a N*T matrix of time series
{
  d = dim(X)
  if (length(d) == 0) return(vectForwardBackward(X,param)) else
  {
    N = d[1]
    result = list()
    llvect = apply(X,1,function(Xvect) vectForwardBackward(Xvect,param)$ll)                          #N-vector,i_th component is the log-likelihood of patch i
    result$ll = sum(llvect)                                                                          #sum on the N patches
    coeffZmat = apply(X,1,function(Xvect) vectForwardBackward(Xvect,param)$coeffZ)                   #3*N-array,i-th column is the lambda(z) vector of patch i
    result$coeffZ = apply(coeffZmat,1,sum)                                                           #sum on the N patches
    coeffQarray = array(apply(X,1,function(Xvect) vectForwardBackward(Xvect,param)$coeffQ),c(3,3,N)) #3*3*N-array,[,,i] is the lambda(z,z') matrix of patch i
    result$coeffQ = apply(coeffQarray,1:2,sum)                                                       #sum on the N patches
    return(result)
  }
}

ForwardBackward = function(X,param)       # X is a list of data vectors for different patches
{
  if (length(dim(X)) == 2) return(MatrixForwardBackward(X,param))
  N = length(X)
  result = list()
  llvect = rep(0,N)
  coeffZmat = matrix(0,3,N)
  coeffQarray = array(0,c(3,3,N))
  for (i in 1:N)
  {
    Xvect = X[[i]]
    llvect[i] = vectForwardBackward(Xvect,param)$ll
    coeffZmat[,i] = vectForwardBackward(Xvect,param)$coeffZ
    coeffQarray[,,i] = vectForwardBackward(Xvect,param)$coeffQ
  }                        
  result$ll = sum(llvect)                                                                                            
  result$coeffZ = apply(coeffZmat,1,sum)                                                           
  result$coeffQ = apply(coeffQarray,1:2,sum)                                                       
  return(result)
}

logLikelihood = function(X,param) ForwardBackward(X,param)$ll


#Maximization functions

p0Maximization = function(coeffZ,p0)
{
  if (p0 == 'free') p0 = sum(coeffZ[2:3])/sum(coeffZ)
  return(p0)
}

gMaximization = function(coeffZ,coeffQ,g)
{
  if (g == 'free') g = (coeffZ[3] + sum(coeffQ[,3]))/(sum(coeffZ[2:3]) + sum(coeffQ[,2:3]))
  return(g)
}

csrMaximization = function(coeffQ,c,s,r)
{
  if (r == 'free') c2 = sum(coeffQ[3,2:3])/sum(coeffQ[3,])      #calculation of c double prime (first step)
  if (s == 'free')                                              #calculation of c prime (first step)
  {
    c1 = sum(coeffQ[2,2:3])/sum(coeffQ[2,])
    if (r == 'free') if(c2 < c1) r = 0
    if (r == 0) c1 = sum(coeffQ[2:3,2:3])/sum(coeffQ[2:3,])  
  }
  if (c == 'free')                                              #estimation of c
  {
    c0 = sum(coeffQ[1,2:3])/sum(coeffQ[1,])                     #calculation of c (first step)
    if (s == 'free') if(c1 < c0) s = 0
    if (s == 0)
    {
      if(r == 0) c0 = sum(coeffQ[,2:3])/sum(coeffQ) else c0 = sum(coeffQ[1:2,2:3])/sum(coeffQ[1:2,])
      c1 = c0
    }
    c = c0
  }
  if (s == 'free') if(c1 < c) s = 0 else s = 1 - (1-c1)/(1-c)   #estimation of s
  c1 = 1- (1-c)*(1-s)
  if (r == 'free') if(c2 < c1) r = 0 else r = 1 - (1-c2)/(1-c1) #estimation of r
  result = list()
  result$c = c
  result$s = s
  result$r = r
  return(result)
}

Maximization = function(coeffZ,coeffQ,model)
{
  result = list()
  result$p0 = p0Maximization(coeffZ,model$p0)
  result$g = gMaximization(coeffZ,coeffQ,model$g)
  csr = csrMaximization(coeffQ,model$c,model$s,model$r)
  result$c = csr$c
  result$s = csr$s
  result$r = csr$r
  return(result)
}


EMiteration = function(X,param,model)
{
  FBresult = ForwardBackward(X,param)
  coeffZ = FBresult$coeffZ
  coeffQ = FBresult$coeffQ
  newParam = Maximization(coeffZ,coeffQ,model)
  EMresult = list()
  EMresult$ll = FBresult$ll
  EMresult$newParam = makeParametersCalculations(newParam)
  return(EMresult)
}

fixToBounds = function(oldParam,param)
{
  newParam = param
  for (parname in c('p0','g','c','s','r'))
  {
    oldPar = oldParam[[parname]]
    newPar = newParam[[parname]]
    if (oldPar > 0.99 & oldPar < newPar) newParam[[parname]] = 1
    if (oldPar < 0.01 & oldPar > newPar) newParam[[parname]] = 0
  }
  return(makeParametersCalculations(newParam))
}

testIdentifiability = function(model)
{
  result = model
  p0 = model$p0
  g = model$g
  c = model$c
  s = model$s
  r = model$r
  alwaysSeedSurvival = FALSE
  if (s == 1 | c == 1 | (g == 1 & r == 1))
  {
    alwaysSeedSurvival = TRUE
    if (r == 'free')
    {
      print('r not identifiable')
      result$r = 1
    }
    if (s == 'free')
    {
      print('s not identifiable')
      result$s = 1
    }
  }
  if (p0 == 1 & alwaysSeedSurvival)
  {
    if (c == 'free') print('c not identifiable')
    result$c = 1
  }
  if (g == 1)
  {
    if (s == 'free' & r == 'free') print('r not identifiable from s')
    result$r = 0
  }
  if (p0 == 'free' & g == 'free' & c == 'free' & s == 'free' & r == 'free')
  {
    #print('too much degrees of liberty, r fixed at 1')
    result$r = 1
  }
  return(result)
}

EMestimation = function(X,nIterations=100,precision = 10^(-6),p0='free',g='free',c='free',s='free',r='free')
{
  model = list()
  model$p0 = p0
  model$g = g
  model$c = c
  model$s = s
  model$r = r
  model = testIdentifiability(model)
  archive = data.frame(matrix(0,nIterations,6))
  names(archive) = c('ll','p0','g','c','s','r')
  param = model
  for (parname in c('p0','g','c','s','r')) if (model[[parname]] == 'free') param[[parname]] = runif(1) 
  param = makeParametersCalculations(param)
  oldLogLikelihood = -Inf
  for (k in 1:nIterations)
  {
    #print(t(param[c('p0','g','c','s','r')]))
    EMresult = EMiteration(X,param,model)
    logLikelihood = EMresult$ll
    archive[k,'ll'] = logLikelihood
    #print(logLikelihood)
    for (parname in c('p0','g','c','s','r')) archive[k,parname] = param[[parname]]
    if(logLikelihood >= oldLogLikelihood & logLikelihood < oldLogLikelihood + precision) break
    #print(k)
    oldLogLikelihood = logLikelihood
    param = fixToBounds(param,EMresult$newParam)
    if (param$c == 1)
    {
      model$c = 1
      model = testIdentifiability(model)
    }
  }
  EMresult = list()
  EMresult$archive = archive[1:k,]
  EMresult$ll = logLikelihood
  EMresult$param = param
  EMresult$nIterations = k - 1
  return(EMresult)
}



