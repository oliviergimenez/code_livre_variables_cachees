
#matrice de transition de Y 
calculateQY = function(param)
{
  result = matrix(0,2,2)
  result[1,1] = 1 - param$c
  result[2,1] = (1 - param$c)*(1 - param$s)*(1 - param$g*param$r)
  result[,2] = 1 - result[,1]
  return(result)
}

#matrice de transition de Zt = (Yt,Xt+1). Les trois états sont dans l'ordre (0,0),(1,0),(1,1).
calculateQZ = function(param)
{
  result = matrix(0,3,3)
  result[1,1] = 1 - param$c
  result[2,1] = (1 - param$c)*(1 - param$s)
  result[3,1] = (1 - param$c)*(1 - param$s)*(1 - param$r)
  result[,2] = (1 - result[,1])*(1 - param$g)
  result[,3] = (1 - result[,1])*param$g
  return(result)
}

#probabilité invariante d'une matrice de transition
calculatePeq = function(Q)
{
  Vmat = eigen(t(Q))$vectors
  #On cherche un vecteur propre dont toutes les composantes sont de même signe et on le renormalise pour que la somme fasse 1
  n = dim(Q)[1]
  for (i in 1:n)
  {
    v = Vmat[,i]
    if ( all(sign(v)>=0) | all(sign(v)<=0) ) return(v/sum(v))
  }
}

# Ajoute à l'objet param ces champs supplémentaires
makeParametersCalculations = function(param)
{
  newparam = param
  newparam$QY = calculateQY(param)
  newparam$QZ = calculateQZ(param)
  #newparam$PZeq = calculatePeq(newparam$QZ)
  #newparam$PYeq = newparam$PZeq[2] + newparam$PZeq[3] #P(Y = 1) à l'équilibre
  newparam$pZ0 = c(1-param$p0, param$p0*(1-param$g), param$p0*param$g)
  return(newparam)
}

makeParametersCalculationsbyC = function(param)
{
  newparam = param
  nCultures = length(param)
  for (k in 1:nCultures) newparam[[k]] = makeParametersCalculations(param[[k]])
  return(newparam)
}

sigmoid = function(x) 1/(1 + exp(-x))

calculateParam = function(hyperParam,lat,lon,vect = FALSE)
{
  result = list()
  for (parname in c('p0','g','c','s'))
  {
    coeffs = hyperParam[,parname]
    linearCombination = coeffs['intercept'] + coeffs['latRatio']*lat + coeffs['lonRatio']*lon
    result[[parname]] = sigmoid(linearCombination)
  }
  result$r = 1
  if(!vect) result = makeParametersCalculations(result)
  return(result)
}

