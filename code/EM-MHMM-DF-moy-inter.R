
library('optimr')
triage= function(c,Y){
  o=dim(Y)
  if(is.null(o)){C=length(Y)
  D= Y
  D[c]=0
  temp= order(D)[2:C]
  res=D[temp]
  }else{  N=o[2]
  C=o[1]
  D=Y
  D[c,]=0
  res=matrix(rep(0,(C-1)*N),(C-1))
  for(i in 1:N){
    temp= order(D[,i])[2:C]
    res[,i]=D[temp,i]
  }
  }
  
  
  
  
  return(res)
  
}

comb = function(n, x) {
  return(factorial(n) / (factorial(x) * factorial(n-x)))
}

## attention le 2 désigne le 0-4 pour le modele interchangeable et le 2 pour moy 1-5 

logitphi = function(beta,K,k,i,I){
  test=comb(K-1,k-1)
  
  test1= exp(-(beta[1]*((i-1)/I)+beta[2]))
  test3=(1/(1+test1))^(K-1) * test1^(K-k)
  result= test* test3
  return(result)
}

funcphi = function(beta,I,K){
  phi=array(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),dim=c(I,K))
  
  
  
  for(k in 1:K){
    for(l in 1:I){
      E=exp(-(beta[1]*(l-1)/I+beta[2]))
      phi[l,k]= comb(K-1,k-1)* (1/(1+E))^(K-1) * E^(K-k)
    }
  }
  
  
  return(phi)
}

funcpi = function(gamma,I){
  pi=rep(0,I)
  
  
  
  
  for(l in 1:I){
    E=exp(-(gamma))
    pi[l]= comb(I-1,l-1)* (1/(1+E))^(I-1) * E^(I-l)
  }
  
  
  return(pi)
}

ZIphi = function(beta,I,K){
  phi=array(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),dim=c(I,K))
  
  if(I==2){
    phi[1,1]=1
  phi[1,2:K]=rep(0,K-1)
    for(k in 1:K){
      l=2
        E=exp(-(beta[1]))
        phi[l,k]= comb(K-1,k-1)* (1/(1+E))^(K-1) * E^(K-k)
      }
    }  else{
  
  for(k in 1:K){
    for(l in 1:I){
      if(l==1){phi[1,1]=1
      phi[1,2:K]=rep(0,K-1)} else{
      E=exp(-(beta[1]*(l-1)/I+beta[2]))
      phi[l,k]= comb(K-1,k-1)* (1/(1+E))^(K-1) * E^(K-k)}
    }
  }
  }
  
  return(phi)
}

logitphi_2 = function(beta,K,k,i,I){
  test=comb(K-1,k-1)
  
  test1= exp(-(beta[1]*(i/I)+beta[2]))
  test3=(1/(1+test1))^(K-1) * test1^(K-k)
  result= test* test3
  return(result)
}
funcphi_2 = function(beta,I,K){
  phi=array(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),dim=c(I,K))
  
  
  
  for(k in 1:K){
    for(l in 1:I){
      E=exp(-(beta[1]*(l)/I+beta[2]))
      phi[l,k]= comb(K-1,k-1)* (1/(1+E))^(K-1) * E^(K-k)
    }
  }
  
  
  return(phi)
}
#1-5
logitA= function(alpha,i,ii,I,K,Y,c,n){
  C=  length(Y[,n])
  
  test=comb(I-1,ii-1)
  
  test1= exp(-(alpha[1]*(i/I)+alpha[2]*(Y[c,n]/K) + alpha[3]*(bij_inter(c,n,Y,K)/bij_inter(c,n,rep(K,C),K)) + alpha[4]))
  test3=(1/(1+test1))^(I-1) * test1^(I-ii)
  result= test* test3
  return(result)
}
#fonction implémenté en premier ! mais constante varie en fonction des autres parmetres 1-5
funcA = function(alpha,I,K,C){
  tot=(bij_inter(1,1,rep(K,C),K))
  A= array(rep(0,I*I*K*tot),dim=c(I,I,K*tot))
  
  
  for(i in 1:tot){
    for(k in 1:K){
      for(j in 1:I){
        for(l in 1:I){
          E=exp(-(alpha[1]*(j)/I+alpha[2]*k/K+ alpha[3]*(i)/tot+alpha[4]))
          A[j,l,(i+((k-1)*tot))]= comb(I-1,l-1)* (1/(1+E))^(I-1) * E^(I-l)
        }
      }
    }
    
  }
  return(A)
}
#0-4
funcA_2 = function(alpha,I,K,C){
  tot=(bij_inter(1,1,rep(K,C),K))
  A= array(rep(0,I*I*K*tot),dim=c(I,I,K*tot))
  
  
  for(i in 1:tot){
    for(k in 1:K){
      for(j in 1:I){
        for(l in 1:I){
          E=exp(-(alpha[1]*(j-1)/I+alpha[2]*(k-1)/K+ alpha[3]*(i-1)/tot+alpha[4]))
          A[j,l,(i+((k-1)*tot))]= comb(I-1,l-1)* (1/(1+E))^(I-1) * E^(I-l)
        }
      }
    }
    
  }
  return(A)
}


bij_inter = function(c,n,Y,K){
  
  DD= dim(Y)
  if(is.null(DD)||DD[2]==1){C= length(Y)
  YY=Y
  YY[c]=1
  if(C==2){res = c+1
  res= res%%C
  if(res==0){return(Y[C])} else{ return(Y[1])}}
  pourtri=YY
  tri1=  triage(c,pourtri)} else{
    N=DD[2]
    C=DD[1]
    if(C==2){res = c+1
    res= res%%C
    if(res==0){return(Y[C,n])} else{ return(Y[1,n])}}
    YY=Y
    YY[c,n]=1
    pourtri=YY[,n]
    tri1=  triage(c,pourtri)}
  
  
  temp=0
  summ=0
  i_c=0
  for(i in (length(tri1)):2){
    temp=temp+1
    if(tri1[temp]!=1){
      
      if(i ==length(tri1)){
        
        
        for(j in K:(K-tri1[1]+2) ){
          
          
          
          summ= comb(j+i-2,i-1)+ summ
        }
      } else { if(tri1[temp] !=tri1[temp-1] ){
        for(j in K:(K-tri1[temp]+2+tri1[temp-1]-1) ){
          
          
          
          summ= comb(j+i-1-tri1[temp-1],i-1)+ summ
        }
      }
        
      }
    }
  }
  
  summ = summ + tri1[length(tri1)]-tri1[length(tri1)-1]+1
  
  
  return(summ)
}

bij_moy = function(c,n,Y,K,C){
  
  D=dim(Y)
  if( is.null(D)) {lis=sum(Y)-Y[c]}
  else{lis=sum(Y[1:C,n])- Y[c,n]}
  result= round(lis/(C-1))
  
 
  etatV=floor(result)
  return(etatV)
}

bij_dist = function(c,n,Y,K,C,dist){
  
  D=dim(Y)
  voisin=which(dist[c,]<70)[which(which(dist[c,]<70)!=c)]
  if( is.null(D)) {
    lis=sum(Y[voisin])}
  else{lis=sum(Y[voisin,n])}
  result= round(lis/length(voisin))
  
  
  
  return(result)
}

field_to_vect_moy = function(c,Y,K){
  D=dim(Y)
  if( is.null(D)) {N=1
  C= length(Y)
  
  
  temp = Y[c]-1
  res= rep(0,N)
  
  
  
  
  for(n in 1:N){
    res[n]=bij_moy(c,n,Y,K,C)
  }
  
  res = temp*K + res
  
  }else{
    C=D[1]
    N= D[2]
    
    temp = Y[c,]-1
    res= rep(0,N)
    
    for(n in 1:N){
      res[n]=bij_moy(c,n,Y,K,C)
    }
    
    res = temp*K + res }
  
  
  return(res)
}

field_to_vect_dist = function(c,Y,K,dist){
  D=dim(Y)
  if( is.null(D)) {N=1
  C= length(Y)
  
  
  temp = Y[c]-1
  res= rep(0,N)
  
  voisin=which(dist[c,]<70)[which(which(dist[c,]<70)!=c)]
  
  
  for(n in 1:N){
    res[n]=round(sum(Y[voisin])/length(voisin))
  
  
  }
  
  res = temp*K + res
  
  }else{
    C=D[1]
    N= D[2]
    voisin=which(dist[c,]<70)[which(which(dist[c,]<70)!=c)]
    temp = Y[c,]-1
    res= rep(0,N)
    
    
    
    for(n in 1:N){
      res[n]=round(sum(Y[voisin,n])/length(voisin))
     
    }
  
    
    res = temp*K + res }
  
  
  return(res)
}

#0-4
logitA_moy= function(alpha,i,ii,I,K,Y,c,n){
  C=  length(Y[,n])
  
  test=comb(I-1,ii-1)
  
  test1= exp(-(alpha[1]*((i-1)/I)+alpha[2]*((Y[c,n]-1)/K) + alpha[3]*((bij_moy(c,n,Y,K,C)-1)/K) + alpha[4]))
  test3=(1/(1+test1))^(I-1) * test1^(I-ii)
  result= test* test3
  return(result)
}


field_to_vect_voisin_inter = function(c,Y,K){
  D=dim(Y)
  if( is.null(D)) {N=1
  C= length(Y)
  
  
  temp = Y[c]-1
  res= rep(0,N)
  
  
  
  
  for(n in 1:N){
    res[n]=bij_inter(c,n,Y,K)
  }
  
  res = temp*bij_inter(c,K,rep(K,C),K) + res
  
  }else{
    C=D[1]
    N= D[2]
    
    temp = Y[c,]-1
    res= rep(0,N)
    
    for(n in 1:N){
      res[n]=bij_inter(c,n,Y,K)
    }
    
    res = temp*bij_inter(c,K,rep(K,C),K) + res }
  
  
  return(res)
}

# 1-5 
logitA_moy2= function(alpha,i,ii,I,K,Y,c,n){
  C=  length(Y[,n])
  
  test=comb(I-1,ii-1)
  
  test1= exp(-(alpha[1]*((i)/I)+alpha[2]*((Y[c,n])/K) + alpha[3]*((bij_moy(c,n,Y,K,C))/K) + alpha[4]))
  test3=(1/(1+test1))^(I-1) * test1^(I-ii)
  result= test* test3
  return(result)
}

#0-4
funcA_moy = function(alpha,I,K,C){
  tot=K
  A= array(rep(0,I*I*K*tot),dim=c(I,I,K*tot))
  
  
  for(i in 1:tot){
    for(k in 1:K){
      for(j in 1:I){
        for(l in 1:I){
          E=exp(-(alpha[1]*(j-1)/I+alpha[2]*(k-1)/K+ alpha[3]*(i-1)/tot+alpha[4]))
          A[j,l,(i+((k-1)*tot))]= comb(I-1,l-1)* (1/(1+E))^(I-1) * E^(I-l)
        }
      }
    }
    
  }
  return(A)
}

#1-5
funcA_moy2 = function(alpha,I,K,C){
  tot=K
  A= array(rep(0,I*I*K*tot),dim=c(I,I,K*tot))
  
  
  for(i in 1:tot){
    for(k in 1:K){
      for(j in 1:I){
        for(l in 1:I){
          E=exp(-(alpha[1]*(j)/I+alpha[2]*(k)/K+ alpha[3]*(i)/tot+alpha[4]))
          A[j,l,(i+((k-1)*tot))]= comb(I-1,l-1)* (1/(1+E))^(I-1) * E^(I-l)
        }
      }
    }
    
  }
  return(A)
}








#EM 
EtapeE_A2_BP= function(Y,etatA2,pi,A2,phi,C,I){
  D=dim(Y)
  if(is.null(D)){
    if(C!=1){N=1} else{N= length(Y)}
    
    
    
    
    
    
    
  } else{
    C=D[1]
    N= D[2]}
  
  
  
  # beug ici converge pas car des fois betaFB = nan ...
  #QQ= array(QQ,dim=c(I,K,bij_inter(1,1,rep(I,C),I))) dans le VEM inter ici on aura du (I,c,n)
  alphaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  betaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  gammaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  xiFB=array(rep(0,C*(N+1)*I*I),dim=c(I,I,C,N+1))
  for(c in 1:C){
    alphaFB[,c,1]=pi
    betaFB[,c,N+1]=rep(1,I)
    for(j in 1:N){
      k=N-j+1
      tempalpha= colSums(alphaFB[,c,j]*phi[,Y[c,j]]* A2[,,etatA2[c,j]])
      alphaFB[,c,j+1] = tempalpha/sum(tempalpha)
    
      tempbeta= rowSums(t(betaFB[,c,k+1]*t(phi[,Y[c,k]]* A2[,,etatA2[c,k]])))
      betaFB[,c,k]= tempbeta/sum(tempbeta)
      
    }
  }
  
  tempvrai=-1
  for(c in 1:C){
    for(j in 1:(N+1)){
      
      temp= alphaFB[,c,j]*betaFB[,c,j]
      
      Vrais=  sum(temp)
      
      
      
      gammaFB[,c,j]=  temp/Vrais
      #comprend pas 
      #test1=recpi*recphi[,2]* Achamp
      #colSums(test1[,,1])
      if(Vrais>tempvrai&& j!=(N+1)){tempvrai=Vrais}
      # xi[i_n-1,i_n,c,n]
      if( j !=1){
        
        xiFB[,,c,j]= t(betaFB[,c,j]*t(alphaFB[,c,j-1]*phi[,Y[c,j-1]]*A2[,,etatA2[c,j-1]]))/(sum(t(betaFB[,c,j]*t(alphaFB[,c,j-1]*phi[,Y[c,j-1]]*A2[,,etatA2[c,j-1]]))  ))  
        
        
      }
      
    }
  }
  
  
  
  result=list(gammaFB,xiFB,tempvrai)
  return(result)
  
  
}

# beaucoup plus lent que non bis ! 
EtapeE_A2_BP_bis= function(Y,etatA2,pi,A2,phi,C){
  D=dim(Y)
  if(is.null(D)){
    if(C!=1){N=1} else{N= length(Y)}

    
  } else{
    C=D[1]
    N= D[2]}
  
  
  
  # beug ici converge pas car des fois betaFB = nan ...
  #QQ= array(QQ,dim=c(I,K,bij_inter(1,1,rep(I,C),I))) dans le VEM inter ici on aura du (I,c,n)
  alphaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  betaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  gammaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  xiFB=array(rep(0,C*(N+1)*I*I),dim=c(I,I,C,N+1))
  for(c in 1:C){
    alphaFB[,c,1]=pi
    betaFB[,c,N+1]=rep(1,I)
    for(j in 1:N){
      k=N-j+1
      
      tempal=0
      tempbet=0
      for(px in 1:I){
        alphaFB[px,c,j+1]= sum(alphaFB[,c,j]*phi[,Y[c,j]]* A2[,px,etatA2[c,j]])
        tempal=tempal + sum(alphaFB[,c,j]*phi[,Y[c,j]]* A2[,px,etatA2[c,j]])
        betaFB[px,c,k]= sum(betaFB[,c,k+1]*phi[px,Y[c,k]]* A2[px,,etatA2[c,k]])
        tempbet=tempbet + sum(betaFB[,c,k+1]*phi[px,Y[c,k]]* A2[px,,etatA2[c,k]])
      }
 
      alphaFB[,c,j+1] = alphaFB[,c,j+1]/tempal
      
     
      betaFB[,c,k]= betaFB[,c,k]/tempbet
      
    }
  }
  
  tempvrai=-1
  for(c in 1:C){
    for(j in 1:(N+1)){
      
      temp= alphaFB[,c,j]*betaFB[,c,j]
      
      Vrais=  sum(temp)
      
      
      
      gammaFB[,c,j]=  temp/Vrais
      #comprend pas 
      #test1=recpi*recphi[,2]* Achamp
      #colSums(test1[,,1])
      if(Vrais>tempvrai&& j!=(N+1)){tempvrai=Vrais}
      # xi[i_n-1,i_n,c,n]
      if( j !=1){
        
        xiFB[,,c,j]= t(betaFB[,c,j]*t(alphaFB[,c,j-1]*phi[,Y[c,j-1]]*A2[,,etatA2[c,j-1]]))/(sum(t(betaFB[,c,j]*t(alphaFB[,c,j-1]*phi[,Y[c,j-1]]*A2[,,etatA2[c,j-1]]))  ))  
        
        
      }
      
    }
  }
  
  
  
  result=list(gammaFB,xiFB,tempvrai)
  return(result)
  
  
}
#etape M
gradphi_FB_gamma= function(gammaFB,Y,K,I){
  o=dim(Y)
  N=o[2]
  C=o[1]
  gradphi_FB_M= function(beta){
    
    
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          test1= exp(-(beta[1]*((i-1)/I)+beta[2]))
          
          
          sumn2= sumn2 + ((Y[c,n]-1)+(1-K)/(1+test1))*gammaFB[i,c,n]
          sumn1= sumn1 + ((Y[c,n]-1)*((i-1)/I)+(1-K)*((i-1)/I)/(1+test1))*gammaFB[i,c,n]
          
          
          
        }
      }
    }
    result1= -(sumn1)
    result2= -(sumn2)
    return(c(result1,result2))}
  return(gradphi_FB_M)
}
Eqphi_FB_gamma= function(gammaFB,Y,K,I){
  o=dim(Y)
  N=o[2]
  C=o[1]
  Eqphi_FB_M= function(beta){
    
    sumn=0
    for(c in 1:C){
      for(n in 1:N){
        for(i in 1:I){
          
          
          sumn= sumn + log(logitphi(beta,K,Y[c,n],i,I))*gammaFB[i,c,n]
        }
      }
    }
    result= -(sumn)
    return(result)}
  return(Eqphi_FB_M)
}
Eqphi_FB_gamma_bis= function(gammaFB,Y,K,I){
  o=dim(Y)
  N=o[2]
  C=o[1]
  
  
  Eqphi_FB_M= function(beta){
    phi= funcphi(beta,I,K)
    sumn=0
    for(c in 1:C){
      for(n in 1:N){
        
        sumn= sumn + sum(log(phi[,Y[c,n]])*gammaFB[,c,n])
        
      }
    }
    result= -(sumn)
    return(result)}
  return(Eqphi_FB_M)
}
#moyenne 0-4
EqA_moy= function(xiFB,Y,K,I){
  o=dim(Y)
  N=o[2]
  C=o[1]
  
  EqA_FBM= function(alpha){
    sumn=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          for(ii in 1 : I){
            sumn= sumn + log(logitA_moy(alpha,i,ii,I,K,Y,c,n))* xiFB[i,ii,c,n+1]
            } } } }
    result= -(sumn)
    return(result)}
  return(EqA_FBM)
}

EqA_moy_bis= function(xiFB,Y,K,I,etatA2){
  o=dim(Y)
  N=o[2]
  C=o[1]
  
  EqA_FBM_bis= function(alpha){
    sumn=0
    A=funcA_moy(alpha,I,K,C)
     for(c in 1:C){
      for(n in 1:N){
       sumn= sumn + sum(diag(log(A[,,etatA2[c,n]])%*% t(xiFB[,,c,n+1])))
      }
    }
    result= -(sumn)
    return(result)}
  return(EqA_FBM_bis)
}

gradA_moy= function(xiFB,Y,K,I,C){
  o=dim(Y)
  N=o[2]
  C=o[1]
  gradAM= function(alpha){
    
    sumo4=0
    sumn4=0
    sumo3=0
    sumn3=0
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          test1= exp(-(alpha[1]*((i-1)/I)+alpha[2]*((Y[c,n]-1)/K) + alpha[3]*((bij_moy(c,n,Y,K,C)-1)/K) +alpha[4]))
          for(ii in 1 : I){
            
            sumn4= sumn4 + ((ii-1)+(1-I)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn3= sumn3 + ((ii-1)*((bij_moy(c,n,Y,K,C)-1)/K)+(1-I)*((bij_moy(c,n,Y,K,C)-1)/K)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn2= sumn2 + ((ii-1)*((Y[c,n]-1)/K)+(1-I)*((Y[c,n]-1)/K)/(1+test1))* xiFB[i,ii,c,n+1]
            sumn1= sumn1 + ((ii-1)*((i-1)/I)+(1-I)*((i-1)/I)/(1+test1))*xiFB[i,ii,c,n+1]


          }
        }
      }
    }
    result1= -(sumn1)
    result2= -(sumn2)
    result3= -(sumn3)
    result4= -(sumn4)
    return(c(result1,result2,result3,result4))}
  return(gradAM)
}
EtapeM_moy = function(gammaFB,xiFB,Y,pi,alpha,beta,K){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  
  local_opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                      "xtol_rel"  = 1.0e-18 )
  opts <- list("algorithm"="NLOPT_LD_AUGLAG","xtol_rel"=1.0e-18,"maxeval"   = 1000,
               "local_opts" = local_opts )
  # epsa=rep(0,4)
  # epsb=rep(0,2)
  #   epsa[1:3]=log(alpha[1:3])
  #   epsa[4]=alpha[4]
  #   epsb[1]=log(beta[1])
  #   epsb[2]=beta[2]
  gradalpha=gradA_moy(xiFB,Y,K,I,C)
  falpha=EqA_moy(xiFB,Y,K,I)
  
  gradbeta=gradphi_FB_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_FB_gamma(gammaFB,Y,K,I)
  
  resuA=nloptr(x0 = alpha, eval_f = falpha,eval_grad_f=gradalpha,lb=c(0,0,0,-100),ub=c(100,100,100,0.01), opts = opts)
  
  resuphi=nloptr(x0 = beta, eval_f = fbeta,eval_grad_f=gradbeta,lb=c(0,-100), opts = opts)
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$solution,resuphi$solution)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}
# best is LBFGSB with 3 min 30 s avec bis ! 
EtapeM_moy_LBFGSB = function(gammaFB,xiFB,Y,pi,alpha,beta,K,etatA2){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  
   # epsa=rep(0,4)
  # epsb=rep(0,2)
  #   epsa[1:3]=log(alpha[1:3])
  #   epsa[4]=alpha[4]
  #   epsb[1]=log(beta[1])
  #   epsb[2]=beta[2]
  
  phi= funcphi(beta,I,K)
  A=funcA_moy(alpha,I,K,C)
  gradalpha=gradA_moy(xiFB,Y,K,I,C)
  falpha=EqA_moy_bis(xiFB,Y,K,I,etatA2)
  
  gradbeta=gradphi_FB_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_FB_gamma_bis(gammaFB,Y,K,I)
  
  resuA=optim(alpha, falpha,gradalpha,method=c("L-BFGS-B"),lower=c(0,0,0,-100),upper=c(75,75,75,0.01))
  
  resuphi=optim( beta, fbeta,gradbeta,method=c("L-BFGS-B"),lower=c(0,-100),upper=c(75,0.01))
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$par,resuphi$par)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}
EtapeM_moy_LBFGSB2= function(gammaFB,xiFB,Y,pi,alpha,beta,K,etatA2){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  
  # epsa=rep(0,4)
  # epsb=rep(0,2)
  #   epsa[1:3]=log(alpha[1:3])
  #   epsa[4]=alpha[4]
  #   epsb[1]=log(beta[1])
  #   epsb[2]=beta[2]
 
  
  gradalpha=gradA_moy(xiFB,Y,K,I,C)
  falpha=EqA_moy_bis(xiFB,Y,K,I,etatA2)
  
  gradbeta=gradphi_FB_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_FB_gamma_bis(gammaFB,Y,K,I)
  
  resuA=optimr(alpha, falpha,gradalpha,method=c("L-BFGS-B"),lower=c(0,0,0,-100),upper=c(75,75,75,0.01))
  
  resuphi=optimr( beta, fbeta,gradbeta,method=c("L-BFGS-B"),lower=c(0,-100),upper=c(75,0.01))
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$par,resuphi$par)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}

logitA_alp= function(alpha,i,ii,I,K,Y,c,n){
  C=  length(Y[,n])
  
  test=comb(I-1,ii-1)
  
  test1= exp(-(alpha[1]*((i-1)/I)+alpha[2]*((Y[c,n]-1)/K) + alpha[3]*((bij_inter(c,n,Y,K)-1)/bij_inter(c,n,rep(K,C),K)) + alpha[4]))
  test3=(1/(1+test1))^(I-1) * test1^(I-ii)
  result= test* test3
  return(result)
}
EqA_inter_bis= function(xiFB,Y,K,I,etatA2){
  EqA_FBM= function(alpha){
    sumn=0
     sumn=0
    A=funcA_2(alpha,I,K,C)
     for(c in 1:C){
      for(n in 1:N){
       sumn= sumn + sum(diag(log(A[,,etatA2[c,n]])%*% t(xiFB[,,c,n+1])))
      }
    }
    result= -(sumn)
    return(result)}
  return(EqA_FBM)
}
gradA_inter= function(xiFB,Y,K,I,C){
  
  gradAM= function(alpha){
    
    sumo4=0
    sumn4=0
    sumo3=0
    sumn3=0
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          test1= exp(-(alpha[1]*((i-1)/I)+alpha[2]*((Y[c,n]-1)/K) + alpha[3]*((bij_inter(c,n,Y,K)-1)/bij_inter(c,n,rep(K,C),K)) +alpha[4]))
          for(ii in 1 : I){
            
            sumn4= sumn4 + ((ii-1)+(1-I)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn3= sumn3 + ((ii-1)*((bij_inter(c,n,Y,K)-1)/bij_inter(c,n,rep(K,C),K))+(1-I)*((bij_inter(c,n,Y,K)-1)/bij_inter(c,n,rep(K,C),K))/(1+test1))*xiFB[i,ii,c,n+1]
            sumn2= sumn2 + ((ii-1)*((Y[c,n]-1)/K)+(1-I)*((Y[c,n]-1)/K)/(1+test1))* xiFB[i,ii,c,n+1]
            sumn1= sumn1 + ((ii-1)*((i-1)/I)+(1-I)*((i-1)/I)/(1+test1))*xiFB[i,ii,c,n+1]
            
            
          }
        }
      }
    }
    result1= -(sumn1)
    result2= -(sumn2)
    result3= -(sumn3)
    result4= -(sumn4)
    return(c(result1,result2,result3,result4))}
  return(gradAM)
}


EtapeM_inter_LBFGSB2= function(gammaFB,xiFB,Y,pi,alpha,beta,K,etatA2){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  
  # epsa=rep(0,4)
  # epsb=rep(0,2)
  #   epsa[1:3]=log(alpha[1:3])
  #   epsa[4]=alpha[4]
  #   epsb[1]=log(beta[1])
  #   epsb[2]=beta[2]
 
  
  gradalpha=gradA_inter(xiFB,Y,K,I,C)
  falpha=EqA_inter_bis(xiFB,Y,K,I,etatA2)
  
  gradbeta=gradphi_FB_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_FB_gamma_bis(gammaFB,Y,K,I)
  
  resuA=optimr(alpha, falpha,gradalpha,method=c("L-BFGS-B"),lower=c(0,0,0,-100),upper=c(70,70,70,0.01))
  
  resuphi=optimr( beta, fbeta,gradbeta,method=c("L-BFGS-B"),lower=c(0,-100),upper=c(70,0.01))
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$par,resuphi$par)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}

EM_alp = function(Y,etatA2,C,K,I,p){
  beta1=rep(0,2)
  alpha1=rep(0,4)
  beta1[1]=runif(1)*5
  beta1[2]=runif(1)*5-5
  gamma1=runif(1)*5
  
  alpha1[1:3]=runif(3)*5
  alpha1[4]=runif(1)*5-5
  
  recgamma=gamma1
  recbeta=beta1
  recalpha=alpha1
  recphi=funcphi(beta1,I,K)
  recpi=funcpi(gamma1,I)
  
  recA=funcA_2(alpha1,I,K,C)
  
  
  fullalpha=matrix(rep(0,4*(p+1)),nrow=4)
  fullbeta=matrix(rep(0,2*(p+1)),nrow=2)
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    
    
    fullalpha[,n]=recalpha
    fullbeta[,n]=recbeta
    fullgamma[,n]=recgamma
    
    
    RES= EtapeE_A2_BP(Y,etatA2,recpi,recA,recphi,C,I)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_inter_LBFGSB2(rho,xi,Y,recpi,recalpha,recbeta,K,etatA2)
    recalpha=RESu$alpha[1:4]
    tempA=RESu$alpha[5]
    recbeta=RESu$beta[1:2]
    tempphi=RESu$beta[3]
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    
    
    recpi=funcpi(recgamma,I)
    recA= funcA_2(recalpha,I,K,C)
    recphi=funcphi(recbeta,I,K)
    if(max(c(max(abs(recbeta -fullbeta[,n])) , max(abs(recalpha-fullalpha[,n])) ,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  fullalpha[,n]=recalpha
  fullbeta[,n]=recbeta
  fullgamma[,n]=recgamma
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  vrphi=funcphi(estimateurs$estimbeta,I,K)
  vrpi=funcpi(estimateurs$estimgamma,I)
  vrA=funcA_2(estimateurs$estimalpha,I,K,C)
  
  AB_FB=EtapeE_A2_BP_LH(Y,vrpi,vrA,vrphi,C)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH(alpha_FB,beta_FB,vrphi,Y,C,N,I)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  Xstate=viterbi_alp(Y,estimateurs$estimalpha,estimateurs$estimbeta,estimateurs$estimgamma,C,K,I,N)
  YN1state=predictFB(alpha_FB,vrphi)
  result=list(ini,estimateurs,full,endn,Vraisemblance,Xstate,YN1state,Y,etatA2)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai','predBG','predYN','Y','Yvoisin')
  
  return(result)
}


EtapeM_moy_nlm= function(gammaFB,xiFB,Y,pi,alpha,beta,K){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  
  # epsa=rep(0,4)
  # epsb=rep(0,2)
  #   epsa[1:3]=log(alpha[1:3])
  #   epsa[4]=alpha[4]
  #   epsb[1]=log(beta[1])
  #   epsb[2]=beta[2]
  gradalpha=gradA_moy(xiFB,Y,K,I,C)
  falpha=EqA_moy_bis(xiFB,Y,K,I)
  
  gradbeta=gradphi_FB_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_FB_gamma_bis(gammaFB,Y,K,I)
  
  resuA=nlm( falpha,alpha)
  
  resuphi=nlm( fbeta,beta)
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$estimate,resuphi$estimate)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}
EtapeM_moy_nlminb= function(gammaFB,xiFB,Y,pi,alpha,beta,K){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  
  # epsa=rep(0,4)
  # epsb=rep(0,2)
  #   epsa[1:3]=log(alpha[1:3])
  #   epsa[4]=alpha[4]
  #   epsb[1]=log(beta[1])
  #   epsb[2]=beta[2]
  gradalpha=gradA_moy(xiFB,Y,K,I,C)
  falpha=EqA_moy_bis(xiFB,Y,K,I)
  
  gradbeta=gradphi_FB_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_FB_gamma_bis(gammaFB,Y,K,I)
  
  resuA=nlminb( alpha,falpha,gradalpha,lower=c(0,0,0,-100),upper=c(100,100,100,0.01))
  
  resuphi=nlminb( beta,fbeta,gradbeta,lower=c(0,-100))
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$estimate,resuphi$estimate)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}




#sans S
EqA_moy_S= function(xiFB,Y,K,I,S,etatA2){
  o=dim(Y)
  N=o[2]
  C=o[1]
  
  EqA_FBM= function(alpha){
    sumn=0
    A=funcA_moy(c(S,alpha[1:3]),I,K,C)

        for(c in 1:C){
         
          for(n in 1:N){
            #ii = in i= in-1
            
            sumn= sumn + sum(diag(log(A[,,etatA2[c,n]])%*% t(xiFB[,,c,n+1])))

          }
        }
    result= -(sumn)
    return(result)}
  return(EqA_FBM)
}
gradA_moy_S= function(xiFB,Y,K,I,C,S){
  o=dim(Y)
  N=o[2]
  C=o[1]
  gradAM= function(alpha){
    
    sumo4=0
    sumn4=0
    sumo3=0
    sumn3=0
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          test1= exp(-(S*((i-1)/I)+alpha[1]*((Y[c,n]-1)/K) + alpha[2]*((bij_moy(c,n,Y,K,C)-1)/K) +alpha[3]))
          for(ii in 1 : I){
            
            sumn4= sumn4 + ((ii-1)+(1-I)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn3= sumn3 + ((ii-1)*((bij_moy(c,n,Y,K,C)-1)/K)+(1-I)*((bij_moy(c,n,Y,K,C)-1)/K)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn2= sumn2 + ((ii-1)*((Y[c,n]-1)/K)+(1-I)*((Y[c,n]-1)/K)/(1+test1))* xiFB[i,ii,c,n+1]
            sumn1= sumn1 + ((ii-1)*((i-1)/I)+(1-I)*((i-1)/I)/(1+test1))*xiFB[i,ii,c,n+1]
            
            
          }
        }
      }
    }
    
    result2= -(sumn2)
    result3= -(sumn3)
    result4= -(sumn4)
    return(c(result2,result3,result4))}
  return(gradAM)
}
EtapeM_moy_S = function(gammaFB,xiFB,Y,pi,alpha,beta,K,etatA2,S){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  alpha_S=c(alpha[2:4])
  

  
  
  gradalpha=gradA_moy_S(xiFB,Y,K,I,C,S)
  falpha=EqA_moy_S(xiFB,Y,K,I,S,etatA2)
  
  gradbeta=gradphi_FB_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_FB_gamma_bis(gammaFB,Y,K,I)
  
  resuA=optimr(alpha_S, falpha,gradalpha,method=c("L-BFGS-B"),lower=c(0,0,-100),upper=c(75,75,0.01))
  
  resuphi=optimr( beta, fbeta,gradbeta,method=c("L-BFGS-B"),lower=c(0,-100),upper=c(75,0.01))
  
  
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$par,resuphi$par)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}

#sans V
EqA_moy_V= function(xiFB,Y,K,I,V,etatA2){
  o=dim(Y)
  N=o[2]
  C=o[1]
  EqA_FBM= function(alpha){
    sumn=0
    A=funcA_moy(c(alpha[1:2],V,alpha[3]),I,K,C)

    for(c in 1:C){
     
      for(n in 1:N){
        #ii = in i= in-1
        
        sumn= sumn + sum(diag(log(A[,,etatA2[c,n]])%*% t(xiFB[,,c,n+1])))
        
        
        
        
      }
    }
  
    result= -(sumn)
    return(result)}
  return(EqA_FBM)
}
gradA_moy_V= function(xiFB,Y,K,I,C,V){
  o=dim(Y)
  N=o[2]
  C=o[1]
  gradAM= function(alpha){
    
    sumo4=0
    sumn4=0
    sumo3=0
    sumn3=0
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          test1= exp(-(alpha[1]*((i-1)/I)+alpha[2]*((Y[c,n]-1)/K) + V*((bij_moy(c,n,Y,K,C)-1)/K) +alpha[3]))
          for(ii in 1 : I){
            
            sumn4= sumn4 + ((ii-1)+(1-I)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn3= sumn3 + ((ii-1)*((bij_moy(c,n,Y,K,C)-1)/K)+(1-I)*((bij_moy(c,n,Y,K,C)-1)/K)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn2= sumn2 + ((ii-1)*((Y[c,n]-1)/K)+(1-I)*((Y[c,n]-1)/K)/(1+test1))* xiFB[i,ii,c,n+1]
            sumn1= sumn1 + ((ii-1)*((i-1)/I)+(1-I)*((i-1)/I)/(1+test1))*xiFB[i,ii,c,n+1]
            
            
          }
        }
      }
    }
    result1= -(sumn1)
    result2= -(sumn2)
    
    result4= -(sumn4)
    return(c(result1,result2,result4))}
  return(gradAM)
}
EtapeM_moy_V = function(gammaFB,xiFB,Y,pi,alpha,beta,K,etatA2,V){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  alpha_V=c(alpha[1:2],alpha[4])
 

  
  
  gradalpha=gradA_moy_V(xiFB,Y,K,I,C,V)
  falpha=EqA_moy_V(xiFB,Y,K,I,V,etatA2)
  
  gradbeta=gradphi_FB_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_FB_gamma_bis(gammaFB,Y,K,I)
  
  resuA=optimr(alpha_V, falpha,gradalpha,method=c("L-BFGS-B"),lower=c(0,0,-100),upper=c(75,75,0.01))
  
  resuphi=optimr( beta, fbeta,gradbeta,method=c("L-BFGS-B"),lower=c(0,-100),upper=c(75,0.01))
  
  
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$par,resuphi$par)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}

#sans D
EqA_moy_D= function(xiFB,Y,K,I,D,etatA2){
  o=dim(Y)
  N=o[2]
  C=o[1]
   EqA_FBM= function(alpha){
    sumn=0
    A=funcA_moy(c(alpha[1],D,alpha[2:3]),I,K,C)
    
    for(c in 1:C){
     
      for(n in 1:N){

        sumn= sumn + sum(diag(log(A[,,etatA2[c,n]])%*% t(xiFB[,,c,n+1])))
  
      }
    }
    result= -(sumn)
    return(result)}
  return(EqA_FBM)
}
gradA_moy_D= function(xiFB,Y,K,I,C,D){
  o=dim(Y)
  N=o[2]
  C=o[1]
  gradAM= function(alpha){
    
    sumo4=0
    sumn4=0
    sumo3=0
    sumn3=0
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          test1= exp(-(alpha[1]*((i-1)/I)+D*((Y[c,n]-1)/K) + alpha[2]*((bij_moy(c,n,Y,K,C)-1)/K) +alpha[3]))
          for(ii in 1 : I){
            
            sumn4= sumn4 + ((ii-1)+(1-I)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn3= sumn3 + ((ii-1)*((bij_moy(c,n,Y,K,C)-1)/K)+(1-I)*((bij_moy(c,n,Y,K,C)-1)/K)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn2= sumn2 + ((ii-1)*((Y[c,n]-1)/K)+(1-I)*((Y[c,n]-1)/K)/(1+test1))* xiFB[i,ii,c,n+1]
            sumn1= sumn1 + ((ii-1)*((i-1)/I)+(1-I)*((i-1)/I)/(1+test1))*xiFB[i,ii,c,n+1]
            
            
          }
        }
      }
    }
    result1= -(sumn1)
    result3= -(sumn3)
    
    result4= -(sumn4)
    return(c(result1,result3,result4))}
  return(gradAM)
}
EtapeM_moy_D = function(gammaFB,xiFB,Y,pi,alpha,beta,K,etatA2,D){
  DI=dim(Y)
  if(is.null(DI)){C=1
  N=length(Y)} else{
    C=DI[1]
    N= DI[2]}
  I= length(pi)
  alpha_D=c(alpha[1:2],alpha[4])

  gradalpha=gradA_moy_D(xiFB,Y,K,I,C,D)
  falpha=EqA_moy_D(xiFB,Y,K,I,D,etatA2)
  
  gradbeta=gradphi_FB_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_FB_gamma_bis(gammaFB,Y,K,I)
  
  resuA=optimr(alpha_D, falpha,gradalpha,method=c("L-BFGS-B"),lower=c(0,0,-100),upper=c(75,75,0.01))
  
  resuphi=optimr( beta, fbeta,gradbeta,method=c("L-BFGS-B"),lower=c(0,-100),upper=c(75,0.01))

  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$par,resuphi$par)
  names(result)= c('gamma','alpha','beta')
  
  return(result)
}

#sans cst
EqA_moy_cst= function(xiFB,Y,K,I,cst,etatA2){
  o=dim(Y)
  N=o[2]
  C=o[1]
   EqA_FBM= function(alpha){
    sumn=0
    A=funcA_moy(c(alpha[1:3],cst),I,K,C)
    
    for(c in 1:C){
      
      for(n in 1:N){
        #ii = in i= in-1
        
        sumn= sumn + sum(diag(log(A[,,etatA2[c,n]])%*% t(xiFB[,,c,n+1])))
        
      }
    }
    result= -(sumn)
    return(result)}
  return(EqA_FBM)
}
gradA_moy_cst= function(xiFB,Y,K,I,C,cst){
  o=dim(Y)
  N=o[2]
  C=o[1]
  gradAM= function(alpha){
    
    sumo4=0
    sumn4=0
    sumo3=0
    sumn3=0
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          test1= exp(-(alpha[1]*((i-1)/I)+alpha[2]*((Y[c,n]-1)/K) + alpha[3]*((bij_moy(c,n,Y,K,C)-1)/K) +cst))
          for(ii in 1 : I){
            
            sumn4= sumn4 + ((ii-1)+(1-I)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn3= sumn3 + ((ii-1)*((bij_moy(c,n,Y,K,C)-1)/K)+(1-I)*((bij_moy(c,n,Y,K,C)-1)/K)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn2= sumn2 + ((ii-1)*((Y[c,n]-1)/K)+(1-I)*((Y[c,n]-1)/K)/(1+test1))* xiFB[i,ii,c,n+1]
            sumn1= sumn1 + ((ii-1)*((i-1)/I)+(1-I)*((i-1)/I)/(1+test1))*xiFB[i,ii,c,n+1]
            
            
          }
        }
      }
    }
    
    result2= -(sumn2)
    result3= -(sumn3)
    result1= -(sumn1)
    return(c(result1,result2,result3))}
  return(gradAM)
}
EtapeM_moy_cst = function(gammaFB,xiFB,Y,pi,alpha,beta,K,etatA2,cst){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  alpha_cst=c(alpha[1:3]) 

  gradalpha=gradA_moy_cst(xiFB,Y,K,I,C,cst)
  falpha=EqA_moy_cst(xiFB,Y,K,I,cst,etatA2)
  
  gradbeta=gradphi_FB_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_FB_gamma_bis(gammaFB,Y,K,I)
  
  resuA=optimr(alpha_cst, falpha,gradalpha,method=c("L-BFGS-B"),lower=c(0,0,0),upper=c(75,75,75))
  
  resuphi=optimr( beta, fbeta,gradbeta,method=c("L-BFGS-B"),lower=c(0,-100),upper=c(75,75,75))
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$par,resuphi$par)
  names(result)= c('gamma','alpha','beta')
  
  return(result)
}

#sans V et D 
EqA_moy_VD= function(xiFB,Y,K,I,V,D,etatA2){
  o=dim(Y)
  N=o[2]
  C=o[1]
  EqA_FBM= function(alpha){
    sumn=0
    A=funcA_moy(c(alpha[1],D,V,alpha[2]),I,K,C)
    
    for(c in 1:C){
      
      for(n in 1:N){
        #ii = in i= in-1
        
        sumn= sumn + sum(diag(log(A[,,etatA2[c,n]])%*% t(xiFB[,,c,n+1])))
        
      }
    }
    result= -(sumn)
    return(result)}
  return(EqA_FBM)
}
gradA_moy_VD= function(xiFB,Y,K,I,C,V,D){
  
  gradAM= function(alpha){
    
    sumo4=0
    sumn4=0
    sumo3=0
    sumn3=0
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          test1= exp(-(alpha[1]*((i-1)/I)+D*((Y[c,n]-1)/K) + V*((bij_moy(c,n,Y,K,C)-1)/K) +alpha[2]))
          for(ii in 1 : I){
            
            sumn4= sumn4 + ((ii-1)+(1-I)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn3= sumn3 + ((ii-1)*((bij_moy(c,n,Y,K,C)-1)/K)+(1-I)*((bij_moy(c,n,Y,K,C)-1)/K)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn2= sumn2 + ((ii-1)*((Y[c,n]-1)/K)+(1-I)*((Y[c,n]-1)/K)/(1+test1))* xiFB[i,ii,c,n+1]
            sumn1= sumn1 + ((ii-1)*((i-1)/I)+(1-I)*((i-1)/I)/(1+test1))*xiFB[i,ii,c,n+1]
            
            
          }
        }
      }
    }
    result1= -(sumn1)
    
    
    result4= -(sumn4)
    return(c(result1,result4))}
  return(gradAM)
}
EtapeM_moy_VD = function(gammaFB,xiFB,Y,pi,alpha,beta,K,V,etatA2,D){
  DI=dim(Y)
  if(is.null(DI)){C=1
  N=length(Y)} else{
    C=DI[1]
    N= DI[2]}
  I= length(pi)
  alpha_VD=c(alpha[1],alpha[4])
  


  gradalpha=gradA_moy_VD(xiFB,Y,K,I,C,V,D)
  falpha=EqA_moy_VD(xiFB,Y,K,I,V,D,etatA2)
  
  gradbeta=gradphi_FB_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_FB_gamma_bis(gammaFB,Y,K,I)
  
  resuA=optimr(alpha_VD, falpha,gradalpha,method=c("L-BFGS-B"),lower=c(0,-100),upper=c(75,0.01))
  
  resuphi=optimr( beta, fbeta,gradbeta,method=c("L-BFGS-B"),lower=c(0,-100),,upper=c(75,0.01))
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$par,resuphi$par)
  names(result)= c('gamma','alpha','beta')
  
  return(result)
}


#culture en plus


#la variable cultd prend la l'anné et N et le champs c pour donnée le paramètre 
#associer a la culture ou le desherbage ou les deux. 
EtapeE_A2_BP_cult= function(Y,etatA2,pi,A2,phi,C,I,cultd){
  D=dim(Y)
  if(is.null(D)){
    if(C!=1){N=1} else{N= length(Y)}
    
    
    
    
    
    
    
  } else{
    C=D[1]
    N= D[2]}
  
  
  
  # beug ici converge pas car des fois betaFB = nan ...
  #QQ= array(QQ,dim=c(I,K,bij_inter(1,1,rep(I,C),I))) dans le VEM inter ici on aura du (I,c,n)
  alphaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  betaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  gammaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  xiFB=array(rep(0,C*(N+1)*I*I),dim=c(I,I,C,N+1))
  for(c in 1:C){
    alphaFB[,c,1]=pi
    betaFB[,c,N+1]=rep(1,I)
    for(j in 1:N){
      k=N-j+1
      tempalpha= colSums(alphaFB[,c,j]*phi[,Y[c,j],cultd[c,j]]* A2[,,etatA2[c,j],cultd[c,j]])
      alphaFB[,c,j+1] = tempalpha/sum(tempalpha)
      
      tempbeta= rowSums(t(betaFB[,c,k+1]*t(phi[,Y[c,k],cultd[c,k]]* A2[,,etatA2[c,k],cultd[c,k]])))
      betaFB[,c,k]= tempbeta/sum(tempbeta)
      
    }
  }
  
  tempvrai=-1
  for(c in 1:C){
    for(j in 1:(N+1)){
      
      temp= alphaFB[,c,j]*betaFB[,c,j]
      
      Vrais=  sum(temp)
      
      
      
      gammaFB[,c,j]=  temp/Vrais
      #comprend pas 
      #test1=recpi*recphi[,2]* Achamp
      #colSums(test1[,,1])
      if(Vrais>tempvrai&& j!=(N+1)){tempvrai=Vrais}
      # xi[i_n-1,i_n,c,n]
      if( j !=1){
        
        xiFB[,,c,j]= t(betaFB[,c,j]*t(alphaFB[,c,j-1]*phi[,Y[c,j-1],cultd[c,j-1]]*A2[,,etatA2[c,j-1],cultd[c,j-1]]))/(sum(t(betaFB[,c,j]*t(alphaFB[,c,j-1]*phi[,Y[c,j-1],cultd[c,j-1]]*A2[,,etatA2[c,j-1],cultd[c,j-1]]))  ))  
        
        
      }
      
    }
  }
  
  
  
  result=list(gammaFB,xiFB,tempvrai)
  return(result)
  
  
}

gradA_moy_cult= function(xiFB,Y,K,I,C,samecult){
  o=dim(Y)
  N=o[2]
  C=o[1]

  gradAM= function(alpha){
    nbrcult=dim(samecult)
    sumo4=0
    sumn4=0
    sumo3=0
    sumn3=0
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(cu in 1:nbrcult[1]){
      c=samecult[cu,1]
      n=samecult[cu,2]
      #ii = in i= in-1
      for(i in 1:I){
        test1= exp(-(alpha[1]*((i-1)/I)+alpha[2]*((Y[c,n]-1)/K) + alpha[3]*((bij_moy(c,n,Y,K,C)-1)/K) +alpha[4]))
        for(ii in 1 : I){
          
          sumn4= sumn4 + ((ii-1)+(1-I)/(1+test1))*xiFB[i,ii,c,n+1]
          sumn3= sumn3 + ((ii-1)*((bij_moy(c,n,Y,K,C)-1)/K)+(1-I)*((bij_moy(c,n,Y,K,C)-1)/K)/(1+test1))*xiFB[i,ii,c,n+1]
          sumn2= sumn2 + ((ii-1)*((Y[c,n]-1)/K)+(1-I)*((Y[c,n]-1)/K)/(1+test1))* xiFB[i,ii,c,n+1]
          sumn1= sumn1 + ((ii-1)*((i-1)/I)+(1-I)*((i-1)/I)/(1+test1))*xiFB[i,ii,c,n+1]
          
          
        }
      }
    }
    result1= -(sumn1)
    result2= -(sumn2)
    result3= -(sumn3)
    result4= -(sumn4)
    return(c(result1,result2,result3,result4))}
  return(gradAM)
}
EqA_moy_bis_cult= function(xiFB,Y,K,I,etatA2,samecult){
  o=dim(Y)
 
  
  EqA_FBM_bis= function(alpha){
    nbrcult=dim(samecult)
    sumn=0
    A=funcA_moy(alpha,I,K,C)
    for(cu in 1:nbrcult[1]){
      c=samecult[cu,1]
      n=samecult[cu,2]
      sumn= sumn + sum(diag(log(A[,,etatA2[c,n]])%*% t(xiFB[,,c,n+1])))
    }
    
    result= -(sumn)
    return(result)}
  return(EqA_FBM_bis)
}

gradphi_FB_gamma_cult= function(gammaFB,Y,K,I,samecult){
  o=dim(Y)
 
  gradphi_FB_M= function(beta){
    
    nbrcult=dim(samecult)
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    for(cu in 1:nbrcult[1]){
      c=samecult[cu,1]
      n=samecult[cu,2]
      #ii = in i= in-1
      for(i in 1:I){
        test1= exp(-(beta[1]*((i-1)/I)+beta[2]))
        
        
        sumn2= sumn2 + ((Y[c,n]-1)+(1-K)/(1+test1))*gammaFB[i,c,n]
        sumn1= sumn1 + ((Y[c,n]-1)*((i-1)/I)+(1-K)*((i-1)/I)/(1+test1))*gammaFB[i,c,n]
        
        
        
      }
    }
    result1= -(sumn1)
    result2= -(sumn2)
    return(c(result1,result2))}
  return(gradphi_FB_M)
}

Eqphi_FB_gamma_bis_cult= function(gammaFB,Y,K,I,samecult){
  o=dim(Y)
  
  
  
  Eqphi_FB_M= function(beta){
    nbrcult=dim(samecult)
    phi= funcphi(beta,I,K)
    sumn=0
    for(cu in 1:nbrcult[1]){
      c=samecult[cu,1]
      n=samecult[cu,2]
      
      sumn= sumn + sum(log(phi[,Y[c,n]])*gammaFB[,c,n])
      
    }
    
    result= -(sumn)
    return(result)}
  return(Eqphi_FB_M)
}

EtapeM_moy_LBFGSB_cult = function(gammaFB,xiFB,Y,pi,alpha,beta,K,etatA2,cultd){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  nbrcult=max(cultd)
  # epsa=rep(0,4)
  # epsb=rep(0,2)
  #   epsa[1:3]=log(alpha[1:3])
  #   epsa[4]=alpha[4]
  #   epsb[1]=log(beta[1])
  #   epsb[2]=beta[2]
  resalpha=array(rep(0,4*nbrcult),dim=c(4,nbrcult))
  resbeta=array(rep(0,2*nbrcult),dim=c(2,nbrcult))
  for(cu in 1 : max(cultd)){
    samecult=which(cu==cultd,arr.ind = TRUE)
    if(length(which(samecult==0))==0){}else{
      gradalpha=gradA_moy_cult(xiFB,Y,K,I,C,samecult)
      falpha=EqA_moy_bis_cult(xiFB,Y,K,I,etatA2,samecult)
      
      gradbeta=gradphi_FB_gamma_cult(gammaFB,Y,K,I,samecult)
      fbeta=Eqphi_FB_gamma_bis_cult(gammaFB,Y,K,I,samecult)
      
      resuA=optim(alpha[,cu], falpha,gradalpha,method=c("L-BFGS-B"),lower=c(0,0,0,-75),upper=c(70,70,70,0.01))
      resalpha[,cu]=resuA$par
      resuphi=optim( beta[,cu], fbeta,gradbeta,method=c("L-BFGS-B"),lower=c(0,-75),upper=c(70,0.01))
      resbeta[,cu]=resuphi$par}
  }
  
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resalpha,resbeta)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}

EM_moy_cult = function(Y,etatA2,C,K,I,p,cultd,newcult){
  o=dim(Y)
  N=o[2]
  C=o[1]
  nbrcult=max(cultd)
  tot=K
  gamma1=runif(1)*5
  recpi=funcpi(gamma1,I)
  
  recA= array(rep(0,I*I*K*tot*nbrcult),dim=c(I,I,K*tot,nbrcult))
  recphi=array(rep(0,I*K*nbrcult),dim=c(I,K,nbrcult))
  recgamma=gamma1
  recbeta=array(rep(0,2*nbrcult),dim=c(2,nbrcult))
  recalpha=array(rep(0,4*nbrcult),dim=c(4,nbrcult))
  
  for(cult in 1:nbrcult){
   
    beta1=rep(0,2)
    alpha1=rep(0,4)
    beta1[1]=runif(1)*5
    beta1[2]=runif(1)*10-5
    
    
    alpha1[1:3]=runif(3)*5
    alpha1[4]=runif(1)*5-5.01
    
    
    recbeta[,cult]=beta1
    recalpha[,cult]=alpha1
    recphi[,,cult]=funcphi(beta1,I,K)
    
    
    recA[,,,cult]=funcA_moy(alpha1,I,K,C)
    
  }


  fullalpha=array(rep(0,4*nbrcult*(p+1)),dim=c(4,nbrcult,p+1))
  fullbeta=array(rep(0,2*nbrcult*(p+1)),dim=c(2,nbrcult,p+1))
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    for(cult in 1:nbrcult){
    
    fullalpha[,cult,n]=recalpha[,cult]
    fullbeta[,cult,n]=recbeta[,cult]
    
    
    }
    
    fullgamma[,n]=recgamma
    if(max(cultd)==min(cultd)){
      RES= EtapeE_A2_BP(Y,etatA2,recpi,recA,recphi,C,I)
    }else{
    RES= EtapeE_A2_BP_cult(Y,etatA2,recpi,recA,recphi,C,I,cultd)}
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_moy_LBFGSB_cult(rho,xi,Y,recpi,recalpha,recbeta,K,etatA2,cultd)


    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    recpi=funcpi(recgamma,I)
    
    maxab=c(0,0)
    for(cult in 1:nbrcult){
      

      recbeta[,cult]=RESu$beta[,cult]
      recalpha[,cult]=RESu$alpha[,cult]
      recphi[,,cult]=funcphi(recbeta[,cult],I,K)
      
      
      recA[,,,cult]=funcA_moy(recalpha[,cult],I,K,C)
      
      
      if(maxab[1]<max(abs(recbeta[,cult] -fullbeta[,cult,n]))){maxab[1]=max(abs(recbeta[,cult] -fullbeta[,cult,n]))}
      if(maxab[2]< max(abs(recalpha[,cult] -fullalpha[,cult,n]))){maxab[2]= max(abs(recalpha[,cult] -fullalpha[,cult,n]))}
      
    }
   
   
    if(max(c(maxab,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  
  
  vrA= array(rep(0,I*I*K*tot*nbrcult),dim=c(I,I,K*tot,nbrcult))
  vrphi=array(rep(0,I*K*nbrcult),dim=c(I,K,nbrcult))
  for(cult in 1:nbrcult){
   
    
    recbeta[,cult]=beta1
    recalpha[,cult]=alpha1
    vrphi[,,cult]=funcphi(estimateurs$estimbeta,I,K)
    
    
    vrA[,,,cult]=funcA_moy(estimateurs$estimalpha,I,K,C)
    
  }
  vrpi=funcpi(estimateurs$estimgamma,I)
  
  AB_FB=EtapeE_moy_LH_cult(Y,vrpi,vrA,vrphi,C,cultd)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH_cult(alpha_FB,beta_FB,vrphi,Y,C,N,I,cultd)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  Xstate=viterbi_moy_cult(Y,estimateurs$estimalpha,estimateurs$estimbeta,estimateurs$estimgamma,C,K,I,N,cultd)
  YN1state=predictFB_cult(alpha_FB,vrphi,newcult)
  result=list(ini,estimateurs,full,endn,Vraisemblance,YN1state,Xstate,Y,etatA2)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai','predYN','predBG','Y','Yvoisin')
  
  return(result)
}

EtapeE_moy_LH_cult= function(Y,pi,A2,phi,C,etatA2,culd){
  
  D=dim(Y)
  if(is.null(D)){
    if(C!=1){N=1} else{N= length(Y)}
  } else{
    C=D[1]
    N= D[2]}
  
  
  
  # beug ici converge pas car des fois betaFB = nan ...
  #QQ= array(QQ,dim=c(I,K,bij_inter(1,1,rep(I,C),I))) dans le VEM inter ici on aura du (I,c,n)
  alphaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  betaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  gammaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  xiFB=array(rep(0,C*(N+1)*I*I),dim=c(I,I,C,N+1))
  for(c in 1:C){
    alphaFB[,c,1]=pi
    betaFB[,c,N+1]=rep(1,I)
    for(j in 1:N){
      k=N-j+1
      tempalpha= colSums(alphaFB[,c,j]*phi[,Y[c,j],culd[c,j]]* A2[,,etatA2[c,j],culd[c,j]])
      alphaFB[,c,j+1] = tempalpha/sum(tempalpha)
      
      tempbeta= rowSums(t(betaFB[,c,k+1]*t(phi[,Y[c,k],culd[c,k]]* A2[,,etatA2[c,k],culd[c,k]])))
      betaFB[,c,k]= tempbeta/sum(tempbeta)
      
    }
  }
  
  
  
  
  result=list(alphaFB,betaFB)
  return(result)
}

viterbi_moy_cult=function(Y,alpha,beta,gamma,C,K,I,N,cultd){
 
  
  nbrcult=max(cultd)
  tot=K
  A2= array(rep(0,I*I*K*tot*nbrcult),dim=c(I,I,K*tot,nbrcult))
  phi=array(rep(0,I*K*nbrcult),dim=c(I,K,nbrcult))
  for(cult in 1:nbrcult){
    
    

    phi[,,cult]=funcphi(beta,I,K)
    
    
    A2[,,,cult]=funcA_moy(alpha,I,K,C)
    
  }
  
 
  pi=funcpi(gamma,I)
  XX=array(rep(0,C*(N+1)),dim=c(C,N+1))
  delta=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  deltaind=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  
  
  
  for(c in 1 :C){
    delta[,c,1]=log(pi)
    etatA2=field_to_vect_moy(c,Y,K)
    for(n in 2:(N+1)){
      
      mad=delta[,c,n-1]+log(phi[,Y[c,n-1],cultd[c,n-1]])+log(A2[,,etatA2[n-1],cultd[c,n-1]])   
      delta[,c,n]=apply(mad,2,max)
      for(i in 1:I){
        deltaind[i,c,n-1]=which.max(mad[,i])
      }
      
    }
    for(i in 1:I){
      deltaind[1,c,N+1]=which.max(delta[,c,N+1])
    }
    
    XX[c,N+1]=deltaind[1,c,N+1]
    for(n in N:1){
      XX[c,n]=deltaind[XX[c,n+1],c,n]
      
    }
    
    
    
  }
  return(XX)
}

predictFB_cult= function(gammaFB,phi,newcult){
  D=dim(gammaFB)
  I=D[1]
  C=D[2]
  N=D[3]-1
  pred= rep(0,3)
  for(c in 1:C){
    sum=0
    for(i in 1:I){
      sum=sum+gammaFB[i,c,N+1]*phi[i,,newcult[c]]
    }
    pred[c]=which.max(sum)
    
    
  }
  return(pred)
}

BP_LLH_cult= function(alphaFB,betaFB,phi,Y,C,N,I,cultd){
  
  
  moy=0
  best = -Inf
  for(c in 1:C){
    for(j in 1: N){
      
      cvoulu=which((1:C)!=c)
      temp= alphaFB[,c,j]*betaFB[,c,j]
      firstpart=  log(sum(temp))
      
      prod2=1
      prod1=0
      for(n in 1:N){
        prodtemp = 1
        prodtemp1=0 
        for(u in cvoulu){
          sumit=0
          for(i in 1:I ){
            sumit=sumit + phi[i,Y[u,n],cultd[u,n]]*alphaFB[i,u,n]/sum(alphaFB[,u,n])
          }
          tempp=1
          prodtemp1= prodtemp1+log(sumit)
        }
        prod1=prod1+prodtemp1
      }
      if(best<firstpart+prod1){best=firstpart+prod1}
      moy= moy + firstpart+prod1
    }
  }
  
  result=list(best,moy/(C*N))
  return(result)
}

#full EM 0-4 + M step 10 times quicker ! 
EM_moy = function(Y,etatA2,C,K,I,p){
  o=dim(Y)
  N=o[2]
  C=o[1]
  beta1=rep(0,2)
  alpha1=rep(0,4)
  beta1[1]=runif(1)*5
  beta1[2]=runif(1)*10-5
  gamma1=runif(1)*5
  
  alpha1[1:3]=runif(3)*5
  alpha1[4]=runif(1)*5-5.01
  
  recgamma=gamma1
  recbeta=beta1
  recalpha=alpha1
  recphi=funcphi(beta1,I,K)
  recpi=funcpi(gamma1,I)
  
  recA=funcA_moy(alpha1,I,K,C)
  
  
  fullalpha=matrix(rep(0,4*(p+1)),nrow=4)
  fullbeta=matrix(rep(0,2*(p+1)),nrow=2)
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    
    
    fullalpha[,n]=recalpha
    fullbeta[,n]=recbeta
    fullgamma[,n]=recgamma
    
    
    RES= EtapeE_A2_BP(Y,etatA2,recpi,recA,recphi,C,I)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_moy_LBFGSB(rho,xi,Y,recpi,recalpha,recbeta,K,etatA2)
    recalpha=RESu$alpha[1:4]
    tempA=RESu$alpha[5]
    recbeta=RESu$beta[1:2]
    tempphi=RESu$beta[3]
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    
    
    recpi=funcpi(recgamma,I)
    recA= funcA_moy(recalpha,I,K,C)
    recphi=funcphi(recbeta,I,K)
    if(max(c(max(abs(recbeta -fullbeta[,n])) , max(abs(recalpha-fullalpha[,n])) ,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  fullalpha[,n]=recalpha
  fullbeta[,n]=recbeta
  fullgamma[,n]=recgamma
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  vrphi=funcphi(estimateurs$estimbeta,I,K)
  vrpi=funcpi(estimateurs$estimgamma,I)
  vrA=funcA_moy(estimateurs$estimalpha,I,K,C)
  
  AB_FB=EtapeE_moy_LH(Y,vrpi,vrA,vrphi,C)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH(alpha_FB,beta_FB,vrphi,Y,C,N,I)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  Xstate=viterbi_moy(Y,estimateurs$estimalpha,estimateurs$estimbeta,estimateurs$estimgamma,C,K,I,N)
  YN1state=predictFB(alpha_FB,vrphi)
  result=list(ini,estimateurs,full,endn,Vraisemblance,YN1state,Xstate,Y,etatA2)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai','predYN','predBG','Y','Yvoisin')
  
  return(result)
}
EM_moy_S = function(Y,etatA2,C,K,I,p,S){
  beta1=rep(0,2)
  alpha1=rep(0,4)
  beta1[1]=runif(1)*5
  beta1[2]=runif(1)*10-5
  gamma1=runif(1)*(-5)
  
  alpha1[2:3]=runif(2)*5
  alpha1[1]=S
  alpha1[4]=runif(1)*5-7
  
  recgamma=gamma1
  recbeta=beta1
  recalpha=alpha1
  recphi=funcphi(beta1,I,K)
  recpi=funcpi(gamma1,I)
  
  recA=funcA_moy(alpha1,I,K,C)
  
  
  fullalpha=matrix(rep(0,4*(p+1)),nrow=4)
  fullbeta=matrix(rep(0,2*(p+1)),nrow=2)
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    
    
    fullalpha[,n]=recalpha
    fullbeta[,n]=recbeta
    fullgamma[,n]=recgamma
    
    
    RES= EtapeE_A2_BP(Y,etatA2,recpi,recA,recphi,C,I)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_moy_S(rho,xi,Y,recpi,recalpha,recbeta,K,etatA2,S)
    recalpha=c(S,RESu$alpha[1:3])
    tempA=RESu$alpha[4]
    recbeta=RESu$beta[1:2]
    tempphi=RESu$beta[3]
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    
    
    recpi=funcpi(recgamma,I)
    recA= funcA_moy(recalpha,I,K,C)
    recphi=funcphi(recbeta,I,K)
    if(max(c(max(abs(recbeta -fullbeta[,n])) , max(abs(recalpha-fullalpha[,n])) ,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  fullalpha[,n]=recalpha
  fullbeta[,n]=recbeta
  fullgamma[,n]=recgamma
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  vrphi=funcphi(estimateurs$estimbeta,I,K)
  vrpi=funcpi(estimateurs$estimgamma,I)
  vrA=funcA_moy(estimateurs$estimalpha,I,K,C)
  
  AB_FB=EtapeE_moy_LH(Y,vrpi,vrA,vrphi,C)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH(alpha_FB,beta_FB,vrphi,Y,C,N,I)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  result=list(ini,estimateurs,full,endn,Vraisemblance)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai')
  
  return(result)
}
EM_moy_V = function(Y,etatA2,C,K,I,p,V){
  beta1=rep(0,2)
  alpha1=rep(0,4)
  beta1[1]=runif(1)*5
  beta1[2]=runif(1)*10-5
  gamma1=runif(1)*(-5)
  
  alpha1[1:2]=runif(2)*5
  alpha1[3]=V
  alpha1[4]=runif(1)*5-7
  
  recgamma=gamma1
  recbeta=beta1
  recalpha=alpha1
  recphi=funcphi(beta1,I,K)
  recpi=funcpi(gamma1,I)
  
  recA=funcA_moy(alpha1,I,K,C)
  
  
  fullalpha=matrix(rep(0,4*(p+1)),nrow=4)
  fullbeta=matrix(rep(0,2*(p+1)),nrow=2)
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    
    
    fullalpha[,n]=recalpha
    fullbeta[,n]=recbeta
    fullgamma[,n]=recgamma
    
    
    RES= EtapeE_A2_BP(Y,etatA2,recpi,recA,recphi,C,I)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_moy_V(rho,xi,Y,recpi,recalpha,recbeta,K,etatA2,V)
    recalpha=c(RESu$alpha[1:2],V,RESu$alpha[3])
    tempA=RESu$alpha[4]
    recbeta=RESu$beta[1:2]
    tempphi=RESu$beta[3]
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    
    
    recpi=funcpi(recgamma,I)
    recA= funcA_moy(recalpha,I,K,C)
    recphi=funcphi(recbeta,I,K)
    if(max(c(max(abs(recbeta -fullbeta[,n])) , max(abs(recalpha-fullalpha[,n])) ,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  fullalpha[,n]=recalpha
  fullbeta[,n]=recbeta
  fullgamma[,n]=recgamma
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  vrphi=funcphi(estimateurs$estimbeta,I,K)
  vrpi=funcpi(estimateurs$estimgamma,I)
  vrA=funcA_moy(estimateurs$estimalpha,I,K,C)
  
  AB_FB=EtapeE_moy_LH(Y,vrpi,vrA,vrphi,C)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH(alpha_FB,beta_FB,vrphi,Y,C,N,I)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  result=list(ini,estimateurs,full,endn,Vraisemblance)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai')
  
  return(result)
}
EM_moy_D = function(Y,etatA2,C,K,I,p,D){
  beta1=rep(0,2)
  alpha1=rep(0,4)
  beta1[1]=runif(1)*5
  beta1[2]=runif(1)*10-5
  gamma1=runif(1)*(-5)
  
  alpha1[1]=runif(1)*5
  alpha1[3]=runif(1)*5
  alpha1[4]=runif(1)*5-7
  alpha1[2]=D
  
  recgamma=gamma1
  recbeta=beta1
  recalpha=alpha1
  recphi=funcphi(beta1,I,K)
  recpi=funcpi(gamma1,I)
  
  recA=funcA_moy(alpha1,I,K,C)
  
  
  fullalpha=matrix(rep(0,4*(p+1)),nrow=4)
  fullbeta=matrix(rep(0,2*(p+1)),nrow=2)
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    
    
    fullalpha[,n]=recalpha
    fullbeta[,n]=recbeta
    fullgamma[,n]=recgamma
    
    
    RES= EtapeE_A2_BP(Y,etatA2,recpi,recA,recphi,C,I)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_moy_D(rho,xi,Y,recpi,recalpha,recbeta,K,etatA2,D)
    recalpha=c(RESu$alpha[1],D,RESu$alpha[2:3])
    tempA=RESu$alpha[4]
    recbeta=RESu$beta[1:2]
    tempphi=RESu$beta[3]
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    
    
    recpi=funcpi(recgamma,I)
    recA= funcA_moy(recalpha,I,K,C)
    recphi=funcphi(recbeta,I,K)
    if(max(c(max(abs(recbeta -fullbeta[,n])) , max(abs(recalpha-fullalpha[,n])) ,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  fullalpha[,n]=recalpha
  fullbeta[,n]=recbeta
  fullgamma[,n]=recgamma
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  vrphi=funcphi(estimateurs$estimbeta,I,K)
  vrpi=funcpi(estimateurs$estimgamma,I)
  vrA=funcA_moy(estimateurs$estimalpha,I,K,C)
  
  AB_FB=EtapeE_moy_LH(Y,vrpi,vrA,vrphi,C)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH(alpha_FB,beta_FB,vrphi,Y,C,N,I)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  result=list(ini,estimateurs,full,endn,Vraisemblance)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai')
  
  return(result)
}
EM_moy_cst = function(Y,etatA2,C,K,I,p,cst){
  beta1=rep(0,2)
  alpha1=rep(0,4)
  beta1[1]=runif(1)*5
  beta1[2]=runif(1)*10-5
  gamma1=runif(1)*(-5)
  
  alpha1[1:3]=runif(3)*5
  alpha1[4]=cst
  recgamma=gamma1
  recbeta=beta1
  recalpha=alpha1
  recphi=funcphi(beta1,I,K)
  recpi=funcpi(gamma1,I)
  
  recA=funcA_moy(alpha1,I,K,C)
  
  
  fullalpha=matrix(rep(0,4*(p+1)),nrow=4)
  fullbeta=matrix(rep(0,2*(p+1)),nrow=2)
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    
    
    fullalpha[,n]=recalpha
    fullbeta[,n]=recbeta
    fullgamma[,n]=recgamma
    
    
    RES= EtapeE_A2_BP(Y,etatA2,recpi,recA,recphi,C,I)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_moy_cst(rho,xi,Y,recpi,recalpha,recbeta,K,etatA2,cst)
    recalpha=c(RESu$alpha[1:3],cst)
    tempA=RESu$alpha[5]
    recbeta=RESu$beta[1:2]
    tempphi=RESu$beta[3]
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    
    
    recpi=funcpi(recgamma,I)
    recA= funcA_moy(recalpha,I,K,C)
    recphi=funcphi(recbeta,I,K)
    if(max(c(max(abs(recbeta -fullbeta[,n])) , max(abs(recalpha-fullalpha[,n])) ,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  fullalpha[,n]=recalpha
  fullbeta[,n]=recbeta
  fullgamma[,n]=recgamma
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  vrphi=funcphi(estimateurs$estimbeta,I,K)
  vrpi=funcpi(estimateurs$estimgamma,I)
  vrA=funcA_moy(estimateurs$estimalpha,I,K,C)
  
  AB_FB=EtapeE_moy_LH(Y,vrpi,vrA,vrphi,C)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH(alpha_FB,beta_FB,vrphi,Y,C,N,I)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  result=list(ini,estimateurs,full,endn,Vraisemblance)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai')
  
  return(result)
}
EM_moy_VD = function(Y,etatA2,C,K,I,p,V,D){
  beta1=rep(0,2)
  alpha1=rep(0,4)
  beta1[1]=runif(1)*5
  beta1[2]=runif(1)*10-5
  gamma1=runif(1)*(-5)
  
  alpha1[1]=runif(1)*5
  alpha1[3]=V
  alpha1[2]=D
  alpha1[4]=runif(1)*5-7
  
  recgamma=gamma1
  recbeta=beta1
  recalpha=alpha1
  recphi=funcphi(beta1,I,K)
  recpi=funcpi(gamma1,I)
  
  recA=funcA_moy(alpha1,I,K,C)
  
  
  fullalpha=matrix(rep(0,4*(p+1)),nrow=4)
  fullbeta=matrix(rep(0,2*(p+1)),nrow=2)
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    
    
    fullalpha[,n]=recalpha
    fullbeta[,n]=recbeta
    fullgamma[,n]=recgamma
    
    
    RES= EtapeE_A2_BP(Y,etatA2,recpi,recA,recphi,C,I)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_moy_VD(rho,xi,Y,recpi,recalpha,recbeta,K,V,etatA2,D)
    recalpha=c(RESu$alpha[1],D,V,RESu$alpha[2])
    tempA=RESu$alpha[5]
    recbeta=RESu$beta[1:2]
    tempphi=RESu$beta[3]
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    
    
    recpi=funcpi(recgamma,I)
    recA= funcA_moy(recalpha,I,K,C)
    recphi=funcphi(recbeta,I,K)
    if(max(c(max(abs(recbeta -fullbeta[,n])) , max(abs(recalpha-fullalpha[,n])) ,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  fullalpha[,n]=recalpha
  fullbeta[,n]=recbeta
  fullgamma[,n]=recgamma
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  vrphi=funcphi(estimateurs$estimbeta,I,K)
  vrpi=funcpi(estimateurs$estimgamma,I)
  vrA=funcA_moy(estimateurs$estimalpha,I,K,C)
  
  AB_FB=EtapeE_moy_LH(Y,vrpi,vrA,vrphi,C)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH(alpha_FB,beta_FB,vrphi,Y,C,N,I)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  result=list(ini,estimateurs,full,endn,Vraisemblance)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai')
  
  return(result)
}

EM_moy_ZI = function(Y,etatA2,C,K,I,p){
  o=dim(Y)
  N=o[2]
  C=o[1]
  beta1=rep(0,2)
  alpha1=rep(0,4)
  beta1[1]=runif(1)*5
  beta1[2]=runif(1)*10-5
  gamma1=runif(1)*5
  
  alpha1[1:3]=runif(3)*5
  alpha1[4]=runif(1)*5-5.01
  
  recgamma=gamma1
  recbeta=beta1
  recalpha=alpha1
  recphi=ZIphi(beta1,I,K)
  recpi=funcpi(gamma1,I)
  
  recA=funcA_moy(alpha1,I,K,C)
  
  
  fullalpha=matrix(rep(0,4*(p+1)),nrow=4)
  fullbeta=matrix(rep(0,2*(p+1)),nrow=2)
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    
    
    fullalpha[,n]=recalpha
    fullbeta[,n]=recbeta
    fullgamma[,n]=recgamma
    
    
    RES= EtapeE_A2_BP(Y,etatA2,recpi,recA,recphi,C,I)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_ZI_LBFGSB(rho,xi,Y,recpi,recalpha,recbeta,K,etatA2)
    recalpha=RESu$alpha[1:4]
    tempA=RESu$alpha[5]
    recbeta=RESu$beta[1:2]
    tempphi=RESu$beta[3]
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    
    
    recpi=funcpi(recgamma,I)
    recA= funcA_moy(recalpha,I,K,C)
    recphi=ZIphi(recbeta,I,K)
    if(max(c(max(abs(recbeta -fullbeta[,n])) , max(abs(recalpha-fullalpha[,n])) ,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  fullalpha[,n]=recalpha
  fullbeta[,n]=recbeta
  fullgamma[,n]=recgamma
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  vrphi=ZIphi(estimateurs$estimbeta,I,K)
  vrpi=funcpi(estimateurs$estimgamma,I)
  vrA=funcA_moy(estimateurs$estimalpha,I,K,C)
  
  AB_FB=EtapeE_moy_LH(Y,vrpi,vrA,vrphi,C)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH(alpha_FB,beta_FB,vrphi,Y,C,N,I)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  Xstate=viterbi_ZImoy(Y,estimateurs$estimalpha,estimateurs$estimbeta,estimateurs$estimgamma,C,K,I,N)
  YN1state=predictFB(alpha_FB,vrphi)
  result=list(ini,estimateurs,full,endn,Vraisemblance,YN1state,Xstate,Y,etatA2)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai','predYN','predBG','Y','Yvoisin')
  
  return(result)
}
EtapeM_ZI_LBFGSB = function(gammaFB,xiFB,Y,pi,alpha,beta,K,etatA2){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  
  # epsa=rep(0,4)
  # epsb=rep(0,2)
  #   epsa[1:3]=log(alpha[1:3])
  #   epsa[4]=alpha[4]
  #   epsb[1]=log(beta[1])
  #   epsb[2]=beta[2]
  
  phi= ZIphi(beta,I,K)
  A=funcA_moy(alpha,I,K,C)
  gradalpha=gradA_moy(xiFB,Y,K,I,C)
  falpha=EqA_moy_bis(xiFB,Y,K,I,etatA2)
  
  gradbeta=gradphi_ZI_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_ZI_gamma_bis(gammaFB,Y,K,I)
  
  resuA=optim(alpha, falpha,gradalpha,method=c("L-BFGS-B"),lower=c(0,0,0,-70),upper=c(70,70,70,0.01))
  
  resuphi=optim( beta, fbeta,gradbeta,method=c("L-BFGS-B"),lower=c(0,-70),upper=c(70,0.01))
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$par,resuphi$par)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}
gradphi_ZI_gamma= function(gammaFB,Y,K,I){
  o=dim(Y)
  N=o[2]
  C=o[1]
  gradphi_FB_M= function(beta){
    
    
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 2:I){
          test1= exp(-(beta[1]*((i-1)/I)+beta[2]))
          
          
          sumn2= sumn2 + ((Y[c,n]-1)+(1-K)/(1+test1))*gammaFB[i,c,n]
          sumn1= sumn1 + ((Y[c,n]-1)*((i-1)/I)+(1-K)*((i-1)/I)/(1+test1))*gammaFB[i,c,n]
          
          
          
        }
      }
    }
    result1= -(sumn1)
    result2= -(sumn2)
    return(c(result1,result2))}
  return(gradphi_FB_M)
}
Eqphi_ZI_gamma_bis= function(gammaFB,Y,K,I){
  o=dim(Y)
  N=o[2]
  C=o[1]
  temp=gammaFB[2:I,,]
  
  
  Eqphi_FB_M= function(beta){
    phi= ZIphi(beta,I,K)
    phi=phi[2:I,]
    sumn=0
    for(c in 1:C){
      for(n in 1:N){
        
        sumn= sumn + sum(log(phi[,Y[c,n]])*temp[,c,n])
        
      }
    }
    result= -(sumn)
    return(result)}
  return(Eqphi_FB_M)
}
viterbi_ZImoy=function(Y,alpha,beta,gamma,C,K,I,N){
  A2=funcA_moy(alpha,I,K,C)
  phi=ZIphi(beta,I,K)
  pi=funcpi(gamma,I)
  XX=array(rep(0,C*(N+1)),dim=c(C,N+1))
  delta=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  deltaind=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  
  
  
  for(c in 1 :C){
    delta[,c,1]=log(pi)
    etatA2=field_to_vect_moy(c,Y,K)
    for(n in 2:(N+1)){
      
      mad=delta[,c,n-1]+log(phi[,Y[c,n-1]])+log(A2[,,etatA2[n-1]])   
      delta[,c,n]=apply(mad,2,max)
      for(i in 1:I){
        deltaind[i,c,n-1]=which.max(mad[,i])
      }
      
    }
    for(i in 1:I){
      deltaind[1,c,N+1]=which.max(delta[,c,N+1])
    }
    
    XX[c,N+1]=deltaind[1,c,N+1]
    for(n in N:1){
      XX[c,n]=deltaind[XX[c,n+1],c,n]
      
    }
    
    
    
  }
  return(XX)
}



gradA_moy_cult_dist_nu0= function(xiFB,Y,K,I,C,samecult,dist,nu0){
  o=dim(Y)
  N=o[2]
  C=o[1]
  
  gradAM= function(alpha){
    nbrcult=dim(samecult)
    sumo4=0
    sumn4=0
    sumo3=0
    sumn3=0
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(cu in 1:nbrcult[1]){
      c=samecult[cu,1]
      n=samecult[cu,2]
      #ii = in i= in-1
      for(i in 1:I){
        test1= exp(-(alpha[1]*((i-1)/I)+alpha[2]*((Y[c,n]-1)/K) + alpha[3]*((bij_dist(c,n,Y,K,C,dist)-1)/K) +nu0))
        for(ii in 1 : I){
          
      
          sumn3= sumn3 + ((ii-1)*((bij_dist(c,n,Y,K,C,dist)-1)/K)+(1-I)*((bij_dist(c,n,Y,K,C,dist)-1)/K)/(1+test1))*xiFB[i,ii,c,n+1]
          sumn2= sumn2 + ((ii-1)*((Y[c,n]-1)/K)+(1-I)*((Y[c,n]-1)/K)/(1+test1))* xiFB[i,ii,c,n+1]
          sumn1= sumn1 + ((ii-1)*((i-1)/I)+(1-I)*((i-1)/I)/(1+test1))*xiFB[i,ii,c,n+1]
          
          
        }
      }
    }
    result1= -(sumn1)
    result2= -(sumn2)
    result3= -(sumn3)

    return(c(result1,result2,result3))}
  return(gradAM)
}
EqA_moy_bis_cult_nu0= function(xiFB,Y,K,I,etatA2,samecult,nu0){
  o=dim(Y)
  
  
  EqA_FBM_bis= function(alpha){
    nbrcult=dim(samecult)
    sumn=0
    A=funcA_moy(c(alpha,nu0),I,K,C)
    for(cu in 1:nbrcult[1]){
      c=samecult[cu,1]
      n=samecult[cu,2]
      sumn= sumn + sum(diag(log(A[,,etatA2[c,n]])%*% t(xiFB[,,c,n+1])))
    }
    
    result= -(sumn)
    return(result)}
  return(EqA_FBM_bis)
}
EtapeM_moy_LBFGSB_cult_dist_nu0 = function(gammaFB,xiFB,Y,pi,alpha,beta,K,etatA2,cultd,dist){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  nbrcult=max(cultd)
  # epsa=rep(0,4)
  # epsb=rep(0,2)
  #   epsa[1:3]=log(alpha[1:3])
  #   epsa[4]=alpha[4]
  #   epsb[1]=log(beta[1])
  #   epsb[2]=beta[2]
  resalpha=array(rep(0,4*nbrcult),dim=c(4,nbrcult))
  resbeta=array(rep(0,2*nbrcult),dim=c(2,nbrcult))
  for(cu in 1 : max(cultd)){
    samecult=which(cu==cultd,arr.ind = TRUE)
    if(length(samecult)==0){}else{
      gradalpha=gradA_moy_cult_dist_nu0(xiFB,Y,K,I,C,samecult,dist,alpha[4,cu])
      falpha=EqA_moy_bis_cult_nu0(xiFB,Y,K,I,etatA2,samecult,alpha[4,cu])
      
      gradbeta=gradphi_FB_gamma_cult_dist(gammaFB,Y,K,I,samecult)
      fbeta=Eqphi_FB_gamma_bis_cult_dist(gammaFB,Y,K,I,samecult)
      
      resuA=optim(alpha[-4,cu], falpha,gradalpha,method=c("L-BFGS-B"),lower=c(0,0,0,-75),upper=c(70,70,70,0.01))
      resalpha[,cu]=c(resuA$par,alpha[4,cu])
      resuphi=optim( beta[,cu], fbeta,gradbeta,method=c("L-BFGS-B"),lower=c(0,-75),upper=c(70,0.01))
      resbeta[,cu]=resuphi$par}
  }
  
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resalpha,resbeta)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}
EM_cult_dist_zi_nu0 = function(Y,etatA2,C,K,I,p,cultd,newcult,dist,nu0){
  o=dim(Y)
  N=o[2]
  C=o[1]
  nbrcult=max(cultd)
  tot=K
  gamma1=runif(1)*5
  recpi=funcpi(gamma1,I)
  
  recA= array(rep(0,I*I*K*tot*nbrcult),dim=c(I,I,K*tot,nbrcult))
  recphi=array(rep(0,I*K*nbrcult),dim=c(I,K,nbrcult))
  recgamma=gamma1
  recbeta=array(rep(0,2*nbrcult),dim=c(2,nbrcult))
  recalpha=array(rep(0,4*nbrcult),dim=c(4,nbrcult))
  
  for(cult in 1:nbrcult){
    
    beta1=rep(0,2)
    alpha1=rep(0,4)
    beta1[1]=runif(1)*5
    beta1[2]=runif(1)*10-5
    
    
    alpha1[1:3]=runif(3)*5
    alpha1[4]=nu0
    
    
    recbeta[,cult]=beta1
    recalpha[,cult]=alpha1
    recphi[,,cult]=ZIphi(beta1,I,K)
    
    recA[,,,cult]=funcA_moy(alpha1,I,K,C)
    
  }
  
  
  fullalpha=array(rep(0,4*nbrcult*(p+1)),dim=c(4,nbrcult,p+1))
  fullbeta=array(rep(0,2*nbrcult*(p+1)),dim=c(2,nbrcult,p+1))
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    for(cult in 1:nbrcult){
      
      fullalpha[,cult,n]=recalpha[,cult]
      fullbeta[,cult,n]=recbeta[,cult]
      
      
    }
    
    fullgamma[,n]=recgamma
    
    RES= EtapeE_A2_BP_cult(Y,etatA2,recpi,recA,recphi,C,I,cultd)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_moy_LBFGSB_cult_dist_nu0(rho,xi,Y,recpi,recalpha,recbeta,K,etatA2,cultd,dist)
    
    
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    recpi=funcpi(recgamma,I)
    
    maxab=c(0,0)
    for(cult in 1:nbrcult){
      
      
      recbeta[,cult]=RESu$beta[,cult]
      recalpha[,cult]=RESu$alpha[,cult]
      recphi[,,cult]=ZIphi(recbeta[,cult],I,K)
      
      
      recA[,,,cult]=funcA_moy(recalpha[,cult],I,K,C)
      
      
      if(maxab[1]<max(abs(recbeta[,cult] -fullbeta[,cult,n]))){maxab[1]=max(abs(recbeta[,cult] -fullbeta[,cult,n]))}
      if(maxab[2]< max(abs(recalpha[,cult] -fullalpha[,cult,n]))){maxab[2]= max(abs(recalpha[,cult] -fullalpha[,cult,n]))}
      
    }
    
    
    if(max(c(maxab,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  
  
  vrA= array(rep(0,I*I*K*tot*nbrcult),dim=c(I,I,K*tot,nbrcult))
  vrphi=array(rep(0,I*K*nbrcult),dim=c(I,K,nbrcult))
  for(cult in 1:nbrcult){
    
    
    
    vrphi[,,cult]=ZIphi(estimateurs$estimbeta[,cult],I,K)
    
    
    vrA[,,,cult]=funcA_moy(estimateurs$estimalpha[,cult],I,K,C)
    
  }
  vrpi=funcpi(estimateurs$estimgamma,I)
  
  AB_FB=EtapeE_moy_LH_cult(Y,vrpi,vrA,vrphi,C,etatA2,cultd)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH_cult(alpha_FB,beta_FB,vrphi,Y,C,N,I,cultd)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  Xstate=viterbi_moy_cult_dist(Y,estimateurs$estimalpha,estimateurs$estimbeta,estimateurs$estimgamma,C,K,I,N,cultd,dist)
  YN1state=predictFB_cult(alpha_FB,vrphi,newcult)
  result=list(ini,estimateurs,full,endn,Vraisemblance,YN1state,Xstate,Y,etatA2)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai','predYN','predBG','Y','Yvoisin')
  
  return(result)
}



gradA_moy_cult_dist= function(xiFB,Y,K,I,C,samecult,dist){
  o=dim(Y)
  N=o[2]
  C=o[1]
  
  gradAM= function(alpha){
    nbrcult=dim(samecult)
    sumo4=0
    sumn4=0
    sumo3=0
    sumn3=0
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(cu in 1:nbrcult[1]){
      c=samecult[cu,1]
      n=samecult[cu,2]
      #ii = in i= in-1
      for(i in 1:I){
        test1= exp(-(alpha[1]*((i-1)/I)+alpha[2]*((Y[c,n]-1)/K) + alpha[3]*((bij_dist(c,n,Y,K,C,dist)-1)/K) +alpha[4]))
        for(ii in 1 : I){
          
          sumn4= sumn4 + ((ii-1)+(1-I)/(1+test1))*xiFB[i,ii,c,n+1]
          sumn3= sumn3 + ((ii-1)*((bij_dist(c,n,Y,K,C,dist)-1)/K)+(1-I)*((bij_dist(c,n,Y,K,C,dist)-1)/K)/(1+test1))*xiFB[i,ii,c,n+1]
          sumn2= sumn2 + ((ii-1)*((Y[c,n]-1)/K)+(1-I)*((Y[c,n]-1)/K)/(1+test1))* xiFB[i,ii,c,n+1]
          sumn1= sumn1 + ((ii-1)*((i-1)/I)+(1-I)*((i-1)/I)/(1+test1))*xiFB[i,ii,c,n+1]
          
          
        }
      }
    }
    result1= -(sumn1)
    result2= -(sumn2)
    result3= -(sumn3)
    result4= -(sumn4)
    return(c(result1,result2,result3,result4))}
  return(gradAM)
}
gradphi_FB_gamma_cult_dist= function(gammaFB,Y,K,I,samecult){
  o=dim(Y)
  
  gradphi_FB_M= function(beta){
    
    nbrcult=dim(samecult)
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    for(cu in 1:nbrcult[1]){
      c=samecult[cu,1]
      n=samecult[cu,2]
      #ii = in i= in-1
      for(i in 2:I){
        test1= exp(-(beta[1]*((i-1)/I)+beta[2]))
        
        
        sumn2= sumn2 + ((Y[c,n]-1)+(1-K)/(1+test1))*gammaFB[i,c,n]
        sumn1= sumn1 + ((Y[c,n]-1)*((i-1)/I)+(1-K)*((i-1)/I)/(1+test1))*gammaFB[i,c,n]
        
        
        
      }
    }
    result1= -(sumn1)
    result2= -(sumn2)
    return(c(result1,result2))}
  return(gradphi_FB_M)
}
Eqphi_FB_gamma_bis_cult_dist= function(gammaFB,Y,K,I,samecult){
  o=dim(Y)
  
  
  
  Eqphi_FB_M= function(beta){
    nbrcult=dim(samecult)
    phi= ZIphi(beta,I,K)
    sumn=0
    for(cu in 1:nbrcult[1]){
      c=samecult[cu,1]
      n=samecult[cu,2]
      beug=log(phi[,Y[c,n]])
      if(length(which(beug==-Inf))!=0){beug[which(beug==-Inf)]=0}
      sumn= sumn + sum(beug*gammaFB[,c,n])
    
    }
    
    result= -(sumn)
    return(result)}
  return(Eqphi_FB_M)
}
EtapeM_moy_LBFGSB_cult_dist = function(gammaFB,xiFB,Y,pi,alpha,beta,K,etatA2,cultd,dist){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  nbrcult=max(cultd)
  # epsa=rep(0,4)
  # epsb=rep(0,2)
  #   epsa[1:3]=log(alpha[1:3])
  #   epsa[4]=alpha[4]
  #   epsb[1]=log(beta[1])
  #   epsb[2]=beta[2]
  resalpha=array(rep(0,4*nbrcult),dim=c(4,nbrcult))
  resbeta=array(rep(0,2*nbrcult),dim=c(2,nbrcult))
  for(cu in 1 : max(cultd)){
    samecult=which(cu==cultd,arr.ind = TRUE)
    if(length(samecult)==0){}else{
      gradalpha=gradA_moy_cult_dist(xiFB,Y,K,I,C,samecult,dist)
      falpha=EqA_moy_bis_cult(xiFB,Y,K,I,etatA2,samecult)
      
      gradbeta=gradphi_FB_gamma_cult_dist(gammaFB,Y,K,I,samecult)
      fbeta=Eqphi_FB_gamma_bis_cult_dist(gammaFB,Y,K,I,samecult)
      
      resuA=optim(alpha[,cu], falpha,gradalpha,method=c("L-BFGS-B"),lower=c(0,0,0,-75),upper=c(70,70,70,0.01))
      resalpha[,cu]=resuA$par
      resuphi=optim( beta[,cu], fbeta,gradbeta,method=c("L-BFGS-B"),lower=c(0,-75),upper=c(70,0.01))
      resbeta[,cu]=resuphi$par}
  }
  
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resalpha,resbeta)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}

EtapeM_moy_LBFGSB_cult_dist_neg = function(gammaFB,xiFB,Y,pi,alpha,beta,K,etatA2,cultd,dist){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  nbrcult=max(cultd)
  # epsa=rep(0,4)
  # epsb=rep(0,2)
  #   epsa[1:3]=log(alpha[1:3])
  #   epsa[4]=alpha[4]
  #   epsb[1]=log(beta[1])
  #   epsb[2]=beta[2]
  resalpha=array(rep(0,4*nbrcult),dim=c(4,nbrcult))
  resbeta=array(rep(0,2*nbrcult),dim=c(2,nbrcult))
  for(cu in 1 : max(cultd)){
    samecult=which(cu==cultd,arr.ind = TRUE)
    if(length(samecult)==0){}else{
      gradalpha=gradA_moy_cult_dist(xiFB,Y,K,I,C,samecult,dist)
      falpha=EqA_moy_bis_cult(xiFB,Y,K,I,etatA2,samecult)
      
      gradbeta=gradphi_FB_gamma_cult_dist(gammaFB,Y,K,I,samecult)
      fbeta=Eqphi_FB_gamma_bis_cult_dist(gammaFB,Y,K,I,samecult)
      
      resuA=optim(alpha[,cu], falpha,gradalpha,method=c("L-BFGS-B"),lower=c(0,-20,0,-75),upper=c(70,70,70,0.01))
      resalpha[,cu]=resuA$par
      resuphi=optim( beta[,cu], fbeta,gradbeta,method=c("L-BFGS-B"),lower=c(0,-75),upper=c(70,0.01))
      resbeta[,cu]=resuphi$par}
  }
  
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resalpha,resbeta)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}
EM_cult_dist_zi_neg = function(Y,etatA2,C,K,I,p,cultd,newcult,dist){
  o=dim(Y)
  N=o[2]
  C=o[1]
  nbrcult=max(cultd)
  tot=K
  gamma1=runif(1)*5
  recpi=funcpi(gamma1,I)
  
  recA= array(rep(0,I*I*K*tot*nbrcult),dim=c(I,I,K*tot,nbrcult))
  recphi=array(rep(0,I*K*nbrcult),dim=c(I,K,nbrcult))
  recgamma=gamma1
  recbeta=array(rep(0,2*nbrcult),dim=c(2,nbrcult))
  recalpha=array(rep(0,4*nbrcult),dim=c(4,nbrcult))
  
  for(cult in 1:nbrcult){
    
    beta1=rep(0,2)
    alpha1=rep(0,4)
    beta1[1]=runif(1)*5
    beta1[2]=runif(1)*10-5
    
    
    alpha1[1:3]=runif(3)*5
    alpha1[4]=runif(1)*5-5.01
    
    
    recbeta[,cult]=beta1
    recalpha[,cult]=alpha1
    recphi[,,cult]=ZIphi(beta1,I,K)
    
    recA[,,,cult]=funcA_moy(alpha1,I,K,C)
    
  }
  
  
  fullalpha=array(rep(0,4*nbrcult*(p+1)),dim=c(4,nbrcult,p+1))
  fullbeta=array(rep(0,2*nbrcult*(p+1)),dim=c(2,nbrcult,p+1))
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    for(cult in 1:nbrcult){
      
      fullalpha[,cult,n]=recalpha[,cult]
      fullbeta[,cult,n]=recbeta[,cult]
      
      
    }
    
    fullgamma[,n]=recgamma
    
    RES= EtapeE_A2_BP_cult(Y,etatA2,recpi,recA,recphi,C,I,cultd)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_moy_LBFGSB_cult_dist_neg(rho,xi,Y,recpi,recalpha,recbeta,K,etatA2,cultd,dist)
    
    
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    recpi=funcpi(recgamma,I)
    
    maxab=c(0,0)
    for(cult in 1:nbrcult){
      
      
      recbeta[,cult]=RESu$beta[,cult]
      recalpha[,cult]=RESu$alpha[,cult]
      recphi[,,cult]=ZIphi(recbeta[,cult],I,K)
      
      
      recA[,,,cult]=funcA_moy(recalpha[,cult],I,K,C)
      
      
      if(maxab[1]<max(abs(recbeta[,cult] -fullbeta[,cult,n]))){maxab[1]=max(abs(recbeta[,cult] -fullbeta[,cult,n]))}
      if(maxab[2]< max(abs(recalpha[,cult] -fullalpha[,cult,n]))){maxab[2]= max(abs(recalpha[,cult] -fullalpha[,cult,n]))}
      
    }
    
    
    if(max(c(maxab,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  
  
  vrA= array(rep(0,I*I*K*tot*nbrcult),dim=c(I,I,K*tot,nbrcult))
  vrphi=array(rep(0,I*K*nbrcult),dim=c(I,K,nbrcult))
  for(cult in 1:nbrcult){
    
    
    
    vrphi[,,cult]=ZIphi(estimateurs$estimbeta[,cult],I,K)
    
    
    vrA[,,,cult]=funcA_moy(estimateurs$estimalpha[,cult],I,K,C)
    
  }
  vrpi=funcpi(estimateurs$estimgamma,I)
  
  AB_FB=EtapeE_moy_LH_cult(Y,vrpi,vrA,vrphi,C,etatA2,cultd)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH_cult(alpha_FB,beta_FB,vrphi,Y,C,N,I,cultd)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  Xstate=viterbi_moy_cult_dist(Y,estimateurs$estimalpha,estimateurs$estimbeta,estimateurs$estimgamma,C,K,I,N,cultd,dist)
  YN1state=predictFB_cult(alpha_FB,vrphi,newcult)
  result=list(ini,estimateurs,full,endn,Vraisemblance,YN1state,Xstate,Y,etatA2)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai','predYN','predBG','Y','Yvoisin')
  
  return(result)
}


EM_cult_dist_zi = function(Y,etatA2,C,K,I,p,cultd,newcult,dist){
  o=dim(Y)
  N=o[2]
  C=o[1]
  nbrcult=max(cultd)
  tot=K
  gamma1=runif(1)*5
  recpi=funcpi(gamma1,I)
  
  recA= array(rep(0,I*I*K*tot*nbrcult),dim=c(I,I,K*tot,nbrcult))
  recphi=array(rep(0,I*K*nbrcult),dim=c(I,K,nbrcult))
  recgamma=gamma1
  recbeta=array(rep(0,2*nbrcult),dim=c(2,nbrcult))
  recalpha=array(rep(0,4*nbrcult),dim=c(4,nbrcult))
  
  for(cult in 1:nbrcult){
    
    beta1=rep(0,2)
    alpha1=rep(0,4)
    beta1[1]=runif(1)*5
    beta1[2]=runif(1)*10-5
    
    
    alpha1[1:3]=runif(3)*5
    alpha1[4]=runif(1)*5-5.01
    
    
    recbeta[,cult]=beta1
    recalpha[,cult]=alpha1
    recphi[,,cult]=ZIphi(beta1,I,K)
    
    recA[,,,cult]=funcA_moy(alpha1,I,K,C)
    
  }
  
  
  fullalpha=array(rep(0,4*nbrcult*(p+1)),dim=c(4,nbrcult,p+1))
  fullbeta=array(rep(0,2*nbrcult*(p+1)),dim=c(2,nbrcult,p+1))
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    for(cult in 1:nbrcult){
      
      fullalpha[,cult,n]=recalpha[,cult]
      fullbeta[,cult,n]=recbeta[,cult]
      
      
    }
    
    fullgamma[,n]=recgamma
    
    RES= EtapeE_A2_BP_cult(Y,etatA2,recpi,recA,recphi,C,I,cultd)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_moy_LBFGSB_cult_dist(rho,xi,Y,recpi,recalpha,recbeta,K,etatA2,cultd,dist)
    
    
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    recpi=funcpi(recgamma,I)
    
    maxab=c(0,0)
    for(cult in 1:nbrcult){
      
      
      recbeta[,cult]=RESu$beta[,cult]
      recalpha[,cult]=RESu$alpha[,cult]
      recphi[,,cult]=ZIphi(recbeta[,cult],I,K)
      
      
      recA[,,,cult]=funcA_moy(recalpha[,cult],I,K,C)
      
      
      if(maxab[1]<max(abs(recbeta[,cult] -fullbeta[,cult,n]))){maxab[1]=max(abs(recbeta[,cult] -fullbeta[,cult,n]))}
      if(maxab[2]< max(abs(recalpha[,cult] -fullalpha[,cult,n]))){maxab[2]= max(abs(recalpha[,cult] -fullalpha[,cult,n]))}
      
    }
    
    
    if(max(c(maxab,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  
  
  vrA= array(rep(0,I*I*K*tot*nbrcult),dim=c(I,I,K*tot,nbrcult))
  vrphi=array(rep(0,I*K*nbrcult),dim=c(I,K,nbrcult))
  for(cult in 1:nbrcult){
    
    
 
    vrphi[,,cult]=ZIphi(estimateurs$estimbeta[,cult],I,K)
    
    
    vrA[,,,cult]=funcA_moy(estimateurs$estimalpha[,cult],I,K,C)
    
  }
  vrpi=funcpi(estimateurs$estimgamma,I)
  
  AB_FB=EtapeE_moy_LH_cult(Y,vrpi,vrA,vrphi,C,etatA2,cultd)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH_cult(alpha_FB,beta_FB,vrphi,Y,C,N,I,cultd)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  Xstate=viterbi_moy_cult_dist(Y,estimateurs$estimalpha,estimateurs$estimbeta,estimateurs$estimgamma,C,K,I,N,cultd,dist)
  YN1state=predictFB_cult(alpha_FB,vrphi,newcult)
  result=list(ini,estimateurs,full,endn,Vraisemblance,YN1state,Xstate,Y,etatA2)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai','predYN','predBG','Y','Yvoisin')
  
  return(result)
}
EM_cult_dist_zi_ini = function(Y,etatA2,C,K,I,p,cultd,newcult,dist,alpha1,beta1,gamma1){
  o=dim(Y)
  N=o[2]
  C=o[1]
  nbrcult=max(cultd)
  tot=K
  gamma1=runif(1)*5
  recpi=funcpi(gamma1,I)
  
  recA= array(rep(0,I*I*K*tot*nbrcult),dim=c(I,I,K*tot,nbrcult))
  recphi=array(rep(0,I*K*nbrcult),dim=c(I,K,nbrcult))
  recgamma=gamma1
  recbeta=array(rep(0,2*nbrcult),dim=c(2,nbrcult))
  recalpha=array(rep(0,4*nbrcult),dim=c(4,nbrcult))
  
  for(cult in 1:nbrcult){
    

    
    
    recbeta[,cult]=beta1[,cult]
    recalpha[,cult]=alpha1[,cult]
    recphi[,,cult]=ZIphi(beta1,I,K)
    
    recA[,,,cult]=funcA_moy(alpha1,I,K,C)
    
  }
  
  
  fullalpha=array(rep(0,4*nbrcult*(p+1)),dim=c(4,nbrcult,p+1))
  fullbeta=array(rep(0,2*nbrcult*(p+1)),dim=c(2,nbrcult,p+1))
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    for(cult in 1:nbrcult){
      
      fullalpha[,cult,n]=recalpha[,cult]
      fullbeta[,cult,n]=recbeta[,cult]
      
      
    }
    
    fullgamma[,n]=recgamma
    
    RES= EtapeE_A2_BP_cult(Y,etatA2,recpi,recA,recphi,C,I,cultd)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_moy_LBFGSB_cult_dist(rho,xi,Y,recpi,recalpha,recbeta,K,etatA2,cultd,dist)
    
    
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    recpi=funcpi(recgamma,I)
    
    maxab=c(0,0)
    for(cult in 1:nbrcult){
      
      
      recbeta[,cult]=RESu$beta[,cult]
      recalpha[,cult]=RESu$alpha[,cult]
      recphi[,,cult]=ZIphi(recbeta[,cult],I,K)
      
      
      recA[,,,cult]=funcA_moy(recalpha[,cult],I,K,C)
      
      
      
      
      if(maxab[1]<max(abs(recbeta[,cult] -fullbeta[,cult,n]))){maxab[1]=max(abs(recbeta[,cult] -fullbeta[,cult,n]))}
      if(maxab[2]< max(abs(recalpha[,cult] -fullalpha[,cult,n]))){maxab[2]= max(abs(recalpha[,cult] -fullalpha[,cult,n]))}
      
    }
    
#     
#     AB_FB=EtapeE_moy_LH_cult(Y,recpi,recA,recphi,C,etatA2,cultd)
#     alpha_FB=AB_FB[[1]]
#     beta_FB=AB_FB[[2]]
#     LLH=BP_LLH_cult(alpha_FB,beta_FB,recphi,Y,C,N,I,cultd)
#     print('LLH =')
#     print(LLH)
    
    
    if(max(c(maxab,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  
  
  vrA= array(rep(0,I*I*K*tot*nbrcult),dim=c(I,I,K*tot,nbrcult))
  vrphi=array(rep(0,I*K*nbrcult),dim=c(I,K,nbrcult))
  for(cult in 1:nbrcult){
    
    
    
    vrphi[,,cult]=ZIphi(estimateurs$estimbeta[,cult],I,K)
    
    
    vrA[,,,cult]=funcA_moy(estimateurs$estimalpha[,cult],I,K,C)
    
  }
  vrpi=funcpi(estimateurs$estimgamma,I)
  
  AB_FB=EtapeE_moy_LH_cult(Y,vrpi,vrA,vrphi,C,etatA2,cultd)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH_cult(alpha_FB,beta_FB,vrphi,Y,C,N,I,cultd)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  Xstate=viterbi_moy_cult_dist(Y,estimateurs$estimalpha,estimateurs$estimbeta,estimateurs$estimgamma,C,K,I,N,cultd,dist)
  YN1state=predictFB_cult(alpha_FB,vrphi,newcult)
  result=list(ini,estimateurs,full,endn,Vraisemblance,YN1state,Xstate,Y,etatA2)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai','predYN','predBG','Y','Yvoisin')
  
  return(result)
}

viterbi_moy_cult_dist=function(Y,alpha,beta,gamma,C,K,I,N,cultd,dist){
  
  
  nbrcult=max(cultd)
  tot=K
  A2= array(rep(0,I*I*K*tot*nbrcult),dim=c(I,I,K*tot,nbrcult))
  phi=array(rep(0,I*K*nbrcult),dim=c(I,K,nbrcult))
  for(cult in 1:nbrcult){
    
    
    
    phi[,,cult]=ZIphi(beta,I,K)
    
    
    A2[,,,cult]=funcA_moy(alpha,I,K,C)
    
  }
  
  
  pi=funcpi(gamma,I)
  XX=array(rep(0,C*(N+1)),dim=c(C,N+1))
  delta=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  deltaind=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  
  
  
  for(c in 1 :C){
    delta[,c,1]=log(pi)
    etatA2=field_to_vect_dist(c,Y,K,dist)
    for(n in 2:(N+1)){
      
      mad=delta[,c,n-1]+log(phi[,Y[c,n-1],cultd[c,n-1]])+log(A2[,,etatA2[n-1],cultd[c,n-1]])   
      delta[,c,n]=apply(mad,2,max)
      for(i in 1:I){
        deltaind[i,c,n-1]=which.max(mad[,i])
      }
      
    }
    for(i in 1:I){
      deltaind[1,c,N+1]=which.max(delta[,c,N+1])
    }
    
    XX[c,N+1]=deltaind[1,c,N+1]
    for(n in N:1){
      XX[c,n]=deltaind[XX[c,n+1],c,n]
      
    }
    
    
    
  }
  return(XX)
}


EqA_moy2= function(xiFB,Y,K,I){
  EqA_FBM= function(alpha){
    sumn=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          for(ii in 1 : I){
            
            sumn= sumn + log(logitA_moy2(alpha,i,ii,I,K,Y,c,n))* xiFB[i,ii,c,n+1]
            
            
          }
        }
      }
    }
    result= -(sumn)
    return(result)}
  return(EqA_FBM)
}
gradA_moy2= function(xiFB,Y,K,I,C){
  
  gradAM= function(alpha){
    
    sumo4=0
    sumn4=0
    sumo3=0
    sumn3=0
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          test1= exp(-(alpha[1]*((i)/I)+alpha[2]*((Y[c,n])/K) + alpha[3]*((bij_moy(c,n,Y,K,C))/K) +alpha[4]))
          for(ii in 1 : I){
            
            sumn4= sumn4 + ((ii-1)+(1-I)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn3= sumn3 + ((ii-1)*((bij_moy(c,n,Y,K,C))/K)+(1-I)*((bij_moy(c,n,Y,K,C))/K)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn2= sumn2 + ((ii-1)*((Y[c,n])/K)+(1-I)*((Y[c,n])/K)/(1+test1))* xiFB[i,ii,c,n+1]
            sumn1= sumn1 + ((ii-1)*((i)/I)+(1-I)*((i)/I)/(1+test1))*xiFB[i,ii,c,n+1]
            
            
          }
        }
      }
    }
    result1= -(sumn1)
    result2= -(sumn2)
    result3= -(sumn3)
    result4= -(sumn4)
    return(c(result1,result2,result3,result4))}
  return(gradAM)
}
EtapeM_moy2 = function(gammaFB,xiFB,Y,pi,alpha,beta,K){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  
  local_opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                      "xtol_rel"  = 1.0e-18 )
  opts <- list("algorithm"="NLOPT_LD_AUGLAG","xtol_rel"=1.0e-18,"maxeval"   = 1000,
               "local_opts" = local_opts )
  # epsa=rep(0,4)
  # epsb=rep(0,2)
  #   epsa[1:3]=log(alpha[1:3])
  #   epsa[4]=alpha[4]
  #   epsb[1]=log(beta[1])
  #   epsb[2]=beta[2]
  gradalpha=gradA_moy2(xiFB,Y,K,I,C)
  falpha=EqA_moy2(xiFB,Y,K,I)
  
  gradbeta=gradphi_FB_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_FB_gamma(gammaFB,Y,K,I)
  
  resuA=nloptr(x0 = alpha, eval_f = falpha,eval_grad_f=gradalpha,lb=c(0,0,0,-100),ub=c(100,100,100,0.01), opts = opts)
  
  resuphi=nloptr(x0 = beta, eval_f = fbeta,eval_grad_f=gradbeta,lb=c(0,-100), opts = opts)
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$solution,resuphi$solution)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}

EqA_moy_V2= function(xiFB,Y,K,I,V){
  EqA_FBM= function(alpha){
    sumn=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          for(ii in 1 : I){
            
            sumn= sumn + log(logitA_moy2(c(alpha[1:2],V,alpha[3]),i,ii,I,K,Y,c,n))* xiFB[i,ii,c,n+1]
            
            
          }
        }
      }
    }
    result= -(sumn)
    return(result)}
  return(EqA_FBM)
}
gradA_moy_V2= function(xiFB,Y,K,I,C,V){
  
  gradAM= function(alpha){
    
    sumo4=0
    sumn4=0
    sumo3=0
    sumn3=0
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          test1= exp(-(alpha[1]*((i)/I)+alpha[2]*((Y[c,n])/K) + V*((bij_moy(c,n,Y,K,C))/K) +alpha[3]))
          for(ii in 1 : I){
            
            sumn4= sumn4 + ((ii-1)+(1-I)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn3= sumn3 + ((ii-1)*((bij_moy(c,n,Y,K,C))/K)+(1-I)*((bij_moy(c,n,Y,K,C))/K)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn2= sumn2 + ((ii-1)*((Y[c,n])/K)+(1-I)*((Y[c,n])/K)/(1+test1))* xiFB[i,ii,c,n+1]
            sumn1= sumn1 + ((ii-1)*((i)/I)+(1-I)*((i/I)/(1+test1)))*xiFB[i,ii,c,n+1]
            
            
          }
        }
      }
    }
    result1= -(sumn1)
    result2= -(sumn2)
    
    result4= -(sumn4)
    return(c(result1,result2,result4))}
  return(gradAM)
}
EtapeM_moy_V2 = function(gammaFB,xiFB,Y,pi,alpha,beta,K,V){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  alpha_V=c(alpha[1:2],alpha[4])
  local_opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                      "xtol_rel"  = 1.0e-18 )
  opts <- list("algorithm"="NLOPT_LD_AUGLAG","xtol_rel"=1.0e-18,"maxeval"   = 1000,
               "local_opts" = local_opts )
  # epsa=rep(0,4)
  # epsb=rep(0,2)
  #   epsa[1:3]=log(alpha[1:3])
  #   epsa[4]=alpha[4]
  #   epsb[1]=log(beta[1])
  #   epsb[2]=beta[2]
  gradalpha=gradA_moy_V2(xiFB,Y,K,I,C,V)
  falpha=EqA_moy_V2(xiFB,Y,K,I,V)
  
  gradbeta=gradphi_FB_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_FB_gamma(gammaFB,Y,K,I)
  
  resuA=nloptr(x0 = alpha_V, eval_f = falpha,eval_grad_f=gradalpha,lb=c(0,0,-100),ub=c(100,100,0.01), opts = opts)
  
  resuphi=nloptr(x0 = beta, eval_f = fbeta,eval_grad_f=gradbeta,lb=c(0,-100), opts = opts)
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$solution,resuphi$solution)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}

EqA_moy_S2= function(xiFB,Y,K,I,S){
  EqA_FBM= function(alpha){
    sumn=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          for(ii in 1 : I){
            
            sumn= sumn + log(logitA_moy2(c(S,alpha[1:3]),i,ii,I,K,Y,c,n))* xiFB[i,ii,c,n+1]
            
            
          }
        }
      }
    }
    result= -(sumn)
    return(result)}
  return(EqA_FBM)
}
gradA_moy_S2= function(xiFB,Y,K,I,C,S){
  
  gradAM= function(alpha){
    
    sumo4=0
    sumn4=0
    sumo3=0
    sumn3=0
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          test1= exp(-(S*((i)/I)+alpha[1]*((Y[c,n])/K) + alpha[2]*((bij_moy(c,n,Y,K,C))/K) +alpha[3]))
          for(ii in 1 : I){
            
            sumn4= sumn4 + ((ii-1)+(1-I)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn3= sumn3 + ((ii-1)*((bij_moy(c,n,Y,K,C))/K)+(1-I)*((bij_moy(c,n,Y,K,C))/K)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn2= sumn2 + ((ii-1)*((Y[c,n])/K)+(1-I)*((Y[c,n])/K)/(1+test1))* xiFB[i,ii,c,n+1]
            sumn1= sumn1 + ((ii-1)*((i)/I)+(1-I)*((i/I)/(1+test1)))*xiFB[i,ii,c,n+1]
            
            
          }
        }
      }
    }
    
    result2= -(sumn2)
    result3= -(sumn3)
    result4= -(sumn4)
    return(c(result2,result3,result4))}
  return(gradAM)
}
EtapeM_moy_S2 = function(gammaFB,xiFB,Y,pi,alpha,beta,K,S){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  alpha_S=c(alpha[2:4])
  local_opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                      "xtol_rel"  = 1.0e-18 )
  opts <- list("algorithm"="NLOPT_LD_AUGLAG","xtol_rel"=1.0e-18,"maxeval"   = 1000,
               "local_opts" = local_opts )
  # epsa=rep(0,4)
  # epsb=rep(0,2)
  #   epsa[1:3]=log(alpha[1:3])
  #   epsa[4]=alpha[4]
  #   epsb[1]=log(beta[1])
  #   epsb[2]=beta[2]
  gradalpha=gradA_moy_S2(xiFB,Y,K,I,C,S)
  falpha=EqA_moy_S2(xiFB,Y,K,I,S)
  
  gradbeta=gradphi_FB_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_FB_gamma(gammaFB,Y,K,I)
  
  resuA=nloptr(x0 = alpha_S, eval_f = falpha,eval_grad_f=gradalpha,lb=c(0,0,-100),ub=c(100,100,0.01), opts = opts)
  
  resuphi=nloptr(x0 = beta, eval_f = fbeta,eval_grad_f=gradbeta,lb=c(0,-100), opts = opts)
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$solution,resuphi$solution)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}


EqA_moy_cst2= function(xiFB,Y,K,I,cst){
  EqA_FBM= function(alpha){
    sumn=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          for(ii in 1 : I){
            
            sumn= sumn + log(logitA_moy2(c(alpha[1:3],cst),i,ii,I,K,Y,c,n))* xiFB[i,ii,c,n+1]
            
            
          }
        }
      }
    }
    result= -(sumn)
    return(result)}
  return(EqA_FBM)
}
gradA_moy_cst2= function(xiFB,Y,K,I,C,cst){
  
  gradAM= function(alpha){
    
    sumo4=0
    sumn4=0
    sumo3=0
    sumn3=0
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          test1= exp(-(alpha[1]*((i)/I)+alpha[2]*((Y[c,n])/K) + alpha[3]*((bij_moy(c,n,Y,K,C))/K) +cst))
          for(ii in 1 : I){
            
            sumn4= sumn4 + ((ii-1)+(1-I)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn3= sumn3 + ((ii-1)*((bij_moy(c,n,Y,K,C))/K)+(1-I)*((bij_moy(c,n,Y,K,C))/K)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn2= sumn2 + ((ii-1)*((Y[c,n])/K)+(1-I)*((Y[c,n])/K)/(1+test1))* xiFB[i,ii,c,n+1]
            sumn1= sumn1 + ((ii-1)*((i)/I)+(1-I)*((i/I)/(1+test1)))*xiFB[i,ii,c,n+1]
            
            
          }
        }
      }
    }
    
    result2= -(sumn2)
    result3= -(sumn3)
    result1= -(sumn1)
    return(c(result1,result2,result3))}
  return(gradAM)
}
EtapeM_moy_cst2 = function(gammaFB,xiFB,Y,pi,alpha,beta,K,cst){
  D=dim(Y)
  if(is.null(D)){C=1
  N=length(Y)} else{
    C=D[1]
    N= D[2]}
  I= length(pi)
  alpha_cst=c(alpha[1:3],cst)
  local_opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                      "xtol_rel"  = 1.0e-18 )
  opts <- list("algorithm"="NLOPT_LD_AUGLAG","xtol_rel"=1.0e-18,"maxeval"   = 1000,
               "local_opts" = local_opts )
  # epsa=rep(0,4)
  # epsb=rep(0,2)
  #   epsa[1:3]=log(alpha[1:3])
  #   epsa[4]=alpha[4]
  #   epsb[1]=log(beta[1])
  #   epsb[2]=beta[2]
  gradalpha=gradA_moy_cst2(xiFB,Y,K,I,C,cst)
  falpha=EqA_moy_cst2(xiFB,Y,K,I,cst)
  
  gradbeta=gradphi_FB_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_FB_gamma(gammaFB,Y,K,I)
  
  resuA=nloptr(x0 = alpha_cst, eval_f = falpha,eval_grad_f=gradalpha,lb=c(0,0,0),ub=c(100,100,100), opts = opts)
  
  resuphi=nloptr(x0 = beta, eval_f = fbeta,eval_grad_f=gradbeta,lb=c(0,-100), opts = opts)
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$solution,resuphi$solution)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}

EqA_moy_D2= function(xiFB,Y,K,I,D){
  EqA_FBM= function(alpha){
    sumn=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          for(ii in 1 : I){
            
            sumn= sumn + log(logitA_moy2(c(alpha[1],D,alpha[2:3]),i,ii,I,K,Y,c,n))* xiFB[i,ii,c,n+1]
            
            
          }
        }
      }
    }
    result= -(sumn)
    return(result)}
  return(EqA_FBM)
}
gradA_moy_D2= function(xiFB,Y,K,I,C,D){
  
  gradAM= function(alpha){
    
    sumo4=0
    sumn4=0
    sumo3=0
    sumn3=0
    sumo2=0
    sumn2=0
    sumo1=0
    sumn1=0
    
    for(c in 1:C){
      for(n in 1:N){
        #ii = in i= in-1
        for(i in 1:I){
          test1= exp(-(alpha[1]*((i)/I)+D*((Y[c,n])/K) + alpha[2]*((bij_moy(c,n,Y,K,C))/K) +alpha[3]))
          for(ii in 1 : I){
            
            sumn4= sumn4 + ((ii-1)+(1-I)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn3= sumn3 + ((ii-1)*((bij_moy(c,n,Y,K,C))/K)+(1-I)*((bij_moy(c,n,Y,K,C))/K)/(1+test1))*xiFB[i,ii,c,n+1]
            sumn2= sumn2 + ((ii-1)*((Y[c,n])/K)+(1-I)*((Y[c,n])/K)/(1+test1))* xiFB[i,ii,c,n+1]
            sumn1= sumn1 + ((ii-1)*((i)/I)+(1-I)*((i/I)/(1+test1)))*xiFB[i,ii,c,n+1]
            
            
          }
        }
      }
    }
    result1= -(sumn1)
    result3= -(sumn3)
    
    result4= -(sumn4)
    return(c(result1,result3,result4))}
  return(gradAM)
}
EtapeM_moy_D2 = function(gammaFB,xiFB,Y,pi,alpha,beta,K,D){
  DI=dim(Y)
  if(is.null(DI)){C=1
  N=length(Y)} else{
    C=DI[1]
    N= DI[2]}
  I= length(pi)
  alpha_D=c(alpha[1:2],alpha[4])
  local_opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                      "xtol_rel"  = 1.0e-18 )
  opts <- list("algorithm"="NLOPT_LD_AUGLAG","xtol_rel"=1.0e-18,"maxeval"   = 1000,
               "local_opts" = local_opts )
  # epsa=rep(0,4)
  # epsb=rep(0,2)
  #   epsa[1:3]=log(alpha[1:3])
  #   epsa[4]=alpha[4]
  #   epsb[1]=log(beta[1])
  #   epsb[2]=beta[2]
  gradalpha=gradA_moy_D2(xiFB,Y,K,I,C,D)
  falpha=EqA_moy_D2(xiFB,Y,K,I,D)
  
  gradbeta=gradphi_FB_gamma(gammaFB,Y,K,I)
  fbeta=Eqphi_FB_gamma(gammaFB,Y,K,I)
  
  resuA=nloptr(x0 = alpha_D, eval_f = falpha,eval_grad_f=gradalpha,lb=c(0,0,-100),ub=c(100,100,0.01), opts = opts)
  
  resuphi=nloptr(x0 = beta, eval_f = fbeta,eval_grad_f=gradbeta,lb=c(0,-100), opts = opts)
  
  
  gamma = - log( ((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1)
  if((sum((1:I)*gammaFB[,,1]) - C ==0)||(((I-1)*C/(sum((1:I)*gammaFB[,,1]) - C)) -1==0)){ gamma = -100}
  
  result=list(gamma,resuA$solution,resuphi$solution)
  names(result)= c('gamma','alpha','beta')
  
  
  return(result)
}




EM_moy2 = function(Y,etatA2,C,K,I,p){
  beta1=rep(0,2)
  alpha1=rep(0,4)
  beta1[1]=runif(1)*5
  beta1[2]=runif(1)*10-5
  gamma1=runif(1)*5
  
  alpha1[1:3]=runif(3)*5
  alpha1[4]=runif(1)*5-5
  
  recgamma=gamma1
  recbeta=beta1
  recalpha=alpha1
  recphi=funcphi(beta1,I,K)
  recpi=funcpi(gamma1,I)
  
  recA=funcA_moy2(alpha1,I,K,C)
  
  
  fullalpha=matrix(rep(0,4*(p+1)),nrow=4)
  fullbeta=matrix(rep(0,2*(p+1)),nrow=2)
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    
    
    fullalpha[,n]=recalpha
    fullbeta[,n]=recbeta
    fullgamma[,n]=recgamma
    
    
    RES= EtapeE_A2_BP(Y,etatA2,recpi,recA,recphi,C,I)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_moy2(rho,xi,Y,recpi,recalpha,recbeta,K)
    recalpha=RESu$alpha[1:4]
    tempA=RESu$alpha[5]
    recbeta=RESu$beta[1:2]
    tempphi=RESu$beta[3]
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    
    
    recpi=funcpi(recgamma,I)
    recA= funcA_moy2(recalpha,I,K,C)
    recphi=funcphi(recbeta,I,K)
    if(max(c(max(abs(recbeta -fullbeta[,n])) , max(abs(recalpha-fullalpha[,n])) ,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  fullalpha[,n]=recalpha
  fullbeta[,n]=recbeta
  fullgamma[,n]=recgamma
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  vrphi=funcphi(estimateurs$estimbeta,I,K)
  vrpi=funcpi(estimateurs$estimgamma,I)
  vrA=funcA_moy2(estimateurs$estimalpha,I,K,C)
  
  AB_FB=EtapeE_moy_LH(Y,vrpi,vrA,vrphi,C)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH(alpha_FB,beta_FB,vrphi,Y,C,N,I)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  Xstate=viterbi_moy(Y,estimateurs$estimalpha,estimateurs$estimbeta,estimateurs$estimgamma,C,K,I,N)
  YN1state=predictFB(alpha_FB,vrphi)
  result=list(ini,estimateurs,full,endn,Vraisemblance,YN1state,Xstate,Y,etatA2)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai','predYN','predBG','Y','Yvoisin')
  
  return(result)
}
EM_moy_S2 = function(Y,etatA2,C,K,I,p,S){
  beta1=rep(0,2)
  alpha1=rep(0,4)
  beta1[1]=runif(1)*5
  beta1[2]=runif(1)*10-5
  gamma1=runif(1)*(-5)
  
  alpha1[2:3]=runif(2)*5
  alpha1[1]=S
  alpha1[4]=runif(1)*5-7
  
  recgamma=gamma1
  recbeta=beta1
  recalpha=alpha1
  recphi=funcphi(beta1,I,K)
  recpi=funcpi(gamma1,I)
  
  recA=funcA_moy2(alpha1,I,K,C)
  
  
  fullalpha=matrix(rep(0,4*(p+1)),nrow=4)
  fullbeta=matrix(rep(0,2*(p+1)),nrow=2)
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    
    
    fullalpha[,n]=recalpha
    fullbeta[,n]=recbeta
    fullgamma[,n]=recgamma
    
    
    RES= EtapeE_A2_BP(Y,etatA2,recpi,recA,recphi,C,I)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_moy_S2(rho,xi,Y,recpi,recalpha,recbeta,K,S)
    recalpha=c(S,RESu$alpha[1:3])
    tempA=RESu$alpha[4]
    recbeta=RESu$beta[1:2]
    tempphi=RESu$beta[3]
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    
    
    recpi=funcpi(recgamma,I)
    recA= funcA_moy2(recalpha,I,K,C)
    recphi=funcphi(recbeta,I,K)
    if(max(c(max(abs(recbeta -fullbeta[,n])) , max(abs(recalpha-fullalpha[,n])) ,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  fullalpha[,n]=recalpha
  fullbeta[,n]=recbeta
  fullgamma[,n]=recgamma
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  vrphi=funcphi(estimateurs$estimbeta,I,K)
  vrpi=funcpi(estimateurs$estimgamma,I)
  vrA=funcA_moy2(estimateurs$estimalpha,I,K,C)
  
  AB_FB=EtapeE_moy_LH(Y,vrpi,vrA,vrphi,C)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH(alpha_FB,beta_FB,vrphi,Y,C,N,I)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  result=list(ini,estimateurs,full,endn,Vraisemblance)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai')
  
  return(result)
}
EM_moy_V2 = function(Y,etatA2,C,K,I,p,V){
  beta1=rep(0,2)
  alpha1=rep(0,4)
  beta1[1]=runif(1)*5
  beta1[2]=runif(1)*10-5
  gamma1=runif(1)*(-5)
  
  alpha1[1:2]=runif(2)*5
  alpha1[3]=V
  alpha1[4]=runif(1)*5-7
  
  recgamma=gamma1
  recbeta=beta1
  recalpha=alpha1
  recphi=funcphi(beta1,I,K)
  recpi=funcpi(gamma1,I)
  
  recA=funcA_moy2(alpha1,I,K,C)
  
  
  fullalpha=matrix(rep(0,4*(p+1)),nrow=4)
  fullbeta=matrix(rep(0,2*(p+1)),nrow=2)
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    
    
    fullalpha[,n]=recalpha
    fullbeta[,n]=recbeta
    fullgamma[,n]=recgamma
    
    
    RES= EtapeE_A2_BP(Y,etatA2,recpi,recA,recphi,C,I)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_moy_V2(rho,xi,Y,recpi,recalpha,recbeta,K,V)
    recalpha=c(RESu$alpha[1:2],V,RESu$alpha[3])
    tempA=RESu$alpha[4]
    recbeta=RESu$beta[1:2]
    tempphi=RESu$beta[3]
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    
    
    recpi=funcpi(recgamma,I)
    recA= funcA_moy2(recalpha,I,K,C)
    recphi=funcphi(recbeta,I,K)
    if(max(c(max(abs(recbeta -fullbeta[,n])) , max(abs(recalpha-fullalpha[,n])) ,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  fullalpha[,n]=recalpha
  fullbeta[,n]=recbeta
  fullgamma[,n]=recgamma
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  vrphi=funcphi(estimateurs$estimbeta,I,K)
  vrpi=funcpi(estimateurs$estimgamma,I)
  vrA=funcA_moy2(estimateurs$estimalpha,I,K,C)
  
  AB_FB=EtapeE_moy_LH(Y,vrpi,vrA,vrphi,C)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH(alpha_FB,beta_FB,vrphi,Y,C,N,I)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  result=list(ini,estimateurs,full,endn,Vraisemblance)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai')
  
  return(result)
}
EM_moy_cst2 = function(Y,etatA2,C,K,I,p,cst){
  beta1=rep(0,2)
  alpha1=rep(0,4)
  beta1[1]=runif(1)*5
  beta1[2]=runif(1)*10-5
  gamma1=runif(1)*(-5)
  
  alpha1[1:2]=runif(2)*5
  alpha1[3]=V
  alpha1[4]=runif(1)*5-7
  
  recgamma=gamma1
  recbeta=beta1
  recalpha=alpha1
  recphi=funcphi(beta1,I,K)
  recpi=funcpi(gamma1,I)
  
  recA=funcA_moy2(alpha1,I,K,C)
  
  
  fullalpha=matrix(rep(0,4*(p+1)),nrow=4)
  fullbeta=matrix(rep(0,2*(p+1)),nrow=2)
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    
    
    fullalpha[,n]=recalpha
    fullbeta[,n]=recbeta
    fullgamma[,n]=recgamma
    
    
    RES= EtapeE_A2_BP(Y,etatA2,recpi,recA,recphi,C,I)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_moy_cst2(rho,xi,Y,recpi,recalpha,recbeta,K,cst)
    recalpha=c(RESu$alpha[1:3],cst)
    tempA=RESu$alpha[5]
    recbeta=RESu$beta[1:2]
    tempphi=RESu$beta[3]
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    
    
    recpi=funcpi(recgamma,I)
    recA= funcA_moy2(recalpha,I,K,C)
    recphi=funcphi(recbeta,I,K)
    if(max(c(max(abs(recbeta -fullbeta[,n])) , max(abs(recalpha-fullalpha[,n])) ,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  fullalpha[,n]=recalpha
  fullbeta[,n]=recbeta
  fullgamma[,n]=recgamma
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  vrphi=funcphi(estimateurs$estimbeta,I,K)
  vrpi=funcpi(estimateurs$estimgamma,I)
  vrA=funcA_moy2(estimateurs$estimalpha,I,K,C)
  
  AB_FB=EtapeE_moy_LH(Y,vrpi,vrA,vrphi,C)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH(alpha_FB,beta_FB,vrphi,Y,C,N,I)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  result=list(ini,estimateurs,full,endn,Vraisemblance)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai')
  
  return(result)
}
EM_moy_D2 = function(Y,etatA2,C,K,I,p,D){
  beta1=rep(0,2)
  alpha1=rep(0,4)
  beta1[1]=runif(1)*5
  beta1[2]=runif(1)*10-5
  gamma1=runif(1)*(-5)
  
  alpha1[1]=runif(1)*5
  alpha1[3]=runif(1)*5
  alpha1[4]=runif(1)*5-7
  alpha1[2]=D
  
  recgamma=gamma1
  recbeta=beta1
  recalpha=alpha1
  recphi=funcphi(beta1,I,K)
  recpi=funcpi(gamma1,I)
  
  recA=funcA_moy2(alpha1,I,K,C)
  
  
  fullalpha=matrix(rep(0,4*(p+1)),nrow=4)
  fullbeta=matrix(rep(0,2*(p+1)),nrow=2)
  fullgamma=matrix(rep(0,(p+1)),nrow=1)
  
  
  
  end=0
  for(n in 1:p){
    
    
    
    
    fullalpha[,n]=recalpha
    fullbeta[,n]=recbeta
    fullgamma[,n]=recgamma
    
    
    RES= EtapeE_A2_BP(Y,etatA2,recpi,recA,recphi,C,I)
    
    
    value = RES[[3]]
    
    temp= value
    rho=RES[[1]]
    xi=RES[[2]]
    
    
    #test de l'élément d'avant ! alpha_it beta_it
    RESu=EtapeM_moy_D2(rho,xi,Y,recpi,recalpha,recbeta,K,D)
    recalpha=c(RESu$alpha[1],D,RESu$alpha[2:3])
    tempA=RESu$alpha[4]
    recbeta=RESu$beta[1:2]
    tempphi=RESu$beta[3]
    print("n =")
    print(n)
    recgamma=RESu$gamma
    
    
    
    recpi=funcpi(recgamma,I)
    recA= funcA_moy2(recalpha,I,K,C)
    recphi=funcphi(recbeta,I,K)
    if(max(c(max(abs(recbeta -fullbeta[,n])) , max(abs(recalpha-fullalpha[,n])) ,max(abs(recgamma-fullgamma[,n]))))<0.01  ){
      print("EM converged ! ")
      break 
    }
    
  }
  
  endn=n
  
  
  
  fullalpha[,n]=recalpha
  fullbeta[,n]=recbeta
  fullgamma[,n]=recgamma
  ini=list(gamma1,alpha1,beta1)
  names(ini)= c('gammainit','alphainit','betainit')
  
  full=list(fullgamma,fullalpha,fullbeta)
  names(full)=c('fullgamma','fullalpha','fullbeta')
  estimateurs= list(recgamma,recalpha,recbeta)
  names(estimateurs)=c('estimgamma','estimalpha','estimbeta')
  vrphi=funcphi(estimateurs$estimbeta,I,K)
  vrpi=funcpi(estimateurs$estimgamma,I)
  vrA=funcA_moy2(estimateurs$estimalpha,I,K,C)
  
  AB_FB=EtapeE_moy_LH(Y,vrpi,vrA,vrphi,C)
  alpha_FB=AB_FB[[1]]
  beta_FB=AB_FB[[2]]
  LLH=BP_LLH(alpha_FB,beta_FB,vrphi,Y,C,N,I)
  Vraisemblance=rep(0,2)
  Vraisemblance[1]=LLH[[1]]
  Vraisemblance[2]=LLH[[2]]
  names(Vraisemblance)=c('best','moyenne')
  result=list(ini,estimateurs,full,endn,Vraisemblance)
  names(result)=c('paraminitaux','paramestim', 'full','nbriter','Vrai')
  
  return(result)
}


viterbi_moy2=function(Y,alpha,beta,gamma,C,K,I,N){
  A2=funcA_moy2(alpha,I,K,C)
  phi=funcphi(beta,I,K)
  pi=funcpi(gamma,I)
  XX=array(rep(0,C*(N+1)),dim=c(C,N+1))
  delta=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  deltaind=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  
  
  
  for(c in 1 :C){
    delta[,c,1]=log(pi)
    etatA2=field_to_vect_moy(c,Y,K)
    for(n in 2:(N+1)){
      
      mad=delta[,c,n-1]+log(phi[,Y[c,n-1]])+log(A2[,,etatA2[n-1]])   
      delta[,c,n]=apply(mad,2,max)
      for(i in 1:I){
        deltaind[i,c,n-1]=which.max(mad[,i])
      }
      
    }
    for(i in 1:I){
      deltaind[1,c,N+1]=which.max(delta[,c,N+1])
    }
    
    XX[c,N+1]=deltaind[1,c,N+1]
    for(n in N:1){
      XX[c,n]=deltaind[XX[c,n+1],c,n]
      
    }
    
    
    
  }
  return(XX)
}
viterbi_moy=function(Y,alpha,beta,gamma,C,K,I,N){
  A2=funcA_moy(alpha,I,K,C)
  phi=funcphi(beta,I,K)
  pi=funcpi(gamma,I)
  XX=array(rep(0,C*(N+1)),dim=c(C,N+1))
  delta=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  deltaind=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  
  
  
  for(c in 1 :C){
    delta[,c,1]=log(pi)
    etatA2=field_to_vect_moy(c,Y,K)
    for(n in 2:(N+1)){
      
      mad=delta[,c,n-1]+log(phi[,Y[c,n-1]])+log(A2[,,etatA2[n-1]])   
      delta[,c,n]=apply(mad,2,max)
      for(i in 1:I){
        deltaind[i,c,n-1]=which.max(mad[,i])
      }
      
    }
    for(i in 1:I){
      deltaind[1,c,N+1]=which.max(delta[,c,N+1])
    }
    
    XX[c,N+1]=deltaind[1,c,N+1]
    for(n in N:1){
      XX[c,n]=deltaind[XX[c,n+1],c,n]
      
    }
    
    
    
  }
  return(XX)
}


viterbi_alp=function(Y,alpha,beta,gamma,C,K,I,N){
  A2=funcA_2(alpha,I,K,C)
  phi=funcphi(beta,I,K)
  pi=funcpi(gamma,I)
  XX=array(rep(0,C*(N+1)),dim=c(C,N+1))
  delta=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  deltaind=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  
  
  
  for(c in 1 :C){
    delta[,c,1]=log(pi)
    etatA2=field_to_vect_voisin_inter(c,Y,K)
    for(n in 2:(N+1)){
      
      mad=delta[,c,n-1]+log(phi[,Y[c,n-1]])+log(A2[,,etatA2[n-1]])   
      delta[,c,n]=apply(mad,2,max)
      for(i in 1:I){
        deltaind[i,c,n-1]=which.max(mad[,i])
      }
      
    }
    for(i in 1:I){
      deltaind[1,c,N+1]=which.max(delta[,c,N+1])
    }
    
    XX[c,N+1]=deltaind[1,c,N+1]
    for(n in N:1){
      XX[c,n]=deltaind[XX[c,n+1],c,n]
      
    }
    
    
    
  }
  return(XX)
}

selection_model_multicoeur_S = function(C,K,I,p,alpha,beta,gamma,N,nbr,S){
  phi=funcphi(beta,I,K)
  pi=funcpi(gamma,I)
  
  Y=NULL
  X=NULL
  
  A2=funcA_moy(alpha,I,K,C)
  
  
  
  
  tempo= bij_moy(1,1,rep(K,C),K,C)
  
  
  X= array(rep(0,C*(N+1)),dim=c(C,N+1))
  Y= array(rep(0,C*(N+1)),dim=c(C,N))
  #creation simulation avec A2 phi pi
  for(i in 1:C){
    f=runif(1)
    if(f<pi[1]){  X[i,1]=1
    } else if(f <pi[1]+ pi[2]){  X[i,1]=2
    } else if(f <pi[1]+ pi[2]+ pi[3]){  X[i,1]=3
    } else if(f <pi[1]+ pi[2]+pi[3]+pi[4]){  X[i,1]=4
    } else {  X[i,1]=5}  
  }
  
  for(i in 1:N){
    for(c in 1:C){
      l= runif(1)
      
      if(l<phi[X[c,i],1]){ Y[c,i]=1
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]){ Y[c,i]=2
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+ phi[X[c,i],3]){ Y[c,i]=3
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+phi[X[c,i],3]+phi[X[c,i],4]){ Y[c,i]=4
      } else { Y[c,i]=5}  
    }
    for(c in 1:C){
      f=field_to_vect_moy(c,Y[,i],K)
      R=  runif(1)
      
      
      if(R<A2[X[c,i],1,f]){ X[c,i+1]=1
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]){ X[c,i+1]=2
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+ A2[X[c,i],3,f]){ X[c,i+1]=3
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+A2[X[c,i],3,f]+A2[X[c,i],4,f]){ X[c,i+1]=4
      } else{ X[c,i+1]=5}  
      
    }
    
    
  }
  
  
  etatA2=Y
  for(c in 1:C){
    etatA2[c,]=field_to_vect_moy(c,Y,K)}
  
  difsim=foreach(i=1:nbr) %dopar% EM_moy(Y,etatA2,C,K,I,p)
  
  best=-1000000
  for(i in 1:nbr){
    if(best<difsim[[i]]$Vrai[2]){result=difsim[[i]]
    best=difsim[[i]]$Vrai[2]}
    
  }
  
  difsim=foreach(i=1:nbr) %dopar% EM_moy_S(Y,etatA2,C,K,I,p,0)
  
  best=-1000000
  for(i in 1:nbr){
    if(best<difsim[[i]]$Vrai[2]){results=difsim[[i]]
    best=difsim[[i]]$Vrai[2]}
    
  }
  
  Bic= log(N*C)*7 - 2*result$Vrai
  Bics= log(N*C)*6 - 2*results$Vrai
  
  ret=list(result$paramestim,result$Vrai,Bic,results$paramestim,Bics,results$Vrai)
  names(ret)=c('para','Vrai','Bic','paras', 'Bics','Vrais')
  
  return(ret)
}

selection_model_multicoeur_V = function(C,K,I,p,alpha,beta,gamma,N,nbr,V){
  phi=funcphi(beta,I,K)
  pi=funcpi(gamma,I)
  
  Y=NULL
  X=NULL
  
  A2=funcA_moy(alpha,I,K,C)
  
  
  
  
  tempo= bij_moy(1,1,rep(K,C),K,C)
  
  
  X= array(rep(0,C*(N+1)),dim=c(C,N+1))
  Y= array(rep(0,C*(N+1)),dim=c(C,N))
  #creation simulation avec A2 phi pi
  for(i in 1:C){
    f=runif(1)
    if(f<pi[1]){  X[i,1]=1
    } else if(f <pi[1]+ pi[2]){  X[i,1]=2
    } else if(f <pi[1]+ pi[2]+ pi[3]){  X[i,1]=3
    } else if(f <pi[1]+ pi[2]+pi[3]+pi[4]){  X[i,1]=4
    } else {  X[i,1]=5}  
  }
  
  for(i in 1:N){
    for(c in 1:C){
      l= runif(1)
      
      if(l<phi[X[c,i],1]){ Y[c,i]=1
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]){ Y[c,i]=2
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+ phi[X[c,i],3]){ Y[c,i]=3
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+phi[X[c,i],3]+phi[X[c,i],4]){ Y[c,i]=4
      } else { Y[c,i]=5}  
    }
    for(c in 1:C){
      f=field_to_vect_moy(c,Y[,i],K)
      R=  runif(1)
      
      
      if(R<A2[X[c,i],1,f]){ X[c,i+1]=1
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]){ X[c,i+1]=2
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+ A2[X[c,i],3,f]){ X[c,i+1]=3
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+A2[X[c,i],3,f]+A2[X[c,i],4,f]){ X[c,i+1]=4
      } else{ X[c,i+1]=5}  
      
    }
    
    
  }
  
  
  etatA2=Y
  for(c in 1:C){
    etatA2[c,]=field_to_vect_moy(c,Y,K)}
  
  difsim=foreach(i=1:nbr) %dopar% EM_moy(Y,etatA2,C,K,I,p)
  
  best=-1000000
  for(i in 1:nbr){
    if(best<difsim[[i]]$Vrai[2]){result=difsim[[i]]
    best=difsim[[i]]$Vrai[2]}
    
  }
  
  difsim=foreach(i=1:nbr) %dopar% EM_moy_V(Y,etatA2,C,K,I,p,0)
  
  best=-1000000
  for(i in 1:nbr){
    if(best<difsim[[i]]$Vrai[2]){results=difsim[[i]]
    best=difsim[[i]]$Vrai[2]}
    
  }
  
  Bic= log(N*C)*7 - 2*result$Vrai
  Bics= log(N*C)*6 - 2*results$Vrai
  
  ret=list(result$paramestim,result$Vrai,Bic,results$paramestim,Bics,results$Vrai)
  names(ret)=c('para','Vrai','Bic','paras', 'Bics','Vrais')
  
  return(ret)
}

selection_model_multicoeur_moy_inter = function(C,K,I,p,alpha,beta,gamma,N,nbr){
  phi=funcphi(beta,I,K)
  pi=funcpi(gamma,I)
  
  Y=NULL
  X=NULL
  
  A2=funcA_moy(alpha,I,K,C)
  
  
  
  
  tempo= bij_moy(1,1,rep(K,C),K,C)
  
  
  X= array(rep(0,C*(N+1)),dim=c(C,N+1))
  Y= array(rep(0,C*(N+1)),dim=c(C,N))
  #creation simulation avec A2 phi pi
  for(i in 1:C){
    f=runif(1)
    if(f<pi[1]){  X[i,1]=1
    } else if(f <pi[1]+ pi[2]){  X[i,1]=2
    } else if(f <pi[1]+ pi[2]+ pi[3]){  X[i,1]=3
    } else if(f <pi[1]+ pi[2]+pi[3]+pi[4]){  X[i,1]=4
    } else {  X[i,1]=5}  
  }
  
  for(i in 1:N){
    for(c in 1:C){
      l= runif(1)
      
      if(l<phi[X[c,i],1]){ Y[c,i]=1
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]){ Y[c,i]=2
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+ phi[X[c,i],3]){ Y[c,i]=3
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+phi[X[c,i],3]+phi[X[c,i],4]){ Y[c,i]=4
      } else { Y[c,i]=5}  
    }
    for(c in 1:C){
      f=field_to_vect_moy(c,Y[,i],K)
      R=  runif(1)
      
      
      if(R<A2[X[c,i],1,f]){ X[c,i+1]=1
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]){ X[c,i+1]=2
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+ A2[X[c,i],3,f]){ X[c,i+1]=3
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+A2[X[c,i],3,f]+A2[X[c,i],4,f]){ X[c,i+1]=4
      } else{ X[c,i+1]=5}  
      
    }
    
    
  }
  
  
  etatAmoymoy=Y
  etatAmoyinter=Y
  
  for(c in 1:C){
    etatAmoymoy[c,]=field_to_vect_moy(c,Y,K)
    etatAmoyinter[c,]=field_to_vect_voisin_inter(c,Y,K)
    }

  
  A21=funcA_2(alpha,I,K,C)
  X1= array(rep(0,C*(N+1)),dim=c(C,N+1))
  Y1= array(rep(0,C*(N+1)),dim=c(C,N))
  #creation simulation avec A2 phi pi
  for(i in 1:C){
    f=runif(1)
    if(f<pi[1]){  X1[i,1]=1
    } else if(f <pi[1]+ pi[2]){  X1[i,1]=2
    } else if(f <pi[1]+ pi[2]+ pi[3]){  X1[i,1]=3
    } else if(f <pi[1]+ pi[2]+pi[3]+pi[4]){  X1[i,1]=4
    } else {  X1[i,1]=5}  
  }
  
  for(i in 1:N){
    for(c in 1:C){
      l= runif(1)
      
      if(l<phi[X1[c,i],1]){ Y1[c,i]=1
      } else if(l <phi[X1[c,i],1]+ phi[X1[c,i],2]){ Y1[c,i]=2
      } else if(l <phi[X1[c,i],1]+ phi[X1[c,i],2]+ phi[X1[c,i],3]){ Y1[c,i]=3
      } else if(l <phi[X1[c,i],1]+ phi[X1[c,i],2]+phi[X1[c,i],3]+phi[X1[c,i],4]){ Y1[c,i]=4
      } else { Y1[c,i]=5}  
    }
    for(c in 1:C){
      f=field_to_vect_voisin_inter(c,Y[,i],K)
      R=  runif(1)
      
      
      if(R<A21[X1[c,i],1,f]){ X1[c,i+1]=1
      } else if(R <A21[X1[c,i],1,f]+ A21[X1[c,i],2,f]){ X1[c,i+1]=2
      } else if(R <A21[X1[c,i],1,f]+ A21[X1[c,i],2,f]+ A21[X1[c,i],3,f]){ X1[c,i+1]=3
      } else if(R <A21[X1[c,i],1,f]+ A21[X1[c,i],2,f]+A21[X1[c,i],3,f]+A21[X1[c,i],4,f]){ X1[c,i+1]=4
      } else{ X1[c,i+1]=5}  
      
    }
    
    
  }
  
  
  etatAinterinter=Y
  etatAintermoy=Y
  for(c in 1:C){
    etatAinterinter[c,]=field_to_vect_voisin_inter(c,Y1,K)
    etatAintermoy[c,]=field_to_vect_moy(c,Y1,K)}
  
  difsimmoy_moy=foreach(i=1:nbr) %dopar% EM_moy(Y,etatAmoymoy,C,K,I,p)
  difsimmoy_inter= foreach(i =1:nbr) %dopar% EM_alp(Y,etatAmoyinter,C,K,I,p)
  
  difsiminter_inter= foreach(i=1:nbr) %dopar% EM_alp(Y1,etatAinterinter,C,K,I,p)
  difsiminter_moy= foreach(i=1:nbr) %dopar%EM_moy(Y1,etatAintermoy,C,K,I,p)
  best=-1000000
  best2=-1000000
  best3=-1000000
  best4=-1000000
  for(i in 1:nbr){
    if(best<difsimmoy_moy[[i]]$Vrai[2]){resultmm=difsimmoy_moy[[i]]
    best=difsimmoy_moy[[i]]$Vrai[2]}
    if(best2<difsimmoy_inter[[i]]$Vrai[2]){resultmi=difsimmoy_inter[[i]]
    best2=difsimmoy_inter[[i]]$Vrai[2]}
    if(best3< difsiminter_inter[[i]]$Vrai[2]){resultii=difsiminter_inter[[i]]
    best3=difsiminter_inter[[i]]$Vrai[2]}
    if(best4< difsiminter_moy[[i]]$Vrai[2]){resultim=difsiminter_moy[[i]]
    best4= difsiminter_moy[[i]]$Vrai[2]}
    
  }
  
  
  
  ret=list(Y,X,Y1,X1,resultmm,resultmi,resultim,resultii,best-best2,best3-best4)
    
  names(ret)=c('Ymoy','Xmoy','Yinter','Xinter','resmm','resmi','resim','resii','selinter','selmoy')
  
  return(ret)
}




sim_EM_moy_multicoeur = function(C,K,I,p,alpha,beta,gamma,N,nbr){
  phi=funcphi(beta,I,K)
  pi=funcpi(gamma,I)
  
  Y=NULL
  X=NULL
  
  A2=funcA_moy(alpha,I,K,C)
  
  
  
  
  tempo= bij_moy(1,1,rep(K,C),K,C)
  
  
  X= array(rep(0,C*(N+1)),dim=c(C,N+1))
  Y= array(rep(0,C*(N+1)),dim=c(C,N))
  #creation simulation avec A2 phi pi
  for(i in 1:C){
    f=runif(1)
    if(f<pi[1]){  X[i,1]=1
    } else if(f <pi[1]+ pi[2]){  X[i,1]=2
    } else if(f <pi[1]+ pi[2]+ pi[3]){  X[i,1]=3
    } else if(f <pi[1]+ pi[2]+pi[3]+pi[4]){  X[i,1]=4
    } else {  X[i,1]=5}  
  }
  
  for(i in 1:N){
    for(c in 1:C){
      l= runif(1)
      
      if(l<phi[X[c,i],1]){ Y[c,i]=1
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]){ Y[c,i]=2
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+ phi[X[c,i],3]){ Y[c,i]=3
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+phi[X[c,i],3]+phi[X[c,i],4]){ Y[c,i]=4
      } else { Y[c,i]=5}  
    }
    for(c in 1:C){
      f=field_to_vect_moy(c,Y[,i],K)
      R=  runif(1)
      
      
      if(R<A2[X[c,i],1,f]){ X[c,i+1]=1
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]){ X[c,i+1]=2
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+ A2[X[c,i],3,f]){ X[c,i+1]=3
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+A2[X[c,i],3,f]+A2[X[c,i],4,f]){ X[c,i+1]=4
      } else{ X[c,i+1]=5}  
      
    }
    
    
  }
  
  
  etatA2=Y
  for(c in 1:C){
    etatA2[c,]=field_to_vect_moy(c,Y,K)}
  
  
  difsim=foreach(i=1:nbr) %dopar% EM_moy(Y,etatA2,C,K,I,p)
  
  best=-1000000
  for(i in 1:nbr){
    if(best<difsim[[i]]$Vrai[2]){result=difsim[[i]]
    best=difsim[[i]]$Vrai[2]}
    
  }
  result$X=X
  return(result)
}

sim_EM_ZI_multicoeur = function(C,K,I,p,alpha,beta,gamma,N,nbr){
  phi=ZIphi(beta,I,K)
  pi=funcpi(gamma,I)
  
  Y=NULL
  X=NULL
  
  A2=funcA_moy(alpha,I,K,C)
  
  
  
  
  tempo= bij_moy(1,1,rep(K,C),K,C)
  
  
  X= array(rep(0,C*(N+1)),dim=c(C,N+1))
  Y= array(rep(0,C*(N+1)),dim=c(C,N))
  #creation simulation avec A2 phi pi
  for(i in 1:C){
    f=runif(1)
    if(f<pi[1]){  X[i,1]=1
    } else if(f <pi[1]+ pi[2]){  X[i,1]=2
    } else if(f <pi[1]+ pi[2]+ pi[3]){  X[i,1]=3
    } else if(f <pi[1]+ pi[2]+pi[3]+pi[4]){  X[i,1]=4
    } else {  X[i,1]=5}  
  }
  
  for(i in 1:N){
    for(c in 1:C){
      l= runif(1)
      
      if(l<phi[X[c,i],1]){ Y[c,i]=1
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]){ Y[c,i]=2
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+ phi[X[c,i],3]){ Y[c,i]=3
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+phi[X[c,i],3]+phi[X[c,i],4]){ Y[c,i]=4
      } else { Y[c,i]=5}  
    }
    for(c in 1:C){
      f=field_to_vect_moy(c,Y[,i],K)
      R=  runif(1)
      
      
      if(R<A2[X[c,i],1,f]){ X[c,i+1]=1
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]){ X[c,i+1]=2
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+ A2[X[c,i],3,f]){ X[c,i+1]=3
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+A2[X[c,i],3,f]+A2[X[c,i],4,f]){ X[c,i+1]=4
      } else{ X[c,i+1]=5}  
      
    }
    
    
  }
  
  
  etatA2=Y
  for(c in 1:C){
    etatA2[c,]=field_to_vect_moy(c,Y,K)}
  
  
  difsim=foreach(i=1:nbr) %dopar% EM_moy_ZI(Y,etatA2,C,K,I,p)
  
  best=-1000000
  for(i in 1:nbr){
    if(best<difsim[[i]]$Vrai[2]){result=difsim[[i]]
    best=difsim[[i]]$Vrai[2]}
    
  }
  result$X=X
  return(result)
}



sim_EM_moy_multicoeur2 = function(C,K,I,p,alpha,beta,gamma,N,nbr){
  phi=funcphi(beta,I,K)
  pi=funcpi(gamma,I)
  
  Y=NULL
  X=NULL
  
  A2=funcA_moy2(alpha,I,K,C)
  
  
  
  
  tempo= bij_moy(1,1,rep(K,C),K,C)
  
  
  X= array(rep(0,C*(N+1)),dim=c(C,N+1))
  Y= array(rep(0,C*(N+1)),dim=c(C,N))
  #creation simulation avec A2 phi pi
  for(i in 1:C){
    f=runif(1)
    if(f<pi[1]){  X[i,1]=1
    } else if(f <pi[1]+ pi[2]){  X[i,1]=2
    } else if(f <pi[1]+ pi[2]+ pi[3]){  X[i,1]=3
    } else if(f <pi[1]+ pi[2]+pi[3]+pi[4]){  X[i,1]=4
    } else {  X[i,1]=5}  
  }
  
  for(i in 1:N){
    for(c in 1:C){
      l= runif(1)
      
      if(l<phi[X[c,i],1]){ Y[c,i]=1
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]){ Y[c,i]=2
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+ phi[X[c,i],3]){ Y[c,i]=3
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+phi[X[c,i],3]+phi[X[c,i],4]){ Y[c,i]=4
      } else { Y[c,i]=5}  
    }
    for(c in 1:C){
      f=field_to_vect_moy(c,Y[,i],K)
      R=  runif(1)
      
      
      if(R<A2[X[c,i],1,f]){ X[c,i+1]=1
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]){ X[c,i+1]=2
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+ A2[X[c,i],3,f]){ X[c,i+1]=3
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+A2[X[c,i],3,f]+A2[X[c,i],4,f]){ X[c,i+1]=4
      } else{ X[c,i+1]=5}  
      
    }
    
    
  }
  
  
  etatA2=Y
  for(c in 1:C){
    etatA2[c,]=field_to_vect_moy(c,Y,K)}
  
  
  difsim=foreach(i=1:nbr) %dopar% EM_moy2(Y,etatA2,C,K,I,p)
  
  best=-1000000
  for(i in 1:nbr){
    if(best<difsim[[nbr]]$Vrai[2]){result=difsim[[nbr]]
    best=difsim[[nbr]]$Vrai[2]}
    
  }
  result=list(result,Y,X)
  
  return(result)
}

Pred_viterbi_EM_moy_multicoeur = function(C,K,I,p,alpha,beta,gamma,nbr,YY,XX){
  o=dim(Y)
  N=o[2]
  C=o[1]
  Y=YY[,(1:(N-1))]
  
  etatA2=Y
  for(c in 1:C){
    etatA2[c,]=field_to_vect_moy(c,Y,K)}
  
  
  difsim=foreach(i=1:nbr) %dopar% EM_moy(Y,etatA2,C,K,I,p)
  
  best=-1000000
  for(i in 1:nbr){
    if(best<difsim[[i]]$Vrai[2]){result=difsim[[i]]
    best=difsim[[i]]$Vrai[2]}
    
  }
  PredictionYN=length(which(result$predYN==YY[,N]))/C
  PredictionX= length(which(XX[,(1:(N))]==result$predBG))/((N)*C)
  results=list(result,PredictionYN,PredictionX)
  names(results)=c('res','YN','X')
  return(results)
}


sim_EM_moy_multicoeur_S_VD = function(C,K,I,p,alpha,beta,gamma,N,nbr){
  phi=funcphi(beta,I,K)
  pi=funcpi(gamma,I)
  
  Y=NULL
  X=NULL
  S=alpha[1]
  V=alpha[3]
  D=alpha[2]
  A2=funcA_moy(alpha,I,K,C)
  
  
  
  
  tempo= bij_moy(1,1,rep(K,C),K,C)
  
  
  X= array(rep(0,C*(N+1)),dim=c(C,N+1))
  Y= array(rep(0,C*(N+1)),dim=c(C,N))
  #creation simulation avec A2 phi pi
  for(i in 1:C){
    f=runif(1)
    if(f<pi[1]){  X[i,1]=1
    } else if(f <pi[1]+ pi[2]){  X[i,1]=2
    } else if(f <pi[1]+ pi[2]+ pi[3]){  X[i,1]=3
    } else if(f <pi[1]+ pi[2]+pi[3]+pi[4]){  X[i,1]=4
    } else {  X[i,1]=5}  
  }
  
  for(i in 1:N){
    for(c in 1:C){
      l= runif(1)
      
      if(l<phi[X[c,i],1]){ Y[c,i]=1
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]){ Y[c,i]=2
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+ phi[X[c,i],3]){ Y[c,i]=3
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+phi[X[c,i],3]+phi[X[c,i],4]){ Y[c,i]=4
      } else { Y[c,i]=5}  
    }
    for(c in 1:C){
      f=field_to_vect_moy(c,Y[,i],K)
      R=  runif(1)
      
      
      if(R<A2[X[c,i],1,f]){ X[c,i+1]=1
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]){ X[c,i+1]=2
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+ A2[X[c,i],3,f]){ X[c,i+1]=3
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+A2[X[c,i],3,f]+A2[X[c,i],4,f]){ X[c,i+1]=4
      } else{ X[c,i+1]=5}  
      
    }
    
    
  }
  
  
  etatA2=Y
  for(c in 1:C){
    etatA2[c,]=field_to_vect_moy(c,Y,K)}
  
  
  difsim=foreach(i=1:nbr) %dopar% EM_moy(Y,etatA2,C,K,I,p)
  difsimS=foreach(i=1:nbr) %dopar% EM_moy_S(Y,etatA2,C,K,I,p,S)
  difsimDV= foreach(i=1:nbr) %dopar% EM_moy_DV(Y,etatA2,C,K,I,p,S)
  result=list(difsim[[1]],difsimS[[1]],difsimDV[[1]],X)
  names(result)=c('fullEM','EM_S','EM_DV','X')
  best=-1000000
  bests=-1000000
  bestdv=-1000000
  for(i in 2:nbr){
    if(best<difsim[[nbr]]$Vrai[2]){result[[1]]=difsim[[nbr]]
    best=difsim[[nbr]]$Vrai[2]}
    if(bests<difsimS[[nbr]]$Vrai[2]){result[[2]]=difsimS[[nbr]]
    bests=difsimS[[nbr]]$Vrai[2]}
    if(bests<difsimDV[[nbr]]$Vrai[2]){result[[3]]=difsimDV[[nbr]]
    bestdv=difsimDV[[nbr]]$Vrai[2]}
    
  }
  
  return(result)
}
sim_EM_moy_multicoeur_S_VD2 = function(C,K,I,p,alpha,beta,gamma,N,nbr){
  phi=funcphi(beta,I,K)
  pi=funcpi(gamma,I)
  
  Y=NULL
  X=NULL
  S=alpha[1]
  V=alpha[3]
  D=alpha[2]
  A2=funcA_moy2(alpha,I,K,C)
  
  
  
  
  tempo= bij_moy(1,1,rep(K,C),K,C)
  
  
  X= array(rep(0,C*(N+1)),dim=c(C,N+1))
  Y= array(rep(0,C*(N+1)),dim=c(C,N))
  #creation simulation avec A2 phi pi
  for(i in 1:C){
    f=runif(1)
    if(f<pi[1]){  X[i,1]=1
    } else if(f <pi[1]+ pi[2]){  X[i,1]=2
    } else if(f <pi[1]+ pi[2]+ pi[3]){  X[i,1]=3
    } else if(f <pi[1]+ pi[2]+pi[3]+pi[4]){  X[i,1]=4
    } else {  X[i,1]=5}  
  }
  
  for(i in 1:N){
    for(c in 1:C){
      l= runif(1)
      
      if(l<phi[X[c,i],1]){ Y[c,i]=1
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]){ Y[c,i]=2
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+ phi[X[c,i],3]){ Y[c,i]=3
      } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+phi[X[c,i],3]+phi[X[c,i],4]){ Y[c,i]=4
      } else { Y[c,i]=5}  
    }
    for(c in 1:C){
      f=field_to_vect_moy(c,Y[,i],K)
      R=  runif(1)
      
      
      if(R<A2[X[c,i],1,f]){ X[c,i+1]=1
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]){ X[c,i+1]=2
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+ A2[X[c,i],3,f]){ X[c,i+1]=3
      } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+A2[X[c,i],3,f]+A2[X[c,i],4,f]){ X[c,i+1]=4
      } else{ X[c,i+1]=5}  
      
    }
    
    
  }
  
  
  etatA2=Y
  for(c in 1:C){
    etatA2[c,]=field_to_vect_moy(c,Y,K)}
  
  
  difsim=foreach(i=1:nbr) %dopar% EM_moy2(Y,etatA2,C,K,I,p)
  difsimS=foreach(i=1:nbr) %dopar% EM_moy_S2(Y,etatA2,C,K,I,p,S)
  difsimDV= foreach(i=1:nbr) %dopar% EM_moy_DV2(Y,etatA2,C,K,I,p,S)
  result=list(difsim[[1]],difsimS[[1]],difsimDV[[1]],X)
  names(result)=c('fullEM','EM_S','EM_DV','X')
  best=-1000000
  bests=-1000000
  bestdv=-1000000
  for(i in 2:nbr){
    if(best<difsim[[nbr]]$Vrai[2]){result[[1]]=difsim[[nbr]]
    best=difsim[[nbr]]$Vrai[2]}
    if(bests<difsimS[[nbr]]$Vrai[2]){result[[2]]=difsimS[[nbr]]
    bests=difsimS[[nbr]]$Vrai[2]}
    if(bests<difsimDV[[nbr]]$Vrai[2]){result[[3]]=difsimDV[[nbr]]
    bestdv=difsimDV[[nbr]]$Vrai[2]}
    
  }
  
  return(result)
}


sim_EM_alp_VRScst = function(C,K,I,p,alpha,beta,gamma,N,alphablocked,nbrcoeur){
  phi=funcphi(beta,I,K)
  pi=funcpi(gamma,I)
  
  S=alphablocked[1]
  V=alphablocked[3]
  cst=alphablocked[4]
  Y=NULL
  X=NULL
  nbr= nbrcoeur
  
  
  registerDoParallel(nbr)
  
  
  
  
  
  sumalpha=rep(0,4)
  sumbeta=rep(0,2)
  sumgamma=rep(0,1)
  summalpha=rep(0,4)
  summbeta=rep(0,2)
  summgamma=rep(0,1)
  sumalphas=rep(0,4)
  sumbetas=rep(0,2)
  sumgammas=rep(0,1)
  summalphas=rep(0,4)
  summbetas=rep(0,2)
  summgammas=rep(0,1)
  
  sumalphacst=rep(0,4)
  sumbetacst=rep(0,2)
  sumgammacst=rep(0,1)
  summalphacst=rep(0,4)
  summbetacst=rep(0,2)
  summgammacst=rep(0,1)
  sumalphav=rep(0,4)
  sumbetav=rep(0,2)
  sumgammav=rep(0,1)
  summalphav=rep(0,4)
  summbetav=rep(0,2)
  summgammav=rep(0,1)
  para_cst=NULL
  para_V=NULL
  para_S=NULL
  para=NULL
  A2=funcA_2(alpha,I,K,C)
  for(j in 1:30){
    
    
    
    
    
    
    
    tempo= bij_inter(1,1,rep(K,C),K)
    
    
    X= array(rep(0,C*(N+1)),dim=c(C,N+1))
    Y= array(rep(0,C*(N+1)),dim=c(C,N))
    #creation simulation avec A2 phi pi
    for(i in 1:C){
      f=runif(1)
      if(f<pi[1]){  X[i,1]=1
      } else if(f <pi[1]+ pi[2]){  X[i,1]=2
      } else if(f <pi[1]+ pi[2]+ pi[3]){  X[i,1]=3
      } else if(f <pi[1]+ pi[2]+pi[3]+pi[4]){  X[i,1]=4
      } else {  X[i,1]=5}  
    }
    
    for(i in 1:N){
      for(c in 1:C){
        l= runif(1)
        
        if(l<phi[X[c,i],1]){ Y[c,i]=1
        } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]){ Y[c,i]=2
        } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+ phi[X[c,i],3]){ Y[c,i]=3
        } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+phi[X[c,i],3]+phi[X[c,i],4]){ Y[c,i]=4
        } else { Y[c,i]=5}  
      }
      for(c in 1:C){
        f=field_to_vect_voisin_inter(c,Y[,i],K)
        R=  runif(1)
        
        
        if(R<A2[X[c,i],1,f]){ X[c,i+1]=1
        } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]){ X[c,i+1]=2
        } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+ A2[X[c,i],3,f]){ X[c,i+1]=3
        } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+A2[X[c,i],3,f]+A2[X[c,i],4,f]){ X[c,i+1]=4
        } else{ X[c,i+1]=5}  
        
      }
      
      
    }
    etatA2=Y
    for(c in 1:C){
      etatA2[c,]=field_to_vect_voisin_inter(c,Y,K)}
    
    bestvrai= rep(-10000,4)
    # names(bestvrai)=c('VraiS','VraiV','Vraicst','Vrai')
    difsim=NULL
    
    difEM_cst=foreach(i=1:nbr) %dopar% EM_alp_cst(Y,etatA2,C,K,I,p,cst)
    difEM_V=foreach(i=1:nbr) %dopar% EM_alp_V(Y,etatA2,C,K,I,p,V)
    difEM_S=foreach(i=1:nbr) %dopar% EM_alp_S(Y,etatA2,C,K,I,p,S)
    difEM=foreach(i=1:nbr) %dopar% EM_alp(Y,etatA2,C,K,I,p)
    
    
    for(i in 1:nbr){
      #ini du meulleur jeu param avec vraisemblance
      if(difEM_cst[[i]]$Vrai[2]>bestvrai[3]){bestvrai[3]=difEM_cst[[i]]$Vrai[2]
      para_cst[j]=list(difEM_cst[[i]]$paramestim)}
      if(difEM_V[[i]]$Vrai[2]>bestvrai[2]){bestvrai[2]=difEM_V[[i]]$Vrai[2]
      para_V[j]=list(difEM_V[[i]]$paramestim)}
      if(difEM_S[[i]]$Vrai[2]>bestvrai[1]){bestvrai[1]=difEM_S[[i]]$Vrai[2]
      para_S[j]=list(difEM_S[[i]]$paramestim)}
      if(difEM[[i]]$Vrai[2]>bestvrai[4]){bestvrai[4]=difEM[[i]]$Vrai[2]
      para[j]=list(difEM[[i]]$paramestim)}
      
    }
    
    #BIAIS et EQM
    sumalpha= sumalpha + para[[j]]$estimalpha
    sumbeta= sumbeta + para[[j]]$estimbeta
    sumgamma =sumgamma + para[[j]]$estimgamma
    summalpha=summalpha+ (alpha - para[[j]]$estimalpha)^2
    summbeta= summbeta +(beta - para[[j]]$estimbeta)^2
    summgamma =summgamma +(gamma- para[j]$estimgamma)^2
    
    sumalphas= sumalphas + para_S[[j]]$estimalpha
    sumbetas= sumbetas +  para_S[[j]]$estimbeta
    sumgammas =sumgammas +  para_S[[j]]$estimgamma
    summalphas=summalphas+ (alpha -  para_S[[j]]$estimalpha)^2
    summbetas= summbetas +(beta - para_S[[j]]$estimbeta)^2
    summgammas =summgammas +(gamma-  para_S[[j]]$estimgamma)^2
    
    
    
    sumalphav= sumalphav + para_V[[j]]$estimalpha
    sumbetav= sumbetav + para_V[[j]]$estimbeta
    sumgammav =sumgammav + para_V[[j]]$estimgamma
    summalphav=summalphav+ (alpha - para_V[[j]]$estimalpha)^2
    summbetav= summbetav +(beta - para_V[[j]]$estimbeta)^2
    summgammav =summgammav +(gamma- para_V[[j]]$estimgamma)^2
    
    sumalphacst= sumalphacst + para_cst[[j]]$estimalpha
    sumbetacst= sumbetacst +  para_cst[[j]]$estimbeta
    sumgammacst =sumgammacst +  para_cst[[j]]$estimgamma
    summalphacst=summalphacst+ (alpha -  para_cst[[j]]$estimalpha)^2
    summbetacst= summbetacst +(beta -  para_cst[[j]]$estimbeta)^2
    summgammacst =summgammacst +(gamma-  para_cst[[j]]$estimgamma)^2
    
    
  }
  
  aBiais=alpha -sumalphas/(j)
  bBiais=beta-sumbetas/(j)
  gBiais=gamma-sumgammas/(j)
  aBiai=alpha -sumalpha/(j)
  bBiai=beta-sumbeta/(j)
  gBiai=gamma-sumgamma/(j)
  aEQM=summalpha/j
  bEQM=summbeta/j
  gEQM=summgamma/j
  aEQMs=summalphas/j
  bEQMs=summbetas/j
  gEQMs=summgammas/j
  aBiaiv=alpha -sumalphav/(j)
  bBiaiv=beta-sumbetav/(j)
  gBiaiv=gamma-sumgammav/(j)
  aBiaicst=alpha -sumalphacst/(j)
  bBiaicst=beta-sumbetacst/(j)
  gBiaicst=gamma-sumgammacst/(j)
  aEQMcst=summalphacst/j
  bEQMcst=summbetacst/j
  gEQMcst=summgammacst/j
  aEQMv=summalphav/j
  bEQMv=summbetav/j
  gEQMv=summgammav/j
  
  result=NULL
  result=list(para_cst,para_V,para_S,para,aBiais,bBiais,gBiais,aBiai,bBiai,gBiai,aEQM,bEQM,gEQM,aEQMs,bEQMs,gEQMs,aBiaiv,bBiaiv,gBiaiv,aBiaicst,bBiaicst,gBiaicst,aEQMcst,bEQMcst,gEQMcst,aEQMv,bEQMv,gEQMv)
  names(result)=c('para_cst','para_V','para_S','para','aBiais','bBiais','gBiais','aBiai','bBiai','gBiai','aEQM','bEQM','gEQM','aEQMs','bEQMs','gEQMs','aBiaiv','bBiaiv','gBiaiv','aBiaicst','bBiaicst','gBiaicst','aEQMcst','bEQMcst','gEQMcst','aEQMv','bEQMv','gEQMv')
  
  save(result,file="selection_model_alpha(7,2,1,-4.5)_N100_C10_fixehySVCST_final.Rdata")
  
  
  return(result)
}

sim_EM_moy_VRScst = function(C,K,I,p,alpha,beta,gamma,N,alphablocked,nbrcoeur){
  phi=funcphi(beta,I,K)
  pi=funcpi(gamma,I)
  
  S=alphablocked[1]
  V=alphablocked[3]
  cst=alphablocked[4]
  Y=NULL
  X=NULL
  nbr= nbrcoeur
  
  
  registerDoParallel(nbr)
  
  
  
  
  
  sumalpha=rep(0,4)
  sumbeta=rep(0,2)
  sumgamma=rep(0,1)
  summalpha=rep(0,4)
  summbeta=rep(0,2)
  summgamma=rep(0,1)
  sumalphas=rep(0,4)
  sumbetas=rep(0,2)
  sumgammas=rep(0,1)
  summalphas=rep(0,4)
  summbetas=rep(0,2)
  summgammas=rep(0,1)
  
  sumalphacst=rep(0,4)
  sumbetacst=rep(0,2)
  sumgammacst=rep(0,1)
  summalphacst=rep(0,4)
  summbetacst=rep(0,2)
  summgammacst=rep(0,1)
  sumalphav=rep(0,4)
  sumbetav=rep(0,2)
  sumgammav=rep(0,1)
  summalphav=rep(0,4)
  summbetav=rep(0,2)
  summgammav=rep(0,1)
  para_cst=NULL
  para_V=NULL
  para_S=NULL
  para=NULL
  A2=funcA_moy(alpha,I,K,C)
  for(j in 1:1){
    
    
    
    
    
    
    
    tempo= bij_moy(1,1,rep(K,C),K,C)
    
    
    X= array(rep(0,C*(N+1)),dim=c(C,N+1))
    Y= array(rep(0,C*(N+1)),dim=c(C,N))
    #creation simulation avec A2 phi pi
    for(i in 1:C){
      f=runif(1)
      if(f<pi[1]){  X[i,1]=1
      } else if(f <pi[1]+ pi[2]){  X[i,1]=2
      } else if(f <pi[1]+ pi[2]+ pi[3]){  X[i,1]=3
      } else if(f <pi[1]+ pi[2]+pi[3]+pi[4]){  X[i,1]=4
      } else {  X[i,1]=5}  
    }
    
    for(i in 1:N){
      for(c in 1:C){
        l= runif(1)
        
        if(l<phi[X[c,i],1]){ Y[c,i]=1
        } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]){ Y[c,i]=2
        } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+ phi[X[c,i],3]){ Y[c,i]=3
        } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+phi[X[c,i],3]+phi[X[c,i],4]){ Y[c,i]=4
        } else { Y[c,i]=5}  
      }
      for(c in 1:C){
        f=field_to_vect_moy(c,Y[,i],K)
        R=  runif(1)
        
        
        if(R<A2[X[c,i],1,f]){ X[c,i+1]=1
        } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]){ X[c,i+1]=2
        } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+ A2[X[c,i],3,f]){ X[c,i+1]=3
        } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+A2[X[c,i],3,f]+A2[X[c,i],4,f]){ X[c,i+1]=4
        } else{ X[c,i+1]=5}  
        
      }
      
      
    }
    etatA2=Y
    for(c in 1:C){
      etatA2[c,]=field_to_vect_moy(c,Y,K)}
    
    bestvrai= rep(-10000,4)
    # names(bestvrai)=c('VraiS','VraiV','Vraicst','Vrai')
    difsim=NULL
    
    difEM_cst=foreach(i=1:nbr) %dopar% EM_moy_cst(Y,etatA2,C,K,I,p,cst)
    difEM_V=foreach(i=1:nbr) %dopar% EM_moy_V(Y,etatA2,C,K,I,p,V)
    difEM_S=foreach(i=1:nbr) %dopar% EM_moy_S(Y,etatA2,C,K,I,p,S)
    difEM=foreach(i=1:nbr) %dopar% EM_moy(Y,etatA2,C,K,I,p)
    
    
    for(i in 1:nbr){
      #ini du meulleur jeu param avec vraisemblance
      if(difEM_cst[[i]]$Vrai[2]>bestvrai[3]){bestvrai[3]=difEM_cst[[i]]$Vrai[2]
      para_cst[j]=list(difEM_cst[[i]]$paramestim)}
      if(difEM_V[[i]]$Vrai[2]>bestvrai[2]){bestvrai[2]=difEM_V[[i]]$Vrai[2]
      para_V[j]=list(difEM_V[[i]]$paramestim)}
      if(difEM_S[[i]]$Vrai[2]>bestvrai[1]){bestvrai[1]=difEM_S[[i]]$Vrai[2]
      para_S[j]=list(difEM_S[[i]]$paramestim)}
      if(difEM[[i]]$Vrai[2]>bestvrai[4]){bestvrai[4]=difEM[[i]]$Vrai[2]
      para[j]=list(difEM[[i]]$paramestim)}
      
    }
    
    #BIAIS et EQM
    sumalpha= sumalpha + para[[j]]$estimalpha
    sumbeta= sumbeta + para[[j]]$estimbeta
    sumgamma =sumgamma + para[[j]]$estimgamma
    summalpha=summalpha+ (alpha - para[[j]]$estimalpha)^2
    summbeta= summbeta +(beta - para[[j]]$estimbeta)^2
    summgamma =summgamma +(gamma- para[[j]]$estimgamma)^2
    
    sumalphas= sumalphas + para_S[[j]]$estimalpha
    sumbetas= sumbetas +  para_S[[j]]$estimbeta
    sumgammas =sumgammas +  para_S[[j]]$estimgamma
    summalphas=summalphas+ (alpha -  para_S[[j]]$estimalpha)^2
    summbetas= summbetas +(beta - para_S[[j]]$estimbeta)^2
    summgammas =summgammas +(gamma-  para_S[[j]]$estimgamma)^2
    
    
    
    sumalphav= sumalphav + para_V[[j]]$estimalpha
    sumbetav= sumbetav + para_V[[j]]$estimbeta
    sumgammav =sumgammav + para_V[[j]]$estimgamma
    summalphav=summalphav+ (alpha - para_V[[j]]$estimalpha)^2
    summbetav= summbetav +(beta - para_V[[j]]$estimbeta)^2
    summgammav =summgammav +(gamma- para_V[[j]]$estimgamma)^2
    
    sumalphacst= sumalphacst + para_cst[[j]]$estimalpha
    sumbetacst= sumbetacst +  para_cst[[j]]$estimbeta
    sumgammacst =sumgammacst +  para_cst[[j]]$estimgamma
    summalphacst=summalphacst+ (alpha -  para_cst[[j]]$estimalpha)^2
    summbetacst= summbetacst +(beta -  para_cst[[j]]$estimbeta)^2
    summgammacst =summgammacst +(gamma-  para_cst[[j]]$estimgamma)^2
    
    
  }
  
  aBiais=alpha -sumalphas/(j)
  bBiais=beta-sumbetas/(j)
  gBiais=gamma-sumgammas/(j)
  aBiai=alpha -sumalpha/(j)
  bBiai=beta-sumbeta/(j)
  gBiai=gamma-sumgamma/(j)
  aEQM=summalpha/j
  bEQM=summbeta/j
  gEQM=summgamma/j
  aEQMs=summalphas/j
  bEQMs=summbetas/j
  gEQMs=summgammas/j
  aBiaiv=alpha -sumalphav/(j)
  bBiaiv=beta-sumbetav/(j)
  gBiaiv=gamma-sumgammav/(j)
  aBiaicst=alpha -sumalphacst/(j)
  bBiaicst=beta-sumbetacst/(j)
  gBiaicst=gamma-sumgammacst/(j)
  aEQMcst=summalphacst/j
  bEQMcst=summbetacst/j
  gEQMcst=summgammacst/j
  aEQMv=summalphav/j
  bEQMv=summbetav/j
  gEQMv=summgammav/j
  
  result=NULL
  result=list(para_cst,para_V,para_S,para,aBiais,bBiais,gBiais,aBiai,bBiai,gBiai,aEQM,bEQM,gEQM,aEQMs,bEQMs,gEQMs,aBiaiv,bBiaiv,gBiaiv,aBiaicst,bBiaicst,gBiaicst,aEQMcst,bEQMcst,gEQMcst,aEQMv,bEQMv,gEQMv)
  names(result)=c('para_cst','para_V','para_S','para','aBiais','bBiais','gBiais','aBiai','bBiai','gBiai','aEQM','bEQM','gEQM','aEQMs','bEQMs','gEQMs','aBiaiv','bBiaiv','gBiaiv','aBiaicst','bBiaicst','gBiaicst','aEQMcst','bEQMcst','gEQMcst','aEQMv','bEQMv','gEQMv')
  
  save(result,file="selection_model_alpha(7,2,1,-4.5)_N100_C10_fixehySVCST_final.Rdata")
  
  
  return(result)
}

predictFB= function(gammaFB,phi){
  D=dim(gammaFB)
  I=D[1]
  C=D[2]
  N=D[3]-1
  pred= rep(0,3)
  for(c in 1:C){
    sum=0
    for(i in 1:I){
      sum=sum+gammaFB[i,c,N+1]*phi[i,]
    }
    pred[c]=which.max(sum)
    
    
  }
  return(pred)
}


EtapeE_moy_LH= function(Y,pi,A2,phi,C){
  D=dim(Y)
  if(is.null(D)){
    if(C!=1){N=1} else{N= length(Y)}

    
  } else{
    C=D[1]
    N= D[2]}

  # beug ici converge pas car des fois betaFB = nan ...
  #QQ= array(QQ,dim=c(I,K,bij_inter(1,1,rep(I,C),I))) dans le VEM inter ici on aura du (I,c,n)
  alphaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  betaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  gammaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  xiFB=array(rep(0,C*(N+1)*I*I),dim=c(I,I,C,N+1))
  for(c in 1:C){
    alphaFB[,c,1]=pi
    betaFB[,c,N+1]=rep(1,I)
    for(j in 1:N){
      k=N-j+1
      etatA2=field_to_vect_moy(c,Y,K)
      tempalpha= colSums(alphaFB[,c,j]*phi[,Y[c,j]]* A2[,,etatA2[j]])
      alphaFB[,c,j+1] = tempalpha/sum(tempalpha)
      #erreur ici ! le beta ne se calcule pas comme ca betaFB[,c,k+1]*phi[,Y[c,k]] 
      #or beta[1,,]*phi[1,].. pas toujours ca peu etre beta[2,,]*phi[1,]
      #
      tempbeta= rowSums(t(betaFB[,c,k+1]*t(phi[,Y[c,k]]* A2[,,etatA2[k]])))
      betaFB[,c,k]= tempbeta/sum(tempbeta)
      
    }
  }
  result=list(alphaFB,betaFB)
  return(result)
}

EtapeE_A2_BP_LH= function(Y,pi,A2,phi,C){
  D=dim(Y)
  if(is.null(D)){
    if(C!=1){N=1} else{N= length(Y)}
    
  } else{
    C=D[1]
    N= D[2]}
  
  # beug ici converge pas car des fois betaFB = nan ...
  #QQ= array(QQ,dim=c(I,K,bij_inter(1,1,rep(I,C),I))) dans le VEM inter ici on aura du (I,c,n)
  alphaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  betaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  gammaFB=array(rep(0,C*(N+1)*I),dim=c(I,C,N+1))
  xiFB=array(rep(0,C*(N+1)*I*I),dim=c(I,I,C,N+1))
  for(c in 1:C){
    alphaFB[,c,1]=pi
    betaFB[,c,N+1]=rep(1,I)
    for(j in 1:N){
      k=N-j+1
      etatA2=field_to_vect_voisin_inter(c,Y,K)
      tempalpha= colSums(alphaFB[,c,j]*phi[,Y[c,j]]* A2[,,etatA2[j]])
      alphaFB[,c,j+1] = tempalpha/sum(tempalpha)
      #erreur ici ! le beta ne se calcule pas comme ca betaFB[,c,k+1]*phi[,Y[c,k]] 
      #or beta[1,,]*phi[1,].. pas toujours ca peu etre beta[2,,]*phi[1,]
      #
      tempbeta= rowSums(t(betaFB[,c,k+1]*t(phi[,Y[c,k]]* A2[,,etatA2[k]])))
      betaFB[,c,k]= tempbeta/sum(tempbeta)
      
    }
  }
  
  
  
  
  
  result=list(alphaFB,betaFB)
  return(result)
  
  
}

BP_LLH= function(alphaFB,betaFB,phi,Y,C,N,I){
  
  
  moy=0
  best = -Inf
  for(c in 1:C){
    for(j in 1: N){
      
      cvoulu=which((1:C)!=c)
      temp= alphaFB[,c,j]*betaFB[,c,j]
      firstpart=  log(sum(temp))

      prod2=1
      prod1=0
      for(n in 1:N){
        prodtemp = 1
        prodtemp1=0 
        for(u in cvoulu){
          sumit=0
          for(i in 1:I ){
            sumit=sumit + phi[i,Y[u,n]]*alphaFB[i,u,n]/sum(alphaFB[,u,n])
            }
          tempp=1
          prodtemp1= prodtemp1+log(sumit)
        }
        prod1=prod1+prodtemp1
      }
      if(best<firstpart+prod1){best=firstpart+prod1}
      moy= moy + firstpart+prod1
    }
  }
  
  result=list(best,moy/(C*N))
  return(result)
}



# hh = histograme du log des data sans les zeros
KULhisto= function(hh){

  
  kull= function(param){
    nbr=sum(hh$counts)
    hh$breaks=round(100*hh$breaks)/100
    dist=hh$breaks[2]- hh$breaks[1]
    j=1
    k=param[1]
    keepind<-list()
    
  densi=rep(0,length(param)+1)
    for(i in param){
      if(i == param[1]){
        ind= which(hh$breaks<i)
        keepind[[1]]=ind
        freq=sum(hh$counts[ind])/nbr
        interval=length(ind)*dist
        densi[1]=freq/interval
      }else{ 
        ind =which((k)<hh$breaks & hh$breaks<i)
        keepind[[j]]=ind
        freq=sum(hh$counts[ind])/nbr
        interval=length(ind)*dist
        densi[j]=freq/interval
        
        k=i
        }
      j=j+1
    }

  ind =which((i)<hh$breaks)
  freq=sum(hh$counts[ind[1:(length(ind)-1)]])/nbr
  interval=(length(ind)-1)*dist
  keepind[[j]]=ind[1:(length(ind)-1)]
  densi[j]=freq/interval
  
  #print("densité")
   # print(densi)
  #  print("parametres")
  # print(param)
  
  sumn=0
  j=1
  for(k in 1 : length(hh$density)){
    if(length(which(keepind[[j]]==k))==0){
      j=j+1
    }
    if(hh$density[k] ==0 || densi[j]==0 ){}else{
    sumn= sumn + dist*hh$density[k]*log(hh$density[k]/densi[j])}
  }
  
  
  
    result= (sumn)
    return(result)}
  return(kull)
}


densit= function(hh,param){
  nbr=sum(hh$counts)
  hh$breaks=floor(100*hh$breaks)/100
  dist=hh$breaks[2]- hh$breaks[1]
  j=1
  k=param[1]
  keepind<-list()
  
  densi=rep(0,length(param)+1)
  for(i in param){
    if(i == param[1]){
      ind= which(hh$breaks<i)
      keepind[[1]]=ind
      freq=sum(hh$counts[ind])/nbr
      interval=length(ind)*dist
      densi[1]=freq/interval
    }else{ 
      ind =which((k)<hh$breaks & hh$breaks<i)
      keepind[[j]]=ind
      freq=sum(hh$counts[ind])/nbr
      interval=length(ind)*dist
      densi[j]=freq/interval
      
      k=i
    }
    j=j+1
  }
  
  ind =which((i)<hh$breaks)
  freq=sum(hh$counts[ind[1:(length(ind)-1)]])/nbr
  interval=(length(ind)-1)*dist
  keepind[[j]]=ind[1:(length(ind)-1)]
  densi[j]=freq/interval
  
  test=0
  for(k in 1:length(keepind)){
   test= test+ length(keepind[[k]])*dist*densi[k]
  }
  result=list(densi,keepind)
  names(result)=c('Densité','interval')
  return(result)
}


classes = function(hh,param){
  hh$breaks=floor(100*hh$breaks)/100
  dist=hh$breaks[2]- hh$breaks[1]
  dist2=dist+0.00001
  fkull=KULhisto(hh)
  
  constraint=array(rep(0,5*4),dim=c(5,4))
  diag(constraint)=rep(-1,4)
  constraint[1,2]=1
  constraint[2,3]=1
  constraint[3,4]=1
  constraint[4,4]=-1
  constraint[5,1]=1
  cires=c(dist2,dist2,dist2,-max(hh$breaks),min(hh$breaks))
  resuA=constrOptim(param, fkull,NULL, ui=constraint, ci=cires )
  
  
  result=resuA$par

  
  
  return(result)
}

classes2 = function(hh,param){
  hh$breaks=floor(100*hh$breaks)/100
  dist=hh$breaks[2]- hh$breaks[1]
  dist2=dist+0.00001
  fkull=KULhisto(hh)
  
  constraint=array(rep(0,4*3),dim=c(4,3))
  diag(constraint)=rep(-1,3)
  constraint[1,2]=1
  constraint[2,3]=1

  constraint[4,3]=1
  cires=c(dist2,dist2,-max(hh$breaks),min(hh$breaks))
  resuA=constrOptim(param, fkull,NULL, ui=constraint, ci=cires )
  
  
  result=resuA$par
  
  
  
  return(result)
}


calcul_de_s_c_g = function(alpha,beta,gamma,cultd,dist,C,N,I,K){
  
  nbrcult=max(cultd)
  freqall=array(rep(0,nbrcult*5*5*5*5),dim=c(nbrcult,5,5,5,5))
  freqX=array(rep(0,I*nbrcult),dim=c(nbrcult,I))
 
  #ici on fait les proba pour comparer avec celle de pluntz
  for(z in 1 : nbrcult){
    pi=funcpi(gamma,I)
      if(nbrcult==1){
        phi=ZIphi(beta,I,K)
        A2=funcA_moy(alpha,I,K,C)
  
      }else{  phi=ZIphi(beta[,z],I,K)
      A2=funcA_moy(alpha[,z],I,K,C) }
      
      for(q in 1:500){
       
        
        X= array(rep(0,C*(N+1)),dim=c(C,N+1))
        Y= array(rep(0,C*(N+1)),dim=c(C,N))
        #creation simulation avec A2
        t=runif(1)
        X[,1]=length(which(cumsum(pi)<t))+1
        
        
        
        
        for(i in 1:N){
          for(c in 1:C){
            l= runif(1)
            
            if(l<phi[X[c,i],1]){ Y[c,i]=1
            } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]){ Y[c,i]=2
            } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+ phi[X[c,i],3]){ Y[c,i]=3
            } else if(l <phi[X[c,i],1]+ phi[X[c,i],2]+phi[X[c,i],3]+phi[X[c,i],4]){ Y[c,i]=4
            } else { Y[c,i]=5}  
          }
          for(c in 1:C){
            f=field_to_vect_dist(c,Y[,i],K,dist)
            R=  runif(1)
            
            
            if(R<A2[X[c,i],1,f]){ X[c,i+1]=1
            } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]){ X[c,i+1]=2
            } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+ A2[X[c,i],3,f]){ X[c,i+1]=3
            } else if(R <A2[X[c,i],1,f]+ A2[X[c,i],2,f]+A2[X[c,i],3,f]+A2[X[c,i],4,f]){ X[c,i+1]=4
            } else{ X[c,i+1]=5}  
            
          }
          
          
        }
        for(i in 1:N){
  
          
          for(c in 1:C){
            freqX[z,X[c,i]]=freqX[z,X[c,i]]+1
            voi=field_to_vect_dist(c,Y[,i],K,dist)%%K
            if(voi==0){voi=5}
            freqall[z,Y[c,i],voi,X[c,i],X[c,i+1]]=freqall[z,Y[c,i],voi,X[c,i],X[c,i+1]]+1
            
          
              
            }
            
            
          }
         
        }}
  
  
 
  
  
  survi=rep(0,nbrcult)
  colon=rep(0,nbrcult)
  colone=rep(0,nbrcult)
  germ=rep(0,nbrcult)

  survibis=rep(0,nbrcult)
  colone2=rep(0,nbrcult)

  

    for(z in 1 : nbrcult){
      if(nbrcult==1){
        colext=funcA_moy(alpha,I,K,C)[1,1,1]
      }else{ colext=funcA_moy(alpha[,z],I,K,C)[1,1,1]}
      
      survi[z]=1-(sum(freqall[z,1,1,-1,1])/sum(freqall[z,1,1,-1,]))/colext
      colon[z]= 1- sum(freqall[z,1,,1,1])/sum(freqall[z,1,,1,])
      colone[z]= 1- (sum(freqall[z,1,-1,1,1])/sum(freqall[z,1,-1,1,]))/colext
      germ[z]=1- sum(freqall[z,1,,-1,])/sum(freqall[z,,,-1,])
     
     
      
    }
  result=list(survi,colon,colone,germ,freqX)
  names(result)<-c('s','cext','cvoi','g','freqX')
  return(result)
  }
  
  
  
