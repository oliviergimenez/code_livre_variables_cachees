#FONCTIONS UTILISEES DANS LE SCRIPT MAIN

library(lavaan)
library(mnormt)
library(car)
library(lavaanPlot)
library(parallel)

normalize <- function(v){

#OBJECTIF : appliquer une transformation puissance sur une variable numérique pour obtenir une distribution proche de la distribution gaussienne

#ARGUMENTS :
#	- v : vecteur numérique correspondant à la variable à transformer

#SORTIE
#	- vecteur numérique correspondant à la variable transformée, centrée et réduite
	
	pT <- powerTransform(v~1) #requiert la bibliothèque "car"
	w <- v^pT$lambda
	return((w-mean(w))/sd(w))
}

extModGen <- function(latVarNam,vecManVarNam,record=FALSE,root="",id=""){

#OBJECTIF : générer un fragment de script en syntaxe compatible avec la bibliothèque lavaan correspondant à fragment de modèle de mesure comprenant un bouquet de variables observées et leur variable latente associée; le modèle inclut des contraintes de standardisation des variables observées, et une contrainte de positivité sur l'effet de la variable latente.

#ARGUMENTS :
#	- latVarNam : chaîne de caractères; nom de la variable latente sous-jacente au bouquet
#	- vecManVar : vecteur de chaînes de caractères; vecManVar[i] est le nom de la i-ème variable observée du bouquet
#	- record : booléen; TRUE commande l'enregistrement dans un fichier texte du fragment de script généré
#	- root : chaîne de caractères ; adresse du répertoire où est enregistré le fichier texte généré quand record = TRUE
#	- id : chaîne de caractères ; nom du fichier texte généré quand record =TRUE

#SORTIE :
#	- CFAMod : chaîne de caractères; fragment de script lavaan 

	nManVar <- length(vecManVarNam)
	latManMod <-  paste(latVarNam,
		' =~ ',
		paste(sapply(vecManVarNam[-nManVar],function(manVarNam) paste('alp_',manVarNam,'_',latVarNam,'*',manVarNam,' +',sep='')),collapse=''),
		' alp_',vecManVarNam[nManVar],'_',latVarNam,'*',vecManVarNam[nManVar],
		sep=''
	)
	resManVar <- paste(sapply(vecManVarNam,function(manVarNam) paste(manVarNam,' ~~ dlt_',manVarNam,'2*',manVarNam,sep='')),collapse=' \n')
	resLatVar <- paste(latVarNam,' ~~ dlt_',latVarNam,'2*',latVarNam,sep='')
	if(length(vecManVarNam)>1) posConstr <- paste(c(
		sapply(vecManVarNam,function(manVarNam) paste('alp_',manVarNam,'_',latVarNam,' > 10^(-3)',sep='')),
		sapply(vecManVarNam,function(manVarNam) paste('alp_',manVarNam,'_',latVarNam,' < 1-10^(-3)',sep='')),
		sapply(vecManVarNam,function(manVarNam) paste('dlt_',manVarNam,'2 == 1 - alp_',manVarNam,'_',latVarNam,'^2',sep=''))
	),collapse='\n')
	if(length(vecManVarNam)==1) posConstr <- paste(
		sapply(vecManVarNam[1],function(manVarNam) paste('alp_',manVarNam,'_',latVarNam,' == 1-10^(-3)',sep='')),
		sapply(vecManVarNam[1],function(manVarNam) paste('dlt_',manVarNam,'2 == 1 - alp_',manVarNam,'_',latVarNam,'^2',sep=''))
	,sep='\n')
	CFAMod <- paste(
		'##Bloc de ',latVarNam,' \n',
		latManMod,'\n',
		'		###Variance residuelle des variables mesurees \n',
		resManVar,'\n',
		'		###Variance de la variable latente \n',
		resLatVar,'\n',
		'		###Contraintes sur les parametres de mesure \n',
		posConstr,'\n',
		sep=''
	)
	if(record) write.table(CFAMod,file=paste(root,"CFAMod_",id,".txt",sep=""),row.names=FALSE,col.names=FALSE)
	return(CFAMod)
}

regBlocGen <- function(lisReg,remAllReg){

#OBJECTIF : générer un fragment de script en syntaxe compatible avec la bibliothèque lavaan correspondant à un modèle relationnel 

#ARGUMENTS :
#	- lisReg : liste de chaînes de caractères; lisReg[[i]] contient l'ensemble des relations (fèches rouges das la charte graphique du chapitre) pointant vers une même variable latente; ces relations sont spécifiées dans une syntaxe similaire à celle des objets formula de la bibliothèque stats de R, utilisée dans les approches statistiques classique telles que lm ou glm : "variableLatenteReponse ~ predicteurLatent1 + predicteurLatent2 + predicteurLatent3" 
#	- remAllReg : variable booléenne; remAllReg = TRUE force toutes les relations spécifiées par lisReg à 0

#SORTIE :
#	- chaîne de caractères; fragment de script lavaan comprenant toutes les relations spécifiées par lisReg et forçant les coefficients beta associés à 0 si remAllReg=TRUE

	if(length(lisReg)>0){
		lisFinDeb <- lapply(lisReg,function(reg){
			vecMemb <- unlist(strsplit(reg,split="~",fixed=TRUE))
			fin <- unlist(strsplit(vecMemb[1],split=" ",fixed=TRUE))
			vecDeb <- unlist(strsplit(vecMemb[2],split="+",fixed=TRUE))
			vecDeb <- sapply(vecDeb,function(deb){
				v <- unlist(strsplit(deb,split=" ",fixed=TRUE))
				indRem <- which(v=="")
				if(length(indRem)>0) v <- v[-indRem]
				return(v)
			})
			return(list(fin=fin,deb=vecDeb))
		})
		vecScrReg <- sapply(lisFinDeb,function(finDeb){
			vecBet <- paste("bet",finDeb$deb,finDeb$fin,sep="_")
			vecMonom <- paste(vecBet,"*",finDeb$deb,sep="")
			membDroite <- paste(vecMonom,collapse=" + ")
			return(paste(finDeb$fin,"~",membDroite,sep=" "))
		})
		if(remAllReg){
			vecScrReg <- c(vecScrReg,
				unlist(sapply(lisFinDeb,function(finDeb){
					return(paste("bet_",finDeb$deb,"_",finDeb$fin," == 0",sep=""))
				}))
			)
		}
		return(paste("##Regressions",paste(vecScrReg,collapse="\n"),sep="\n"))
	}
	return("")
}

covBlocGen <- function(lisReg,vecNamLV,force=NULL){

#OBJECTIF : générer un fragment de script en syntaxe compatible avec la bibliothèque lavaan correspondant aux corrélations libres (paramètres gamma) entre les variables latentes exogènes du modèle

#ARGUMENTS :
#	- lisReg : liste de chaînes de caractères; lisReg[[i]] contient l'ensemble des relations (fèches rouges das la charte graphique du chapitre) pointant vers une même variable latente; est utilisée ici pour identifier les variables latentes exogènes;
#	- vecNamLV : vecteur de chaînes de caractères; vecNamLV[i] est la i-ème variable latente du MES à construire. Remarque : les variables citées dans vecNamLV ne sont pas nécessairement toutes citées dans lisReg; Remarque : les ordinations de lisReg et vecNamLV sont indépendantes;
#	- force : matrice symétrique de nombres réels; force[i,j] contient une valeur que l'on souhaite fixer comme contrainte sur la corrélation entre les variables rownames(force)[i] et colnames(force)[j]; la diagonale de cette matrice n'est pas lue.

#SORTIE :
#	- chaîne de caractères; fragment de script lavaan comprenant la spécification de toutes les corrélations libres entre variables latentes exogènes dans le modèle, ainsi que d'éventuelles contraintes numériques sur certaines d'entre elles.


	vecLab <- vecNamLV
	DAG <- matrix(0,length(vecLab),length(vecLab))
	rownames(DAG) <- vecLab
	colnames(DAG) <- vecLab
	if(length(lisReg)>0){
		lisFinDeb <- lapply(lisReg,function(reg){
			vecMemb <- unlist(strsplit(reg,split="~",fixed=TRUE))
			fin <- unlist(strsplit(vecMemb[1],split=" ",fixed=TRUE))
			vecDeb <- unlist(strsplit(vecMemb[2],split="+",fixed=TRUE))
			vecDeb <- sapply(vecDeb,function(deb){
				v <- unlist(strsplit(deb,split=" ",fixed=TRUE))
				indRem <- which(v=="")
				if(length(indRem)>0) v <- v[-indRem]
				return(v)
			})
			return(list(fin=fin,deb=vecDeb))
		})
		# vecLab <- unique(unlist(lisFinDeb))
		if(is.na(sum(match(unique(unlist(lisFinDeb)),vecNamLV)))) return(NA)
		for(i in 1:length(lisReg)){
			fin <- lisFinDeb[[i]]$fin
			vDeb <- lisFinDeb[[i]]$deb
			DAG[fin,vDeb]<-1
		}
		nbPath <- DAG
		for(i in 1:(length(vecLab)-1)) nbPath <- nbPath%*%DAG
		if(sum(nbPath)>0) return(NA)
	} 
	vecIndExo <- which(apply(DAG,1,sum)==0)
	vecLabExo <- vecLab[vecIndExo]
	if(length(vecIndExo)==1) return("##Pas de correlations")
	return(paste(
		c(
			"##Correlations",
			sapply(1:(length(vecIndExo)-1),function(i){
				paste(
					sapply((i+1):length(vecIndExo),function(j){
						paste(vecLabExo[i]," ~~ gam_",vecLabExo[i],"_",vecLabExo[j],"*",vecLabExo[j],sep="")
					}),
					collapse="\n"
				)
			}),
			sapply(1:(length(vecIndExo)-1),function(i){
				paste(
					sapply((i+1):length(vecIndExo),function(j){
						if(!is.null(force)){
							k <- match(vecLabExo[i],rownames(force))
							l <- match(vecLabExo[j],rownames(force))
							if(!is.na(k+l)){
								return(	paste("gam_",vecLabExo[i],"_",vecLabExo[j]," == ",force[k,l],sep=""))
							}
						}
						return(paste("gam_",vecLabExo[i],"_",vecLabExo[j],"^2 < ","dlt_",vecLabExo[i],"2*dlt_",vecLabExo[j],"2",sep=""))
					}),
					collapse="\n"
				)
			})
		),
		collapse="\n"
	))
}

sigGen <- function(lisReg,vecNamLV,remAllReg){

#OBJECTIF : générer l'expression développée des variances des variables latentes du modèle en fonction des variances résiduelles (paramètres delta), des relations entre variables latentes (paramètres beta) et des corrélations libres entre variables latentes exogènes (paramètres gamma), par application récursive de l'équation (5) du chapitre.

#ARGUMENTS :
#	- lisReg : liste de chaînes de caractères; lisReg[[i]] contient l'ensemble des relations (fèches rouges das la charte graphique du chapitre) pointant vers une même variable latente; est utilisée ici pour construire le graphe dirigé acyclique associé au modèle relationnel;
#	- vecNamLV : vecteur de chaînes de caractères; vecNamLV[i] est la i-ème variable latente du MES à construire. Remarque : les variables citées dans vecNamLV ne sont pas nécessairement toutes citées dans lisReg; Remarque : les ordinations de lisReg et vecNamLV sont indépendantes;
#	- remAllReg : variable booléenne; remAllReg = TRUE force toutes les relations spécifiées par lisReg à 0

#SORTIE :
#	- vecteur de chaîne de caractères; le i-eme élément du vecteur contient un fragment de script lavaan donnant l'expression de la variance de la variable latente vecNamLV[i] en fonction des coefficients beta, gamma et delta du modèle. Il y a donc une correspondance terme à terme entre vecNamLV et la sortie.

	if(length(lisReg)>0){
		lisFinDeb <- lapply(lisReg,function(reg){
			vecMemb <- unlist(strsplit(reg,split="~",fixed=TRUE))
			fin <- unlist(strsplit(vecMemb[1],split=" ",fixed=TRUE))
			vecDeb <- unlist(strsplit(vecMemb[2],split="+",fixed=TRUE))
			vecDeb <- sapply(vecDeb,function(deb){
				v <- unlist(strsplit(deb,split=" ",fixed=TRUE))
				indRem <- which(v=="")
				if(length(indRem)>0) v <- v[-indRem]
				return(v)
			})
			return(list(fin=fin,deb=vecDeb))
		})
		# vecLab <- unique(unlist(lisFinDeb))
		if(is.na(sum(match(unique(unlist(lisFinDeb)),vecNamLV)))) return(NA)
		vecLab <- vecNamLV
		
		# Construction du graphe des beta
		DAG <- matrix(0,length(vecLab),length(vecLab))
		rownames(DAG) <- vecLab
		colnames(DAG) <- vecLab
		for(i in 1:length(lisReg)){
			fin <- lisFinDeb[[i]]$fin
			vDeb <- lisFinDeb[[i]]$deb
			DAG[fin,vDeb]<-1
		}
		
		# Vérification du caractère acyclique du graphe
		nbPath <- DAG
		for(i in 1:(length(vecLab)-1)) nbPath <- nbPath%*%DAG
		if(sum(nbPath)>0) return(NA) 
		
		# Calcul de la plus longue chaine d'ancêtre pour chaque variable latente
		vLMax <- rep(0,length(vecLab))
		vIni <- rep(0,length(vecLab))
		vIni[which(apply(DAG,1,sum)==0)] <- 1
		vActu <- vIni
		for(i in 1:(length(vecLab)-1)){
			vActu <- (DAG%*%vActu>0)*1
			vLMax <- vLMax+vActu		
		}

		# Calcul récursif des covariances symboliques
		
		matCovStr <- matrix(0,length(vecLab),length(vecLab))
		rownames(matCovStr) <- vecLab
		colnames(matCovStr) <- vecLab
		vecInd0 <- which(vLMax==0)
		for(i in vecInd0){
			for(j in vecInd0){
				if(i==j){
					matCovStr[i,j] <- paste("dlt_",rownames(DAG)[i],"2",sep="")
				}
				if(i!=j) matCovStr[i,j] <- paste("gam",
					rownames(DAG)[min(c(i,j))],
					rownames(DAG)[max(c(i,j))],
					sep="_"
				)	
			}
		}
		for(l in 1:(max(vLMax))){
			vecIndPrev <- which((vLMax<=(l-1))*(vLMax>0)==1)
			vecInd <- which(vLMax==l)
			for(i in vecInd){
				for(j in vecInd0){
					vTerm <- NULL 
					vAntI <- which(DAG[i,]>0)
					for(k in vAntI){
						if(!remAllReg){
							vTerm <- c(
								vTerm,
								paste(
									"bet_",rownames(DAG)[k],"_",rownames(DAG)[i],"*",
									"(",matCovStr[k,j],")",
									sep=""
								)
							)
						}
					}
					term <- paste(vTerm,collapse=" + ")
					matCovStr[i,j] <- term
					matCovStr[j,i] <- term
				}
				for(j in c(vecIndPrev,vecInd)){
					vTerm <- NULL 
					vAntI <- which(DAG[i,]>0)
					vAntJ <- which(DAG[j,]>0)
					for(k in vAntI){
						for(l in vAntJ){
							if(!remAllReg){
								vTerm <- c(
									vTerm,
									paste(
										"bet_",rownames(DAG)[k],"_",rownames(DAG)[i],"*",
										"bet_",rownames(DAG)[l],"_",rownames(DAG)[j],"*",
										"(",matCovStr[k,l],")",
										sep=""
									)
								)
							}
						}
					}
					if(i==j) term <- paste(c(vTerm,paste("dlt_",rownames(DAG)[i],"2",sep="")),collapse =" + ")
					else term <- paste(vTerm,collapse=" + ")
					matCovStr[i,j] <- term
					if(vLMax[j]<=(l-1)) matCovStr[j,i] <- term			
				}
			}	
		}
		return(diag(matCovStr))
	}else{
		return(paste("dlt_",vecNamLV,"2",sep=""))
	}
}

SEMGen <- function(lisBloc,lisReg,covBloc=TRUE,stdConstr=TRUE,remAllReg=FALSE, record=FALSE,force=NULL,root="./",id="test"){

#OBJECTIF : générer le script intégral d'un SEM en syntaxe compatible avec la bibliothèque lavaan de R à partir des spécifications utilisateur d'un modèle de mesure, un modèle relationnel et d'une série de contraintes sur les variances et covariances entre variables latentes

#ARGUMENTS:
#	- lisBloc : liste de vecteurs de chaines de caractères; names(lisBloc)[i] est le nom d'une variable latente, et lisBloc[[i]] est un vecteur de noms de variables observées rattachées en bouquet à cette variable latente
#	- lisReg : liste de chaînes de caractères; lisReg[[i]] contient l'ensemble des relations (fèches rouges dans la charte graphique du chapitre) pointant vers une même variable latente; ces relations sont spécifiées dans une syntaxe similaire à celle des objets formula de la bibliothèque stats de R, utilisée dans les approches statistiques classique telles que lm ou glm : "variableLatenteReponse ~ predicteurLatent1 + predicteurLatent2 + predicteurLatent3" 
#	- covBloc : variable booléenne; si covBloc=TRUE (défaut) alors les corrélations entre variables exogènes sont toutes autorisées et associées à un terme de type gamma, sinon elles sont toutes forcées à 0 (hypothèse d'indépendance entre variables exogènes) ;
#	- stdConstr : variable booléenne; si stdConstr=TRUE (défaut) alors les variables latentes du modèles sont toutes standardisées et le faisceau de contraintes sur les paramètres est généré en conséquence, sinon les variances sont laissées libres;
#	- remAllReg : variable booléenne; remAllReg = TRUE force toutes les relations spécifiées par lisReg (paramètres beta) à 0;
#	- record : variable booléenne; record =TRUE commande l'enregistrement dans un fichier texte du script généré
#	- force : matrice symétrique définie positive de nombres réels; force[i,j] contient une valeur que l'on souhaite fixer comme contrainte sur la corrélation (paramètre gamma) entre les variables exogène rownames(force)[i] et colnames(force)[j]; la diagonale de cette matrice n'est pas lue.
#	- root : chaîne de caractères ; adresse du répertoire où est enregistré le fichier texte généré quand record = TRUE
#	- id : chaîne de caractères ; nom du fichier texte généré quand record =TRUE;

#SORTIE
#	- chaîne de caractère contenant un script de MES compatible avec la syntaxe lavaan

	vecNamLV <- NULL
	if(length(lisReg)>0){
		lisFinDeb <- lapply(lisReg,function(reg){
			vecMemb <- unlist(strsplit(reg,split="~",fixed=TRUE))
			fin <- unlist(strsplit(vecMemb[1],split=" ",fixed=TRUE))
			vecDeb <- unlist(strsplit(vecMemb[2],split="+",fixed=TRUE))
			vecDeb <- sapply(vecDeb,function(deb){
				v <- unlist(strsplit(deb,split=" ",fixed=TRUE))
				indRem <- which(v=="")
				if(length(indRem)>0) v <- v[-indRem]
				return(v)
			})
			return(list(fin=fin,deb=vecDeb))
		})
		vecNamLV <- c(vecNamLV,unique(unlist(lisFinDeb)))
	}
	vecNamLV <- unique(c(vecNamLV,names(lisBloc)))
	vecModMes <- sapply((1:length(lisBloc)),function(iBloc){
		return(extModGen(names(lisBloc)[iBloc],lisBloc[[iBloc]]))
	})
	modIntern <- paste(
		"#Modele relationnel",
		regBlocGen(lisReg,remAllReg),
		sep="\n\n"
	)
	if(covBloc){
		modIntern <- paste(modIntern,
			covBlocGen(lisReg,vecNamLV,force),
			sep="\n\n"
		)
	}
	if(stdConstr){
		if(!is.null(force)){
			# constr[match(rownames(force),vecNamLV)] <- diag(force)
			modIntern <- paste(modIntern,
				paste(c("##Contraintes de forçage sur certaines variances variables latentes exogènes",
					paste(
						"dlt_",
						rownames(force),
						"2 == ",
						diag(force),
						sep=""
					)
				),collapse="\n"),
				sep="\n\n"
			)
			indNotForced <- which(is.na(match(vecNamLV,rownames(force))))
			if(length(indNotForced)>0){
				modIntern <- paste(modIntern,
					paste(c("##Contraintes de standardisation des variables latentes non forcées",
						paste(
							"(",
							sigGen(lisReg,vecNamLV,remAllReg)[indNotForced],
							" - 1)^2 + ((dlt_",
							vecNamLV[indNotForced],
							"2-10^(-3))-sqrt((dlt_",
							vecNamLV[indNotForced],
							"2-10^(-3))^2))^2 < 10^(-6)",
							sep=""
						)
					),collapse="\n"),
					sep="\n\n"
				)
			}
		}else{
			modIntern <- paste(modIntern,
				paste(c("##Contraintes de standardisation des variables latentes",
					paste(
						"(",
						sigGen(lisReg,vecNamLV,remAllReg),
						" - 1)^2 + ((dlt_",
						vecNamLV,
						"2-10^(-3))-sqrt((dlt_",
						vecNamLV,
						"2-10^(-3))^2))^2 < 10^(-6)",
						sep=""
					)
				),collapse="\n"),
				sep="\n\n"
			)
		}
	}
	semScript <- paste(c("#Modele de mesure",vecModMes,modIntern),collapse="\n\n") 
	if(record) cat(semScript,file=paste(root,id,".txt",sep=""))
	return(semScript)
}

rem1Reg <- function(lisReg){

#OBJECTIF : générer tous les modèles relationnels amputés d'une relation (i.e. une flèche rouge) à partir d'un modèle relationnel de référence

#ARGUMENT :
#	- lisReg : liste de chaînes de caractères; lisReg[[i]] contient l'ensemble des relations (fèches rouges dans la charte graphique du chapitre) pointant vers une même variable latente; ces relations sont spécifiées dans une syntaxe similaire à celle des objets formula de la bibliothèque stats de R, utilisée dans les approches statistiques classique telles que lm ou glm : "variableLatenteReponse ~ predicteurLatent1 + predicteurLatent2 + predicteurLatent3" 

#SORTIE
#	- lisLisReg : liste de listes de chaines de caractère; lisLisReg[[i]] contient un objet similaire à lisReg mais où la i-ème flèche (en comptant dans l'ordre des éléments de lisReg et des prédicteurs au sein de ces éléments) a été soustraite

	lisLisReg <- NULL
	if(length(lisReg)>0){
		lisFinDeb <- lapply(lisReg,function(reg){
			vecMemb <- unlist(strsplit(reg,split="~",fixed=TRUE))
			fin <- unlist(strsplit(vecMemb[1],split=" ",fixed=TRUE))
			vecDeb <- unlist(strsplit(vecMemb[2],split="+",fixed=TRUE))
			vecDeb <- sapply(vecDeb,function(deb){
				v <- unlist(strsplit(deb,split=" ",fixed=TRUE))
				indRem <- which(v=="")
				if(length(indRem)>0) v <- v[-indRem]
				return(v)
			})
			return(list(fin=fin,deb=vecDeb))
		})
		for(i in 1:length(lisFinDeb)){
			if(length(lisFinDeb[[i]]$deb)==1){
				if(length(lisFinDeb)>1){
					blocReg2 <- lapply((1:length(lisFinDeb))[-i],function(indReg){
						reg <- lisFinDeb[[indReg]]
						return(paste(reg$fin,paste(reg$deb,collapse=" + "),sep=" ~ "))
					})
					lisLisReg <- c(lisLisReg,list(blocReg2))
				}
			}else{
				for(j in 1:length(lisFinDeb[[i]]$deb)){
					lisFinDeb2 <- lisFinDeb
					lisFinDeb2[[i]]$deb <- lisFinDeb2[[i]]$deb[-j]
					blocReg2 <- lapply(lisFinDeb2,function(reg){
						return(paste(reg$fin,paste(reg$deb,collapse=" + "),sep=" ~ "))
					})
					lisLisReg <- c(lisLisReg,list(blocReg2))
				}
			}
			
		}
	}
	return(lisLisReg)
}
