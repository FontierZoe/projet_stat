
rm(list=ls())



library ( MASS )
############### CLASSIFICATION D'IMAGES DE CHATS ET DE CHIENS ########################
############################# Gautier Appert #########################################
# Groupe : Alexandre Junius  - Xavier Van Poelvoorde - Marie Guegain

# 1. DECOUVERTE DE LA BASE DE DONNEES ET REDUCTION DE LA DIMENSION ----
#QUESTION 1: Importer les bases de donnees ----
setwd("/Users/zoefontier/Desktop/ensae 2A/stats")
load("Xtest.RData")
load("Xtrain.RData")
load("Ytrain.RData")
load("Ytest.RData")

# Conversion de quelques images en matrices 200x200
Image1 <- matrix(Xtest[1,], nrow=200, ncol=200)
Image2 <- matrix(Xtest[2,], nrow=200, ncol=200)
Image3 <- matrix(Xtest[3,], nrow=200, ncol=200)
Image30 <- matrix(Xtest[13,], nrow=200, ncol=200)

Image1bis <-matrix(Xtrain[1,], nrow=200, ncol=200)
Image2bis <-matrix(Xtrain[2,], nrow=200, ncol=200)
Image3bis <-matrix(Xtrain[3,], nrow=200, ncol=200)

#Fonction de rotation pour les images
rotate = function(x) t(apply(x, 2, rev))

# Affichage de quelques images de Xtest
image(rotate(Image1), col=grey(seq(0,1, length = 256)))
image(rotate(Image2), col=grey(seq(0,1, length = 256)))
image(rotate(Image3), col=grey(seq(0,1, length = 256)))
image(rotate(Image30), col=grey(seq(0,1, length = 256)))

rasterImage(Image1,0,0,dim(Image1)[2],dim(Image1)[1])

# Affichage de quelques images de Xtrain
image(rotate(Image1bis), col=grey(seq(0,1, length = 256)))
image(rotate(Image2bis), col=grey(seq(0,1, length = 256)))
image(rotate(Image3bis), col=grey(seq(0,1, length = 256)))


#QUESTION 2 : Réduction de la dimension des donnees par une ACP ----

#Concatenation des bases Xtrain et Xtest
Xtot = rbind(Xtrain,Xtest)

#Centrage des vecteurs colonnes 
X = scale(Xtot,scale=F)

#ACP en utilisant une decomposition en valeur singulière (SVD) de X. 
svdres = svd(X, nu=5, nv=15)
#On ne retient que les 30 premieres composantes principales

str(svdres)
#Part de variance expliquee par 30 composantes principales ? 
L=svdres$d^2
pct=cumsum(L)/sum(L)
plot(pct,type='l')

pct[15]
sum(pct[1:30])
#On obtient 80,4%

#Composantes principales C
C=X%*%svdres$v

#Decoupage de la base en Ctrain et Ctest
Ctrain = C[1:315,]
Ctest = C[316:363,]




indice1=which(Ytest==1)
indice0=which(Ytest==0)

plot(Ctest[indice1,1],Ctest[indice1,2],col="blue",xlim=c(-200,200),ylim=c(-100,100))
points(Ctest[indice0,1],Ctest[indice0,2],col="red")

     

############################################################
## Function that computes Maximum Likelihood estimators---##
############################################################



computeML=function(C,Y){
  
  indice1=which(Y==1)
  indice0=which(Y==0)
  
  ## piHat
  piHat= length(indice1)/length(Y)
  
  ## muHat
   mu1Hat= colMeans(C[indice1,])
   mu0Hat= colMeans(C[indice0,])
  
  ## SigmaHat
   Cmu1=apply(C[indice1,],1,FUN=function(x){x-mu1Hat})
   Sigma1Hat = Cmu1%*%t(Cmu1)/length(indice1)
   
   Cmu0=apply(C[indice0,],1,FUN=function(x){x-mu0Hat})
   Sigma0Hat = Cmu0%*%t(Cmu0)/length(indice0)
   
   return(list(piHat=piHat,mu1Hat=mu1Hat, mu0Hat= mu0Hat, Sigma1Hat=Sigma1Hat, Sigma0Hat=Sigma0Hat))
  
}



computeLogRatio=function(cNew,pi,mu1,mu0,Sigma1,Sigma0){
  return(
  -(1/2)*log(det(Sigma1))-(1/2)*(cNew-mu1)%*%
    solve(Sigma1)%*%(cNew-mu1)
  +log(pi)+(1/2)*log(det(Sigma0))
  +(1/2)*(cNew-mu0)%*%solve(Sigma0)%*%(cNew-mu0)-log(1-pi)
  )
}



computePred=function(C,pi,mu1,mu0,Sigma1,Sigma0){
  return(ifelse(apply(C,1,FUN=function(x){ computeLogRatio(x,pi,mu1,mu0,Sigma1,Sigma0) })>0,1,0))
  }


ML=computeML(Ctrain,Ytrain)
QDA.Pred = computePred(Ctest,ML$piHat,ML$mu1Hat,ML$mu0Hat,ML$Sigma1Hat,ML$Sigma0Hat)

sum(abs(Ytest-QDA.Pred))/length(Ytest)


indice1=which(QDA.Pred==1)
indice0=which(QDA.Pred==0)

plot(Ctest[indice1,1],Ctest[indice1,2],col="blue",xlim=c(-200,200),ylim=c(-200,100),pch=19)
points(Ctest[indice0,1],Ctest[indice0,2],col="red",pch=19)

indexDiff=which(QDA.Pred != Ytest)


# Conversion de quelques images en matrices 200x200
Image1 <- matrix(Xtest[indexDiff[6],], nrow=200, ncol=200)

# Affichage de quelques images de Xtest
image(rotate(Image1), col=grey(seq(0,1, length = 256)))

rasterImage(Image1,0,0,dim(Image1)[2],dim(Image1)[1])

# Affichage de quelques images de Xtrain
image(rotate(Image1bis), col=grey(seq(0,1, length = 256)))


indexDiff=which(QDA.Pred != Ytest)




length(ML$mu1Hat)
dim(Ctrain)
a=Ctrain-ML$mu1Hat
dim(a)

b=sweep(Ctrain,2,ML$mu1Hat)

a==b
a=matrix(1:12,nrow=3,ncol=4)
a-1:4

t(apply(a,1,FUN=function(x){x-1:4}))

computeML()





lda.fit = lda ( Ctrain,Ytrain )
lda.fit


plot ( lda.fit )


Ctest = data.frame(Ctest)
head(Ctest)
dim(Ctest)

lda.pred = predict (lda.fit ,Ctest)

lda.pred$class


qda.fit = qda (Ctrain,Ytrain )
str(qda.fit)
ML$mu1Hat
qda.fit$ldet
qda.fit$means
qda.pred=predict (qda.fit , as.data.frame(Ctest))$class



aa=as.data.frame(Ctest)

head(aa)
class(aa)

sum(abs(Ytest-as.numeric(as.character(qda.pred))))/length(Ytest)
qda.pred


logit=glm(Ytrain~Ctrain, family = "binomial")
summary(logit)



# 2. MODELE LOGISTIQUE ET PENALITE RIDGE ----

#Question 6: Coder une fonction qui calcule le maximum de vraisemblance penalise ----

#Importation de la library Matrix pour sa fonction norm qui permet de calculer la norme2 d'un vecteur :
library(Matrix)
#Exemple : norm(matrix(1:3, ncol = 1), "F")

#Fonction pi qui permettra d'ecrire le vecteur PI
pi=function(x){
  exp(x)/(1+exp(x))
}

#Fonction pi(1-pi) qui permettra d'écrire W 
pi2=function(x){
  exp(x)/((1+exp(x))^2)
}

#Algorithme de Newton-Raphson: 
newton <- function(C, y, lambda, eps){
  beta_old=rep(1,30)
  beta_new=rep(0,30)                                                   #Valeur initiale de Beta
  while (norm(as.matrix(beta_new-beta_old), "F")>eps) {         
    beta_old<-beta_new
    c_beta = C%*%beta_old                                     
    PI=sapply(c_beta,FUN=pi)                                         #vecteur PI
    W=diag(sapply(c_beta,FUN=pi2))                                   #Création de la matrice diagonale W
    z=solve(W)%*%(y-PI)+C%*%beta_old                                   #Création de z
    beta_new= solve(t(C)%*%W%*%C+2*lambda*diag(1,30))%*%t(C)%*%W%*%z   #Calcul de la nouvelle valeur de Beta
  }
  beta_new
}  



# 3. VALIDATION CROISEE ET PREDICTION DES LABELS SUR LA BASE TEST ----

#QUESTION 8: Validation croisee ----

#Creation des echantillons de la base de donnees
library(caret)
echantillon <-createFolds(1:nrow(Ctrain),k=15)

#Calcul de prediction des yi
prediction = function(base,beta) {                                  
  c_beta <- base%*%beta
  PI <- sapply(c_beta, FUN=pi)
  y <- rep(0,nrow(C))
  y <- ifelse(PI>=0.5,1,0)                   #selon la regle de decision                    
  y                                                     
}

#Support de lambda 
support <- c(0,100,1000,2000,10000,15000,45000)

#Fonction pour calculer la bonne valeur de lambda parmi celles du support
lambda_max <- function(support,echantillon) {
  vector_lambda=rep(0,length(support))                
  vector_beta=matrix(data=0,length(echantillon),ncol(Ctrain))
  vector_MSE=rep(0,length(echantillon))
  for (i in 1:length(support)) {
    for (j in 1:length(echantillon)) {
      Ctrain_VC = Ctrain[unlist(echantillon[-j]),]                        
      Ytrain_VC=Ytrain[unlist(echantillon[-j])]                            
      Ctrain_init = Ctrain[unlist(echantillon[j]),]                       
      Ytrain_init=Ytrain[unlist(echantillon[j])] 
      vector_prediction = rep(0,nrow(Ctrain_init))                           
      vector_beta[j,]=newton(Ctrain_VC,Ytrain_VC,support[i],0.0001)
      vector_prediction=prediction(Ctrain_init,vector_beta[j,])
      vector_MSE[j]=(1/length(Ytrain_init))*norm(as.matrix(vector_prediction) - as.matrix(Ytrain_init),"1")
    }
    vector_lambda[i]=(1/length(echantillon))*sum(vector_MSE)                
  }
  k=which.min(vector_lambda)                                 
  support[k]
}


#Valeur de lambda apres validation croisee:
lambda <- lambda_max(support,echantillon)                                 
#On obtient lambda = 45 000
#Il faut donc choisir le modele avec penalite

#QUESTION 9: Prediction sur toute la base ----

#Estimation du modele logistique sur toute la base de donnees Ctrain
Beta=newton(Ctrain,Ytrain,lambda, 0.0001)

#Prediction des labels sur la base Ctest
Prediction_C_test=prediction(Ctest,Beta)

#Fonction calculant le pourcentage d'erreurs
error <- function(x,y) { 
  E=0
  for (i in 1:length(x)) {
    if (x[i]!=y[i]) {
      E=E+1
    }
  }
  E/length(x)
}


#Calcul de l'erreur
error (Ytest,lda.pred$class)
#On trouve 12,5% d'erreurs (rejet si > 50%)

