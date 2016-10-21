#######################################################################################
## Name:         mlevCSrealMFVB.R (Demo for Insight Interview)
## Purpose :     Real time Bayesian semiparametric mixed effects regression analysis
##               via mean field variational Bayes
## Dataset:      New South Wales Perinatal Data Collection 
## Author :      Cathy Y. Y. Lee
## Last updated: 14 SEP 2015
## R Version:    3.2.0 (2015-04-16)                                                             
## Input data files:  CS2010.csv                                                     
## Output data files: ---
## Required R packages: MASS, magic, Matrix, mgcv, maptools, ggmap, RgoogleMaps
########################################################################################

setwd("/Users/cathylee/Documents/UTS/PhD/Programs/R/R.scripts")

# Clear R memory:

rm(list=ls())
set.seed(12345)

# Load required libraries:

libname <- c("MASS","magic","Matrix","mgcv","maptools","ggmap","RgoogleMaps","lattice")
lapply(libname, require, character.only=T)

# Set flags:

userSavedFits <- T;    
doBatchFit <- T;       doRealTimeAddOn <- T
plotMFVB <- T;         printLogML <- F
createPDF <- F;        userwait <- F

# Set iteration size values:

maxIterONL <- 1;       maxIter <- 200
   
# Set colours and plotting size parameters:

ptCol <- "dodgerblue"
lineColMFVB <- "blue"
lineColTruth <- "grey40"
lwdVal <- 2; cexVal <- 1.5

numPlots <- 20
output.path <- "/Users/cathylee/Documents/UTS/PhD/Programs/R/R.outputs/SimResults/Chpt.Six.SimResults/"
MFVB.plots <- vector(numPlots, mode='list')

# Set hyperparameters:

nuVal <- 2;   A.u <- 1e5 
A.R <- 1e5;   sigsq.beta <- 1e5

# Load user-written functions:

source("twoLevelMFVB.r")
source("rankprofile.r")
source("addAlpha.r")
source("bernFunctions.r")

# Read data from a polygon shapefile into a SpatialPolygonsDataFrame object:

NSW.SLA06.map <- readShapePoly("NSWSLA06.shp",IDvar="SLA_CODE06", 
                               proj4string=CRS("+proj=longlat +ellps=clrk66"))

# Retrieves spatial bounding box from spatial data:

#bb <- bbox(NSW.SLA06.map)

# Query the Google server for a static map, defined primarily by its lat/lon range:

#MyMap <- GetMap.bbox(bb[1,], bb[2,], destfile = "Australia.png",GRAYSCALE=FALSE)
#saveRDS(MyMap,file=paste(getwd(),"/MyMap.rds",sep=""))

# Read in caesarean section data from the New South Wales Perinatal Data Collection:

CSdata <- read.table("/Users/cathylee/Documents/UTS/PhD/Programs/R/R.data/Data4Thesis/CS2010.csv",
                     header=T,sep=",")
head(CSdata)

# Extract latitudes and longitudes of hospital locations from Google Maps API:

#hoslatlon <- geocode(as.character(unique(CSdata$hospital)))
hoslatlon <- read.table("hoslatlon.csv",header=T,sep=",")

# Standardise maternal age:

age <- CSdata$age
meanAge <- mean(age); sdAge <- sd(age)

# Set up response and predictor variables, and dimension parameters for the full data analysis:

yMax <- CSdata$cs_birth               # Response variable: caesarean section (Yes or No)
x1Max <- CSdata$pricare               # Predictor: private obstetric care (Yes or No)
x2Max <- (age - meanAge)/sdAge        # Predictor: maternal age (continuous, standardised)
x3Max <- CSdata$sla06res              # Statistical Local Area of residence 
idnumMax <- CSdata$hosID              # Hospital ID
 
mMax <- length(unique(idnumMax))      # Number of groups/clusters (i.e. hospitals)
nVecMax <- as.vector(table(idnumMax)) # Sample size within hospitals
numObsMax <- sum(nVecMax)             # Total number of observations

currStt <- 1                          # Indicator for the random effects 
reBlockIndsMax <- vector("list", length=mMax)
for (i in 1:mMax)
{         
   currEnd <- currStt + nVecMax[i] - 1
   reBlockIndsMax[i] <- list(currStt:currEnd)
   currStt <- currEnd + 1  
}

# Obtain spline basis function design matrices for the full data analysis::

numIntKnots <- 25
a <- 1.01*min(x2Max)-0.01*max(x2Max); b <- 1.01*max(x2Max)-0.01*min(x2Max)
intKnots <-  quantile(unique(x2Max),seq(0,1,length=numIntKnots+2)[-c(1,numIntKnots+2)])
ZGmax <- ZOSull(x2Max,intKnots=intKnots,range.x=c(a,b))
ncZG <- ncol(ZGmax); L <- length(ncZG)

# Obtain constant matrices and dimension variables for the full data analysis::

XRmax <- cbind(rep(1,numObsMax),x1Max); XGmax <- x2Max
CGmax <- cbind(XRmax,XGmax,ZGmax);      Xmax <- cbind(XRmax,XGmax)
ZRmax <- matrix(0,numObsMax,2*mMax)
for (i in 1:mMax)
{
   indsCurr <- (1:numObsMax)[idnumMax==i]
   ZRmax[indsCurr,c((2*i-1),2*i)] <- cbind(rep(1,nVecMax[i]),x1Max[indsCurr])
}

ncXR <- ncol(XRmax); ncCG <- ncol(CGmax)
ncX <- ncol(Xmax)

# Panel plot of caesarean data: each panel respresent a hospital:

dev.new(height=6.5,width=10.5)
source("plotPanel.r")
wait()

# Set up response and predictor variables, and dimension parameters for the
# batch and real-time data analysis:

propWarm <- 0.95
nVecWarm <- round(propWarm*nVecMax)
nVecReal <- nVecMax - nVecWarm
numObsWarm <- sum(nVecWarm)

indsWarm <- NULL; indsWarmList <- vector("list",mMax)
indsRealList <- vector("list",mMax)
for (i in 1:mMax)
{
   indsWarm <- c(indsWarm,reBlockIndsMax[[i]][1:nVecWarm[i]])
   indsWarmList[[i]] <- reBlockIndsMax[[i]][1:nVecWarm[i]]
   indsRealList[[i]] <- setdiff(reBlockIndsMax[[i]],indsWarmList[[i]])
}
indsReal <- unlist(indsRealList)

y <- yMax[indsWarm]
x1 <- x1Max[indsWarm]
x2 <- x2Max[indsWarm]
x3 <- x3Max[indsWarm]

reBlockInds <- vector("list",length=mMax)
currStt <- 1
for (i in 1:mMax)
{
   currEnd <- currStt + length(indsWarmList[[i]]) - 1
   reBlockInds[i] <- list(currStt:currEnd)
   currStt <- currEnd + 1
}

# Feed XRwarm, XGwarm and Xwarm etc into the batch MFVB Algorithm to obtain
# the starting values for q-density moments:

if(doBatchFit|doRealTimeAddOn)
{
   XR <- cbind(rep(1,numObsWarm),x1)
   XG <- x2
   X <- cbind(XR,XG)
   ZG <- ZGmax[indsWarm,]
   CG <- cbind(XR,XG,ZG)
   ZR <- ZRmax[indsWarm,]
   sla <- x3
}

# Perform batch-based run:

if (doBatchFit)
{
   m <- mMax
   nVec <- nVecWarm
   numObs <- numObsWarm
   
   if (!userSavedFits)
   {
      MFVBbatchFit <- twoLevelMFVB(y,XG,XR,ZG,ZR=NULL,reBlockInds,ncZG,
                                   responseType="Bernoulli",
                                   doStreamlined=TRUE,
                                   maxIter=200,useMatForDg=TRUE)
      saveRDS(MFVBbatchFit,file=paste(getwd(),"/MFVBbatchFit.rds",sep=""))
   }
   
   if (userSavedFits)
      MFVBbatchFit <- readRDS(paste(getwd(),"/MFVBbatchFit.rds",sep=""))
}

if (doRealTimeAddOn)
{
   # Set up grid for plotting later on:
    
   ng <- 101
   x1g <- rep(mean(x1Max),ng)
   x2g <- seq(a,b,length=ng)
   ZGg <- ZOSull(x2g,intKnots=intKnots,range.x=c(a,b))
   Xg <- cbind(rep(1,ng),x1g,x2g) 
   CGg <- cbind(Xg,ZGg)
       
   # Initialize parameters with batch values:

   nVec <- nVecWarm                                      # Sample size within each hospital
   mu.q.betauG <- MFVBbatchFit$mu.q.betauG               # Regression and spline coefficients 
   Sigma.q.betauG <- MFVBbatchFit$Sigma.q.betauG         # Covariance matrix for regression and spline coefficients
   Gm <- MFVBbatchFit$Gm                                 # G matrix
   Hm <- MFVBbatchFit$Hm                                 # H matrix 
   mu.q.uRm <- MFVBbatchFit$mu.q.uRm                     # Random intercept and slope coefficients
   Sigma.q.uRm <- MFVBbatchFit$Sigma.q.uRm               # Covariance matrix for random intercepts and slopes
   xiVec <- MFVBbatchFit$xiVec                           # Variational parameters
   M.q.inv.SigmaR <- MFVBbatchFit$M.q.inv.SigmaR  
   mu.q.recip.aR <- MFVBbatchFit$mu.q.recip.aR           # Auxilliary parameters
   mu.q.recip.sigsq.u <- MFVBbatchFit$mu.q.recip.sigsq.u # Smooting parameter for spline
   rankMat <- NULL
   
   indsCurrList <- vector("list",length=mMax)
   for (i in 1:mMax)
      indsCurrList[[i]] <- indsWarmList[[i]]

   yCurrList <- vector("list",length=mMax)
   x2CurrList <- vector("list",length=mMax)
   XRCurrList <- vector("list",length=mMax)
   CGcurrList <- vector("list",length=mMax)
   xiCurrList <- vector("list",length=mMax)
   slaCurrList <- vector("list",length=mMax)
   for (i in 1:mMax)
   {
       yCurrList[[i]] <- y[reBlockInds[[i]]]
       x2CurrList[[i]] <- x2[reBlockInds[[i]]]
       XRCurrList[[i]] <- XR[reBlockInds[[i]],]
       CGcurrList[[i]] <- CG[reBlockInds[[i]],]
       xiCurrList[[i]] <- xiVec[reBlockInds[[i]]]
       slaCurrList[[i]] <- sla[reBlockInds[[i]]]
   }

   CGTy <- crossprod(CG,y-0.5)
   ZRTy <- vector("list", length=mMax)
   for (i in 1:mMax)
      ZRTy[[i]] <- crossprod(XR[reBlockInds[[i]],],y[reBlockInds[[i]]]-0.5)

   #dev.new(height=6.5,width=10.5)
   blankPlot()
   addWhitePoly()
   text(0.5,0.7,"WELCOME TO REAL TIME \n VARIATIONAL BAYES ANALYSIS!",cex=1.5,col="navy")
   text(0.5,0.5,paste("Batch Phase (Initalisation): No. hospitals = ",mMax),cex=1.5,col="navy")
   #text(0.5,0.5,paste("Batch: No. hospitals = ",mMax,", No. obs = ",numObsWarm),cex=1.5,col="navy")
   wait()
   
   ii <- 0
   while (ii < mMax) # Loop over each hospital:
   {
      # Read in new data:

      ii <- ii + 1; n <- 0
      
      #if (ii > 1) dev.new(height=6.5,width=11\0.5)
      blankPlot()
      par(mfrow=c(1,1))
      addWhitePoly()
      text(0.5,0.5,paste("Commencing update for \n Hospital ID = ",ii),cex=2,col="DeepPink")
      Sys.sleep(1)
      
      while (n < sum(nVecReal[ii])) # Loop over new records within each hospital:
      {
         if (ii <= 3)
         {
            addWhitePoly()
            text(0.5,0.7,"Percentage of iterations completed is",cex=1,col="blue")
            text(0.5,0.4,round(n*100/nVecReal[ii]),cex=7,col="blue")
            Sys.sleep(0.1)
         }
         
         # Update dimension parameters and data-related vectors and matrices:
         
         n <- n + 1
         nVec[ii] <- nVec[ii] + 1
         numObs <- sum(nVec)
      
         indsCurrList[[ii]] <- sort(c(indsCurrList[[ii]],indsRealList[[ii]][n]))
         reBlockInds <- vector("list",length=mMax)
         currStt <- 1
         for (i in 1:mMax)
         {
            currEnd <- currStt + length(indsCurrList[[i]]) - 1
            reBlockInds[i] <- list(currStt:currEnd)
            currStt <- currEnd + 1
         }

         yNew <- yMax[indsRealList[[ii]][n]]
         x2New <- x2Max[indsRealList[[ii]][n]]         
         Xnew <- Xmax[indsRealList[[ii]][n],]
         XRnew <- XRmax[indsRealList[[ii]][n],]
         ZGnew <- ZGmax[indsRealList[[ii]][n],]
         CGnew <- c(Xnew,ZGnew)
         slaNew <- x3Max[indsRealList[[ii]][n]]         

         xiNewSqd <- diag(matrix(CGnew,1,ncCG)%*%(Sigma.q.betauG
                          + tcrossprod(mu.q.betauG))%*%matrix(CGnew,ncCG,1))
         EsqMatCurr <- (-Sigma.q.betauG)%*%Gm%*%Hm + tcrossprod(mu.q.betauG,mu.q.uRm)
         xiNewSqd <- (xiNewSqd + 2*diag(matrix(CGnew,1,ncCG)%*%EsqMatCurr%*%matrix(XRnew,2,1)))
         EsqMatCurr <- Sigma.q.uRm + tcrossprod(mu.q.uRm)
         xiNewSqd <- (xiNewSqd + diag(matrix(XRnew,1,2)%*%EsqMatCurr%*%matrix(XRnew,2,1))) 
         xiNew <- sqrt(xiNewSqd)
         
         yCurrList[[ii]] <- c(yCurrList[[ii]],yNew)
         x2CurrList[[ii]] <- c(x2CurrList[[ii]],x2New)
         XRCurrList[[ii]] <- rbind(XRCurrList[[ii]],XRnew)
         CGcurrList[[ii]] <- rbind(CGcurrList[[ii]],CGnew)
         xiCurrList[[ii]] <- c(xiCurrList[[ii]],xiNew)
         slaCurrList[[ii]] <- c(as.character(slaCurrList[[ii]]),as.character(slaNew))
    
         y <- NULL;  x2 <- NULL;    XR <- NULL
         CG <- NULL; xiVec <- NULL; sla <- NULL
         for (i in 1:m)
         {
            y <- c(y,yCurrList[[i]])
            x2 <- c(x2,x2CurrList[[i]])
            XR <- rbind(XR,XRCurrList[[i]])
            CG <- rbind(CG,CGcurrList[[i]])
            xiVec <- c(xiVec,xiCurrList[[i]])
            sla <- c(sla,as.character(slaCurrList[[i]]))
         }
         
         CGTy <- CGTy + as.matrix(CGnew*(yNew-0.5))
         ZRTy[[ii]] <- ZRTy[[ii]] + as.matrix(XRnew*(yNew-0.5))
         
         # Create lists of matrices required for MFVB:
            
         G <- vector("list",length=mMax) 
         H <- vector("list",length=mMax)

         # Mean field variational Bayes algorithmic updates:
         
         itnum <- 0; converged <- FALSE         
         while (!converged) 
         {
            itnum <- itnum + 1
         
            wtVec <- 2*lambda(xiVec)
            phiVec <- phi(xiVec)
            
            ridgeVec <- rep((1/sigsq.beta),ncX)
            for (ell in 1:L)
               ridgeVec <- c(ridgeVec,rep(mu.q.recip.sigsq.u[ell],ncZG[ell]))
                      
            # Update q*(beta,u) parameters:
         
            sVec <- rep(0,ncCG); Smat <- matrix(0,ncCG,ncCG)
            for (i in 1:m)
            {            
               G[[i]] <- (crossprod(CG[reBlockInds[[i]],],
                          wtVec[reBlockInds[[i]]]*XR[reBlockInds[[i]],]))
               H[[i]] <- (solve(crossprod(XR[reBlockInds[[i]],],
                          wtVec[reBlockInds[[i]]]*XR[reBlockInds[[i]],]) + M.q.inv.SigmaR))
               sVec <- sVec + as.vector(G[[i]]%*%H[[i]]%*%ZRTy[[i]])
               Smat <- Smat + G[[i]]%*%H[[i]]%*%t(G[[i]])
            }
            
            Sigma.q.betauG <- (solve(crossprod(CG,wtVec*CG) + diag(ridgeVec) - Smat))
            mu.q.betauG <- as.vector(Sigma.q.betauG%*%(CGTy - sVec))
            
            Sigma.q.uR <- vector("list",length=m)
            mu.q.uR <- vector("list",length=m)
            
            for (i in 1:m)
            {
               Sigma.q.uR[[i]] <- (H[[i]] + H[[i]]%*%crossprod(G[[i]],
                                   (Sigma.q.betauG%*%G[[i]]%*%H[[i]])))
               mu.q.uR[[i]] <- as.vector(H[[i]]%*%(ZRTy[[i]] -
                                         crossprod(G[[i]],mu.q.betauG)))
            }
                                 
            # Update xiVec parameter:
                
            xiSqd <- diag(CG%*%(Sigma.q.betauG + tcrossprod(mu.q.betauG))%*%t(CG))
            for (i in 1:m)
            {
               EsqMatCurr <- (-Sigma.q.betauG)%*%G[[i]]%*%H[[i]]
                   + tcrossprod(mu.q.betauG,mu.q.uR[[i]])
               xiSqd[reBlockInds[[i]]] <- (xiSqd[reBlockInds[[i]]] 
                   + 2*diag(CG[reBlockInds[[i]],]%*%EsqMatCurr%*%t(XR[reBlockInds[[i]],])))
               EsqMatCurr <- Sigma.q.uR[[i]] + tcrossprod(mu.q.uR[[i]])
               xiSqd[reBlockInds[[i]]] <- (xiSqd[reBlockInds[[i]]] 
                   + diag(XR[reBlockInds[[i]],]%*%EsqMatCurr%*%t(XR[reBlockInds[[i]],])))
            }
            xiVec <- sqrt(xiSqd)

            # Reassign G and H matrices:
      
            Gm <- G[[m]]; Hm <- H[[m]]
            mu.q.uRm <- mu.q.uR[[m]]
            Sigam.q.uRm <- Sigma.q.uR[[m]]
            xiCurrList <- vector("list",length=m)
            for (i in 1:m)
               xiCurrList[[i]] <- xiVec[reBlockInds[[i]]]

            # Update q*(a.R) parameters:
            
            B.q.a.R <- nuVal*diag(M.q.inv.SigmaR) + (1/A.R^2)
            mu.q.recip.a.R <- (0.5*(nuVal + ncXR))/B.q.a.R
            
            # Update q*(SigmaR^{-1}) parameters:
            
            B.q.SigmaR <- 2*nuVal*diag(mu.q.recip.a.R)
            indsStt <- ncCG + 1
            for (i in 1:m)
               B.q.SigmaR <- B.q.SigmaR + Sigma.q.uR[[i]] + tcrossprod(mu.q.uR[[i]])
            
            M.q.inv.SigmaR <- (nuVal + m + ncXR - 1)*solve(B.q.SigmaR)
            
            # Update q*(a.u) parameters:
            
            B.q.a.u <- mu.q.recip.sigsq.u + (1/A.u^2)
            mu.q.recip.a.u <- 1/B.q.a.u
            
            # Update q*(sigsq.u) parameters:
               
            indsStt <- ncX+1; A.q.sigsq.u <- NULL; B.q.sigsq.u <- NULL
            trSpline <- NULL
            for (ell in 1:L)
            {
               indsEnd <- indsStt + ncZG[ell] - 1; inds <- indsStt:indsEnd
               A.q.sigsq.u[ell] <- ncZG[ell] + 1
               B.q.sigsq.u[ell] <- (2*mu.q.recip.a.u[ell] + sum(mu.q.betauG[inds]^2)
                                    + sum(diag(Sigma.q.betauG[inds,inds])))
               indsStt <- indsEnd + 1                                              
            }  
            mu.q.recip.sigsq.u <- A.q.sigsq.u/B.q.sigsq.u
            
            if (itnum >= maxIterONL) converged <- TRUE
         }
      }

      # MFVB estimate for the spline:
      
      fhatg <- as.vector(CGg%*%mu.q.betauG)
      sdhatg <- as.vector(sqrt(diag(CGg%*%Sigma.q.betauG%*%t(CGg))))
      lowg <- fhatg - qnorm(0.975)*sdhatg
      uppg <- fhatg + qnorm(0.975)*sdhatg

      if (ii == 1)
      {    
         addWhitePoly()
         text(0.5,0.7,"Graphic outputs for \n hospital benchmarking...",cex=2,col="green3")
         Sys.sleep(1)
      }
            
      layout(matrix(c(1,2,3,4),2,2,byrow=T),height=c(4,4),width=c(4,4))
      
      # Figure 1: Mean field variational Bayes approximate posterior means of
      # the regression function and the dashed curves are pointwise 95% credible
      # sets. The blue open circles are the batch data (jittered) and the red
      # solid circules are the real-time data:
      
      par(mar=c(4,5,0,0),las=1)
      plot(0,0,type="n",bty="l",xlim=range(x2g),ylim=c(-0.05,0.45),
           xlab="Maternal age (years)",ylab="Probability of caesarean",
           cex.lab=cexVal,cex.axis=cexVal,main="",cex.main=cexVal,axes=F)
      axis(2,at=seq(0,0.4,0.1),label=seq(0,0.4,0.1),cex.axis=cexVal)
      axis(1,at=seq(a,b,length=8),label=seq(20,34,length=8),cex.axis=cexVal)
      points(x2[y>0],runif(length(x2[y>0]),0.4,0.45),col=ptCol,cex=0.5)
      points(x2[y==0],runif(length(x2[y==0]),-0.05,0),col=ptCol,cex=0.5)
     
      for (k in 1:nVecReal[[ii]]) # Plot new data:
      {
         if (rev(yCurrList[[ii]])[1:nVecReal[ii]][k]==0)
            points(rev(x2CurrList[[ii]])[1:nVecReal[ii]][k],
                   rev(yCurrList[[ii]])[1:nVecReal[ii]][k],
                   col="red",pch=16,cex=cexVal)
         if (rev(yCurrList[[ii]])[1:nVecReal[ii]][k]==1)
            points(rev(x2CurrList[[ii]])[1:nVecReal[ii]][k],
                   rev(yCurrList[[ii]])[1:nVecReal[ii]][k]-0.55,
                   col="red",pch=16,cex=cexVal)
      }

      # Plot full data estimate:
      
      MFVBfullFit <- readRDS(paste(getwd(),"/MFVBfullFit.rds",sep=""))
      mu.q.betauG <- MFVBfullFit$mu.q.betauG
      Sigma.q.betauG <- MFVBfullFit$Sigma.q.betauG
      fFullg <- as.vector(CGg%*%mu.q.betauG)
      sdFullg <- as.vector(sqrt(diag(CGg%*%Sigma.q.betauG%*%t(CGg))))
      lowFullg <- fFullg - qnorm(0.975)*sdFullg
      uppFullg <- fFullg + qnorm(0.975)*sdFullg
      polygon(c(x2g,rev(x2g)),c(inv.logit(lowFullg),rev(inv.logit(uppFullg))),
              col="grey90",border=FALSE)

      # Plot real-time estimate:
      
      lines(x2g,inv.logit(fhatg),col=lineColMFVB,lwd=lwdVal,lty=1)
      lines(x2g,inv.logit(lowg),type="l",col=lineColMFVB,lty=2,lwd=lwdVal)
      lines(x2g,inv.logit(uppg),type="l",col=lineColMFVB,lty=2,lwd=lwdVal)
      legend(-0.3,0.2,legend=c("New data","Full data estimate","Real-time estimate"),
             col=c("red","grey90",lineColMFVB),pch=c(19,15,NA),lty=c(NA,NA,2),bty="n",
             cex=cexVal,lwd=c(NA,NA,lwdVal))

      # Figure 2: Plot of Statistical Local Areas of residence for women attending
      # the hospital of interest:
               
      slaFreq <- data.frame(table(slaCurrList[[ii]]))
      slaFreq$Freq <- slaFreq[,"Freq"]/nVec[[ii]]
      slaFreq <- slaFreq[order(-slaFreq[,"Freq"]),]
      slaFreq$Cumu <- cumsum(slaFreq[,"Freq"])
      slaFreq$flag <- as.numeric(slaFreq$Cumu<0.95)
      if (nrow(slaFreq) < 5) slaFreq$flag <- 1
      HosSLA <- as.character(slaFreq[slaFreq$flag==1,1])
      NSW.SLA06.map@data$SLA_4DIGIT<- substr(NSW.SLA06.map@data$SLA_5DIGIT,2,5)
      
      # Subset the polygon shapefile to including the relevant SLAs:
      
      HOS.SLA06 <- NSW.SLA06.map[NSW.SLA06.map@data$SLA_4DIGIT%in%HosSLA,]
      HOSSLA06.shp <- SpatialPolygons(HOS.SLA06@polygons,proj4string=HOS.SLA06@proj4string)

      # Overlay the hospital location and relevant SLA polygons on background image of 
      # map of New South Wales, Australia:

      MyMap <- readRDS(paste(getwd(),"/MyMap.rds",sep=""))
      PlotOnStaticMap(MyMap,lat=hoslatlon[ii,2],lon=hoslatlon[ii,1],
                      cex=cexVal,pch=20,mar=c(4,5,0,1),col="red",add=F)
      PlotPolysOnStaticMap(MyMap,HOSSLA06.shp,col=add.alpha("#507415",alpha=.4))

      # Figure 3: Caterpillar plot of hospital-specific deviation from the overall
      # intercept for benchmarking. Assign hospitals into 1) higher-than-avergae;
      # 2) no different from average; or 3) lower-than-average:
      
      randInt <- matrix(0,mMax,5); 
      for (i in 1:mMax)
      {
          randInt[i,1:4] <- c(i,mu.q.uR[[i]][1],
                      mu.q.uR[[i]][1]-1.96*Sigma.q.uR[[1]][1,1],
                      mu.q.uR[[i]][1]+1.96*Sigma.q.uR[[1]][1,1])
          if (i == ii) randInt[i,5] <- 1
      }
      randInt <- cbind(1:mMax,randInt[order(randInt[,2]),])
      colnames(randInt) <- c("index","hosID","int","low","upp","flag")

      par(mar=c(4,5,0,0),las=1)
      plot(randInt[,"index"],randInt[,"int"],pch=19,cex.main=cexVal,
           ylim=range(c(randInt[,"low"],randInt[,"upp"])),
           cex.axis=cexVal,cex.lab=cexVal,cex=1,xaxt="n",axes=F,
           ylab="Hospital deviation from average",
           xlab="Hospitals ordered by increasing deviation")
      abline(h=0,col="dodgerblue",lty=2,lwd=lwdVal)
      axis(1,randInt[,"index"],randInt[,"upp"],label=randInt[,"hosID"],cex.axis=0.8)
      axis(2,at=seq(-0.4,0.6,0.2),label=seq(-0.4,0.6,0.2),cex.axis=cexVal)
      segments(randInt[randInt[,"flag"]==0,"index"],randInt[randInt[,"flag"]==0,"low"],
               randInt[randInt[,"flag"]==0,"index"],randInt[randInt[,"flag"]==0,"upp"],
               lwd=lwdVal)
      segments(randInt[randInt[,"flag"]==1,"index"],randInt[randInt[,"flag"]==1,"low"],
               randInt[randInt[,"flag"]==1,"index"],randInt[randInt[,"flag"]==1,"upp"],
               lwd=lwdVal,col="red")
      points(randInt[randInt[,"flag"]==0,"index"],randInt[randInt[,"flag"]==0,"int"],
             col="black",pch=19)
      points(randInt[randInt[,"flag"]==1,"index"],randInt[randInt[,"flag"]==1,"int"],
             col="red",pch=19)
      legend("topleft",legend=c(paste("Hospital ID is",ii),paste("Sample size is",numObs)),
              cex=cexVal,text.col=rep("forestgreen",2),bty="n")
      legend("bottomright",legend=c(">0 : Higher than average", "=0: No different from average",
            "<0: Lower than average"),cex=cexVal,text.col=rep("navyblue",3),bty="n")

      # Figure 4: Heat map to display changes in hospital ranking after updates of
      # each hospital:

      randInt <- randInt[order(randInt[,"hosID"]),]
      rankMat <- cbind(rankMat,randInt[,"index"])
      rankMat2 <- data.frame(1:mMax,randInt[,"flag"],rankMat)
      colnames(rankMat2) <- "clusters"
      rank.profile(rankMat2)
      Sys.sleep(0.1)
      
      MFVB.plots[[ii]] <- recordPlot()
      if (!userwait) {dev.flush(); Sys.sleep(0.1)}
      if (ii<=3) wait()
      if (createPDF) dev.off()
   }
}   

if (createPDF) 
{
   pdf(paste(output.path,"mlevCSrealMFVB.pdf",sep=""),width=12)
   for (MFVB.plot in MFVB.plots)
   {
      replayPlot(MFVB.plot)
   }
   dev.off()
}   
   
######### End of mlevCSrealMFVB.R ##########


