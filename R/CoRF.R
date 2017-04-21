
#' @title Calculate the oob-auc of a rfsrc object
#' @description Calculate the oob-auc of a rfsrc object (from randomForestSRC) using rfsrc$predicticed.oob and the binary response variable.
#' @param rf A rfsrc object.
#' @return Returns the oob-auc.
#' @export

AUC.randomForest <- function(rf)
{
    r <- rank(rf$predicted.oob[,2])
    rd <- mean(r[rf$class == levels(rf$class)[2]])
    nd <- sum(rf$class == levels(rf$class)[2])
    nnd <- length(rf$class) - nd
    return((rd - nd/2 - 0.5)/nnd)
}

#' @title Calculate the auc of RF predictions
#' @description Calculate the auc of a RF using a set of RF prediction. Primarly useful for cross-validated predictions or predictions made on independent data.
#' @param PredResp A vector of predictions (fraction of votes) and the predicted response vector binded together, e.g. cbind(Prediction,Response).
#' @return Returns the auc.
#' @export

AUC.randomForest.Pred <- function(PredResp) 	#input: cbind(Pred,Resp)
{
    if(all(sort(unique(PredResp[,2]))==c(1,2))) PredResp[,2] <- PredResp[,2]-1
    if(!all(unique(sort(PredResp[,2]))==c(0,1))) print("unknown input for AUC.randomForest.Pred")

    r <- rank(PredResp[,1])
    rd <- mean(r[PredResp[,2] == 1])
    nd <- sum(PredResp[,2] == 1)
    nnd <- length(PredResp[,2]) - nd
    return((rd - nd/2 - 0.5)/nnd)
}


#' @title CoRF: improved high-dimensional prediction with the RF by the use of co-data.
#' @description This fits a RF guided by co-data (CoRF).
#' @param Y A response variable.
#' @param X The primary set of variables (n rows, and p colums).
#' @param CoData A data.frame containing the co-data, p rows and one columns per set of co-data.
#' @param CoDataModelText Optionally a text string containing the specification of the co-data model. Not needed if CoDataRelation is specified.
#' @param CoDataRelation The to be fitted relationship in the co-data model, e.g. linear, (monotome) increasing, (monotome) decreasing. Alternatively, choose one of the scam smooth contrained contructs (mpd, mpi, mdcv, mdcx, micv, micx, cv, cx, see ??shape.constrained.smooth.terms).
#' @param ScreenCoData Boolean that indicates whether or not co-data selection step should be conducted Default is TRUE.
#' @param ScreenPvalThreshold The threshold value used in the co-data selection step (when ScreenCoData is FALSE). If TRUE it Defaults to 0.05/CD, where CD is the number of co-data sets.
#' @param TuneGamma If TRUE this sets GammaSeq to c(1/3, 2/3, 0.9, 1, 1.1, 1.2, 1.3). CoRF is fitted for each value of gamma, and returns results for all refited forests. Default is FALSE, in which case GammaSeq = 1.
#' @param GammaSeq Specifies the sequence of gamma values. Default is to only use a gamma of 1.
#' @param BaseRF Optionally use an earlier fitted base RF (i.e. uniform sampling probabilities).
#' @param ForestsInOutput Whether or not to save the various rfsrc objects. These objects are needed to make further predictions. Defaults to TRUE. Optionally set to FALSE if the CoRF objects become too big and are not needed.
#' @param setseed seed used to fit the RF.
#' @param importance default set to "none". Not needed to fit CoRF and very computational expensive.
#' @param nodesize Sets the minimal node size, set at recommended default for CoRF (2).
#' @param ntree The number of trees used in CoRF, default set to 2000. Convergences improves for a larger of number trees. If set too low the co-data model will give a poor fit.
#' @return Returns an CoRF object with the following components.
#' \item{InitialCall}{Details of the statements used when calling CoRF.}
#' \item{SamplingProbabilities}{The sampling probabilities used in refitting. A matrix of size p x number of gammas.}
#' \item{ScreenPval}{The screening p-values. Each p-value is the result of glm fit per type of co-data (univariable) to the variables used. For monotome increasing/decreasing relationships p-value is one-sided.}
#' \item{SavedForests}{A list of the fitted RFs. Element [[1]] contains the base RF, subsequent elements contain the refitted RFs. See attr(SavedForests[[1]],"WhichForest").}
#' \item{CoDataModel}{Fit of the co-data model.}
#' \item{ResultPerGamma}{Result overview of the fit of the base RF and of the refitted RFs.}
#' @section Reference: CoRF paper
#' @author authors
#' @examples
#' #---Run CoRF:
#' #load(LNM_Example)
#' #CoDataRelation <- c("increasing","linear","decreasing")
#' #CoRF_Fit <- CoRF(Y=RespTrain,X=TrainData,CoData=CoDataTrain,CoDataRelation=CoDataRelation,ntree=2000)
#'
#' #---These elements contains the rfsrc objects
#' #CoRF_Fit$SavedForests contains the rfsrc objects
#' #CoRF_Fit$SavedForests[[1]]  #The base RF
#' #CoRF_Fit$SavedForests[[2]]  #Subsequent numbers contain CoRF fits
#'
#' #---Overview of results, per gamma. The first row is the base RF:
#' #CoRF_Fit$ResultPerGamma
#'
#' #---SamplingProbabilities used to refit:
#' #CoRF_Fit$SamplingProbabilities
#'
#' #---Co-data model:
#' #CoRF_Fit$CoDataModel
#'
#' #---Plot fit of the co-data model:
#' #plot(CoRF_Fit$CoData$Corrs,plogis(predict(CoRF_Fit$CoDataModel)))
#'
#' #---The second option to run CoRF is through specifying the co-data model with CoDataModelText
#' #Example to run CoRF through CoDataModelText
#'
#' #---first take care of the missing values.
#' #CoData$pvalsVUmc[is.na(CoData$pvalsVUmc)] <- mean(CoData$pvalsVUmc,na.rm=TRUE)
#'
#' #---specify the co-data model:
#' #CoDataModelText <- "~ s(Corrs,k=25,bs=\"mpi\",m=2)+RoepmanGenes+s(pvalsVUmc,k=25,bs=\"mpd\",m=2)"
#' #CoRF_Fit <- CoRF(Y=RespTrain,X=TrainData,CoData=CoDataTrain,CoDataModelText=CoDataModelText)
#'
#' #---The third way to run CoRF is directly through randomForestSRC and scam (or glm).
#' #DF <- data.frame(Ydf=RespTrain,Xdf=TrainData)
#' #Forest <- rfsrc(Ydf ~ .,data=DF,ntree=2000,var.used="all.trees",importance=c("none"),nodesize=2,seed=1)
#' #CoDataTrain$pvalsVUmc[is.na(CoDataTrain$pvalsVUmc)] <- mean(CoDataTrain$pvalsVUmc,na.rm=TRUE)
#' #CoDataModell <- scam(VarUsed/sum(VarUsed)~  s(Corrs,k=25,bs="mpi",m=2)+RoepmanGenes+s(pvalsVUmc,k=25,bs="mpd",m=2),data=CoDataTrain,family=quasibinomial)
#' #preds <- as.numeric(plogis(predict(CoDataModell)))
#' #P <- length(preds)
#' #preds2 <- pmax(preds2-1/P,0)
#' #Mtry <- ceiling(sqrt(sum(preds2!=0)))
#' #ReffitedCoRF <- rfsrc(Ydf ~ .,data=DF,ntree=2000,var.used="all.trees",importance=c("none"),xvar.wt=preds2,mtry=Mtry,nodesize=2,setseed=1)


#' @export

CoRF <- function(Y,X,CoData,CoDataModelText=NULL,CoDataRelation=NULL,ScreenCoData=FALSE,ScreenPvalThreshold=NULL,TuneGamma=FALSE,GammaSeq=NULL,GammaNVar=NULL,BaseRF=NULL,ForestsInOutput=TRUE,setseed=1,importance=c("none"),nodesize=2,ntree=2000,...)
{
    Y <- factor(Y)
    DF <- data.frame(Ydf=Y,Xdf=X)
    P <- dim(X)[2]
    N <- dim(X)[1]
    CD <- dim(CoData)[2]

    if(is.null(ScreenPvalThreshold)) ScreenPvalThreshold <- 0.10/CD       #divide by CD for bonferroni correction
    if(class(CoData)=="matrix") CoData <- data.frame(CoData)
    if(is.null(colnames(CoData))) colnames(CoData) <- gsub(" ","",paste("CoData",1:CD))
    if(!is.null(colnames(CoData))) if(sum(is.na(colnames(CoData)))>0) colnames(CoData) <- gsub(" ","",paste("CoData",1:CD))
    if(is.null(setseed)) setseed <- sample(1E7)
    if(is.na(setseed)) setseed <- sample(1E7)
    if(!is.null(GammaSeq)) if(sum(is.na(GammaSeq))>0) stop("NA's in GammaSeq")
    if(!is.null(GammaNVar)) if(sum(is.na(GammaNVar))>0) stop("NA's in GammaNVar")

    mtry <- ceiling(sqrt(P))
    xvar.wt <- rep(1/P,times=P)

    if(is.null(CoDataModelText))
    {
        if(is.null(CoDataRelation))
        {
            print("No vector with CoDataRelations and no CoDataModelText found, fitting linear effects with glm")
            CoDataRelation <- rep("linear",times=CD)
        } else {
            CoDataRelation <- replace(CoDataRelation,CoDataRelation=="increasing","mpi")
            CoDataRelation <- replace(CoDataRelation,CoDataRelation=="decreasing","mpd")
        }
    }


    #####################################################################################################################
    #--Train base forest
    #####################################################################################################################


    if(is.null(BaseRF)==TRUE)
    {
        print("train base RF")

        #Forest <- randomForestSRC::rfsrc(Ydf ~ .,data=DF,importance=importance,nodesize=nodesize,ntree=ntree,var.used=var.used,xvar.wt=xvar.wt,mtry=mtry,...)
#ptm <- proc.time()
        Forest <- randomForestSRC::rfsrc(Ydf ~ .,data=DF,importance=importance,nodesize=nodesize,ntree=ntree,xvar.wt=NULL,mtry=mtry,var.used="all.trees",seed=setseed)
#print(proc.time() - ptm)
        VarUsed <- as.numeric(Forest$var.used)

        print("Base RF fitted")

    } else {
        Forest <- BaseRF
        VarUsed <- as.numeric(Forest$var.used)
    }


    ######################################################################################################################
    #--Train Co-data model
    ######################################################################################################################


    if(!is.null(CoDataModelText))
    {
        modelPart1 <- c("VarUsed/sum(VarUsed)")
        model <- formula(paste(modelPart1,CoDataModelText))

        if(sum(is.na(CoData))==0)
        {
            model <- formula(paste(modelPart1,CoDataModelText))

            print("Fitting co-data model  ")
            print(paste(modelPart1,CoDataModelText))
            match <- nchar(CoDataModelText)-stringdist::stringdist(CoDataModelText,"s(")

            if(match>=2)
            {
                print("fitting with scam")

                CoDataModel <- scam::scam(model,data=CoData,family=quasibinomial)
                CoDataNA <- CoData

                print("  ")
                print("scam co-data model fitted in")
                print(CoDataModel$CPU.time)
                print("  ")
            } else {

                print("fitting with glm")

                CoDataModel <- glm(model,data=CoData,family=quasibinomial)
            }
        }
        if(sum(is.na(CoData))>0)
        {
            print("When CoDataModelText is specified, the user needs to take of NA's, see example")
            stop("Codata contains NA's")
        }
    }


    ######################################################################################################################
    ###Screen and select co-data, takes only place if CoDataModelText is Null
    ######################################################################################################################


    if(ScreenCoData==TRUE & is.null(CoDataModelText))
    {
        ScreenPval <- numeric(CD)+NA

        for(cd in 1:CD)
        {
            #cd <- 3

            Subset <- !is.na(CoData[,cd])
            modelPart1 <- c("VarUsed[Subset]/sum(VarUsed[Subset]) ~")
            textmodel <- paste(colnames(CoData)[cd])
            model <- formula(paste(modelPart1,textmodel))
            ScreenModel <- glm(model,data=CoData[Subset,],family=quasibinomial)

            if(class(CoData[,cd]) == "numeric" | class(CoData[,cd]) == "integer")
            {
               if(CoDataRelation[cd] == "mpd" | CoDataRelation[cd] == "mdcx" | CoDataRelation[cd] == "mdcv") ScreenPval[cd] <- pt(summary(ScreenModel)$coeff[2,3],df=summary(ScreenModel)$df.residual)
               if(CoDataRelation[cd] == "mpi" | CoDataRelation[cd] == "micx" | CoDataRelation[cd] == "micv") ScreenPval[cd] <- 1-pt(summary(ScreenModel)$coeff[2,3],df=summary(ScreenModel)$df.residual)
               if(CoDataRelation[cd] == "linear" | CoDataRelation[cd] == "cx" | CoDataRelation[cd] == "cv") ScreenPval[cd] <- summary(ScreenModel)$coeff[2,4]
            }
            if(class(CoData[,cd]) == "factor")
            {
               ScreenPval[cd] <- anova(ScreenModel, test="Chisq")[2,5]
            }
       }

        OmitCoData <- ScreenPval>ScreenPvalThreshold
        CoDataNA <- CoData[,!OmitCoData,drop=FALSE]

        print("removing co-data:")
        for(h in 1:sum(OmitCoData)) print(paste(colnames(CoData)[OmitCoData][h],"- screenings pvalue",round(ScreenPval[OmitCoData],3)[h]))

        CoDataRelationNA <- CoDataRelation[!OmitCoData]
        names(ScreenPval) <- colnames(CoData)
    }

    if(ScreenCoData==FALSE | !is.null(CoDataModelText))
    {
        CoDataNA <- CoData
        ScreenPval <- NULL
    }


    ######################################################################################################################
    ###Create models when CoDataModelText is Null
    ###Also takes care of missing values
    ######################################################################################################################


    CD2 <- dim(CoDataNA)[2]
    if(CD2==0) print("no codata selected")

    if(is.null(CoDataModelText) & CD2>0)
    {
        CollecText <- c()

        for(cd in 1:CD2)
        {
            #cd <- 1
            N_NA <- sum(is.na(CoDataNA[,cd]))

            if(CoDataRelationNA[cd] == "linear")
            {
                if(N_NA==0)        #No missing values
                {
                    textmodel <- paste(colnames(CoDataNA)[cd])
                }

                if(N_NA>0 & N_NA<10)  #Less than 10 NAs, just impute the mean
                {
                    print(paste("Found less than 10 missing values in:",colnames(CoDataNA)[cd]))
                    print("Imputing with mean")

                    MissingCD <- is.na(CoDataNA[,cd])
                    CoDataNA[MissingCD,cd] <- mean(CoDataNA[!MissingCD,cd])
                    textmodel <- paste(colnames(CoDataNA)[cd])
                }

                if(N_NA>=10)      #More than 10 NAs, model NA's with extra factor
                {
                    print(paste("Found at least 10 missing values in:",colnames(CoDataNA)[cd]))
                    print("Adding extra factor to co-data model to model these missing values")

                    textNA <- gsub(" ","",paste(colnames(CoDataNA)[cd],"Missing"))
                    MissingCD <- is.na(CoDataNA[,cd])
                    CoDataNA[MissingCD,cd] <- mean(CoDataNA[!MissingCD,cd])
                    CoDataNA <- cbind(CoDataNA,factor(MissingCD))
                    colnames(CoDataNA)[dim(CoDataNA)[2]] <- textNA

                    textmodel <- paste(colnames(CoDataNA)[cd])
                    textmodel <- gsub(" ","",textmodel)
                    textmodel <- c(textmodel,textNA)
                }
            }

            if(CoDataRelationNA[cd] == "mpi" | CoDataRelationNA[cd] == "mpd")
            {
                if(N_NA==0)      #No missing values
                {
                    textmodel <- paste("s(",colnames(CoDataNA)[cd],",k=25,bs=\"",CoDataRelationNA[cd],"\",m=2)")
                    textmodel <- gsub(" ","",textmodel)
                }

                if(N_NA>0 & N_NA<10)   #Less than 10 NAs, just impute the mean
                {
                    print(paste("Found less than 10 missing values in:",colnames(CoDataNA)[cd]))
                    print("Imputing with mean")

                    MissingCD <- is.na(CoDataNA[,cd])
                    CoDataNA[MissingCD,cd] <- mean(CoDataNA[!MissingCD,cd])

                    textmodel <- paste("s(",colnames(CoDataNA)[cd],",k=25,bs=\"",CoDataRelationNA[cd],"\",m=2)")
                    textmodel <- gsub(" ","",textmodel)
                }

                if(N_NA>=10)       #More than 10 NAs, model NA's with extra factor
                {
                    print(paste("Found at least 10 missing values in:",colnames(CoDataNA)[cd]))
                    print("Adding extra factor to co-data model to model these missing values")

                    textNA <- gsub(" ","",paste(colnames(CoDataNA)[cd],"Missing"))
                    MissingCD <- is.na(CoDataNA[,cd])
                    CoDataNA[MissingCD,cd] <- mean(CoDataNA[!MissingCD,cd])
                    CoDataNA <- cbind(CoDataNA,factor(MissingCD))
                    colnames(CoDataNA)[dim(CoDataNA)[2]] <- textNA

                    textmodel <- paste("s(",colnames(CoDataNA)[cd],",k=25,bs=\"",CoDataRelationNA[cd],"\",m=2,by=",textNA,")")
                    textmodel <- gsub(" ","",textmodel)
                }
            }

            CollecText <- c(CollecText,textmodel)
        }

        ModelText <- paste(" ~",paste(CollecText,collapse="+"))
        modelPart1 <- c("VarUsed/sum(VarUsed)")
        model <- formula(paste(modelPart1,ModelText))

        print("Fitting co-data model")
        print(paste(modelPart1,ModelText))

        fitGLM <- all(CoDataRelationNA == "linear")
        if(fitGLM)
        {
            CoDataModel <- glm(model,data=CoDataNA,family=quasibinomial)
        }
        if(!fitGLM)
        {
            CoDataModel <- scam::scam(model,data=CoDataNA,family=quasibinomial)

            print("  ")
            print("scam co-data model fitted in")
            print(CoDataModel$CPU.time)
            print("  ")
        }
    }

    ##################################################################################################
    ###Get sampling probs and check if all values for gamma make sense..
    ##################################################################################################

    if(CD2>0)
    {
        preds0 <- plogis(predict(CoDataModel))

        if(is.null(GammaSeq))
        {
            if(is.null(GammaNVar))
            {
                if(TuneGamma) GammaSeq <- c(1/3,2/3,0.9,1,1.1,1.2,1.3)
                if(!TuneGamma) GammaSeq <- 1
            } else {
                quants <- quantile(preds0,probs=rev(c(1-GammaNVar/P)))
                GammaSeq <- quants/(1/P)
            }
        }

        NSteps <- length(GammaSeq)
        predsCand <- matrix(nrow=P,ncol=NSteps)

        maxprob <- mean(sort(unique(preds0),decreasing=TRUE)[1:2])
        UpperGamma <- max(1,maxprob/(1/P))

        RemoveGamma <- numeric(NSteps)+0
        RemoveGamma[GammaSeq > UpperGamma] <- 1

        if(sum(RemoveGamma==0)>0)
        {
            for(s in 1:NSteps)
            {
                tmp <- pmax(preds0-GammaSeq[s]*1/P,0)

                if(sum(tmp!=0)>0 & RemoveGamma[s]==0)
                {
                   predsCand[,s] <- tmp/sum(tmp)
                   if(s > 1) if(sum(round(predsCand[,s-1],10)!=round(predsCand[,s],10))==0) RemoveGamma[s] <- 1
                } else {
                   RemoveGamma[s] <- 1
                }
            }

            predsCand2 <- predsCand[,RemoveGamma==0,drop=FALSE]
            GammaSeq2 <- GammaSeq[RemoveGamma==0]
            NSteps2 <- length(GammaSeq2)
            NUsed <- colSums(predsCand2!=0)
        }

        if(sum(RemoveGamma==0)==0)
        {
            print(paste("Maximum gamma value is",round(UpperGamma,3)))
            print("Gamma(s) supplied higher than maximum")
            print("Only base RF returned")
            NSteps2 <- 0
        }
        if(sum(RemoveGamma==1)>1)
        {
            print(paste("Maximum gamma value is",round(UpperGamma,3)))
            print(paste("Removing gamma values",paste(round(GammaSeq[RemoveGamma==1],1),collapse=",")))
        }
    }

    if(CD2==0) NSteps2 <- 0


    ##################################################################################################
    ###Initialise variables for output saving, and some output saving
    ##################################################################################################


    SavedForests <- list()
    ActiveVariables=Correct=Errors=SaveAUC=Brier=UsedVariables <- numeric(NSteps2+1)+NA

    SavedForests[[1]] <- Forest
    attr(SavedForests[[1]],"WhichForest") <- "Base Forest"
    if(nlevels(Y)==2)
    {
        SaveAUC[1] <- AUC.randomForest(Forest)
        tab <- table(Y,Forest$predicted.oob[,2]>0.5)
        Errors[1] <- 1-sum(diag(tab))/sum(tab)
        Correct[1] <- sum(diag(tab))
        YBin <- as.numeric(Y)-1
        Brier[1] <- mean((Forest$predicted.oob[,2]-YBin)^2)
    }

    ActiveVariables[1] <- P
    UsedVariables[1] <- sum(as.numeric(Forest$var.used)!=0)

    print("  ")
    print("-----------------------------------------")
    print("  ")
    print("Base RF:")
    if(nlevels(Y)==2)
    {
        printout <- data.frame(round(SaveAUC[1],3),round(Errors[1],3),round(Brier[1],3),Correct[1],ActiveVariables[1],UsedVariables[1])
        names(printout) <- c("AUC","Error","Brier","Correct","ActiveVariables","UsedVariables")
        print(printout)
    }
    print("  ")
    print("-----------------------------------------")
    print("  ")


    ##################################################################################################
    ###Retrain the forest with different values of lambda
    ##################################################################################################

    EnterNextLoop <- TRUE
    if(CD2>0) if(sum(RemoveGamma==0)>0)
    {
        text1 <- "train random forest with gamma "
        text2 <- "gamma "
        savetext <- character(NSteps2)
        ActiveVariables[2:(NSteps2+1)] <- NUsed

        for(f in 1:NSteps2)
        {
            print(paste(text1,round(GammaSeq2[f],3)))
            savetext[f] <- paste(text2,round(GammaSeq2[f],3))
            preds <- predsCand2[,f]

            #SavedForests[[f+1]] <- randomForestSRC::rfsrc(Ydf ~ .,data=DF,importance=importance,nodesize=nodesize,ntree=ntree,var.used=var.used,xvar.wt=preds,mtry=ceiling(sqrt(NUsed[f])),...)
#ptm <- proc.time()
            SavedForests[[f+1]] <- randomForestSRC::rfsrc(Ydf ~ .,data=DF,importance=importance,nodesize=nodesize,ntree=ntree,var.used="all.trees",xvar.wt=preds,mtry=ceiling(sqrt(NUsed[f])),seed=setseed)
#print(proc.time() - ptm)
            attr(SavedForests[[f+1]],"WhichForest") <- paste("Forest with gamma",GammaSeq2[f])

            if(nlevels(Y)==2)
            {
                tab <- table(Y,SavedForests[[f+1]]$predicted.oob[,2]>0.5)
                SaveAUC[f+1] <- AUC.randomForest(SavedForests[[f+1]])
                Errors[f+1] <- 1-sum(diag(tab))/sum(tab)
                Correct[f+1] <- sum(diag(tab))
                Brier[f+1] <- mean((SavedForests[[f+1]]$predicted.oob[,2]-YBin)^2)
            }

            UsedVariables[f+1] <- sum(as.numeric(SavedForests[[f+1]]$var.used)!=0)

            if(nlevels(Y)==2)
            {
                printout <- data.frame(round(SaveAUC[f+1],3),round(Errors[f+1],3),round(Brier[f+1],3),Correct[f+1],ActiveVariables[f+1],UsedVariables[f+1])
                names(printout) <- c("AUC","Error","Brier","Correct","ActiveVariables","UsedVariables")
                print(printout)
            }

            print("  ")
            print("-----------------------------------------")
            print("  ")
        }

        ################################################################################################################################
        ###Save outputs
        #################################################################################################################################

        SamplingProbabilities <- predsCand2
        PredictedProbs <- preds0
        colnames(SamplingProbabilities) <- savetext
        outMat <- data.frame(Gamma=c(NA,GammaSeq2),AUC=SaveAUC,OOB=Errors,NCorrect=Correct,NFalse=N-Correct,Brier=Brier,ActiveVariables=ActiveVariables,UsedVariables=UsedVariables)

        EnterNextLoop <- FALSE
    }

    if(EnterNextLoop)
    {
        print("Returning only the base RF")
        SamplingProbabilities <- NA
        CoDataModel <- NA
        outMat <- data.frame(Gamma=NA,AUC=SaveAUC,OOB=Errors,NCorrect=Correct,NFalse=N-Correct,Brier=Brier,ActiveVariables=ActiveVariables,UsedVariables=UsedVariables)
    }

    if(!ForestsInOutput) SavedForests <- NULL
    InitialCall <- list(CoData=CoData,Y=Y,CoDataModelText=CoDataModelText,CoDataRelation=CoDataRelation,ScreenCoData=ScreenCoData,TuneGamma=TuneGamma,GammaSeq=GammaSeq,GammaNVar=GammaNVar,ntree=ntree,nodesize=nodesize,setseed=setseed)
    return(list(InitialCall=InitialCall,SamplingProbabilities=SamplingProbabilities,ScreenPval=ScreenPval,SavedForests=SavedForests,CoDataModel=CoDataModel,ResultPerGamma=outMat))
}


#' @title n-fold cross-validation of CoRF
#' @description Conducts a n-fold cross-validation of a CoRF object to obtain an unbiased estimate of the predictive error. Primarly for when gamma is tuned.
#' @param CoRFObject A CoRF object.
#' @param nfold The cross-validation fold.
#' @param SeedCVO The seed used to draw the cross-validation folds.
#' @param GammaSeq Set gamma sequence to be used in the cross-validation. By default uses the same sequence as used to fit the CoRF object.
#' @param ntree Set the number of trees used per RF in the cross-validation. By default uses the same number of trees as used to fit the CoRF object.
#' @return Returns An object with the result of the cross-validation, containng the following components.
#' \item{foldid}{The cross-validation folds used}
#' \item{OptGammaPerFold}{When more than one gamma was used, this returns per fold the gamma with lowest brier/auc.}
#' \item{Predictions}{The cross-validated predictions, for the base RF, the CoRF with gamma=1, and optionally the predictions for the CoRF selected by brier or auc.}
#' \item{cvResult}{The result of the cross-validation, e.g. cross-validated auc, brier, and accuracy.}
#' @export


cvCoRF <- function(CoRFObject,nfold=10,SeedCVO=123,GammaSeq=NULL,ntree=NULL,SelectForest=TRUE)
{
    X <- CoRFObject$SavedForests[[1]]$xvar
    Y <- CoRFObject$InitialCall$Y
    N <- dim(X)[1]

    set.seed(SeedCVO)
    cvo <- c(rep(1:nfold,each=floor(N/nfold)),sample(1:nfold,N-floor(N/nfold)*nfold,replace=FALSE))[sample(1:N,N,replace=FALSE)]

    Predictions <- matrix(nrow=N,ncol=4)

    if(is.null(GammaSeq)) GammaSeq <- CoRFObject$InitialCall$GammaSeq
    if(is.null(ntree)) ntree <- CoRFObject$InitialCall$ntree
    OptGammaPerFold <- matrix(nrow=nfold,ncol=2)
    cvResultPerGamma=cvPreds <- list()

    for(f in 1:nfold)
    {
        #f <- 1

        print("##############################################################################")
        print("nfold")
        print(f)

        index <- which(cvo==f)

        XTrain <- X[-index,]
        YTrain <- Y[-index]

        XTest <- X[index,,drop=FALSE]
        YTest <- Y[index]

        cvCoRFObject <- CoRF(Y=YTrain,X=XTrain,CoData=CoRFObject$InitialCall$CoData,CoDataModelText=CoRFObject$InitialCall$CoDataModelText,
        CoDataRelation=CoRFObject$InitialCall$CoDataRelation,ScreenCoData=CoRFObject$InitialCall$ScreenCoData,TuneGamma=CoRFObject$InitialCall$TuneGamma,
        GammaSeq=GammaSeq,nodesize=CoRFObject$InitialCall$nodesize,ntree=ntree,setseed=CoRFObject$InitialCall$setseed)

        DFTest <- data.frame(Ydf=factor(YTest),Xdf=XTest)


        ####################################################################################################


        predBase <- predict(cvCoRFObject$SavedForests[[1]],newdata=DFTest,importance=c("none"))$predicted
        rownames(predBase) <- which(cvo==f)
        if(f==1) cvPreds[[1]] <- predBase
        if(f>1) cvPreds[[1]] <- rbind(cvPreds[[1]],predBase)

        cvResultPerGamma[[f]] <- cvCoRFObject$ResultPerGamma

        if(f==nfold) attr(cvPreds[[1]],"WhichForest") <- "Base Forest"

        for(g in 1:length(GammaSeq))
        {
            whichgamma <- (cvCoRFObject$ResultPerGamma$Gamma[-1]==GammaSeq[g])

            if(sum(whichgamma)==1)
            {
                predGamma <- predict(cvCoRFObject$SavedForests[[which(whichgamma)+1]],newdata=DFTest,importance=c("none"))$predicted
                rownames(predGamma) <- which(cvo==f)
            } else {
                predGamma <- predBase
                predGamma[] <- NA
            }

            if(f==1) cvPreds[[g+1]] <- predGamma

            if(f>1) cvPreds[[g+1]] <- rbind(cvPreds[[g+1]],predGamma)

            if(f==nfold) attr(cvPreds[[g+1]],"WhichForest") <- paste("Gamma",round(GammaSeq[g],2))
        }


        ####################################################################################################


        if(SelectForest & nlevels(Y)==2)
        {
            WhichMaxAUC <- which.max(cvCoRFObject$ResultPerGamma$AUC)
            WhichMinBrier <- which.min(cvCoRFObject$ResultPerGamma$Brier)
            OptGammaPerFold[f,] <- cbind(cvCoRFObject$ResultPerGamma$Gamma[WhichMaxAUC],cvCoRFObject$ResultPerGamma$Gamma[WhichMinBrier])

            currentfold <- (sum(cvo<f)+1):sum(cvo<=f)

            if(f==1)
            {
                PredMaxAUC <- cvPreds[[WhichMaxAUC]][currentfold,]
                PredMinBrier <- cvPreds[[WhichMinBrier]][currentfold,]
            }
            if(f>1)
            {
                PredMaxAUC <- rbind(PredMaxAUC,cvPreds[[WhichMaxAUC]][currentfold,])
                PredMinBrier <- rbind(PredMinBrier,cvPreds[[WhichMinBrier]][currentfold,])
            }
        }
    }


    ####################################################################################################


    if(SelectForest & nlevels(Y)==2)
    {
        cvPreds[[length(cvPreds)+1]] <- PredMaxAUC
        attr(cvPreds[[length(cvPreds)]],"WhichForest") <- paste("Max AUC")
        cvPreds[[length(cvPreds)+1]] <- PredMinBrier
        attr(cvPreds[[length(cvPreds)]],"WhichForest") <- paste("Min Brier")
    }

    for(g in 1:length(cvPreds))
    {
        o <- order(as.numeric(rownames(cvPreds[[g]])))
        tmp <- attr(cvPreds[[g]],"WhichForest")
        cvPreds[[g]] <- cvPreds[[g]][o,]
        attr(cvPreds[[g]],"WhichForest") <- tmp
    }

    cvAUC=cvBrier=cvError <- numeric()

    if(nlevels(Y)==2)
    {
        YBin <- as.numeric(Y)-1

        for(g in 1:length(cvPreds))
        {
            if(sum(is.na(cvPreds[[g]]))==0)
            {
                cvAUC[g] <- AUC.randomForest.Pred(cbind(cvPreds[[g]][,2],YBin))
                cvBrier[g] <- mean((cvPreds[[g]][,2]-YBin)^2)
                tabBase <- table(Y,cvPreds[[g]][,2]>0.5)
                cvError[g] <- 1-sum(diag(tabBase))/sum(tabBase)
            }
            if(sum(is.na(cvPreds[[g]]))>0)
            {
                cvAUC[g]=cvBrier[g]=cvError[g] <- NA
            }
            names(cvAUC)[g]=names(cvBrier)[g]=names(cvError)[g] <- attr(cvPreds[[g]],"WhichForest")
        }
    }

    if(nlevels(Y)==2)
    {
        cvResult <- cbind(cvBrier,cvError,cvAUC)
        out <- list(foldid=cvo,OptGammaPerFold=OptGammaPerFold,Y=Y,cvPreds=cvPreds,cvResult=cvResult)
    }

    if(nlevels(Y)>2) out <- list(foldid=cvo,Y=Y,cvPreds=cvPreds,cvResultPerGamma=cvResultPerGamma)
    return(out)
}
