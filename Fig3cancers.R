library(cowplot)

library(ggplot2)
library(grid)
library(gridExtra)
library(scales)

# Define Support Functions ------------------------------------------------




Convert1Year <- function(Vect) {
  as.integer(trimws(strsplit(
    x = Vect, split = "years", fixed = TRUE
  )))
}

GetTick <- function(Min, Max) {
  LowTick <- 10 ^ floor(log10(Min))
  LowTick <- floor(Min / LowTick) * LowTick

  HighTick <- 10 ^ floor(log10(Max))
  HighTick <- ceiling(Max / HighTick) * HighTick

  Pwr <- floor(log10(LowTick))
  Idx <- LowTick / (10 ^ Pwr)

  Ticks <- NULL
  TickLab <- NULL

  while (Idx * 10 ^ Pwr <= HighTick) {
    Ticks <- c(Ticks, Idx * 10 ^ Pwr)

    if (Idx == 1) {
      TickLab <- c(TickLab, bquote("10" ^  ~ .(Pwr)))
    } else {
      TickLab <- c(TickLab, "")
    }

    if (Idx < 9) {
      Idx <- Idx + 1
    } else {
      Idx <- 1
      Pwr <- Pwr + 1
    }

  }

  # TickLab[1] <- LowTick
  # TickLab[length(TickLab)] <- HighTick

  AllInfo <- list()
  AllInfo$Ticks <- Ticks
  AllInfo$TickLab <- as.vector(TickLab)

  return(AllInfo)

}





k_mid <- 0.044


# Other support functions -------------------------------------------------


Fit2ParModel <- function(FittedData, a=NULL, b=NULL, k_mid=NULL, Rounds=10){
  
  StartList <- list(a = runif(1), b = runif(1), k_mid = 0.05)
  
  StartList <- StartList[c(is.null(a), is.null(b), is.null(k_mid))]
  LowVect <- rep(1e-10, times=length(StartList))
  
  out2_m <-
    nls(
      log(y) ~ log(a / (exp(
        b * exp(-k_mid * x)
      ) - 1)),
      data = FittedData,
      start = StartList,
      algorithm = "port",
      lower = LowVect,
      control = nls.control(
        maxiter = 1000,
        tol = 1e-8,
        warnOnly = TRUE
      )
    )
  
  PEst <- summary(out2_m)$parameters[, "Estimate"]
  names(PEst) <- rownames(summary(out2_m)$parameters)
  StartList[names(PEst)] <- PEst
  
  if(Rounds > 0){
    for (i in 1:Rounds) {
      
      out2_m <-
        nls(
          log(y) ~ log(a / (exp(
            b * exp(-k_mid * x)
          ) - 1)),
          data = FittedData,
          start = StartList,
          lower = LowVect,
          algorithm = "port",
          control = nls.control(
            maxiter = 1000,
            tol = 1e-8/i,
            warnOnly = TRUE
          )
        )
      
      PEst <- summary(out2_m)$parameters[, "Estimate"]
      names(PEst) <- rownames(summary(out2_m)$parameters)
      StartList[names(PEst)] <- PEst
      
    }
  }
  
  R2Mb <- 1 - sum((summary(out2_m)$residuals)^2)/sum((log(FittedData$y)-mean(log(FittedData$y)))^2)
  
  PVals <- c(NA, NA, NA)
  names(PVals) <- c("a", "b", "k_mid")
  PVals[names(PEst)] <- PEst
  PVals[is.na(PVals)] <- c(a, b, k_mid)
  
  ReturnVect <- c(PVals, AIC(out2_m), R2Mb, deviance(out2_m), sum(summary(out2_m)$residuals^2))
  
  names(ReturnVect) <- c("a", "b", "k", "AIC", "R^2", "Dev", "Resid")
  
  return(ReturnVect)
  
}










Fit1ParModel <- function(FittedData, a=NULL, k_mid=NULL, Rounds=10){
  
  StartList <- list(a = runif(1), k_mid = 0.05)
  
  StartList <- StartList[c(is.null(a), is.null(k_mid))]
  LowVect <- rep(1e-10, times=length(StartList))
  
  out1_m <-
    nls(
      log(y) ~ log(a * exp(k_mid * x)),
      data = FittedData,
      start = StartList,
      lower = LowVect,
      algorithm = "port",
      control = nls.control(
        maxiter = 1000,
        tol = 1e-8,
        warnOnly = TRUE
      )
    )
  
  PEst <- summary(out1_m)$parameters[, "Estimate"]
  names(PEst) <- rownames(summary(out1_m)$parameters)
  StartList[names(PEst)] <- PEst
  
  for (i in 1:Rounds) {
    
    
    out1_m <-
      nls(
        log(y) ~ log(a * exp(k_mid * x)),
        data = FittedData,
        start = StartList,
        lower = LowVect,
        algorithm = "port",
        control = nls.control(
          maxiter = 1000,
          tol = 1e-8,
          warnOnly = TRUE
        )
      )
    
    PEst <- summary(out1_m)$parameters[, "Estimate"]
    names(PEst) <- rownames(summary(out1_m)$parameters)
    StartList[names(PEst)] <- PEst
    
  }
  
  
  R1Mb <- 1 - sum((summary(out1_m)$residuals)^2)/sum((log(FittedData$y)-mean(log(FittedData$y)))^2)
  
  PVals <- c(NA, NA)
  names(PVals) <- c("a", "k_mid")
  PVals[names(PEst)] <- PEst
  PVals[is.na(PVals)] <- c(a, k_mid)
  
  ReturnVect <- c(PVals, AIC(out1_m), R1Mb, deviance(out1_m), sum(summary(out1_m)$residuals^2))
  
  names(ReturnVect) <- c("a", "k", "AIC", "R^2", "Dev", "Resid")
  
  return(ReturnVect)
  
}

TryFitting <- function(Model, FittedData, a = NULL, b = NULL, k_mid = NULL,
                       t_0 = NULL, t_0_low = -Inf, t_0_high = Inf,
                       Rounds = 20, Retires = 10){
  
  AllFits <- NULL
  
  if(Model == "A"){
    for(i in 1:Retires){
      TR <- try(ModelPars <- Fit1ParModel(FittedData = FittedData, a = a, k_mid = k_mid, Rounds = Rounds))
      if(length(TR)!=1){
        AllFits <- rbind(AllFits, ModelPars)
      }
    }
  }
  
  
  if(Model == "B"){
    for(i in 1:Retires){
      TR <- try(ModelPars <- Fit2ParModel(FittedData = FittedData, a = a, b = b, k_mid = k_mid, Rounds = Rounds))
      if(length(TR)!=1){
        AllFits <- rbind(AllFits, ModelPars)
      }
    }
  }
  
  
  if(Model == "C"){
    for(i in 1:Retires){
      TR <- try(ModelPars <- FitPLModel(FittedData = FittedData, Rounds = Rounds))
      if(length(TR)!=1){
        AllFits <- rbind(AllFits, ModelPars)
      }
    }
  }
  
  
  if(Model == "TM"){
    for(i in 1:Retires){
      if(is.null(k_mid)){
        TR <- try(ModelPars <- FitTimeModel(FittedData = FittedData, a = a, t_0 = t_0, k = NULL,
                                            t_0_low = t_0_low, t_0_high = t_0_high, Rounds = Rounds))
      } else{
        TR <- try(ModelPars <- FitTimeModel(FittedData = FittedData, a = a, t_0 = t_0, k = -k_mid,
                                            t_0_low = t_0_low, t_0_high = t_0_high, Rounds = Rounds))
      }
      if(length(TR)!=1){
        AllFits <- rbind(AllFits, ModelPars)
      }
    }
  }
  
  return(AllFits[which.min(AllFits[,"Resid"]),])
  
}



Predict.Model <- function(Model, Pars, NewX){
  
  if(Model == "A"){
    Model.Pred <- Pars["a"] * exp(Pars["k"] * NewX)
  }
  
  if(Model == "B"){
    Model.Pred <- Pars["a"] / (exp(Pars["b"] * exp(-Pars["k"] * NewX)) - 1)
  }
  
  if(Model == "C"){
    Model.Pred <- exp(Pars["b"])*NewX^Pars["a"]
  }
  
  if(Model == "TM"){
    Model.Pred <- Pars["a"]/(exp(exp(Pars["k"]*(NewX-Pars["t_0"])))-1)
  }
  
  return(Model.Pred)
  
}


GetTick <- function(Min, Max) {
  LowTick <- 10 ^ floor(log10(Min))
  LowTick <- floor(Min / LowTick) * LowTick
  
  HighTick <- 10 ^ floor(log10(Max))
  HighTick <- ceiling(Max / HighTick) * HighTick
  
  Pwr <- floor(log10(LowTick))
  Idx <- LowTick / (10 ^ Pwr)
  
  Ticks <- NULL
  TickLab <- NULL
  
  while (Idx * 10 ^ Pwr <= HighTick) {
    Ticks <- c(Ticks, Idx * 10 ^ Pwr)
    
    if (Idx == 1) {
      TickLab <- c(TickLab, bquote("10" ^  ~ .(Pwr)))
    } else {
      TickLab <- c(TickLab, "")
    }
    
    if (Idx < 9) {
      Idx <- Idx + 1
    } else {
      Idx <- 1
      Pwr <- Pwr + 1
    }
    
  }
  
  # TickLab[1] <- LowTick
  # TickLab[length(TickLab)] <- HighTick
  
  AllInfo <- list()
  AllInfo$Ticks <- Ticks
  AllInfo$TickLab <- as.vector(TickLab)
  
  return(AllInfo)
  
}


Model.To.Use <- "B"






# Process Infection Data --------------------------------------------------




setwd("~/Downloads")

library(xlsx)


Dis1 <- read.xlsx("Diseases2.xlsx", sheetIndex = 1,stringsAsFactors=FALSE)

Dis1$Age=as.numeric(as.character(Dis1$Age))
Dis1$Inc=as.numeric(as.character(Dis1$Inc))
Dis1$LowerCI=as.numeric(as.character(Dis1$LowerCI))
Dis1$UpperCI=as.numeric(as.character(Dis1$UpperCI))
Dis1$Pop=as.numeric(as.character(Dis1$Pop))


InfToPlot=c(  
  # "MRSA"           ,             "West Nile Virus Disease"   ,           
            # "Streptococcus pneumoniae"   ,  "group B Streptococcus "    
           # ,    "Legionellosis"                 
            # ,"Haemophilus influenzae "  ,
                     )






# Process Cancer Data -----------------------------------------------------

library("SEER2R")

AllData <- read.SeerStat(DICfileName = "Incidence-Full-1Y-18.dic", TXTfileName = "Incidence-Full-1Y-18.txt", UseVarLabelsInData = TRUE)

# unique(AllData$`Site_recode_ICDO3/WHO_2008`)

AllData$AgeAdjusted_Rate <- as.numeric(levels(AllData$AgeAdjusted_Rate))[AllData$AgeAdjusted_Rate]
AllData$Lower_Confidence_Interval <- as.numeric(levels(AllData$Lower_Confidence_Interval))[AllData$Lower_Confidence_Interval]
AllData$Upper_Confidence_Interval <- as.numeric(levels(AllData$Upper_Confidence_Interval))[AllData$Upper_Confidence_Interval]



SumData <- NULL
PlotList <- list()

PlotList <- sapply(trimws(unique(AllData$`Site_recode_ICDO3/WHO_2008`)),function(x) NULL)

# uniscalif ---------------------------------------------------------------

k_mid <- 0.044


Newmybdf=data.frame()

Newmy2bdf=data.frame()


T1 <- read.csv("suppT1temp.csv", header=TRUE, stringsAsFactors=FALSE)

ToPlot=T1[,1][(T1$"IM.II.R.2" >= sort(T1$"IM.II.R.2",decreasing = TRUE)[20])]
ToPlot=T1[,1][(T1$"IM.II.R.2" >= sort(T1$"IM.II.R.2",decreasing = TRUE)[30])&(T1$"IM.II.AIC" <= sort(T1$"IM.II.AIC",decreasing = FALSE)[30])]

ToPlot=T1[,1][(T1$"IM.II.AIC" <= sort(T1$"IM.II.AIC",decreasing = FALSE)[20])]

ToPlot=sort(ToPlot,decreasing = TRUE)


# ToPlot <- c("Chronic Myeloid Leukemia",
#                "Soft Tissue including Heart",
#                "Brain",
#                "Colon and Rectum",
#                "Lip",
#                "Stomach",
#                "Digestive System",
#                "Salivary Gland",
#                "Descending Colon",
#                # "Rectum and Rectosigmoid Junction",
#                # "Intrahepatic Bile Duct",
#                # "Large Intestine, NOS",
#                # "Miscellaneous",
#                # "Mesothelioma",
#                "Gallbladder",
#                # "Aleukemic, Subleukemic and NOS",
#                # "Other Acute Leukemia",
#                # "Acute Monocytic Leukemia",
#                # "Other Lymphocytic Leukemia",
#                # "Leukemia",
#                "Myeloma",
#                "Non-Hodgkin Lymphoma",
#                "NHL - Extranodal",
#                "Other Urinary Organs",
#                "Vagina",
#                "Vulva",
#                "Urinary System",
#                
#                "Nose, Nasal Cavity and Middle Ear",
#                "Ovary",
#                "Ureter"
#                # "Lymphoma"
#                # "Penis"
# )
# 
# ToPlot <- c("Chronic Myeloid Leukemia",
#             "Soft Tissue including Heart",
#             "Brain",
#             "Colon and Rectum",
#             "Lip",
#             "Stomach",
#             "Digestive System",
#             "Salivary Gland",
#             "Brain and other nervous system",
#             # "Rectum and Rectosigmoid Junction",
#             # "Intrahepatic Bile Duct",
#             # "Large Intestine, NOS",
#             # "Miscellaneous",
#             # "Mesothelioma",
#             "Gallbladder",
#             # "Aleukemic, Subleukemic and NOS",
#             # "Other Acute Leukemia",
#             # "Acute Monocytic Leukemia",
#             # "Other Lymphocytic Leukemia",
#             # "Leukemia",
#             "Lymphoma",
#             "Non-Hodgkin Lymphoma",
#             "NHL - Extranodal",
#             "ONHL - nodal",
#             "Vagina",
#             "Vulva",
#             "Urinary System",
#             
#             "Melanoma of the skin",
#             "Rectum",
#             "Rectum and Rectosigmoid Junction"
#             # "Lymphoma"
#             # "Penis"
# )
# 
# 
# ToPlot <- c(
#   # "Rectum and Rectosigmoid Junction",
#                "Soft Tissue including Heart",
#                "Brain",
#                "Colon and Rectum",
#                "Lip",
#                "Stomach",
#                "Digestive System",
#                "Brain and Other Nervous System",
#                
#                # "Rectum",
#                # "Rectum and Rectosigmoid Junction",
#                # "Intrahepatic Bile Duct",
#                # "Large Intestine, NOS",
#                # "Miscellaneous",
#                # "Mesothelioma",
#                "Gallbladder",
#                # "Aleukemic, Subleukemic and NOS",
#                # "Other Acute Leukemia",
#                # "Acute Monocytic Leukemia",
#                # "Other Lymphocytic Leukemia",
#                # "Leukemia",
#                "Sigmoid Colon",
#                "Non-Hodgkin Lymphoma",
#                "NHL - Extranodal",
#                "Colon excluding Rectum",
#                
#                # "Esophagus",
#                "Vulva",
#                "Urinary System",
#                
#                # "Kidney and Renal Pelvis",
#                "Rectosigmoid Junction",
#                "Ureter",
#                "Esophagus"   ,
#                "Rectum and Rectosigmoid Junction",
#                "Rectum"           ,
#                "Kidney and Renal Pelvis"  
#                # "Lymphoma"
#                # "Penis"
# )

AllData$`Site_recode_ICDO3/WHO_2008`=trimws(AllData$`Site_recode_ICDO3/WHO_2008`)

# Newmybdf=rbind(Newmybdf,data.frame(cbind(-38:123,(exp(1)-1)/(exp(exp(-k_mid*(-38:123)))-1),"Theory")))
# colnames(Newmybdf)=c("Converted.Age"  ,  "AgeAdjusted_Rate", "type")

# Newmy2bdf=rbind(Newmy2bdf,data.frame(cbind(NA,NA,"Theory")))
# colnames(Newmy2bdf)=c("Converted.Age"  ,  "AgeAdjusted_Rate", "type")


# }


# for(Disease in InfToPlot){
  
#   
#   print(Disease)
#   
#   
#   XFil <- filter(Dis1,Name==Disease)$Age
#   YFil <- filter(Dis1,Name==Disease)$Inc
#   
#   Tks <-
#     GetTick(Min = (10 ^ -1.2) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)),
#             Max = (10 ^ 1.9) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)))
#   
#   ToFil <- YFil == 0 | XFil < 16 | is.na(XFil) | is.na(YFil)
#   
#   XFil <- XFil[!ToFil]
#   YFil <- YFil[!ToFil]
#   
#   FittedData <- data.frame(cbind(XFil, YFil))
#   colnames(FittedData) <- c("x", "y")
#   
#   ModData <- TryFitting(Model = Model.To.Use, FittedData = FittedData, k_mid = k_mid,
#                         Rounds = 4, Retires = 100)
#   
#   
#   # Model.Pred <- Predict.Model(Model = Model.To.Use, Pars = ModData, NewX = 18:92)
#   
#   # Model.Pred[XFil<12 | XFil>83] <- NA
#   
#   t0 <- log(abs(ModData["b"]))/k_mid
#   
#   DataSubset <- data.frame(cbind(XFil, YFil,Disease))
#   
#   colnames(DataSubset) <- c("Converted.Age", "AgeAdjusted_Rate","type")
#   
#   # inc0 = a2_m/(exp(b2_m*exp(-k_mid*t0))-1)
#   
#   if(t0<0){t0=-50}
#   
#   inc0=ModData["a"]/(exp(ModData["b"]*exp(-k_mid*(t0+100)))-1)/((exp(1)-1)/(exp(exp(-k_mid*(100)))-1))
#   
#   Newmydf=DataSubset
#   
#   Newmydf$Converted.Age=as.numeric(as.character(DataSubset$Converted.Age)) - t0
#   Newmydf$AgeAdjusted_Rate=as.numeric(as.character(Newmydf$AgeAdjusted_Rate))/inc0
#   Newmydf[as.numeric(as.character(Newmydf$AgeAdjusted_Rate))<= 0]=NA
#   Newmydf[as.numeric(as.character(Newmydf$Converted.Age))<=  17-t0]=NA
#   
#   
#   Newmybdf=rbind(Newmybdf,Newmydf)
#   
# 
#   Newmy2df=DataSubset
#   
#   Newmy2df$AgeAdjusted_Rate[as.numeric(as.character(Newmy2df$Converted.Age)) < 17] <- NA
#   
#   
#   Newmy2bdf=rbind(Newmy2bdf,Newmy2df)
#   
#   
# }


for(Cancer in ToPlot){
  
  # Cancer <- unique(AllData$`Site_recode_ICDO3/WHO_2008`)[25]
  
  if(Cancer == "All Sites"){
    next()
  }
  
  
  DataSubset <- AllData[AllData$`Site_recode_ICDO3/WHO_2008` == Cancer &
                          AllData$Race_recode_W_B_AI_API == "White" &
                          AllData$Sex == "Male and female",]
  
  if(trimws(DataSubset[1,]$`Site_recode_ICDO3/WHO_2008`)%in% ToPlot){
    
    
    print(Cancer)
    
    
    Tks <- GetTick(Min = max(c(0.1, min(DataSubset$Lower_Confidence_Interval, na.rm = TRUE))),
                   Max = max(DataSubset$Upper_Confidence_Interval, na.rm = TRUE))
    
    Converted.Age <- Convert1Year(DataSubset$Age_recode_with_single_ages_and_85)
    Converted.Age[Converted.Age==84] <- NA
    
    DataSubset <- cbind(DataSubset, Converted.Age)
    
    XFil <- DataSubset$Converted.Age
    YFil <- DataSubset$AgeAdjusted_Rate
    
    DataSubset <- DataSubset[-which(is.na(XFil)),]
    
    ToFil <- YFil == 0 | XFil < 18 | is.na(XFil) | is.na(YFil)
    
    XFil <- XFil[!ToFil]
    YFil <- YFil[!ToFil]
    
    out2_m <-
      nls(
        log(YFil) ~ log(a / (exp(
          b * exp(-k_mid * XFil)
        ) - 1)),
        start = list(a = 1, b = 1),
        algorithm = "port",
        lower = c(1e-10, 1e-10),
        control = nls.control(
          maxiter = 1000,
          tol = 1e-8,
          warnOnly = TRUE
        )
      )
    a2_m <- summary(out2_m)$parameters["a", "Estimate"]
    b2_m <- summary(out2_m)$parameters["b", "Estimate"]
    
    for (i in 1:10) {
      out2_m <-
        nls(
          log(YFil) ~ log(a / (exp(
            b * exp(-k_mid * XFil)
          ) - 1)),
          start = list(a = a2_m, b = b2_m),
          lower = c(1e-10, 1e-10),
          algorithm = "port",
          control = nls.control(
            maxiter = 1000,
            tol = 1e-8/i,
            warnOnly = TRUE
          )
        )
      a2_m <- summary(out2_m)$parameters["a", "Estimate"]
      b2_m <- summary(out2_m)$parameters["b", "Estimate"]
      
    }
    
    ModelB.Pred <- c(a2_m / (exp(b2_m * exp(-k_mid * DataSubset$Converted.Age)) - 1))
    ModelB.Pred[DataSubset$Converted.Age<18 & DataSubset$Converted.Age>83] <- NA
    DataSubset <- cbind(DataSubset, ModelB.Pred)
    # 
    # 
    #   out1_m <-
    #     nls(
    #       log(YFil) ~ log(a * exp(k_mid * XFil)),
    #       start = list(a = 0.1),
    #       lower = c(1e-10),
    #       algorithm = "port",
    #       control = nls.control(
    #         maxiter = 1000,
    #         tol = 1e-8,
    #         warnOnly = TRUE
    #       )
    #     )
    #   a1_m <- summary(out1_m)$parameters["a", "Estimate"]
    # 
    # 
    #   for (i in 1:10) {
    #     out1_m <-
    #       nls(
    #         log(YFil) ~ log(a * exp(k_mid * XFil)),
    #         start = list(a = a1_m),
    #         lower = c(1e-10),
    #         algorithm = "port",
    #         control = nls.control(
    #           maxiter = 1000,
    #           tol = 1e-8/i,
    #           warnOnly = TRUE
    #         )
    
    
    
    
    #       )
    #     a1_m <- summary(out1_m)$parameters["a", "Estimate"]
    # 
    #   }
    # 
    #   ModelA.Pred <- a1_m * (exp(k_mid * DataSubset$Converted.Age))
    #   ModelA.Pred[DataSubset$Converted.Age<18 & DataSubset$Converted.Age>83] <- NA
    #   DataSubset <- cbind(DataSubset, ModelA.Pred)
    # 
    #   out_pl <-
    #     nls(
    #       log(YFil) ~ a*log(XFil) + b,
    #       start = list(a = 1, b = 1),
    #       algorithm = "port",
    #       control = nls.control(
    #         maxiter = 1000,
    #         tol = 1e-8,
    #         warnOnly = TRUE
    #       )
    #     )
    #   a_pl <- summary(out_pl)$parameters["a", "Estimate"]
    #   b_pl <- summary(out_pl)$parameters["b", "Estimate"]
    # 
    # 
    #   for (i in 1:10) {
    #     out_pl <-
    #       nls(
    #         log(YFil) ~ a*log(XFil) + b,
    #         start = list(a = a_pl, b = b_pl),
    #         algorithm = "port",
    #         control = nls.control(
    #           maxiter = 1000,
    #           tol = 1e-8/i,
    #           warnOnly = TRUE
    #         )
    #       )
    #     a_pl <- summary(out_pl)$parameters["a", "Estimate"]
    #     b_pl <- summary(out_pl)$parameters["b", "Estimate"]
    # 
    #   }
    # 
    #   ModelC.Pred <- exp(b_pl)*(DataSubset$Converted.Age)^a_pl
    #   ModelC.Pred[DataSubset$Converted.Age<18 & DataSubset$Converted.Age>83] <- NA
    #   DataSubset <- cbind(DataSubset, ModelC.Pred)
    # 
    # 
    #   R2Ma <- 1 - sum((summary(out1_m)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
    R2Mb <- 1 - sum((summary(out2_m)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
    #   R2Mc <- 1 - sum((summary(out_pl)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
    # 
    #   SumData <- rbind(SumData, c(
    #     Cancer, k_mid,
    #     a1_m, AIC(out1_m), R2Ma,
    #     sum(summary(out1_m)$residuals^2),
    #     a2_m, b2_m, AIC(out2_m), R2Mb,
    #     sum(summary(out2_m)$residuals^2),
    #     a_pl, b_pl, AIC(out_pl), R2Mc,
    #     sum(summary(out_pl)$residuals^2)))
    # 
    #   BgkCol <- "yellow"
    # 
    #   if(R2Mb > 0.9 & R2Mc > 0.9 & AIC(out2_m)>AIC(out_pl)){
    #     BgkCol <- "yellow"
    #   }
    # 
    #   if(R2Mb < 0.85 & R2Mc < 0.85){
    #     BgkCol <- "gray"
    #   }
    # 
    #   if(R2Mb > 0.9 & R2Mc > 0.9 & AIC(out2_m)<AIC(out_pl)){
    #     BgkCol <- "yellow"
    #   }
    # 
    
    
    # 
    # 
    #   p <-
    #     ggplot(
    #       data = DataSubset,
    #       mapping = aes(
    #         x = Converted.Age,
    #         y = AgeAdjusted_Rate
    #       )
    #     )
    # 
    #   p <- p +  theme(
    #     plot.title=element_text(face="bold", size=25),
    #     axis.text=element_text(size=25),
    #     panel.grid.major = element_blank(),
    #     panel.grid.minor = element_blank(),
    #     panel.background = element_rect(fill = adjustcolor(BgkCol, alpha.f = 0.1)),
    #     plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), "cm")
    #   ) + labs(x = NULL, y = NULL) +
    #     ggtitle(trimws(Cancer)) +
    #   geom_point(color="blue")
    # 
    #   p <- p + scale_y_log10(
    #     breaks = Tks$Ticks,
    #     labels = Tks$TickLab,
    #     limits = range(Tks$Ticks)
    #   ) + scale_x_continuous(
    #     breaks = c(0, 40, 80),
    #     limits = c(0, 85)
    #     )
    # 
    #   p <- p  + geom_pointrange(color="blue",
    #       aes(
    #         ymin = Lower_Confidence_Interval,
    #         ymax = Upper_Confidence_Interval
    #       )
    #     )
    # 
    #   p <- p + geom_line(
    #     mapping = aes(x = Converted.Age, y = ModelB.Pred),
    #     color = adjustcolor("red", alpha.f = 1),
    #     size = 3
    #   ) + geom_line(
    #     mapping = aes(x = Converted.Age, y = ModelA.Pred),
    #     color = adjustcolor("red", alpha.f = 0.3),
    #     size = 2
    #   ) + geom_line(
    #     mapping = aes(x = Converted.Age, y = ModelC.Pred),
    #     color = adjustcolor("green", alpha.f = 0.8),
    #     size = 2
    #   ) 
    # 
    #   PlotList[[trimws(Cancer)]] <- p
    # 
    #   # print(length(PlotList))
    # 
    #   print(p)
    #   
    
    # if(R2Mb>0.98){
    
    t0 <- log(b2_m)/k_mid
    
    
    # inc0 = a2_m/(exp(b2_m*exp(-k_mid*t0))-1)
    
    if(t0< -50){t0=-50}
    
    inc0=a2_m/(exp(b2_m*exp(-k_mid*(t0+100)))-1)/((exp(1)-1)/(exp(exp(-k_mid*(100)))-1))
    
    
    Newmydf=data.frame(cbind(floor(DataSubset$Converted.Age-t0),DataSubset$AgeAdjusted_Rate/inc0,DataSubset$`Site_recode_ICDO3/WHO_2008`))
    
    Newmydf$X2[as.numeric(as.character(Newmydf$X2)) <= 0] <- NA
    
    Newmydf$X2[as.numeric(as.character(Newmydf$X1)) < 17-t0] <- NA
    
    colnames(Newmydf) <- c("Converted.Age", "AgeAdjusted_Rate","type")
    
    Newmybdf=rbind(Newmybdf,Newmydf)
    
    
    
    Newmy2df=data.frame(cbind(floor(DataSubset$Converted.Age),DataSubset$AgeAdjusted_Rate,DataSubset$`Site_recode_ICDO3/WHO_2008`))
    
    Newmy2df$X2[as.numeric(as.character(Newmy2df$X1)) < 17] <- NA
    
    colnames(Newmy2df) <- c("Converted.Age", "AgeAdjusted_Rate","type")
    
    Newmy2bdf=rbind(Newmy2bdf,Newmy2df)
    
    
    
  }
}

Newmybdf=cbind(Newmybdf,"Both")



colnames(Newmybdf)=c("time","risk","type","gender")

Theorydf=data.frame(cbind(-41:141,(exp(1)-1)/(exp(exp(-k_mid*(-41:141)))-1),"Theory - IM-II"))
colnames(Theorydf)=c("time","risk","type")

Theorydf2=data.frame(cbind(-41:141,(exp(1)-1)*exp(k_mid*(-41:141)),"Theory - IM-I"))
colnames(Theorydf2)=c("time","risk","type")


Newmyplot <-
 ggplot(Newmybdf,aes(x=as.numeric(as.character(time)),y=as.numeric(as.character(risk)),colour=type,group=type))+geom_line()+scale_y_log10( labels = comma,breaks=c(0.1,1,10,100))+xlim(-43,133)

Newmyplot = Newmyplot+geom_line(data=Theorydf,aes(x=as.numeric(as.character(time)),y=as.numeric(as.character(risk))),linetype = 2,size=1)
Newmyplot = Newmyplot+geom_line(data=Theorydf2,aes(x=as.numeric(as.character(time)),y=as.numeric(as.character(risk))),linetype = 2,size=1)

Newmyplot

colnames(Newmy2bdf)=c("time","risk","type")

Newmy2plot <-
  ggplot(Newmy2bdf,aes(x=as.numeric(as.character(time)),y=as.numeric(as.character(risk)),colour=as.character(type)))+geom_line(size=1)+scale_y_log10()+xlim(17,85)

Newmy2plot



ggsave(filename = "Fig4h2.ggplot.AIC20.pdf", plot = Newmyplot, width = 8, height = 8, scale = 2)
ggsave(filename = "Fig4g2.ggplot.AIC20.pdf", plot = Newmy2plot, width = 8, height = 8, scale = 2)


Expos <- c("Chronic Myeloid Leukemia",
            "Soft Tissue including Heart",
            "Brain",
            # "Colon and Rectum",
            # "Lip",
            # "Stomach",
            # "Digestive System",
            "Salivary Gland",
            # "Descending Colon",
            # "Gallbladder",
#              "Myeloma",
#             "Non-Hodgkin Lymphoma",
#             "NHL - Extranodal",
            # "Vagina",
            # "Vulva",
            # "Urinary System",
            # "Nose, Nasal Cavity and Middle Ear",
            # "Ovary",
            # "Ureter",
            "MRSA",
            "West Nile Virus Disease",     
"Streptococcus pneumoniae"   ,  "group B Streptococcus "    
            )
  



Newmy3bdf=filter(Newmy2bdf,type %in% Expos)

Newmy3plot <-
  ggplot(Newmy3bdf,aes(x=as.numeric(as.character(time)),y=as.numeric(as.character(risk)),colour=as.character(type)))+geom_line()+scale_y_log10()+xlim(17,85)

Newmy3plot

Newmy4bdf=filter(Newmy2bdf,!(type %in% Expos))

Newmy4plot <-
  ggplot(Newmy4bdf,aes(x=as.numeric(as.character(time)),y=as.numeric(as.character(risk)),colour=as.character(type)))+geom_line()+scale_y_log10()+xlim(17,85)

Newmy4plot

ggsave(filename = "Fig4g3.ggplot.pdf", plot = Newmy3plot, width = 8, height = 8, scale = 2)
ggsave(filename = "Fig4g4.ggplot.pdf", plot = Newmy4plot, width = 8, height = 8, scale = 2)

# male and female ---------------------------------------------------------

k_mid=0.038

Newmybdf=data.frame()


# 
# for(Disease in c("West Nile Virus Disease Female")){
#   
#   
#   print(Disease)
#   
#   
#   XFil <- filter(Dis1,Name==Disease)$Age
#   YFil <- filter(Dis1,Name==Disease)$Inc
#   
#   Tks <-
#     GetTick(Min = (10 ^ -1.2) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)),
#             Max = (10 ^ 1.9) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)))
#   
#   ToFil <- YFil == 0 | XFil < 16 | is.na(XFil) | is.na(YFil)
#   
#   XFil <- XFil[!ToFil]
#   YFil <- YFil[!ToFil]
#   
#   FittedData <- data.frame(cbind(XFil, YFil))
#   colnames(FittedData) <- c("x", "y")
#   
#   ModData <- TryFitting(Model = Model.To.Use, FittedData = FittedData, k_mid = k_mid,
#                         Rounds = 4, Retires = 100)
#   
#   
#   # Model.Pred <- Predict.Model(Model = Model.To.Use, Pars = ModData, NewX = 18:92)
#   
#   # Model.Pred[XFil<12 | XFil>83] <- NA
#   
#   t0 <- log(abs(ModData["b"]))/k_mid
#   
#   DataSubset <- data.frame(cbind(XFil, YFil,Disease))
#   
#   colnames(DataSubset) <- c("Converted.Age", "AgeAdjusted_Rate","type")
#   
#   # inc0 = a2_m/(exp(b2_m*exp(-k_mid*t0))-1)
#   
#   if(t0<0){t0=-50}
#   
#   inc0=ModData["a"]/(exp(ModData["b"]*exp(-k_mid*(t0+100)))-1)/((exp(1)-1)/(exp(exp(-k_mid*(100)))-1))
#   
#   Newmydf=DataSubset
#   
#   Newmydf$Converted.Age=as.numeric(as.character(DataSubset$Converted.Age)) - t0
#   Newmydf$AgeAdjusted_Rate=as.numeric(as.character(Newmydf$AgeAdjusted_Rate))/inc0
#   Newmydf[as.numeric(as.character(Newmydf$AgeAdjusted_Rate))<= 0]=NA
#   Newmydf[as.numeric(as.character(Newmydf$Converted.Age))<=  17-t0]=NA
#   
#   
#   
#   if(Disease=="West Nile Virus Disease Male"){ Newmydf=cbind(Newmydf,"Male")}
#   if(Disease=="West Nile Virus Disease Female"){ Newmydf=cbind(Newmydf,"Female")}
#   colnames(Newmydf)=c("time","risk","type","gender")
#   
#   
#   Newmybdf=rbind(Newmybdf,Newmydf)
#   
#   
# }

for(Cancer in ToPlot){
  
  # Cancer <- unique(AllData$`Site_recode_ICDO3/WHO_2008`)[25]
  
  if(Cancer == "All Sites"){
    next()
  }
  
  
  DataSubset <- AllData[AllData$`Site_recode_ICDO3/WHO_2008` == Cancer &
                          AllData$Race_recode_W_B_AI_API == "White" &
                          AllData$Sex == "  Female",]
  
  if(sum(DataSubset$AgeAdjusted_Rate, na.rm = TRUE)==0){
    
    SumData <- rbind(SumData, c(Cancer, NA, NA))
    next
  }
  
  if(trimws(DataSubset[1,]$`Site_recode_ICDO3/WHO_2008`)%in% ToPlot){
    
    
    print(Cancer)
    
    
    Tks <- GetTick(Min = max(c(0.1, min(DataSubset$Lower_Confidence_Interval, na.rm = TRUE))),
                   Max = max(DataSubset$Upper_Confidence_Interval, na.rm = TRUE))
    
    Converted.Age <- Convert1Year(DataSubset$Age_recode_with_single_ages_and_85)
    Converted.Age[Converted.Age==84] <- NA
    
    DataSubset <- cbind(DataSubset, Converted.Age)
    
    XFil <- DataSubset$Converted.Age
    YFil <- DataSubset$AgeAdjusted_Rate
    
    DataSubset <- DataSubset[-which(is.na(XFil)),]
    
    ToFil <- YFil == 0 | XFil < 18 | is.na(XFil) | is.na(YFil)
    
    XFil <- XFil[!ToFil]
    YFil <- YFil[!ToFil]
    
    out2_m <-
      nls(
        log(YFil) ~ log(a / (exp(
          b * exp(-k_mid * XFil)
        ) - 1)),
        start = list(a = 1, b = 1),
        algorithm = "port",
        lower = c(1e-10, 1e-10),
        control = nls.control(
          maxiter = 1000,
          tol = 1e-8,
          warnOnly = TRUE
        )
      )
    a2_m <- summary(out2_m)$parameters["a", "Estimate"]
    b2_m <- summary(out2_m)$parameters["b", "Estimate"]
    
    for (i in 1:10) {
      out2_m <-
        nls(
          log(YFil) ~ log(a / (exp(
            b * exp(-k_mid * XFil)
          ) - 1)),
          start = list(a = a2_m, b = b2_m),
          lower = c(1e-10, 1e-10),
          algorithm = "port",
          control = nls.control(
            maxiter = 1000,
            tol = 1e-8/i,
            warnOnly = TRUE
          )
        )
      a2_m <- summary(out2_m)$parameters["a", "Estimate"]
      b2_m <- summary(out2_m)$parameters["b", "Estimate"]
      
    }
    
#     ModData.F <- TryFitting(Model = "B", FittedData = FittedData,
#                             Rounds = 4, Retires = 50)
# 
#     a2_m <- ModData.F["a"]
#     b2_m <- ModData.F["b"]
    
    ModelB.Pred <- c(a2_m / (exp(b2_m * exp(-k_mid * DataSubset$Converted.Age)) - 1))
    ModelB.Pred[DataSubset$Converted.Age<18 & DataSubset$Converted.Age>83] <- NA
    DataSubset <- cbind(DataSubset, ModelB.Pred)
    

    # 
    # 
    #   out1_m <-
    #     nls(
    #       log(YFil) ~ log(a * exp(k_mid * XFil)),
    #       start = list(a = 0.1),
    #       lower = c(1e-10),
    #       algorithm = "port",
    #       control = nls.control(
    #         maxiter = 1000,
    #         tol = 1e-8,
    #         warnOnly = TRUE
    #       )
    #     )
    #   a1_m <- summary(out1_m)$parameters["a", "Estimate"]
    # 
    # 
    #   for (i in 1:10) {
    #     out1_m <-
    #       nls(
    #         log(YFil) ~ log(a * exp(k_mid * XFil)),
    #         start = list(a = a1_m),
    #         lower = c(1e-10),
    #         algorithm = "port",
    #         control = nls.control(
    #           maxiter = 1000,
    #           tol = 1e-8/i,
    #           warnOnly = TRUE
    #         )
    #       )
    #     a1_m <- summary(out1_m)$parameters["a", "Estimate"]
    # 
    #   }
    # 
    #   ModelA.Pred <- a1_m * (exp(k_mid * DataSubset$Converted.Age))
    #   ModelA.Pred[DataSubset$Converted.Age<18 & DataSubset$Converted.Age>83] <- NA
    #   DataSubset <- cbind(DataSubset, ModelA.Pred)
    # 
    #   out_pl <-
    #     nls(
    #       log(YFil) ~ a*log(XFil) + b,
    #       start = list(a = 1, b = 1),
    #       algorithm = "port",
    #       control = nls.control(
    #         maxiter = 1000,
    #         tol = 1e-8,
    #         warnOnly = TRUE
    #       )
    #     )
    #   a_pl <- summary(out_pl)$parameters["a", "Estimate"]
    #   b_pl <- summary(out_pl)$parameters["b", "Estimate"]
    # 
    # 
    #   for (i in 1:10) {
    #     out_pl <-
    #       nls(
    #         log(YFil) ~ a*log(XFil) + b,
    #         start = list(a = a_pl, b = b_pl),
    #         algorithm = "port",
    #         control = nls.control(
    #           maxiter = 1000,
    #           tol = 1e-8/i,
    #           warnOnly = TRUE
    #         )
    #       )
    #     a_pl <- summary(out_pl)$parameters["a", "Estimate"]
    #     b_pl <- summary(out_pl)$parameters["b", "Estimate"]
    # 
    #   }
    # 
    #   ModelC.Pred <- exp(b_pl)*(DataSubset$Converted.Age)^a_pl
    #   ModelC.Pred[DataSubset$Converted.Age<18 & DataSubset$Converted.Age>83] <- NA
    #   DataSubset <- cbind(DataSubset, ModelC.Pred)
    # 
    # 
    #   R2Ma <- 1 - sum((summary(out1_m)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
    R2Mb <- 1 - sum((summary(out2_m)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
    #   R2Mc <- 1 - sum((summary(out_pl)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
    # 
    #   SumData <- rbind(SumData, c(
    #     Cancer, k_mid,
    #     a1_m, AIC(out1_m), R2Ma,
    #     sum(summary(out1_m)$residuals^2),
    #     a2_m, b2_m, AIC(out2_m), R2Mb,
    #     sum(summary(out2_m)$residuals^2),
    #     a_pl, b_pl, AIC(out_pl), R2Mc,
    #     sum(summary(out_pl)$residuals^2)))
    # 
    #   BgkCol <- "yellow"
    # 
    #   if(R2Mb > 0.9 & R2Mc > 0.9 & AIC(out2_m)>AIC(out_pl)){
    #     BgkCol <- "yellow"
    #   }
    # 
    #   if(R2Mb < 0.85 & R2Mc < 0.85){
    #     BgkCol <- "gray"
    #   }
    # 
    #   if(R2Mb > 0.9 & R2Mc > 0.9 & AIC(out2_m)<AIC(out_pl)){
    #     BgkCol <- "yellow"
    #   }
    # 
    
    
    # 
    # 
    #   p <-
    #     ggplot(
    #       data = DataSubset,
    #       mapping = aes(
    #         x = Converted.Age,
    #         y = AgeAdjusted_Rate
    #       )
    #     )
    # 
    #   p <- p +  theme(
    #     plot.title=element_text(face="bold", size=25),
    #     axis.text=element_text(size=25),
    #     panel.grid.major = element_blank(),
    #     panel.grid.minor = element_blank(),
    #     panel.background = element_rect(fill = adjustcolor(BgkCol, alpha.f = 0.1)),
    #     plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), "cm")
    #   ) + labs(x = NULL, y = NULL) +
    #     ggtitle(trimws(Cancer)) +
    #   geom_point(color="blue")
    # 
    #   p <- p + scale_y_log10(
    #     breaks = Tks$Ticks,
    #     labels = Tks$TickLab,
    #     limits = range(Tks$Ticks)
    #   ) + scale_x_continuous(
    #     breaks = c(0, 40, 80),
    #     limits = c(0, 85)
    #     )
    # 
    #   p <- p  + geom_pointrange(color="blue",
    #       aes(
    #         ymin = Lower_Confidence_Interval,
    #         ymax = Upper_Confidence_Interval
    #       )
    #     )
    # 
    #   p <- p + geom_line(
    #     mapping = aes(x = Converted.Age, y = ModelB.Pred),
    #     color = adjustcolor("red", alpha.f = 1),
    #     size = 3
    #   ) + geom_line(
    #     mapping = aes(x = Converted.Age, y = ModelA.Pred),
    #     color = adjustcolor("red", alpha.f = 0.3),
    #     size = 2
    #   ) + geom_line(
    #     mapping = aes(x = Converted.Age, y = ModelC.Pred),
    #     color = adjustcolor("green", alpha.f = 0.8),
    #     size = 2
    #   ) 
    # 
    #   PlotList[[trimws(Cancer)]] <- p
    # 
    #   # print(length(PlotList))
    # 
    #   print(p)
    #   
    
    # if(R2Mb>0.98){
    
    t0 <- log(b2_m)/k_mid
    
    
    # inc0 = a2_m/(exp(b2_m*exp(-k_mid*t0))-1)
    
    if(t0< -50){t0=-50}
    
    inc0=a2_m/(exp(b2_m*exp(-k_mid*(t0+100)))-1)/((exp(1)-1)/(exp(exp(-k_mid*(100)))-1))
    
    
    Newmydf=data.frame(cbind(floor(DataSubset$Converted.Age-t0),DataSubset$AgeAdjusted_Rate/inc0,DataSubset$`Site_recode_ICDO3/WHO_2008`))
    
    Newmydf$X2[as.numeric(as.character(Newmydf$X2)) <= 0] <- NA
    
    Newmydf$X2[as.numeric(as.character(Newmydf$X1)) < 17-t0] <- NA
    Newmydf=cbind(Newmydf,"Female")
    colnames(Newmydf)=c("time","risk","type","gender")
    
    
    Newmybdf=rbind(Newmybdf,Newmydf)
    
    
    
    
    
    
    
  }
}


k_mid=0.046

# 
# for(Disease in c("West Nile Virus Disease Male")){
#   
#   
#   print(Disease)
#   
#   
#   XFil <- filter(Dis1,Name==Disease)$Age
#   YFil <- filter(Dis1,Name==Disease)$Inc
#   
#   Tks <-
#     GetTick(Min = (10 ^ -1.2) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)),
#             Max = (10 ^ 1.9) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)))
#   
#   ToFil <- YFil == 0 | XFil < 16 | is.na(XFil) | is.na(YFil)
#   
#   XFil <- XFil[!ToFil]
#   YFil <- YFil[!ToFil]
#   
#   FittedData <- data.frame(cbind(XFil, YFil))
#   colnames(FittedData) <- c("x", "y")
#   
#   ModData <- TryFitting(Model = Model.To.Use, FittedData = FittedData, k_mid = k_mid,
#                         Rounds = 4, Retires = 100)
#   
#   
#   # Model.Pred <- Predict.Model(Model = Model.To.Use, Pars = ModData, NewX = 18:92)
#   
#   # Model.Pred[XFil<12 | XFil>83] <- NA
#   
#   t0 <- log(abs(ModData["b"]))/k_mid
#   
#   DataSubset <- data.frame(cbind(XFil, YFil,Disease))
#   
#   colnames(DataSubset) <- c("Converted.Age", "AgeAdjusted_Rate","type")
#   
#   # inc0 = a2_m/(exp(b2_m*exp(-k_mid*t0))-1)
#   
#   if(t0<0){t0=-50}
#   
#   inc0=ModData["a"]/(exp(ModData["b"]*exp(-k_mid*(t0+100)))-1)/((exp(1)-1)/(exp(exp(-k_mid*(100)))-1))
#   
#   Newmydf=DataSubset
#   
#   Newmydf$Converted.Age=as.numeric(as.character(DataSubset$Converted.Age)) - t0
#   Newmydf$AgeAdjusted_Rate=as.numeric(as.character(Newmydf$AgeAdjusted_Rate))/inc0
#   Newmydf[as.numeric(as.character(Newmydf$AgeAdjusted_Rate))<= 0]=NA
#   Newmydf[as.numeric(as.character(Newmydf$Converted.Age))<=  17-t0]=NA
#   
#   
#   
#   if(Disease=="West Nile Virus Disease Male"){ Newmydf=cbind(Newmydf,"Male")}
#   if(Disease=="West Nile Virus Disease Female"){ Newmydf=cbind(Newmydf,"Female")}
#   colnames(Newmydf)=c("time","risk","type","gender")
#   
#   
#   Newmybdf=rbind(Newmybdf,Newmydf)
#   
#   
# }


for(Cancer in ToPlot){
  
  # Cancer <- unique(AllData$`Site_recode_ICDO3/WHO_2008`)[25]
  
  if(Cancer == "All Sites"){
    next()
  }
  
  
  DataSubset <- AllData[AllData$`Site_recode_ICDO3/WHO_2008` == Cancer &
                          AllData$Race_recode_W_B_AI_API == "White" &
                          AllData$Sex == "  Male",]
  
  if(sum(DataSubset$AgeAdjusted_Rate, na.rm = TRUE)==0){
    
    SumData <- rbind(SumData, c(Cancer, NA, NA))
    next
  }
  
  if(trimws(DataSubset[1,]$`Site_recode_ICDO3/WHO_2008`)%in% ToPlot){
    
    
    print(Cancer)
    
    
    Tks <- GetTick(Min = max(c(0.1, min(DataSubset$Lower_Confidence_Interval, na.rm = TRUE))),
                   Max = max(DataSubset$Upper_Confidence_Interval, na.rm = TRUE))
    
    Converted.Age <- Convert1Year(DataSubset$Age_recode_with_single_ages_and_85)
    Converted.Age[Converted.Age==84] <- NA
    
    DataSubset <- cbind(DataSubset, Converted.Age)
    
    XFil <- DataSubset$Converted.Age
    YFil <- DataSubset$AgeAdjusted_Rate
    
    DataSubset <- DataSubset[-which(is.na(XFil)),]
    
    ToFil <- YFil == 0 | XFil < 18 | is.na(XFil) | is.na(YFil)
    
    XFil <- XFil[!ToFil]
    YFil <- YFil[!ToFil]
    
    out2_m <-
      nls(
        log(YFil) ~ log(a / (exp(
          b * exp(-k_mid * XFil)
        ) - 1)),
        start = list(a = 1, b = 1),
        algorithm = "port",
        lower = c(1e-10, 1e-10),
        control = nls.control(
          maxiter = 1000,
          tol = 1e-8,
          warnOnly = TRUE
        )
      )
    a2_m <- summary(out2_m)$parameters["a", "Estimate"]
    b2_m <- summary(out2_m)$parameters["b", "Estimate"]
    
    for (i in 1:10) {
      out2_m <-
        nls(
          log(YFil) ~ log(a / (exp(
            b * exp(-k_mid * XFil)
          ) - 1)),
          start = list(a = a2_m, b = b2_m),
          lower = c(1e-10, 1e-10),
          algorithm = "port",
          control = nls.control(
            maxiter = 1000,
            tol = 1e-8/i,
            warnOnly = TRUE
          )
        )
      a2_m <- summary(out2_m)$parameters["a", "Estimate"]
      b2_m <- summary(out2_m)$parameters["b", "Estimate"]
      
    }
    
#     
#     ModData.M <- TryFitting(Model = "B", FittedData = FittedData,
#                             Rounds = 4, Retires = 50)
#     
#     a2_m <- ModData.M["a"]
#     b2_m <- ModData.M["b"]
    
    
    ModelB.Pred <- c(a2_m / (exp(b2_m * exp(-k_mid * DataSubset$Converted.Age)) - 1))
    ModelB.Pred[DataSubset$Converted.Age<18 & DataSubset$Converted.Age>83] <- NA
    DataSubset <- cbind(DataSubset, ModelB.Pred)
    # 
    # 
    #   out1_m <-
    #     nls(
    #       log(YFil) ~ log(a * exp(k_mid * XFil)),
    #       start = list(a = 0.1),
    #       lower = c(1e-10),
    #       algorithm = "port",
    #       control = nls.control(
    #         maxiter = 1000,
    #         tol = 1e-8,
    #         warnOnly = TRUE
    #       )
    #     )
    #   a1_m <- summary(out1_m)$parameters["a", "Estimate"]
    # 
    # 
    #   for (i in 1:10) {
    #     out1_m <-
    #       nls(
    #         log(YFil) ~ log(a * exp(k_mid * XFil)),
    #         start = list(a = a1_m),
    #         lower = c(1e-10),
    #         algorithm = "port",
    #         control = nls.control(
    #           maxiter = 1000,
    #           tol = 1e-8/i,
    #           warnOnly = TRUE
    #         )
    #       )
    #     a1_m <- summary(out1_m)$parameters["a", "Estimate"]
    # 
    #   }
    # 
    #   ModelA.Pred <- a1_m * (exp(k_mid * DataSubset$Converted.Age))
    #   ModelA.Pred[DataSubset$Converted.Age<18 & DataSubset$Converted.Age>83] <- NA
    #   DataSubset <- cbind(DataSubset, ModelA.Pred)
    # 
    #   out_pl <-
    #     nls(
    #       log(YFil) ~ a*log(XFil) + b,
    #       start = list(a = 1, b = 1),
    #       algorithm = "port",
    #       control = nls.control(
    #         maxiter = 1000,
    #         tol = 1e-8,
    #         warnOnly = TRUE
    #       )
    #     )
    #   a_pl <- summary(out_pl)$parameters["a", "Estimate"]
    #   b_pl <- summary(out_pl)$parameters["b", "Estimate"]
    # 
    # 
    #   for (i in 1:10) {
    #     out_pl <-
    #       nls(
    #         log(YFil) ~ a*log(XFil) + b,
    #         start = list(a = a_pl, b = b_pl),
    #         algorithm = "port",
    #         control = nls.control(
    #           maxiter = 1000,
    #           tol = 1e-8/i,
    #           warnOnly = TRUE
    #         )
    #       )
    #     a_pl <- summary(out_pl)$parameters["a", "Estimate"]
    #     b_pl <- summary(out_pl)$parameters["b", "Estimate"]
    # 
    #   }
    # 
    #   ModelC.Pred <- exp(b_pl)*(DataSubset$Converted.Age)^a_pl
    #   ModelC.Pred[DataSubset$Converted.Age<18 & DataSubset$Converted.Age>83] <- NA
    #   DataSubset <- cbind(DataSubset, ModelC.Pred)
    # 
    # 
    #   R2Ma <- 1 - sum((summary(out1_m)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
    R2Mb <- 1 - sum((summary(out2_m)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
    #   R2Mc <- 1 - sum((summary(out_pl)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
    # 
    #   SumData <- rbind(SumData, c(
    #     Cancer, k_mid,
    #     a1_m, AIC(out1_m), R2Ma,
    #     sum(summary(out1_m)$residuals^2),
    #     a2_m, b2_m, AIC(out2_m), R2Mb,
    #     sum(summary(out2_m)$residuals^2),
    #     a_pl, b_pl, AIC(out_pl), R2Mc,
    #     sum(summary(out_pl)$residuals^2)))
    # 
    #   BgkCol <- "yellow"
    # 
    #   if(R2Mb > 0.9 & R2Mc > 0.9 & AIC(out2_m)>AIC(out_pl)){
    #     BgkCol <- "yellow"
    #   }
    # 
    #   if(R2Mb < 0.85 & R2Mc < 0.85){
    #     BgkCol <- "gray"
    #   }
    # 
    #   if(R2Mb > 0.9 & R2Mc > 0.9 & AIC(out2_m)<AIC(out_pl)){
    #     BgkCol <- "yellow"
    #   }
    # 
    
    
    # 
    # 
    #   p <-
    #     ggplot(
    #       data = DataSubset,
    #       mapping = aes(
    #         x = Converted.Age,
    #         y = AgeAdjusted_Rate
    #       )
    #     )
    # 
    #   p <- p +  theme(
    #     plot.title=element_text(face="bold", size=25),
    #     axis.text=element_text(size=25),
    #     panel.grid.major = element_blank(),
    #     panel.grid.minor = element_blank(),
    #     panel.background = element_rect(fill = adjustcolor(BgkCol, alpha.f = 0.1)),
    #     plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), "cm")
    #   ) + labs(x = NULL, y = NULL) +
    #     ggtitle(trimws(Cancer)) +
    #   geom_point(color="blue")
    # 
    #   p <- p + scale_y_log10(
    #     breaks = Tks$Ticks,
    #     labels = Tks$TickLab,
    #     limits = range(Tks$Ticks)
    #   ) + scale_x_continuous(
    #     breaks = c(0, 40, 80),
    #     limits = c(0, 85)
    #     )
    # 
    #   p <- p  + geom_pointrange(color="blue",
    #       aes(
    #         ymin = Lower_Confidence_Interval,
    #         ymax = Upper_Confidence_Interval
    #       )
    #     )
    # 
    #   p <- p + geom_line(
    #     mapping = aes(x = Converted.Age, y = ModelB.Pred),
    #     color = adjustcolor("red", alpha.f = 1),
    #     size = 3
    #   ) + geom_line(
    #     mapping = aes(x = Converted.Age, y = ModelA.Pred),
    #     color = adjustcolor("red", alpha.f = 0.3),
    #     size = 2
    #   ) + geom_line(
    #     mapping = aes(x = Converted.Age, y = ModelC.Pred),
    #     color = adjustcolor("green", alpha.f = 0.8),
    #     size = 2
    #   ) 
    # 
    #   PlotList[[trimws(Cancer)]] <- p
    # 
    #   # print(length(PlotList))
    # 
    #   print(p)
    #   
    
    # if(R2Mb>0.98){
    
    t0 <- log(b2_m)/k_mid
    
    
    # inc0 = a2_m/(exp(b2_m*exp(-k_mid*t0))-1)
    
    if(t0< -50){t0=-50}
    
    inc0=a2_m/(exp(b2_m*exp(-k_mid*(t0+100)))-1)/((exp(1)-1)/(exp(exp(-k_mid*(100)))-1))
    
    
    Newmydf=data.frame(cbind(floor(DataSubset$Converted.Age-t0),DataSubset$AgeAdjusted_Rate/inc0,DataSubset$`Site_recode_ICDO3/WHO_2008`))
    
    Newmydf$X2[as.numeric(as.character(Newmydf$X2)) <= 0] <- NA
    
    Newmydf$X2[as.numeric(as.character(Newmydf$X1)) < 17-t0] <- NA
    Newmydf=cbind(Newmydf,"Male")
    colnames(Newmydf)=c("time","risk","type","gender")
    
    
    Newmybdf=rbind(Newmybdf,Newmydf)
    
    
    
    
    
    
    
  }
}


# plot --------------------------------------------------------------------

thEnd=135

Theorydf=data.frame(rbind(
  cbind(-36:thEnd,(exp(1)-1)/(exp(exp(-k_mid*(-36:thEnd)))-1),"Theory - IM-II","Male"),
  cbind(-49:thEnd,(exp(1)-1)/(exp(exp(-0.038*(-49:thEnd)))-1),"Theory - IM-II","Female")
)
)
colnames(Theorydf)=c("time","risk","type","gender")

Theorydf2=data.frame(rbind(
  cbind(-36:thEnd,(exp(1)-1)*exp(k_mid*(-36:thEnd)),"Theory - IM-I","Male"),
  cbind(-49:thEnd,(exp(1)-1)*exp(0.038*(-49:thEnd)),"Theory - IM-I","Female")
)
)
colnames(Theorydf2)=c("time","risk","type","gender")




Newmyplot <-
  ggplot(Newmybdf,aes(x=as.numeric(as.character(time)),y=as.numeric(as.character(risk)),colour=type,group=interaction(gender,type)))+geom_line()+scale_y_log10( labels = comma,breaks=c(0.1,1,10,100))+xlim(-51,thEnd)

Newmyplot = Newmyplot+geom_line(data=Theorydf,aes(x=as.numeric(as.character(time)),y=as.numeric(as.character(risk))),linetype = 2,size=1)
Newmyplot = Newmyplot+geom_line(data=Theorydf2,aes(x=as.numeric(as.character(time)),y=as.numeric(as.character(risk))),linetype = 2,size=1)

Newmyplot

ggsave(filename = "Fig4i.ggplot.R230.pdf", plot = Newmyplot, width = 8, height = 8, scale = 2)



# AllCancers --------------------------------------------------------------

k_mid=0.044

SumData <- NULL
PlotList <- list()

PlotList <- sapply(trimws(unique(AllData$`Site_recode_ICDO3/WHO_2008`)),function(x) NULL)


for(Cancer in unique(AllData$`Site_recode_ICDO3/WHO_2008`)){
  
  # Cancer <- unique(AllData$`Site_recode_ICDO3/WHO_2008`)[25]
  
  if(Cancer == "All Sites"){
    next()
  }
  
  
  DataSubset <- AllData[AllData$`Site_recode_ICDO3/WHO_2008` == Cancer &
                          AllData$Race_recode_W_B_AI_API == "White" &
                          AllData$Sex == "Male and female",]
  
  if(trimws(DataSubset[1,]$`Site_recode_ICDO3/WHO_2008`)%in% ToPlot){
    
    
    print(Cancer)
    
    
    Tks <- GetTick(Min = max(c(0.1, min(DataSubset$Lower_Confidence_Interval, na.rm = TRUE))),
                   Max = max(DataSubset$Upper_Confidence_Interval, na.rm = TRUE))
    
    Converted.Age <- Convert1Year(DataSubset$Age_recode_with_single_ages_and_85)
    Converted.Age[Converted.Age==84] <- NA
    
    DataSubset <- cbind(DataSubset, Converted.Age)
    
    XFil <- DataSubset$Converted.Age
    YFil <- DataSubset$AgeAdjusted_Rate
    
    DataSubset <- DataSubset[-which(is.na(XFil)),]
    
    ToFil <- YFil == 0 | XFil < 18 | is.na(XFil) | is.na(YFil)
    
    XFil <- XFil[!ToFil]
    YFil <- YFil[!ToFil]
    
    out2_m <-
      nls(
        log(YFil) ~ log(a / (exp(
          b * exp(-k_mid * XFil)
        ) - 1)),
        start = list(a = 1, b = 1),
        algorithm = "port",
        lower = c(1e-10, 1e-10),
        control = nls.control(
          maxiter = 1000,
          tol = 1e-8,
          warnOnly = TRUE
        )
      )
    a2_m <- summary(out2_m)$parameters["a", "Estimate"]
    b2_m <- summary(out2_m)$parameters["b", "Estimate"]
    
    for (i in 1:10) {
      out2_m <-
        nls(
          log(YFil) ~ log(a / (exp(
            b * exp(-k_mid * XFil)
          ) - 1)),
          start = list(a = a2_m, b = b2_m),
          lower = c(1e-10, 1e-10),
          algorithm = "port",
          control = nls.control(
            maxiter = 1000,
            tol = 1e-8/i,
            warnOnly = TRUE
          )
        )
      a2_m <- summary(out2_m)$parameters["a", "Estimate"]
      b2_m <- summary(out2_m)$parameters["b", "Estimate"]
      
    }
    
    ModelB.Pred <- c(a2_m / (exp(b2_m * exp(-k_mid * DataSubset$Converted.Age)) - 1))
    ModelB.Pred[DataSubset$Converted.Age<18 | DataSubset$Converted.Age>83] <- NA
    DataSubset <- cbind(DataSubset, ModelB.Pred)
    
    
      out1_m <-
        nls(
          log(YFil) ~ log(a * exp(k_mid * XFil)),
          start = list(a = 0.1),
          lower = c(1e-10),
          algorithm = "port",
          control = nls.control(
            maxiter = 1000,
            tol = 1e-8,
            warnOnly = TRUE
          )
        )
      a1_m <- summary(out1_m)$parameters["a", "Estimate"]
    
    
      for (i in 1:10) {
        out1_m <-
          nls(
            log(YFil) ~ log(a * exp(k_mid * XFil)),
            start = list(a = a1_m),
            lower = c(1e-10),
            algorithm = "port",
            control = nls.control(
              maxiter = 1000,
              tol = 1e-8/i,
              warnOnly = TRUE
            )
          )
        a1_m <- summary(out1_m)$parameters["a", "Estimate"]
    
      }
    
      ModelA.Pred <- a1_m * (exp(k_mid * DataSubset$Converted.Age))
      ModelA.Pred[DataSubset$Converted.Age<18 | DataSubset$Converted.Age>83] <- NA
      DataSubset <- cbind(DataSubset, ModelA.Pred)
    
      out_pl <-
        nls(
          log(YFil) ~ a*log(XFil) + b,
          start = list(a = 1, b = 1),
          algorithm = "port",
          control = nls.control(
            maxiter = 1000,
            tol = 1e-8,
            warnOnly = TRUE
          )
        )
      a_pl <- summary(out_pl)$parameters["a", "Estimate"]
      b_pl <- summary(out_pl)$parameters["b", "Estimate"]
    
    
      for (i in 1:10) {
        out_pl <-
          nls(
            log(YFil) ~ a*log(XFil) + b,
            start = list(a = a_pl, b = b_pl),
            algorithm = "port",
            control = nls.control(
              maxiter = 1000,
              tol = 1e-8/i,
              warnOnly = TRUE
            )
          )
        a_pl <- summary(out_pl)$parameters["a", "Estimate"]
        b_pl <- summary(out_pl)$parameters["b", "Estimate"]
    
      }
    
      ModelC.Pred <- exp(b_pl)*(DataSubset$Converted.Age)^a_pl
      ModelC.Pred[DataSubset$Converted.Age<18 | DataSubset$Converted.Age>83] <- NA
      DataSubset <- cbind(DataSubset, ModelC.Pred)
    
    
      R2Ma <- 1 - sum((summary(out1_m)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
    R2Mb <- 1 - sum((summary(out2_m)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
      R2Mc <- 1 - sum((summary(out_pl)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
    
      SumData <- rbind(SumData, c(
        Cancer, k_mid,
        a1_m, AIC(out1_m), R2Ma,
        sum(summary(out1_m)$residuals^2),
        a2_m, b2_m, AIC(out2_m), R2Mb,
        sum(summary(out2_m)$residuals^2),
        a_pl, b_pl, AIC(out_pl), R2Mc,
        sum(summary(out_pl)$residuals^2)))
    
      BgkCol <- "white"
    
#       if(R2Mb > 0.9 & R2Mc > 0.9 & AIC(out2_m)>AIC(out_pl)){
#         BgkCol <- "yellow"
#       }
#     
#       if(R2Mb < 0.85 & R2Mc < 0.85){
#         BgkCol <- "gray"
#       }
#     
#       if(R2Mb > 0.9 & R2Mc > 0.9 & AIC(out2_m)<AIC(out_pl)){
#         BgkCol <- "yellow"
#       }
#     
    
    
    
    
      p <-
        ggplot(
          data = DataSubset,
          mapping = aes(
            x = Converted.Age,
            y = AgeAdjusted_Rate
          )
        )
    
      p <- p +  theme(
        plot.title=element_text(face="bold", size=25),
        axis.text=element_text(size=25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = adjustcolor(BgkCol, alpha.f = 0.1)),
        plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), "cm")
      ) + labs(x = NULL, y = NULL) +
        ggtitle(trimws(Cancer)) +
      geom_point(color="blue")
    
      p <- p + scale_y_log10(
        breaks = Tks$Ticks,
        labels = Tks$TickLab,
        limits = range(Tks$Ticks)
      ) + scale_x_continuous(
        breaks = c(0, 45, 90),
        limits = c(0, 90)
        )
    
      p <- p  + geom_pointrange(color="blue",
          aes(
            ymin = Lower_Confidence_Interval,
            ymax = Upper_Confidence_Interval
          )
        )
    
      p <- p + geom_line(
        mapping = aes(x = Converted.Age, y = ModelB.Pred),
        color = adjustcolor("red", alpha.f = 1),
        size = 3
#       ) + geom_line(
#         mapping = aes(x = Converted.Age, y = ModelA.Pred),
#         color = adjustcolor("red", alpha.f = 0.3),
#         size = 2
      ) + geom_line(
        mapping = aes(x = Converted.Age, y = ModelC.Pred),
        color = adjustcolor("green", alpha.f = 0.8),
        size = 2
      ) 
      
    
      PlotList[[trimws(Cancer)]] <- p
    
    
    
    
  }
}



NewPlotList <- list(
  PlotList$`Chronic Myeloid Leukemia`,
  PlotList$`Soft Tissue including Heart`,
  PlotList$Brain,
  PlotList$`Colon and Rectum`,
                    PlotList$Lip,
                    PlotList$Stomach)

ml=plot_grid(plotlist=NewPlotList, ncol = 3)

ggsave(filename = "Fig3a.pdf", plot = ml, width = 8,  scale = 2)


