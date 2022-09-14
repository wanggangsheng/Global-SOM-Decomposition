library(pacman)
p_load(ggplot2,readxl,openxlsx,ggsci,ggmap,maps,tidyverse,magrittr,dplyr,forestplot,grid,cowplot,ggpubr,
       lattice,ggridges,raster,rasterVis,RColorBrewer,caret,gbm,pdp,Hmisc,GGally)

#=============
R2 <- function(obs,sim){
  obs_mean <- mean(obs)
  R2 <-  1-sum((sim-obs)^2)/sum((obs-obs_mean)^2)
  R2
}
RMSE <- function(obs,sim){
  n <- length(obs)
  MSE <- (sum((sim-obs)^2))/n
  RMSE <- MSE^0.5
  RMSE
}
AIC <- function(obs,sim,p){
  ssr <- sum((sim-obs)^2)
  n <- length(obs)
  AIC <- n*log(ssr/n) + 2*p                              #p is the number of parameters, n is the length of observations
  AIC
}
#=============

#create new folder, save all results
ifelse(!dir.exists("globalresults_final"), dir.create("globalresults_final"), FALSE)
IndexsumData <- read.csv("arranged.csv")
inputData <- read.xlsx("dataForML.xlsx")
#define font sizes
font_size1 <- 6    #legend font size
font_size2 <- 7    #text font size
font_size3 <- 8    #title font size
anno_size <- 2.5   #annotate text size
#redord statistical values
STATVALUES <- data.frame(min = c(0), first_Qu = c(0), median = c(0), 
                         mean = c(0), third_Qu = c(0), max = c(0), 
                         CI95_lower = c(0), CI95_higher = c(0),
                         CI90_lower = c(0), CI90_higher = c(0),
                         CI50_lower = c(0), CI50_higher = c(0))
#1.contrast the distribution of parameters among different pool-model and field/incubation condition
experimentType <- c("incubation","field_experiment")
modelType <- paste(c("one","two","three"), rep("-pool",3), sep = "")
typicalModel <- c("ANIMO","DAISY","CLM","DAYCENT")
incuStartPosition <- which(inputData$labORfield == experimentType[1])
fieldStartPosition <- which(inputData$labORfield == experimentType[2])
onepoolStartPosition <- which(inputData$MODELS == modelType[1])
twopoolStartPosition <- which(inputData$MODELS == modelType[2])
threepoolStartPosition <- which(inputData$MODELS == modelType[3])
k_onepoolIncu <- inputData$`k1(day-1)`[(onepoolStartPosition+1):(twopoolStartPosition[1]-1)]
k_twopoolIncu <- c(inputData$`k1(day-1)`[(twopoolStartPosition[1]+1):(threepoolStartPosition[1]-1)],
                   inputData$`k2(day-1)`[(twopoolStartPosition[1]+1):(threepoolStartPosition[1]-1)])
k_threepoolIncu <- c(inputData$`k1(day-1)`[(threepoolStartPosition[1]+1):(twopoolStartPosition[2]-1)],
                     inputData$`k2(day-1)`[(threepoolStartPosition[1]+1):(twopoolStartPosition[2]-1)],
                     inputData$`k3(day-1)`[(threepoolStartPosition[1]+1):(twopoolStartPosition[2]-1)])
k_twopoolField <- c(inputData$`k1(day-1)`[(twopoolStartPosition[2]+1):(threepoolStartPosition[2]-1)],
                    inputData$`k2(day-1)`[(twopoolStartPosition[2]+1):(threepoolStartPosition[2]-1)])
k_threepoolField <- c(inputData$`k1(day-1)`[(threepoolStartPosition[2]+1):nrow(inputData)],
                      inputData$`k2(day-1)`[(threepoolStartPosition[2]+1):nrow(inputData)],
                      inputData$`k3(day-1)`[(threepoolStartPosition[2]+1):nrow(inputData)])
#set k value equals 0 to NA
k_onepoolIncu[k_onepoolIncu == 0] <- NA
k_twopoolIncu[k_twopoolIncu == 0] <- NA
k_threepoolIncu[k_threepoolIncu == 0] <- NA
k_twopoolField[k_twopoolField == 0] <- NA
k_threepoolField[k_threepoolField == 0] <- NA
f_onepoolIncu <- inputData$`f1(%)`[(onepoolStartPosition+1):(twopoolStartPosition[1]-1)]
f_twopoolIncu <- c(inputData$`f1(%)`[(twopoolStartPosition[1]+1):(threepoolStartPosition[1]-1)],
                   inputData$`f2(%)`[(twopoolStartPosition[1]+1):(threepoolStartPosition[1]-1)])
f_threepoolIncu <- c(inputData$`f1(%)`[(threepoolStartPosition[1]+1):(twopoolStartPosition[2]-1)],
                     inputData$`f2(%)`[(threepoolStartPosition[1]+1):(twopoolStartPosition[2]-1)],
                     inputData$`f3(%)`[(threepoolStartPosition[1]+1):(twopoolStartPosition[2]-1)])
f_twopoolField <- c(inputData$`f1(%)`[(twopoolStartPosition[2]+1):(threepoolStartPosition[2]-1)],
                    inputData$`f2(%)`[(twopoolStartPosition[2]+1):(threepoolStartPosition[2]-1)])
f_threepoolField <- c(inputData$`f1(%)`[(threepoolStartPosition[2]+1):nrow(inputData)],
                      inputData$`f2(%)`[(threepoolStartPosition[2]+1):nrow(inputData)],
                      inputData$`f3(%)`[(threepoolStartPosition[2]+1):nrow(inputData)])
k_valuesIncu <- c(k_onepoolIncu, k_twopoolIncu, k_threepoolIncu)
k_valuesField <- c(k_twopoolField, k_threepoolField)
f_valuesIncu <- c(f_onepoolIncu, f_twopoolIncu, f_threepoolIncu)
f_valuesField <- c(f_twopoolField, f_threepoolField)
k_incuVsField <- data.frame(rateConstant = c(k_valuesIncu, k_valuesField),
                            model = c(rep("one-pool", length(k_onepoolIncu)),
                                      rep("two-pool", each = length(k_twopoolIncu)),
                                      rep("three-pool", each = length(k_threepoolIncu)),
                                      rep("two-pool", each = length(k_twopoolField)),
                                      rep("three-pool", each = length(k_threepoolField))),
                            pools = c(rep("fast_onepool", length(k_onepoolIncu)),
                                      rep(paste0(c("fast","slow"),"_twopool"), each = length(k_twopoolIncu)/2),
                                      rep(paste0(c("fast","slow","passive"),"_threepool"), each = length(k_threepoolIncu)/3),
                                      rep(paste0(c("fast","slow"),"_twopool"), each = length(k_twopoolField)/2),
                                      rep(paste0(c("fast","slow","passive"),"_threepool"), each = length(k_threepoolField)/3)),
                            experiment = c(rep("incubation", length(k_valuesIncu)),
                                           rep("field_experiment", length(k_valuesField))))
f_incuVsField <- data.frame(poolSize = c(f_valuesIncu, f_valuesField),
                            model = c(rep("one-pool", length(f_onepoolIncu)),
                                      rep("two-pool", each = length(f_twopoolIncu)),
                                      rep("three-pool", each = length(f_threepoolIncu)),
                                      rep("two-pool", each = length(f_twopoolField)),
                                      rep("three-pool", each = length(f_threepoolField))),
                            pools = c(rep("fast_onepool", length(f_onepoolIncu)),
                                      rep(paste0(c("fast","slow"),"_twopool"), each = length(f_twopoolIncu)/2),
                                      rep(paste0(c("fast","slow","passive"),"_threepool"), each = length(f_threepoolIncu)/3),
                                      rep(paste0(c("fast","slow"),"_twopool"), each = length(f_twopoolField)/2),
                                      rep(paste0(c("fast","slow","passive"),"_threepool"), each = length(f_threepoolField)/3)),
                            experiment = c(rep("incubation", length(f_valuesIncu)),
                                           rep("field_experiment", length(f_valuesField))))
#1.1plot k/f among defferent pool-model
kYticks <- c(expression("10"^-6*""),expression("10"^-4*""),expression("10"^-2*""),1)
fYticks <- c(0,25,50,75,100)
k_incuOnepool <- subset(k_incuVsField, experiment == "incubation" & model == "one-pool")
k_incuTwopool <- subset(k_incuVsField, experiment == "incubation" & model == "two-pool")
k_incuThreepool <- subset(k_incuVsField, experiment == "incubation" & model == "three-pool")
k_incuThreepool$pools <- factor(k_incuThreepool$pools,
                                levels=paste0(c("fast","slow","passive"),"_threepool")) 
f_incuOnepool <- subset(f_incuVsField, experiment == "incubation" & model == "one-pool")
f_incuTwopool <- subset(f_incuVsField, experiment == "incubation" & model == "two-pool")
f_incuThreepool <- subset(f_incuVsField, experiment == "incubation" & model == "three-pool")
f_incuThreepool$pools <- factor(f_incuThreepool$pools,
                                levels=paste0(c("fast","slow","passive"),"_threepool"))
kDataList <- list(k_incuTwopool,k_incuThreepool)
kYName <- c(expression("Reference decomposition rate (d"^-1*")"), NULL)
fDataList <- list(f_incuTwopool,f_incuThreepool)
fYName <- c(expression("Relative pool size (%)"), NULL, NULL)
XName <- c("Two-pool model (M2)","Three-pool model (M3)")
XLabels_K <- list(c("k1","k2"),c("k1","k2","k3"))
XLabels_F <- list(c("f1","f2"),c("f1","f2","f3"))
kYLabels <- list(kYticks, NULL)
fYLabels <- list(fYticks, NULL)
P_kList <- list()
P_fList <- list()
for(i in 1:length(kDataList)){
  P_kList[[i]] <- ggplot(kDataList[[i]], aes(x = pools, y = log10(rateConstant), fill = pools)) +
    geom_violin(trim = FALSE, na.rm = T, scale = "width", color = "white") +
    geom_boxplot(na.rm = T, width = 0.2, color="black" ,size = 1) +
    scale_y_continuous(name = kYName[i], breaks = seq(-6, 0, 2), limits = c(-6,0), labels = kYLabels[[i]])+
    scale_x_discrete(name = XName[i], labels = XLabels_K[[i]]) +
    scale_fill_manual(values = c("#56B4E9","#E69F00","#73BF00")) +
    theme_bw() +
    theme(plot.title = element_text(size = font_size3, face =  "bold"),
          axis.title = element_text(size = font_size3),
          axis.text = element_text(size = font_size2),
          legend.position = "none",
          panel.grid = element_blank(),
          text = element_text(family = font_used_inall))
}
for(i in 1:length(fDataList)){
  P_fList[[i]] <- ggplot(fDataList[[i]], aes(x = pools, y = poolSize, fill = pools)) +
    geom_violin(trim = FALSE, na.rm = T, scale = "width", color = "white") +
    geom_boxplot(na.rm = T, width = 0.2, color="black" ,size = 1) +
    scale_y_continuous(name = fYName[i], breaks = seq(0, 100, 25), limits = c(0,100), labels = fYLabels[[i]])+
    scale_x_discrete(name = XName[i], labels = XLabels_F[[i]]) +
    scale_fill_manual(values = c("#56B4E9","#E69F00","#73BF00")) +
    theme_bw() +
    theme(plot.title = element_text(size = font_size3, face =  "bold"),
          axis.title = element_text(size = font_size3),
          axis.text = element_text(size = font_size2),
          legend.position = "none",
          plot.margin = unit(c(0.2,0.2,0.2,0.45),"cm"), 
          panel.grid = element_blank(),                 
          text = element_text(family = font_used_inall))
}
P_F <- ggarrange(P_fList[[1]],P_fList[[2]],
                 labels = c("d","e"),
                 font.label = list(size = font_size3, family = font_used_inall))
#===============================================================================
#1.2.1plot k among different experiment conditions
kYticks2 <- c(expression("10"^-7*""),expression("10"^-5*""),expression("10"^-3*""),expression("10"^-1*""),1)
k_typicalTwopool <- data.frame(rateConstant = c(inputData$`k1(day-1)`[which(inputData$MODELS == typicalModel[1])],
                                                inputData$`k2(day-1)`[which(inputData$MODELS == typicalModel[1])],
                                                inputData$`k1(day-1)`[which(inputData$MODELS == typicalModel[2])],
                                                inputData$`k2(day-1)`[which(inputData$MODELS == typicalModel[2])]),
                               Model = rep(c(typicalModel[1],typicalModel[2]), each = 2),
                               pools = rep(paste0(c("fast","slow"),"_twopool"),2))
k_typicalTwopool$pools <- factor(k_typicalTwopool$pools,
                                 levels = paste0(c("fast","slow"),"_twopool"))  
k_typicalThreeopool <- data.frame(rateConstant = c(inputData$`k1(day-1)`[which(inputData$MODELS == typicalModel[3])],
                                                   inputData$`k2(day-1)`[which(inputData$MODELS == typicalModel[3])],
                                                   inputData$`k3(day-1)`[which(inputData$MODELS == typicalModel[3])],
                                                   inputData$`k1(day-1)`[which(inputData$MODELS == typicalModel[4])],
                                                   inputData$`k2(day-1)`[which(inputData$MODELS == typicalModel[4])],
                                                   inputData$`k3(day-1)`[which(inputData$MODELS == typicalModel[4])]),
                                  Model = rep(c(typicalModel[3],typicalModel[4]), each = 3),
                                  pools = rep(paste0(c("fast","slow","passive"),"_threepool"),2))
k_typicalThreeopool$pools <- factor(k_typicalThreeopool$pools,
                                    levels = paste0(c("fast","slow","passive"),"_threepool"))
k_typicalThreeopool$Model[which(k_typicalThreeopool$Model == "CLM")] <- "CLMcn"
Pk_incuVSFieldTwopool <- P_kList[[1]] +  
  scale_y_continuous(name = kYName[1], breaks = seq(-7, 1, 2), limits = c(-7,1), labels = kYticks2) +
  geom_point(data = k_typicalTwopool, mapping = aes(x = pools, y = log10(rateConstant), shape = Model), size = 4) +
  scale_shape_manual(values = c(0, 1)) +
  theme(legend.position = c(0.20,0.20),
        legend.title = element_text(size = font_size3),
        legend.text = element_text(size = font_size2),
        plot.margin = unit(c(0.2,0.2,0.2,0.45),"cm")) +            
  guides(fill = F)                                    
Pk_incuVSFieldThreepool <- P_kList[[2]] +  
  scale_y_continuous(name = NULL, breaks = seq(-7, 1, 2), limits = c(-7,1), labels = NULL) +
  geom_point(data = k_typicalThreeopool, mapping = aes(x = pools, y = log10(rateConstant), shape = Model), size = 4) +
  scale_shape_manual(values = c(2, 3)) +
  theme(legend.position = c(0.20,0.20),
        legend.title = element_text(size = font_size3),
        legend.text = element_text(size = font_size2),
        plot.margin = unit(c(0.2,0.2,0.2,0.6),"cm")) +
  guides(fill = F)
P_K <- ggarrange(Pk_incuVSFieldTwopool,Pk_incuVSFieldThreepool,
                 labels = c("b","c"),
                 font.label = list(size = font_size3, family = font_used_inall))
#=====================================================================================================
#===============================================================================
#plot sampling position
xbrks <- seq(-180, 192, 60)
ybrks <- seq(-90, 90, 30)
xlbls <- unlist(lapply(xbrks, function(x) 
  ifelse(x < 0, gsub("[-]","",paste(x, "буW")), ifelse(x > 0, paste(x, "буE"), x))))
ylbls <- unlist(lapply(ybrks, function(x) 
  ifelse(x < 0, gsub("[-]","",paste(x, "буS")), ifelse(x > 0, paste(x, "буN"), x))))
mapworld <- borders("world",colour = "gray50",fill="white")
mp <- ggplot(inputData) + mapworld + 
  geom_point(aes(x = `LON(E)`, y = `LAT(N)`), color="red") +
  scale_y_continuous(name = "Latitude", breaks = ybrks, limits = c(-90,90), expand = c(0,0), labels = ylbls) +
  scale_x_continuous(name = "Longitude", breaks = xbrks, limits = c(-180,192), expand = c(0.01,0.01), labels = xlbls) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = font_size3),
        axis.text = element_text(size = font_size2),
        legend.position = "none",
        text = element_text(family = font_used_inall))
######################################################################
##########################  Figure 1   ###############################
ggarrange(mp,P_K,P_F,ncol=1, nrow=3, labels = c("a",NA,NA),
          font.label = list(size = font_size3, family = font_used_inall))
######################   save as pdf file  ###########################
ggsave("Figure 1.pdf", width = 14, height = 20, unit = "cm", dpi = 600)

#Kruskal-Wallis Test for k1,k2,f1,f2(two-pool vs. three-pool)
twopool_k1 <- k_incuTwopool$rateConstant[which(k_incuTwopool$pools == "fast_twopool")]
twopool_k2 <- k_incuTwopool$rateConstant[which(k_incuTwopool$pools == "slow_twopool")]
twopool_f1 <- f_incuTwopool$poolSize[which(f_incuTwopool$pools == "fast_twopool")]
twopool_f2 <- f_incuTwopool$poolSize[which(f_incuTwopool$pools == "slow_twopool")]
threepool_k1 <- k_incuThreepool$rateConstant[which(k_incuThreepool$pools == "fast_threepool")]
threepool_k2 <- k_incuThreepool$rateConstant[which(k_incuThreepool$pools == "slow_threepool")]
threepool_k3 <- k_incuThreepool$rateConstant[which(k_incuThreepool$pools == "passive_threepool")]
threepool_f1 <- f_incuThreepool$poolSize[which(f_incuThreepool$pools == "fast_threepool")]
threepool_f2 <- f_incuThreepool$poolSize[which(f_incuThreepool$pools == "slow_threepool")]
threepool_f3 <- f_incuThreepool$poolSize[which(f_incuThreepool$pools == "passive_threepool")]
KWdatak1 <- data.frame(x = c(twopool_k1,threepool_k1), g = factor(rep(1:2,c(length(twopool_k1),length(threepool_k1)))))
KWdatak2 <- data.frame(x = c(twopool_k2,threepool_k2), g = factor(rep(1:2,c(length(twopool_k2),length(threepool_k2)))))
KWdataf1 <- data.frame(x = c(twopool_f1,threepool_f1), g = factor(rep(1:2,c(length(twopool_f1),length(threepool_f1)))))
KWdataf2 <- data.frame(x = c(twopool_f2,threepool_f2), g = factor(rep(1:2,c(length(twopool_f2),length(threepool_f2)))))
KWDataList <- list(KWdatak1, KWdatak2, KWdataf1, KWdataf2)
#save KW-Test results(p-values)
KWTestResult <- data.frame(parameter = c("k1","k2","f1","f2"), p_value = c(rep(0,4)))
for(i in 1:length(KWDataList)){
  KWTest <- kruskal.test(x~g, data = KWDataList[[i]])
  KWTestResult$p_value[i] <- KWTest$p.value
}
write.xlsx(KWTestResult,"KWTest.xlsx")
kineticValue <- list(k_onepoolIncu,twopool_k1,twopool_k2,twopool_f1,twopool_f2,
                     threepool_k1,threepool_k2,threepool_k3,threepool_f1,threepool_f2,threepool_f3)
for(i in 1:length(kineticValue)){
  datatmp <- na.omit(kineticValue[[i]])
  STATVALUES[i+6,] <- c(summary(datatmp),quantile(datatmp, probs = c(0.025,0.975)),
                        quantile(datatmp, probs = c(0.05,0.95)),
                        quantile(datatmp, probs = c(0.25,0.75)))
}
rownames(STATVALUES) <- c("R2cum","R2noncum","R2all","RMSEcum","RMSEnoncum","RMSEall","obs_onepool_k",
                          paste0("obs_twopool_",c("k1","k2","f1","f2")),
                          paste0("obs_threepool_",c("k1","k2","k3","f1","f2","f3")))
write.csv(STATVALUES,"STATVALUES.csv")

#===============================================================================
for(dataForm in 0:1){
  saveDifData <- paste0("globalresults_final/",c("initialDataForm","log10DataForm"))
  saveDifPrameter <- c(paste0("globalresults_final/initialDataForm/modelwith_",c(5,6,7),"preds"),
                       paste0("globalresults_final/log10DataForm/modelwith_",c(5,6,7),"preds"))
  dirs <- c(saveDifData, saveDifPrameter)
  for(i in 1:length(dirs)){
    ifelse(!dir.exists(dirs[i]), dir.create(dirs[i]), FALSE)
  }
  #three predictors condition
  predictors1 <- c("MAP","MAT","sand","clay","pH")
  predictors2 <- c("MAP","MAT","sand","clay","temperature","pH")
  predictors3 <- c("MAP","MAT","sand","clay","temperature","pH","elevation")
  predictors_List <- list(predictors1, predictors2, predictors3)
  outcomes <- c(paste0("k",1:3,"(day-1)"),paste0("f",1:3,"(%)"))
  outcome_variables <- c(paste0("k",1:3),paste0("f",1:3))
  varImp_PlotXname <- "Relative importance (%)"
  bestModels <- list()        
  model_list <- list()        
  varImpContratList <- list()
  pdpContrastList <- list()
  n_loop <- length(outcomes) 
  n_loop_ml <- 100              #number of loop for machine learning
  for(j in 1:length(predictors_List)){
    set.seed(2021)           
    i <- 0      
    count <- 0                
    n_picture <- 0          
    while(i < n_loop){
      i <- i + 1
      if(count == 0){
        solutions_InputData <- inputData[twopoolStartPosition[1]:threepoolStartPosition[1],]
        if(i == 3){i <- i + 1}
        if(i == 6){ 
          count  <-  count + 1
          i <- 1
        }
      }
      if(count == 1){
        solutions_InputData <- inputData[threepoolStartPosition[1]:twopoolStartPosition[2],]
      }
      n_picture <- n_picture + 1
      #read predictors and corresponding outcome data
      predData <- solutions_InputData[,match("LAT(N)",colnames(inputData)):match("ELEVATION(m)",colnames(inputData))]
      if(j == 1){predData <- predData[,-c(13,12,11,10,6,8,1,2)]}   
      if(j == 2){predData <- predData[,-c(13,12,11,10,6,1,2)]}     
      if(j == 3){predData <- predData[,-c(12,11,10,6,1,2)]}        
      colnames(predData) <- predictors_List[[j]]
      outDataInit <- solutions_InputData[,outcomes[i]]
      #remove rows containing NA
      noNArows <- which(rowSums(is.na(predData)) == 0 & !is.na(outDataInit))
      predData <- predData[noNArows,]
      outDataInit <- outDataInit[noNArows]
      outData <- outDataInit
      if(dataForm == 1){
        outData <- log10(outData)
      }
      #train/test edition(75% for training,25% for testing)
      inTrain <- createDataPartition(outData, p = 3/4, list = F)
      trainx <- predData[inTrain,]
      testx <- predData[-inTrain,]
      trainy <- outData[inTrain]
      testy <- outData[-inTrain]
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      #1.choose gbm for machine learning
      fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 3,
                                 returnResamp = "all")
      gbmBestFit <- train(predData, outData, 
                          method = "gbm", 
                          trControl = fitControl, 
                          verbose = F, 
                          metric = "Rsquared")
      pred_rf <- predict(gbmBestFit, predData)
      pred_rf <- data.frame(pred = pred_rf, obs = outData)
      R2_gbmBestFit <- R2(pred_rf$obs,pred_rf$pred)
      #train/test edition gbm
      gbmFitBestTrainTest <- train(trainx, trainy, 
                                   method = "gbm", 
                                   trControl = fitControl, 
                                   verbose = F, 
                                   metric = "Rsquared")
      pred_rf <- predict(gbmFitBestTrainTest, testx)
      pred_rf <- data.frame(pred = pred_rf, obs = testy)
      R2_gbmFitBestTrainTest <- R2(pred_rf$obs,pred_rf$pred)
      for(k in 1:n_loop_ml){
        gbmFit <- train(predData, outData, 
                        method = "gbm", 
                        trControl = fitControl, 
                        verbose = F, 
                        metric = "Rsquared")
        pred_rf <- predict(gbmFit, predData)
        pred_rf <- data.frame(pred = pred_rf, obs = outData)
        R2_gbmFit <- R2(pred_rf$obs,pred_rf$pred)
        if(R2_gbmFit > R2_gbmBestFit){
          gbmBestFit <- gbmFit
          R2_gbmBestFit <- R2_gbmFit
        }
        #train/test edition gbm
        gbmFitTrainTest <- train(trainx, trainy, 
                                 method = "gbm", 
                                 trControl = fitControl, 
                                 verbose = F, 
                                 metric = "Rsquared")
        pred_rf <- predict(gbmFitTrainTest, testx)
        pred_rf <- data.frame(pred = pred_rf, obs = testy)
        R2_gbmFitTrainTest <- R2(pred_rf$obs,pred_rf$pred)
        if(R2_gbmFitTrainTest > R2_gbmFitBestTrainTest){
          gbmFitBestTrainTest <- gbmFitTrainTest
          R2_gbmFitBestTrainTest <- R2_gbmFitTrainTest
        }
      }
      colnames(stepDateFrame) <- c("Y", predictors_List[[j]])
      lmFit <- lm(Y ~., stepDateFrame)                 
      lmFitStep <- step(lmFit, trace = F)                
      model_list[[n_picture]] <- list(gbmBestFit,gbmFitBestTrainTest,lmFit,lmFitStep)
	  bestModels[[j]] <- model_list
  }
  save(bestModels, file = "bestModels.Rdata")
}

#set seed to ensure the reproducibility
set.seed(999)
breaks <- list(c(-4,-3,-2,-1,0),c(-5,-4,-3,-2),c(0,20,40,60),c(40,60,80,100),
               c(-3,-2,-1),c(-4,-3,-2,-1),c(-6,-5,-4,-3),
               c(0,2.5,5,7.5,10),c(0,25,50,75))
ticks <- list(c(expression("10"^-4),expression("10"^-3),expression("10"^-2),expression("10"^-1),1),
              c(expression("10"^-5),expression("10"^-4),expression("10"^-3),expression("10"^-2)),
              c(0,20,40,60),
              c(40,60,80,100),
              c(expression("10"^-3),expression("10"^-2),expression("10"^-1)),
              c(expression("10"^-4),expression("10"^-3),expression("10"^-2),expression("10"^-1)),
              c(expression("10"^-6),expression("10"^-5),expression("10"^-4),expression("10"^-3)),
              c(0,2.5,5,7.5,10),
              c(0,25,50,75))
predVSobs_PlotXname <- c(expression("M2-k1 observed (d"^-1*")"),expression("M2-k2 observed (d"^-1*")"),
                         "M2-f1 observed (%)","M2-f2 observed (%)",
                         expression("M3-k1 observed (d"^-1*")"),expression("M3-k2 observed (d"^-1*")"),
                         expression("M3-k3 observed (d"^-1*")"),
                         "M3-f1 observed (%)","M3-f2 observed (%)","M3-f3 observed (%)")
predVSobs_PlotYname <- c(expression("M2-k1 predicted (d"^-1*")"),expression("M2-k2 predicted (d"^-1*")"),
                         "M2-f1 predicted (%)","M2-f2 predicted (%)",
                         expression("M3-k1 predicted (d"^-1*")"),expression("M3-k2 predicted (d"^-1*")"),
                         expression("M3-k3 predicted (d"^-1*")"),
                         "M3-f1 predicted (%)","M3-f2 predicted (%)","M3-f3 predicted (%)")
varImp_PlotXname <- "Relative importance (%)"
varImp_PlotYname <- c(paste0("M2-k",1:2),paste0("M2-f",1:2),
                      paste0("M3-k",1:3),paste0("M3-f",1:3))
gbmVsLmPlotName <- c(paste0("gbmVsLm_two_k",1:2,".png"), paste0("gbmVsLm_two_f",1:2,".png"),
                     paste0("gbmVsLm_three_k",1:3,".png"), paste0("gbmVsLm_three_f",1:3,".png"))
obsVsPredResults <- c(paste0("obsVsPred_M2-k",1:2),paste0("obsVsPred_M2-f",1:2),
                      paste0("obsVsPred_M3-k",1:3),paste0("obsVsPred_M3-f",1:3))
obsVsPredResults <- paste0(obsVsPredResults,".xlsx")
for(k in 1:2){    #the optimal model is k=1 and j=3 (initial dataform and 7preds)
  modelDir <- paste0(getwd(),"/",saveDifData[k])
  INDEX <- read.csv(paste0(modelDir,"/INDEX_2.csv")
  obsVsPredValues <- list()
  for(j in 1:3){
    bestModels_str <- load(paste0(modelDir,"/bestModels.Rdata"))
    bestModels <- eval(parse(text = bestModels_str))
    model_multiPred <- bestModels[[j]]
    colors <- c("#156077","#f46f20")
    P_gbmVsLmList <- list()
    P_gbmvarImpList <- list()
    P_gbmpdpList <- list()
    obsVsPredValuesTmp <- data.frame()
    for(i in 1:(length(model_multiPred)-1)){
      model_gbm <- model_multiPred[[i]][[1]]
      pred_gbm <- predict(model_gbm, model_gbm$trainingData[,-length(model_gbm$trainingData)])
      pred_rfgbm <- data.frame(pred = pred_gbm, obs = model_gbm$trainingData[,length(model_gbm$trainingData)])
      pred_rfgbm$pred[pred_rfgbm$pred < summary(pred_rfgbm$obs)[1]] <- summary(pred_rfgbm$obs)[1]
      model_lm <- model_multiPred[[i]][[3]]
      pred_lm <- predict(model_lm, model_lm$model[,-1])
      pred_rflm <- data.frame(pred = pred_lm, obs = model_lm$model[,1])
      pred_rflm$pred[pred_rflm$pred < summary(pred_rflm$obs)[1]] <- summary(pred_rflm$obs)[1]
      obsVsPredValuesTmp <- data.frame(obs = pred_rfgbm$obs,
                                         pred_gbm = pred_rfgbm$pred,
                                         pred_lm = pred_rflm$pred)
      if(j == 1){
        obsVsPredValues[[i]] <- obsVsPredValuesTmp
      }
      else{
        obsVsPredValues[[i]] <- cbind(obsVsPredValues[[i]],obsVsPredValuesTmp)
      }
      if(j == 3){
        write.xlsx(obsVsPredValues[[i]],obsVsPredResults[[i]])
        file.copy(obsVsPredResults[[i]], modelDir, overwrite = T)
        file.remove(obsVsPredResults[[i]])
      }
      gbmVsLmData <- data.frame(pred = c(pred_rfgbm$pred,pred_rflm$pred),
                                obs = c(pred_rfgbm$obs,pred_rflm$obs),
                                model = c(rep("gbm",nrow(pred_rfgbm)),rep("lm",nrow(pred_rflm))))
      #change k values to log10 form for plotting
      if(k == 1 & i %in% c(1,2,5,6,7)){
        gbmVsLmData$pred <- log10(gbmVsLmData$pred)
        gbmVsLmData$obs <- log10(gbmVsLmData$obs)
        pred_rfgbm$pred <- log10(pred_rfgbm$pred)
        pred_rfgbm$obs <- log10(pred_rfgbm$obs)
      }
      gbmVsLmData$model <- factor(gbmVsLmData$model, levels = c("gbm","lm"))
      P_gbmVsLmList[[i]] <- ggplot(gbmVsLmData, aes(y = pred, x = obs, color = model, shape = model)) +
        geom_point() +
        scale_color_manual(values = colors) +
        geom_abline(intercept = 0, slope = 1, size = 0.5, color = "grey") +     # 1:1-Line
        scale_shape_manual(values = c(16,17)) +
        scale_x_continuous(name = predVSobs_PlotXname[i], 
                           limits = c(min(pred_rfgbm$obs),
                                      max(pred_rfgbm$obs)+(max(pred_rfgbm$obs)-min(pred_rfgbm$obs))/5),
                           breaks = breaks[[i]], labels = ticks[[i]]) +
        scale_y_continuous(name = predVSobs_PlotYname[i], 
                           limits = c(min(pred_rfgbm$obs),
                                      max(pred_rfgbm$obs)+(max(pred_rfgbm$obs)-min(pred_rfgbm$obs))/5),
                           breaks = breaks[[i]], labels = ticks[[i]]) +
        theme_bw() +
        annotate("text", label = bquote(R ^ 2 * " = " * .(round(INDEX[5*(i-1)+1,6*(j-1)+5],2))),
                 x = min(pred_rfgbm$obs),
                 y = max(pred_rfgbm$obs) + (max(pred_rfgbm$obs)-min(pred_rfgbm$obs))/6,
                 hjust = 0, size = anno_size, family = font_used_inall, color = colors[1]) +
        annotate("text", label = paste("RMSEn=",round(INDEX[5*(i-1)+1,6*(j-1)+6],2)),
                 x =  min(pred_rfgbm$obs),
                 y =  max(pred_rfgbm$obs)+(max(pred_rfgbm$obs)-min(pred_rfgbm$obs))/15,
                 hjust = 0, size = anno_size, family = font_used_inall, color = colors[1]) +
        annotate("text", label = bquote(R ^ 2 * " = " * .(round(INDEX[5*(i-1)+4,6*(j-1)+5],2))),
                 x =  min(pred_rfgbm$obs)+(max(pred_rfgbm$obs)-min(pred_rfgbm$obs))/1.9,
                 y =  max(pred_rfgbm$obs) + (max(pred_rfgbm$obs)-min(pred_rfgbm$obs))/6,
                 hjust = 0, size = anno_size, family = font_used_inall, color = colors[2]) +
        annotate("text", label = paste("RMSEn=",round(INDEX[5*(i-1)+4,6*(j-1)+6],2)),
                 x =  min(pred_rfgbm$obs)+(max(pred_rfgbm$obs)-min(pred_rfgbm$obs))/1.9,
                 y =  max(pred_rfgbm$obs)+(max(pred_rfgbm$obs)-min(pred_rfgbm$obs))/15,
                 hjust = 0, size = anno_size, family = font_used_inall, color = colors[2]) +
        theme(axis.text = element_text(size = font_size2),
              axis.title = element_text(size = font_size3),
              panel.grid = element_blank(),
              legend.direction = "horizontal",
              legend.background = element_rect(color = "black"),
              legend.key = element_blank(),
              legend.position = c(0.2,0.94),
              legend.text = element_text(size = font_size2),
              legend.title = element_text(size = font_size3),
              text = element_text(family = font_used_inall))
      P_gbmVsLmList[[i]]
      #relative importance
      mlvariable_importance <- varImp(model_gbm)
      mlvarImpDataFull <- data.frame(imp = mlvariable_importance$importance$Overall,
                                     var = rownames(mlvariable_importance$importance))
      mlvarImpDataFull$imp <- mlvarImpDataFull$imp/sum(mlvarImpDataFull$imp)*100   #normalization
      mlvarImpDataFull$imp <- round(mlvarImpDataFull$imp,1)   
      if(j == 1){
        mlvarImpDataFull$imp[5] <- 100 - sum(mlvarImpDataFull$imp[-5])  #insure the sum equal to 100
        mlvarImpDataFull$imp[5] <- round(mlvarImpDataFull$imp[5],1)
        mlvarImpDataFull$var <- c("MAP","MAT","Sand","Clay","pH")
      }
      if(j == 2){
        mlvarImpDataFull$imp[6] <- 100 - sum(mlvarImpDataFull$imp[-6])  #insure the sum equal to 100
        mlvarImpDataFull$imp[6] <- round(mlvarImpDataFull$imp[6],1)
        mlvarImpDataFull$var <- c("MAP","MAT","Sand","Clay","IncT","pH")
      }
      if(j == 3){
        mlvarImpDataFull$imp[7] <- 100 - sum(mlvarImpDataFull$imp[-7])  #insure the sum equal to 100
        mlvarImpDataFull$imp[7] <- round(mlvarImpDataFull$imp[7],1)
        mlvarImpDataFull$var <- c("MAP","MAT","Sand","Clay","IncT","pH","Elev")
      }
      P_gbmvarImpList[[i]] <- ggplot(mlvarImpDataFull, aes(imp, reorder(var, imp), fill = var)) + 
        geom_bar(stat = "identity") +
        geom_text(aes(x = imp + 3, y = var, label = imp), 
                  family = font_used_inall, size = anno_size) +
        theme_bw() +
        scale_x_continuous(limit = c(0, max(mlvarImpDataFull$imp) + 8), expand = c(0.0, 0.01)) + 
        labs(x = varImp_PlotXname, y = varImp_PlotYname[i]) +
        theme(axis.text = element_text(size = font_size2),
              axis.title = element_text(size = font_size3),
              legend.position = "none",
              panel.grid = element_blank(),
              text = element_text(family = font_used_inall))
      }
    }
    #combination chart of gbm VS. lm
######################################################################
##########################  Figure 2   ###############################
    ggarrange(P_gbmVsLmList[[1]],P_gbmVsLmList[[2]],P_gbmVsLmList[[3]],
              P_gbmVsLmList[[5]],P_gbmVsLmList[[6]],P_gbmVsLmList[[7]],
              P_gbmVsLmList[[8]],P_gbmVsLmList[[9]],P_legend_GbmVsLm,
              ncol=3, nrow=3, labels = c(letters[1:8]), 
              font.label = list(size = font_size3, family = font_used_inall),
              common.legend = TRUE, legend = "none")
    ######################   save as pdf file ########################
    ggsave("Figure 2.pdf", width = 17, height = 18, unit = "cm", dpi = 600)
#combination chart of variable relative importance
######################################################################
##########################  Figure 3   ###############################
    ggarrange(P_gbmvarImpList[[1]],P_gbmvarImpList[[2]],P_labels,
              P_gbmvarImpList[[5]],P_gbmvarImpList[[6]],P_gbmvarImpList[[7]],
              ncol=3, nrow=2, labels = c("a","b",NA,"c","d","e"), 
              font.label = list(size = font_size3, family = font_used_inall),
              common.legend = TRUE, legend = "none")
    ######################   save as pdf file ########################
    ggsave("Figure 3.pdf", width = 17, height = 18, unit = "cm", dpi = 600)
  }
}

#Read the monthly average precipitation and sum up to get MAP
globalDataset <- paste0(getwd(),"/globalDatasets")
MMP_dir <- dir(paste0(globalDataset,"/wc2.1_30s_prec"), full.names = T)
MMP <- stack(MMP_dir)
MAP <- sum(MMP, na.rm = T)
MMT_dir <- dir(paste0(globalDataset,"/wc2.1_30s_tavg"), full.names = T)
MMT <- stack(MMT_dir)                                 
MAT <- mean(MMT, na.rm = T)
rm(MMP, MMT)
SAND <- raster(paste0(globalDataset,"/SOIL_1km/SAND_CONTENT/SNDPPT_M_sl1_1km_ll.tif"))
CLAY <- raster(paste0(globalDataset,"/SOIL_1km/CLAY_CONTENT/CLYPPT_M_sl1_1km_ll.tif"))
SILT <- raster(paste0(globalDataset,"/SOIL_1km/SILT_CONTENT/SLTPPT_M_sl1_1km_ll.tif"))
PH <- raster(paste0(globalDataset,"/SOIL_1km/PH_H2O/PHIHOX_M_sl1_1km_ll.tif"))
PH <- PH/10
#read DEM data
ELEVATION <- raster(paste0(globalDataset,"/DEM/globalDEM.tif"))
AOI <- SAND@extent                                      #Get the spatial range of SAND data
MAP_clip <- crop(MAP, AOI)
MAT_clip <- crop(MAT, AOI)
ELEVATION_clip <- crop(ELEVATION, AOI)
TEMPERATURE <- CLAY
TEMPERATURE[TEMPERATURE >= 0] <- 20
pred1 <- stack(MAP_clip, MAT_clip, SAND, CLAY, PH)
pred2 <- stack(MAP_clip, MAT_clip, SAND, CLAY, TEMPERATURE, PH)
pred3 <- stack(MAP_clip, MAT_clip, SAND, CLAY, TEMPERATURE, PH, ELEVATION_clip)
names(pred1) <- predictors1
names(pred2) <- predictors2
names(pred3) <- predictors3
global_pred <- list(pred1, pred2, pred3)
P_mean_predVarWithY <- list()
P_predVar <- list()
rm(MAP, MAT, SAND, CLAY, PH, MAP_clip, MAT_clip)
for(k in 1:2){
  modelDir <- paste0(getwd(),"/",saveDifData[k])
  for(j in 1:3){
    bestModels_str <- load(paste0(modelDir,"/bestModels.Rdata")) 
    bestModels <- eval(parse(text = bestModels_str))
    model_multiPred <- bestModels[[j]]
    globalOutcomes <- list()
    for(i in 5:9){
      bestModel <- bestModels[[j]][[i]][[1]]
      globalOutcomes[[i-4]] <- raster::predict(object = global_pred[[j]], model = bestModel)
    }
    save(globalOutcomes, file = "global_predictios.Rdata")                
    file.copy("global_predictios.Rdata", saveDifPrameter[3*(k-1)+j], overwrite = T)
    file.remove("global_predictios.Rdata")
    tifResultNames <- c(paste0("k",1:3,"_1km.tif"), paste0("f",1:2,"_1km.tif"))
    for(i in 1:5){
      writeRaster(globalOutcomes[[i]], filename = tifResultNames[i], datatype='FLT4S')
      file.copy(tifResultNames[i], saveDifPrameter[3*(k-1)+j], overwrite = T)
      file.remove(tifResultNames[i])
    }
  }
}

#append the summary results of global predictions to the obs summary results
STATVALUES <- read.csv("STATVALUES.csv")
P_medianPredWithY <- list()
P_red <- list()
tifResultNames <- c(paste0("k",1:3,"_1km.tif"), paste0("f",1:2,"_1km.tif"))
legendNames <- c(expression("M3-k1 (d"^-1*")  "),expression("M3-k2 (d"^-1*")  "),expression("M3-k3 (d"^-1*")  "),
                 paste0("M3-f",1:3," (%)  "))
legendBreaks <- list(c(6e-4,0.1,0.2,0.3),
                     c(9e-5,0.01,0.02,0.035),
                     c(5e-7,5e-5,1e-4,1.3e-4),
                     c(0.03,2,4,6.4),
                     c(0.17,20,40,60),
                     c(30,50,70,90,100))
legendLables <- list(c("6.0e-4","1.0e-1","2.0e-1","3.0e-1"),
                     c("9.0e-5","1.0e-2","2.0e-2","3.5e-2"),
                     c("5.0e-7","5.0e-5","1.0e-4","1.3e-4"),
                     c("0.03","2","4","6.4"),
                     c("0.17","20","40","67.9"),
                     c("30","50","70","90","100"))
xbrks <-  c(seq(-180, 190, 60))
ybrks <-  seq(-60, 90, 30)
xlbls <- unlist(lapply(xbrks, function(x) 
  ifelse(x < 0, gsub("[-]","",paste(x, "буW")), ifelse(x > 0, paste(x, "буE"), x))))
ylbls <- unlist(lapply(ybrks, function(x) 
  ifelse(x < 0, gsub("[-]","",paste(x, "буS")), ifelse(x > 0, paste(x, "буN"), x))))
outPicNames <- c("k_combined.png","f_combined.png")
for(k in 1:1){  #only plot the best model result(initial dataform, 7preds)
  modelDir <- paste0(getwd(),"/",saveDifData[k])
  SUMVALUES <- read.csv(paste0(modelDir,"/SUMVALUES.csv"))
  obsMin <- SUMVALUES$min_obs_s3[c(17,21,25,29,33,37)] 
  obsMax <- SUMVALUES$max_obs_s3[c(17,21,25,29,33,37)]  
  for(j in 3:3){
    for(i in 1:6){ #Figure 4(combine k1,k2,k3), Figure S3(combined f1,f2,f3)
      #Read the output tif file with a resolution of (1/120)degree
      if(i < 6){
        predNames <- paste0(saveDifPrameter[3*(k-1)+j], "/", tifResultNames[i])
        pred <- raster(predNames)
        }
      if(i == 6){  #f3 = 100-f1-f2
        predf1 <- raster(paste0(saveDifPrameter[3*(k-1)+j], "/", tifResultNames[4]))
        predf2 <- raster(paste0(saveDifPrameter[3*(k-1)+j], "/", tifResultNames[5]))
        pred <- 100-(predf1+predf2)
      }
      #Resolution adjusted to 0.5degree by means of aggregation and transferred to data frame
      pred <- aggregate(pred, fact = 60)
      predDataFrame <- rasterToPoints(pred) %>% as.data.frame()
      colnames(predDataFrame) <- c("x","y","value")
      predDataFrame$value[predDataFrame$value < obsMin[i]] <- obsMin[i]
      predDataFrame$value[predDataFrame$value > obsMax[i]] <- obsMax[i]
      #set legend breaks in graph of parameter global distribution
      min_globalPred <- min(predDataFrame$value)
      max_globalPred <- max(predDataFrame$value)
      STATVALUES[i+17,-1] <- c(summary(predDataFrame$value),
                               quantile(predDataFrame$value, probs = c(0.025,0.975)),
                               quantile(predDataFrame$value, probs = c(0.05,0.95)),
                               quantile(predDataFrame$value, probs = c(0.25,0.75)))
      #reverse Log form data to original value
      if(k == 2){
        predDataFrame$value <- 10^predDataFrame$value 
      }
      #The median values with 50%CI of the predicted results changed with latitude were counted and plotted
      predDataFrame %>% group_by(y) %>% summarise(median = median(value),
                                                  CI_low = quantile(value,0.25),
                                                  CI_high = quantile(value,0.75)) -> median_predDataFrame
      max_median <- max(median_predDataFrame$median)
      max_median_sig <- signif(max_median,1)      
      xbrk_medianPredWithY <- c(0,max_median_sig/2,max_median_sig)  
      xlm_medianPredWithY <- c(0,max(median_predDataFrame$CI_high))
      P_medianPredWithY[[i]] <- ggplot(median_predDataFrame) + 
        geom_ribbon(aes(y = y, xmin = CI_low, xmax = CI_high), 
                    orientation = "y", fill = "grey") +
        geom_line(aes(x = median, y = y), orientation = "y", size = 0.5) +      
        scale_x_continuous(name = legendNames[i], breaks = xbrk_medianPredWithY, 
                           limits = xlm_medianPredWithY, expand = c(0,max_median/7)) +
        scale_y_continuous(name = NULL, breaks = ybrks, limits = c(-60, 90), labels = ylbls) +
        theme_bw() + 
        theme(axis.text = element_text(size = font_size2),
              axis.title = element_text(size = font_size3),
              panel.grid = element_blank(),
              plot.margin = unit(c(0.2,0.2,0.2,0.7),"cm"),         
              text = element_text(family = font_used_inall))
      #Mapping the world and overlaying the dependent variable prediction results
      WorldData <- ggplot2::map_data("world")
      WorldData %>% filter(region != "Antarctica") -> WorldData   
      WorldData <- fortify(WorldData)
      P_red[[i]] <- ggplot(predDataFrame) +
        geom_map(data = WorldData, map = WorldData,
                 aes(x = long, y = lat, group = group, map_id = region),
                 fill = "grey", colour = "grey", size = 0.5) +
        geom_tile(aes(x = x, y = y, fill = value)) +
        scale_fill_distiller(palette = "Spectral", breaks = legendBreaks[[i]], 
                             labels = legendLables[[i]]) +  
        labs(fill = legendNames[i]) + xlab("Longitude") + ylab("Latitude") +
        theme_bw() +
        scale_x_continuous(breaks = xbrks, limits = c(-180, 194), labels = xlbls, expand = c(0.01, 0.01)) +
        scale_y_continuous(breaks = ybrks, limits = c(-60, 90), labels = ylbls) +
        theme(legend.position = "right",
              legend.margin = margin(0.1,0.1,0.1,0.02,"cm"),       
              legend.key.width = unit(0.2,"cm"),                    
              legend.key.height = unit(0.8,"cm"),                  
              legend.text = element_text(size = font_size2),        
              legend.title = element_blank(),             
              legend.background = element_blank(),
              axis.text = element_text(size = font_size2),
              axis.title = element_text(size = font_size3),
              panel.grid = element_blank(),
              text = element_text(family = font_used_inall)        
              )
    }
######################################################################
##########################  Figure 4   ###############################
    ggarrange(ggarrange(P_red[[1]],P_medianPredWithY[[1]],widths = c(0.7,0.3),
                        labels = c("a","b"),                       
                        font.label = list(size = font_size3, family = font_used_inall)),  
              ggarrange(P_red[[2]],P_medianPredWithY[[2]],widths = c(0.7,0.3),
                        labels = c("c","d"),
                        font.label = list(size = font_size3, family = font_used_inall)),
              ggarrange(P_red[[3]],P_medianPredWithY[[3]],widths = c(0.7,0.3),
                        labels = c("e","f"),
                        font.label = list(size = font_size3, family = font_used_inall)),
              ncol = 1, nrow = 3)
    ######################   save as pdf file ########################
    ggsave("Figure 4.pdf", width = 17, height = 18, unit = "cm", dpi = 600) 
    }
}

#save the summary results of global predictions
STATVALUES[c(18:23),1] <- c(paste0("globalpred_threepool_",c("k1","k2","k3","f1","f2","f3")))
write.csv(STATVALUES,"STATVALUES.csv")
