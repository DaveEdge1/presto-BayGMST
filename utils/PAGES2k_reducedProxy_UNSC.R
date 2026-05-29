# Written by Wan He, Dec 2016
# modified by Julien Emile-Geay, Feb 21 2017
# - eliminated "eigenvalue" selection bc results were garbage
# - changed PCA to PCR
# - only fill in composite up from tMin onwards
rm(list=ls())
library(dplyr)
library(fda)
library(ggplot2)
library(reshape2)
library(cowplot)
library(grid)
library(glmnet)
library(parallel)
library(matrixStats)
library(dr)
#library(edrGraphicalTools)
library(stringr)
library(spls)
library(superpc)
library(tidyr)
# library(ggmap)   # PReSto template patch: unused — was bloating Docker image with geo deps
# library(maptools)
# library(maps)    # PReSto template patch: unused

library(config)
library(yaml)
library(here)

cfg <- read_yaml(here("config.yml"))

# TO-DOS/ISSUES:
#   1. Breaks down if t1 != 1, make dynamic.
#   2. Fix grabbing desired proxies by type

set.seed(2018)
# specify options to treat the data
t1 <- cfg$partition_years$t1
t2 <- cfg$partition_years$t2
t3 <- 2000#cfg$partition_years$t3
tStart = 1 #define start year (remember: the Common Era does not have a year 0).
tEnd   = 2000 #define end year for the analysis
tce  = tStart:tEnd
nce = length(tce)

#load the temperature data
temp <- read.csv(here(cfg$folder_paths$instr_temp_path))
colnames(temp) <- c("year","T","l95","u95")
print(temp)

#load the proxy data
proxydata =read.csv(here('data','PAGES2K_proxy_matrix_screened_1900-2000.csv'), header=TRUE, sep=",")
#View(proxydata)

#load metadata
metadataproxy =read.csv(here('data','PAGES2K_proxy_metadata_screened_1900-2000.csv'), header=TRUE, sep=",")
#View(metadataproxy)

### GRAB THE DESIRED PROXY TYPE FOR RP
# ptype may be a single string or a list (from PReSto multi-select).
# yaml::read_yaml returns a plain character for a scalar and a list for a
# YAML sequence; unlist() normalises both to a character vector.
ptype_sel <- as.character(unlist(cfg$ptype))
if (length(ptype_sel) == 1 && ptype_sel == 'ALL') {
  NULL # nothing to do — keep every archive type
} else {
  desired_proxies <- metadataproxy[sub("\\..*$", "", metadataproxy$ptype) %in% ptype_sel, ]
  # Fail loudly rather than silently building an empty proxy matrix: an empty
  # selection means the requested ptype(s) don't exist in the metadata (e.g. a
  # LiPD archiveType that wasn't normalised, or a typo in user_config.yml).
  if (nrow(desired_proxies) == 0) {
    stop(sprintf(
      "No proxies match ptype selection {%s}. Available archive types: {%s}.",
      paste(ptype_sel, collapse = ", "),
      paste(sort(unique(sub("\\..*$", "", metadataproxy$ptype))), collapse = ", ")
    ))
  }
  proxydata <- proxydata[, colnames(proxydata) %in% c('year', desired_proxies$X), drop = FALSE]
}

#View(proxydata)
############################

t=proxydata[,1]
proxy=as.data.frame(proxydata[,-1])[which(t%in%tce),]
np=dim(proxy)[2]
#Get the years of the first available data for each proxy
yearMin= matrix(1, 1, np) 
for (i in 1:np) {
  yearMin[i]=min(which(is.finite(proxy[,i])))
}
#ny=dim(proxy)[1]
ny=tEnd-min(yearMin)+1
tStart <- min(yearMin)

#####HISTOGRAM FIRST YEAR#############
datahist <- data.frame(Year = c(yearMin))
histogram_first <- ggplot(data = datahist)+
  geom_histogram(mapping = aes(x = Year),bins = 10,fill='white',color='black')+
  ylab('')+
  theme_bw()+
  theme(axis.title.x = element_text(size=18),axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))
histogram_first


# define analysis options
RP_style = cfg$rp_method  #possible choices: "PCR", "LASSO", "SIR", "SPLS" #### EDIT THIS YML!!!
print(paste0('Method to compute RP: ', RP_style))
## DEFINE FUNCTIONS USED LATER
#define function %!in% which is same as ~ismember
'%!in%' <- function(x,y)!('%in%'(y,x))

#define rfind fun ction as similar to find in matlab
rfind <- function(x)seq(along=x)[as.logical(x)] 

# scale the proxy records to unit variance and zero mean
proxy_scaled = scale(proxy)
#View(proxy_scaled)

# NEW
proxy_filled <- apply(proxy_scaled, 2, function(x) {
  seen_first <- cumsum(!is.na(x)) > 0
  x[is.na(x) & seen_first] <- 0
  x
})
rownames(proxy_filled) <- rownames(proxy_scaled)
colnames(proxy_filled) <- colnames(proxy_scaled)
#View(proxy_filled)

# count the number of available datapoints in each record
navl = matrix(1, 1, np) 
for (i in 1:np){  
  navl[i] = sum(is.finite(proxy_filled[,i]))
}
#View(navl)

##  MAKE A TIME-DEPENDENT COMPOSITE
chunk = 250 # define  segment length
ns = ny/chunk  # number of segments
calib = c(t2:t3)#Define calibration year
RP = matrix(NA,ny,ns)#Initialise the reduced proxies
nprox= matrix(0,1,ns)#Initialise nprox that is later used to store number of proxies used in each segmentation
temp_lastc <- temp[which(temp[,1]%in%calib),2] #last century temperature

#evalues=matrix(1,1,npcs)#Initialise a matrix to store the eigenvalues of PCs
coefflist <- list()
#npcs <- c(11,13,9,11,7,12,10,6) #first PC whose adjusted R2 > 70% THIS WAS HARDCODED LIKE THIS ORIGINALLY
target_adj_r2 <- 0.70
limitup <- c(1,2,2,2,2,1,1,2) ##Marginal Dimension Tests (SIR)
nslicesv <- c(5,7)
proport <- c(1,1,1,1,1,0.35,0.6,0.75)
pZerov <- c(30,30,30,30,30,30,30,30)
#tsholdv <- c(0.53,1.16,0.36,2.62,1.16)
tsholdv <- NULL
npcs2 <- c(3,3,3,3,3,3,3,3)



for (k in 1:ns){   # loop over segments
  show(k)
  tMin = tStart+(k-1)*chunk
  segProxies = rfind(yearMin <=tMin) 
  timeSpan  = rfind(tce>=tMin)
  nt = length(timeSpan)
  one=matrix(1,nt,1)#Initialise a ny*1 one matrix that can be later appended to PC matrix
  proxy_finite= proxy_filled[timeSpan,segProxies]
  cleanNANs <- apply(X = proxy_finite,MARGIN = 2,function(x) sum(is.na(x))/dim(proxy_finite)[1])<0.05
  proxy_finite <- proxy_finite[,cleanNANs]
  time_calib <- rfind(timeSpan %in% calib)
  proxy_finite_calib= proxy_finite[time_calib,]
  cleanNANs_calib <- apply(X = proxy_finite_calib,MARGIN = 2,function(x) sum(is.na(x))/dim(proxy_finite_calib)[1])<0.05
  proxy_finite <- proxy_finite[,cleanNANs_calib]
  proxy_finite[!is.finite(proxy_finite)]=0
  nprox[k] = dim(proxy_finite)[2]
  
  if (RP_style == "PCR") {
    pca <- prcomp(proxy_finite, center = FALSE, scale. = FALSE)
    calib_idx <- rfind(timeSpan %in% calib)
    if (length(temp_lastc) != length(calib_idx)) {
      stop("Calibration temperature and proxy calibration rows do not have the same length.")
    }
    npcs_max <- min(
      ncol(pca$x),
      length(temp_lastc) - 2
    )
    adj_r2 <- sapply(1:npcs_max, function(j) {
      pc_calib <- pca$x[calib_idx, 1:j, drop = FALSE]
      summary(lm(temp_lastc ~ pc_calib))$adj.r.squared
    })
    hits <- which(adj_r2 >= target_adj_r2)
    if (length(hits) > 0) {
      npcs_help <- hits[1]
    } else {
      npcs_help <- which.max(adj_r2)
    }
    show(paste("Segment", k, "using", npcs_help, "PCs; adj R2 =", round(adj_r2[npcs_help], 3)))
    pc_calib <- pca$x[calib_idx, 1:npcs_help, drop = FALSE]
    pc_all <- pca$x[, 1:npcs_help, drop = FALSE]
    a <- lm(temp_lastc ~ pc_calib)
    coeff <- as.numeric(a$coefficients)
    RP[timeSpan, k] <- cbind(one, pc_all) %*% coeff
  }
  if (RP_style=="SPCR") {
    datapc <- list(x=t(proxy_finite[rfind(timeSpan %in% calib),]),y=temp_lastc)
    trainobj <- superpc.train(datapc,type='regression')
    cvobj <- superpc.cv(trainobj,datapc)
    superpc.plotcv(cvobj)
    tsholdv[k] <- cvobj$thresholds[which.max(cvobj$scor[1,])]
    datatest <- list(x=t(proxy_finite))
    fitobj <- superpc.predict(trainobj, datapc, datatest, threshold=tsholdv[k], 
                              n.components=npcs2[k], prediction.type="continuous")
    datapred <- data.frame(temp_lastc,fitobj$v.pred[rfind(timeSpan %in% calib),])
    modelopred <- lm(temp_lastc~.,data = datapred)
    RP[timeSpan,k]   <- predict(modelopred,newdata = data.frame(fitobj$v.pred))
  }
  
  if (RP_style == "LASSO") {
    # Performs a lasso regularization technique to avoid the use of PCA in multiple regression
    Xmatrix <- proxy_finite[rfind(timeSpan %in% calib), ]
    cv.lasso <- cv.glmnet(x = Xmatrix, y = temp_lastc, alpha = 1, parallel = TRUE)
    coeff <- coef(cv.lasso, s = "lambda.min")
    coefflist[[k]] <- coeff
    RP[timeSpan, k] <- predict(cv.lasso, newx = proxy_finite, s = "lambda.min")
  }
  if (RP_style=="SIR"){
    ### TYLER: NEED TO FIGURE OUT HOW TO INSTALL edrGraphicalTools!!
    #Performs a Sliced Inverse Regression (SIR)
    Xmatrix <- proxy_finite[rfind(timeSpan %in% calib),]
    show(dim(Xmatrix)[2])
    if(k >= 6){
      selectModelCSS <- edrSelec(Y = temp_lastc,X = Xmatrix,H = 8,K = 1,
                               method = 'CSS',pZero = nprox[5],NZero = 10000,zeta = 0.1)
      show(max(selectModelCSS$vectSqCor))
      #indexselect <- selectModelCSS$scoreVar>quantile(selectModelCSS$scoreVar,proport[k])
      indexselect <- sort(selectModelCSS$scoreVar,decreasing = T,index.return=T)$ix<=nprox[5]
      show(sum(indexselect))
      Xmatrix <- Xmatrix[,indexselect]
    }
    SIRdr <- dr(temp_lastc ~ Xmatrix,nslices=8)
    namesafterSIR <- names(SIRdr$evectors[,1])%>%str_replace('Xmatrix','')
    Xmatrix <- Xmatrix[,namesafterSIR]
    xtemp <- as.matrix(Xmatrix %*% SIRdr$evectors[,1:(limitup[k])]) 
    modeltemp <- lm(temp_lastc~xtemp)
    xtemptot <- as.matrix(proxy_finite[,namesafterSIR] %*% SIRdr$evectors[,1:(limitup[k])])
    colnames(xtemptot) <- names(modeltemp$coefficients)[(1:(limitup[k]))+1]
    RP[timeSpan,k] <- cbind(one,xtemptot)%*%modeltemp$coefficients
  }
    
  if (RP_style=="SPLS"){
    Xmatrix <- proxy_finite[rfind(timeSpan %in% calib),]
    show(dim(Xmatrix)[2])
    Kint <- seq(1,min(dim(Xmatrix)[2],floor(0.9*dim(Xmatrix)[1]-1)))
    #Kint <- seq(1,70)
    colVars(Xmatrix)
    cvspls <- cv.spls(x = Xmatrix,y = temp_lastc,K = Kint,eta = seq(0.1,0.9,0.1))
    SPLSmodel <- spls(x = Xmatrix,y = temp_lastc,K = cvspls$K.opt,eta = cvspls$eta.opt)
    RP[timeSpan,k] <- predict.spls(SPLSmodel,newx = proxy_finite)
  }
}

colnames(RP) <- paste0('RP',seq(1,8))
RP <- data.frame(Year=1:2000,RP)

RPn <- RP %>% gather(key = 'RPNumber',value = 'Value',RP1:RP8)

plotcombined <- ggplot(data = RPn,mapping = aes(x = Year,y = Value))+
  geom_line(mapping = aes(color=RPNumber))+
  scale_color_brewer(palette = 4,type = 'div')+
  ggtitle(RP_style)+
  theme_bw()
plotcombined

#save(RPn,file=paste0('Plot_RP',RP_style,'.RData'))
#dataset <- 'All'
#save(RP,file = paste0('./results/RPs/RP_new_',dataset,'_',RP_style,'.RData'))


###SINGLE RP SERIES 
RPind <- data.frame(
  Year = RP$Year,
  RP1 = apply(RP[, -1], 1, function(x) mean(x, na.rm = TRUE))
)
write.csv(RPind, here("data", "RPind.csv"), row.names = FALSE)
fig_path <- file.path(cfg$folder_paths$figures_dir, "RP_ts.png")
png(filename = fig_path, width = 1200, height = 700, res = 150)
plot(RPind$Year, RPind$RP1, type = "l",
     xlab = "Year", ylab = "RP1")
dev.off()
head(RPind)





