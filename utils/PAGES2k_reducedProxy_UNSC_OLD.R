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
library(edrGraphicalTools)
library(stringr)
library(spls)
library(superpc)
library(tidyr)
library(ggmap)
library(maptools)
library(maps)

set.seed(2018)
# specify options to treat the data
tStart = 1 #define start year (remember: the Common Era does not have a year 0).
tEnd   = 2000 #define end year for the analysis
tce  = tStart:tEnd
nce = length(tce)

#load the temperature data
#temp <- read.table('./data/had4_krig_ama_v2_0_0.txt', header=FALSE, sep="") #last century temperature
temp <- read.table('./data/HadCRUT.4.4.0.0.ama_ns_avg_1850-2015.txt', header=FALSE, sep="") #last century temperature ## HARMONIZE WITH THE YML FILE!!!!

#load the proxy data
proxydata =read.table('./data/proxy_ama_2.0.0.txt', header=TRUE, sep="") ## HARMONIZE WITH THE YML FILE!!!!
#proxydata =read.table('./data/proxy_ama_2.0.0_PAGES-crit-regional+FDR.txt', header=TRUE, sep="")
#proxydata =read.table('./data/proxy_ama_2.0.0_No_tree.txt', header=TRUE, sep="")

#load metadata
metadataproxy <- read.csv(file = 'data/metadata_2.0.0.csv',header = T) ### NOT NECCESSARY FOR US ATM!

#####PROXY MAP###########
metadataproxy2 <- as.data.frame(cbind(colnames(metadataproxy)[-1],t(metadataproxy[,-1])))
rownames(metadataproxy2) <- NULL
colnames(metadataproxy2) <- c('Name',as.vector(metadataproxy[,1]))
metadataproxy2$geo_latitude <- as.numeric(as.character(metadataproxy2$geo_latitude))
metadataproxy2$geo_longitude <- as.numeric(as.character(metadataproxy2$geo_longitude))

simbolos <- c('bivalve'=0,'borehole'=1,'coral'=2,'glacier ice'=3,'hybrid'=4,
              'lake sediment'=5,'marine sediment'=6,'sclerosponge'=7,
              'speleothem'=8,'tree'=9)

colores <- c('bivalve'='#F2D715','borehole'='#BBB87E','coral'='#FB9C05',
             'glacier ice'='#89B8E2','hybrid'='#2B8FCE','lake sediment'='#57649C',
             'marine sediment'='#955421','sclerosponge'='#EA132B',
             'speleothem'='#F81F6C','tree'='#489268')

#scale_color_brewer(type = 'div',palette = 'Spectral')
#guides(color=guide_legend(title='Type of Proxy:'))+
mapWorld <- borders("world", colour="gray50", fill="#F6F6F6") 
plotproxies <- ggplot(data = metadataproxy2)+mapWorld+
  geom_point(aes(x=geo_longitude, y=geo_latitude,color=archiveType,shape=archiveType), size=3,stroke=1)+
  scale_color_manual(name='Type of Proxy:',values = colores)+
  scale_shape_manual(name='Type of Proxy:',values = simbolos)+
  xlab('Longitude')+ylab('Latitude')+
  theme_bw()+
  theme(axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),
        legend.title = element_text(size=14),legend.text = element_text(size=12))
############################

t=proxydata[,1]
proxy=as.data.frame(proxydata[,-1])[which(t%in%tce),]
ny=dim(proxy)[1]
np=dim(proxy)[2]
#Get the years of the first available data for each proxy
yearMin= matrix(1, 1, np) 
for (i in 1:np) {
  yearMin[i]=min(which(is.finite(proxy[,i])))
}


#####HISTOGRAM FIRST YEAR#############
datahist <- data.frame(Year = c(yearMin))
histogram_first <- ggplot(data = datahist)+
  geom_histogram(mapping = aes(x = Year),bins = 10,fill='white',color='black')+
  ylab('')+
  theme_bw()+
  theme(axis.title.x = element_text(size=18),axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))
########AREA PLOT##############################
etiquetas <- as.character(metadataproxy2$archiveType)
proxy2 <- proxy
for(x in 1:257){
  proxy2[!is.na(proxy2[,x]),x] <- etiquetas[x]
}
#scale_fill_brewer(type = 'div',palette = 'Spectral')+
proxy2 <- as.data.frame(t(proxy2))
proxy2 <- proxy2 %>% mutate(Name=rownames(proxy2))
proxy2 <- proxy2 %>% gather(key = 'Year',value = 'Type',-Name) %>%
  filter(!is.na(Type))
taproxies <- proxy2 %>% dplyr::select(-Name) %>% mutate(Year=as.numeric(Year)) %>% 
  group_by(Year,Type) %>% summarise(Number=n()) %>% arrange(Year,Type)
plottaproxies <- ggplot(data = taproxies, mapping =aes(x = Year,y = Number,fill=Type))+
  geom_area()+
  scale_fill_manual(values = colores)+
  theme_bw()+  theme(axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),
                     legend.title = element_text(size=14),legend.text = element_text(size=12),
                     legend.position = 'none')
#Combination map-area plot:
#combinedmap_area <- plot_grid(plotproxies,plottaproxies,nrow=2,
#                              rel_widths = c(1,0.45),rel_heights = c(1,0.45))

totalproxies <- taproxies %>% group_by(Year) %>% summarise(total=sum(Number))
# define analysis options
RP_style="SIR" #possible choices: "PCR", "LASSO", "SIR", "SPLS"
## DEFINE FUNCTIONS USED LATER
#define function %!in% which is same as ~ismember
'%!in%' <- function(x,y)!('%in%'(y,x))

#define rfind fun ction as similar to find in matlab
rfind <- function(x)seq(along=x)[as.logical(x)] 


# scale the proxy records to unit variance and zero mean
proxy_scaled = scale(proxy)

# count the number of available datapoints in each record
navl = matrix(1, 1, np) 
for (i in 1:np)
{  
  navl[i] = sum(is.finite(proxy_scaled[,i]))
}

##  MAKE A TIME-DEPENDENT COMPOSITE
chunk = 250 # define  segment length
ns = ny/chunk  # number of segments
calib = c(1900:2000)#Define calibration year
RP = matrix(NA,ny,ns)#Initialise the reduced proxies
nprox= matrix(0,1,ns)#Initialise nprox that is later used to store number of proxies used in each segmentation
temp_lastc <- temp[which(temp[,1]%in%calib),2] #last century temperature

#evalues=matrix(1,1,npcs)#Initialise a matrix to store the eigenvalues of PCs
coefflist <- list()
npcs <- c(11,13,9,11,7,12,10,6) #first PC whose adjusted R2 > 70%
limitup <- c(1,2,2,2,2,1,1,2) ##Marginal Dimension Tests (SIR)
nslicesv <- c(5,7)
proport <- c(1,1,1,1,1,0.35,0.6,0.75)
pZerov <- c(30,30,30,30,30,30,30,30)
#tsholdv <- c(0.53,1.16,0.36,2.62,1.16)
tsholdv <- NULL
npcs2 <- c(3,3,3,3,3,3,3,3)

for (k in 1:ns)
{   # loop over segments
  show(k)
  tMin = tStart+(k-1)*chunk
  segProxies = rfind(yearMin <=tMin) 
  timeSpan  = rfind(tce>=tMin)
  nt = length(timeSpan)
  one=matrix(1,nt,1)#Initialise a ny*1 one matrix that can be later appended to PC matrix
  proxy_finite= proxy_scaled[timeSpan,segProxies]
  cleanNANs <- apply(X = proxy_finite,MARGIN = 2,function(x) sum(is.na(x))/dim(proxy_finite)[1])<0.05
  proxy_finite <- proxy_finite[,cleanNANs]
  time_calib <- rfind(timeSpan %in% calib)
  proxy_finite_calib= proxy_finite[time_calib,]
  cleanNANs_calib <- apply(X = proxy_finite_calib,MARGIN = 2,function(x) sum(is.na(x))/dim(proxy_finite_calib)[1])<0.05
  proxy_finite <- proxy_finite[,cleanNANs_calib]
  proxy_finite[!is.finite(proxy_finite)]=0
  nprox[k] = dim(proxy_finite)[2]
  
  if (RP_style=="PCR") { 
      pca=prcomp(proxy_finite)
      a <- lm(temp_lastc~pca$x[rfind(timeSpan %in% calib),1:npcs[k]])
      coeff <- as.numeric(a$coefficients)
      RP[timeSpan,k]   = cbind(one,pca$x[,1:npcs[k]])%*%coeff
  }
  if (RP_style=="sPCR") {
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
  
    
  if (RP_style=="LASSO"){
    #Performs a Lasso regularization technique to avoid the use of PCA in multiple regression
    Xmatrix <- proxy_finite[rfind(timeSpan %in% calib),]
    cv.lasso <- cv.glmnet(x = Xmatrix,y = temp_lastc,alpha=1,parallel = T)
    coeff <- coef.cv.glmnet(cv.lasso,s = 'lambda.min')
    coefflist[[k]] <- coeff
    RP[timeSpan,k] <- predict(cv.lasso,proxy_finite,s = 'lambda.min')
  }
  if (RP_style=="SIR"){
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

save(RPn,file=paste0('Plot_RP',RP_style,'.RData'))
dataset <- 'All'
save(RP,file = paste0('./results/RPs/RP_new_',dataset,'_',RP_style,'.RData'))


###SINGLE SERIES 
RPind <- apply(RP[,-1], 1, function(x) mean(x,na.rm=T))
save(RPind,file = paste0('./results/RPs/RPind_',dataset,'_',RP_style,'.RData'))
plot(RPind,type='l')
