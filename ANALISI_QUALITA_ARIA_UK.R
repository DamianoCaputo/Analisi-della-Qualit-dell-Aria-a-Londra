#LIBRERIE---------------
library(openair)
library(readr)
library(worldmet)
library(tidyverse)
library(GGally)
library(openairmaps)
library(FactoMineR)
library(FactoInvestigate)
library(Factoshiny)
library(factoextra)
library(pls)
library(gamlss)
library(ISLR)
library(dplyr)
library(tidyr)
library(mhsmm)
library(zoo)
library(lmerTest)
library(MixGHD)
library(mclust)
library(gridExtra)

#DATASET----------------
uk_2015_20 <- importAURN(site = "kc1", year = 2015:2020, meta = TRUE)
uk_15_20 <- uk_2015_20

##SALVA FILE----
write.csv(uk_2015_20,file ="uk_air_quality_2015_20.csv")

##PULIZIA DATI----
uk_15_20 <- uk_15_20[,-c(12,13,14,15,21)]

na_count <-sapply(uk_15_20, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)

##DATI NORMALIZZATI----
logT_dta <- log(uk_15_20[,4:11]+1)

logT_dta$pm2.5[logT_dta$pm2.5 == -Inf] <- 0
logT_dta$pm10[logT_dta$pm10 == -Inf] <- 0

#ANALISI GENERALE----
##SUMMARY PLOT-----
summaryPlot(uk_15_20[,-c(15,16)], type="density")

##TIME SERIES------ 
timeVariation(uk_15_20,pollutant="co")
timeVariation(uk_15_20,pollutant="nox")
timeVariation(uk_15_20,pollutant="no2")
timeVariation(uk_15_20,pollutant="no")
timeVariation(uk_15_20,pollutant="o3")
timeVariation(uk_15_20,pollutant="so2")
timeVariation(uk_15_20,pollutant="pm2.5")
timeVariation(uk_15_20,pollutant="pm10")

##TIME SERIES TOTALI------
timeVariation(uk_15_20,pollutant=c("co","nox","no2","no","o3",
                                  "so2","pm2.5","pm10"),type="year")

timePlot(uk_15_20, pollutant = c("co","nox","no2","no","o3","so2",
                                 "pm2.5","pm10"),
         avg.time = "year", normalise = "1/1/1998",lwd = 4, lty = 1,
         group = TRUE, ylim = c(0, 200))
##CALENDAR PLOT-----
#SO2
calendarPlot(uk_15_20, pollutant = "so2", year = 2015)
calendarPlot(uk_15_20, pollutant = "so2", year = 2016)
calendarPlot(uk_15_20, pollutant = "so2", year = 2017)
calendarPlot(uk_15_20, pollutant = "so2", year = 2018)
calendarPlot(uk_15_20, pollutant = "so2", year = 2019)
calendarPlot(uk_15_20, pollutant = "so2", year = 2020)

#O3
calendarPlot(uk_15_20, pollutant = "o3", year = 2015)
calendarPlot(uk_15_20, pollutant = "o3", year = 2016)
calendarPlot(uk_15_20, pollutant = "o3", year = 2017)
calendarPlot(uk_15_20, pollutant = "o3", year = 2018)
calendarPlot(uk_15_20, pollutant = "o3", year = 2019)
calendarPlot(uk_15_20, pollutant = "o3", year = 2020)

##CORRELAZIONI------
dta <- uk_15_20[,-c(1,2,3,15,16)]

corr <- corPlot(dta,main="Correlazioni")

ggpairs(dta)+theme_bw()

##WIND ROSE-----
windRose(uk_15_20,type = "year",pollutant="pm10")
windRose(uk_15_20, type = "year",pollutant="pm2.5")
windRose(uk_15_20, type = "year",pollutant="co")
windRose(uk_15_20, type = "year",pollutant="no")
windRose(uk_15_20, type = "year",pollutant="nox")
windRose(uk_15_20, type = "year",pollutant="so2")
windRose(uk_15_20, type = "year",pollutant="no2")
windRose(uk_15_20, type = "year",pollutant="o3")

##POLAR PLOT------
#CONCENTRAZIONE DEGLI INQUINANTI IN BASE AI VENTI 

polarPlot(uk_15_20,pollutant="co",x="ws",wd="wd",type="year")
polarPlot(uk_15_20,pollutant="nox",x="ws",wd="wd",type="year")
polarPlot(uk_15_20,pollutant="no2",x="ws",wd="wd",type="year")
polarPlot(uk_15_20,pollutant="no",x="ws",wd="wd", type="year")
polarPlot(uk_15_20,pollutant="o3",x="ws",wd="wd", type="year")
polarPlot(uk_15_20,pollutant="so2",x="ws",wd="wd", type="year")
polarPlot(uk_15_20,pollutant="pm2.5",x="ws",wd="wd", type="year")
polarPlot(uk_15_20,pollutant="pm10",x="ws",wd="wd", type="year")


##MAPPA INTERATTIVA-----
polarMap(uk_15_20,pollutant="o3",latitude = "latitude",
         longitude = "longitude",type = "site",iconWidth=300,
         iconHeight=300,fig.width=8,fig.height=8)

##TREND ANALYSIS------
#Usa il bootstrap 
##CLASSICO
TheilSen(uk_15_20, pollutant = "co", 
         ylab = "co", 
         deseason = TRUE,
         date.format = "%Y",main="Trend Del CO")

TheilSen(uk_15_20, pollutant = "nox", 
         ylab = "nox", 
         deseason = TRUE,
         date.format = "%Y",main="Trend Degli NOx")

TheilSen(uk_15_20, pollutant = "no", 
         ylab = "no", 
         deseason = TRUE,
         date.format = "%Y",main="Trend Degli NO")

TheilSen(uk_15_20, pollutant = "no2", 
         ylab = "no2", 
         deseason = TRUE,
         date.format = "%Y",main="Trend Del NO2")

TheilSen(uk_15_20, pollutant = "o3", 
         ylab = "o3", 
         deseason = TRUE,
         date.format = "%Y",main="Trend Dell'Ozono")

TheilSen(uk_15_20, pollutant = "so2", 
         ylab = "so2", 
         deseason = TRUE,
         date.format = "%Y",main="Trend Dell'Anidride Solforosa")

TheilSen(uk_15_20, pollutant = "pm10", 
         ylab = "pm10", 
         deseason = TRUE,
         date.format = "%Y",main="Trend del PM10")

TheilSen(uk_15_20, pollutant = "pm2.5", 
         ylab = "pm2.5", 
         deseason = TRUE,
         date.format = "%Y",main="Trend del PM2.5")
##SMOOTH TREND
smoothTrend(uk_15_20, pollutant = c("co", "nox", "no2"), 
            type = c("season"),
            date.breaks = 3, lty = 0)

smoothTrend(uk_15_20, pollutant = c( 'no',"o3", "so2"), 
            type = c("season"),
            date.breaks = 3, lty = 0)

smoothTrend(uk_15_20, pollutant = c( 'pm10',"pm2.5"), 
            type = c("season"),
            date.breaks = 3, lty = 0)
#PCA----
dta_pca <- PCA(uk_15_20[,4:14])
Investigate(dta_pca)

pca_uk5 <- PCA(uk_15_20[, -c(15,16)],quali.sup = c(1,2,3),
               quanti.sup = c(12,13,14), graph = F)
plot.PCA(pca_uk5,choix='var',col.quanti.sup='#0000FF')
plot.PCA(pca_uk5,invisible=c('quali','ind.sup'),label =c('ind'))
# Rapp. dimensioni dati 
fviz_eig(pca_uk5, addlabels = TRUE)#Rapp. Dimensioni  ok 

#Qualità dalle rappresentazioni 
fviz_cos2(pca_uk5, choice = "var") #Qualita delle rappresentazioni ok

#ordine delle variabili per dimensioni 
corrplot(pca_uk5$var$cos2, is.corr=FALSE)#ordine delle dimensioni ok

fviz_pca_var(pca_uk5, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlappin
)#Importanza aper cos2
#pca con le quali suo e quanti sup

##CLUSTERING PCA-----
#Clustering delle varibili 
# Create a grouping variable using kmeans
# Create 3 groups of variables (centers = 3)
set.seed(123)
res.km <- kmeans(pca_uk5$var$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_var(pca_uk5, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster")

#PCR----
#Principal Component Regression
#con senso a lvl chimico 
set.seed(12345)
uk_pcr <- uk_15_20[,c(4:11)]

pcr_co <- pcr(co~., data=uk_pcr,scale=TRUE, validation = "CV")
summary(pcr_co)
validationplot(pcr_co)#Validation plot
validationplot(pcr_co,val.type="MSEP")#MSEP validation 
validationplot(pcr_co,val.type="R2") #R2 validation 
predplot(pcr_co)#predizioni delle variabili 
coefplot(pcr_co)#coefficenti di regressione 

pcr_nox <- pcr(nox~., data=uk_pcr,scale=TRUE, validation = "CV")
summary(pcr_nox)
validationplot(pcr_nox)#Validation plot
validationplot(pcr_nox,val.type="MSEP")#MSEP validation 
validationplot(pcr_nox,val.type="R2") #R2 validation 
predplot(pcr_nox)#predizioni delle variabili 
coefplot(pcr_nox)#coefficenti di regressione 

pcr_no2 <- pcr(no2~., data=uk_pcr,scale=TRUE, validation = "CV")
summary(pcr_no2)
validationplot(pcr_no2)#Validation plot
validationplot(pcr_no2,val.type="MSEP")#MSEP validation 
validationplot(pcr_no2,val.type="R2") #R2 validation 
predplot(pcr_no2)#predizioni delle variabili 
coefplot(pcr_no2)#coefficenti di regressione 

pcr_no <- pcr(no~., data=uk_pcr,scale=TRUE, validation = "CV")
summary(pcr_no)
validationplot(pcr_no)
validationplot(pcr_no,val.type="MSEP")
validationplot(pcr_no,val.type="R2")
predplot(pcr_no)
coefplot(pcr_no)

pcr_o3 <- pcr(o3~., data=uk_pcr,scale=TRUE, validation = "CV")
summary(pcr_o3)
validationplot(pcr_o3)#Validation plot
validationplot(pcr_o3,val.type="MSEP")#MSEP validation 
validationplot(pcr_o3,val.type="R2") #R2 validation 
predplot(pcr_o3)#predizioni delle variabili 
coefplot(pcr_o3)#coefficenti di regressione 

pcr_so2 <- pcr(so2~., data=uk_pcr,scale=TRUE, validation = "CV")
summary(pcr_so2)
validationplot(pcr_so2)#Validation plot
validationplot(pcr_so2,val.type="MSEP")#MSEP validation 
validationplot(pcr_so2,val.type="R2") #R2 validation 
predplot(pcr_so2)#predizioni delle variabili 
coefplot(pcr_so2)#coefficenti di regressione 

pcr_pm10 <- pcr(pm10~., data=uk_pcr,scale=TRUE, validation = "CV")
summary(pcr_pm10)
validationplot(pcr_pm10)#Validation plot
validationplot(pcr_pm10,val.type="MSEP")#MSEP validation 
validationplot(pcr_pm10,val.type="R2") #R2 validation 
predplot(pcr_pm10)#predizioni delle variabili 
coefplot(pcr_pm10)#coefficenti di regressione 


pcr_pm2.5 <- pcr(pm2.5~., data=uk_pcr,scale=TRUE, validation = "CV")
summary(pcr_pm2.5)
validationplot(pcr_pm2.5)#Validation plot
validationplot(pcr_pm2.5,val.type="MSEP")#MSEP validation 
validationplot(pcr_pm2.5,val.type="R2") #R2 validation 
predplot(pcr_pm2.5)#predizioni delle variabili 
coefplot(pcr_pm2.5)#coefficenti di regressione

#REGRESSIONE----
dta_xNA <- na.omit(uk_15_20)
model_norm_AT <- gamlss( air_temp ~ co+nox+no2+no+o3+so2+pm10+pm2.5,
                         data=dta_xNA,family = NO)

model_norm2_AT <- gamlss( air_temp ~ co+o3+so2+pm10+pm2.5,
                          data=dta_xNA,family = NO)

model_norm_3 <- gamlss(ws~co+o3+pm10+pm2.5,data=dta_xNA,family = NO)

plot(model_norm_AT)
plot(model_norm2_AT)
plot(model_norm_3)

#HMM-----
##k=2-----
K=2
start.val <- hmmspec(init = rep(1/K, K),
                     trans = matrix(1/K, nrow = K, ncol = K),
                     parms.emis = list(mu = c(0.1,0.3),
                                       sigma=c(1,1,1)),
                     dens.emis = dnorm.hsmm)
hmm_CO_K2 <- hmmfit(logT_dta$co, start.val, mstep = mstep.norm)

plot(hmm_CO_K2$loglik,type="b",ylab="Log-likelihood",xlab="Iteration")

states_CO_K2 <- hmm_CO_K2$yhat

plot(logT_dta$co,col=states_CO_K2)

##CO-----
plot(logT_dta$co)
abline(h=c(0.1,0.25,0.4))

K=3
start.val <- hmmspec(init = rep(1/K, K),
                     trans = matrix(1/K, nrow = K, ncol = K),
                     parms.emis = list(mu = c(0.1,0.25,0.4),
                                       sigma=c(1,1,1)),
                     dens.emis = dnorm.hsmm)
hmm_CO_K3 <- hmmfit(logT_dta$co, start.val, mstep = mstep.norm)

plot(hmm_CO_K3$loglik,type="b",ylab="Log-likelihood",xlab="Iteration")

states_CO_K3 <- hmm_CO_K3$yhat

plot(logT_dta$co,col=states_CO_K3)

##NOx-----
plot(logT_dta$nox)
abline(h=c(2.5,3.5,4.5))

K=3
start.val <- hmmspec(init = rep(1/K, K),
                     trans = matrix(1/K, nrow = K, ncol = K),
                     parms.emis = list(mu =c(2.5,3.5,4.5),
                                       sigma=c(1,1,1)),
                     dens.emis = dnorm.hsmm)
hmm_NOx_K3 <- hmmfit(logT_dta$nox, start.val, mstep = mstep.norm)

plot(hmm_NOx_K3$loglik,type="b",ylab="Log-likelihood",xlab="Iteration")

states_NOx_K3 <- hmm_NOx_K3$yhat

plot(logT_dta$nox,col=states_NOx_K3)

##NO2-----
plot(logT_dta$no2)
abline(h=c(2,3,4))
K=3
start.val <- hmmspec(init = rep(1/K, K),
                     trans = matrix(1/K, nrow = K, ncol = K),
                     parms.emis = list(mu = c(2,3,4),sigma=c(1,1,1)),
                     dens.emis = dnorm.hsmm)
hmm_NO2_K3 <- hmmfit(logT_dta$no2, start.val, mstep = mstep.norm)

plot(hmm_NO2_K3$loglik,type="b",ylab="Log-likelihood",xlab="Iteration")

states_NO2_K3 <- hmm_NO2_K3$yhat

plot(logT_dta$no2,col=states_NO2_K3)

##NO----
plot(logT_dta$no)
abline(h=c(0.5,1.5,2.5))
K=3
start.val <- hmmspec(init = rep(1/K, K),
                     trans = matrix(1/K, nrow = K, ncol = K),
                     parms.emis = list(mu = c(0.5,1.5,2.5),
                                       sigma=c(1,1,1)),
                     dens.emis = dnorm.hsmm)
hmm_NO_K3 <- hmmfit(logT_dta$no, start.val, mstep = mstep.norm)

plot(hmm_NO_K3$loglik,type="b",ylab="Log-likelihood",xlab="Iteration")

states_NO_K3 <- hmm_NO_K3$yhat

plot(logT_dta$no,col=states_NO_K3)

##O3----
plot(logT_dta$o3)
abline(h=c(2.5,3.5,4.2))

K=3
start.val <- hmmspec(init = rep(1/K, K),
                     trans = matrix(1/K, nrow = K, ncol = K),
                     parms.emis = list(mu =c(2.5,3.5,4.2),
                                       sigma=c(1,1,1)),
                     dens.emis = dnorm.hsmm)
hmm_O3_K3 <- hmmfit(logT_dta$o3, start.val, mstep = mstep.norm)

plot(hmm_O3_K3$loglik,type="b",ylab="Log-likelihood",xlab="Iteration")

states_O3_K3 <- hmm_O3_K3$yhat

plot(logT_dta$o3,col=states_O3_K3)

##SO2----
plot(logT_dta$so2)
abline(h=c(0,1,2))

K=3
start.val <- hmmspec(init = rep(1/K, K),
                     trans = matrix(1/K, nrow = K, ncol = K),
                     parms.emis = list(mu = c(0,1,2),sigma=c(1,1,1)),
                     dens.emis = dnorm.hsmm)
hmm_SO2_K3 <- hmmfit(logT_dta$so2, start.val, mstep = mstep.norm)

plot(hmm_SO2_K3$loglik,type="b",ylab="Log-likelihood",xlab="Iteration")

states_SO2_K3 <- hmm_SO2_K3$yhat

plot(logT_dta$so2,col=states_SO2_K3)

##PM10----
plot(logT_dta$pm10)
abline(h=c(2,3,4))

K=3
start.val <- hmmspec(init = rep(1/K, K),
                     trans = matrix(1/K, nrow = K, ncol = K),
                     parms.emis = list(mu = c(2,3,4),sigma=c(1,1,1)),
                     dens.emis = dnorm.hsmm)
hmm_PM10_K3 <- hmmfit(logT_dta$pm10, start.val, mstep = mstep.norm)

plot(hmm_PM10_K3$loglik,type="b",ylab="Log-likelihood",xlab="Iteration")

states_PM10_K3 <- hmm_PM10_K3$yhat

plot(logT_dta$pm10,col=states_PM10_K3)

##PM2.5----
plot(logT_dta$pm2.5)
abline(h=c(1,2,3))

K=3
start.val <- hmmspec(init = rep(1/K, K),
                     trans = matrix(1/K, nrow = K, ncol = K),
                     parms.emis = list(mu = c(1,2,3),sigma=c(1,1,1)),
                     dens.emis = dnorm.hsmm)
hmm_PM25_K3 <- hmmfit(logT_dta$pm2.5, start.val, mstep = mstep.norm)

plot(hmm_PM25_K3$loglik,type="b",ylab="Log-likelihood",xlab="Iteration")

states_PM25_K3 <- hmm_PM25_K3$yhat

plot(logT_dta$pm2.5,col=states_PM25_K3)

##BOOTSTRAP O3-------
B <- 200 # replicates
mu.boot_O3_K3 <- matrix(NA,3,B)
sigma.boot_O3_K3 <- matrix(NA,3,B)
for (b in 1:B)
{
  true.par <- hmmspec(init = hmm_O3_K3$model$init,
                      trans = hmm_O3_K3$model$transition,
                      parms.emis = list(mu = hmm_O3_K3$model$parms.emission$mu,
                                        sigma=hmm_O3_K3$model$parms.emission$sigma),
                      dens.emis = dnorm.hsmm)
  train <- simulate(true.par, nsim = 107, seed = 1234, rand.emis = rnorm.hsmm)
  mod.boot_O3_K3 <- hmmfit(train, true.par, mstep = mstep.norm)
  mu.boot_O3_K3[,b] <- mod.boot_O3_K3$model$parms.emission$mu
  sigma.boot_O3_K3[,b] <- mod.boot_O3_K3$model$parms.emission$sigma
}
apply(mu.boot_O3_K3,1,mean)
apply(sigma.boot_O3_K3,1,mean)

hmm_O3_K3$model$parms.emission

#HMM MULTIVARIATI-------------
##FUNZIONI-----
dmvnorm.hsmm <- function(x, j, model) 
{
  dmvnorm(x,mean = model$parms.emission$mu[[j]],
          sigma = model$parms.emission$sigma[[j]])
}

mstep.mvnorm <- function(x, wt) {
  emission <- list(mu = list(), sigma = list())
  for(i in 1:ncol(wt)) {
    tmp <- cov.wt(x, wt[, i])
    emission$mu[[i]] <- tmp$center
    emission$sigma[[i]] <- tmp$cov
  }
  emission
}
##TOTALI-----
dta_xNA <- logT_dta
dta_xNA <- na.aggregate(dta_xNA)

K=3
start.val <- hmmspec(init = rep(1/K, K),
                     trans = matrix(1/K, nrow = K, ncol = K),
                     parms.emis = list(mu = list(c(0.1,2.5,2,0.5,2.5,
                                                   0,2,1),
                                                 c(0.25,3.5,3,1.5,
                                                   3.5,1,3,2),
                                                 c(0.4,4.5,4,2.5,
                                                   4.2,2,4,3)),
                                       sigma=list(diag(8), diag(8), 
                                                  diag(8))),
                     dens.emis = dmvnorm.hsmm)

hmm_mv_k3 <- hmmfit(matrix(unlist(dta_xNA),ncol=8), start.val, 
                    mstep = mstep.mvnorm)

states_mv_k3 <- hmm_mv_k3$yhat

plot(dta_xNA$no2,col=states_mv_k3)
plot(dta_xNA$so2,col=states_mv_k3)
plot(dta_xNA$pm10,col=states_mv_k3)
plot(dta_xNA$o3,col=states_mv_k3)
plot(dta_xNA$co,col=states_mv_k3)
plot(dta_xNA$no,col=states_mv_k3)
plot(dta_xNA$nox,col=states_mv_k3)
plot(dta_xNA$PM2.5,col=states_mv_k3)






