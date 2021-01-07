# DEMAND PEAK MODELS #

#### 0. Libraries ####

library(mgcViz)
library(corrplot)
library(electBook)
library(RhpcBLASctl); blas_set_num_threads(1) # Optional
library(stringr)
library(magrittr)
library(matrixStats)
library(lubridate)
library(dplyr)
library(mvnfast)
library(rstudioapi)
library(gridExtra)
source("99_Utility_v2.R")

# Working directory to be changed adequately
setwd("C:/Users/yvenn/OneDrive/Bureau/Thèse/10. Articles/3. GitHub/Multi-resolution-forecasting/1. Sample Code") 

#### 1. Data Prepartion ####

data_list = prep()

datMxLearn = data_list[[1]] # low-resolution
datMxTest = data_list[[2]] # low-resolution
DataLearn = data_list[[3]] # high-resolution
DataTest = data_list[[4]] # high-resolution

#### 2. Model Fitting ####

#### Persistence baseline ####

dd = datMxTest %>% filter(ymd >= "2015-07-01")

# Peak Load Last Year
mape(dd$load,dd$peak24) # 4.38
rmse(dd$load,dd$peak24) # 0.343
mae(dd$load,dd$peak24) # 0.230

#### Multi-resolution ####
## GEV

MR_gev_function= function(df){
  return(gamV(list(load ~ dow + s(toy, k = 20) +
                          ti(matTem, matInt, k = c(10, 15), mc = c(TRUE, FALSE)) +
                          ti(matTem95, matInt, k = c(5, 5), mc = c(TRUE, FALSE)) +
                          ti(matLag, matInt,k = c(10, 10), mc = c(TRUE, FALSE)), 
                        ~ 1, ~ 1),
                   data = df, family = gevlss))}

MR_rolling_pred = rolling_origin(na.omit(datMxLearn),MR_gev_function, family='gevlss')

# save(file = "2. Pred Signals/1. GAMs/MR_rolling_pred_more_knots_factor.RData", MR_rolling_pred)

## Scaled-T

MR_scat_function= function(df){
  return(gamV(load ~ dow + s(toy, k = 20) +
                     ti(matTem, matInt, k = c(15, 10), mc = c(TRUE, FALSE)) +
                     ti(matTem95, matInt, k = c(5, 5), mc = c(TRUE, FALSE)) +
                     ti(matLag, matInt,k = c(10, 10), mc = c(TRUE, FALSE)),
              data = df, family = scat))}

MR_scat_rolling_pred = rolling_origin(na.omit(datMxLearn),MR_scat_function, family='scat')

# save(file = "2. Pred Signals/1. GAMs/MR_scat_rolling_pred_factor.RData", MR_scat_rolling_pred)

## Gaussian

MR_gauss_function= function(df){
  return(bamV(load ~ dow + s(toy, k = 20) +
                     ti(matTem, matInt, k = c(15, 10), mc = c(TRUE, FALSE)) +
                     ti(matTem95, matInt, k = c(5, 5), mc = c(TRUE, FALSE)) +
                     ti(matLag, matInt, k = c(10, 10), mc = c(TRUE, FALSE)),
              data = df, family = gaussian, aGam=list(discrete = T)))}

MR_rolling_pred = rolling_origin(datMxLearn,MR_gauss_function, family="gaussian")

# save(file = "2. Pred Signals/1. GAMs/MR_gauss_rolling_pred_more_knots_factor.RData", MR_rolling_pred)

#### Low-resolution ####

## Gaussian

low_resolution_gauss = function(df){
  return(bamV(load ~ dow + s(tod24) + s(toy, k= 20) +
                s(peak24, k=20) + s(tempMax, k=20) + s(temp95Max, k=20) + s(tempMin, k=20) + 
                s(temp95Min, k=20),
              data = df, aGam=list(discrete=T)))}

## Scaled-T

low_resolution_scat = function(df){
  return(gamV(load ~ dow + s(tod24) + s(toy, k= 20) +
                s(peak24, k=20) + s(tempMax, k=20) + s(temp95Max, k=20) + s(tempMin, k=20) + 
                s(temp95Min, k=20),
              data = df, family='scat'))}

## GEV

low_resolution_gev = function(df){
  return(gamV(list(load ~ dow + s(tod24) + s(toy, k= 20) +
                s(peak24, k=20) + s(tempMax, k=20) + s(temp95Max, k=20) + s(tempMin, k=20) + 
                s(temp95Min, k=20),
                ~dow + s(tod24) + s(toy, k= 20) +
                  s(peak24, k=20) + s(tempMax, k=20) + s(temp95Max, k=20) + s(tempMin, k=20) + 
                  s(temp95Min, k=20),
                ~1),
              data = df, family=gevlss))}

gauss_rolling_pred = rolling_origin(datMxLearn, low_resolution_gauss, family="gaussian")

scat_rolling_pred = rolling_origin(datMxLearn, low_resolution_scat, family="scat")

gev_rolling_pred = rolling_origin(datMxLearn, low_resolution_gev, family="gevlss")

# save(file = "2. Pred Signals/1. GAMs/gauss_rolling_pred_more_knots_factor.RData", gauss_rolling_pred)

# save(file = "2. Pred Signals/1. GAMs/scat_rolling_pred_more_knots_factor.RData", scat_rolling_pred)

# save(file = "2. Pred Signals/1. GAMs/gev_rolling_pred_more_knots_factor.RData", gev_rolling_pred)

#### High-resolution ####

high_res_function = function(df){
  return(bamV(load ~ factor_dow + factor_tod + 
                s(load24, k=20) +
                s(toy, k = 50) +
                s(temp, k = 20) +
                s(temp95, k = 20) +
                ti(temp, tod, k = c(10, 10)) +
                ti(temp95, tod, k = c(10, 10)) +
                ti(load24, tod, k = c(10, 10)) +
                ti(toy,tod, k = c(20, 20)),
              data = df, family=gaussian, aGam=list(discrete = T)))
}

high_res_pred = rolling_origin(DataLearn %>% mutate(factor_dow = factor(dow),factor_tod=factor(tod)),high_res_function, family="gaussian")

# save(file = "2. Pred Signals/1. GAMs/high_res_rolling_bam.RData", high_res_pred)

#### 3. Results ####

## High resolution

# Gauss
load("2. Pred Signals/1. GAMs/high_res_rolling_bam.RData")
test_and_pred = DataTest %>%
  mutate(pred = high_res_pred[[1]]) %>%
  select(ymd,load,pred) %>%
  group_by(ymd) %>%
  summarise_all(max) %>%
  ungroup(.) %>%
  mutate(ym = paste(year(ymd),"-",month(ymd), sep=""))
last_year(test_and_pred) # 2.43, 0.130, 0.155 
HR_gauss = bind_cols(metrics_over_time(datMxTest %>% mutate(pred = test_and_pred$pred), model='HR-Gauss', res='HR'),
                     AIC_metric_GAM(mod_list = high_res_pred))

# FCNN
HRFCNN = as.matrix(read.csv("2. Pred Signals/2. NNs/1. FC/high-resolution_FC_rolling_pred.csv",header=F))
test_and_pred = DataTest %>%
  mutate(pred = HRFCNN) %>%
  select(ymd,load,pred) %>%
  group_by(ymd) %>%
  summarise_all(max) %>%
  ungroup(.) %>%
  mutate(ym = paste(year(ymd),"-",month(ymd), sep=""))
last_year(test_and_pred) ### 1.43%, 0.0777, 0.0998
HR_FCNN = bind_cols(metrics_over_time(datMxTest %>% mutate(pred = test_and_pred$pred), model='HR-FCNN', res='HR'),
                    AIC_metric_NN("2. Pred Signals/2. NNs/1. FC/fitted/high-resolution_FC_rolling_","high",3051))

##### Low resolution

# FCNN
LRFCNN = as.matrix(read.csv("2. Pred Signals/2. NNs/1. FC/low-resolution_FC_rolling_pred.csv",header=F))
test_and_pred = datMxTest %>% mutate(pred = LRFCNN)
last_year(test_and_pred) # 2.11, 0.112, 0.144
LR_FCNN = bind_cols(metrics_over_time(datMxTest %>% mutate(pred = LRFCNN), model='LR-FCNN', res='LR'),
                    AIC_metric_NN("2. Pred Signals/2. NNs/1. FC/fitted/low-resolution_FC_rolling_","low",1751))

# Gauss
load("2. Pred Signals/1. GAMs/gauss_rolling_pred_more_knots_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = gauss_rolling_pred[[1]])
last_year(test_and_pred) # 2.26, 0.123, 0.144
LR_gauss = bind_cols(metrics_over_time(datMxTest %>% mutate(pred = gauss_rolling_pred[[1]]), model='LR-Gauss', res='LR'),
                     AIC_metric_GAM(mod_list = gauss_rolling_pred))

# Scat
load("2. Pred Signals/1. GAMs/scat_rolling_pred_more_knots_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = scat_rolling_pred[[1]])
last_year(test_and_pred) # 1.92, 0.105, 0.129
LR_scat = bind_cols(metrics_over_time(datMxTest %>% mutate(pred = scat_rolling_pred[[1]]), model='LR-Scat', res='LR'),
                    AIC_metric_GAM(mod_list = scat_rolling_pred))

# Gev
load("2. Pred Signals/1. GAMs/gev_rolling_pred_more_knots_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = gev_rolling_pred[[1]])
last_year(test_and_pred) # 2.67, 0.145, 0.169
LR_gev = bind_cols(metrics_over_time(datMxTest %>% mutate(pred = gev_rolling_pred[[1]]), model='LR-Gev', res='LR'),
                   AIC_metric_GAM(mod_list = gev_rolling_pred))

##### Multi-resolution

# CNN
MRCNN = apply(bind_cols(lapply(list.files(path = "2. Pred Signals/2. NNs/2. Hybrid/2. Without Trend/", pattern="*.csv", full.names=T), read.delim, header=F)),
              1,
              mean)

#MRCNN=as.matrix(read.csv("2. Pred Signals/2. NNs/2. Hybrid/multi-resolution_CNN_rolling_pred_1.csv",header=F))
test_and_pred = datMxTest %>% mutate(pred = MRCNN)
last_year(test_and_pred) # 1.56, 0.0844, 0.105
MR_CNN = bind_cols(metrics_over_time(datMxTest %>% mutate(pred = MRCNN), model='MR-CNN', res='MR'),
                   AIC_metric_NN("2. Pred Signals/2. NNs/2. Hybrid/fitted/multi-resolution_CNN_rolling_","multi",7315))

# Gauss
load("2. Pred Signals/1. GAMs/MR_gauss_rolling_pred_more_knots_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = MR_rolling_pred[[1]])
last_year(test_and_pred) # 1.42, 0.0765, 0.0963
MR_gauss = bind_cols(metrics_over_time(datMxTest %>% mutate(pred = MR_rolling_pred[[1]]), model='MR-Gauss', res='MR'),
                     AIC_metric_GAM(mod_list = MR_rolling_pred))

# Scat
load("2. Pred Signals/1. GAMs/MR_scat_rolling_pred_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = MR_scat_rolling_pred[[1]])
last_year(test_and_pred) # 1.41, 0.0755, 0.0959
MR_scat = bind_cols(metrics_over_time(datMxTest %>% mutate(pred = MR_scat_rolling_pred[[1]]), model='MR-Scat', res='MR'),
                    AIC_metric_GAM(mod_list = MR_scat_rolling_pred))

# Gev
load("2. Pred Signals/1. GAMs/MR_rolling_pred_more_knots_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = MR_rolling_pred[[1]])
last_year(test_and_pred) # 1.53, 0.0819, 0.103
MR_gev = bind_cols(metrics_over_time(datMxTest %>% mutate(pred = MR_rolling_pred[[1]]), model='MR-Gev', res='MR'),
                   AIC_metric_GAM(mod_list = MR_rolling_pred))

## Plots

all = bind_rows(MR_gev, MR_scat, MR_gauss, MR_CNN,
                LR_gev, LR_scat, LR_gauss, LR_FCNN,
                HR_gauss, HR_FCNN)

plot_metrics(all,"mape")
plot_metrics(all,"mae")
plot_metrics(all,"rmse")

gam_consistents = bind_rows(MR_gev, MR_scat, MR_gauss,
                            LR_gev, LR_scat, LR_gauss)

plot_metrics(gam_consistents,"AIC")

######### END