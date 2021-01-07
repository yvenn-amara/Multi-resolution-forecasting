# INSTANT OF PEAK MODELS #

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

# Instant of Peak Last Year
acc(dd$tod + 1,dd$tod24 + 1) # 55.9
mape(dd$tod  + 1,dd$tod24  + 1) # 8.71
rmse(dd$tod  + 1,dd$tod24  + 1) # 5.36
mae(dd$tod  + 1,dd$tod24  + 1) # 2.49

#### Multi-resolution ####

ocat_function = function(df){
  return(gamV(todFrom1 ~ dow + s(toy, k= 20) + s(tod24) +
                ti(matTem, matInt, k = c(15, 10), mc = c(TRUE, FALSE)) +
                ti(matTem95, matInt, k = c(5, 5), mc = c(TRUE, FALSE)) +
                ti(matLag, matInt, k = c(10, 10), mc = c(TRUE, FALSE)), 
              data = df, family = ocat(R = 48)))}

ocat_rolling_pred = rolling_origin(datMxLearn,ocat_function, family="ocat")
# save(file = "2. Pred Signals/1. GAMs/ocat_rolling_pred_more_knots.RData", ocat_rolling_pred)

#### Low-resolution ####

low_resolution_ocat = function(df){
  return(gamV(todFrom1 ~ dow + s(tod24) + s(toy, k= 20) +
                s(peak24, k=20) + s(tempMax, k=20) + s(temp95Max, k=20) + s(tempMin, k=20) + 
                s(temp95Min, k=20),
              data = df, family = ocat(R = 48)))}

lr_ocat_rolling_pred = rolling_origin(datMxLearn, low_resolution_ocat, family="ocat")

# save(file = "2. Pred Signals/1. GAMs/low_resolution_ocat.RData", lr_ocat_rolling_pred)

#### High-resolution ####
# Fitted in the DEMAND PEAK MODELS file

#### 3. Results ####

#### High-resolution ####

## Gaussian
load("2. Pred Signals/1. GAMs/high_res_rolling_bam.RData")
test_and_pred = high_res_to_peak(DataTest, high_res_pred[[1]])
last_year_instant(test_and_pred) #Acc: 17.4 %, MAPE:  7.99, RMSE: 4.81 MW, MAE 2.48 MW 
HR_gauss_metrics = metrics_over_time_instant(datMxTest %>% mutate(pred = test_and_pred$pred), model='HR-Gauss', res='HR')


# FCNN
HRFCNN = as.vector(as.matrix(read.csv("2. Pred Signals/2. NNs/1. FC/high-resolution_FC_rolling_pred.csv",header=F)))
test_and_pred = high_res_to_peak(DataTest, HRFCNN)
last_year_instant(test_and_pred) #Acc: 10.3 %, MAPE: 8.07 %, RMSE: 4.63 MW, MAE: 2.53 MW 
HRFCNN_metrics = metrics_over_time_instant(datMxTest %>% mutate(pred = test_and_pred$pred), model='HR-FCNN', res='HR')

#### Low-resolution ####

# Ocat
load("2. Pred Signals/1. GAMs/low_resolution_ocat.RData")
test_and_pred = datMxTest %>%
  mutate(pred = lr_ocat_rolling_pred[[1]])
last_year_instant(test_and_pred) #Acc: 53.8 %, MAPE: 7.24 %, RMSE, 4.22 MW, MAE: 2.11 MW
LR_ocat_metrics = bind_cols(metrics_over_time_instant(datMxTest %>% mutate(pred = test_and_pred$pred), model='LR-Ocat', res='LR'),
                            AIC_metric_GAM_instant(lr_ocat_rolling_pred))

# FCNN
LRFCNN = to_classes(round(as.matrix(read.csv("2. Pred Signals/2. NNs/1. FC/low-resolution_FC_rolling_pred_Instant.csv",header=F))))
test_and_pred = datMxTest %>% mutate(pred = LRFCNN)
last_year_instant(test_and_pred) #Acc: 60 %, MAPE: 6.45 %, RMSE: 4.4 MW, MAE: 1.94 MW
LRFCNN_metrics = metrics_over_time_instant(datMxTest %>% mutate(pred = test_and_pred$pred), model='LR-FCNN', res='LR')

#### Multi-resolution ####

# Ocat
load("2. Pred Signals/1. GAMs/ocat_rolling_pred_more_knots.RData")
test_and_pred = datMxTest %>%
  mutate(pred = ocat_rolling_pred[[1]])
last_year_instant(test_and_pred) #Acc: 53.2%, MAPE: 6.83%, RMSE: 4.08 MW, MAE 2.01 MW
MR_ocat_metrics = bind_cols(metrics_over_time_instant(datMxTest %>% mutate(pred = test_and_pred$pred), model='MR-Ocat', res='MR'),
                            AIC_metric_GAM_instant(ocat_rolling_pred))

# MRCNN
MRCNN= to_classes(round(as.matrix(read.csv("2. Pred Signals/2. NNs/2. Hybrid/multi-resolution_CNN_rolling_pred_Instant_1.csv",header=F))))
test_and_pred = datMxTest %>% mutate(pred = MRCNN)
last_year_instant(test_and_pred) #Acc: 60.9 %, MAPE: 5.82 %, RMSE: 3.85 MW, MAE: 1.7 MW 
MRCNN_metrics = metrics_over_time_instant(datMxTest %>% mutate(pred = test_and_pred$pred), model='MR-CNN', res='MR')


## Plots
all = bind_rows(HRFCNN_metrics,LRFCNN_metrics,MRCNN_metrics,
                HR_gauss_metrics,LR_ocat_metrics,MR_ocat_metrics)

plot_metrics(all,"acc")
plot_metrics(all,"mae")
plot_metrics(all,"rmse")

######### END