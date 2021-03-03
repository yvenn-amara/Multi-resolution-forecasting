# SIGNIFICANCE TESTING #

#### 0. Libraries and Data Prep ####

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
library(forecast)
library(multDM)
library(formattable)

# Working directory to be changed adequately
setwd("C:/Users/yvenn/OneDrive/Bureau/Thèse/10. Articles/3. GitHub/Multi-resolution-forecasting/1. Sample Code")
source("99_Utility_v3.R")

data_list = prep()

datMxLearn = data_list[[1]] # low-resolution
datMxTest = data_list[[2]] # low-resolution
DataLearn = data_list[[3]] # high-resolution
DataTest = data_list[[4]] # high-resolution

#### 1. Importing Prediction signals ####

filter_test = "2015-07-01"

#### DP ####

## High resolution
# ARIMA
load("2. Pred Signals/0. ARIMA/hr_arima.RData")
test_and_pred = DataTest %>%
  mutate(pred = HR_arima[[1]]) %>%
  select(ymd,load,pred) %>%
  group_by(ymd) %>%
  summarise_all(max) %>%
  ungroup(.) %>%
  mutate(ym = paste(year(ymd),"-",month(ymd), sep="")) %>% 
  filter(ymd >= filter_test)
DP_HR_arima = test_and_pred$pred
DP_HR_arima_error = test_and_pred$load - test_and_pred$pred

# Gauss
load("2. Pred Signals/1. GAMs/high_res_rolling_bam.RData")
test_and_pred = DataTest %>%
  mutate(pred = high_res_pred[[1]]) %>%
  select(ymd,load,pred) %>%
  group_by(ymd) %>%
  summarise_all(max) %>%
  ungroup(.) %>%
  mutate(ym = paste(year(ymd),"-",month(ymd), sep="")) %>% 
  filter(ymd >= filter_test)
DP_HR_gauss_error =  test_and_pred$load - test_and_pred$pred
DP_HR_gauss = test_and_pred$pred

# FCNN
HRFCNN = as.matrix(read.csv("2. Pred Signals/2. NNs/1. FC/high-resolution_FC_rolling_pred.csv",header=F))
test_and_pred = DataTest %>%
  mutate(pred = HRFCNN) %>%
  select(ymd,load,pred) %>%
  group_by(ymd) %>%
  summarise_all(max) %>%
  ungroup(.) %>%
  mutate(ym = paste(year(ymd),"-",month(ymd), sep="")) %>% 
  filter(ymd >= filter_test)
DP_HRFCNN_error =  test_and_pred$load - test_and_pred$pred
DP_HRFCNN = test_and_pred$pred

##### Low resolution
# ARIMA
load("2. Pred Signals/0. ARIMA/lr_arima.RData")
test_and_pred = datMxTest %>% mutate(pred = LR_arima[[1]]) %>% 
  filter(ymd >= filter_test)
DP_LR_arima_error = test_and_pred$load - test_and_pred$pred
DP_LR_arima = test_and_pred$pred

# FCNN
LRFCNN = as.matrix(read.csv("2. Pred Signals/2. NNs/1. FC/low-resolution_FC_rolling_pred.csv",header=F))
test_and_pred = datMxTest %>% mutate(pred = LRFCNN) %>% 
  filter(ymd >= filter_test)
DP_LRFCNN_error =  test_and_pred$load - test_and_pred$pred
DP_LRFCNN = test_and_pred$pred

# Gauss
load("2. Pred Signals/1. GAMs/gauss_rolling_pred_more_knots_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = gauss_rolling_pred[[1]]) %>% 
  filter(ymd >= filter_test)
DP_LR_gauss_error =  test_and_pred$load - test_and_pred$pred
DP_LR_gauss = test_and_pred$pred

# Scat
load("2. Pred Signals/1. GAMs/scat_rolling_pred_more_knots_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = scat_rolling_pred[[1]]) %>% 
  filter(ymd >= filter_test)
DP_LR_scat_error =  test_and_pred$load - test_and_pred$pred
DP_LR_scat = test_and_pred$pred

# Gev
load("2. Pred Signals/1. GAMs/gev_rolling_pred_more_knots_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = gev_rolling_pred[[1]]) %>% 
  filter(ymd >= filter_test)
DP_LR_gev_error =  test_and_pred$load - test_and_pred$pred
DP_LR_gev = test_and_pred$pred

##### Multi-resolution

# CNN
# MRCNN = apply(bind_cols(lapply(list.files(path = "2. Pred Signals/2. NNs/2. Hybrid/2. Without Trend/", pattern="*.csv", full.names=T), read.delim, header=F)),
#               1,
#               mean)

MRCNN=as.matrix(read.csv("2. Pred Signals/2. NNs/2. Hybrid/multi-resolution_CNN_rolling_pred_1.csv",header=F))
test_and_pred = datMxTest %>% mutate(pred = MRCNN) %>% 
  filter(ymd >= filter_test)
DP_MRCNN_error =  test_and_pred$load - test_and_pred$pred
DP_MRCNN = test_and_pred$pred

# Gauss
load("2. Pred Signals/1. GAMs/MR_gauss_rolling_pred_more_knots_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = MR_rolling_pred[[1]]) %>% 
  filter(ymd >= filter_test)
DP_MR_gauss_error =  test_and_pred$load - test_and_pred$pred
DP_MR_gauss = test_and_pred$pred

# Scat
load("2. Pred Signals/1. GAMs/MR_scat_rolling_pred_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = MR_scat_rolling_pred[[1]]) %>% 
  filter(ymd >= filter_test)
DP_MR_scat_error =  test_and_pred$load - test_and_pred$pred
DP_MR_scat = test_and_pred$pred

# Gev
load("2. Pred Signals/1. GAMs/MR_rolling_pred_more_knots_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = MR_rolling_pred[[1]]) %>% 
  filter(ymd >= filter_test)
DP_MR_gev_error =  test_and_pred$load - test_and_pred$pred
DP_MR_gev = test_and_pred$pred


#### 2. Diebold-Mariano Tests ####

### MultDM package

all_DP = c("DP_HR_arima","DP_HR_gauss","DP_HRFCNN",
        "DP_LR_arima","DP_LR_gauss","DP_LR_scat","DP_LR_gev","DP_LRFCNN",
        "DP_MR_gauss","DP_MR_scat","DP_MR_gev","DP_MRCNN"
)

DP_significance = tibble(model=all_DP, DP_HR_arima=NA,DP_HR_gauss=NA,DP_HRFCNN=NA,DP_LR_arima=NA,DP_LR_gauss=NA,DP_LR_scat=NA,
                         DP_LR_gev=NA,DP_LRFCNN=NA,DP_MR_gauss=NA,DP_MR_scat=NA,DP_MR_gev=NA,DP_MRCNN=NA) %>%
               column_to_rownames(var="model")

## DP

for (i in 1:12){
  for (j in 1:12){
    if (j>i){
      DP_significance[i,j] = round(DM.test(eval(as.name(all_DP[i])),eval(as.name(all_DP[j])),(datMxTest %>% filter(ymd >= filter_test))$load ,h=1,
                                           loss.type="AE",c=FALSE,H1="same")$p.value,3)
      DP_significance[j,i] = DP_significance[i,j]
    }
  }
}

significance_formatter = formatter("span", 
                                    style = x ~ style(color = ifelse(x > 0.05, "red", 
                                                              "green")))

formattable(DP_significance %>%
              rownames_to_column(var='model'),list(DP_HR_arima = significance_formatter,	
                                                   DP_HR_gauss	= significance_formatter,
                                                   DP_HRFCNN	= significance_formatter,
                                                   DP_LR_arima	= significance_formatter,
                                                   DP_LR_gauss	= significance_formatter,
                                                   DP_LR_scat	= significance_formatter,
                                                   DP_LR_gev	= significance_formatter,
                                                   DP_LRFCNN = significance_formatter,	
                                                   DP_MR_gauss	= significance_formatter, 
                                                   DP_MR_scat = significance_formatter,	
                                                   DP_MR_gev = significance_formatter,
                                                   DP_MRCNN = significance_formatter
              ))


###### Metrics Boxplots ######

lastyear_DP = datMxTest %>% filter(ymd >= filter_test) %>% select(load) %>% as.matrix(.)
lastyear_IP = datMxTest %>% filter(ymd >= filter_test) %>% select(todFrom1) %>% as.matrix(.)

all_DP_gams = c("DP_HR_gauss",
                "DP_LR_scat",
                "DP_MR_scat")

all_DP_NNs = c("DP_HRFCNN",
               "DP_LRFCNN",
               "DP_MRCNN")

pdf(file = "./3. Plots/DP_GAM_AE_BOX.pdf", width = 11, height = 8)
metrics_box(all_DP_gams, change_unit=FALSE, loss=loss_mae,loss_name="Absolute Errors [MW]")
dev.off()

pdf(file = "./3. Plots/DP_GAM_SE_BOX.pdf", width = 11, height = 8)
metrics_box(all_DP_gams, change_unit=FALSE,loss=loss_mse,loss_name="Squared Errors [MW²]")
dev.off()

pdf(file = "./3. Plots/DP_NN_AE_BOX.pdf", width = 11, height = 8)
metrics_box(all_DP_NNs, change_unit=FALSE, loss=loss_mae,loss_name="Absolute Errors [MW]")
dev.off()

pdf(file = "./3. Plots/DP_NN_SE_BOX.pdf", width = 11, height = 8)
metrics_box(all_DP_NNs, change_unit=FALSE,loss=loss_mse,loss_name="Squared Errors [MW²]")
dev.off()

##### END