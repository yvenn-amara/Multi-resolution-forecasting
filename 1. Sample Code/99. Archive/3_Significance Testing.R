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
  filter(ymd >= "2015-07-01")
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
  filter(ymd >= "2015-07-01")
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
  filter(ymd >= "2015-07-01")
DP_HRFCNN_error =  test_and_pred$load - test_and_pred$pred
DP_HRFCNN = test_and_pred$pred

##### Low resolution
# ARIMA
load("2. Pred Signals/0. ARIMA/lr_arima.RData")
test_and_pred = datMxTest %>% mutate(pred = LR_arima[[1]]) %>% 
  filter(ymd >= "2015-07-01")
DP_LR_arima_error = test_and_pred$load - test_and_pred$pred
DP_LR_arima = test_and_pred$pred

# FCNN
LRFCNN = as.matrix(read.csv("2. Pred Signals/2. NNs/1. FC/low-resolution_FC_rolling_pred.csv",header=F))
test_and_pred = datMxTest %>% mutate(pred = LRFCNN) %>% 
  filter(ymd >= "2015-07-01")
DP_LRFCNN_error =  test_and_pred$load - test_and_pred$pred
DP_LRFCNN = test_and_pred$pred

# Gauss
load("2. Pred Signals/1. GAMs/gauss_rolling_pred_more_knots_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = gauss_rolling_pred[[1]]) %>% 
  filter(ymd >= "2015-07-01")
DP_LR_gauss_error =  test_and_pred$load - test_and_pred$pred
DP_LR_gauss = test_and_pred$pred

# Scat
load("2. Pred Signals/1. GAMs/scat_rolling_pred_more_knots_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = scat_rolling_pred[[1]]) %>% 
  filter(ymd >= "2015-07-01")
DP_LR_scat_error =  test_and_pred$load - test_and_pred$pred
DP_LR_scat = test_and_pred$pred

# Gev
load("2. Pred Signals/1. GAMs/gev_rolling_pred_more_knots_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = gev_rolling_pred[[1]]) %>% 
  filter(ymd >= "2015-07-01")
DP_LR_gev_error =  test_and_pred$load - test_and_pred$pred
DP_LR_gev = test_and_pred$pred

##### Multi-resolution

# CNN
# MRCNN = apply(bind_cols(lapply(list.files(path = "2. Pred Signals/2. NNs/2. Hybrid/2. Without Trend/", pattern="*.csv", full.names=T), read.delim, header=F)),
#               1,
#               mean)

MRCNN=as.matrix(read.csv("2. Pred Signals/2. NNs/2. Hybrid/multi-resolution_CNN_rolling_pred_1.csv",header=F))
test_and_pred = datMxTest %>% mutate(pred = MRCNN) %>% 
  filter(ymd >= "2015-07-01")
DP_MRCNN_error =  test_and_pred$load - test_and_pred$pred
DP_MRCNN = test_and_pred$pred

# Gauss
load("2. Pred Signals/1. GAMs/MR_gauss_rolling_pred_more_knots_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = MR_rolling_pred[[1]]) %>% 
  filter(ymd >= "2015-07-01")
DP_MR_gauss_error =  test_and_pred$load - test_and_pred$pred
DP_MR_gauss = test_and_pred$pred

# Scat
load("2. Pred Signals/1. GAMs/MR_scat_rolling_pred_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = MR_scat_rolling_pred[[1]]) %>% 
  filter(ymd >= "2015-07-01")
DP_MR_scat_error =  test_and_pred$load - test_and_pred$pred
DP_MR_scat = test_and_pred$pred

# Gev
load("2. Pred Signals/1. GAMs/MR_rolling_pred_more_knots_factor.RData")
test_and_pred = datMxTest %>% mutate(pred = MR_rolling_pred[[1]]) %>% 
  filter(ymd >= "2015-07-01")
DP_MR_gev_error =  test_and_pred$load - test_and_pred$pred
DP_MR_gev = test_and_pred$pred


#### IP ####

#### High-resolution ####

## Gaussian
load("2. Pred Signals/1. GAMs/high_res_rolling_bam.RData")
test_and_pred = high_res_to_peak(DataTest, high_res_pred[[1]]) %>% 
  filter(ymd >= "2015-07-01")
IP_HR_gauss = test_and_pred$pred
IP_HR_gauss_error = test_and_pred$todFrom1 - test_and_pred$pred

# FCNN
IP_HRFCNN = as.vector(as.matrix(read.csv("2. Pred Signals/2. NNs/1. FC/high-resolution_FC_rolling_pred.csv",header=F)))
test_and_pred = high_res_to_peak(DataTest, IP_HRFCNN) %>% 
  filter(ymd >= "2015-07-01")
IP_HRFCNN = test_and_pred$pred
IP_HRFCNN_error = test_and_pred$todFrom1 - test_and_pred$pred

#### Low-resolution ####

# Ocat
load("2. Pred Signals/1. GAMs/low_resolution_ocat.RData")
test_and_pred = datMxTest %>%
  mutate(pred = lr_ocat_rolling_pred[[1]]) %>% 
  filter(ymd >= "2015-07-01")
IP_LR_ocat = test_and_pred$pred
IP_LR_ocat_error = test_and_pred$todFrom1 - test_and_pred$pred

# FCNN
IP_LRFCNN = to_classes(round(as.matrix(read.csv("2. Pred Signals/2. NNs/1. FC/low-resolution_FC_rolling_pred_Instant.csv",header=F))))
test_and_pred = datMxTest %>% mutate(pred = IP_LRFCNN) %>% 
  filter(ymd >= "2015-07-01")
IP_LRFCNN = test_and_pred$pred
IP_LRFCNN_error = test_and_pred$todFrom1 - test_and_pred$pred

#### Multi-resolution ####

# Ocat
load("2. Pred Signals/1. GAMs/ocat_rolling_pred_more_knots.RData")
test_and_pred = datMxTest %>%
  mutate(pred = ocat_rolling_pred[[1]]) %>% 
  filter(ymd >= "2015-07-01")
IP_MR_ocat = test_and_pred$pred
IP_MR_ocat_error = test_and_pred$todFrom1 - test_and_pred$pred

# MRCNN
IP_MRCNN = to_classes(round(as.matrix(read.csv("2. Pred Signals/2. NNs/2. Hybrid/multi-resolution_CNN_rolling_pred_Instant_1.csv",header=F))))
test_and_pred = datMxTest %>% mutate(pred = IP_MRCNN) %>% 
  filter(ymd >= "2015-07-01")
IP_MRCNN = test_and_pred$pred
IP_MRCNN_error = test_and_pred$todFrom1 - test_and_pred$pred


#### 3. Diebold-Mariano Tests ####

### MultDM package

all_DP = c("DP_HR_arima","DP_HR_gauss","DP_HRFCNN",
        "DP_LR_arima","DP_LR_gauss","DP_LR_scat","DP_LR_gev","DP_LRFCNN",
        "DP_MR_gauss","DP_MR_scat","DP_MR_gev","DP_MRCNN"
)

all_IP = c("IP_HR_gauss","IP_HRFCNN",
           "IP_LR_ocat","IP_LRFCNN",
           "IP_MR_ocat","IP_MRCNN"
)

# all_resid = c("HR_arima_error","HR_gauss_error","HRFCNN_error",
#         "LR_arima_error","LR_gauss_error","LR_scat_error","LR_gev_error","LRFCNN_error",
#         "MR_gauss_error","MR_scat_error","MR_gev_error","MRCNN_error"
# )

DP_significance = tibble(model=all_DP, DP_HR_arima=NA,DP_HR_gauss=NA,DP_HRFCNN=NA,DP_LR_arima=NA,DP_LR_gauss=NA,DP_LR_scat=NA,
                         DP_LR_gev=NA,DP_LRFCNN=NA,DP_MR_gauss=NA,DP_MR_scat=NA,DP_MR_gev=NA,DP_MRCNN=NA) %>%
               column_to_rownames(var="model")

IP_significance = tibble(model=all_IP, IP_HR_gauss=NA,IP_HRFCNN=NA,IP_LR_ocat=NA,IP_LRFCNN=NA,IP_MR_ocat=NA,IP_MRCNN=NA) %>%
  column_to_rownames(var="model")

## DP

for (i in 1:12){
  for (j in 1:12){
    if (j>i){
      DP_significance[i,j] = round(DM.test(eval(as.name(all_DP[i])),eval(as.name(all_DP[j])),(datMxTest %>% filter(ymd >= "2015-07-01"))$load ,h=1,
                                           loss.type="SE",c=FALSE,H1="same")$p.value,3)
      # significance[i,j] = dm.test(eval(as.name(all_resid[i])),eval(as.name(all_resid[j])))$p.value
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

## IP

for (i in 1:6){
  for (j in 1:6){
    if (j>i){
      IP_significance[i,j] = round(DM.test(eval(as.name(all_IP[i])),eval(as.name(all_IP[j])),(datMxTest  %>% filter(ymd >= "2015-07-01"))$todFrom1,h=1,
                                           loss.type="SE",c=FALSE,H1="same")$p.value,3)
      IP_significance[j,i] = IP_significance[i,j]
    }
  }
}

formattable(IP_significance %>%
              rownames_to_column(var='model'),list(IP_HR_gauss = significance_formatter,	
                                                   IP_HRFCNN	= significance_formatter,
                                                   IP_LR_ocat	= significance_formatter,
                                                   IP_LRFCNN	= significance_formatter,
                                                   IP_MR_ocat	= significance_formatter,
                                                   IP_MRCNN	= significance_formatter
              ))

##### Block Bootstrap approach ######


IP_significance = tibble(model=all_IP, IP_HR_gauss=NA,IP_HRFCNN=NA,IP_LR_ocat=NA,IP_LRFCNN=NA,IP_MR_ocat=NA,IP_MRCNN=NA) %>%
  column_to_rownames(var="model")


loss = function(resid){return(abs(resid))}
# loss = function(resid){return(as.numeric(abs(resid) > 0))}
# loss = function(resid){return(resid^2)}


for (i in 1:6){
  for (j in 1:6){
    if (j>i){
      
      mean_difference = mean( loss(IP_MRCNN_error) - loss(IP_LRFCNN_error), trim=trim )
      
      library(boot)
      out = tsboot(loss(eval(as.name(all_IP[i]))) - loss(eval(as.name(all_IP[j]))), statistic = function(x) mean(x, trim=trim), R = 500, l = round(n/24), sim = "fixed")
      
      mean_sample_centered = out$t-mean(out$t)
      upper_bound = quantile(mean_sample_centered,0.05)
      lower_bound = quantile(mean_sample_centered,0.95)
      
      IP_significance[i,j] = between(mean_difference,upper_bound,lower_bound)
      IP_significance[j,i] = IP_significance[i,j]
    }
  }
}





n = length(IP_MRCNN_error)

# Sign loss (but can be abs or square loss or anything else)
#loss = function(x) sign(x)


trim = 0

loss = function(resid){return(abs(resid))}
# loss = function(resid){return(as.numeric(abs(resid) > 0))}
# loss = function(resid){return(resid^2)}

# These are our autocorrelated differences in losses
plot( loss(IP_MRCNN_error) - loss(IP_LRFCNN_error), pch= 16 )

# First forecast is better, but we need uncertainty estimate
mean_difference = mean( loss(IP_MRCNN_error) - loss(IP_LRFCNN_error), trim=trim )

# Get uncertainty using block-bootstrap
library(boot)
out = tsboot(loss(IP_MRCNN_error) - loss(IP_LRFCNN_error), statistic = function(x) mean(x, trim=trim), R = 500, l = round(n/24), sim = "fixed")

# var = out$t
# var_centered = var - mean(var)

mean_sample_centered = out$t-mean(out$t)
upper_bound = quantile(mean_sample_centered,0.05)
lower_bound = quantile(mean_sample_centered,0.95)
between(mean_difference,upper_bound,lower_bound)

#var_centered_sign = -var_centered
 
hist(mean_sample_centered)
abline(v = mean_difference, col = 2, lty = 2)
abline(v = quantile(mean_sample_centered,0.05), col = 3, lty = 2)
abline(v = quantile(mean_sample_centered,0.95), col = 3, lty = 2)

between(mean_difference )
quantile(var_centered,0.05)
quantile(var_centered,0.95)




######

###### Metrics Boxplots ######

lastyear_DP = datMxTest %>% filter(ymd >= "2015-07-01") %>% select(load) %>% as.matrix(.)
lastyear_IP = datMxTest %>% filter(ymd >= "2015-07-01") %>% select(todFrom1) %>% as.matrix(.)

all_DP_gams = c("DP_HR_gauss",
                "DP_LR_scat",
                "DP_MR_scat")

all_DP_NNs = c("DP_HRFCNN",
               "DP_LRFCNN",
               "DP_MRCNN")

all_IP_gams = c("IP_HR_gauss",
                "IP_LR_ocat",
                "IP_MR_ocat")

all_IP_NNs = c("IP_HRFCNN",
               "IP_LRFCNN",
               "IP_MRCNN")

metrics_box = function(model_names,loss,loss_name,output="DP",change_unit=FALSE){
  losses = NULL
  
  if(output == 'DP'){
    for(i in 1:length(model_names)){
      losses = rbind(losses,
                     tibble(Model=rep(model_names[i],length(lastyear_DP)),
                            loss_mae=c(loss(lastyear_DP,eval(as.name(model_names[i])),change_unit))  ))
    }
  }
  else if(output == 'IP'){
    for(i in 1:length(model_names)){
      losses = rbind(losses,
                     tibble(Model=rep(model_names[i],length(lastyear_IP)),
                            loss_mae=c(loss(lastyear_IP,eval(as.name(model_names[i])),change_unit))  ))
    }
  }
  

  
  ggplot(data = losses, aes(x=Model, y=loss_mae)) + 
    geom_boxplot(aes(fill=Model)) +
    theme_Publication() +
    theme(legend.position='none') +
    ylab(loss_name)
  
}


histo_box = function(model_names,output="DP"){
  losses = NULL
  
  if(output == 'DP'){
    for(i in 1:length(model_names)){
      losses = rbind(losses,
                     tibble(model=rep(model_names[i],length(lastyear_DP)),
                            loss_mae=c(loss_mae(lastyear_DP,eval(as.name(model_names[i]))))  ))
    }
  }
  else if(output == 'IP'){
    for(i in 1:length(model_names)){
      losses = rbind(losses,
                     tibble(model=rep(model_names[i],length(lastyear_IP)),
                            loss_mae=c(loss_mae(lastyear_IP,eval(as.name(model_names[i]))))  ))
    }
  }
  
  
  
  ggplot(data = losses, aes(x=loss_mae)) + 
    geom_histogram(aes(fill=model)) +
    theme_Publication() +
    theme(legend.position='bottom')
  
}


metrics_box(all_DP_gams, change_unit=FALSE, loss=loss_mae,loss_name="Absolute Errors [MW]")
metrics_box(all_DP_gams, change_unit=FALSE,loss=loss_mse,loss_name="Squared Errors [MW²]")
metrics_box(all_DP_gams, change_unit=FALSE,loss=loss_mape,loss_name="MAPE [%]")


metrics_box(all_DP_NNs, change_unit=FALSE, loss=loss_mae,loss_name="Absolute Errors [MW]")
metrics_box(all_DP_NNs, change_unit=FALSE,loss=loss_mse,loss_name="Squared Errors [MW²]")
metrics_box(all_DP_NNs, change_unit=FALSE,loss=loss_mape,loss_name="MAPE [%]")


metrics_box(all_IP_gams,change_unit=FALSE, loss=loss_mae,loss_name="Absolute Errors [half-hour]",output="IP")
metrics_box(all_IP_gams,change_unit=FALSE, loss=loss_mse,loss_name="Squared Errors [half-hour²]",output="IP")

metrics_box(all_IP_NNs,change_unit=FALSE, loss=loss_mae,loss_name="Absolute Errors [half-hour]",output="IP")
metrics_box(all_IP_NNs,change_unit=FALSE, loss=loss_mse,loss_name="Squared Errors [half-hour²]",output="IP")


histo_box(all_DP_gams)
histo_box(all_DP_NNs)
histo_box(all_IP_gams,output="IP")
histo_box(all_IP_NNs,output="IP")



loss_table = function(model_names,output="DP"){
  losses = NULL
  
  if(output == 'DP'){
    for(i in 1:length(model_names)){
      losses = rbind(losses,
                     tibble(model=rep(model_names[i],length(lastyear_DP)),
                            loss_mae=c(loss_mae(lastyear_DP,eval(as.name(model_names[i]))))  ))
    }
  }
  else if(output == 'IP'){
    for(i in 1:length(model_names)){
      losses = rbind(losses,
                     tibble(model=rep(model_names[i],length(lastyear_IP)),
                            loss_mae=c(loss_mae(lastyear_IP,eval(as.name(model_names[i]))))  ))
    }
  }
  return(losses)
}


histo_box_pairs = function(loss_table,model_names,output="DP"){
  require(gridExtra)
  if(output=="IP"){
    plots=list()
    counter = 0
    for (i in 1:6){
      for (j in 1:6){
        if(i > j){
          counter= counter +1
          first = (loss_table %>% filter(model==model_names[i]))$loss_mae
          second = (loss_table %>% filter(model==model_names[j]))$loss_mae
          for_plot = tibble(difference=first-second)
          
          plots[[counter]] = ggplot(data=for_plot,aes(x=difference)) + 
            geom_histogram() +
            theme_Publication() +
            theme(legend.position='bottom') +
            ggtitle(paste(model_names[i],"-",model_names[j]))
        }
      }

    }
    grid.arrange(plots[[1]], plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],plots[[7]],plots[[8]],
                 plots[[9]],plots[[10]],plots[[11]],plots[[12]],plots[[13]],plots[[14]],plots[[15]],ncol=4)

  }
  
  if(output=="DP"){
    plots=list()
    counter = 0
    for (i in 1:6){
      for (j in 1:6){
        if(i > j){
          counter= counter +1
          first = (loss_table %>% filter(model==model_names[i]))$loss_mae
          second = (loss_table %>% filter(model==model_names[j]))$loss_mae
          for_plot = tibble(difference=first-second)
          
          plots[[counter]] = ggplot(data=for_plot,aes(x=difference)) + 
            geom_histogram() +
            theme_Publication() +
            theme(legend.position='bottom') +
            ggtitle(paste(model_names[i],"-",model_names[j]))
        }
      }
      
      grid.arrange(plots[[1]], plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],plots[[7]],plots[[8]],
                   plots[[9]],plots[[10]],plots[[11]],plots[[12]],plots[[13]],plots[[14]],plots[[15]],ncol=4)
      
      
    }

    
  }
  
}

a = loss_table(all_IP,output="IP")
histo_box_pairs(a,output="IP")


all_DP_relevant = c("DP_HRFCNN",
                    "DP_LR_gauss","DP_LR_scat","DP_LRFCNN",
                    "DP_MR_scat","DP_MRCNN"
)

a = loss_table(all_DP_relevant,output="DP")
histo_box_pairs(a,all_DP_relevant,output="IP")
