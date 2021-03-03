library(tidyverse)
`%notin%` <- Negate(`%in%`)

#### 0. Data Preparation ####

prep = function(){
Data = read.csv("1. Data/UK.csv") %>%
  select(date, Annee, LOAD, wM, wM_s95, Posan, Tendance, JourSemaine, Instant, LOAD.48) %>%
  rename("year"="Annee",
         "load"="LOAD",
         "temp"="wM",
         "temp95"="wM_s95",
         "toy"="Posan",
         "timeCount"="Tendance",
         "dow"="JourSemaine",
         "tod"="Instant",
         "load24"="LOAD.48") %>%
  na.omit(.)

# Discarding days where some of the 48 infra-day time steps are missing
Data$ymd <- as.Date(Data$date)
good <- which( table(Data$ymd) == 48 )
Data <- Data[Data$ymd %in% as.Date(names(table(Data$ymd))[good]), ]

# Keeping only full years starting from 1st July 2011
Data = Data %>%
  filter(ymd >= as.Date('2011-07-01'))

# Adding the year-month field for the rolling-origin training
Data = Data %>%
  mutate(ym = paste(year(ymd),"-",month(ymd), sep=""))

# Create covariates for smaller (daily) dataset
datMx = Data %>% group_by(ymd) %>% 
  summarize(ym = first(ym),
            dow = first(dow), 
            toy = first(toy), 
            timeCount = first(timeCount), 
            tod = tod[which.max(load)], 
            load = max(load), 
            tempMax = max(temp), 
            temp95Max = max(temp95), 
            tempMin = min(temp), 
            temp95Min = min(temp95)) %>%
  mutate(dow = as.factor(dow))

datMx <- datMx %>% 
  mutate(todFrom1 = tod + 1, 
         peak24 = lag(load), 
         tod24 = lag(tod))

# Create daily training and testing sets
datMxLearn <- datMx
datMxTest <- datMx %>% filter(ymd >= as.Date("2012-07-01"))

# Create by-instant training and testing sets
DataLearn <- Data
DataTest <- Data %>% filter(ymd >= as.Date("2012-07-01"))

### Create matrices of temperature, smooth temp, instant and lagged demand for functional effects
# Training set
datMxLearn$matTem <- matrix(DataLearn$temp, nrow(datMxLearn), 48, byrow = TRUE)
datMxLearn$matTem95 <- matrix(DataLearn$temp95, nrow(datMxLearn), 48, byrow = TRUE)
datMxLearn$matInt <- matrix(DataLearn$tod, nrow(datMxLearn), 48, byrow = TRUE)
datMxLearn$matLag <- matrix(DataLearn$load24, nrow(datMxLearn), 48, byrow = TRUE)

# Testing set
datMxTest$matTem <- matrix(DataTest$temp, nrow(datMxTest), 48, byrow = TRUE)
datMxTest$matTem95 <- matrix(DataTest$temp95, nrow(datMxTest), 48, byrow = TRUE)
datMxTest$matInt <- matrix(DataTest$tod, nrow(datMxTest), 48, byrow = TRUE)
datMxTest$matLag <- matrix(DataTest$load24, nrow(datMxTest), 48, byrow = TRUE)


## Datasets saved for training
# write.csv(datMxLearn, file="1. Data/low-resolution_train.csv", row.names = FALSE) # low-resolution and multi-resolution
# write.csv(datMxTest, file="1. Data/low-resolution_test.csv", row.names = FALSE) # low-resolution and multi-resolution
# write.csv(DataLearn, file="1. Data/high-resolution_train.csv", row.names = FALSE) # high-resolution
# write.csv(DataTest, file="1. Data/high-resolution_test.csv", row.names = FALSE) # high-resolution

## Save for NNrank model
# ordinal_learn = to_ordinal(datMxLearn)
# ordinal_test = to_ordinal(datMxTest)
# write.csv(ordinal_learn, file="1. Data/ordinal_learn.csv", row.names = FALSE)
# write.csv(ordinal_test, file="1. Data/ordinal_test.csv", row.names = FALSE)

return(list(datMxLearn,datMxTest,DataLearn,DataTest))

}

#### 1. Training Framework ####

arima_rolling = function(df){
  ym_unique = unique(df$ym)
  
  pred_signal = NULL
  models = list()

  for(i in 12:(length(ym_unique)-1)){
    print(paste("Iteration",i-11))
    
    load_train = df %>%
      filter(ym %in% ym_unique[1:i]) %>%
      select(load) %>%
      as.matrix(.)
    
    model = auto.arima(ts(load_train,frequency=1))
    models[[i-11]] = model
    
    refit = Arima(df %>%
                    filter(ym %in% ym_unique[1:(i+1)]) %>%
                    select(load) %>%
                    as.matrix(.),
                  model=model)
    
    pred_signal = c(pred_signal,refit$fitted[(nrow(load_train)+1):nrow(refit$fitted)]) 
    
  }
  
  return(list(pred_signal,models))
  
}


rolling_origin = function(df, model_function, family = 'gauss'){
  
  ym_unique = unique(df$ym)
  
  pred_signal = NULL
  models = list()
  
  
    for(i in 12:(length(ym_unique)-1)){
      
      print(paste("Iteration",i-11))
      
      filtered_df = df %>%
        filter(ym %in% ym_unique[1:i])
      
      model = model_function(filtered_df)
      
      models[[i-11]] = model
      
      if (family == 'gevlss'){
        
        
        lss = predict(model, df %>% filter(ym == ym_unique[i+1]))
        
        
        mu = as.vector(lss[,1])
        sigma = exp(as.vector(lss[,2]))
        xi = as.vector(lss[,3])
        gamma = -digamma(1)
        
        
        if (xi[1] == 0){
          pred = mu + sigma*gamma
        }
        
        else if (xi[1] < 1){
          pred = mu + sigma*(gamma(1-xi)-1)/xi
        }
        
        else if (xi[1] >= 1){
          pred = NA
        }
        pred_signal = c(pred_signal,pred)
      }
      
      else if (family == 'elf'){
        
        temp = NULL
        
        for (j in 1:length(model)){
          temp = bind_cols(temp, tibble(j=predict(model[[j]], df %>% filter(ym == ym_unique[i+1])))) 
        }
        
        pred_signal = bind_rows(pred_signal, temp)
        
      }
      
      else if (family == 'ocat'){
        
        pred = predict(model, df %>% filter(ym == ym_unique[i+1]), type="response")
        
        temp = NULL
        for (k in 1:nrow(pred)){
          temp = c(temp, which.max(pred[k,]))
        }
        
        pred_signal = c(pred_signal, temp)
        
        
      }
      
      else{
        pred_signal = c(pred_signal,predict(model, df %>% filter(ym == ym_unique[i+1]))) 
      }
      
      }

  return(list(pred_signal,models))
}

#### 2. Ordinal output for NNrank model (IP)  ####

## Classes to ordinal (1 ...1,0,...0)
to_ordinal = function(df){
  ordinal_output = NULL
  for (i in 1:nrow(df)){
    ordinal_output = rbind(ordinal_output,c(rep(1,df$todFrom1[i]),rep(0,48-df$todFrom1[i])))
  }
  
  return(df %>% cbind(.,as.data.frame(ordinal_output)))
}

## Ordinal (1 ...1,0,...0) to classes
to_classes = function(df){
  classes_output = NULL
  for (i in 1:nrow(df)){
    classes_output = c(classes_output,sum(df[i,]))
  }
  
  return(classes_output)
  
}

#### 3. Predictions ####

# Calculates the expected value of a GEV random variable from the estimated location, scale and shape
gev_pred = function(lss){
  
  mu = as.vector(lss[,1])
  sigma = exp(as.vector(lss[,2]))
  xi = as.vector(lss[,3])
  gamma = -digamma(1)
  
  
  if (xi[1] == 0){
    pred = mu + sigma*gamma
  }
  
  else if (xi[1] < 1){
    pred = mu + sigma*(gamma(1-xi)-1)/xi
  }
  
  else if (xi[1] >= 1){
    pred = NA
  }
  
  return(pred)
}

# Estimated DP from the high_res NN
high_res_to_max=function(DataLearn,high_res_signal){
  
  
  pred = DataLearn[1:length(high_res_signal),] %>%
    mutate(pred = high_res_signal) %>%
    select(ymd,pred) %>%
    group_by(ymd) %>%
    summarise_all(max) %>%
    select(pred)%>%
    as.matrix(.)
  return(pred)
  
}

# Estimated IP from the high_res NN
high_res_to_peak = function(DataTest,high_res_pred){
  demand_pred = DataTest %>%
    mutate(pred = high_res_pred) %>%
    select(ymd,tod,pred)
  
  peak_pred = DataTest %>%
    mutate(pred = high_res_pred) %>%
    select(ymd,pred) %>%
    group_by(ymd) %>%
    summarise_all(max) %>%
    ungroup(.)
  
  test_and_pred = inner_join(demand_pred,peak_pred,by=c("pred","ymd")) %>%
    select(ymd,tod) %>%
    rename("pred"="tod") %>%
    inner_join(.,datMxTest %>% select(ymd,todFrom1),by='ymd') %>%
    mutate(ym = paste(year(ymd),"-",month(ymd), sep=""))
  
  return(test_and_pred)
}

#### 4. Performance Metrics ####

loss_mape = function(y_true,y_pred){
  return(100*(abs((y_true-y_pred)/y_true)))
}

mape = function(y_true,y_pred){
  return(100*mean(abs((y_true-y_pred)/y_true)))
}

sd_mape = function(y_true,y_pred){
  l_bar = mape(y_true,y_pred)
  return(sqrt(mean((100*(abs((y_true-y_pred)/y_true))-l_bar)^2)))
}

mae = function(y_true,y_pred){
  return(mean(abs(y_true-y_pred)))
}

sd_mae = function(y_true,y_pred){
  l_bar = mae(y_true,y_pred)
  return(sqrt(mean((abs(y_true-y_pred)-l_bar)^2)))
}

rmse = function(y_true,y_pred){
  return(sqrt(mean((y_true-y_pred)^2)))
}



acc = function(y_true,y_pred){
  temp = (y_true == y_pred)
  return(100*sum(temp)/length(temp))
}

AIC_metric_GAM = function(mod_list){
  
  df_AIC = tibble(window = 1:48, AIC=0)
  
  for(i in 1:48){
    observed = mod_list[[2]][[i]]$y
    fitted = mod_list[[2]][[i]]$fitted.values
    parameters = length(mod_list[[2]][[i]]$coefficients)
    df_AIC$AIC[i]=AIC(mod_list[[2]][[i]])
  }
  return(df_AIC)
}

AIC_metric_GAM_instant = function(mod_list){
  
  df_AIC = tibble(window = 1:48, AIC=0)
  
  for(i in 1:48){
    df_AIC$AIC[i]=AIC(mod_list[[2]][[i]])
  }
  return(df_AIC)
}

metrics_over_time = function(df,model,res){
  ym_unique = unique(df$ym)
  
  metrics = tibble(ym = ym_unique, ymd = parse_date_time(ym_unique, orders="%Y-%m"),mape=0,rmse=0,mae=0,
                                                                                    rolling_mape=0,rolling_rmse=0,rolling_mae=0)
  
  for(i in 1:length(ym_unique)){
    
    filtered = df %>%
      filter(ym == ym_unique[i])
    
    metrics$mape[i]=mape(filtered$load,filtered$pred)
    metrics$rmse[i]=rmse(filtered$load,filtered$pred)
    metrics$mae[i]=mae(filtered$load,filtered$pred)
    
    rolling_filtered = df %>%
      filter(ym %in% ym_unique[1:i])      
    
    metrics$rolling_mape[i]=mape(rolling_filtered$load,rolling_filtered$pred)
    metrics$rolling_rmse[i]=rmse(rolling_filtered$load,rolling_filtered$pred)
    metrics$rolling_mae[i]=mae(rolling_filtered$load,rolling_filtered$pred)
    
  }
  
  metrics$Model = model
  metrics$Resolution = res
  
  return(metrics)
}

metrics_over_time_instant = function(df,model,res){
  ym_unique = unique(df$ym)
  
  metrics = tibble(ym = ym_unique, ymd = parse_date_time(ym_unique, orders="%Y-%m"),mape=0,rmse=0,acc=0,mae=0,
                                                                                    rolling_mape=0,rolling_rmse=0,rolling_acc=0,rolling_mae=0)
  
  for(i in 1:length(ym_unique)){
    
    filtered = df %>%
      filter(ym == ym_unique[i])
    
    metrics$mape[i]=mape(filtered$todFrom1,filtered$pred)
    metrics$rmse[i]=rmse(filtered$todFrom1,filtered$pred)
    metrics$acc[i]=acc(filtered$todFrom1,filtered$pred)
    metrics$mae[i]=mae(filtered$todFrom1,filtered$pred)
    
    rolling_filtered = df %>%
      filter(ym %in% ym_unique[1:i])      
    
    metrics$rolling_mape[i]=mape(rolling_filtered$todFrom1,rolling_filtered$pred)
    metrics$rolling_rmse[i]=rmse(rolling_filtered$todFrom1,rolling_filtered$pred)
    metrics$rolling_acc[i]=acc(rolling_filtered$todFrom1,rolling_filtered$pred)
    metrics$rolling_mae[i]=mae(rolling_filtered$todFrom1,rolling_filtered$pred)
    
  }
  
  metrics$Model = model
  metrics$Resolution = res
  
  return(metrics)
}
 
last_year = function(df){
  df = df %>%
    filter(ymd >= "2015-07-01")
  
  print(paste("MAPE:",mape(df$load,df$pred), "sd:", sd_mape(df$load,df$pred)))
  print(paste("MAE:",mae(df$load,df$pred), "sd:", sd_mae(df$load,df$pred)))
  print(paste("RMSE:",rmse(df$load,df$pred)))
}

last_year_sd = function(df){
  df = df %>%
    filter(ymd >= "2015-07-01")
  
  print(paste("MAPE-sd:",sd(df$mape)))
  
  
  
  print(paste("MAE-sd:",sd(df$mae)))
  print(paste("RMSE-sd:",sd(df$rmse)))
}

last_year_instant = function(df){
  df = df %>%
    filter(ymd >= "2015-07-01")
  
  print(paste("Accuracy:",acc(df$todFrom1,df$pred)))
  print(paste("MAPE:",mape(df$todFrom1,df$pred)))
  print(paste("RMSE:",rmse(df$todFrom1,df$pred)))
  print(paste("MAE:",mae(df$todFrom1,df$pred)))
}

last_year_instant_sd = function(df){
  df = df %>%
    filter(ymd >= "2015-07-01")
  
  print(paste("Accuracy:",sd(df$acc)))
  print(paste("MAPE:",sd(df$mape)))
  print(paste("MAE:",sd(df$mae)))
  print(paste("RMSE:",sd(df$rmse)))
}

##### 5. Plots 

plot_metrics=function(dataset, metric){
  
  cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000","#994F00")
  
  if (metric=="mape") {
    ggplot(dataset, aes(ymd, group=Model, colour=Model, shape=Model)) + 
      geom_line(aes(y=rolling_mape, linetype = Model)) +
      geom_point(aes(y=rolling_mape),size=2) +
      xlab('Month') +
      ylab('Rolling MAPE') +
      theme_Publication() +
      # theme(legend.title=element_text(size=30),
      #       legend.position = c(.95, .95),
      #       legend.justification = c("right", "top"),
      #       legend.box.just = "right",
      #       legend.margin = margin(6, 6, 6, 6),
      #       legend.text=element_text(size=24),
      #       axis.title = element_text(size=24),
      #       axis.text = element_text(size=24),
      #       axis.ticks = element_line(colour = "black", size = 0.2),
      #       panel.grid.major = element_line(colour = "black", size = 0.2,linetype="dashed"),
      #       panel.grid.minor = element_blank()) +
      # scale_colour_manual(values=cbPalette) +
      scale_shape_manual(values=c(15,15,15,16,16,16,17,17,18)) +
      scale_y_continuous(breaks = round(seq(0, max(dataset$rolling_mape), by = 0.5),1)) +
      scale_linetype_manual(values=c("solid","longdash","dashed","solid","longdash","dashed","solid","longdash","solid"))
    
    
  } else if (metric=="rmse"){
    ggplot(dataset, aes(ymd, group=Model, colour=Model)) + 
      geom_line(aes(y=rolling_rmse, linetype = Model), size = 2.5) +
      geom_point(aes(y=rolling_rmse), size = 3) +
      xlab('Month') +
      ylab('Rolling RMSE') +
      theme_minimal() +
      theme(legend.title=element_text(size=30),
            legend.text=element_text(size=24),
            axis.title = element_text(size=24),
            axis.text = element_text(size=24)) +
      scale_colour_manual(values=cbPalette) #+
      #scale_shape_manual(values=seq(0,15))
  } 
  # else if (metric=="mis-AIC"){
  #   ggplot(dataset, aes(window, group=Model, colour=Resolution)) + 
  #     geom_line(aes(y=mis_AIC), size = 1.5) +
  #     geom_point(aes(y=mis_AIC), size = 3) +
  #     xlab('Window') +
  #     ylab('mis-AIC') +
  #     theme_minimal() +
  #     theme(legend.title=element_text(size=30),
  #           legend.text=element_text(size=24),
  #           axis.title = element_text(size=24),
  #           axis.text = element_text(size=24)) +
  #     scale_colour_manual(values=cbPalette) #+
  #     #scale_shape_manual(values=c(4,8,15,16,17,18,21,22,3,42))
  # }
  
  else if (metric=="AIC"){
    ggplot(dataset, aes(ymd, group=Model, colour=Model)) + 
      geom_line(aes(y=AIC), size = 1.5) +
      geom_point(aes(y=AIC), size = 3) +
      xlab('Month') +
      ylab('AIC') +
      theme_minimal() +
      theme(legend.title=element_text(size=30),
            legend.text=element_text(size=24),
            axis.title = element_text(size=24),
            axis.text = element_text(size=24)) +
      scale_colour_manual(values=cbPalette) #+
      #scale_shape_manual(values=c(4,8,15,16,17,18,21,22,3,42))
  }
  
  else if (metric=="acc"){
    ggplot(dataset, aes(ymd, group=Model, colour=Model)) + 
      geom_line(aes(y=rolling_acc), size = 1.5) +
      geom_point(aes(y=rolling_acc), size = 3) +
      xlab('Month') +
      ylab('Rolling Accuracy') +
      theme_minimal() +
      theme(legend.title=element_text(size=30),
            legend.text=element_text(size=24),
            axis.title = element_text(size=24),
            axis.text = element_text(size=24)) +
      scale_colour_manual(values=cbPalette) #+
      #scale_shape_manual(values=c(4,8,15,16,17,18,21,22,3,42))
  }
  
  else if (metric=="mae"){
    ggplot(dataset, aes(ymd, group=Model, colour=Model)) + 
      geom_line(aes(y=rolling_mae), size = 1.5) +
      geom_point(aes(y=rolling_mae), size = 3) +
      xlab('Month') +
      ylab('Rolling MAE') +
      theme_minimal() +
      theme(legend.title=element_text(size=30),
            legend.text=element_text(size=24),
            axis.title = element_text(size=24),
            axis.text = element_text(size=24)) +
      scale_colour_manual(values=cbPalette) #+
      #scale_shape_manual(values=c(4,8,15,16,17,18,21,22,3,42))
  }
  
  
  
}

#### Theme inspired from https://github.com/koundy/ggplot_theme_Publication/blob/master/R/ggplot_theme_Publication.R#L7
theme_Publication = function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            #panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = c(.95, .95),
            legend.justification = c("right", "top"),
            legend.box.just = "right",
            #legend.position = "bottom",
            #legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            legend.key.size= unit(1, "cm"),
            legend.margin = margin(6, 6, 6, 6),
            legend.title = element_text(face="italic"),
            panel.grid.major = element_line(colour = "black", size = 0.2,linetype="dashed"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


######### END