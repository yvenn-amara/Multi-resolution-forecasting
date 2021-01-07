# 0. Set Session

import numpy as np
import keras
import tensorflow as tf
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
import copy
from datetime import datetime as dt
from keras.models import Sequential, Model
from keras.layers import Input, Dense, LeakyReLU, BatchNormalization, Dropout, Activation
from keras.optimizers import RMSprop, Adam, Nadam, Adagrad, Adamax
from keras.callbacks import EarlyStopping, ModelCheckpoint, TensorBoard, ReduceLROnPlateau
from keras.initializers import glorot_uniform
from keras.activations import relu, elu
from tensorflow.keras.models import load_model

np.random.seed(24)

# 1. Functions

#Get dummies
def dummies(df, colname):
    
    """ Takes a dataframe and transforms one of its column into multiple dummy variables (0-1).
    It returns a dataframe containing only the dummy variables creating, thus it requires 
    a concatenation with the orginial data to be used for training a model """
    
    X = pd.get_dummies(df[colname],
                         prefix = colname, 
                         prefix_sep = '_', 
                         drop_first=True)
    return(X)

def prep(df):
    
    df_test = df[df['dateDay'] >= dt.date(dt.strptime('2012-07-01', '%Y-%m-%d'))]

    # Feed Forward Input
    x_ffw_train = copy.deepcopy(pd.concat([df[features],
                       dummies(df, 'dow'),
                       dummies(df, 'tod24'),
                       df['dateDay']
                      ],
                      axis = 1))
    x_ffw_test = x_ffw_train[x_ffw_train['dateDay'] >= dt.date(dt.strptime('2012-07-01', '%Y-%m-%d'))]

    # Output
    y_train = copy.deepcopy(df[output_cols].values)
    
    y_test = copy.deepcopy(df_test[output_cols].values)
    
    return(x_ffw_train.drop("dateDay", axis = 1).values, y_train, x_ffw_test.drop("dateDay", axis = 1).values, y_test)

# 2. Data Preparation

rolling_train = pd.read_csv('1. Data/ordinal_learn.csv')
rolling_train['ymd'] = pd.to_datetime(rolling_train['ymd'])
rolling_train['dateDay'] = [dt.date(d) for d in rolling_train['ymd']]

rolling_test = pd.read_csv('1. Data/ordinal_test.csv')
rolling_test['ymd'] = pd.to_datetime(rolling_test['ymd'])
rolling_test['dateDay'] = [dt.date(d) for d in rolling_test['ymd']]

output_cols = [col for col in rolling_train if col.startswith('V')]

features = ["toy", 
            #"trend",
            "peak24",
            "tempMax",
            "temp95Max",
            "tempMin",
            "temp95Min"
            ]
selection = features + ["tod24","dow", 'dateDay', 'ymd', 'ym'] + output_cols

rolling_train = rolling_train[selection].dropna()
rolling_test = rolling_test[selection].dropna()

all_data = prep(rolling_train)
training, holdout = (all_data[0],all_data[1]), (all_data[2],all_data[3]) 

print(training[0].shape)
print(training[1].shape)

print(holdout[0].shape)
print(holdout[1].shape)

x_scaler = MinMaxScaler()

ym_unique = rolling_train['ym'].unique()
len(ym_unique)

# Fitting only the scaler on the first year of data
filtered_train = rolling_train[rolling_train['ym'].isin(ym_unique[:12])] 
index_train = len(filtered_train)

x_scaler.fit(training[0][:index_train])

training_scaled = (x_scaler.transform(training[0]),
                   training[1]
                  )

holdout_scaled = (x_scaler.transform(holdout[0]),
                  holdout[1]
                  )

# 3. Modelling

seed = 24

input_ffw = Input(shape=(training_scaled[0].shape[1],))

ffw = Dense(50,use_bias=True,activation=None,kernel_initializer=glorot_uniform(seed=seed))(input_ffw)
ffw = BatchNormalization()(ffw)
ffw = Dropout(0.1)(ffw)
ffw = Activation('relu')(ffw)
ffw = Dense(48,use_bias=True,activation=None,kernel_initializer=glorot_uniform(seed=seed))(ffw)
ffw = Activation('sigmoid')(ffw)

model = Model(inputs=input_ffw, outputs=ffw)
model.summary()

# Select Optimizer
optimizer = Nadam(lr=0.001)

# Compile Model
model.compile(loss='mse', optimizer=optimizer)
model.save('./weights/low-resolution_FC_rolling_Instant')

# Rolling origin training

pred_signal = []
fitted_signal = []
for i in range(len(ym_unique)-12):
    print("================ Iteration",i+1, "================")
    model = load_model('./weights/low-resolution_FC_rolling_Instant')
    filtered_train = rolling_train[rolling_train['ym'].isin(ym_unique[:i+12])] #starting to train with the full first year
    index_train = len(filtered_train)
    
    filtered_test = rolling_train[rolling_train['ym'] == ym_unique[i+12]]
    index_test = len(filtered_test)
    
    # train
    history = model.fit(x=training_scaled[0][:index_train+1], 
                  y=training_scaled[1][:index_train+1], 
                  batch_size=1024,
                  epochs=2000, 
                  verbose=0
                 )

    # fitted signal
    fitted_signal.append(list(model.predict(training_scaled[0][:index_train])))

    # pred
    pred_signal.append(model.predict(training_scaled[0][index_train:index_train+index_test]))


pred_signal = np.vstack(pred_signal)

for i in range(48):
    path = '2. Pred Signals/2. NNs/1. FC/fitted/low-resolution_FC_rolling_Instant_' + str(i) + '.csv'    
    np.savetxt(path, np.vstack(fitted_signal[i]), delimiter=',')

np.savetxt('2. Pred Signals/2. NNs/1. FC/low-resolution_FC_rolling_pred_Instant.csv', pred_signal, delimiter=',')