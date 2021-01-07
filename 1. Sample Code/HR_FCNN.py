# 0. Set Session

import numpy as np
np.random.seed(24)

##Set GPU number and tensorflow configuratoin\n",
#import os
#os.environ['CUDA_VISIBLE_DEVICE'] = '0'
import tensorflow as tf
#tf.set_random_seed(24)
import keras
#config = tf.ConfigProto()
#config.gpu_options.allow_growth=True
#sess = tf.Session(config=config)

#%matplotlib inline
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
import copy
import warnings
from datetime import datetime as dt
warnings.filterwarnings(action='once')

from keras.models import Sequential, Model
from keras.layers import Input, Dense, LeakyReLU, BatchNormalization, Dropout, Activation
from keras.optimizers import RMSprop, Adam, Nadam, Adagrad, Adamax
from keras.callbacks import EarlyStopping, ModelCheckpoint, TensorBoard, ReduceLROnPlateau
from keras.initializers import glorot_uniform
from keras.activations import relu, elu
from tensorflow.keras.models import load_model

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

def mape(y_true, y_pred): 
    y_true, y_pred = np.array(y_true), np.array(y_pred)
    return np.mean(np.abs((y_true - y_pred) / y_true)) * 100
def rmse(y_true, y_pred): 
    y_true, y_pred = np.array(y_true), np.array(y_pred)
    return np.sqrt(np.mean((y_true - y_pred)**2))

def fit_plot(y_true,y_pred):
    fig = plt.figure(1, figsize=(18, 10))
    plt.plot(y_true)
    plt.plot(y_pred)
    plt.title('Fitted vs Observed', fontsize=18)
    plt.ylabel('Load [MW]', fontsize=18)
    plt.xlabel('Index', fontsize=18)
    plt.legend(['Observed', 'Fitted'], loc='upper right', fontsize=18)
    plt.show()
    print("MAPE:",mape(y_true,y_pred))
    print("RMSE:",rmse(y_true,y_pred))

def prep(df):
    
    df_test = df[df['dateDay'] >= dt.date(dt.strptime('2012-07-01', '%Y-%m-%d'))]

    # Feed Forward Input
    x_ffw_train = copy.deepcopy(pd.concat([df[features],
                       dummies(df, 'dow'),
                       dummies(df, 'tod'),
                       df['dateDay']
                      ],
                      axis = 1))
    x_ffw_test = x_ffw_train[x_ffw_train['dateDay'] >= dt.date(dt.strptime('2012-07-01', '%Y-%m-%d'))]

    # Output
    y_train = copy.deepcopy(df['load'].values)
    y_train.shape = (y_train.shape[0],1)
    
    y_test = copy.deepcopy(df_test['load'].values)
    y_test.shape = (y_test.shape[0],1)
    
    return(x_ffw_train.drop("dateDay", axis = 1).values, y_train, x_ffw_test.drop("dateDay", axis = 1).values, y_test)

# 2. Data Preparation

rolling_train = pd.read_csv('1. Data/high-resolution_train.csv')
rolling_train['ymd'] = pd.to_datetime(rolling_train['ymd'])
rolling_train['dateDay'] = [dt.date(d) for d in rolling_train['ymd']]

rolling_test = pd.read_csv('1. Data/high-resolution_test.csv')
rolling_test['ymd'] = pd.to_datetime(rolling_test['ymd'])
rolling_test['dateDay'] = [dt.date(d) for d in rolling_test['ymd']]

features = ["toy",
            "temp",
            "temp95",
            "load24"
            ]
selection = features + ["tod","dow", "load", 'dateDay', 'ymd', 'ym']

all_data = prep(rolling_train)
training, holdout = (all_data[0],all_data[1]), (all_data[2],all_data[3]) 

print(training[0].shape)
print(training[1].shape)

print(holdout[0].shape)
print(holdout[1].shape)

## Scaling Data

x_scaler, y_scaler = MinMaxScaler(), MinMaxScaler()

ym_unique = rolling_train['ym'].unique()
len(ym_unique)

# Fitting only the scaler only on the first year of data
filtered_train = rolling_train[rolling_train['ym'].isin(ym_unique[:12])] 
index_train = len(filtered_train)

x_scaler.fit(training[0][:index_train])
y_scaler.fit(training[1][:index_train])

training_scaled = (x_scaler.transform(training[0]),
                   y_scaler.transform(training[1])
                  )

holdout_scaled = (x_scaler.transform(holdout[0]),
                  y_scaler.transform(holdout[1])
                  )


## Modelling

seed = 24

input_ffw = Input(shape=(training_scaled[0].shape[1],))

ffw = Dense(50,use_bias=True,activation=None,kernel_initializer=glorot_uniform(seed=seed))(input_ffw)
ffw = BatchNormalization()(ffw)
ffw = Dropout(0.1)(ffw)
ffw = Activation('relu')(ffw)

ffw = Dense(1,use_bias=True,activation=None,kernel_initializer=glorot_uniform(seed=seed))(ffw)
ffw = Activation('relu')(ffw)

model = Model(inputs=input_ffw, outputs=ffw)

model.summary()

# Select Optimizer
optimizer = Nadam(lr=0.001)

# Compile Model
model.compile(loss='mse', optimizer=optimizer)
model.save('./weights/high-resolution_FC_rolling')

# Rolling origin training
fitted_signal = []
pred_signal = []
for i in range(len(ym_unique)-12):
#for i in range(1):
    print("================ Iteration",i+1, "================")
    model = load_model('./weights/high-resolution_FC_rolling')
    filtered_train = rolling_train[rolling_train['ym'].isin(ym_unique[:i+12])] #starting to train with the full first year
    index_train = len(filtered_train)
    
    filtered_test = rolling_train[rolling_train['ym'] == ym_unique[i+12]]
    index_test = len(filtered_test)
    
    # train
    history = model.fit(x=training_scaled[0][:index_train], 
                  y=training_scaled[1][:index_train], 
                  batch_size=1024,
                  epochs=2000,
                  verbose=0
                 )

    # fitted signal
    fitted_signal.append(list(y_scaler.inverse_transform(model.predict(training_scaled[0][:index_train]))))

    # pred
    pred_signal.append(y_scaler.inverse_transform(model.predict(training_scaled[0][index_train:index_train+index_test])))
    
pred_signal = np.vstack(pred_signal)

fit_plot(holdout[1],pred_signal)


for i in range(48):
    path = '2. Pred Signals/2. NNs/1. FC/fitted/high-resolution_FC_rolling_' + str(i) + '.csv'    
    np.savetxt(path, np.vstack(fitted_signal[i]), delimiter=',')

np.savetxt('2. Pred Signals/2. NNs/1. FC/high-resolution_FC_rolling_pred.csv', pred_signal, delimiter=',')