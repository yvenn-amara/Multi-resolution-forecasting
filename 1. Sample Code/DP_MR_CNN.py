
# 0. Set Session

import numpy as np
import tensorflow as tf
import keras
import matplotlib.pyplot as plt
import pandas as pd
import os
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
import copy
from datetime import datetime as dt
from keras.models import Sequential, Model
from keras.layers import Input, Activation, Dense, Conv1D, Flatten, Dropout, LeakyReLU, BatchNormalization, concatenate, MaxPooling1D
from keras.optimizers import RMSprop, Adam, Nadam, Adagrad
from keras.callbacks import EarlyStopping, ModelCheckpoint, TensorBoard, ReduceLROnPlateau
from keras.initializers import glorot_uniform
from tensorflow.keras.models import load_model
np.random.seed(24)

# 1. Functions

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

def prep(data):
        
    # FFW Input
    
    x_ffw = pd.concat([data['toy'],
                       #data['trend'],
                       dummies(data,'dow')],
                      axis = 1)
    x_ffw = x_ffw.values
    
    # CNN inputs
    x_inst = data[[col for col in rolling_train if col.startswith('matInt')]].values
    x_tem = data[[col for col in rolling_train if col.startswith('matTem.')]].values
    x_tem95 = data[[col for col in rolling_train if col.startswith('matTem95')]].values
    x_lag = data[[col for col in rolling_train if col.startswith('matLag')]].values
    
    ## Scaling data before stacking
    x_ffw_scaler,x_inst_scaler,x_tem_scaler,x_tem95_scaler,x_lag_scaler=MinMaxScaler(),MinMaxScaler(),MinMaxScaler(),MinMaxScaler(),MinMaxScaler()
    
    # Fitting scalers only on first year of data
    
    filtered_train = rolling_train[rolling_train['ym'].isin(ym_unique[:12])] 
    index_train = len(filtered_train)
    
    x_ffw_scaler.fit(x_ffw[:index_train])
    x_inst_scaler.fit(x_inst[:index_train])
    x_tem_scaler.fit(x_tem[:index_train])
    x_tem95_scaler.fit(x_tem95[:index_train])
    x_lag_scaler.fit(x_lag[:index_train])
    
    # Scale the whole data
    
    x_ffw   = x_ffw_scaler.transform(x_ffw)
    x_inst  = x_inst_scaler.transform(x_inst)
    x_tem   = x_tem_scaler.transform(x_tem)
    x_tem95 = x_tem95_scaler.transform(x_tem95)
    x_lag   = x_lag_scaler.transform(x_lag)
    
    # Stack the cnn inputs
    x_cnn_1 = np.stack((x_tem,x_inst), axis = 2)
    x_cnn_2 = np.stack((x_tem95,x_inst), axis = 2)
    x_cnn_3 = np.stack((x_lag,x_inst), axis = 2)

    # Output
    y_scaler = MinMaxScaler()
    y = data['load'].values
    
    y.shape = (y.shape[0],1)
    
    y_scaler.fit(y[:index_train])
    y = y_scaler.transform(y)
    
    
    tables = x_ffw,x_cnn_1,x_cnn_2,x_cnn_3, y
    scalers = x_ffw_scaler,x_inst_scaler,x_tem_scaler,x_tem95_scaler,x_lag_scaler,y_scaler

    return(tables,scalers)

def conv1D_simple(input_cnn, num_filters):
    cnn = Conv1D(kernel_size = 3, filters = num_filters,
             kernel_initializer=glorot_uniform(seed=seed),
             padding='same')(input_cnn)
    cnn = BatchNormalization()(cnn)
    cnn = Activation('relu')(cnn)
    cnn = MaxPooling1D(pool_size=2)(cnn)
    #cnn = Dropout(0.1)(cnn) # back to none if it does not work
    
    return cnn

def conv1D_block(input_cnn, num_filters):
    cnn = Conv1D(kernel_size = 3, filters = num_filters,
             kernel_initializer=glorot_uniform(seed=seed),
             padding='same')(input_cnn)
    cnn = BatchNormalization()(cnn)
    cnn = Activation('relu')(cnn)
    
    cnn = Conv1D(kernel_size = 3, filters = num_filters,
             kernel_initializer=glorot_uniform(seed=seed),
             padding='same')(cnn)
    cnn = BatchNormalization()(cnn)
    cnn = Activation('relu')(cnn)
    cnn = MaxPooling1D(pool_size=2)(cnn)
    #cnn = Dropout(0.2)(cnn)
    
    return cnn

# 2. Data Preparation

rolling_train = pd.read_csv('1. Data/low-resolution_train.csv')
rolling_train['ymd'] = pd.to_datetime(rolling_train['ymd'])
rolling_train['dateDay'] = [dt.date(d) for d in rolling_train['ymd']]

rolling_test = pd.read_csv('1. Data/low-resolution_test.csv')
rolling_test['ymd'] = pd.to_datetime(rolling_test['ymd'])
rolling_test['dateDay'] = [dt.date(d) for d in rolling_test['ymd']]

features = [col for col in rolling_train if col.startswith('matLag') or col.startswith('matTem') or col.startswith('matInt')]
selection = features + ['dateDay', 'dow', 'toy','load','ymd','ym','trend']

rolling_train = rolling_train[selection].dropna()
rolling_test = rolling_test[selection].dropna()

ym_unique = rolling_train['ym'].unique()

training, train_scalers = prep(rolling_train)

print(training[0].shape)
print(training[1].shape)
print(training[2].shape)
print(training[3].shape)
print(training[4].shape)

# 3. Modelling

seed = 24

input_ffw = Input(shape=(training[0].shape[1],))
input_cnn_1 = Input(shape=(training[1].shape[1],training[1].shape[2]))
input_cnn_2 = Input(shape=(training[2].shape[1],training[2].shape[2]))
input_cnn_3 = Input(shape=(training[3].shape[1],training[3].shape[2]))

ffw = Dense(30,use_bias=True,activation=None,kernel_initializer=glorot_uniform(seed=seed))(input_ffw)
ffw = BatchNormalization()(ffw)
ffw = Dropout(0.1)(ffw)
ffw = Activation('relu')(ffw)

cnn_1 = conv1D_simple(input_cnn_1, 8)
cnn_1 = conv1D_simple(cnn_1, 16)
cnn_1 = conv1D_simple(cnn_1, 32)
cnn_1 = Flatten()(cnn_1)
cnn_1 = Dropout(0.1)(cnn_1)

cnn_2 = conv1D_simple(input_cnn_2, 8)
cnn_2 = conv1D_simple(cnn_2, 16)
cnn_2 = conv1D_simple(cnn_2, 32)
cnn_2 = Flatten()(cnn_2)
cnn_2 = Dropout(0.1)(cnn_2)

cnn_3 = conv1D_simple(input_cnn_3, 8)
cnn_3 = conv1D_simple(cnn_3, 16)
cnn_3 = conv1D_simple(cnn_3, 32)
cnn_3 = Flatten()(cnn_3)
cnn_3 = Dropout(0.1)(cnn_3)

x = concatenate([ffw, cnn_1, cnn_2, cnn_3])
#x = Dense(300,use_bias=True,kernel_initializer=glorot_uniform(seed=seed))(x)
#x = BatchNormalization()(x)
#x = Activation('relu')(x)
#x = Dropout(0.4)(x)
x = Dense(1,use_bias=True,activation = None,kernel_initializer=glorot_uniform(seed=seed))(x)
x = Activation('relu')(x)

model = Model(inputs=[input_ffw, input_cnn_1,input_cnn_2,input_cnn_3], outputs=x)

optimizer = Nadam(lr=1e-4)

model.compile(loss='mse', optimizer=optimizer, metrics=['mae'])

model.summary()

# Save model
model.save('weights/multi-resolution_CNN_rolling')

# Rolling origin training
fitted_signal = []
pred_signal = []
for i in range(len(ym_unique)-12):
    print("================ Iteration",i+1, "================")
    model = load_model('Final_Hybrid_rolling')
    filtered_train = rolling_train[rolling_train['ym'].isin(ym_unique[:i+12])] #starting to train with the full first year
    index_train = len(filtered_train)
    
    filtered_test = rolling_train[rolling_train['ym'] == ym_unique[i+12]] #starting to train with the full first year
    index_test = len(filtered_test)
    
    # train
    history = model.fit(x=[training[0][:index_train],training[1][:index_train],training[2][:index_train],training[3][:index_train]],
                        y=training[4][:index_train],
                        batch_size=16, # back to 16 if it does not work
                        epochs=300,# back to 300 if it does not work
                        verbose = 0)

    # fitted signal
    fitted_signal.append(list(train_scalers[5].inverse_transform(model.predict([training[0][:index_train],
                                                          training[1][:index_train],
                                                          training[2][:index_train],
                                                          training[3][:index_train]]))))

    # pred signal
    pred_signal.append(train_scalers[5].inverse_transform(model.predict([training[0][index_train:index_train+index_test],
                                                          training[1][index_train:index_train+index_test],
                                                          training[2][index_train:index_train+index_test],
                                                          training[3][index_train:index_train+index_test]])))
pred_signal = np.vstack(pred_signal)

for i in range(48):
    path = '2. Pred Signals/2. NNs/2. Hybrid/fitted/multi-resolution_CNN_rolling_' + str(i) + '.csv'    
    np.savetxt(path, np.vstack(fitted_signal[i]), delimiter=',')


np.savetxt('2. Pred Signals/2. NNs/2. Hybrid/multi-resolution_CNN_rolling_pred_1.csv', pred_signal, delimiter=',')