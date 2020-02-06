#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# matlab engine is used to execute pqmodel function in pqmodel.m
import matlab.engine
import pandas as pd
import numpy as np

# python wrapper for function found in pqmodel.m
# see that function for detailed documentation about its inner workings
# returns
def pqmodel(num_samples=10, samp_freq=16000, fund_freq=50, num_cycles=10, amplitude=1):
    '''
    This is a python wrapper for function found in pqmodel.m
    (See https://data.mendeley.com/datasets/6kmkk9bjdx/1)
    Returns the matlab matrix as a numpy array.
    See that function for detailed documentation about its inner workings, but
    I'll copy the details about input and ouput here:

    %Input Parameters
    %   ns: Number of total signals/samples to generate of each class (default
    %   value: 10). Range of possible values: 1-1000000
    %   fs: Sampling frequency of the signals(default value: 16kHz). Range
    %   of possible values:200-30000Hz
    %   f: Fundamental frequency of the electrical signal (default value:
    %   50Hz). Range of possible values:40-100Hz
    %   n: Number of total cycles (periods) of the fundamental frequency
    %   contained in each sample (default value: 10). Range of possible
    %   values:3-100
    %   A: Amplitude (default value: per unit). Range of possible
    %   values:0.1-400000
    %
    %Output
    %   out: Matrix of dimensions(Number of samples per class, Number of points per signal, Number of classes)
    '''
    eng = matlab.engine.start_matlab()
    # Converting to floats because matlab doesn't like the ints.
    samples = eng.pqmodel(float(num_samples), float(samp_freq), float(fund_freq), float(num_cycles), float(amplitude))
    return np.asarray(samples)

def pqmodel_df(num_samples=10, samp_freq=16000, fund_freq=50, num_cycles=10, amplitude=1):
    '''
    This function returns the output of the above function to a pandas DataFrame.
    The DataFrame gives the id, class number, and samples of each wave.
    The data can be easily written to a csv from here using: .to_csv('filename.csv', index=False)
    '''
    # Call matlab function
    m = pqmodel(num_samples, samp_freq, fund_freq, num_cycles, amplitude)

    num_classes = np.shape(m)[2]

    df = pd.DataFrame()
    for class_num in range(1, num_classes + 1):
        classDF = pd.DataFrame(m[:,:,class_num-1], columns=[f's{i}' for i in range(np.shape(m)[1])])
        classDF['class'] = class_num
        df = df.append(classDF)

    df.reset_index(drop=True, inplace=True)
    df['id'] = [i for i in range(len(df))]
    df = df[['id'] + ['class'] + df.columns[:-2].tolist()]

    return df
