import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd

testfile='../data/data_set_ALL_AML_independent.csv'
trainfile='../data/data_set_ALL_AML_train.csv'
patient_cancer='../data/actual.csv'

train = pd.read_csv(trainfile).transpose()
test = pd.read_csv(testfile).transpose()
patient_cancer = pd.read_csv(patient_cancer).transpose()

print train.transpose()