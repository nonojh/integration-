from numpy import *
from pandas import *
import time
import csv
import numpy
import pandas as pd
import math
import random

from sklearn.metrics import roc_auc_score

from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import cross_val_score
from scipy.optimize import fsolve 


import csv
import numpy
import pandas as pd
import math
import random

raw_num = 25332

Distance_cfdScore = '/home/nonojh/crisprCas9/SSL/data/cfdScore_distance.txt'
Distance_MITwebsite = '/home/nonojh/crisprCas9/SSL/data/MITwebsite_distance.txt'
Distance_MITscores = '/home/nonojh/crisprCas9/SSL/data/MITscores_distance.txt'
Distance_CropitScores = '/home/nonojh/crisprCas9/SSL/data/CropitScores_distance.txt'
Distance_CCTopSores = '/home/nonojh/crisprCas9/SSL/data/CCTopSores_distance.txt'
Distance_scores = '/home/nonojh/crisprCas9/SSL/data/scores_distance.txt'
Distance_phylop = '/home/nonojh/crisprCas9/SSL/data/phylop_distance.txt'
Distance_Promoter = '/home/nonojh/crisprCas9/SSL/data/Promoter_distance.txt'
Distance_Enhancer = '/home/nonojh/crisprCas9/SSL/data/Enhancer_distance.txt'
Distance_segwayPromoter = '/home/nonojh/crisprCas9/SSL/data/segwayPromoter_distance.txt'
Distance_segwayEnhancer = '/home/nonojh/crisprCas9/SSL/data/segwayEnhancer_distance.txt'

Data_crispr = '/home/nonojh/crisprCas9/SSL/data/Crispr_SSL_pre.txt'


sigma_cfd = 0.001
sigma_Mitweb = 0.02
sigma_Mit = 0.034
sigma_Cropit = 62.5


Distance_path2='/home/nonojh/crisprCas9/SSL/data/distance_method2.txt'
Distance_path3='/home/nonojh/crisprCas9/SSL/data/distance_method3.txt'

Distance_path=[Distance_cfdScore,Distance_MITwebsite,Distance_MITscores,Distance_CropitScores,
            Distance_CCTopSores,Distance_scores,Distance_phylop,Distance_Promoter,Distance_Enhancer,
               Distance_segwayPromoter,Distance_segwayEnhancer]


sigma_list=[0.001,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]

Features_list=['cfdScore','MITwebsite','MITscores','CropitScores','CCTopSores','scores','phylop','Promoter','Enhancer','segwayPromoter','segwayEnhancer']


dist_train = []

#read input data
csvfile = file(Distance_path3,'rb')
reader = csv.reader(csvfile)
for row in reader:
    dist_train.append(row)
csvfile.close()

L = []
Index = mat(numpy.eye(raw_num))
D = mat(numpy.eye(raw_num))
for i in range(0,raw_num):
    d = 0.0
    for j in range(0,raw_num):
        d = d+float(dist_train[i][j])
    D[i,i] = d
L = D-mat(dist_train).astype('float64')

max_auc = 0
#X_train, X_test, y_train, y_test = train_test_split(L_list[Data_flag], o_true, test_size=0.3, random_state=42)

l_5fold=[]
l = range(0,raw_num)
l_train,l_test = train_test_split(l,test_size=0.2)
l_5fold.append(l_test)
l0_train,l0_test = train_test_split(l_train,test_size=0.5)
l1_train,l1_test = train_test_split(l0_train,test_size=0.5)
l2_train,l2_test = train_test_split(l0_test,test_size=0.5)
l_5fold.append(l1_train)
l_5fold.append(l1_test)
l_5fold.append(l2_train)
l_5fold.append(l2_test)

def func(x):
    u0 = float(x)
    f1= Index+u0*L
    y=mat(o_train)
    f2 = f1.I*y
    return (f-y).T*(f-y)+u0*f2.T*L*f2

for i in range(0,5):
    l_test = l_5fold[i]
    l_train = []
    for j in range(0,5):
        if i==j:
            continue
        else:
            for k in range(0,len(l_5fold[j])):
                l_train.append(l_5fold[j][k] )

    o_train = []
    for j in range(0,raw_num):
        if j in l_train:
            o_train.append(o_true[j])
            elif j in l_test:
                o_train.append([0.0])

    u = fsolve(func,-1)
    f = Index+u*L.I*mat(o_train)
    max_auc += roc_auc_score(o_true,f)
    break

wcsvfile = file('/home/nonojh/crisprCas9/SSL/data/result.txt','wb')
writer = csv.writer(wcsvfile)
writer.writerow(round(max_auc,3))
wcsvfile.close()
