from numpy import *
from pandas import *
import csv
import os
import numpy
import pandas as pd
import math
import random

from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import cross_val_score
from scipy.optimize import fsolve 

Crispr_data2= '/home/nonojh/crisprCas9/SSL/equal/data/Crispr_Data2.txt'
Crispr_positive='/home/nonojh/crisprCas9/SSL/equal/data/Crispr_SSL_positive.txt'
Crispr_negative='/home/nonojh/crisprCas9/SSL/equal/data/Crispr_SSL_negative'


Name_data_list = ['cfdScore','MITwebsite','MITscores','CropitScores','CCTopSores','PhastCons','phylop','Promoter',
                  'Enhancer','segwayPromoter','segwayEnhancer']
Distance_path='/home/nonojh/crisprCas9/SSL/equal/distance/'

Result = '/home/nonojh/crisprCas9/SSL/equal/result/'
input_list=[2,3,4,5,6,7,8,9,10,11,12]

for index0 in range(0,166):
    index=0
    dist_train = [[],[],[],[],[],[],[],[],[],[],[]]
    for i in range(0,11):
        #read input data
        filename=Distance_path+str(index)+'/Distance_feature'+str(i)+'.txt'
        csvfile = file(filename,'rb')
        reader_crispr = csv.reader(csvfile)
        for row in reader_crispr:
            dist_train[i].append(row)
        csvfile.close()
    num=len(row)

    L_list = [[],[],[],[],[],[],[],[],[],[],[]]
    for i in range(0,11):
        Index = mat(numpy.eye(num))
        D = mat(numpy.eye(num))
        for i1 in range(0,num):
            d = 0.0
            for j1 in range(0,num):
                d = d+float(dist_train[i][i1][j1])
            D[i1,i1] = d
            L = D-mat(dist_train[i]).astype('float64')
            L_list[i] = L


    ##Read  Outcome Data
    o_true = []
    filename=Distance_path+str(index)+'/Crispr_Data'+str(index)+'.txt'
    csvfile = file(filename,'rb')
    reader = csv.reader(csvfile)
    for row in reader:
        o_true.append([float(row[1])])
    csvfile.close()

    l_5fold=[]
    l = range(0,num)
    l_train,l_test = train_test_split(l,test_size=0.2)
    l_5fold.append(l_test)
    l0_train,l0_test = train_test_split(l_train,test_size=0.5)
    l1_train,l1_test = train_test_split(l0_train,test_size=0.5)
    l2_train,l2_test = train_test_split(l0_test,test_size=0.5)
    l_5fold.append(l1_train)
    l_5fold.append(l1_test)
    l_5fold.append(l2_train)
    l_5fold.append(l2_test)


    directory_name=Result+str(index)
    if not os.path.exists(directory_name):
        os.mkdir(directory_name)

    for Data_flag in range(0,11):
        filename=directory_name+'/result_feature'+str(Data_flag)+'.txt'
        wcsvfile = file(filename,'wb')
        writer = csv.writer(wcsvfile)
        max_auc = 0
        best_u = 0
            
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
            for j in range(0,num):
                if j in l_train:
                    o_train.append(o_true[j])
                elif j in l_test:
                    o_train.append([0.0])
    
            auc=-1
            u=0
            count=0
            for u0 in numpy.arange(0,1.001,0.001):
                f1= Index+u0*L_list[Data_flag]
                f = mat(f1.I)*mat(o_train)
                tmp=roc_auc_score(o_true,f)
                if (auc<tmp):
                    auc=tmp
                    u=u0
    
            best_u +=u
            max_auc +=auc
            writer.writerow([round(u,6),round(auc,6)])
        
        writer.writerow([round(best_u/5,6),round(max_auc/5,6)])
        wcsvfile.close()
    break
