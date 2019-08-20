from numpy import *
from pandas import *
import csv
import numpy
import pandas as pd
import math
import random

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

features_list=['cfdScore','MITwebsite','MITscores','CropitScores','CCTopSores','scores','phylop','Promoter','Enhancer','segwayPromoter','segwayEnhancer']

input_list=[2,3,4,5,6,7,8,9,10,11,12]
feature_range=[]
raw_num=0
for i in range(0,11):
    min_value=10000
    max_value=-10000
    csvfile = file(Data_crispr,'rb')
    reader = csv.reader(csvfile)
    for row in reader:
        if i==0:
            raw_num=raw_num+1
        if row[input_list[i]]=='NA':
            current = 0
        else:
            current=float(row[input_list[i]])
        if current > max_value:
            max_value = current
        if current < min_value:
            min_value = current
    csvfile.close()
    
    feature_range.append(array([min_value,max_value,max_value-min_value]))


def WeightGaussian(a,b,sig):
    distance = round(math.exp(-numpy.sum(numpy.square(array(a) - array(b)))/numpy.square(sig)),8)
    return distance


def Knn_GetDistance(data,num,k,flag):
    wcsvfile = file(Distance_path3,'wb')
    writer = csv.writer(wcsvfile)
    
    for i in range(0,num):
        A = data[i]
        arr=[0]*num
        for j in range(0,num):
            if i==j:
                arr[j]==0
            else:
                B = data[j]
                temp = WeightGaussian(A,B,1)
                arr[j] = temp

        k_arr = [0]*num

        A = list(arr)
        A_sort = argsort(A)
        B = A_sort[:-(k+1):-1]
        for j in B:
            k_arr[j] = A[j]

        writer.writerow(k_arr)

    wcsvfile.close()


k=100
i_train = []
csvfile = file(Data_crispr,'rb')
reader = csv.reader(csvfile)
for row in reader:
    tmp = []
    for j in range(2,13):
        if row[j]=="NA":
            tmp.append(0)
        else:
            tmp.append(float((float(row[j])-feature_range[j-2][0])/feature_range[j-2][2])) #normalization
    i_train.append(tmp)
csvfile.close()
wcsvfile = file('/home/nonojh/crisprCas9/SSL/data/Crispr_SSL_pre2.txt','wb')
writer = csv.writer(wcsvfile)
writer.writerows(i_train)
wcsvfile.close()
Knn_GetDistance(i_train,raw_num,k,0)
