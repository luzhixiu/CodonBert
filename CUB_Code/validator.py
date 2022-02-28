#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 14:22:44 2019

@author: lu
"""

# This code is used to analyze the connection, difference corelation of two datasets
import scipy.stats as stats
import scipy.stats as ss
import matplotlib.pyplot as plt
import numpy as np
import random
import math
from numbers import Number

#testing functions
#random.seed(0)
#list1=random.sample(range(100), 100)
#list2=random.sample(range(100), 100)




def testCorelation(x,y,corelationFunction):
    import scipy.stats as ss
    if "pearson" in corelationFunction:
#        print ("p value %f"%stats.pearsonr(x, y)[1])
        return  stats.pearsonr(x, y)[0]
    elif "spearman" in corelationFunction:
        return stats.spearmanr(x,y,nan_policy="omit")[0]
    elif "kendall" in corelationFunction:
        return stats.kendalltau(x,y,nan_policy="omit")[0]

def plotCorelation(x,y,xLabel="list1",yLabel="list2",logScale="no",showCorelation="yes"):
    plt.figure(dpi=300)
    if(logScale=="yes"):
        x=logify(x)
        y=logify(y)
    if "yes" in showCorelation:
        plt.text(0.1,0.68,"R =: %0.4f"%(testCorelation(x,y,"pearson")))
    x=np.array(x)
    y=np.array(y)
    
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    b=estimate_coef(x,y)
    y_pred=b[0]+b[1]*x
    plt.plot(x, y_pred, color = "g",linewidth=1)
    plt.scatter(x,y,s=0.5)
    plt.show()

   

def testTopAndBotAgreeMent(x,y):
    listSize=len(x)
    print ("Total genes samples %d"%listSize)
    selectPercentList=[0.1,0.5,1,5,10,15,20,25]
    disagreeIndexList=[]
    for selectPercentage in selectPercentList:
        disagreeIndexList=[]
        botSelectCutOff=int(listSize*0.01*selectPercentage)
        topSelectCutOff=int(listSize*0.01*(100-selectPercentage))
        list1=x
        list2=y
        rankList1=ss.rankdata(list1,method="min")
        rankList2=ss.rankdata(list2,method="min")
        # options are: min, max,averaage, dense,ordinal
    #    print rankList1
    #    print rankList2    
        selectRange=int(listSize*0.01*selectPercentage)
        if selectRange==0:
            print ("provided list is too small")
            selectRange=0.000000000001
        print ("Selected Samples: %d"%(selectRange))
        
        cnt=0.0
    
    #    print "Top cutoff: %d"%topSelectCutOff
        for i in range(listSize):
            if rankList1[i]>topSelectCutOff and rankList2[i]>topSelectCutOff:
                cnt+=1
        percentage=cnt/selectRange*100
        print ("Top %s%%: Intersect: %d (%0.2f%%)"%(selectPercentage,cnt,percentage))
    
    
        cnt=0.0
    #    print "Bot cutoff: %d"%botSelectCutOff
        for i in range(listSize):
            if rankList1[i]<botSelectCutOff and rankList2[i]<botSelectCutOff:
                cnt+=1
            else:
                disagreeIndexList.append(i)
        percentage=cnt/selectRange*100
        print ("Bot %s%%: Intersect: %d (%0.2f%%)"%(selectPercentage,cnt,percentage))
        print()
    print ("Used min as the ranking Method (allows ties), ranking list [0,2,3,2] will give [1,2,4,2]")
    
    
def logify(myList):
#    print (myList)
    logList=[]
    for i in myList:
        if i<=0:
            i=0.0001# treat 0 as a really small number so log of it makes sense
        if not isinstance(i, Number):
            i=0.0001
        logList.append(math.log(float(i)))
        print (logList)
#        logList.append(math.log(i+1))
    return logList

def validate(x,y,corelationFunction="pearson",logScale="no",xLabel="x",yLabel="y",showCorelation="yes"):
    if(logScale=="yes"):
        x=logify(x)
        y=logify(y)
#    testTopAndBotAgreeMent(x,y)
    print ("corelation is: %s "%testCorelation(x,y,corelationFunction))
    plotCorelation(x,y,xLabel=xLabel,yLabel=yLabel)
    



#====================private helper functions
def estimate_coef(x, y): 
    # number of observations/points 
    n = np.size(x) 
    # mean of x and y vector 
    m_x, m_y = np.mean(x), np.mean(y) 
    # calculating cross-deviation and deviation about x 
    SS_xy = np.sum(y*x) - n*m_y*m_x 
    SS_xx = np.sum(x*x) - n*m_x*m_x 
    # calculating regression coefficients 
    b_1 = SS_xy / SS_xx 
    b_0 = m_y - b_1*m_x 
    return(b_0, b_1)  
    
