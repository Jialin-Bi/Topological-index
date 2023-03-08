# -*- coding: utf-8 -*-
"""
Created on Mon May 16 20:59:45 2022

@author: bijialin
"""

import pandas as pd
import os
from scipy.stats import pearsonr
import math
import numpy as np
import operator
import random
from scipy import stats
import csv
import matplotlib.pyplot as plt


def pearson(x,y):
   x = pd.Series(x)
   y = pd.Series(y)
   r = x.corr(y, method="pearson")
   return r

def pearson_stat(x,y):
    x = pd.Series(x)
    y = pd.Series(y)
    r,p = stats.pearsonr(x, y)
    return r,p


def spearman(x,y):
   x = pd.Series(x)
   y = pd.Series(y)
   r = x.corr(y, method="spearman")
   return r

def spearman_stat(x,y):
    x = pd.Series(x)
    y = pd.Series(y)
    r,p = stats.spearmanr(x, y)
    return r,p

if __name__ == '__main__':
    
    data=pd.read_excel("./6M0J-RBD1-index.xlsx")
    data_key=data.keys()
    data_value=data.values
    
    Note2=open('name.txt',mode='r')
    file_name= eval(Note2.read())
    Note2.close()
    
    #delta delta G
    data_dict={}   
    for i in range(len(data_value[:,1])):
        data_dict[str(data_value[i,4])+'_'+str(data_value[i,5])+'_'+str(data_value[i,6])]=data_value[i,8]
        
    name_keys=[]
    for na in file_name:
        name_keys.append(na[-7:])
        
    resid1='C'#'C'
    resid2='C'#'C
    area='_mutation'    
    # ad='_low'#'_up' or '_low' or '_c'
    ad='_low'
    dim='_A1_'#'_A0_'

    
    #mutation index
    mutation_dict={}
    mutation_dict2={}
    site_dict={}
    flag=1
    # ename='6M0J_E_A_411_D'

    for ename in file_name:
        result_name=ename+'_'+resid1+'_'+resid1+area+dim+'index.csv'
        file = pd.read_csv('./result/mutation_result'+dim+resid1+'_'+resid1+ad+'/'+result_name,index_col=0)
        file_col=file.columns
        file_index=file.index#centrality
        file_value=file.values#flitration
        key=ename[-7:]
        col=list(file_col)
        index=list(file_index)
        if flag ==1:
            flag=0
            for m in range(len(col)):
                for f in range(len(index)):
                    fla=index[f]
                    mutation_dict[col[m]+'_'+str(fla)]={}
                    mutation_dict[col[m]+'_'+str(fla)][key]=file_value[f,m]
                mutation_dict[col[m]+'_sum']={}
                mutation_dict[col[m]+'_sum'][key]=sum(file_value[:,m])    
        else:
            for m in range(len(col)):
                for f in range(len(index)):
                    fla=index[f]
                    mutation_dict[col[m]+'_'+str(fla)][key]=file_value[f,m]
                mutation_dict[col[m]+'_sum'][key]=sum(file_value[:,m])  


    flag=1
    for ename in range(333,527):
        result_name='6M0J_origin_'+str(ename)+'_'+resid1+'_'+resid2+area+dim+'index.csv'
        file=pd.read_csv('./result/origin_mutation_result'+dim+resid1+'_'+resid2+ad+'/'+result_name,index_col=0)
        file_col=file.columns
        file_index=file.index
        file_value=file.values
        
        key=str(ename)
        col=list(file_col)
        index=list(file_index)
        if flag ==1:
            flag=0
            for m in range(len(col)):
                for f in range(len(index)):
                    fla=index[f]
                    mutation_dict2[col[m]+'_'+str(fla)]={}
                    mutation_dict2[col[m]+'_'+str(fla)][key]=file_value[f,m]
                mutation_dict2[col[m]+'_sum']={}
                mutation_dict2[col[m]+'_sum'][key]=sum(file_value[:,m])    
        else:
            for m in range(len(col)):
                for f in range(len(index)):
                    fla=index[f]
                    mutation_dict2[col[m]+'_'+str(fla)][key]=file_value[f,m]
                mutation_dict2[col[m]+'_sum'][key]=sum(file_value[:,m])  
        

    
    '''
    print("***********************************************************************")
    
    yuansu='D' #'P'#'Charged'#'V'
    result_name='mutation_result'+dim+resid1+'_'+resid2+ad+'fenzu_add_'+yuansu+'.csv'
    # result_name='mutation_result'+dim+resid1+'_'+resid2+ad+'fenzu_add2.csv'
    File1=open('./'+result_name,"w+",newline='') 
    writing_File1=csv.writer(File1)
    
    
    Charged={'R','H','K','D','E'}
    Polar={'S','T','N','Q'}
    Hydrophobic={'A','V','I','L','M','F','Y','W'}
    Special={'C','G','P'}
    # one={'D','G','N','P','R','V','E','H','K','T'}
    one={yuansu}
    fenzu=one
    # fenzu=Alanine
    
    method=list(mutation_dict.keys())
    count=0
    method_result=dict()#pearson
    method2_result=dict()#spearman
    method1_result_all=dict()#pearson
    method2_result_all=dict()#spearman
    for key in method:#
        print(key)
        x=[]
        y=[]       
        for key2 in name_keys:
            name=key2.split('_')
            if name[0] in fenzu:
                x.append(data_dict[key2])
                val=mutation_dict[key][key2]+mutation_dict2[key][key2[2:5]]
                y.append(np.nan_to_num(val))
        p,P=pearson_stat(x,y)
        s,P2=spearman_stat(x,y)
        print('pearson= %f，p = %f'%(p,P))
        print('spearman = %f，s = %f'%(s,P2))
        method_result[key]=p
        method2_result[key]=s
        method1_result_all[key]=p
        method2_result_all[key]=s
        count+=1
        if count ==len(index)+1:
            writing_File1.writerow(['method']+list(method_result.keys()))  
            writing_File1.writerow(['P']+list(method_result.values()))   
            writing_File1.writerow(['S']+list(method2_result.values()))  
            method_result=dict()
            method2_result=dict()
            count=0            
          
    method1=list(method1_result_all.values())
    method2=list(method2_result_all.values())
    writing_File1.writerow([max(np.nan_to_num(method1))]+[min(np.nan_to_num(method1))])  
    writing_File1.writerow([max(np.nan_to_num(method2))]+[min(np.nan_to_num(method2))])  
    # writing_File1.writerow([max(method1_result_all.values())]+[min(method1_result_all.values())])     
    # writing_File1.writerow([max(method2_result_all.values())]+[min(method2_result_all.values())])        
    File1.close()
    
    '''
    
    
    print("***********************************************************************")
    yuansu_all=['R','H','K','D','E','S','T','N','Q','V','I','L','F','Y','W','C','G','P']
    for yuansu in yuansu_all[:]:
        #correlation
        result_name='mutation_result'+dim+resid1+'_'+resid2+ad+'fenzu_add_'+yuansu+'.csv'
        File1=open('./'+result_name,"w+",newline='') 
        writing_File1=csv.writer(File1)
        
        
        Charged={'R','H','K','D','E'}
        Polar={'S','T','N','Q'}
        Hydrophobic={'A','V','I','L','M','F','Y','W'}
        Special={'C','G','P'}
        one={yuansu}
        fenzu=one
        # fenzu=Alanine
        
        #correlation
        method=list(mutation_dict.keys())
        count=0
        method_result=dict()#pearson
        method2_result=dict()#spearman
        method1_result_all=dict()#pearson
        method2_result_all=dict()#spearman
        for key in method:
            print(key)
            x=[]
            y=[]       
            for key2 in name_keys:
                name=key2.split('_')
                if name[0] in fenzu:
                    x.append(data_dict[key2])
                    val=mutation_dict[key][key2]+mutation_dict2[key][key2[2:5]]
                    y.append(np.nan_to_num(val))
            p,P=pearson_stat(x,y)
            s,P2=spearman_stat(x,y)
            print('pearson= %f，p = %f'%(p,P))
            print('spearman = %f，s = %f'%(s,P2))
            method_result[key]=p
            method2_result[key]=s
            method1_result_all[key]=p
            method2_result_all[key]=s
            count+=1
            if count ==len(index)+1:
                writing_File1.writerow(['method']+list(method_result.keys()))  
                writing_File1.writerow(['P']+list(method_result.values()))   
                writing_File1.writerow(['S']+list(method2_result.values()))  
                method_result=dict()
                method2_result=dict()
                count=0            
              
        method1=list(method1_result_all.values())
        method2=list(method2_result_all.values())
        writing_File1.writerow([max(np.nan_to_num(method1))]+[min(np.nan_to_num(method1))])  
        writing_File1.writerow([max(np.nan_to_num(method2))]+[min(np.nan_to_num(method2))])  
        File1.close()
