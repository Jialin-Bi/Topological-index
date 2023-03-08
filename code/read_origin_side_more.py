# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 19:53:06 2022

@author: bijialin
"""
#Read the coordinates of the corresponding elements in the mutation pdb and store them

import pandas as pd
import numpy as np
import csv

def readpdb(name):
    pdbfile = open(name)
    lines = pdbfile.read().splitlines()
    PRO = []
    for i in range(0,len(lines)):
        line = lines[i]
        if line[0:4] == 'ATOM':
            AtomType = line[13:14]
            chain = line[21:22]
            resid = line[22:26].strip()
            if AtomType in pro_ele_list and chain in ['A','E']:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                PRO.append([chain,resid,AtomType,x,y,z])
                continue
            else:
                continue
        else:
            continue 
    return PRO	


if __name__ == "__main__":
    resid1='C'#'O'
    area='origin_mutation'#binding
    pro_ele_list = [resid1]
    
    Note2=open('name.txt',mode='r')
    file_name= eval(Note2.read())
    Note2.close()
    
    if area =='origin_mutation':
        pdbname='6M0J_m.pdb'

    ename1='0'
    for ename in file_name[:]: 
        if ename[-5:-2]!=ename1:
            ename1=ename[-5:-2]
            print(ename)
            mut_mutation='./6m0j-rbd1-pdb-file/'+ename+'/pdbfile/'+pdbname
            pro_mut_m = readpdb(mut_mutation)
    

            ##storage excel document
            data=pd.DataFrame(pro_mut_m)
            data.columns=['chain','num','atom','X','Y','Z']
            writer=pd.ExcelWriter("./"+area+"_coodinates_"+resid1+'_'+resid1+"/6M0J_"+ename[-5:-2]+'_'+resid1+'_'+resid1+"_"+area+"_coodinates.xlsx")
            data.to_excel(writer) 
            writer.save()
        