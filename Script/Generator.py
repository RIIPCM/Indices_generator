# -*- coding: utf-8 -*-
"""
Created on Fri Oct 01 10:57:25 2022

@authors: Bartosz Kowal (1), Krzysztof Smalara(1), Jiri Mazurek(2)
(1) Rzeszów University of Technology, (2) Silesian University in Opava

#License
Smalara, K., Kowal, B., Mazurek, J. (2022). Inconsistency matrix generator with Saaty scale with several indices. https://github.com/RIIPCM/Indices_generator
"""

import numpy as np
import scipy.linalg as la
import random
import csv
from datetime import datetime
import time

def create_matrix(size):    
    matrix = [[0 for x in range(size)] for y in range(size)]
    
    for n in range(size):
        for m in range(size):
            if n==m:
                matrix[n][m]=1
    
    return matrix


def return_max_eigenvalue(n_matrix):

    nump_table = np.array(n_matrix).astype(np.float64)
    evals, evecs = la.eig(nump_table)
    evals = evals.real
    eigen_value_max = max(evals)
    

    return eigen_value_max


def get_RI(matrix):
    size=len(matrix)
    if size == 3:
        RI=0.5247
    elif size ==4:
        RI=0.8815
    elif size ==5:
        RI=1.1086     
    elif size ==6:
        RI=1.2479
    elif size ==7:
        RI=1.3417     
    elif size ==8:
        RI=1.4056
    elif size ==9:
        RI=1.4499    
    elif size ==10:
        RI=1.4854
    elif size ==11:
        RI=1.5140   
    elif size ==12:
        RI=1.5365           
    elif size ==13:
        RI=1.5551    
    elif size ==14:
        RI=1.5713
    else:
        return 0
    return RI         
        
        

def return_CR(matrix):
    size=len(matrix)
    CI = (return_max_eigenvalue(matrix) - size)/(size-1)
    RI = get_RI(matrix)
    CR = CI / RI
    return CR


def fill_matrix(matrix):
    
    choice_list = [1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3, 1/2, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    
    size=len(matrix)
    
    for n in range(size):
        for m in range(size):
            if n > m:
                matrix[n][m]=random.choice(choice_list)
    for n in range(0,size):
        for m in range(0,size):
            if n < m:
                matrix[n][m]=1/matrix[m][n]         
                
    return matrix



###Indicies
def dynamic_vectors(matrix):
  nump_table=np.array(matrix).astype(np.float64)
  evals, evecs = la.eig(nump_table)
  abs_real_evecs = np.absolute(evecs.real)
  vector_sum=0
  list_of_vecs=[]
  for i in range(len(matrix)):
    vector_sum=vector_sum+abs_real_evecs[i][0]
  for i in range(len(matrix)):
    list_of_vecs.append(abs_real_evecs[i][0]/vector_sum)
  return list_of_vecs


def Pelaez_Lamata(matrix):
  PLI=0
  for i in range(len(matrix)-2):
    for j in range(len(matrix)-1):
      for k in range(len(matrix)):
        if j>i:
          if k>j:
            PLI=PLI+ ( 
                (matrix[i][k]/(matrix[i][j]*matrix[j][k])) 
                + 
                ((matrix[i][j]*matrix[j][k])/matrix[i][k])
                -2)
  n=len(matrix)
  PLI=PLI*(6/(n*(n-1)*(n-2)))
  del matrix
  return PLI



def Tirads_geometric_consistency(matrix):
  TGCI=0
  n=len(matrix)
  for k in range(len(matrix)):
    for j in range(len(matrix)):
      for i in range(len(matrix)):
        if k>j:
          if j>i:
            TGCI=TGCI+(np.log(matrix[i][j]*matrix[j][k]*matrix[k][i]))*(np.log(matrix[i][j]*matrix[j][k]*matrix[k][i]))
  TGCI=2*TGCI
  TGCI=TGCI/(n*(n-1)*(n-2))
  return TGCI


def Koczkodaj_index(matrix):
  min_list=[]
  for i in range(len(matrix)):
    #print(i)
    for j in range(len(matrix)):
      #print(j)
      for k in range(len(matrix)):
        #print(k)
        min_a=abs(1-matrix[i][j]/(matrix[i][k]*matrix[k][j]))
        min_b=abs(1-(matrix[i][k]*matrix[k][j])/matrix[i][j])
        minimum=min(min_a,min_b)
        min_list.append(minimum)
  del matrix

  return(max(min_list))

####
####saving
####
def save_to_csv(matrix):
    
    filename= datetime.now()
    filename= str(filename.hour)+"-"+str(filename.minute)+"-"+str(filename.second)+"-"+str(filename.microsecond)+".csv"
    
    CR = str(return_CR(matrix))
    CR_str = "CR: "+ CR
    KI = str(Koczkodaj_index(matrix))
    KI_str = "KI: " + KI
    PLI = str(Pelaez_Lamata(matrix))
    PLI_str = "PLI: " + PLI
    TGCI = str(Tirads_geometric_consistency(matrix))
    TGCI_str = "T-GCI: "+ TGCI
    
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(matrix)
        writer.writerow([CR_str])
        writer.writerow([KI_str])
        writer.writerow([PLI_str])
        writer.writerow([TGCI_str])
    return filename,CR,KI,PLI,TGCI



def main_loop(size,loops):
    output_list=[]
    labels = ['filename','CR','KI','PLI','T-GCI']

    for n in range(loops):
        matrix = create_matrix(size)
        matrix = fill_matrix(matrix)
        output_list.append(save_to_csv(matrix))
        ##sleep because loop is too fast
        time.sleep(0.1)
        
    filename='output-n'+str(len(matrix))+'.csv'
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(labels)
        writer.writerows(output_list)
    return 0
	

#first parameter is size of matrix
#second parameter is number of generated matrices
main_loop(3,50)
