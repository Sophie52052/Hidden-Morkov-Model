# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 09:24:39 2020

@author: Sophie
"""

import numpy as np
import math
from Bio import SeqIO 
from functools import reduce
import gzip
import time

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Data Preprocess
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def un_gz(file_name):
    """ungz zip file"""
    
    f_name = file_name.replace(".gz", "")
    g_file = gzip.GzipFile(file_name)
    open(f_name, "wb+").write(g_file.read())

print("GRCh38_latest_genomic.fna.gz 解壓縮中 ...")
un_gz("GRCh38_latest_genomic.fna.gz")


fin = open('GRCh38_latest_genomic.fna', 'r') 
print("讀取GRCh38...")
fout = open('data.txt', 'w') 
fout2 = open('data2.txt', 'w') 
item="NC_000006.12 Homo sapiens chromosome 6, GRCh38.p13 Primary Assembly"
print("正在尋找chromosome 6...")
for record in SeqIO.parse(fin,'fasta'): 
    #for item in my_list: 
     if item == record.description: 
      fout.write(">NC_000006.12:100001-200000 Homo sapiens chromosome 6, GRCh38.p13 Primary Assembly\n")
      fout2.write(">NC_000006.12:400001-40050 Homo sapiens chromosome 6, GRCh38.p13 Primary Assembly\n")
      #fout.write(">" + record.id + "\n")
      S= str(record.seq[100000:200000])
      S2=str(record.seq[300000:402000])
      fout.write(S.lower() + "\n") 
      fout2.write(S2.lower() + "\n") 

fin.close() 
fout.close()
fout2.close()


#test.txt
with open("data.txt", "r") as file:
    lines = file.readlines()
#print(lines)
lines.remove(lines[0])
space = ''
for line in lines:
    space += line.strip()#delet \n
per = space.split()#分割字串
seq = ''.join(per)#connect
#print(space)
#print(seq)
seq_list=[]
for i in range(len(seq)):
    per=seq[i]
    seq_list.append(per)
    
#test.txt
with open("data2.txt", "r") as file:
    lines2 = file.readlines()
#print(lines)
lines2.remove(lines2[0])
space = ''
for line in lines2:
    space += line.strip()#delet \n
per2 = space.split()#分割字串
seq2 = ''.join(per2)#connect
#print(space)
#print(seq)
seq_list2=[]
for i in range(len(seq2)):
    per=seq2[i]
    seq_list2.append(per)

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Order 0
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
from guppy import hpy
hp = hpy()

t0 = time.time()
total_base0=len(seq)

A_num=seq.count('a')
T_num=seq.count('t')
C_num=seq.count('c')
G_num=seq.count('g')
A_pc=A_num/total_base0
T_pc=T_num/total_base0
C_pc=C_num/total_base0
G_pc=G_num/total_base0

value_list = ["a", "t", "c", "g",]
probability = [A_pc, T_pc, C_pc, G_pc]
#origin
#order0_pro=math.pow(A_pc,A_num)*math.pow(T_pc,T_num)*math.pow(C_pc,C_num)*math.pow(G_pc,G_num)

#log base 2
order0_pro=math.log(A_pc,2)*A_num +math.log(T_pc,2)*T_num+math.log(C_pc,2)*C_num+math.log(G_pc,2)*G_num
t1 = time.time()
print("Markov Chain :")
print("Order 0 : " + str(order0_pro) +"  Time : "+str(t1-t0))

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Order 1
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
t0 = time.time()
#slide
order1=[]
total_base1=len(seq)-1 
for i in range(0,len(seq)-1):
    per=seq[i]+seq[i+1]
    order1.append(per)
    
#base_o1=["AA","AT","AC","AG","TA","TT","TC","TG",
#       "CA","CT","CC","CG","GA","GT","GC","GG"]
dict_o1 = {}
for key in order1:
    dict_o1[key] = dict_o1.get(key, 0) + 1/total_base1 #已是機率
#print (dict_o1["TT"])
    

    
states = ["A","T","C","G"]
base_o1=[["aa","at","ac","ag"],["ta","tt","tc","tg"],
         ["ca","ct","cc","cg"],["ga","gt","gc","gg"]]

base_o1m=[[dict_o1["aa"],dict_o1["at"],dict_o1["ac"],dict_o1["ag"]],
          [dict_o1["ta"],dict_o1["tt"],dict_o1["tc"],dict_o1["tg"]],
          [dict_o1["ca"],dict_o1["ct"],dict_o1["cc"],dict_o1["cg"]],
          [dict_o1["ga"],dict_o1["gt"],dict_o1["gc"],dict_o1["gg"]]]

base_o1m=[[(base_o1m[0][0])/(base_o1m[0][0]+base_o1m[0][1]+base_o1m[0][2]+base_o1m[0][3]),
             (base_o1m[0][1])/(base_o1m[0][0]+base_o1m[0][1]+base_o1m[0][2]+base_o1m[0][3]),
             (base_o1m[0][2])/(base_o1m[0][0]+base_o1m[0][1]+base_o1m[0][2]+base_o1m[0][3]),
             (base_o1m[0][3])/(base_o1m[0][0]+base_o1m[0][1]+base_o1m[0][2]+base_o1m[0][3])],
           [(base_o1m[1][0])/(base_o1m[1][0]+base_o1m[1][1]+base_o1m[1][2]+base_o1m[1][3]),
           (base_o1m[1][1])/(base_o1m[1][0]+base_o1m[1][1]+base_o1m[1][2]+base_o1m[1][3]),
           (base_o1m[1][2])/(base_o1m[1][0]+base_o1m[1][1]+base_o1m[1][2]+base_o1m[1][3]),
           (base_o1m[1][3])/(base_o1m[1][0]+base_o1m[1][1]+base_o1m[1][2]+base_o1m[1][3])],
          [(base_o1m[2][0])/(base_o1m[2][0]+base_o1m[2][1]+base_o1m[2][2]+base_o1m[2][3]),
           (base_o1m[2][1])/(base_o1m[2][0]+base_o1m[2][1]+base_o1m[2][2]+base_o1m[2][3]),
           (base_o1m[2][2])/(base_o1m[2][0]+base_o1m[2][1]+base_o1m[2][2]+base_o1m[2][3]),
           (base_o1m[2][3])/(base_o1m[2][0]+base_o1m[2][1]+base_o1m[2][2]+base_o1m[2][3])],
          [(base_o1m[3][0])/(base_o1m[3][0]+base_o1m[3][1]+base_o1m[3][2]+base_o1m[3][3]),
           (base_o1m[3][1])/(base_o1m[3][0]+base_o1m[3][1]+base_o1m[3][2]+base_o1m[3][3]),
           (base_o1m[3][2])/(base_o1m[3][0]+base_o1m[3][1]+base_o1m[3][2]+base_o1m[3][3]),
           (base_o1m[3][3])/(base_o1m[3][0]+base_o1m[3][1]+base_o1m[3][2]+base_o1m[3][3])]]

i = 0

if seq[0]=='a':
    prob=math.log(A_pc,2)
elif seq[0]=='t':
    prob=math.log(T_pc,2)
elif seq[0]=='c':
    prob=math.log(C_pc,2)
elif seq[0]=='g':
    prob=math.log(G_pc,2)


while i != total_base1:
    this_base=seq[i]
    change=seq[i]+seq[i+1]

    if this_base == "a":       
        if change == "aa":
            prob = prob + math.log(base_o1m[0][0],2)       
        elif change == "at":
            prob = prob + math.log(base_o1m[0][1],2)
        elif change == "ac":
            prob = prob + math.log(base_o1m[0][2],2)
        else:
            prob = prob + math.log(base_o1m[0][3],2)
    elif this_base == "t":      
        if change == "ta":
            prob = prob + math.log(base_o1m[1][0],2)       
        elif change == "tt":
            prob = prob + math.log(base_o1m[1][1],2)
        elif change == "tc":
            prob = prob + math.log(base_o1m[1][2],2)
        else:
            prob = prob + math.log(base_o1m[1][3],2)
    elif this_base == "c":      
        if change == "ca":
            prob = prob + math.log(base_o1m[2][0],2)    
        elif change == "ct":
            prob = prob + math.log(base_o1m[2][1],2)
        elif change == "cc":
            prob = prob + math.log(base_o1m[2][2],2)
        else:
            prob = prob + math.log(base_o1m[2][3],2)          
    elif this_base == "g":        
        if change == "ga":
            prob = prob + math.log(base_o1m[3][0],2)        
        elif change == "gt":
            prob = prob + math.log(base_o1m[3][1],2)         
        elif change == "gc":
            prob = prob + math.log(base_o1m[3][2],2)
        else:
            prob = prob + math.log(base_o1m[3][3],2)
           
    i += 1
t1 = time.time()
order1_pro  = prob
print("order 1 : " + str(order1_pro)+"  Time : "+str(t1-t0))

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Order 2
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
t0 = time.time()
order2=[]
total_base2=len(seq)-2 
for i in range(0,len(seq)-2):
    per=seq[i]+seq[i+1]+seq[i+2]
    order2.append(per)  

dict_o2 = {}
for key in order2:
    dict_o2[key] = dict_o2.get(key, 0) + 1/total_base2 #已是機率

base_o2=[[["aaa","aat","aac","aag"],["ata","att","atc","atg"],["aca","act","acc","acg"],["aga","agt","agc","agg"]],
          [["taa","tat","tac","tag"],["tta","ttt","ttc","ttg"],["tca","tct","tcc","tcg"],["tga","tgt","tgc","tgg"]],
          [["caa","cat","cac","cag"],["cta","ctt","ctc","ctg"],["cca","cct","ccc","ccg"],["cga","cgt","cgc","cgg"]],
          [["gaa","gat","gac","gag"],["gta","gtt","gtc","gtg"],["gca","gct","gcc","gcg"],["gga","ggt","ggc","ggg"]]]

base_o2m=[[[dict_o2["aaa"],dict_o2["aat"],dict_o2["aac"],dict_o2["aag"]],
           [dict_o2["ata"],dict_o2["att"],dict_o2["atc"],dict_o2["atg"]],
           [dict_o2["aca"],dict_o2["act"],dict_o2["acc"],dict_o2["acg"]],
           [dict_o2["aga"],dict_o2["agt"],dict_o2["agc"],dict_o2["agg"]]],
          [[dict_o2["taa"],dict_o2["tat"],dict_o2["tac"],dict_o2["tag"]],
           [dict_o2["tta"],dict_o2["ttt"],dict_o2["ttc"],dict_o2["ttg"]],
           [dict_o2["tca"],dict_o2["tct"],dict_o2["tcc"],dict_o2["tcg"]],
           [dict_o2["tga"],dict_o2["tgt"],dict_o2["tgc"],dict_o2["tgg"]]],
          [[dict_o2["caa"],dict_o2["cat"],dict_o2["cac"],dict_o2["cag"]],
           [dict_o2["cta"],dict_o2["ctt"],dict_o2["ctc"],dict_o2["ctg"]],
           [dict_o2["cca"],dict_o2["cct"],dict_o2["ccc"],dict_o2["ccg"]],
           [dict_o2["cga"],dict_o2["cgt"],dict_o2["cgc"],dict_o2["cgg"]]],
          [[dict_o2["gaa"],dict_o2["gat"],dict_o2["gac"],dict_o2["gag"]],
           [dict_o2["gta"],dict_o2["gtt"],dict_o2["gtc"],dict_o2["gtg"]],
           [dict_o2["gca"],dict_o2["gct"],dict_o2["gcc"],dict_o2["gcg"]],
           [dict_o2["gga"],dict_o2["ggt"],dict_o2["ggc"],dict_o2["ggg"]]]]

base_o1m=[[(base_o1m[0][0])/(base_o1m[0][0]+base_o1m[0][1]+base_o1m[0][2]+base_o1m[0][3]),
             (base_o1m[0][1])/(base_o1m[0][0]+base_o1m[0][1]+base_o1m[0][2]+base_o1m[0][3]),
             (base_o1m[0][2])/(base_o1m[0][0]+base_o1m[0][1]+base_o1m[0][2]+base_o1m[0][3]),
             (base_o1m[0][3])/(base_o1m[0][0]+base_o1m[0][1]+base_o1m[0][2]+base_o1m[0][3])],
           [(base_o1m[1][0])/(base_o1m[1][0]+base_o1m[1][1]+base_o1m[1][2]+base_o1m[1][3]),
           (base_o1m[1][1])/(base_o1m[1][0]+base_o1m[1][1]+base_o1m[1][2]+base_o1m[1][3]),
           (base_o1m[1][2])/(base_o1m[1][0]+base_o1m[1][1]+base_o1m[1][2]+base_o1m[1][3]),
           (base_o1m[1][3])/(base_o1m[1][0]+base_o1m[1][1]+base_o1m[1][2]+base_o1m[1][3])],
          [(base_o1m[2][0])/(base_o1m[2][0]+base_o1m[2][1]+base_o1m[2][2]+base_o1m[2][3]),
           (base_o1m[2][1])/(base_o1m[2][0]+base_o1m[2][1]+base_o1m[2][2]+base_o1m[2][3]),
           (base_o1m[2][2])/(base_o1m[2][0]+base_o1m[2][1]+base_o1m[2][2]+base_o1m[2][3]),
           (base_o1m[2][3])/(base_o1m[2][0]+base_o1m[2][1]+base_o1m[2][2]+base_o1m[2][3])],
          [(base_o1m[3][0])/(base_o1m[3][0]+base_o1m[3][1]+base_o1m[3][2]+base_o1m[3][3]),
           (base_o1m[3][1])/(base_o1m[3][0]+base_o1m[3][1]+base_o1m[3][2]+base_o1m[3][3]),
           (base_o1m[3][2])/(base_o1m[3][0]+base_o1m[3][1]+base_o1m[3][2]+base_o1m[3][3]),
           (base_o1m[3][3])/(base_o1m[3][0]+base_o1m[3][1]+base_o1m[3][2]+base_o1m[3][3])]]

          
i=0  
base_o2m_4=[[[]for i in range(4)] for i in range(4)]  
  
for i in range(len(base_o2m)):
    for j in range(len(base_o2m[i])) :
        total=0
        for k in range(len(base_o2m[i][j])) :
            total+=base_o2m[i][j][k]
        #print(total)
        for k in range(len(base_o2m[i][j])) :
            prob=base_o2m[i][j][k]/total
            base_o2m_4[i][j].append(prob)
            #print(prob)         
          


prob=0
this_base=seq[0]
change=seq[0]+seq[1]
#i=0,1
if this_base == "a":      
    if change == "aa":
        prob = math.log(A_pc,2) + math.log(base_o1m[0][0],2)       
    elif change == "at":
        prob = math.log(A_pc,2) + math.log(base_o1m[0][1],2)
    elif change == "ac":
        prob = math.log(A_pc,2) + math.log(base_o1m[0][2],2)
    else:
        prob = math.log(A_pc,2) + math.log(base_o1m[0][3],2)
elif this_base == "t":      
    if change == "ta":
        prob = math.log(T_pc,2) + math.log(base_o1m[1][0],2)       
    elif change == "tt":
        prob = math.log(T_pc,2) + math.log(base_o1m[1][1],2)
    elif change == "tc":
        prob = math.log(T_pc,2) + math.log(base_o1m[1][2],2)
    else:
        prob = math.log(T_pc,2) + math.log(base_o1m[1][3],2)
elif this_base == "c":      
    if change == "ca":
        prob = math.log(C_pc,2) + math.log(base_o1m[2][0],2)    
    elif change == "ct":
        prob = math.log(C_pc,2) + math.log(base_o1m[2][1],2)
    elif change == "cc":
        prob = math.log(C_pc,2) + math.log(base_o1m[2][2],2)
    else:
        prob = math.log(C_pc,2) + math.log(base_o1m[2][3],2)          
elif this_base == "g":        
    if change == "ga":
        prob = math.log(G_pc,2) + math.log(base_o1m[3][0],2)        
    elif change == "gt":
        prob = math.log(G_pc,2) + math.log(base_o1m[3][1],2)         
    elif change == "gc":
        prob = math.log(G_pc,2) + math.log(base_o1m[3][2],2)
    else:
        prob = math.log(G_pc,2) + math.log(base_o1m[3][3],2)


i = 0
while i != total_base2:
    this_base2=seq[i]
    this_base2_1=seq[i+1]
    change=seq[i+2]
    if this_base2 == "a":
        if this_base2_1 == "a":
            if change == "a":
                prob = prob + math.log(base_o2m_4[0][0][0],2)
            elif change == "t":
                prob = prob + math.log(base_o2m_4[0][0][1],2)
            elif change == "c":
                prob = prob + math.log(base_o2m_4[0][0][2],2)
                    
            else:
                prob = prob + math.log(base_o2m_4[0][0][3],2)
            
        elif this_base2_1 == "t":
            if change == "a":
                prob = prob + math.log(base_o2m_4[0][1][0],2)
            elif change == "t":
                prob = prob + math.log(base_o2m_4[0][1][1],2)
            elif change == "c":
                prob = prob + math.log(base_o2m_4[0][1][2],2)
            else:
                prob = prob + math.log(base_o2m_4[0][1][3],2)
            
        elif this_base2_1 == "c":
            if change == "a":
                prob = prob + math.log(base_o2m_4[0][2][0],2)
            elif change == "t":
                prob = prob + math.log(base_o2m_4[0][2][1],2)
            elif change == "c":
                prob = prob + math.log(base_o2m_4[0][2][2],2)
            else:
                prob = prob + math.log(base_o2m_4[0][2][3],2)         
            
        else:
            if change == "a":
                prob = prob + math.log(base_o2m_4[0][3][0],2)
            elif change == "t":
                prob = prob + math.log(base_o2m_4[0][3][1],2)
            elif change == "c":
                prob = prob + math.log(base_o2m_4[0][3][2],2)
            else:
                prob = prob + math.log(base_o2m_4[0][3][3],2)             

    elif this_base2 == "t":
        if this_base2_1 == "a":
            if change == "a":
                prob = prob + math.log(base_o2m_4[1][0][0],2)
            elif change == "t":
                prob = prob + math.log(base_o2m_4[1][0][1],2)
            elif change == "c":
                prob = prob + math.log(base_o2m_4[1][0][2],2)
                    
            else:
                prob = prob + math.log(base_o2m_4[1][0][3],2)
            
        elif this_base2_1 == "t":
            if change == "a":
                prob = prob + math.log(base_o2m_4[1][1][0],2)
            elif change == "t":
                prob = prob + math.log(base_o2m_4[1][1][1],2)
            elif change == "c":
                prob = prob + math.log(base_o2m_4[1][1][2],2)
            else:
                prob = prob + math.log(base_o2m_4[1][1][3],2)
            
        elif this_base2_1 == "c":
            if change == "a":
                prob = prob + math.log(base_o2m_4[1][2][0],2)
            elif change == "t":
                prob = prob + math.log(base_o2m_4[1][2][1],2)
            elif change == "c":
                prob = prob + math.log(base_o2m_4[1][2][2],2)
            else:
                prob = prob + math.log(base_o2m_4[1][2][3],2)           
            
        else:
            if change == "a":
                prob = prob + math.log(base_o2m_4[1][3][0],2)
            elif change == "t":
                prob = prob + math.log(base_o2m_4[1][3][1],2)
            elif change == "c":
                prob = prob + math.log(base_o2m_4[1][3][2],2)
            else:
                prob = prob + math.log(base_o2m_4[1][3][3],2)  

    elif this_base2 == "c":
        if this_base2_1 == "a":
            if change == "a":
                prob = prob + math.log(base_o2m_4[2][0][0],2)
            elif change == "t":
                prob = prob + math.log(base_o2m_4[2][0][1],2)
            elif change == "c":
                prob = prob + math.log(base_o2m_4[2][0][2],2)
                    
            else:
                prob = prob + math.log(base_o2m_4[2][0][3],2)
            
        elif this_base2_1 == "t":
            if change == "a":
                prob = prob + math.log(base_o2m_4[2][1][0],2)
            elif change == "t":
                prob = prob + math.log(base_o2m_4[2][1][1],2)
            elif change == "c":
                prob = prob + math.log(base_o2m_4[2][1][2],2)
            else:
                prob = prob + math.log(base_o2m_4[2][1][3],2)
            
        elif this_base2_1 == "c":
            if change == "a":
                prob = prob + math.log(base_o2m_4[2][2][0],2)
            elif change == "t":
                prob = prob + math.log(base_o2m_4[2][2][1],2)
            elif change == "c":
                prob = prob + math.log(base_o2m_4[2][2][2],2)
            else:
                prob = prob + math.log(base_o2m_4[2][2][3],2)          
            
        else:
            if change == "a":
                prob = prob + math.log(base_o2m_4[2][3][0],2)
            elif change == "t":
                prob = prob + math.log(base_o2m_4[2][3][1],2)
            elif change == "c":
                prob = prob + math.log(base_o2m_4[2][3][2],2)
            else:
                prob = prob + math.log(base_o2m_4[2][3][3],2)

    else :
        if this_base2_1 == "a":
            if change == "a":
                prob = prob + math.log(base_o2m_4[3][0][0],2)
            elif change == "t":
                prob = prob + math.log(base_o2m_4[3][0][1],2)
            elif change == "c":
                prob = prob + math.log(base_o2m_4[3][0][2],2)
                    
            else:
                prob = prob + math.log(base_o2m_4[3][0][3],2)
            
        elif this_base2_1 == "t":
            if change == "a":
                prob = prob + math.log(base_o2m_4[3][1][0],2)
            elif change == "t":
                prob = prob + math.log(base_o2m_4[3][1][1],2)
            elif change == "c":
                prob = prob + math.log(base_o2m_4[3][1][2],2)
            else:
                prob = prob + math.log(base_o2m_4[3][1][3],2)
            
        elif this_base2_1 == "c":
            if change == "a":
                prob = prob + math.log(base_o2m_4[3][2][0],2)
            elif change == "t":
                prob = prob + math.log(base_o2m_4[3][2][1],2)
            elif change == "c":
                prob = prob + math.log(base_o2m_4[3][2][2],2)
            else:
                prob = prob + math.log(base_o2m_4[3][2][3],2)         
            
        else:
            if change == "a":
                prob = prob + math.log(base_o2m_4[3][3][0],2)
            elif change == "t":
                prob = prob + math.log(base_o2m_4[3][3][1],2)
            elif change == "c":
                prob = prob + math.log(base_o2m_4[3][3][2],2)
            else:
                prob = prob + math.log(base_o2m_4[3][3][3],2)

    i += 1
#print(prob) #有少歐
    
t1 = time.time()      
order2_pro  = prob
print("order 2 : " + str(order2_pro)+"  Time : "+str(t1-t0))


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Hidden Morkov model
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
t0 = time.time()
num_seq=[]

for i in range(0, len(seq_list)):
    if(seq_list[i]=='a'):
        num_seq.append(0)
    elif(seq_list[i]=='t'):
        num_seq.append(1)
    elif(seq_list[i]=='c'):
        num_seq.append(2)
    elif(seq_list[i]=='g'):
        num_seq.append(3)
        
num_seq2=[]

for i in range(0, len(seq_list2)):
    if(seq_list2[i]=='a'):
        num_seq2.append(0)
    elif(seq_list2[i]=='t'):
        num_seq2.append(1)
    elif(seq_list2[i]=='c'):
        num_seq2.append(2)
    elif(seq_list2[i]=='g'):
        num_seq2.append(3)
        
class HMM:
    
    def __init__(self, transition_pro, emission_pro, initial_pro):
        self.transition_matix = np.array(transition_pro)     # 狀態轉移 NxN
        self.emissio_matix = np.array(emission_pro)     # 每個狀態產生各個base (N狀態 x M種base)
        self.initial_pro = np.array(initial_pro)   # 初始狀態，1xN
        self.state = self.transition_matix.shape[0]  #狀態
        self.per = self.emissio_matix.shape[1]  #(ATCG)

    def printresult(self):
        print("Initial Probability :\n", self.initial_pro)
        print("transition Probability :")
        for i in range(self.state):
            print(" ",self.transition_matix[i, :])
        print("emission_pro:")
        for i in range(self.state):
            print(" ",self.emissio_matix[i, :])

    def Forwardfun(self, seq_len, seq, forward, total):
        total[0] = 0.0

        # 第一個值
        for i in range(self.state):
            forward[0, i] = self.initial_pro[i] * self.emissio_matix[i, seq[0]] #計算alpha[index, i]計算在index個值從各個狀態產生的機率
            total[0] += forward[0, i] #total[index]是第index個值在各個狀態的機率總和
            
        #prob取log的另一個方法
        for i in range(self.state):
            forward[0, i] /= total[0]

        # 第二個以後
        for t in range(1, seq_len): #(1 ~ T-1)
            total[t] = 0.0
            for j in range(self.state): 
                p = 0.0  # p是前一個base在各個狀態到j狀態的所有機率
                for i in range(self.state):
                    p += forward[t-1, i] * self.transition_matix[i, j]
                forward[t, j] = p * self.emissio_matix[j, seq[t]]  # 現在base的機率總和 =上個base的所有機率*現在狀態產生這個base的機率
                total[t] += forward[t, j]
                
            for j in range(self.state):
                forward[t, j] /= total[t]

    def Backwardfun(self, seq_len, seq, backward, total):
        # 最後一個
        for i in range(self.state): # T-1
            backward[seq_len - 1, i] = 1.0

        # 除了最後一個以外的 (T-2 ~ 0)
        for t in range(seq_len - 2, -1, -1):
            for i in range(self.state):
                p2 = 0.0
                for j in range(self.state):
                    #現在這個base從狀態i轉到狀態j，產生t+1的base，再乘上backward
                    p2 += self.transition_matix[i, j] * self.emissio_matix[j, seq[t + 1]] * backward[t + 1, j]
                backward[t, i] = p2 / total[t + 1]

    # Learning 1 : gamma=((forward*backward)/total)
    # gamma[sequence, Nstate]
    def learning_fxb(self, seq_len, forward, backward, gamma):
        for t in range(seq_len): #0~T-1
            learning_fxb_total = 0.0 #總和
            for j in range(self.state):
                gamma[t, j] = forward[t, j] * backward[t, j]
                learning_fxb_total += gamma[t, j]
                
            for i in range(self.state):
                gamma[t, i] = gamma[t, i] / learning_fxb_total  

    # xi[sequence, Nstate, Nstate]
    def learning2_Xi(self, seq_len, seq, forward, backward, gamma, xi):
        for t in range(seq_len - 1): #0~T-2
            learning2_xi_total = 0.0
            for i in range(self.state):
                for j in range(self.state):
                    xi[t, i, j] = forward[t, i] * self.transition_matix[i, j] * self.emissio_matix[j, seq[t + 1]] * backward[t + 1, j]
                    learning2_xi_total += xi[t, i, j]

            for i in range(self.state):
                for j in range(self.state):
                    xi[t, i, j] = xi[t, i, j] / learning2_xi_total
        

    # Baum-Welch算法 ， HMM={A,B,initial_pro,N,M}
    def EM(self, seq_len, seq, forward, backward, gamma):
        
        run = 0 

        #  初始化
        xi = np.zeros((seq_len, self.state, self.state))           
        initial_pro = np.zeros((seq_len), np.float)                 
        total_gamma = np.zeros((self.state), np.float)                
        total_xi_state = np.zeros((self.state, self.state), np.float)      
        total_xi_state2base = np.zeros((self.state, self.per), np.float)      
        total = np.zeros((seq_len), np.float)               

        while True:
            # E_step
            self.Forwardfun(seq_len, seq, forward, total)
            self.Backwardfun(seq_len, seq, backward, total)
            self.learning_fxb(seq_len, forward, backward, gamma)
            self.learning2_Xi(seq_len, seq, forward, backward, gamma, xi)

            for i in range(self.state):
                initial_pro[i] += gamma[0, i]   
                for t in range(seq_len):
                    total_gamma[i] += gamma[t, i]  

                # total_xi_state:狀態轉狀態的機率和
                for j in range(self.state):
                    for t in range(seq_len - 1):
                        total_xi_state[i, j] += xi[t, i, j]  

                # total_xi_state2base:在i狀態產生鹼基k的機率和
                for k in range(self.per):
                    for t in range(seq_len):    
                        if seq[t] == k:  
                             total_xi_state2base[i, k] += gamma[t, i]  # i狀態產生k(base)的機率和

            # M_step 重新計算transition_pro(A)和emission_pro(B)
            for i in range(self.state):
                self.initial_pro[i] = initial_pro[i] 
                for j in range(self.state):
                    self.transition_matix[i, j] = total_xi_state[i, j] / total_gamma[i] #重新計算transition_pro
                    total_xi_state[i, j] = 0.0  #state i to j

                for k in range(self.per):
                    self.emissio_matix[i, k] = total_xi_state2base[i, k] / total_gamma[i] #重新計算emission_pro
                    total_xi_state2base[i, k] = 0.0  #state i to generate base(k)

                initial_pro[i] = total_gamma[i] = 0.0

            # Iteration
            if run==2:
                print("迭代次數：", run)
                break;
            else:
                run += 1
                continue
                
    def calculate(self, seq):
        t3 = time.time()
        state=np.array([0, 1, 2])
        max_prob = 0
        max_temp = -1

        # 第一個
        for j in state:
            if(max_temp < self.initial_pro[j] * self.emissio_matix[j, seq[0]]):
                max_temp = self.initial_pro[j] * self.emissio_matix[j, seq[0]]
                now_state = j
        max_prob += math.log(max_temp, 2)
   
        for i in range(1, len(seq)):
            max_temp = -1
            for k in state:
                if(max_temp < self.transition_matix[now_state, k] * self.emissio_matix[k, seq[i]]):
                    max_temp = self.transition_matix[now_state, k] * self.emissio_matix[k, seq[i]]
                    now_state = k
            max_prob += math.log(max_temp, 2)
            
      
        if len(seq) < 100002 :
            t1 = time.time()
            print("100001-200000 Maximun Probability : "+ str(max_prob) +"  Time : "+str(t1-t0))
        else :
            t2 = time.time()
            print("300001-402000 Maximun Probability : "+ str(max_prob)+"  Time : "+str(t2-t3))
        
        
        
if __name__ == "__main__":
    print("計算HMM中...")
    print("Hidden Markov Model")
    
    value_list = ["a", "t", "c", "g",]
    probability = [A_pc, T_pc, C_pc, G_pc]
    
    base_o1m=[[dict_o1["aa"],dict_o1["at"],dict_o1["ac"],dict_o1["ag"]],
          [dict_o1["ta"],dict_o1["tt"],dict_o1["tc"],dict_o1["tg"]],
          [dict_o1["ca"],dict_o1["ct"],dict_o1["cc"],dict_o1["cg"]],
          [dict_o1["ga"],dict_o1["gt"],dict_o1["gc"],dict_o1["gg"]]]
    
    State_int=[base_o1m[0][0]+base_o1m[0][1]+base_o1m[0][2]+base_o1m[0][3],
           base_o1m[1][0]+base_o1m[1][1]+base_o1m[1][2]+base_o1m[1][3],
           base_o1m[2][0]+base_o1m[2][1]+base_o1m[2][2]+base_o1m[2][3],
           base_o1m[3][0]+base_o1m[3][1]+base_o1m[3][2]+base_o1m[3][3]]
    
    base_o1m_4=[[(base_o1m[0][0])/(base_o1m[0][0]+base_o1m[0][1]+base_o1m[0][2]+base_o1m[0][3]),
             (base_o1m[0][1])/(base_o1m[0][0]+base_o1m[0][1]+base_o1m[0][2]+base_o1m[0][3]),
             (base_o1m[0][2])/(base_o1m[0][0]+base_o1m[0][1]+base_o1m[0][2]+base_o1m[0][3]),
             (base_o1m[0][3])/(base_o1m[0][0]+base_o1m[0][1]+base_o1m[0][2]+base_o1m[0][3])],
           [(base_o1m[1][0])/(base_o1m[1][0]+base_o1m[1][1]+base_o1m[1][2]+base_o1m[1][3]),
           (base_o1m[1][1])/(base_o1m[1][0]+base_o1m[1][1]+base_o1m[1][2]+base_o1m[1][3]),
           (base_o1m[1][2])/(base_o1m[1][0]+base_o1m[1][1]+base_o1m[1][2]+base_o1m[1][3]),
           (base_o1m[1][3])/(base_o1m[1][0]+base_o1m[1][1]+base_o1m[1][2]+base_o1m[1][3])],
          [(base_o1m[2][0])/(base_o1m[2][0]+base_o1m[2][1]+base_o1m[2][2]+base_o1m[2][3]),
           (base_o1m[2][1])/(base_o1m[2][0]+base_o1m[2][1]+base_o1m[2][2]+base_o1m[2][3]),
           (base_o1m[2][2])/(base_o1m[2][0]+base_o1m[2][1]+base_o1m[2][2]+base_o1m[2][3]),
           (base_o1m[2][3])/(base_o1m[2][0]+base_o1m[2][1]+base_o1m[2][2]+base_o1m[2][3])],
          [(base_o1m[3][0])/(base_o1m[3][0]+base_o1m[3][1]+base_o1m[3][2]+base_o1m[3][3]),
           (base_o1m[3][1])/(base_o1m[3][0]+base_o1m[3][1]+base_o1m[3][2]+base_o1m[3][3]),
           (base_o1m[3][2])/(base_o1m[3][0]+base_o1m[3][1]+base_o1m[3][2]+base_o1m[3][3]),
           (base_o1m[3][3])/(base_o1m[3][0]+base_o1m[3][1]+base_o1m[3][2]+base_o1m[3][3])]]
    
   
    transition_pro = [[base_o1m_4[0][0], base_o1m_4[0][1], base_o1m_4[0][2],base_o1m_4[0][3]],
                      [base_o1m_4[1][0], base_o1m_4[1][1], base_o1m_4[1][2],base_o1m_4[1][3]],
                      [base_o1m_4[2][0], base_o1m_4[2][1], base_o1m_4[2][2],base_o1m_4[2][3]],
                      [base_o1m_4[3][0], base_o1m_4[3][1], base_o1m_4[3][2],base_o1m_4[3][3]],]
    
    
    emission_pro = [[A_pc, T_pc, C_pc, G_pc], [A_pc, T_pc, C_pc, G_pc],[A_pc, T_pc, C_pc, G_pc],[A_pc, T_pc, C_pc, G_pc],]
    
    initial_pro = [State_int[0],State_int[1], State_int[2],State_int[3]]
    
    hmm = HMM(transition_pro, emission_pro, initial_pro)    # state=4； per=4

    seq = num_seq
    seq_len = len(seq)
    seq2 = num_seq2
    seq_len2 = len(seq2)
    forward = np.zeros((seq_len, hmm.state), np.float)   
    backward = np.zeros((seq_len, hmm.state), np.float)
    gamma = np.zeros((seq_len, hmm.state), np.float)
    hmm.EM(seq_len, seq, forward, backward, gamma)
    hmm.printresult()
    hmm.calculate(seq)
    
    hmm.calculate(seq2)
    
hp.heap()








