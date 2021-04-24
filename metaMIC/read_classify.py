#!/usr/bin/env python

import pandas as pd
import numpy as np
import multiprocessing
import argparse,operator,os,random,sys,time
import random,subprocess
import pysam
import collections
import warnings
import math

def parseargs():
    parser=argparse.ArgumentParser(description="Calculate read features")
    parser.add_argument('--bam',help='index bam file for alignment')
    parser.add_argument('--output',help='output directory for MetaREA')  
    parser.add_argument("--mlen",type=int, default=5000,help='minimum contig length [default: 5000bp]')  
    args=parser.parse_args()
    return args
             
def window_read_cal(samfile,ref,lens):
    read_count={"start_pos":[],"proper_count":[],"inversion_count":[],"clipped_count":[],"supplementary_count":[],"read_count":[]}
    for i in range(300,lens,100):
        start=i
        end=i+100
        read_count["start_pos"].append(i) 
        proper_num=0
        inversion_num=0
        clipped_num=0
        supplementary_num=0 
        readcount=0    
        for read in samfile.fetch(ref,start,end):
            proper_num+=(read.rnext==read.tid and read.is_proper_pair)
            inversion_num+=(read.rnext==read.tid and read.is_paired and ((read.is_reverse)+(read.mate_is_reverse)) != 1)
            clipped_num+=(read.get_cigar_stats()[0][4] > 20)
            supplementary_num+=(read.is_supplementary and read.get_cigar_stats()[0][5]>20) 
            readcount+=1                                    
        read_count["read_count"].append(readcount)
        read_count["proper_count"].append(proper_num)
        read_count["inversion_count"].append(inversion_num)
        read_count["clipped_count"].append(clipped_num)
        read_count["supplementary_count"].append(supplementary_num)
        if (lens-end) < 300:
            break
    return read_count                         
                 
def read_cal(args,samfile):
    if os.path.exists(args.output+"/temp/read_feature/read_feature.txt"):
        return 0            
    references=samfile.references
    lengths=samfile.lengths
    read_store={"contig":[],"start_pos":[],"read_count":[],"proper_read_count":[],
    "inversion_read_count":[],"clipped_read_count":[],"supplementary_read_count":[],"length":[]} 
    i=0
    for ref,lens in zip(references,lengths):
        print(i)
        i+=1
        if lens < args.mlen:
            continue 
        read_count=window_read_cal(samfile,ref,lens)     
        read_store["contig"].extend([ref]*len(read_count["start_pos"]))     
        read_store["start_pos"].extend(read_count["start_pos"]) 
        read_store["read_count"].extend(read_count["read_count"])
        read_store["proper_read_count"].extend(read_count["proper_count"])
        read_store["inversion_read_count"].extend(read_count["inversion_count"])
        read_store["clipped_read_count"].extend(read_count["clipped_count"])
        read_store["supplementary_read_count"].extend(read_count["supplementary_count"])   
        read_store["length"].extend([lens]*len(read_count["start_pos"]))
    data=pd.DataFrame(read_store)  
    data.to_csv(args.output+"/temp/read_feature/read_feature.txt",sep="\t")                                                                     

def main():
    args=parseargs()
    warnings.filterwarnings("ignore")
    if not os.path.isdir(args.output+"/temp/read_feature"):
        os.system("mkdir -p "+ args.output+"/temp/read_feature")  
    samfile=pysam.AlignmentFile(args.bam,"rb")   
    read_cal(args,samfile)          
                     
if __name__=="__main__":
    main()    
