#!/usr/bin/env python

import pandas as pd
import numpy as np 
import math
import sys
import os
import argparse
import warnings
import pysam
import collections
import multiprocessing

def parseargs():
    parser=argparse.ArgumentParser(description="Calculate coverage/fragment coverage of assemblies")
    parser.add_argument('--bam',help='index bam file for alignment')
    parser.add_argument('--output',help='output directory for MetaREA results')  
    parser.add_argument('--bedtools',default='bedtools',help="path to bedtools")
    parser.add_argument("--mlen",type=int, default=5000,help='minimum contig length [default: 5000bp]')
    args=parser.parse_args()
    return args
        
def window_coverage_cal(args,coverage):
    """
    Using sliding window approach to smooth out features
    """               
    cov = {}
    moving_sum_array = []
    moving_dev_array = []
    pos  = []
    for i in range(300,len(coverage),100):
        start=i
        end = i + 100
        moving_sum_array.append(np.mean(coverage[start:end]))
        moving_dev_array.append(np.sqrt(np.var(coverage[start:end]))/np.mean(coverage[start:end]))
        pos.append(start)
        if len(coverage) - end <= 300:
            break
    cov["coverage"] = moving_sum_array
    cov["deviation"] = moving_dev_array      
    cov["pos"]=pos
    return cov                         
           
def depthparse(args):
    if not os.path.exists(os.path.join(args.output, "temp/coverage/contigs.depth")):
        command=args.bedtools + ' genomecov -ibam ' + args.bam + " -d > " + os.path.join(args.output, "temp/coverage/contigs.depth")
        os.system(command)
    file = open(os.path.join(args.output, "temp/coverage/contigs.depth"),"r")
    lines = file.readlines()      
    depth_dict={"position":[],"depth":[]}
    previous_contig=0   
    cov_dict={"contig":[],"start_pos":[],"normalized_coverage":[],"normalized_deviation":[],"mean_coverage":[]}            
    for line in lines:  
        contig=line.strip().split()[0]
        if previous_contig==0:
            previous_contig=contig
        if contig != previous_contig:
            if len(depth_dict["position"]) < args.mlen:
                depth_dict={"position":[],"depth":[]}    
            else:       
                coverage=list(depth_dict['depth'])  
                cov=window_coverage_cal(args,coverage) 
                mean_cov=mean_coverage(coverage)                                    
                cov_dict["normalized_coverage"].extend(cov['coverage']/mean_cov) 
                cov_dict["normalized_deviation"].extend(cov['deviation']) 
                cov_dict["mean_coverage"].extend([mean_cov]*len(cov['coverage']))
                cov_dict["contig"].extend([contig]*len(cov['coverage']))   
                cov_dict["start_pos"].extend(cov["pos"])  
                depth_dict={"position":[],"depth":[]}                              
        depth_dict["position"].append(int(line.strip().split()[1]))
        depth_dict["depth"].append(int(line.strip().split()[2]))   
        previous_contig=contig             
    data=pd.DataFrame(cov_dict)
    os.makedirs(os.path.join(args.output, "temp/coverage/"))
    data.to_csv(os.path.join(args.output, "temp/coverage/coverage.txt"),sep="\t")
    return data  
                  
def mean_coverage(coverage):
    """
    Using sliding window approach with 300bp window size to calculate mean coverage of contigs
    """
    moving_sum_array=[]
    for i in range(300,len(coverage),150):
        start=i
        end=start+300
        moving_sum_array.append(sum(coverage[start:end])/len(coverage[start:end]))
        if (len(coverage) - end) < 300:
            break
    return np.mean(moving_sum_array)            
                                                                                                                                                                                                                                                                                                                                                                                                      
def main():
    args=parseargs()
    warnings.filterwarnings("ignore")    
    samfile=pysam.AlignmentFile(args.bam,"rb")
    depthparse(args)       
        
if __name__=='__main__':
    main()    
