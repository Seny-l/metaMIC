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
                  
def fragment_distribution(samfile):
    all_reads=samfile.fetch()
    size_freq=collections.defaultdict(int)
    for read in all_reads:
        if read.rnext == read.tid and read.is_paired:
            size = abs(read.isize)
            size_freq[size]+=1                           
    return size_freq  
    
def FragMAD(freq):
    """
    calculate median and median absolute deviation fragment size distribution
    """
    all_size=[]
    for key,value in freq.items():
        all_size.extend([key]*int(value))
    median_size = np.median(all_size)
    residuals = abs(np.array(all_size)-median_size)
    mad_size = 1.4826*np.median(residuals)        
    return median_size,mad_size  
    
def fragment_coverage_cal(reads,mu,dev,length):
    """
    calculate fragment coverage per contig
    """      
    frag_coverage = np.array([0]*length)
    for read in reads:
        if read.rnext == read.tid and read.is_proper_pair:
            size = abs(read.isize)
            if (mu - 3*dev <= size <= mu + 3*dev):
                if read.next_reference_start < read.reference_start:
                    start=min(read.next_reference_start,read.reference_start,read.reference_end)
                    end = start + size                 
                    frag_coverage[start:end] += 1                                                         
    return frag_coverage                               
                                                                                                
def window_read_cal(reads,mu,dev):
    read_dict={"start_pos":[],"read_count":[],"proper_read_count":[],"inversion_read_count":[],"clipped_read_count":[],
    "supplementary_read_count":[],"discordant_size_count":[],"discordant_loc_count":[]}
    read_temp={"num_read":0,"num_proper":0,"num_inversion":0,"num_clipped":0,"num_supplementary":0,"num_discordant_size":0,
    "num_discordant_loc":0}
    pos=0
    for read in reads:
        new_pos=math.floor((read.reference_start-300)/100)*100+300
        if read.reference_start < 300:
            continue 
        if pos == 0:
            pos = new_pos 
        elif new_pos != pos:
            read_dict["start_pos"].append(pos)
            read_dict["read_count"].append(read_temp["num_read"])
            read_dict["proper_read_count"].append(read_temp["num_proper"])
            read_dict["inversion_read_count"].append(read_temp["num_inversion"])
            read_dict["clipped_read_count"].append(read_temp["num_clipped"])
            read_dict["supplementary_read_count"].append(read_temp["num_supplementary"])
            read_dict["discordant_size_count"].append(read_temp["num_discordant_size"])
            read_dict["discordant_loc_count"].append(read_temp["num_discordant_loc"])
            read_temp={"num_read":0,"num_proper":0,"num_inversion":0,"num_clipped":0,"num_supplementary":0,"num_discordant_size":0,"num_discordant_loc":0}
            pos = new_pos
        read_temp["num_read"] += 1  
        if read.is_paired:                                                                            
            if read.rnext==read.tid:
                if read.is_proper_pair: 
                    read_temp["num_proper"] += 1                                                                                                                                                                                                                                    
                if (read.is_reverse+read.mate_is_reverse) != 1:
                    read_temp["num_inversion"] += 1 
                if not mu-3*dev<= abs(read.isize) <= mu+3*dev:
                    read_temp["num_discordant_size"] += 1                                          
            else:
                read_temp["num_discordant_loc"] += 1                                                                                   
        if read.get_cigar_stats()[0][4] > 20:
            read_temp["num_clipped"] += 1
        if (read.is_supplementary and read.get_cigar_stats()[0][5]>20):
            read_temp["num_supplementary"] += 1
    return read_dict                                          
    
    
def window_frag_cal(coverage):
    """
    Using sliding window approach to smooth out features
    """               
    coverage=np.array(coverage)
    cov = {"pos":[],"coverage":[],"deviation":[]}
    for i in range(300,len(coverage),100):
        start=i
        end = i + 100
        cov["coverage"].append(np.mean(coverage[start:end]))
        cov["deviation"].append(np.sqrt(np.var(coverage[start:end]))/np.mean(coverage[start:end]))
        cov["pos"].append(start)
        if len(coverage) - end <= 300:
            break
    return cov     
                                                                                                                                        
def read_cal(args,mu,dev):
    if os.path.exists(args.output+"/temp/read_feature/read_feature.txt"):
        return 0            
    samfile=pysam.AlignmentFile(args.bam,"rb") 
    references=samfile.references
    lengths = samfile.lengths
    read_dicts={"contig":[],"start_pos":[],"read_count":[],"proper_read_count":[],"inversion_read_count":[],
    "clipped_read_count":[],"supplementary_read_count":[],"discordant_size_count":[],"discordant_loc_count":[],"length":[]}  
    for ref, lens in zip(references, lengths):
        if lens < args.mlen:
            continue
        contig_reads = samfile.fetch(ref)  
        read_dict=window_read_cal(contig_reads,mu,dev) 
        read_dicts["start_pos"].extend(read_dict["start_pos"])   
        read_dicts["contig"].extend([ref]*len(read_dict["start_pos"]))   
        read_dicts["read_count"].extend(read_dict["read_count"])
        read_dicts["proper_read_count"].extend(read_dict["proper_read_count"]) 
        read_dicts["inversion_read_count"].extend(read_dict["inversion_read_count"])
        read_dicts["clipped_read_count"].extend(read_dict["clipped_read_count"])
        read_dicts["supplementary_read_count"].extend(read_dict["supplementary_read_count"])  
        read_dicts["discordant_size_count"].extend(read_dict["discordant_size_count"])  
        read_dicts["discordant_loc_count"].extend(read_dict["discordant_loc_count"])  
        read_dicts["length"].extend([lens]*len(read_dict["start_pos"]))
    data=pd.DataFrame(read_dicts)  
    data.to_csv(args.output+"/temp/read_feature/read_feature.txt",sep="\t")  
    
def fragment_cal(args,mu,dev):
    if os.path.exists(args.output+"/temp/coverage/fragment_coverage.txt"):
        return 0
    samfile=pysam.AlignmentFile(args.bam,"rb") 
    references=samfile.references
    lengths = samfile.lengths
    frag_dict={"contig":[],"start_pos":[],"normalized_fragment_coverage":[],"normalized_fragment_deviation":[]}
    for ref, lens in zip(references,lengths):   
        if lens < args.mlen:
            continue
        reads = samfile.fetch(ref)
        frag_coverage=fragment_coverage_cal(reads,mu,dev,lens)  
        fragcov = window_frag_cal(frag_coverage)   
        frag_dict["contig"].extend([ref]*len(fragcov['pos']))       
        frag_dict["start_pos"].extend(fragcov["pos"])        
        frag_dict["normalized_fragment_coverage"].extend(fragcov["coverage"]/np.mean(fragcov["coverage"]))
        frag_dict["normalized_fragment_deviation"].extend(fragcov["deviation"])
    data=pd.DataFrame(frag_dict)                                 
    data.to_csv(args.output+"/temp/coverage/fragment_coverage.txt",sep="\t")
                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
def main():
    args=parseargs()
    warnings.filterwarnings("ignore")
    if not os.path.isdir(args.output+"/temp/read_feature"):
        os.system("mkdir -p "+ args.output+"/temp/read_feature")  
    if not os.path.isdir(args.output+"/temp/coverage"):
        os.system("mkdir -p "+ args.output+"/temp/coverage")          
    samfile=pysam.AlignmentFile(args.bam,"rb") 
    size_freq = fragment_distribution(samfile) 
    mu,dev = FragMAD(size_freq)          
    fragment_cal(args,mu,dev)   
    pool=[multiprocessing.Process(target=read_cal, args=(args,mu,dev,)),
    multiprocessing.Process(target=fragment_cal,args=(args,mu,dev,))]
    for t in pool:
        t.start()
    for t in pool:
        t.join()        
                                         
if __name__=="__main__":
    main() 
