#!/usr/bin/env python

import pandas as pd
import numpy as np
import math
import sys,os,argparse,warnings
import pysam
import collections
import re

def parseargs():
    parser = argparse.ArgumentParser(description="calculate pileup features")
    parser.add_argument("--pileup",help="path to pileup file")
    parser.add_argument('--bam',help='index bam file for alignment')
    parser.add_argument("--output",help="output directory for MetaREA results")
    parser.add_argument("--mlen",type=int, default=5000,help='minimum contig length [default: 5000bp]')
    args = parser.parse_args()
    return args
    
def contig_pool(samfile):
    contig_len={}
    for (ref,lens) in zip(samfile.references,samfile.lengths):
        contig_len[ref]=lens
    return contig_len        
         
def pileup_window_cal(pileup_dict):
    window_dict={"contig":[],"start_pos":[],"correct_portion":[],"ambiguous_portion":[],"disagree_portion":[],
    "deletion_portion":[],"insert_portion":[],"coverage":[],"deviation":[]}       
    for i in range(300,len(pileup_dict['correct']),100):
        start=i
        end=i+100
        total=np.sum(pileup_dict['depth'][start:end])
        window_dict["contig"].append(pileup_dict["contig"][0])
        window_dict["start_pos"].append(start)
        window_dict["correct_portion"].append(np.sum(pileup_dict['correct'][start:end])/total)
        window_dict["ambiguous_portion"].append(np.sum(pileup_dict["ambiguous"][start:end])/total)        
        window_dict["insert_portion"].append(np.sum(pileup_dict['insert'][start:end])/total)    
        window_dict["deletion_portion"].append(np.sum(pileup_dict['deletion'][start:end])/total)   
        window_dict["disagree_portion"].append(np.sum(pileup_dict['disagree'][start:end])/total) 
        window_dict["coverage"].append(np.mean(pileup_dict["depth"][start:end]))
        window_dict["deviation"].append(np.sqrt(np.var(pileup_dict["depth"][start:end]))/np.mean(pileup_dict["depth"][start:end]))        
        if len(pileup_dict['correct']) - (i+100) <= 300:
                break
    return window_dict                                                        
            
def pileupfile_parse(args,samfile):
    """
    process pileup file
    """    
    samfile=pysam.AlignmentFile(args.bam,"rb") 
    contig_len=contig_pool(samfile)
    prev_contig=None    
    pileup_dict={"contig":[],"correct":[],"ambiguous":[],"insert":[],
        "deletion":[],"disagree":[],"depth":[]} 
    window_pileup_dict={"contig":[],"start_pos":[],"correct_portion":[],"ambiguous_portion":[],"disagree_portion":[],
    "deletion_portion":[],"insert_portion":[],"normalized_coverage":[],"normalized_deviation":[],"mean_coverage":[]}  
    for line in open(args.pileup,"r"):
        record = line.strip().split('\t')
        if contig_len[record[0]] < args.mlen:
            continue            
        if prev_contig is None:
            prev_contig=record[0]          
        if record[0] !=prev_contig:  
            window_data=pileup_window_cal(pileup_dict) 
            mean_cov=np.mean(window_data["coverage"])
            window_pileup_dict["contig"].extend(window_data["contig"]) 
            window_pileup_dict["start_pos"].extend(window_data["start_pos"]) 
            window_pileup_dict["correct_portion"].extend(window_data["correct_portion"]) 
            window_pileup_dict["ambiguous_portion"].extend(window_data["ambiguous_portion"]) 
            window_pileup_dict["disagree_portion"].extend(window_data["disagree_portion"]) 
            window_pileup_dict["deletion_portion"].extend(window_data["deletion_portion"])   
            window_pileup_dict["insert_portion"].extend(window_data["insert_portion"])     
            window_pileup_dict["normalized_coverage"].extend(window_data["coverage"]/mean_cov)
            window_pileup_dict["normalized_deviation"].extend(window_data["deviation"])
            window_pileup_dict["mean_coverage"].extend([mean_cov]*len(window_data["start_pos"]))                                           
            pileup_dict={"contig":[],"correct":[],"ambiguous":[],"insert":[],"deletion":[],"disagree":[],"depth":[]}
            prev_contig = record[0]                
        pileup_dict['contig'].append(record[0])
        match_detail=record[4]        
        pileup_dict['correct'].append(match_detail.count('.')+match_detail.count(','))
        pileup_dict['ambiguous'].append(match_detail.count('*'))
        pileup_dict['insert'].append(match_detail.count("+"))
        pileup_dict['deletion'].append(match_detail.count("-"))
        pileup_dict['depth'].append(int(record[3]))
        st = ''.join(re.split('[\+|\-][0-9]+[ATCGatcg]+',match_detail))
        numd = st.count('a')+st.count('A')+st.count('t')+st.count('T')+st.count('c')+st.count('C')+st.count('g')+st.count('G')
        pileup_dict['disagree'].append(numd)
    if not os.path.exists(args.output+"/temp/pileup"):
        os.system("mkdir -p "+ args.output+"/temp/pileup")            
    data=pd.DataFrame(window_pileup_dict)                                                       
    data.to_csv(args.output+"/temp/pileup/pileup_feature.txt",sep="\t") 
    return data
                                                
def main():
    args=parseargs()
    warnings.filterwarnings("ignore")       
    data = pileupfile_parse(args)  
                                                         
if __name__=='__main__':
    main()    
