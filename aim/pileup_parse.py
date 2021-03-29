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
    parser.add_argument("--output",help="output directory for MetaREA results")
    parser.add_argument("--mlen",type=int, default=5000,help='minimum contig length [default: 5000bp]')    
    args = parser.parse_args()
    return args
    
def pileup_window_cal(pileup_dict):
    window_dict={"contig":[],"start_pos":[],"correct_portion":[],"ambiguous_portion":[],"disagree_portion":[],
    "deletion_portion":[],"insert_portion":[]}        
    for i in range(300,len(pileup_dict['position']),100):
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
        if len(pileup_dict['position']) - (i+100) <= 300:
                break
    return window_dict                        
            
def file_parse(args):
    """
    process pileup file
    """ 
    chunk=0
    prev_contig=None    
    pileup_dict={"contig":[],"position":[],"correct":[],"ambiguous":[],"insert":[],
        "deletion":[],"disagree":[],"depth":[]} 
    window_dict={"contig":[],"start_pos":[],"correct_portion":[],"ambiguous_portion":[],"disagree_portion":[],
    "deletion_portion":[],"insert_portion":[]} 
    for line in open(args.pileup,"r"):
        record = line.strip().split('\t')
        if prev_contig is None:
            prev_contig=record[0]
        if record[0] !=prev_contig:
            if len(pileup_dict['position']) < args.mlen:
                pileup_dict={"contig":[],"position":[],"correct":[],"ambiguous":[],"insert":[],
                "deletion":[],"disagree":[],"depth":[]}
            else:                                
                window_data=pileup_window_cal(pileup_dict) 
                window_dict["contig"].extend(window_data["contig"]) 
                window_dict["start_pos"].extend(window_data["start_pos"]) 
                window_dict["correct_portion"].extend(window_data["correct_portion"]) 
                window_dict["ambiguous_portion"].extend(window_data["ambiguous_portion"]) 
                window_dict["disagree_portion"].extend(window_data["disagree_portion"]) 
                window_dict["deletion_portion"].extend(window_data["deletion_portion"])   
                window_dict["insert_portion"].extend(window_data["insert_portion"])                                   
                pileup_dict={"contig":[],"position":[],"correct":[],"ambiguous":[],"insert":[],
                "deletion":[],"disagree":[],"depth":[]}
        pileup_dict['contig'].append(record[0])
        pileup_dict['position'].append(record[1])
        match_detail=record[4]
        pileup_dict['correct'].append(match_detail.count('.')+match_detail.count(','))
        pileup_dict['ambiguous'].append(match_detail.count('*'))
        pileup_dict['insert'].append(match_detail.count("+"))
        pileup_dict['deletion'].append(match_detail.count("-"))
        pileup_dict['depth'].append(int(record[3]))
        st = ''.join(re.split('[\+|\-][0-9]+[ATCGatcg]+',match_detail))
        m_list = np.array([st.count('a')+st.count('A'),st.count('t')+st.count('T'),st.count('c')+st.count('C'),st.count('g')+st.count('G')])
        pileup_dict['disagree'].append(np.sum(m_list))
        prev_contig=record[0]
    if not os.path.exists(args.output+"/temp/pileup"):
        os.system("mkdir -p "+ args.output+"/temp/pileup")
    data=pd.DataFrame(window_dict)                                                 
    data.to_csv(args.output+"/temp/pileup/pileup_feature.txt",sep="\t")  
    return data
                                                
def main():
    args=parseargs()
    warnings.filterwarnings("ignore")    
    data=file_parse(args)  
                                                 
if __name__=='__main__':
    main()   
