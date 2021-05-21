#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse,operator,os,random,sys,time
import pysam
import re
import collections
import warnings
import math

def parseargs():
    parser=argparse.ArgumentParser(description="Calculate the number of read breakpoints per contig per base")
    parser.add_argument('--bam',help='index bam file for alignment')
    parser.add_argument('--output',help='output directory for MetaREA')  
    parser.add_argument("--mlen",type=int, default=5000,help='minimum contig length [default: 5000bp]')  
    args=parser.parse_args()
    return args
        
def read_breakpoint_per_contig(samfile,ref,lens):   
    reads=samfile.fetch(contig=ref)      
    break_count={"breakcount":np.array([0]*lens),"readcount":np.array([0]*lens)}
    for read in reads:
        ref_end=read.reference_end
        ref_start=read.reference_start
        read_start=read.query_alignment_start
        read_end=read.query_alignment_end
        break_count["readcount"][ref_start:ref_end] += 1
        if read.is_supplementary: 
            if re.match('^([0-9]+H)',read.cigarstring):
                break_count["breakcount"][read.get_blocks()[0][0]]+=1
            else:
                if len(read.get_blocks())==1:
                    break_count["breakcount"][read.get_blocks()[0][1]-1]+=1    
                else:
                    break_count["breakcount"][read.get_blocks()[-1][1]-1]+=1                     
        if read.get_cigar_stats()[0][4]>0:
            if re.match('^([0-9]+S)',read.cigarstring):
                break_count["breakcount"][read.get_blocks()[0][0]]+=1
            if (read.cigarstring).endswith('S'): 
                if len(read.get_blocks())==1:
                    break_count["breakcount"][read.get_blocks()[0][1]-1]+=1 
                else:                
                    break_count["breakcount"][read.get_blocks()[-1][1]-1]+=1    
    data=pd.DataFrame(break_count)      
    data['position']=data.index+1
    data['contig']=ref  
    data=data.loc[data['breakcount']>0,]
    return data
                                            
def window_break_cal(data):
    data['start_pos']=[math.floor(x)*100+300 for x in (data['position']-300)/100]
    data=data.loc[data['start_pos']>=300,]
    data['read_breakpoint_ratio']=data['read_breakpoint_count']/data['read_count']
    data['index']=data['contig']+'_'+[str(int(x)) for x in data['start_pos']]
    grouped=data.groupby(['index'])
    read_break_ratio=pd.DataFrame(grouped['read_breakpoint_ratio'].max()) 
    read_break_ratio['contig']=['_'.join(x.split("_")[:-1]) for x in read_break_ratio.index]
    read_break_ratio['start_pos']=[int(x.split("_")[-1]) for x in read_break_ratio.index]   
    read_break_ratio.index=range(read_break_ratio.shape[0])
    return read_break_ratio     
    
def read_breakpoint_cal(args):
    if os.path.exists(os.path.join(args.output, "temp/read_breakpoint/read_breakpoint_per_window.txt")):
        return 0
    if os.path.exists(os.path.join(args.output, "temp/read_breakpoint/read_breakpoint_per_base.txt")):
        read_breakpoint_data = pd.read_csv(os.path.join(args.output, "temp/read_breakpoint/read_breakpoint_per_base.txt"),sep="\t",index_col=0)
        window_read_breakpoint_data=window_break_cal(read_breakpoint_data)  
        window_read_breakpoint_data.to_csv(os.path.join(args.output, "temp/read_breakpoint/read_breakpoint_per_window.txt"),sep="\t")
        return 0                        
    samfile=pysam.AlignmentFile(args.bam,"rb") 
    references=samfile.references
    lengths=samfile.lengths
    read_breakpoint_pool={"contig":[],"position":[],"read_breakpoint_count":[],"read_count":[]}
    for ref,lens in zip(references,lengths):
        if lens < args.mlen:
            continue
        contig_break_data=read_breakpoint_per_contig(samfile,ref,lens)       
        if contig_break_data.shape[0]>0:
            read_breakpoint_pool["read_breakpoint_count"].extend(list(contig_break_data['breakcount']))
            read_breakpoint_pool["read_count"].extend(list(contig_break_data['readcount'])) 
            read_breakpoint_pool["contig"].extend([ref]*contig_break_data.shape[0])
            read_breakpoint_pool["position"].extend(list(contig_break_data['position'])) 
    read_breakpoint_data=pd.DataFrame(read_breakpoint_pool)    
    read_breakpoint_data.to_csv(os.path.join(args.output, "temp/read_breakpoint/read_breakpoint_per_base.txt"),sep="\t")
    window_read_breakpoint_data=window_break_cal(read_breakpoint_data)                        
    window_read_breakpoint_data.to_csv(os.path.join(args.output, "temp/read_breakpoint/read_breakpoint_per_window.txt"),sep="\t")
    
def main():
    args=parseargs()
    warnings.filterwarnings("ignore")
    if not os.path.isdir(os.path.join(args.output, "temp/read_breakpoint")):
        os.makedirs(os.path.join(args.output, "temp/read_breakpoint"))
    read_breakpoint_cal(args)  
                         
if __name__=="__main__":
    main()    
