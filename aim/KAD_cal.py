#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse,operator,os,random,sys,time
import multiprocessing
import pysam
import re
import math
import collections
import warnings
from Bio import SeqIO

def parseargs():
    parser=argparse.ArgumentParser(description="Calculate KAD values")
    parser.add_argument('--bam',help='index bam file for alignment')
    parser.add_argument("--contig",help="input assembly fasta file")
    parser.add_argument('--output',help='output directory for MetaREA')  
    parser.add_argument("--thread",type=int, default=8,help='length of k-mer')
    parser.add_argument('--jellyfish',default='jellyfish',help='path to jellyfish')  
    parser.add_argument("--mlen",type=int, default=5000,help='minimum contig length [default: 5000bp]')
    args=parser.parse_args()
    return args
    
#args.bam='/home1/Laisenying/Data-analysis/phage/quality/train/temp/sam/contigs.filter.sort.bam'  
#args.contig='/home1/Laisenying/Data-analysis/phage/quality/train/temp/contig/filtered_contigs.fa'
#args.output='/home1/Laisenying/Data-analysis/phage/quality/train'
#args.jellyfish='jellyfish'  
    
def KAD(args,contig,file):
    if os.path.exists(args.output+"/temp/KAD/KAD_data/"+str(contig)+".KAD"):
        return 0
    contig_file=args.output+"/temp/split/contigs/"+file+".fa"
    read_file=args.output+"/temp/split/reads/"+"'"+str(contig)+".read.fa"+"'"
    # kmer count
    outputdir=args.output+"/temp/KAD/temp"
    contig_command1=' '.join([args.jellyfish,"count -m 25 -o",outputdir+"/"+"'"+str(contig)+".jf"+"'","-s 100M -t 8",contig_file])
    contig_command2=' '.join([args.jellyfish,"dump -c -t -o",outputdir+"/"+"'"+str(contig)+"_count.txt"+"'",outputdir+"/"+"'"+str(contig)+".jf"+"'"])
    os.system(contig_command1) 
    os.system(contig_command2)         
    read_command1=' '.join([args.jellyfish,"count -m 25 -o",outputdir+"/"+"'"+str(contig)+".read.jf"+"'","-s 100M -t 8",read_file])
    read_command2=' '.join([args.jellyfish,"dump -c -t -o",outputdir+"/"+"'"+str(contig)+"_count.read.txt"+"'",outputdir+"/"+"'"+str(contig)+".read.jf"+"'"])
    os.system(read_command1)  
    os.system(read_command2) 
    assembly_kmer=pd.read_csv(args.output+"/temp/KAD/temp/"+str(contig)+"_count.txt",sep="\t",header=None) 
    assembly_kmer.index=assembly_kmer[0]
    try:
        read_kmer=pd.read_csv(args.output+"/temp/KAD/temp/"+str(contig)+"_count.read.txt",sep="\t",header=None)   
        read_kmer.index=read_kmer[0]            
    except:
        # zero reads mapped to contig
        return 0  
    shared_kmer=set(assembly_kmer.loc[assembly_kmer[1]==1,0]).intersection(read_kmer.index)     
    if len(shared_kmer)==0: 
        kmer_depth=pd.value_counts(read_kmer.loc[read_kmer[1]>5,1]).index[0]
    else:                
        kmer_depth=pd.value_counts(read_kmer.loc[shared_kmer,][1]).index[0]        
    assembly_kmer.columns=['k-mer','assembly_count']
    read_kmer.columns=['k-mer','read_count']
    assembly_kmer.index=range(assembly_kmer.shape[0])
    read_kmer.index=range(read_kmer.shape[0])
    kmer_result=pd.merge(assembly_kmer,read_kmer,how='outer')
    kmer_result=kmer_result.fillna(0)
    kmer_result['KAD']=np.log2((kmer_result['read_count']+kmer_depth)/(kmer_depth*(kmer_result['assembly_count']+1)))
    kmer_result.loc[(kmer_result['read_count']==1)*(kmer_result['assembly_count']==0),'KAD']=np.nan
    kmer_result=kmer_result.loc[kmer_result['KAD']==kmer_result['KAD'],]     
    kmer_result.to_csv(args.output+"/temp/KAD/KAD_data/"+str(contig)+".KAD",sep="\t")  
    
def seq_parse(args):
    input=SeqIO.parse(args.contig,"fasta")
    contig_seqs={}
    for record in input:
        if len(record.seq) >= args.mlen:
            contig_seqs[record.id]=str(record.seq)
    return contig_seqs  
    
def kmer_parse(seq,data):
    seq_kmer={"position":[],"KAD":[]}
    for i in range(len(seq)):
        if seq[i:(i+25)] in data.index:
            seq_kmer["KAD"].append(data.loc[seq[i:(i+25)],"KAD"])
            seq_kmer["position"].append(i+1)         
        if (i+25)>=len(seq):
            break
    return seq_kmer                         

def KAD_window_cal(data):
    data=pd.DataFrame(data)
    data['start_pos']=[math.floor(x)*100+300 for x in (data['position']-300)/100] 
    data=data.loc[data['start_pos']>=300,]    
    data['abs_KAD']=np.abs(data['KAD'])  
    data['statu']=(data['abs_KAD']>0.5)+0  
    grouped=data.groupby(['start_pos'])
    KAD_mean=pd.DataFrame(grouped['abs_KAD'].mean())
    KAD_abnormal_count=pd.DataFrame(grouped['statu'].sum())                    
    KAD_dev=pd.DataFrame(np.sqrt(grouped['KAD'].var()))
    KAD_window_data=pd.concat([KAD_mean,KAD_abnormal_count/100,KAD_dev],axis=1)
    KAD_window_data.columns=['mean_KAD','abnormal_KAD_ratio','dev_KAD']
    KAD_window_data['start_pos']=[int(x) for x in KAD_window_data.index]
    KAD_window_data.index=range(KAD_window_data.shape[0])   
    return KAD_window_data              
    
def KAD_feature(args):
    seq_data=seq_parse(args)    
    KAD_dict = {"contig":[],'start_pos':[],'mean_KAD':[], 'abnormal_KAD_ratio':[], 'dev_KAD':[]}
    for contig,seq in seq_data.items():
        if os.path.exists(args.output+"/temp/KAD/KAD_data/"+contig+".KAD"):
            try:
                KAD_data=pd.read_csv(args.output+"/temp/KAD/KAD_data/"+contig+".KAD",index_col=0,sep="\t")                             
                KAD_data=KAD_data.drop_duplicates(['k-mer'])     
            except:
                continue
            KAD_data.index = KAD_data['k-mer']
            seq_kmer=kmer_parse(seq,KAD_data)      
            KAD_window=KAD_window_cal(seq_kmer)            
            KAD_dict["contig"].extend([contig]*KAD_window.shape[0])
            KAD_dict["start_pos"].extend(list(KAD_window["start_pos"]))
            KAD_dict["mean_KAD"].extend(list(KAD_window["mean_KAD"]))  
            KAD_dict["abnormal_KAD_ratio"].extend(list(KAD_window["abnormal_KAD_ratio"]))  
            KAD_dict["dev_KAD"].extend(list(KAD_window["dev_KAD"]))                           
    KAD_window_data = pd.DataFrame(KAD_dict)                    
    KAD_window_data.to_csv(args.output+"/temp/KAD/KAD_window_data.txt",sep="\t")
                                 
                                               
    
def main():
    args=parseargs()
    warnings.filterwarnings("ignore")
    contig_data=pd.read_csv(args.output+"/temp/split/contig_name.txt",header=None)
    split_data=pd.read_csv(args.output+"/temp/split/split_file_name.txt",header=None)
    data=pd.concat([contig_data,split_data],axis=1)
    data.columns=['contig','file']
    data.index=data['contig']
    os.system("mkdir -p "+args.output+'/temp/KAD/temp')
    os.system("mkdir -p "+args.output+'/temp/KAD/KAD_data')  
    pool=multiprocessing.Pool(processes=args.thread)
    for contig in list(np.unique(data['contig'])):  
        try:
            file=data.loc[contig,'file']
            t = pool.apply_async(func=KAD,args=(args,contig,file,))   
        except:
            continue            
    pool.close()
    pool.join()   
    KAD_feature(args) 
                        
if __name__=='__main__':
    main() 
