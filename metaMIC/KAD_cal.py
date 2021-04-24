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

base_path="/".join(sys.argv[0].split("/")[:-1])


def parseargs():
    parser=argparse.ArgumentParser(description="Calculate KAD values")
    parser.add_argument('--bam',help='index bam file for alignment')
    parser.add_argument("--contig",help="input assembly fasta file")
    parser.add_argument('--output',help='output directory for MetaREA')  
    parser.add_argument("--thread",type=int, default=8,help='length of k-mer')
    parser.add_argument("--samtools",default='samtools',help='path to samtools')
    parser.add_argument('--jellyfish',default='jellyfish',help='path to jellyfish')  
    parser.add_argument("--mlen",type=int, default=5000,help='minimum contig length [default: 5000bp]')
    args=parser.parse_args()
    return args

def contig_pool(samfile):
    contig_len={}
    for (ref,lens) in zip(samfile.references,samfile.lengths):
        contig_len[ref]=lens
    return contig_len        
        
def split_sam(args):
    split_command=' '.join(['sh',os.path.join(base_path,"split_sam.sh"),args.contig,args.bam,args.output,args.samtools])            
    os.system(split_command) 
    
def seq_parse(args):
    input=SeqIO.parse(args.contig,"fasta")
    contig_seqs={}
    for record in input:
        if len(record.seq) >= args.mlen:
            contig_seqs[record.id]=str(record.seq)
    return contig_seqs  
    
def kmer_parse(seq,pool):
    seq_kmer={"position":[],"KAD":[]}
    for i in range(len(seq)):
        if seq[i:(i+25)] in pool:
            seq_kmer["KAD"].append(pool[seq[i:(i+25)]])
            seq_kmer["position"].append(i+1)         
        if (i+25) >= len(seq):
            break
    return seq_kmer                        
    
def KAD_window_cal(seq_kmer):
    KAD_window_dict={"start_pos":[],"mean_KAD":[],"abnormal_KAD_ratio":[],"dev_KAD":[]}
    for i in range(300,len(seq_kmer['position']),100):
        KAD_window_dict["start_pos"].append(i)
        mean_KAD = np.mean(np.abs(seq_kmer['KAD'][i:i+100]))        
        KAD_window_dict["mean_KAD"].append(mean_KAD) 
        KAD_window_dict["abnormal_KAD_ratio"].append(np.sum(np.abs(seq_kmer['KAD'][i:i+100])>0.5)/100) 
        KAD_window_dict["dev_KAD"].append(np.sqrt(np.var(np.abs(seq_kmer['KAD'][i:i+100]))))               
    return KAD_window_dict                                 
    
def KAD_feature(args):
    seq_data=seq_parse(args)    
    KAD_dict = {"contig":[],'start_pos':[],'mean_KAD':[], 'abnormal_KAD_ratio':[], 'dev_KAD':[]}
    for contig,seq in seq_data.items():
        if len(seq) < args.mlen:
            continue
        if os.path.exists(args.output+"/temp/KAD/KAD_data/"+contig+".KAD"):
            try:
                KAD_data=pd.read_csv(args.output+"/temp/KAD/KAD_data/"+contig+".KAD",index_col=0,sep="\t")                             
                KAD_data=KAD_data.drop_duplicates(['k-mer'])     
            except:
                continue
            KAD_data.index = KAD_data['k-mer']
            KAD_pool=KAD_data.loc[:,'KAD'].to_dict()
            seq_kmer=kmer_parse(seq,KAD_pool)                              
            KAD_window=KAD_window_cal(seq_kmer)            
            KAD_dict["contig"].extend([contig]*len(KAD_window['start_pos']))
            KAD_dict["start_pos"].extend(KAD_window['start_pos'])
            KAD_dict["mean_KAD"].extend(KAD_window["mean_KAD"])  
            KAD_dict["abnormal_KAD_ratio"].extend(KAD_window["abnormal_KAD_ratio"])  
            KAD_dict["dev_KAD"].extend(KAD_window["dev_KAD"])  
    return KAD_dict                                     
    
                                                                
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
    kmer_result.loc[:,['k-mer','KAD']].to_csv(args.output+"/temp/KAD/KAD_data/"+str(contig)+".KAD",sep="\t")  

def KAD_cal(args):
    if os.path.exists(args.output+"/temp/KAD/KAD_window_data.txt"):
        return 0   
    split_sam(args)              
    contig_data=pd.read_csv(args.output+"/temp/split/contig_name.txt",header=None)
    split_data=pd.read_csv(args.output+"/temp/split/split_file_name.txt",header=None)
    data=pd.concat([contig_data,split_data],axis=1)
    data.columns=['contig','file']
    data.index=data['contig']
    contig_file=data.loc[:,'file'].to_dict()
    os.system("mkdir -p "+args.output+'/temp/KAD/temp')
    os.system("mkdir -p "+args.output+'/temp/KAD/KAD_data')  
    pool=multiprocessing.Pool(processes=args.thread)
    samfile=pysam.AlignmentFile(args.bam,"rb") 
    contig_len = contig_pool(samfile)
    for contig,file in contig_file.items():  
        if contig_len[contig] < args.mlen:
            continue
        try:
            t = pool.apply_async(func=KAD,args=(args,contig,file,))   
        except:
            continue            
    pool.close()
    pool.join()   
    KAD_dict = KAD_feature(args)
    KAD_window_data = pd.DataFrame(KAD_dict)                    
    KAD_window_data.to_csv(args.output+"/temp/KAD/KAD_window_data.txt",sep="\t")                                
                         
                                                   
def main():
    args=parseargs()
    warnings.filterwarnings("ignore")
    KAD_cal(args)
                            
if __name__=='__main__':
    main()  
