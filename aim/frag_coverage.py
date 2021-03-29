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
    parser.add_argument("--mlen",type=int, default=5000,help='minimum contig length')
    args=parser.parse_args()
    return args
    
def fragment_distribution(args,samfile):
    all_reads=samfile.fetch()
    references=samfile.references
    lengths=samfile.lengths
    size_freq=collections.defaultdict(int)
    pool={}
    contig_pool={}
    discordant_pool={}
    for ref,lens in zip(references,lengths):
        if lens < args.mlen:
            continue
        contig_pool[ref]=[]
        discordant_pool[ref]=[]
    for read in all_reads:
        contig = samfile.get_reference_name(read.tid) 
        if contig not in contig_pool:
            continue
        if read.rnext == read.tid:
            if read.qname in pool and pool[read.qname][0] == read.tid:
                mate_read = pool[read.qname]
                size=abs(max(mate_read[1] + mate_read[2] - read.pos, read.pos+read.rlen-mate_read[1]))
                size_freq[size]+=1
                contig_pool[contig].append([min(mate_read[1],read.pos),size])
            else:
                pool[read.qname] = (read.tid,read.pos,read.rlen)  
        else:
            contig1=samfile.get_reference_name(read.tid)
            contig2=samfile.get_reference_name(read.rnext)
            discordant_pool[contig1].append([read.pos,contig2])                
    return size_freq, contig_pool, discordant_pool   
    
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
    
def discordant_loc_count(args,discordant_pool):
    """
    Read pairs with mates mapped to different contigs
    """
    discs={"contig":[],"start_pos":[],"discordant_loc_count":[]}
    for key,value in discordant_pool.items():
        if len(discordant_pool[key]) == 0:
            continue
        subdata = pd.DataFrame(value)  
        subdata['start_pos']=[math.floor((x-300)/100)*100+300 for x in subdata[0]]       
        subdata=subdata.loc[subdata['start_pos']>=300,] 
        if subdata.shape[0]==0:
            continue    
        subdata['count']=1
        grouped=subdata.groupby(['start_pos'])
        data=pd.DataFrame(grouped[1].value_counts())     
        data['start_pos'] = [x[0] for x in data.index]         
        data['target'] = [x[1] for x in data.index] 
        data.index=range(data.shape[0])
        grouped=data.groupby(['start_pos'])     
        data=pd.DataFrame(grouped[1].max()) 
        positions=list(data.index)    
        discs["contig"].extend([key]*len(positions))
        discs["start_pos"].extend(positions)
        discs["discordant_loc_count"].extend(list(data[1]))
    data=pd.DataFrame(discs)
    os.system("mkdir -p "+ args.output+"/temp/read_feature/")
    data.to_csv(args.output+"/temp/read_feature/discordant_loc_feature.txt",sep='\t')
    
def discordant_size_count(args,contig_frag,mu,dev):
    """
    Read pairs with mates too for or too close to each other
    """
    discs={"contig":[],"start_pos":[],"discordant_size_count":[]}
    for key,value in contig_frag.items():
        if len(contig_frag[key])==0:
            continue
        subdata = pd.DataFrame(value)            
        subdata['start_pos']=[math.floor((x-300)/100)*100+300 for x in subdata[0]]
        subdata['count']=1
        subdata=subdata.loc[subdata['start_pos']>=300,] 
        if subdata.shape[0]==0:
            continue  
        subdata['statu'] = (subdata[1] < (mu - 3*dev)) + (subdata[1] > mu + 3*dev) 
        grouped = subdata.groupby(['start_pos'])    
        counts = grouped['statu'].sum()+0        
        discs["contig"].extend([key]*counts.shape[0])      
        discs["start_pos"].extend(list(counts.index))
        discs["discordant_size_count"].extend(list(counts))
    data = pd.DataFrame(discs)
    os.system("mkdir -p "+ args.output+"/temp/read_feature/")
    data.to_csv(args.output+"/temp/read_feature/discordant_size_feature.txt",sep='\t')
    
def fragment_coverage_per_contig(args,length,reads,mu,dev):
    """
    calculate fragment coverage per contig
    """    
    frag_coverage = [0]*length
    for read in reads:
        if read.cigarstring!=None:
            if read.rnext == read.tid:
                if (mu - 3*dev <= abs(read.isize) <= mu + 3*dev) and read.is_proper_pair and read.next_reference_start < read.reference_start:
                    if read.next_reference_start + read.rlen < read.reference_start:
                        for i in range(read.next_reference_start + read.rlen, read.reference_start):
                            frag_coverage[i]+=1
                    else:
                        for i in range(read.next_reference_start,read.reference_end):
                            frag_coverage[i]+=1
    frag_dict = window_coverage_cal(args,frag_coverage)
    return frag_dict     
    
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
                        
def fragment_coverage(args,samfile,mu,dev):
    """
    Calculate fragment coverage
    """   
    refs = samfile.references
    lengths = samfile.lengths
    norm_fragcov=[]
    norm_fragdev=[]
    position=[]
    contigs=[]
    for contig,length in zip(refs,lengths):   
        if length < args.mlen:
            continue       
        reads = samfile.fetch(contig,0,length)
        fragcov=fragment_coverage_per_contig(args,length,reads,mu,dev)
        contigs.extend([contig]*len(fragcov['pos']))
        norm_fragcov.extend(fragcov['coverage']/np.mean(fragcov['coverage']))
        norm_fragdev.extend(fragcov['deviation'])
        position.extend(fragcov['pos'])
    data=pd.DataFrame({"contig":contigs,"start_pos":position,
    "normalized_fragment_coverage":norm_fragcov,"normalized_fragment_deviation":norm_fragdev})    
    os.system("mkdir -p " + args.output+"/temp/coverage/")   
    data.to_csv(args.output+"/temp/coverage/fragment_coverage.txt",sep="\t")
               
def contig_index(data):
    """
    return the index of contigs
    """
    contig_to_pos={}
    data.index=range(data.shape[0])
    data['pos'] = data.index
    grouped = data.groupby(['contig'])  
    end_index=pd.DataFrame(grouped['pos'].max())     
    start_index=pd.DataFrame(grouped['pos'].min())
    for contig in start_index.index:
        start=start_index.loc[contig,'pos']
        end=end_index.loc[contig,"pos"]
        pos=range(start,end+1) 
        contig_to_pos[contig]=pos
    return contig_to_pos                    
                 
         
                                                                                                                                                                                                                                                                                                                                                                                        
def main():
    args=parseargs()
    warnings.filterwarnings("ignore")    
    samfile=pysam.AlignmentFile(args.bam,"rb")
    os.system("mkdir -p "+args.output+"/temp/coverage")
    size_freq,contig_pool,discordant_pool=fragment_distribution(args,samfile)
    mu,dev = FragMAD(size_freq)
    pool=[multiprocessing.Process(target=discordant_size_count, args=(args,contig_pool,mu,dev,)),
    multiprocessing.Process(target=discordant_loc_count,args=(args,discordant_pool,)),
    multiprocessing.Process(target=fragment_coverage,args=(args,samfile,mu,dev,))]
    for t in pool:
        t.start()
    for t in pool:
        t.join()        
        
if __name__=='__main__':
    main()
