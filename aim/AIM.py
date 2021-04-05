#!/usr/bin/env python

import pandas as pd
import numpy as np
import multiprocessing
from optparse import OptionParser  
import argparse, operator, os, random, sys, time
import random, subprocess
import pysam
import collections
import re
import math
from Bio import SeqIO
import joblib
import warnings
from sklearn.ensemble import IsolationForest
from scipy import stats
import logging

base_path="/".join(sys.argv[0].split("/")[:-1])

def get_opts():
    parser=OptionParser()
    parser.add_option("-1","--r1",dest="read1",help="paired-end #1 fasta/q files") 
    parser.add_option("-2","--r2",dest="read2",help="paired-end #2 fasta/q files") 
    parser.add_option("-p","--r",dest="read",help="smart pairing (ignoring #2 fasta/q)") 
    parser.add_option("-c","--contig",dest="assemblies",help="input assembly fasta file")
    parser.add_option("--bam",dest="bamfile",help="index bam file for alignment") 
    parser.add_option("-a","--assembler",dest="assembler",default="MEGAHIT",help="The assembler-specific model or user-trained model used for assembled fasta file [MEGAHIT/IDBA_UD/Train]") 
    parser.add_option("-o","--output",dest="output",help='output directory for AIM results')
    parser.add_option("-m","--mode", dest="mode",default='meta',help="Applied to single/metagenomic assemblies [meta/single]")
    parser.add_option("-t","--threads",dest="threads",type=int, default=8,help='Maximum number of threads [default: 8]')
    parser.add_option("-l","--mlen",dest="min_length",type=int, default=5000,help='Minimum contig length [default: 5000bp]') 
    parser.add_option("-s","--slen",dest="split_length",type=int, default=1000,help='Minimum length to be broken [default: 1000bp]') 
    parser.add_option("--pileup", dest="pileup",default='pileup',help="path to pileup file")
    parser.add_option("--samtools", dest="samtools",default='samtools',help="path to samtools")
    parser.add_option("--bwa", dest="bwa",default='bwa',help="path to bwa") 
    parser.add_option("--bedtools", dest="bedtools",default='bedtools',help="path to bedtools")
    parser.add_option("--jellyfish",dest="jellyfish",default="jellyfish",help="path to jellyfish")  
    parser.add_option("--train",dest="train",action="store_true",help='Training on user-specific datasets')  
    parser.add_option("--label",dest="label",help='Misassembled label for training assemblies') 
    parser.add_option("--no-breakpoints",dest="no_break",action="store_true",help='Do not locate possible breakpoints')    
    parser.add_option("--no-correct",dest="no_correct",action="store_true",help='Do not break misassembled contigs at breakpoints')
    (options, args) = parser.parse_args()
    return (options, args)  
       


def filter_contig(options):
    """
    remove assemblies with length smaller than required minimal length
    """
    input=SeqIO.parse(options.assemblies,"fasta")
    input=(record for record in input if len(record.seq)>1000)
    os.system("mkdir -p "+ options.output+"/temp/contig")
    SeqIO.write(input,options.output+"/temp/contig/filtered_contigs.fa","fasta")
    options.assemblies=options.output+"/temp/contig/filtered_contigs.fa"
    return options
    
def mapping(options):
    """
    mapping reads to assemblies
    """
    command_line=options.bwa + ' index ' + options.assemblies
    os.system(command_line)
    os.system("mkdir -p "+ options.output+"/temp/sam/")
    output_bam=options.output+"/temp/sam/contigs.filter.sort.bam"
    options.bamfile=output_bam
    if options.read is None:
        command_bwa=options.bwa+' mem -a -t ' + str(options.threads) + ' ' + options.assemblies + ' ' + options.read1 + ' ' + options.read2 + ' | ' + options.samtools + ' view -h -q 10 -m 50 -F 4 -b | '+ options.samtools + ' sort > ' + output_bam           
    else:
        command_bwa=options.bwa+' mem -a -p -t ' + str(options.threads) + ' ' + options.assemblies + ' ' + options.read + ' | ' + options.samtools + ' view -h -q 10 -m 50 -F 4 -b | '+ options.samtools + ' sort > ' + output_bam       
    os.system(command_bwa)
    command_index=options.samtools+' index '+output_bam
    os.system(command_index)  
    return options
    
def extract_feature(options):
    extract_command = ' '.join(['python',os.path.join(base_path,"extract.py"),'--bam',options.bamfile,'--contig',options.assemblies,'--output',options.output,'--mlen',str(options.min_length),'--pileup',options.pileup,'--samtools',options.samtools,'--jellyfish',options.jellyfish,"--thread",str(options.threads)]) 
    os.system(extract_command)                                   
    
    
def cov_thread_cal(data):
    up = np.percentile(np.array(data),95)
    low = np.percentile(np.array(data),5)
    return up,low
        

def contig_fea_generate(data):
    # transfer window based features to contig based feature
    grouped=data.groupby(['contig'])   
    mean_data=grouped.mean()   
    mean_data=mean_data.loc[:,mean_data.columns!='start_pos'] 
    data['count']=1
    data['proper_status']=data['proper_read_ratio']<=0.9
    data['clipped_status']=data['clipped_read_ratio']>=0.1
    data['supplementary_status']=data['supplementary_read_ratio']>=0.1
    data['inversion_status']=data['inversion_read_ratio']>=0.1  
    data['discordant_loc_status']=data['discordant_loc_ratio']>=0.1 
    data['discordant_size_status']=data['discordant_size_ratio']>=0.1  
    data['KAD_status']=data['mean_KAD']>0.75
    data['disagree_status']=data['disagree_portion']>0.1
    data['coverage_status']=(data['normalized_coverage']>cov_thread_cal(list(data['normalized_coverage']))[0])+(data['normalized_coverage']<cov_thread_cal(list(data['normalized_coverage']))[1])
    data['fragment_status']=(data['normalized_fragment_coverage']>cov_thread_cal(list(data['normalized_fragment_coverage']))[0])+(data['normalized_fragment_coverage']<cov_thread_cal(list(data['normalized_fragment_coverage']))[1])    
    data['deviation_status']=data['normalized_deviation']>cov_thread_cal(list(data['normalized_deviation']))[0]
    data['fragment_deviation_status']=data['normalized_fragment_deviation']>cov_thread_cal(list(data['normalized_fragment_deviation']))[0]    
    grouped=data.groupby(['contig'])   
    width_data=pd.concat([grouped['proper_status'].sum()/grouped['count'].sum(),grouped['clipped_status'].sum()/grouped['count'].sum(),grouped['supplementary_status'].sum()/grouped['count'].sum(),grouped['inversion_status'].sum()/grouped['count'].sum(),grouped['discordant_loc_status'].sum()/grouped['count'].sum(),grouped['discordant_size_status'].sum()/grouped['count'].sum(),grouped['KAD_status'].sum()/grouped['count'].sum(),grouped['disagree_status'].sum()/grouped['count'].sum(),grouped['coverage_status'].sum()/grouped['count'].sum(),grouped['fragment_status'].sum()/grouped['count'].sum(),grouped['deviation_status'].sum()/grouped['count'].sum(),grouped['fragment_deviation_status'].sum()/grouped['count'].sum()],axis=1)
    width_data.columns=['proper_read_width','clipped_read_width','supplementary_read_width','inversion_read_width','discordant_loc_width','discordant_size_width','KAD_width','disagree_width','coverage_width','fragment_width','deviation_width','fragment_deviation_width']
    dev_data=pd.concat([pd.DataFrame(np.sqrt(grouped['normalized_coverage'].var())),pd.DataFrame(np.sqrt(grouped['normalized_fragment_coverage'].var()))],axis=1)
    dev_data.columns=['window_cov_dev','window_frag_cov_dev']
    break_data=pd.DataFrame(grouped['read_breakpoint_ratio'].max())
    break_data.columns=['read_breakpoint_max']
    contig_data=pd.concat([mean_data,width_data,dev_data,break_data],axis=1) 
    features=['coverage_width', 'deviation_width','normalized_deviation', 'window_cov_dev', 'fragment_width', 'fragment_deviation_width',
    'normalized_fragment_deviation', 'window_frag_cov_dev', 'proper_read_ratio', 
    'clipped_read_ratio', 'supplementary_read_ratio', 'inversion_read_ratio', 'discordant_loc_ratio', 
    'discordant_size_ratio', 'read_breakpoint_ratio', 'proper_read_width', 'clipped_read_width', 
    'supplementary_read_width', 'inversion_read_width', 'discordant_loc_width', 'discordant_size_width',
     'read_breakpoint_max', 'disagree_width', 'correct_portion', 'ambiguous_portion', 'insert_portion', 
     'deletion_portion', 'disagree_portion', 'mean_KAD', 'abnormal_KAD_ratio', 'dev_KAD', 'KAD_width', 'coverage_diff','length'] 
    contig_data=contig_data.loc[:,features] 
    contig_data=contig_data.fillna(0)
    contig_data=contig_data.replace(np.inf,0)
    return contig_data
                 
                
def cal_feature(options):   
    read_feature = pd.read_csv(options.output+"/temp/read_feature/read_feature.txt",sep="\t",index_col=0)
    read_feature['proper_read_ratio'] = read_feature['proper_read_count']/read_feature['read_count']
    read_feature['inversion_read_ratio'] = read_feature['inversion_read_count']/read_feature['proper_read_count']
    read_feature['clipped_read_ratio'] = read_feature['clipped_read_count']/read_feature['proper_read_count']
    read_feature['supplementary_read_ratio'] = read_feature['supplementary_read_count']/read_feature['proper_read_count']
    read_feature['discordant_loc_ratio'] = read_feature['discordant_loc_count']/read_feature['proper_read_count']
    read_feature['discordant_size_ratio'] = read_feature['discordant_size_count']/read_feature['proper_read_count']
    window_read_data = read_feature.loc[:,['contig','start_pos','proper_read_ratio','inversion_read_ratio','clipped_read_ratio','supplementary_read_ratio','discordant_loc_ratio','discordant_size_ratio','length']]    
    frag_coverage = pd.read_csv(options.output+"/temp/coverage/fragment_coverage.txt",sep="\t",index_col=0)    
    pileup_feature = pd.read_csv(options.output+"/temp/pileup/pileup_feature.txt",sep="\t",index_col=0)
    KAD_feature = pd.read_csv(options.output+"/temp/KAD/KAD_window_data.txt",sep="\t",index_col=0)            
    breakpoint_data = pd.read_csv(options.output+"/temp/read_breakpoint/read_breakpoint_per_window.txt",sep="\t",index_col=0)
    window_data = pd.merge(window_read_data,frag_coverage,on=['contig','start_pos'])
    window_data = pd.merge(window_data,pileup_feature,on=['contig','start_pos'])
    window_data = pd.merge(window_data,KAD_feature,on=['contig','start_pos'])
    window_data = pd.merge(window_data,breakpoint_data,on=['contig','start_pos'],how='left')
    window_data = window_data.fillna(0)    
    window_data['coverage_diff'] = window_data['normalized_coverage'] - window_data['normalized_fragment_coverage']
    os.system("mkdir -p "+options.output+"/feature_matrix")    
    window_data = window_data.loc[window_data['mean_coverage']>5,]
    window_data.to_csv(options.output+"/feature_matrix/window_fea_matrix.txt",sep="\t")
    contig_data = contig_fea_generate(window_data)          
    contig_data.to_csv(options.output+"/feature_matrix/contig_fea_matrix.txt",sep="\t")      
    return contig_data 
    
def check_feature(options):
    if not os.path.getsize(options.output+"/temp/read_feature/read_feature.txt"):
        logging.error("Can not find file: " + options.output+"/temp/read_feature/read_feature.txt")
    if not os.path.getsize(options.output+"/temp/coverage/fragment_coverage.txt"):
        logging.error("Can not find file: " + options.output+"/temp/coverage/fragment_coverage.txt")    
    if not os.path.getsize(options.output+"/temp/read_breakpoint/read_breakpoint_per_window.txt"): 
        logging.error("Can not find file: " + options.output+"/temp/read_breakpoint/read_breakpoint_per_window.txt")
    if not os.path.getsize(options.output+"/temp/pileup/pileup_feature.txt"): 
        logging.error("Can not find file: " + options.output+"/temp/pileup/pileup_feature.txt")        
    if not os.path.getsize(options.output+"/temp/KAD/KAD_window_data.txt"):                                          
        logging.error("Can not find file: " + options.output+"/temp/KAD/KAD_window_data.txt")       

def predict(options,data):
    # Generate MetaREA contig scores
    min_length=options.min_length
    test_data = data.loc[data['length']>min_length,data.columns!='length']
    if options.assembler=='MEGAHIT':
        model_path=os.path.join(base_path,"dataset/MEGAHIT/model6")
    elif options.assembler=='IDBA_UD':
        model_path=os.path.join(base_path,"model/IDBA_UD") 
    elif options.assembler=='Train':
        model_path = os.path.join(base_path,"train_model")                
    else:
        model_path=os.path.join(base_path,"model/merge")        
    score=pd.DataFrame(np.zeros([test_data.shape[0],10])) 
    score.index=test_data.index   
    for i in range(10):                                                                                  
        rf=joblib.load(model_path+"/RF"+str(i)+'.pkl')
        pro=pd.DataFrame(rf.predict_proba(test_data))  
        pro.index=test_data.index
        score.loc[pro.index,i] = pro[1]
    score['aim_contig_score']=score.mean(axis=1)
    score=score.loc[:,['aim_contig_score']]
    score['length']=data.loc[score.index,'length']  
    score.to_csv(options.output+"/aim_contig_score.txt",sep="\t")    
    return score  

def findcut(score):
    # find score cutoff to identify candidate misassembled contigs
    score_cut=np.percentile(score['aim_contig_score'],98)
    return score_cut

def contig_pos(data):
    contig_to_pos = {}
    data.index=range(data.shape[0])
    data['pos']=data.index
    grouped=data.groupby(['contig'])
    contig_end=pd.DataFrame(grouped['pos'].max())
    contig_start=pd.DataFrame(grouped['pos'].min())
    for contig in contig_start.index:
        start = contig_start.loc[contig,'pos']
        end = contig_end.loc[contig,'pos']
        poss = range(start,end+1)
        contig_to_pos[contig] = poss
    return contig_to_pos       
    
def Isolation_forest(options,window_data):
    features=['correct_portion', 'ambiguous_portion',
       'disagree_portion', 'deletion_portion', 'insert_portion', 'mean_KAD',
       'abnormal_KAD_ratio', 'dev_KAD', 'normalized_fragment_coverage',
       'normalized_fragment_deviation', 'normalized_coverage',
       'normalized_deviation', 'mean_coverage', 'coverage_diff',
       'read_breakpoint_ratio', 'proper_read_ratio', 'inversion_read_ratio',
       'clipped_read_ratio', 'supplementary_read_ratio',
       'discordant_loc_ratio', 'discordant_size_ratio']
    Xdata=window_data.loc[:,features]    
    n_samples=Xdata.shape[0]    
    if options.mode=='meta':
        iso_parameter={"outlier_fraction":0.001,"outlier_threshold":0.01}
    else:
        iso_parameter={"outlier_fraction":0.0001,"outlier_threshold":0.001}                  
    # fit the model
    clf=IsolationForest(contamination=iso_parameter["outlier_fraction"])
    clf.fit(Xdata)
    score_pred=clf.decision_function(Xdata)
    Xdata['outlier_score']=-score_pred
    threshold=stats.scoreatpercentile(-score_pred,100*(1-iso_parameter["outlier_threshold"]))                
    Xdata['outlier_thred']=threshold
    score_pred_data=Xdata.loc[:,['outlier_score','outlier_thred','read_breakpoint_ratio']]
    score_pred_data['start_pos']=[int(x.split("_")[-1]) for x in score_pred_data.index]
    score_pred_data['contig']=['_'.join(x.split("_")[:-1]) for x in score_pred_data.index]
    return score_pred_data                      
        
def meta_breakpoint_detect(options,filter_score):
    # identify misassembly breakpoints on metagenomic misassembled contigs
    if not os.path.exists(options.output+"/feature_matrix/window_fea_matrix.txt"):
        logging.error("Missing window feature matrix file:" + options.output+"/feature_matrix/window_fea_matrix.txt")
        exit(-1)
    if not os.path.exists(options.output+"/temp/read_breakpoint/read_breakpoint_per_base.txt"):
        logging.error("Missing read breakpoint file:"+ options.output+"/temp/read_breakpoint/read_breakpoint_per_base.txt")
        exit(-1)            
    window_data=pd.read_csv(options.output+"/feature_matrix/window_fea_matrix.txt",sep="\t",index_col=0)
    window_data.index=window_data['contig']
    window_data=window_data.loc[set(filter_score.index).intersection(window_data.index),]      
    window_data=window_data.fillna(0)
    window_data=window_data.replace(np.inf,0)
    window_data.index=window_data['contig']+"_"+[str(int(x)) for x in window_data['start_pos']]       
    score_pred_data = Isolation_forest(options,window_data)     
    score_pred_data.to_csv(options.output+"/outlier_score.txt",sep="\t")            
    contig_to_pos=contig_pos(score_pred_data) 
    outlier_data=pd.DataFrame(columns=['contig','start_pos','outlier_score','outlier_thred','read_breakpoint_ratio'])
    for contig in np.unique(score_pred_data['contig']):
        contig_outlier_score=score_pred_data.iloc[contig_to_pos[contig],]
        contig_outlier_score=contig_outlier_score.iloc[np.argsort(-contig_outlier_score['outlier_score']),]
        contig_outlier_score=contig_outlier_score.iloc[:1,]
        contig_outlier_score=contig_outlier_score.loc[:,['contig','start_pos','outlier_score','outlier_thred','read_breakpoint_ratio']]
        outlier_data=pd.concat([outlier_data,contig_outlier_score])
    read_breakpoint=pd.read_csv(options.output+"/temp/read_breakpoint/read_breakpoint_per_base.txt",sep="\t",index_col=0)
    read_breakpoint['start_pos']=[int(x)*100+300 for x in list((read_breakpoint['position']-300)/100)]
    result = pd.merge(read_breakpoint,outlier_data,on=['contig','start_pos'],how='right')
    result['read_breakpoint_count']=result['read_breakpoint_count'].fillna(0)
    result['position']=result['position'].fillna(0)
    result=result.iloc[np.argsort(-result['read_breakpoint_count']),]
    result=result.drop_duplicates(['contig'])
    result.loc[result['position']==0,'position']=result.loc[result['position']==0,'start_pos']+50
    result['length']=list(filter_score.loc[result['contig'],'length'])
    result['aim_contig_score']=list(filter_score.loc[result['contig'],'aim_contig_score'])
    breakpoint_result=result.loc[:,['contig', 'position', 'read_breakpoint_count', 'read_count',
        'read_breakpoint_ratio', 'outlier_score',"outlier_thred",'aim_contig_score',"length"]]
    breakpoint_result.columns=['contig', 'misassembly_breakpoint', 'read_breakpoint_count', 'read_count',
       'read_breakpoint_ratio','outlier_score',"outlier_thred",'aim_contig_score',"contig_length"]       
    breakpoint_result.to_csv(options.output+"/mis_breakpoint.txt",sep="\t")   
    return breakpoint_result     
        
def single_breakpoint_detect(options):
    # identify misassembly breakpoints on single genome assemblies
    if not os.path.exists(options.output+"/feature_matrix/window_fea_matrix.txt"):
        logging.error("Missing window feature matrix file:" + options.output+"/feature_matrix/window_fea_matrix.txt")
        exit(-1)
    if not os.path.exists(options.output+"/temp/read_breakpoint/read_breakpoint_per_base.txt"):
        logging.error("Missing read breakpoint file: " + options.output+"/temp/read_breakpoint/read_breakpoint_per_base.txt")
        exit(-1)             
    window_data=pd.read_csv(options.output+"/feature_matrix/window_fea_matrix.txt",sep="\t",index_col=0)
    window_data=window_data.fillna(0)
    window_data=window_data.replace(np.inf,0)
    window_data.index=window_data['contig']+"_"+[str(int(x)) for x in window_data['start_pos']]                    
    score_pred_data = Isolation_forest(options,window_data)    
    score_pred_data.to_csv(options.output+"/outlier_score.txt")                     
    read_breakpoint=pd.read_csv(options.output+"/temp/read_breakpoint/read_breakpoint_per_base.txt",sep="\t",index_col=0)
    read_breakpoint['start_pos']=[int(x)*100+300 for x in list((read_breakpoint['position']-300)/100)]
    score_pred_data=score_pred_data.loc[score_pred_data['outlier_score']>score_pred_data['outlier_thred'],]
    score_pred_data.index=range(score_pred_data.shape[0])
    result = pd.merge(read_breakpoint,score_pred_data,on=['contig','start_pos'])
    result = result.iloc[np.argsort(-result['read_breakpoint_count']),]
    result['id']=result['contig']+"_"+[str(int(x)) for x in result['start_pos']]
    result = result.drop_duplicates(['id'])
    result = result.loc[result['read_breakpoint_count']/result['read_count']>0.2,]
    result = result.loc[result['read_breakpoint_count']>5,]
    contig_len=window_data.loc[:,['contig','length']].drop_duplicates()
    contig_len.index=contig_len['contig']
    result['contig_length']=list(contig_len.loc[result['contig'],'length'])
    breakpoint_result = result.loc[:,["contig","position","read_breakpoint_count","read_count","read_breakpoint_ratio","outlier_score","outlier_thred","contig_length"]]
    breakpoint_result.columns=["contig","misassembly_breakpoint","read_breakpoint_count","read_count","read_breakpoint_ratio","outlier_score","outlier_thred","contig_length"]
    return breakpoint_result
                                        
def corrected_contig(options,breakpoint_result):
    breakpoint_result=breakpoint_result.loc[breakpoint_result['misassembly_breakpoint']>options.split_length,]  
    breakpoint_result=breakpoint_result.loc[(breakpoint_result['contig_length']-breakpoint_result['misassembly_breakpoint'])>options.split_length,]
    breakpoint_result=breakpoint_result.loc[breakpoint_result['read_breakpoint_count']>5,]
    breakpoint_result=breakpoint_result.loc[breakpoint_result['read_breakpoint_ratio']>0.1,]
    corrected_contig_file=options.output+"/aim_corrected_contigs.fa"
    corrected_file=open(corrected_contig_file,"w")
    original_file=options.assemblies 
    input=SeqIO.parse(original_file,"fasta")
    breakcontigs=list(np.unique(breakpoint_result['contig']))
    for record in input:
        if record.id in breakcontigs:
            corrected_file.write(">"+record.id+"_1\n")
            breakpoint=int(list(breakpoint_result.loc[breakpoint_result['contig']==record.id,'misassembly_breakpoint'])[0])
            corrected_file.write(str(record.seq[:breakpoint])+"\n")
            corrected_file.write(">"+record.id+"_2\n")
            corrected_file.write(str(record.seq[breakpoint:])+"\n")   
        else:
            corrected_file.write(">"+record.id+"\n")   
            corrected_file.write(str(record.seq)+"\n")  
    corrected_file.close()  
    return breakpoint_result     
    
def train_model(options,data):
    if not os.path.exists(options.label):
        logging.error("Please provide train labels for contigs [0: no misassembly,1: misassembly]")
        exit(-1)    
    train_commond = ' '.join(['python3',os.path.join(base_path,"train.py"),'--data',options.output+"/feature_matrix/contig_fea_matrix.txt",'--label',options.label,"--thread",str(options.threads)])                            
                                                                                                                                
                                                                                                                                                                                                                                                      
def pipeline(options):
    ##### Feature extraction #########
    print("Step1: Feature extraction")
    extract_feature(options)   
    check_feature(options)       
    contig_data=cal_feature(options)    
    # Identify misassembled contigs 
    if options.train:
        # Training on user provided datasets
        print("Training models")
        localtime = time.asctime(time.localtime(time.time()))
        print("Time:"+localtime)
        train_model(options,contig_data)
        print("Finished training")               
        return 0                                                                                
    if options.mode=='meta':
        # Applicaton of AMI to metagenomics
        print("Step2: Identify misassembled contigs")
        localtime = time.asctime(time.localtime(time.time()))
        print("Time:"+localtime)        
        score=predict(options,contig_data)                           
        # Misassembly breakpoint identification 
        if options.no_break:
            print("Finished")            
            return 0
        print("Step3: Misassembly Breakpoint Identification") 
        score_cut=findcut(score)
        filter_score=score.loc[score['aim_contig_score']>score_cut,]
        breakpoint_result = meta_breakpoint_detect(options,filter_score)
        localtime = time.asctime(time.localtime(time.time()))
        print("Time:"+localtime)        
        # Break misassembled contigs at misassembled breakpoints 
        if options.no_correct:
            print("Finished")            
            return 0
        print("Step4: Correcting misassembled contigs...")    
        localtime = time.asctime(time.localtime(time.time()))                           
        correct_result=corrected_contig(options,breakpoint_result)
        print("A total of " + str(correct_result.shape[0]) + " misassembled contigs are corrected")
        localtime = time.asctime(time.localtime(time.time()))
        print("Time:"+localtime)
        print("Finished")    
        return 0
    else:
        # Applicaton of AIM to Isolates
        print("Step 2: Misassembly Breakpoint Identification")   
        localtime = time.asctime(time.localtime(time.time()))  
        print("Time:"+localtime)
        breakpoint_result =  single_breakpoint_detect(options)                     
        breakpoint_result.to_csv(options.output+"/mis_breakpoint.txt",sep="\t")
        localtime = time.asctime(time.localtime(time.time()))  
        print("Time:"+localtime)
        if options.no_correct:
            print("Finished")            
            return 0
        print("Step3: Correcting misassembled contigs")
        breakpoint_result = breakpoint_result.loc[breakpoint_result['misassembly_breakpoint'] > options.split_length,] 
        breakpoint_result = breakpoint_result.loc[(breakpoint_result['contig_length']-breakpoint_result['misassembly_breakpoint']) > options.split_length,]
        breakpoint_result=breakpoint_redult.iloc[np.argsort(-breakpoint_result['outlier_score']),]
        breakpoint_result = breakpoint_result.drop_duplicates(['contig'])                                   
        correct_result=corrected_contig(options,breakpoint_result)
        print("A total of " + str(correct_result.shape[0]) + " misassembled contigs are corrected")
        localtime = time.asctime(time.localtime(time.time()))  
        print("Time:"+localtime)
        print("Finished")   
        return 0         
        
def bamindex(options):                                          
    command_index=options.samtools+' index '+ options.bamfile
    os.system(command_index)     
    return options               
                                             
def main():
    (options, args)=get_opts() 
    # removing contigs with length smaller than required minimal length
    warnings.filterwarnings("ignore")   
    options.output = os.path.abspath(options.output)
    if not os.path.exists(options.output):
        os.system("mkdir -p "+ options.output)         
    localtime = time.asctime(time.localtime(time.time()))
    print("Time: "+ localtime)
    if not os.path.exists(options.assemblies):
        logging.error("Can not find contigs: " + options.assemblies)
        exit(-1)  
    else:
        if os.path.exists(options.assemblies):
            if os.path.exists(options.output+"/temp/contig/filtered_contigs.fa"):
                options.assemblies = options.output+"/temp/contig/filtered_contigs.fa"
            else:
                options=filter_contig(options)                                                                 
    if options.bamfile is None:
        if os.path.exists(options.output+"/temp/sam/contigs.filter.sort.bam"):
            options.bamfile = options.output+"/temp/sam/contigs.filter.sort.bam"
        else: 
            print("Mapping paired-end reads to assemblies")                         
            options=mapping(options) 
    elif os.path.exists(options.bamfile): 
        options=bamindex(options)                             
    if not os.path.exists(options.bamfile):
        logging.error("Can not find bamfile: "+ options.bamfile)
        exit(-1)                                                  
    pipeline(options)      
    
if __name__=='__main__':
    main()     
