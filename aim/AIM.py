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

base_path="/".join(sys.argv[0].split("/")[:-1])

def get_opts():
    parser=OptionParser()
    parser.add_option("-1","--r1",dest="read1",help="paired-end #1 fasta/q files") 
    parser.add_option("-2","--r2",dest="read2",help="paired-end #2 fasta/q files") 
    parser.add_option("-p","--r",dest="read",help="smart pairing (ignoring #2 fasta/q)") 
    parser.add_option("-c","--contig",dest="assemblies",help="input assembly fasta file")
    parser.add_option("--bam",dest="bamfile",help="index bam file for alignment") 
    parser.add_option("-a","--assembler",dest="assembler",default="MEGAHIT",help="The assembler used for assembly fasta file[MEGAHIT/IDBA_UD]") 
    parser.add_option("-o","--output",dest="output",help='output directory for MetaREA results')
    parser.add_option("-t","--threads",dest="threads",type=int, default=8,help='Maximum number of threads [default: 8]')
    parser.add_option("-l","--mlen",dest="min_length",type=int, default=5000,help='Minimum contig length [default: 5000bp]') 
    parser.add_option("-s","--slen",dest="split_length",type=int, default=1000,help='Minimum length to be broken [default: 1000bp]') 
    parser.add_option("--samtools", dest="samtools",default='samtools',help="path to samtools")
    parser.add_option("--bwa", dest="bwa",default='bwa',help="path to bwa") 
    parser.add_option("--bedtools", dest="bedtools",default='bedtools',help="path to bedtools")
    parser.add_option("--jellyfish",dest="jellyfish",default="jellyfish",help="path to jellyfish")  
    parser.add_option("--no-breakpoints",dest="no_break",action="store_true",help='Do not locate possible breakpoints')    
    parser.add_option("--no-correct",dest="no_correct",action="store_true",help='Do not break misassembled contigs at breakpoints')
    parser.add_option("-m","--mode", dest="mode",default='meta',help="Applied to single/metagenomic assemblies [meta/single]")
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
    
def fragcoverage_cal(options):
    # calculate fragment coverage for assemblies
    if not os.path.exists(options.output+"/temp/coverage/fragment_coverage.txt"):        
        coverage_command=' '.join(['python3',os.path.join(base_path,"frag_coverage.py"),
        '--bam',options.bamfile,'--output',options.output,'--mlen', str(options.min_length)])
        os.system(coverage_command) 
    
def read_feature(options):
    # extract features containing in mapped reads
    if not os.path.exists(options.output+"/temp/read_feature/read_feature.txt"):
        print("calculate paired read features...")    
        read_command=' '.join(['python3',os.path.join(base_path,"read_classify.py"),'--bam',options.bamfile,'--output',options.output,'--mlen',str(options.min_length)]) 
        os.system(read_command)                 
    
def coverage_cal(options):
    # calculate coverage for assemblies
    if not os.path.exists(options.output+"/temp/coverage/contigs.depth"):
        command=options.bedtools + ' genomecov -ibam ' + options.bamfile + " -d > " + options.output+"/temp/coverage/contigs.depth"
        os.system(command)
    if not os.path.exists(options.output+"/temp/coverage/coverage.txt"):        
        coverage_command=' '.join(['python3',os.path.join(base_path,"coverage.py"),
        '--bam',options.bamfile,'--output',options.output,'--bedtools',options.bedtools, "--mlen",str(options.min_length)])
        os.system(coverage_command) 
           
def pileup(options):
    print("calculate nucleotide variants...")
    # extract pileup features
    bam=options.bamfile
    pileup_file=options.output+"/temp/pileup/contigs_pipelup.out"
    if not os.path.exists(pileup_file):
        os.system("mkdir -p " + options.output+"/temp/pileup")
        pileup_command=' '.join([options.samtools,'mpileup -C 50 -A -f',options.assemblies,bam," | awk","'","$3 !=","\"N\"","'",">",pileup_file])                     
        os.system(pileup_command)  
    if not os.path.exists(options.output+"/temp/pileup/pileup_feature.txt"):
        pileup_command=' '.join(['python3',os.path.join(base_path,"pileup_parse.py"),"--pileup",pileup_file,"--output",options.output, "--mlen",str(options.min_length)])
        os.system(pileup_command)  
    
def read_breakpoint(options):
    # count the number of read breakpoints per base/window  
    print("calculate read breakpoints...")    
    if not os.path.exists(options.output+"/temp/read_breakpoint/read_breakpoint_per_window.txt"):
        read_breakpoint_command=' '.join(['python3',os.path.join(base_path,"read_breakpoint.py"),"--bam",options.bamfile,"--output",options.output,"--mlen",str(options.min_length)])          
        os.system(read_breakpoint_command)  
    
def KAD_cal(options):
    # split contigs and extract mapped reads
    if not os.path.exists(options.output+"/temp/KAD/KAD_window_data.txt"):
        split_command=' '.join(['sh',os.path.join(base_path,"split_sam.sh"),options.assemblies,options.output,options.samtools])            
        os.system(split_command)           
        KAD_command=' '.join(['python3',os.path.join(base_path,"KAD_cal.py"),'--bam',options.bamfile,'--contig',options.assemblies,'--output',options.output,'--thread',str(options.threads),'--jellyfish',options.jellyfish]) 
        os.system(KAD_command)

def read_feature_ratio(data):
    data['proper_read_ratio']=data['proper_read_count']/data['read_count'] 
    data['inversion_read_ratio']=data['inversion_read_count']/data['proper_read_count']   
    data['clipped_read_ratio']=data['clipped_read_count']/data['proper_read_count']
    data['supplementary_read_ratio']=data['supplementary_read_count']/data['proper_read_count']
    data['discordant_loc_ratio']=data['discordant_loc_count']/data['proper_read_count']
    data['discordant_size_ratio']=data['discordant_size_count']/data['proper_read_count']
    data = data.fillna(0)
    data=data.loc[:,['contig','start_pos','proper_read_ratio','inversion_read_ratio','clipped_read_ratio',
    'supplementary_read_ratio','discordant_loc_ratio','discordant_size_ratio','length']]
    return data
    
def cov_thread_cal(data):
    up = np.percentile(np.array(data),95)
    low = np.percentile(np.array(data),5)
    return up,low
        

def contig_fea(data):
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
    features=['coverage_width', 'deviation_width', 'normalized_coverage', 
    'normalized_deviation', 'window_cov_dev', 'fragment_width', 'fragment_deviation_width',
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
                 
                
def data_generation(options):   
    KAD_data=pd.read_csv(options.output+"/temp/KAD/KAD_window_data.txt",sep="\t",index_col=0)
    coverage=pd.read_csv(options.output+"/temp/coverage/coverage.txt",sep="\t",index_col=0)
    fragment_coverage=pd.read_csv(options.output+"/temp/coverage/fragment_coverage.txt",sep="\t",index_col=0)
    coverage_data=pd.merge(fragment_coverage,coverage,on=['contig','start_pos'])
    read_discordant_loc=pd.read_csv(options.output+"/temp/read_feature/discordant_loc_feature.txt",sep="\t",index_col=0)   
    read_discordant_size=pd.read_csv(options.output+"/temp/read_feature/discordant_size_feature.txt",sep="\t",index_col=0)    
    read_feature=pd.read_csv(options.output+"/temp/read_feature/read_feature.txt",sep="\t",index_col=0) 
    read_data=pd.merge(read_discordant_loc,read_discordant_size,on=['contig','start_pos'],how='right') 
    read_data=pd.merge(read_feature,read_data,how='left',on=['contig','start_pos'])  
    read_data=read_data.fillna(0) 
    read_data=read_feature_ratio(read_data)    
    pileup_data=pd.read_csv(options.output+"/temp/pileup/pileup_feature.txt",sep="\t",index_col=0)  
    window_read_breakpoint=pd.read_csv(options.output+"/temp/read_breakpoint/read_breakpoint_per_window.txt",sep="\t",index_col=0)
    window_data=pd.merge(pileup_data,KAD_data,on=['contig','start_pos'])
    window_data=pd.merge(window_data,coverage_data,on=['contig','start_pos'])  
    window_data['coverage_diff']=window_data['normalized_coverage']-window_data['normalized_fragment_coverage'] 
    window_data=pd.merge(window_data,window_read_breakpoint,on=['contig','start_pos'],how='left')
    window_data=window_data.fillna(0)
    window_data=pd.merge(window_data,read_data)
    window_data=window_data.loc[window_data['mean_coverage']>5,]
    os.system("mkdir -p "+options.output+"/feature_matrix")
    window_data.to_csv(options.output+"/feature_matrix/window_fea_matrix.txt",sep="\t")
    contig_data=contig_fea(window_data)
    contig_data.to_csv(options.output+"/feature_matrix/contig_fea_matrix.txt",sep="\t")    
    return contig_data

def predict(options,data):
    # Generate MetaREA contig scores
    min_length=options.min_length
    test_data = data.loc[data['length']>min_length,data.columns!='length']
    if options.assembler=='MEGAHIT':
        model_path=os.path.join(base_path,"model/MEGAHIT")
    elif options.assembler=='IDBA_UD':
        model_path=os.path.join(base_path,"model/IDBA_UD") 
    else:
        model_path=os.path.join(base_path,"model/merge")        
    score=pd.DataFrame(np.zeros([test_data.shape[0],10])) 
    score.index=test_data.index   
    #model_path='/home1/Laisenying/Data-analysis/projects/assembly_quality/MetaREA/model/MEGAHIT'            
    for i in range(10):                                                                                  
        rf=joblib.load(model_path+"/RF"+str(i)+'.pkl')
        pro=pd.DataFrame(rf.predict_proba(test_data))  
        pro.index=test_data.index
        score.loc[pro.index,i] = pro[1]
    score['MetaREA_contig_score']=score.mean(axis=1)
    score=score.loc[:,['MetaREA_contig_score']]
    score['length']=data.loc[score.index,'length']  
    score.to_csv(options.output+"/MetaREA_contig_score.txt",sep="\t")    
    return score  

def findcut(score):
    # find score cutoff to identify candidate misassembled contigs
    score_cut=np.percentile(score['MetaREA_contig_score'],98)
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
        outlier_fraction=0.001
    else:
        outlier_fraction=0.0001                    
    # fit the model
    clf=IsolationForest(contamination=outlier_fraction)
    clf.fit(Xdata)
    score_pred=clf.decision_function(Xdata)
    Xdata['outlier_score']=-score_pred
    if options.mode=='meta':
        threshold=stats.scoreatpercentile(-score_pred,100*(1-0.01))
    else:
        threshold=stats.scoreatpercentile(-score_pred,100*(1-0.001))                
    Xdata['outlier_thred']=threshold
    score_pred_data=Xdata.loc[:,['outlier_score','outlier_thred','read_breakpoint_ratio']]
    score_pred_data['start_pos']=[int(x.split("_")[-1]) for x in score_pred_data.index]
    score_pred_data['contig']=['_'.join(x.split("_")[:-1]) for x in score_pred_data.index]
    return score_pred_data                      
        
def meta_breakpoint_detect(options,filter_score):
    # identify misassembly breakpoints on metagenomic misassembled contigs
    if not os.path.exists(options.output+"/feature_matrix/window_fea_matrix.txt"):
        print("Missing window feature matrix file")
        exit(-1)
    if not os.path.exists(options.output+"/temp/read_breakpoint/read_breakpoint_per_base.txt"):
        print("Missing read breakpoint file")
        exit(-1)            
    window_data=pd.read_csv(options.output+"/feature_matrix/window_fea_matrix.txt",sep="\t",index_col=0)
    window_data.index=window_data['contig']
    window_data=window_data.loc[set(filter_score.index).intersection(window_data.index),]      
    window_data=window_data.fillna(0)
    window_data=window_data.replace(np.inf,0)
    window_data.index=window_data['contig']+"_"+[str(int(x)) for x in window_data['start_pos']]       
    score_pred_data = Isolation_forest(options,window_data)     
    score_pred_data.to_csv(options.output+"/outlier_score.txt")            
    contig_to_pos=contig_pos(score_pred_data) 
    outlier_data=pd.DataFrame(columns=['contig','start_pos','outlier_score','outlier_thred','read_breakpoint_ratio'])
    for contig in np.unique(score_pred_data['contig']):
        contig_outlier_score=score_pred_data.iloc[contig_to_pos[contig],]
        contig_outlier_score=contig_outlier_score.iloc[np.argsort(-contig_outlier_score['outlier_score']),]
        contig_outlier_score=contig_outlier_score.iloc[:2,]
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
    result['MetaREA_contig_score']=list(filter_score.loc[result['contig'],'MetaREA_contig_score'])
    breakpoint_result=result.loc[:,['contig', 'position', 'read_breakpoint_count', 'read_count',
        'read_breakpoint_ratio', 'outlier_score',"outlier_thred",'MetaREA_contig_score',"length"]]
    breakpoint_result.columns=['contig', 'misassembly_breakpoint', 'read_breakpoint_count', 'read_count',
       'read_breakpoint_ratio','outlier_score',"outlier_thred",'MetaREA_contig_score',"contig_length"]       
    breakpoint_result.to_csv(options.output+"/mis_breakpoint.txt",sep="\t")   
    return breakpoint_result     
    
    

def single_breakpoint_detect(options):
    # identify misassembly breakpoints on single genome assemblies
    if not os.path.exists(options.output+"/feature_matrix/window_fea_matrix.txt"):
        print("Missing window feature matrix file")
        exit(-1)
    if not os.path.exists(options.output+"/temp/read_breakpoint/read_breakpoint_per_base.txt"):
        print("Missing read breakpoint file")
        exit(-1)             
    window_data=pd.read_csv(options.output+"/feature_matrix/window_fea_matrix.txt",sep="\t",index_col=0)
    window_data=window_data.fillna(0)
    window_data=window_data.replace(np.inf,0)
    window_data.index=window_data['contig']+"_"+[str(int(x)) for x in window_data['start_pos']]                    
    score_pred_data = Isolation_forest(options,window_data)    
    score_pred_data.to_csv(options.ouput+"/outlier_score.txt")                     
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
    corrected_contig_file=options.output+"/corrected_contigs.fa"
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
                                                                                                                                                                                                                                                      
def pipeline(options):
    ##### Feature extraction #########
    print("############################################################################")
    print("Step1: Feature extraction")
    print("############################################################################")
    print("calculate coverage/fragment coverage for assemblies...")        
    fragcoverage_cal(options) 
    pool=multiprocessing.Pool(processes=options.threads)
    pool = [multiprocessing.Process(target=coverage_cal,args=(options,)),
        multiprocessing.Process(target=pileup,args=(options,)),
        multiprocessing.Process(target=read_breakpoint,args=(options,)),
        multiprocessing.Process(target=read_feature,args=(options,))]
    for t in pool:
        t.start()
    for t in pool:
        t.join()   
    print("Calculate KAD values...") 
    KAD_cal(options)  
    ######### Matrix generation ############################################
    print("Generate matrixs...") 
    contig_data=data_generation(options)    
    ######### Identify misassembled contigs ################################
    if options.mode=='meta':
        print("############################################################################")
        print("Step2: Identify misassembled contigs")
        print("############################################################################")
        score=predict(options,contig_data)            
    ######### Misassembly breakpoint identification ########################
        if not options.no_break:
            print("############################################################################")
            print("Step3: Misassembly Breakpoint Identification")
            print("############################################################################")       
            score_cut=findcut(score)
            filter_score=score.loc[score['MetaREA_contig_score']>score_cut,]
            breakpoint_result = meta_breakpoint_detect(options,filter_score)
    ######### Break misassembled contigs at misassembled breakpoints ########
            if not options.no_correct:
                print("############################################################################")
                print("Step4: Correcting misassembled contigs...")
                print("############################################################################")                                     
                correct_result=corrected_contig(options,breakpoint_result)
                print("A total of " + str(correct_result.shape[0]) + " misassembled contigs are corrected")
                print("################### Finishing MetaREA ######################")    
    elif options.mode=="single":
        print("############################################################################")
        print("Step 2: Misassembly Breakpoint Identification")
        print("############################################################################")       
        breakpoint_result =  single_breakpoint_detect(options)     
        breakpoint_result = breakpoint_result.loc[breakpoint_result['misassembly_breakpoint']>1000,] 
        breakpoint_result = breakpoint_result.loc[(breakpoint_result['contig_length']-breakpoint_result['misassembly_breakpoint'])>1000,]
        breakpoint_result.to_csv(options.output+"/mis_breakpoint.txt",sep="\t")
        if not options.no_correct:
            print("############################################################################")
            print("Step3: Correcting misassembled contigs...")
            print("############################################################################")  
            breakpoint_result=breakpoint_redult.iloc[np.argsort(-breakpoint_result['outlier_score']),]
            breakpoint_result = breakpoint_result.drop_duplicates(['contig'])                                   
            correct_result=corrected_contig(options,breakpoint_result)
        
def bamfilter(options):
    output_bam = output_bam=options.output+"/temp/sam/contigs.filter.sort.bam"
    command_bam = options.samtools + ' view -h -q 10 -m 50 -F 4 -b ' + options.bamfile + " | " + options.samtools + " sort " + ' > ' + output_bam        
    os.system(command_bam)                                       
    command_index=options.samtools+' index '+output_bam
    os.system(command_index)  
    options.bamfile = output_bam                  
                                             
def main():
    (options, args)=get_opts() 
    # removing contigs with length smaller than required minimal length
    warnings.filterwarnings("ignore")   
    if not os.path.exists(options.output):
        options.output = os.path.abspath(options.output)
        os.system("mkdir -p "+options.output)         

    if not os.path.exists(options.assemblies):
        print("Can not find "+options.assemblies)
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
            print("############################################################################")
            print("Mapping paired-end reads to assemblies")
            print("############################################################################")                   
            options=mapping(options) 
    elif os.path.exists(options.bamfile): 
        if not os.path.exists(options.output+"/temp/sam/contigs.filter.sort.bam"):
            options=bamfilter(options)
        else:
            options.bamfile = options.output+"/temp/sam/contigs.filter.sort.bam"                          
    if not os.path.exists(options.bamfile):
        print("Can not find bamfile: "+options.bamfile)
        exit(-1)
    if not os.path.exists(options.assemblies):
        print("Can not find assemblies: "+options.assemblies)
        exit(-1)                                                    
    pipeline(options)      
    
if __name__=='__main__':
    main()  
