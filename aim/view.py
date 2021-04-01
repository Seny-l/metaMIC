import pandas as pd
import numpy as np
import argparse,operator,os,random,sys,time
import warnings
from bokeh.plotting import figure, output_file, show,save
from bokeh.layouts import column
from bokeh.layouts import gridplot

def parseargs():
    parser=argparse.ArgumentParser(description="Visualization of contig quality")
    parser.add_argument("--contig",help="contig name for visualization")
    parser.add_argument('--output',help='output directory for AIM result')  
    parser.add_argument("--width",type=int, default=1000,help='with for visualization')
    args=parser.parse_args()
    return args
        
    
def visualize(args):  
    data = pd.read_csv(args.output+"/feature_matrix/window_fea_matrix.txt",sep="\t",index_col=0)
    data.index=data['contig']
    contig_data = data.loc[args.contig,]
    contig_data = contig_data.loc[:,['start_pos','inversion_read_ratio','clipped_read_ratio','supplementary_read_ratio',
    'discordant_loc_ratio','discordant_size_ratio','disagree_portion','normalized_fragment_coverage','mean_KAD','normalized_coverage']]    
    outlier_data = pd.read_csv(args.output+"/outlier_score.txt",sep="\t",index_col=0)
    outlier_data['start_pos'] = [int(x.split("_")[-1]) for x in outlier_data.index]
    outlier_data['contig'] = ['_'.join(x.split("_")[:-1]) for x in outlier_data.index]
    outlier_data.index=outlier_data['contig']
    if args.contig in outlier_data.index:
        contig_outlier = outlier_data.loc[args.contig,] 
        contig_outlier = contig_outlier.loc[:,['start_pos','outlier_score']] 
        contig_data = pd.merge(contig_outlier,contig_data,on=['contig','start_pos'])
    breakdata = pd.read_csv(args.output+"/temp/read_breakpoint/read_breakpoint_per_base.txt",sep="\t",index_col=0)
    breakdata.index=breakdata['contig']
    contig_break = breakdata.loc[args.contig,['contig','position','read_breakpoint_count']]
    contig_break.index=range(contig_break.shape[0])
    null_data = pd.DataFrame({"position":range(1,np.max(contig_data['start_pos']),1)})
    contig_break=pd.merge(pd.DataFrame({"position":range(1,np.max(contig_data['start_pos']),1)}),contig_break,on="position",how='left')
    contig_break['read_breakpoint_count'] = contig_break['read_breakpoint_count'].fillna(0)
    os.system("mkdir -p " + args.output+"/html")
    output_file(args.output+"/html/"+str(args.contig)+'.html')          
    p1 = figure(plot_width=args.width, plot_height=150,title="fragment coverage")
    p1.line(list(contig_data['start_pos']), list(contig_data['normalized_fragment_coverage']), color='navy', alpha=1,line_width=2)
    p2 = figure(plot_width=args.width, plot_height=120,title="inversion reads")
    p2.vbar(x=list(contig_data['start_pos']),width=0.9,bottom=0,top=list(contig_data['inversion_read_ratio']), color='darkgreen', alpha=0.8)        
    p3 = figure(plot_width=args.width, plot_height=120,title="clipped reads")
    p3.vbar(x=list(contig_data['start_pos']),width=0.9,bottom=0,top=list(contig_data['clipped_read_ratio']), color='darkgreen', alpha=0.8)    
    p4 = figure(plot_width=args.width, plot_height=120,title="discordant reads(paired reads with abnormal distance)")
    p4.vbar(x=list(contig_data['start_pos']),width=0.9,bottom=0,top=list(contig_data['discordant_size_ratio']), color='darkgreen', alpha=0.8)    
    p5 = figure(plot_width=args.width, plot_height=120,title="discordant reads(paired reads mapped to different contigs)")
    p5.vbar(x=list(contig_data['start_pos']),width=0.9,bottom=0,top=list(contig_data['discordant_loc_ratio']), color='darkgreen', alpha=0.8)    
    p6 = figure(plot_width=args.width, plot_height=120,title="supplementary reads")
    p6.vbar(x=list(contig_data['start_pos']),width=0.9,bottom=0,top=list(contig_data['supplementary_read_ratio']), color='darkgreen', alpha=0.8)    
    p7 = figure(plot_width=args.width, plot_height=120,title="disagreement bases")
    p7.vbar(x=list(contig_data['start_pos']),width=0.9,bottom=0,top=list(contig_data['disagree_portion']), color='darkorange', alpha=0.8)    
    p8 = figure(plot_width=args.width, plot_height=120,title="KAD")
    p8.vbar(x=list(contig_data['start_pos']),width=0.9,bottom=0,top=list(contig_data['mean_KAD']), color='darkred', alpha=0.8)    
    p9 = figure(plot_width=args.width, plot_height=120,title="read breakpoints")
    p9.vbar(x=list(contig_break['position']),width=0.9,bottom=0,top=list(contig_break['read_breakpoint_count']), color='darkred', alpha=0.8)
    if "outlier_score" in contig_data.columns:
        p10 = figure(plot_width=args.width, plot_height=120,title="outlier score")
        p10.line(list(contig_data['start_pos']), list(contig_data['outlier_score']), color='darkblue', alpha=1,line_width=2)    
        grid=gridplot([[p1], [p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10]])
    else:
        grid=gridplot([[p1], [p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9]])                
    save(grid)   
    
def main():
    args=parseargs()
    warnings.filterwarnings("ignore")
    visualize(args)
                            
if __name__=='__main__':
    main()  
