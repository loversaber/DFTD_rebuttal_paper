from sys import argv
import random
import pandas as pd
import os
from multiprocessing import Pool
import numpy as np
import glob
from scipy import stats

def original_samples():
    d_pop_true_samples={}
    actual_sample_dir="/nfs/users/nfs_q/ql4/lustre_ql4/work/pro_nc_paper/pro360samples/final_vcf/cal_MAF_proportion_3focal/permutation_test/origin_st_ed_samples"
    actual_ALL_sample_dir="/nfs/users/nfs_q/ql4/lustre_ql4/work/pro_nc_paper/pro360samples/final_vcf/fast_phase/unphased/time_series/wwp_Ne"
    for pop in d_pop_years_st_ed:
        d_pop_true_samples[pop]={}
        st_file=f"{actual_sample_dir}/{pop}_start.samples"
        ed_file=f"{actual_sample_dir}/{pop}_end.samples"
        df_st=pd.read_csv(st_file,sep="\t",header=None,names=["name1","name2"])
        df_ed=pd.read_csv(ed_file,sep="\t",header=None,names=["name1","name2"])
        ALL_samples=f"{actual_ALL_sample_dir}/{pop}_ALL.samples"
        df_ALL=pd.read_csv(ALL_samples,sep="\t",header=0)
        st_samples=df_st["name1"].to_list()
        ed_samples=df_ed["name1"].to_list()
        all_samples=df_ALL["sample_2"].to_list()#st_samples+ed_samples
        d_pop_true_samples[pop]["st"]=st_samples
        d_pop_true_samples[pop]["ed"]=ed_samples
        d_pop_true_samples[pop]["ALL"]=all_samples
    return d_pop_true_samples

def write_samples(samples_file,samples_list):
    df=pd.DataFrame(samples_list)
    #print(df)
    df[1]=df[0]
    df.to_csv(samples_file,index=False,header=False,sep="\t")

def split_3pop_vcf_randomly(run,samples_list0):
    samples_list=samples_list0.copy()
    d_pop_sted_samples={}#store the random sampling
    for pop in d_pop_years_st_ed:
        d_pop_sted_samples[pop]={}
        pop_sampls_list=d_pop_true_samples[pop]["ALL"]
        random.shuffle(pop_sampls_list)

        st=d_pop_years_st_ed[pop]["st"]#num
        #st_samples=d_pop_sted_samples[pop]["st"]=d_pop_true_samples[pop]["st"]
        st_samples=d_pop_sted_samples[pop]["st"]=pop_sampls_list[:st]
        d_pop_sted_samples[pop]["st_index"]=[n for n,i in enumerate(samples_list0) if i in d_pop_sted_samples[pop]["st"]]

        ed=d_pop_years_st_ed[pop]["ed"]
        #ed_samples=d_pop_sted_samples[pop]["ed"]=d_pop_true_samples[pop]["ed"]
        ed_samples=d_pop_sted_samples[pop]["ed"]=pop_sampls_list[st:st+ed]
        d_pop_sted_samples[pop]["ed_index"]=[n for n,i in enumerate(samples_list0) if i in d_pop_sted_samples[pop]["ed"]]

        all_samples=d_pop_sted_samples[pop]["ALL"]=d_pop_true_samples[pop]["ALL"]#d_pop_sted_samples[pop]["st"]+d_pop_sted_samples[pop]["ed"]
        d_pop_sted_samples[pop]["ALL_index"]=[n for n,i in enumerate(samples_list0) if i in d_pop_sted_samples[pop]["ALL"]]
        ###write the samples to the samples file which will be used in prune
        print(pop,len(st_samples),len(ed_samples),st,ed)
        all_samples_file=pop+"_all.samples_"+run
        st_samples_file=pop+"_start.samples_"+run
        ed_samples_file=pop+"_end.samples_"+run
        d_pop_sted_samples[pop]["st_file"]=st_samples_file
        d_pop_sted_samples[pop]["ed_file"]=ed_samples_file
        pop_samples_files=[st_samples_file,ed_samples_file,all_samples_file]
        pop_samples_lists=[st_samples,ed_samples,all_samples]
        for pop_samples_file,pop_samples_list in zip(pop_samples_files,pop_samples_lists):
            write_samples(pop_samples_file,pop_samples_list)
    return d_pop_sted_samples#store the two time point samples and samples' index

def cal_GT_rate_MAF(tp_GT,threshold_rate):
    gt=[i for i in tp_GT if i!="./."]
    gt_rate=len(gt)/len(tp_GT)#will increase more variant for gt13 in Nara
    #gt_rate=round(len(gt)/len(tp_GT),6)
    if gt_rate>=threshold_rate:
        k=[]
        for i in gt:
            c=i.split("/")
            k.extend(c)
        n_B=len([j for j in k if j=="1"])
        n_Total=len(gt)*2
        maf=n_B/n_Total#round(n_B/n_Total,6)
        syb=1
    else:
        syb=0
        maf=0
        n_Total=0
        n_B=0
    return syb,n_B,n_Total,maf,gt_rate

def write_pop_samples_vcf_file(d_vcf,d_MAF,d_gt,d_allgt13_vcf,d_pop_sted_samples,line):#store the GT1/3 and maf of snp in two time point
    line_list=line.strip().split("\t")
    vcf_col=line_list[:9]
    maf_col1=line_list[:5]
    samples_col=line_list[9:]#0/0:0.0,-3.05,-16.3:0:31:3:0
    for pop in d_pop_sted_samples:
        pop_vcf=d_vcf[pop]
        pop_allgt13_vcf=d_allgt13_vcf[pop]
        pop_maf_handle=d_MAF[pop]
        st_index=d_pop_sted_samples[pop]["st_index"]
        ed_index=d_pop_sted_samples[pop]["ed_index"]
        pop_all_samples_index_for_tp=st_index+ed_index#!!!!!!!this is a mistake if calculating RSB
        pop_all_samples_index_for_allGT=d_pop_sted_samples[pop]["ALL_index"]

        st_GT=[samples_col[i].split(":")[0] for i in st_index ]
        ed_GT=[samples_col[i].split(":")[0] for i in ed_index ]
        st_rate_threshold=d_rate_threshold[pop]["st"]
        ed_rate_threshold=d_rate_threshold[pop]["ed"]
        st_mark,n_B1,n_Total1,maf1,gt_rate1=cal_GT_rate_MAF(st_GT,st_rate_threshold)
        ed_mark,n_B2,n_Total2,maf2,gt_rate2=cal_GT_rate_MAF(ed_GT,ed_rate_threshold)
        #gt_rate_handle=d_gt[pop]
        #gt_line=[str(gt_rate1),str(gt_rate2)]
        #gt_rate_handle.write("\t".join(maf_col1+gt_line)+"\n")
        ####RSB based on maf >=0.05
        all_GT=[samples_col[i].split(":")[0] for i in pop_all_samples_index_for_allGT]
        ccc=len(all_GT)
        #print(f"{pop}_all_GT13:{ccc}")
        all_mark,n_B0,n_Total0,maf0,gt_rate0=cal_GT_rate_MAF(all_GT,1/3)
        #write gt rate to file
        gt_rate_handle=d_gt[pop]
        gt_line=[str(gt_rate1),str(gt_rate2),str(gt_rate0)]
        gt_rate_handle.write("\t".join(maf_col1+gt_line)+"\n")
        if all_mark==1 and maf1>=0.05 and maf2>=0.05:
            #print(f"{pop}_ALLGT13_yes_{gt_rate0}_{maf1}_{maf2}")
            all_samples_allgt13_info=[samples_col[i] for i in pop_all_samples_index_for_allGT]
            pop_allgt13_vcf.write("\t".join(vcf_col+all_samples_allgt13_info)+"\n")
        #else:
            #print(f"{pop}_ALLGT13_no_{gt_rate0}_{maf1}_{maf2}")
        ####
        if st_mark==1 and ed_mark==1:
            all_samples_gt_info=[samples_col[i] for i in pop_all_samples_index_for_tp]
            pop_vcf.write("\t".join(vcf_col+all_samples_gt_info)+"\n")
            maf_line=[str(n_B1),str(n_Total1),str(maf1),str(n_B2),str(n_Total2),str(maf2)]
            pop_maf_handle.write("\t".join(maf_col1+maf_line)+"\n")


def write_header_to_pop_vcf(d_vcf_handle,line):
    for pop in d_vcf_handle:#.keys():
        pop_vcf_file=d_vcf_handle[pop]
        pop_vcf_file.write(line)


def write_header_to_pop_allgt13vcf(d_vcf_handle):
    for pop in d_vcf_handle:
        d_vcf_handle[pop].write('##fileformat=VCFv4.0\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

def pro_vcf(run,d_vcf,d_vcf_file,d_MAF,d_MAF_file,d_gt,d_gt_file,d_allgt13_vcf,d_allgt13_vcf_file):
    with open(vcf_file,"r") as f:
        for line in f:
            if line.startswith("#CHROM"):
                line_list=line.strip().split("\t")
                vcf_col=line_list[:9]
                maf_col1=line_list[:5]
                samples_list=line_list[9:]
                d_pop_sted_samples=split_3pop_vcf_randomly(run,samples_list)
                #d_pop_sted_samples=generate_origin_sted_samples(samples_list)
                for pop in d_vcf:
                    pop_vcf_file=d_vcf[pop]
                    pop_maf_handle=d_MAF[pop]
                    gt_rate_handle=d_gt[pop]
                    allgt13_vcf_handle=d_allgt13_vcf[pop]
                    st_index=d_pop_sted_samples[pop]["st_index"]
                    ed_index=d_pop_sted_samples[pop]["ed_index"]
                    pop_all_samples_list_for_tp=[samples_list[i] for i in st_index+ed_index]
                    pop_all_samples_list_for_allGT13=[samples_list[i] for i in d_pop_sted_samples[pop]["ALL_index"]]
                    pop_vcf_file.write("\t".join(vcf_col+pop_all_samples_list_for_tp)+"\n")
                    allgt13_vcf_handle.write("\t".join(vcf_col+pop_all_samples_list_for_allGT13)+"\n")
                    pop_maf_handle.write("\t".join(maf_col1+["nB_st","nTotal_st","MAF_st","nB_ed","nTotal_ed","MAF_ed"])+"\n")
                    gt_rate_handle.write("\t".join(maf_col1+["gtRate_st","gtRate_ed","gtRate_ALL"])+"\n")

            elif line.startswith("##"):
                write_header_to_pop_vcf(d_vcf,line)
            else:
                write_pop_samples_vcf_file(d_vcf,d_MAF,d_gt,d_allgt13_vcf,d_pop_sted_samples,line)
    return d_pop_sted_samples


def pro_prune(run,pop,d_vcf_file,d_MAF_file,d_pop_sted_samples):
    print("#Plink Prune:")
    pop_vcf_file=d_vcf_file[pop]
    pop_maf_file=d_MAF_file[pop]
    cmd="plink --double-id --indep-pairwise 50 20 0.99 --keep {:s}_all.samples_{:s} --out {:s}_pruned_{:s} --vcf {:s}".format(pop,run,pop,run,pop_vcf_file)
    os.system(cmd)
    df_pruned_in_snps=pd.read_csv(pop+"_pruned_"+run+".prune.in",header=None)
    df_pruned_in_snps=df_pruned_in_snps.copy()
    df_pruned_in_snps=df_pruned_in_snps.sample(n=d_pop_nSNPs[pop])#select same number of SNPs like the original run
    snp_needed=df_pruned_in_snps[0].to_list()
    st_samples=d_pop_sted_samples[pop]["st"]
    ed_samples=d_pop_sted_samples[pop]["ed"]
    df_snp_MAF=pd.read_csv(pop_maf_file,sep="\t",header=0)
    df_snp_MAF=df_snp_MAF.copy()
    df_snp_MAF=df_snp_MAF[df_snp_MAF["ID"].isin(snp_needed)]
    df_snp_MAF["Change"]=(df_snp_MAF["MAF_st"]-df_snp_MAF["MAF_ed"]).abs()
    df_snp_MAF=df_snp_MAF.sort_values(["Change"],ascending=False).reset_index(drop=True)
    n=2.5
    q25=df_snp_MAF.head(int(len(df_snp_MAF)*(n/100)))
    pop_q25=pop+"_rank25.txt_"+run
    q25.to_csv(pop_q25,sep="\t",header=True,index=False,float_format='%.6f')
    pruned_in_file="{:s}_pruned_{:s}.prune.in".format(pop,run)
    #num_lines1=sum(1 for line in open(pruned_in_file))
    num_lines1=len(df_pruned_in_snps)
    print(pruned_in_file,num_lines1)
    return q25,pop_q25,num_lines1,df_snp_MAF

def get_rank25(run,df,pop):#merge rank25 bed
    df=df.copy()
    df["#CHROM"]=df["#CHROM"].astype(str)
    df["start"]=df["POS"]-100000
    df["end"]=df["POS"]+100000
    #df["start"].values[[df["start"].values<0]]=0
    df[["start"]]=df[["start"]].clip(lower=0)
    df_bed=df[["#CHROM","start","end"]]
    df_bed=df_bed.sort_values(by=["#CHROM","start"],ascending=True)
    bed_file=pop+".bed_"+run
    bed_merged_file=pop+".merged.bed_"+run#merge some region with overlap
    df_bed.to_csv(bed_file,sep="\t",header=False,index=False)
    cmd_bed="bedtools merge -i %s > %s"%(bed_file,bed_merged_file)
    os.system(cmd_bed)
    return df,bed_merged_file

def cal_intersection(run,d_vcf_file,d_MAF_file,d_pop_sted_samples):
    print("#Merge rank25 in 3 pops:")
    l=[]
    d_num_lines2={}
    d_pop_maf={}
    d_pop_maf_change25={}
    for pop in ["Freycinet","Narawntapu","WWP"]:
        df25,rank25file,n_pruned,df_snp_MAF=pro_prune(run,pop,d_vcf_file,d_MAF_file,d_pop_sted_samples)
        d_pop_maf_change25[pop]=df25
        d_pop_maf[pop]=df_snp_MAF
        d_num_lines2[pop]=n_pruned
        df25_interval,pop_bed_merged_file=get_rank25(run,df25,pop)
        l.append(pop_bed_merged_file)
    cmd_bed1="bedtools intersect -a {:s} -b {:s} > frey_nara_intersection.bed_{:s}".format(l[0],l[1],run)
    cmd_bed2="bedtools intersect -a frey_nara_intersection.bed_{:s} -b {:s} > frey_nara_wwp_intersection.bed_{:s}".format(run,l[2],run)
    os.system(cmd_bed1)
    os.system(cmd_bed2)
    df_inter=pd.read_csv("frey_nara_wwp_intersection.bed_{:s}".format(run),sep="\t",header=None,names=["CHR","start","end"])
    df_inter=df_inter.copy()
    df_inter["length"]=df_inter["end"]-df_inter["start"]+1
    total_length=df_inter["length"].sum()
    n_inter=len(df_inter)
    return df_inter,n_inter,total_length,d_num_lines2,d_pop_maf,d_pop_maf_change25

def global_var(run):
    d_vcf={}
    d_vcf_file={}
    d_MAF={}
    d_MAF_file={}
    d_gt={}
    d_gt_file={}
    d_allgt13_vcf={}
    d_allgt13_vcf_file={}
    for pop1,pop2 in zip(["Freycinet","Narawntapu","WWP"],["Frey","Nara","WWP"]):
        pop2_vcf_file=pop2+".gt13.vcf_"+run#two time points 1/3
        d_vcf[pop1]=open(pop2_vcf_file,"w")
        d_vcf_file[pop1]=pop2_vcf_file
        pop2_maf_file=pop2+".gt13.maf_"+run#store MAF
        d_MAF[pop1]=open(pop2_maf_file,"w")
        d_MAF_file[pop1]=pop2_maf_file
        pop2_gt_file=pop2+".gtRate_"+run#genotyping rate
        d_gt[pop1]=open(pop2_gt_file,"w")
        d_gt_file[pop1]=pop2_gt_file
        pop2_allgt13_vcf_file=pop2+".gtAll13.vcf_"+run
        d_allgt13_vcf[pop1]=open(pop2_allgt13_vcf_file,"w")
        d_allgt13_vcf_file[pop1]=pop2_allgt13_vcf_file
    write_header_to_pop_allgt13vcf(d_allgt13_vcf)#write manual 2 header lines
    return d_vcf,d_vcf_file,d_MAF,d_MAF_file,d_gt,d_gt_file,d_allgt13_vcf,d_allgt13_vcf_file

########################Cal Rsb value 
def make_tp_vcf(allgt13_vcf_file,chr_i,out_prefix,keep_samples):
    vcftools="/software/CGP/modules/installs/vcftools/vcftools-0.1.16/bin/vcftools"
    cmd ="{:s} --vcf {:s} --recode --chr {:s} --out {:s} --keep {:s}".format(vcftools,allgt13_vcf_file,chr_i,out_prefix,keep_samples)
    os.system(cmd)
    opt_tp_vcf=out_prefix+".recode.vcf"
    return opt_tp_vcf

def cal_rsb(run,pop,d_allgt13_vcf_file,d_pop_sted_samples):
    st_file=d_pop_sted_samples[pop]["st_file"]
    ed_file=d_pop_sted_samples[pop]["ed_file"]
    allGT13_dir="/nfs/users/nfs_q/ql4/lustre_ql4/work/pro_nc_paper/pro360samples/final_vcf/fast_phase/unphased/new_calc_RSB/setMAF005_RSB/allGT13_MAF005_VCF/"
    allgt13_vcf_file=f"{allGT13_dir}{pop}.gt13_maf005.vcf"#d_allgt13_vcf_file[pop]#pop is "Freycinet","Narawntapu","WWP"
    vcftools="/software/CGP/modules/installs/vcftools/vcftools-0.1.16/bin/vcftools"
    rscript="/nfs/users/nfs_q/ql4/anaconda3/envs/R3.6/bin/Rscript"
    l_df_rsb=[]
    for chr_i in range(1,7):
        chr_i=str(chr_i)
        opt_prefix=pop+"_"+chr_i+"_"+run
        pop_chr_rsb=opt_prefix+".rsb"
        pre_prefix=opt_prefix+"_pre"
        post_prefix=opt_prefix+"_post"
        opt_pre_vcf =make_tp_vcf(allgt13_vcf_file,chr_i,pre_prefix,st_file)
        opt_post_vcf=make_tp_vcf(allgt13_vcf_file,chr_i,post_prefix,ed_file)
        cmd_rsb="{:s} cal_Rsb.R {:s} {:s} {:s}".format(rscript,opt_post_vcf,opt_pre_vcf,pop_chr_rsb)
        os.system(cmd_rsb)
        df_rsb=pd.read_csv(pop_chr_rsb,sep="\t",header=0)
        #   NA should be removed or this column be str
        df_rsb["RSB_POST_PRE"]=df_rsb["RSB_POST_PRE"].str.strip()
        df_rsb=df_rsb[df_rsb["RSB_POST_PRE"]!="NA"]
        #####IMPORTANT#22.05.12
        df_rsb=df_rsb[df_rsb["RSB_POST_PRE"]==df_rsb["RSB_POST_PRE"]]
        df_rsb["RSB_POST_PRE"]=df_rsb["RSB_POST_PRE"].astype(float)
        df_rsb["LOGPVALUE"]=df_rsb["LOGPVALUE"].astype(float)
        l_df_rsb.append(df_rsb)
    df_all_rsb=pd.concat(l_df_rsb)
    df_all_rsb=df_all_rsb.sort_values(["RSB_POST_PRE","LOGPVALUE"],ascending=[False,False]).reset_index(drop=True)
    n=2.5
    rsb25=df_all_rsb.head(int(len(df_all_rsb)*(n/100)))
    return rsb25,df_all_rsb

########Calculating the Composite pvalue
def get_quantile(df,col1):
    df=df.copy()
    value_list1=df[col1].to_list()
    df[col1+"_q"]=(stats.rankdata(value_list1, 'min'))/len(value_list1)
    return df

def split_genome(df,col1,col2):
    print("#Composite Pvalue:")
    df=df.copy()
    l_df=[]
    df=get_quantile(df,col1)
    df=get_quantile(df,col2)
    for chr_i in df["CHR"].unique():
        df_chr=df[df["CHR"]==chr_i]
        df_chr=df_chr.copy()
        df_chr=df_chr.sort_values(by=["POSITION"])
        pos_list=df_chr["POSITION"].to_list()
        st=df_chr["POSITION"].iloc[0]
        ed=df_chr["POSITION"].iloc[-1]
        print("MaxMin",st,ed)
        if st!=ed:#not only one pos
            index_st=st
            step=100000#100kb
            nrange=np.ceil((ed-st)/step)+1#should plus 1 to contain the end
            index_ed=index_st+nrange*step
            arange_array=np.arange(index_st,index_ed,step)
            #print(chr_i,"yes",arange_array)
        else:
            index_st=st
            step=100000#100kb
            index_ed=index_st+100001
            arange_array=np.arange(index_st,index_ed,step)
            print(chr_i,"no")
        print(chr_i)
        for index,i in enumerate(arange_array[:-1]):
            n_window=index+1
            start=i
            if index!=len(arange_array)-1:
                end=arange_array[index+1]
            else:
                end=arange_array[-1]+step
            df_window=df_chr[(df_chr["POSITION"]>=start)&(df_chr["POSITION"]<=end)]
            df_window=df_window.copy()
            if not df_window.empty:
                #print(df_window[["CHR","POSITION",col1,col1+"_q",col2,col2+"_q"]])
                s1=len(df_window[df_window[col1]==df_window[col1]])
                s2=len(df_window[df_window[col2]==df_window[col2]])
                df_window["Nwindow"]=n_window
                df_window["start"]=start
                df_window["end"]=end
                df_window["S1"]=s1
                df_window["S2"]=s2
                if s1!=0 and s2!=0:
                    max1=df_window[col1].max()
                    max1_q=df_window[df_window[col1]==max1][col1+"_q"].values[0]
                    max2=df_window[col2].max()
                    max2_q=df_window[df_window[col2]==max2][col2+"_q"].values[0]
                    #print(max1,max2,s1,s2)
                    if max1_q==1:
                        max1_q=0.99999
                    if max2_q==1:
                        max2_q=0.99999
                    Pi1=1-(max1_q)**s1
                    Pi2=1-(max2_q)**s2
                    df_window["pval1"]=Pi1
                    df_window["pval2"]=Pi2
                    #score=(-2)*(np.log(Pi1)+np.log(Pi2))
                    #pvalue=chi2.sf(score,4)
                    stats2,pval2=stats.combine_pvalues([Pi1,Pi2])
                    log_pval2=-np.log10(pval2)
                    pval2=f"{pval2:.6f}"
                    df_window["Cpvalue"]=pval2
                    #log_pval2=f"{log_pval2:.6f}"
                    df_window["LogCpvalue"]=log_pval2
                else:
                    df_window["Cpvalue"]=1
                    df_window["LogCpvalue"]=0
                    if s1==0:
                        print(f"#No {col1}")
                    else:
                        print(f"#No {col2}")
                l_df.append(df_window)
    df_end0=pd.concat(l_df)
    n_p=3
    #df_end0.head(int(len(df_end0)*(n_p/100)))
    df_composite_pvalue03=df_end0[df_end0["LogCpvalue"]>df_end0["LogCpvalue"].quantile(0.97)]
    return df_composite_pvalue03,df_end0

def get_inter(df_rank,df_inter):    
    n_pop=0
    for chr_i in df_rank["CHR"].unique():
        df_inter_tmp=df_inter[df_inter["CHR"]==chr_i]
        if not df_inter_tmp.empty:
            lows=df_inter_tmp["start"].to_numpy()
            ups=df_inter_tmp["end"].to_numpy()
            df_tmp=df_rank[df_rank["CHR"]==chr_i]
            pos_ok=[i for i in df_tmp["POSITION"].to_list() if in_range(i,lows,ups)]
            n_pop+=len(pos_ok)
    return n_pop

def get_inter2(df_rank,df_inter):
    l=[]
    for chr_i in df_inter["CHR"].unique():#not finished (should cal for each region)
        df_inter_tmp=df_inter[df_inter["CHR"]==chr_i]
        lows=df_inter_tmp["start"].to_numpy()
        ups=df_inter_tmp["end"].to_numpy()
        df_tmp=df_rank[df_rank["CHR"]==chr_i]
        pos_ok=[i for i in df_tmp["POSITION"].to_list() if in_range(i,lows,ups)]
        df_tmp=df_tmp[df_tmp["POSITION"].isin(pos_ok)]["Cpvalue"].mean()
    return n_pop

def main_pro(run):
    run=str(run)
    d_vcf,d_vcf_file,d_MAF,d_MAF_file,d_gt,d_gt_file,d_allgt13_vcf,d_allgt13_vcf_file=global_var(run)
    d_pop_sted_samples=pro_vcf(run,d_vcf,d_vcf_file,d_MAF,d_MAF_file,d_gt,d_gt_file,d_allgt13_vcf,d_allgt13_vcf_file)
    for pop in d_vcf:
        d_vcf[pop].close()
        d_MAF[pop].close()
        d_gt[pop].close()
        d_allgt13_vcf[pop].close()
    d_num_lines={}
    for pop in d_vcf_file:
        num_lines=sum(1 for line in open(d_vcf_file[pop]))-48
        d_num_lines[pop]=num_lines
        print(d_vcf_file[pop],num_lines)
    df_inter,n_inter,toatl_len,d_num_lines2,d_pop_maf,d_pop_maf_change25=cal_intersection(run,d_vcf_file,d_MAF_file,d_pop_sted_samples)#intersections for the 3 pop 
    #l_rsb25=[]
    ####
    print("###Rsb Calculation:")
    #n_inter_rsb25=0
    df_inter_maf25={}
    d_inter_rsb25={}
    d_inter_cp03={}
    for pop in d_allgt13_vcf_file:
        df_rsb25,df_all_rsb=cal_rsb(run,pop,d_allgt13_vcf_file,d_pop_sted_samples)
        df_all_maf=d_pop_maf[pop]#for composite analysis #only pruned SNPs
        df_all_maf=df_all_maf.rename(columns={"#CHROM":"CHR","POS":"POSITION"})
        df_maf_rsb=pd.merge(df_all_maf,df_all_rsb,on=["CHR","POSITION"],how="left")
        df_composite_pvalue03,df_composite_pvalue=split_genome(df_maf_rsb,"Change","RSB_POST_PRE")
        ###for MAF change;check how many snps in intersections
        df_maf25=d_pop_maf_change25[pop]
        df_maf25=df_maf25.rename(columns={"#CHROM":"CHR","POS":"POSITION"})
        df_inter_maf25[pop]=get_inter(df_maf25,df_inter)
        ###for Rsb interval
        d_inter_rsb25[pop]=get_inter(df_rsb25,df_inter)
        ###for Composite pvalue
        d_inter_cp03[pop]=get_inter(df_composite_pvalue03,df_inter)
        #d_inter_cp_mean=get_inter2(df_composite_pvalue,df_inter)#for Cpvaluein the region

    patterns=(f"*_{run}_*.log","*_{:s}".format(run),"*_pruned_{:s}*".format(run),"*_{:s}_p*.vcf".format(run),"*_{:s}.rsb".format(run))
    for pat in patterns:
        for i in glob.glob(pat):
            print(i)
            os.remove(i)
    return run,n_inter,toatl_len,d_num_lines["Freycinet"],d_num_lines["Narawntapu"],d_num_lines["WWP"],d_num_lines2["Freycinet"],d_num_lines2["Narawntapu"],d_num_lines2["WWP"],d_inter_rsb25["Freycinet"],d_inter_rsb25["Narawntapu"],d_inter_rsb25["WWP"],d_inter_cp03["Freycinet"],d_inter_cp03["Narawntapu"],d_inter_cp03["WWP"],df_inter_maf25["Freycinet"],df_inter_maf25["Narawntapu"],df_inter_maf25["WWP"]#n_inter_rsb25

def in_range(x,lows,ups):
    return np.any((x<=ups)&(x>=lows))

if __name__=="__main__":
    d_pop_nSNPs={"Freycinet":4028,"Narawntapu":37407,"WWP":6611}
    vcf_header=open("vcf_header.txt","w")
    vcf_header.write('##fileformat=VCFv4.0\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    vcf_header.close()
    vcf_file="/nfs/users/nfs_q/ql4/lustre_ql4/work/pro_nc_paper/pro360samples/final_vcf/cal_MAF_proportion_3focal/sorted.vcf"
    d_pop_years_st_ed={"Freycinet":{"st":29,"ed":20},\
                 "Narawntapu":{"st":53,"ed":27},\
                 "WWP":{"st":21,"ed":43}}
    d_rate_threshold={"Freycinet":{"st":1/3,"ed":1/2},\
                "Narawntapu":{"st":1/3,"ed":1/3},\
                "WWP":{"st":1/2,"ed":1/3}}
    d_pop_true_samples=original_samples()
    final_file=open("final_inter_count_MAF005.txt","w")
    final_file.write("Run\tN_inter\tLength\tFrey\tNara\tWWP\tFrey_prune\tNara_prune\tWWP_prune\tFrey_Rsb\tNara_Rsb\tWWP_Rsb\tFrey_Cp\tNara_Cp\tWWP_Cp\tFrey_Afc\tNara_Afc\tWWP_Afc\n")
    p=Pool(8)#6
    multi_res=[p.apply_async(main_pro,(run,))for run in range(1000)]
    results=[[res.get()[0],res.get()[1],res.get()[2],res.get()[3],res.get()[4],res.get()[5],res.get()[6],res.get()[7],res.get()[8],res.get()[9],res.get()[10],res.get()[11],res.get()[12],res.get()[13],res.get()[14],res.get()[15],res.get()[16],res.get()[17]] for res in multi_res]
    for item in results:
        final_file.write(f"{item[0]}\t{item[1]}\t{item[2]}\t{item[3]}\t{item[4]}\t{item[5]}\t{item[6]}\t{item[7]}\t{item[8]}\t{item[9]}\t{item[10]}\t{item[11]}\t{item[12]}\t{item[13]}\t{item[14]}\t{item[15]}\t{item[16]}\t{item[17]}\n")
    final_file.close()
