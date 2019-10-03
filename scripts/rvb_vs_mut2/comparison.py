#!/usr/bin/env python

import allel
import pandas as pd
import os
import subprocess
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import prody



"""
Рисуем различные диаграммы Венна
"""



def hapcall_s1_s2(patient):
        
#     fig=plt.figure()
    d={}
    for i in range(1,3):
        path_to_vcfs = f"/media/emir/Storage/LINUX/gatk/gatk_source/Patient_{patient}/file_samples"
#         vcf_name = f"{patient}_t_{i}.bam.g.vcf"
        vcf_name = f"gvcf_list_{i}.list.raw_snp.vcf"
        d["s"+str(i)]= allel.vcf_to_dataframe(os.path.join(path_to_vcfs, vcf_name),
                                              fields="*",
                                              alt_number=1)
    merged = pd.merge(d["s1"], d["s2"], left_on=["CHROM","POS"], 
                      right_on=["CHROM","POS"], how="inner")
    venn2(subsets=(d["s1"].shape[0],d["s2"].shape[0],merged.shape[0]),
          set_labels=(f"Mut2_pat{patient}s1",f"Mut2_pat{patient}s2"))
    plt.savefig(f"/media/emir/Storage/Cancer/mutect/output/pat{patient}/pics/mut2_s1_s2.png",
                dpi=250)
#     plt.clf()



def mutect2_s1_s2(patient):
        
#     fig=plt.figure()
    d={}
    for i in range(1,3):
        path_to_mut_files = f"/media/emir/Storage/Cancer/mutect/output/pat{patient}/s{i}"
        mutect2_name = "mutect2.pkl"
        d["s"+str(i)] = pd.read_pickle(os.path.join(path_to_mut_files,mutect2_name))

    merged = pd.merge(d["s1"], d["s2"], left_on=["CHROM","POS"], 
                      right_on=["CHROM","POS"], how="inner")

    venn2(subsets=(d["s1"].shape[0],d["s2"].shape[0],merged.shape[0]),
          set_labels=(f"Mut2_pat{patient}s1",f"Mut2_pat{patient}s2"))
#     plt.savefig(f"/media/emir/Storage/Cancer/mutect/output/pat{patient}/pics/mut2_s1_s2.png",
#                 dpi=250)
#     plt.clf()
#     return fig

    

def rvboost_s1_s2(patient):
    
#     fig=plt.figure()
    d={}
    for i in range(1,3):
        path_to_rvb_files = f"/media/emir/Storage/Cancer/mutect/rvboost/{patient}s{i}"
        rvboost_name = f"intersect.{patient}_s{i}_0.0.vcf"
        d["s"+str(i)] = allel.vcf_to_dataframe(os.path.join(path_to_rvb_files, rvboost_name),
                                          fields=["CHROM", "POS"],
                                          alt_number=1)
    merged = pd.merge(d["s1"], d["s2"], left_on=["CHROM","POS"], 
                      right_on=["CHROM","POS"], how="inner")
    venn2(subsets=(d["s1"].shape[0],d["s2"].shape[0],merged.shape[0]),
          set_labels=(f"RVB_pat{patient}s1",f"RVB_pat{patient}s2"))
    plt.savefig(f"/media/emir/Storage/Cancer/mutect/output/pat{patient}/pics/rvb_s1_s2.png",
                dpi=250)
    plt.clf()
#     return fig
    

"""    
Делаем таблицу.
Листы вида: 
all_muts_s{sample} -- все найденные mutect2 мутации 
mut2_exclusive_muts_s{sample} -- нашлись только mutect2 
commn_muts_for_rvb_mut2_s{sample} -- общие для выдачи mutect2 и rvboost
"""    
   
    
    
def pfam_annotate(data):
    
    data["PFAM_domains"] = "-"
    data["Link_to_pfam"] = "-"
    for index,row in enumerate(data.iterrows()):
        pfam_domains=[]
        uniprot_id = row[1][18].split(",")
        aa_pos=int(row[1][12])
        for seq in uniprot_id:
            if seq == "-":
                continue
            else:
                pfam_url = "https://pfam.xfam.org/protein/"
                full_url = os.path.join(pfam_url,seq)
                try:
                    f=prody.searchPfam(seq)
                    for i in f.items():
                        start_pos = int(i[1]["locations"][0]["start"])
                        end_pos = int(i[1]["locations"][0]["end"])
                        if aa_pos >= start_pos and aa_pos <= end_pos:
                            pfam_domains.append(str(i[1]["id"])+":"+\
                                                str(start_pos)+"-"+str(end_pos))
                except Exception:
                    continue
        if len(pfam_domains) == 0:
            continue
        else:
            data.loc[index,"PFAM_domains"] = ",".join(pfam_domains)
            data.loc[index,"Link_to_pfam"] = full_url
    return data



def save_mutect2(*args):
    
    patient=args[0]
    sample=args[1]
    
    path_to_mutect2 = f"/media/emir/Storage/Cancer/mutect/output/pat{patient}/s{sample}"
    mutect2_name = f"output_transcr_new_predicted_dbsnp.vcf"
    mutect2_wo_ann=allel.vcf_to_dataframe(os.path.join(path_to_mutect2, mutect2_name),
                                          fields=["numalt"], alt_number=1)
    alt_len = max(mutect2_wo_ann["numalt"])
    col_list=["CHROM","POS","ID","REF","DP","FILTER_PASS","ANN_Annotation",
             "ANN_Annotation_Impact","ANN_Gene_Name","ANN_Gene_ID",
             "ANN_HGVS_c","ANN_HGVS_p","ANN_AA_pos"]
    triple_list=[]
    alt_cols=["ALT","dbNSFP_Polyphen2_HVAR_score","dbNSFP_SIFT_score","dbNSFP_MetaLR_score",
           "dbNSFP_Polyphen2_HDIV_score","dbNSFP_Uniprot_acc","dbNSFP_CADD_phred",
           "dbNSFP_Polyphen2_HDIV_pred","dbNSFP_MutationTaster_score","dbNSFP_SIFT_pred",
           "dbNSFP_MutationTaster_pred","dbNSFP_Polyphen2_HVAR_pred","dbNSFP_MetaLR_pred"]
    for j in alt_cols:
        for i in range(1,alt_len+1):
            triple_list.append(j+f"_{i}")
            col_list.append(j+f"_{i}")
    mutect2_wo_ann=allel.vcf_to_dataframe(os.path.join(path_to_mutect2, mutect2_name),
                                          fields="*", alt_number=alt_len,
                                         exclude_fields="ANN")
    mutect2_w_ann = allel.vcf_to_dataframe(os.path.join(path_to_mutect2, mutect2_name), 
                                           fields="ANN", alt_number=alt_len,
                                           transformers=allel.ANNTransformer())
    mutect2=pd.concat([mutect2_wo_ann,mutect2_w_ann],axis=1)
    mutect2=mutect2[col_list]
    mutect2.fillna("-",inplace=True)
    for i in alt_cols:
        mutect2[i] = mutect2[i+"_1"].map(str)+","+mutect2[i+"_2"].map(str)+","+ \
        mutect2[i+"_3"].map(str)
    mutect2.drop(triple_list,axis=1,inplace=True)
    return mutect2
    mutect2=pfam_annotate(mutect2)
    mutect2.to_pickle(os.path.join(path_to_mutect2,"mutect2.pkl"))
    
    
    
def rvboost_mutect2(*args):
    
    patient=args[0]
    sample=args[1]
    
    path_to_mutect2 = f"/media/emir/Storage/Cancer/mutect/output/pat{patient}/s{sample}"
    mutect2_name = "mutect2.pkl"
    mutect2=pd.read_pickle(os.path.join(path_to_mutect2,mutect2_name))
#     mutect2.query("ANN_Annotation_Impact == \"HIGH\" | ANN_Annotation_Impact == \"MODERATE\"",inplace=True)
    
    try:
        path_to_rvboost = f"/media/emir/Storage/Cancer/mutect/rvboost/{patient}s{sample}"
        rvboost_name = f"intersect.{patient}_s{sample}_0.0.vcf"
        rvboost = allel.vcf_to_dataframe(os.path.join(path_to_rvboost, rvboost_name), 
                                                fields=["CHROM","POS"], alt_number=1)
#         rvboost.query("SNPEFF_IMPACT == \"HIGH\" | SNPEFF_IMPACT == \"MODERATE\"",inplace=True)
        rvb_label=f"RVB_pat{patient}s{sample}"
        subtract = pd.merge(mutect2, rvboost, left_on=["CHROM","POS"], 
                            right_on=["CHROM","POS"], how="outer", indicator=True)
        subtract.query("_merge == \"left_only\"",inplace=True)
        merged = pd.merge(mutect2, rvboost, left_on=["CHROM","POS"],
                          right_on=["CHROM","POS"], how="inner")     
    except OSError:
        path_to_rvboost = f"/media/emir/Storage/Cancer/mutect/rvboost/{patient}_s1s2"
        rvboost_name = f"{patient}.rvboost.s1s2.filtered.tsv"
        rvboost = pd.read_csv(os.path.join(path_to_rvboost, rvboost_name), sep="\t",
                             usecols=["CHROM","START_POSITION"])
        rvb_label=f"RVB_pat{patient}s1s2"
        merged = pd.merge(mutect2, rvboost, left_on=["CHROM","POS"], 
                          right_on=["CHROM","START_POSITION"], how="inner")
        subtract = pd.merge(mutect2, rvboost, left_on=["CHROM","POS"], 
                            right_on=["CHROM","START_POSITION"], how="outer", indicator=True)
        subtract.query("_merge == \"left_only\"",inplace=True)
    
    venn2(subsets=(mutect2.shape[0],rvboost.shape[0],merged.shape[0]),
          set_labels=(f"Mut2_pat{patient}s{sample}",rvb_label))
    plt.savefig(f"/media/emir/Storage/Cancer/mutect/output/pat{patient}/pics/rvb_mut2_s{sample}.png", dpi=250)
    plt.clf()
    
    merged.reset_index(drop=True,inplace=True)
    mutect2.reset_index(drop=True,inplace=True)
    subtract.reset_index(drop=True,inplace=True)
    path_to_xlsx = "/media/emir/Storage/Cancer/mutect/output"
    
    try:
        if args[2]:
            try:
                with pd.ExcelWriter(os.path.join(path_to_xlsx,f"pat{patient}.xlsx"), 
                                    mode="a") as writer:
                    mutect2.to_excel(writer, sheet_name=f"all_muts_s{sample}")
#                     subtract.iloc[:,:subtract.shape[1]-1].to_excel(writer, sheet_name=f"mut2_exclusive_muts_s{sample}")
#                     merged.to_excel(writer, sheet_name=f"commn_muts_for_rvb_mut2_s{sample}")
                    subtract.iloc[:,:subtract.shape[1]-2].to_excel(writer, sheet_name=f"mut2_exclusive_muts_s{sample}")
                    merged.iloc[:,:merged.shape[1]-1].to_excel(writer, sheet_name=f"commn_muts_for_rvb_mut2_s{sample}")
            except OSError:
                with pd.ExcelWriter(os.path.join(path_to_xlsx, f"pat{patient}.xlsx")) as writer:
                    mutect2.to_excel(writer, sheet_name=f"all_muts_s{sample}")
#                     subtract.iloc[:,:subtract.shape[1]-1].to_excel(writer, sheet_name=f"mut2_exclusive_muts_s{sample}")
#                     merged.to_excel(writer, sheet_name=f"commn_muts_for_rvb_mut2_s{sample}")
                    subtract.iloc[:,:subtract.shape[1]-2].to_excel(writer, sheet_name=f"mut2_exclusive_muts_s{sample}")
                    merged.iloc[:,:merged.shape[1]-1].to_excel(writer, sheet_name=f"commn_muts_for_rvb_mut2_s{sample}")
    except IndexError:
        return None
    
    
    
"""
Фильтруем мутации
"""
    
    
    
def filter_merged_muts(*args):
    
    patient=args[0]
    sample=args[1]
    
    df = compare_to_cosmic(patient,sample)

#     df.query("NLOD > TLOD",inplace=True)
    df.query("ANN_Annotation_Impact == \"HIGH\" | ANN_Annotation_Impact == \"MODERATE\"", inplace=True)
    col_list=df.columns.tolist()
    df_new=pd.DataFrame(columns=col_list)
    
    for i in df.iterrows(): 
        cadd=i[1][19].split(",")
        hdiv=i[1][20].split(",")
        sift=i[1][22].split(",")
        taster=i[1][23].split(",")
        hvar=i[1][24].split(",")
        metalr=i[1][25].split(",")
        for a,b,c,d,e,f in zip(cadd,hdiv,sift,taster,hvar,metalr):
            if a == "-":
                continue
            if float(a)>=30.0 or b=="D" or c=="D" or d=="D" or d=="A" or e=="D" or f=="D":
                df_new.loc[i[0]]=i[1]
                break
                
#     df.query("dbNSFP_CADD_phred > 30 | dbNSFP_MutationTaster_pred == \"D\" | dbNSFP_Polyphen2_HDIV_pred == \"D\" | dbNSFP_Polyphen2_HVAR_pred == \"D\" | dbNSFP_MutationTaster_pred == \"A\" | dbNSFP_SIFT_pred == \"D\" | dbNSFP_MetaLR_pred == \"D\"", inplace=True)
    df = df_new.reset_index(drop=True)

    filtered_snps = "filtered_snps.xlsx"
    path_to_mutect2 = "/media/emir/Storage/Cancer/mutect/output"
    path_to_file = os.path.join(path_to_mutect2,filtered_snps)

    try:
        with pd.ExcelWriter(path_to_file, mode="a") as writer:
            df.to_excel(writer, sheet_name=f"pat{patient}s{sample}")
    except OSError:
        with pd.ExcelWriter(path_to_file) as writer:
            df.to_excel(writer, sheet_name=f"pat{patient}s{sample}")


            
"""
Сравниваем результаты mutect2 с кодирующими космиковскими мутациями
"""            
            
    
            
# path_to_cosmic = "/media/emir/Storage/Cancer/mutect/annotation"
# cosmic_name = "CosmicCodingMuts.vcf.gz"
# cosmic = allel.vcf_to_dataframe(os.path.join(path_to_cosmic, cosmic_name), 
#                                 fields="*", alt_number=1)
# cosmic["CHROM"] = "chr"+cosmic["CHROM"]
# cosmic=cosmic[~cosmic.GENE.str.contains("ENST")]
# cosmic.to_pickle(os.path.join(path_to_cosmic, "cosmic.pkl"))

def compare_to_cosmic(*args):
    
    path_to_cosmic = "/media/emir/Storage/Cancer/mutect/annotation"
    cosmic = pd.read_pickle(os.path.join(path_to_cosmic, "cosmic.pkl"))
    cosmic=cosmic[["CHROM", "POS"]]
    
    patient=args[0]
    sample=args[1]
    
    # args[2] -- rvboost
    # no args[2] -- mutect2
    
    try:
        path_to_rvboost = f"/media/emir/Storage/Cancer/mutect/rvboost/{patient}_s1s2"
        rvboost_name = f"{patient}.rvboost.s1s2.filtered.tsv"
        rvboost = pd.read_csv(os.path.join(path_to_rvboost, rvboost_name), sep="\t",
                             usecols=["CHROM","START_POSITION"])
        path_to_mutect2 = f"/media/emir/Storage/Cancer/mutect/output/pat{patient}/s{sample}"
        mutect2_name = "mutect2.pkl"
        mutect2=pd.read_pickle(os.path.join(path_to_mutect2,mutect2_name))
        subtract = pd.merge(mutect2, rvboost, left_on=["CHROM","POS"], 
                            right_on=["CHROM","START_POSITION"], how="outer", indicator=True)
        subtract.query("_merge == \"left_only\"",inplace=True)
        merged = pd.merge(subtract, cosmic, left_on=["CHROM","POS"],
                          right_on=["CHROM","POS"], how="inner")
        merged.query("_merge == \"left_only\"",inplace=True)
        merged.reset_index(drop=True,inplace=True)
        return merged.iloc[:,:merged.shape[1]-2]
    except OSError:
        path_to_rvboost = f"/media/emir/Storage/Cancer/mutect/rvboost/{patient}s{sample}"
        rvboost_name = f"intersect.{patient}_s{sample}_0.0.vcf"
        rvboost = allel.vcf_to_dataframe(os.path.join(path_to_rvboost, rvboost_name), 
                                                fields=["CHROM","POS"], alt_number=1)
        path_to_mutect2 = f"/media/emir/Storage/Cancer/mutect/output/pat{patient}/s{sample}"
        mutect2_name = "mutect2.pkl"
        mutect2=pd.read_pickle(os.path.join(path_to_mutect2,mutect2_name))
        subtract = pd.merge(mutect2, rvboost, left_on=["CHROM","POS"], 
                            right_on=["CHROM","POS"], how="outer", indicator=True)
        subtract.query("_merge == \"left_only\"",inplace=True)
        merged = pd.merge(subtract, cosmic, left_on=["CHROM","POS"],
                          right_on=["CHROM","POS"], how="inner")
        merged.query("_merge == \"left_only\"",inplace=True)
        merged.reset_index(drop=True,inplace=True)
        return merged.iloc[:,:merged.shape[1]-2]
    
    
    
"""
Сколько мутаций общие для космиковских rvboost и mutect2
"""
   
    
    
def compare_cosmic(*args):
    
    patient=args[0]
    sample=args[1]

    path_to_rvboost = f"/media/emir/Storage/Cancer/mutect/rvboost/{patient}s{sample}"
    rvboost_name = f"intersect.{patient}_s{sample}_0.0.vcf"
    rvboost = allel.vcf_to_dataframe(os.path.join(path_to_rvboost, rvboost_name), 
                                            fields="*", alt_number=1)
    rvboost.query("SNPEFF_IMPACT == \"HIGH\" | SNPEFF_IMPACT == \"MODERATE\"",inplace=True)
    merged_rvb = pd.merge(rvboost, cosmic, left_on=["CHROM","POS"],
                          right_on=["CHROM","POS"], how="inner")

    path_to_mutect2 = f"/media/emir/Storage/Cancer/mutect/output/pat{patient}/s{sample}"
    mutect2_name = "mutect2.pkl"
    mutect2=pd.read_pickle(os.path.join(path_to_mutect2,mutect2_name))
    mutect2.query("ANN_Annotation_Impact == \"HIGH\" | ANN_Annotation_Impact == \"MODERATE\"", inplace=True)

    merged_m2 = pd.merge(mutect2, cosmic, left_on=["CHROM","POS"],
                        right_on=["CHROM","POS"], how="inner")

    merged = pd.merge(merged_rvb, merged_m2, left_on=["CHROM","POS"],
                        right_on=["CHROM","POS"], how="inner")
    venn2(subsets=(merged_m2.shape[0],merged_rvb.shape[0],merged.shape[0]),
          set_labels=(f"Mut2_pat{patient}s{sample}",f"rvb_pat{patient}s{sample}"))
#     plt.clf()
#     return mutect2.shape[0], merged.shape[0]
             
            
            
"""
Сохраняем все картинки
"""
    

    
def save_pics(patient):
    
    mutect2_s1_s2(patient)
    rvboost_s1_s2(patient)
    [rvboost_mutect2(patient,sample) for sample in range(1,3)]

    
    
"""
Показываем картинки
"""    
    
    
    
def show_pics(patient):
    
    path_to_pics=f"/media/emir/Storage/Cancer/mutect/output/pat{patient}/pics"
    titles=["mut2", "rvb", "rvb vs mut2 (s1)", "rvb vs mut2 (s2)"]
    names=["mut2_s1_s2.png", "rvb_s1_s2.png", "rvb_mut2_s1.png", "rvb_mut2_s2.png"]
    fig=plt.figure(figsize=(18,12))
    columns = 2
    rows = 2
    for i in range(1, columns*rows+1):
        img = plt.imread(os.path.join(path_to_pics, names[i-1]))
        ax=fig.add_subplot(rows, columns, i)
        ax.set_title(titles[i-1])
        ax.set_xticks([]); ax.set_yticks([])
        plt.imshow(img)

