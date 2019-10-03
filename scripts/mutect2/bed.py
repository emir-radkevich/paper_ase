#!/usr/bin/env python

import allel
import pandas as pd
import os
from tqdm import tqdm_notebook as tqdm
import numpy as np
from pyensembl import EnsemblRelease
# в конде:
# pyensembl install --release 75 --species human


data = EnsemblRelease(75)

global counter
counter=0

def make_bed_file(*args):
    
    patient=args[0]
    a=np.empty([1,3])
    global bed_pat
    global counter
    if counter == 0: bed_pat=patient
    
    try:
        try:
            sample=args[1]
            path_to_files = f"/media/emir/Storage/Cancer/mutect/rvboost/{patient}s{sample}"
            file_name = f"intersect.{patient}_s{sample}_0.0.vcf"
            output_name=f"genes_pat{patient}s{sample}_ensembl.bed"
            callset=allel.read_vcf(os.path.join(path_to_files, file_name),
                                   fields=["CHROM","POS"])
            for chrom,pos in tqdm(zip(callset["variants/CHROM"], callset["variants/POS"])):
                gene_obj=data.genes_at_locus(contig=chrom, position=pos.item())
                for obj in gene_obj:
                    bed_row=["chr"+str(obj.contig),obj.start,obj.end]
                    a=np.vstack([a,bed_row])

        except FileNotFoundError:
            if bed_pat == patient and counter != 0:
                return
            else:
                path_to_files = f"/media/emir/Storage/Cancer/mutect/rvboost/{patient}_s1s2"
                file_name = f"{patient}.rvboost.s1s2.filtered.tsv"
                output_name=f"genes_pat{patient}_s1s2_ensembl.bed"
                callset = np.loadtxt(fname=os.path.join(path_to_files,file_name),
                                     delimiter="\t", skiprows=1,usecols=(1,2),dtype=str)
                for chrom_pos in tqdm(callset):
                    contig=str(chrom_pos[0]).split("chr")[1]
                    gene_obj=data.genes_at_locus(contig=contig, 
                                                 position=int(chrom_pos[1]))
                    for obj in gene_obj:
                        bed_row=["chr"+str(obj.contig),obj.start,obj.end]
                        a=np.vstack([a,bed_row])
    except OSError:
        return
        
        
    a=np.delete(a,(0),axis=0)
    df=pd.DataFrame(data=a[0:,0:], columns=["CHROM","START","END"])
    df.drop_duplicates(inplace=True)
    
    arr=df.to_numpy()
    arr[:,1:]=arr[:,1:].astype(int)
    arr

    d={}
    for i in arr:
        if i[0] not in d.keys():
            d[i[0]]=[[i[1],i[2]]]
        else:
            d[i[0]].append([i[1],i[2]])

    with open(os.path.join(path_to_files, output_name),
              "w") as genes_str:        
        for k,v in d.items():
            v.sort(key=lambda interval: interval[0])
            merged = [v[0]]
            for current in v:
                previous = merged[-1]
                if current[0] <= previous[1]:
                    previous[1] = max(previous[1], current[1])
                else:
                    merged.append(current)
            for i in merged:
                genes_str.write(k+"\t"+str(i[0])+"\t"+str(i[1])+"\n")

    counter+=1
    bed_pat=patient
                

