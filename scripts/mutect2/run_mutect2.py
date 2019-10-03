#!/usr/bin/env python


import subprocess
import os

def search_mutations(*args):
    
    patient=args[0]
    sample=args[1]
   
    try:
        if args[2]:
            result=subprocess.run("BAM=/media/emir/Seagate/data/patients; \
            DIR_FILES=/media/emir/Storage/LINUX/gatk; \
            DIR_OUT=/media/emir/Storage/Cancer/mutect; \
            GATK=$DIR_OUT/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar; \
            java -jar $GATK \
            GetSampleName \
            -I $BAM/Patient_{0}/Tumor_{1}/reorder.sorted.bam \
            -O $BAM/Patient_{0}/Tumor_{1}/name; \
            name_tumor=$(less $BAM/Patient_{0}/Tumor_{1}/name); \
            rm $BAM/Patient_{0}/Tumor_{1}/name; \
            java -jar $GATK \
            GetSampleName \
            -I $BAM/Patient_{0}/Non-Tumor_{1}/reorder.sorted.bam \
            -O $BAM/Patient_{0}/Non-Tumor_{1}/name; \
            name_normal=$(less $BAM/Patient_{0}/Non-Tumor_{1}/name); \
            rm $BAM/Patient_{0}/Non-Tumor_{1}/name; \
            java -jar $GATK \
            Mutect2 \
            -R $DIR_FILES/ucsc_hg19.fasta \
            -I $BAM/Patient_{0}/Tumor_{1}/reorder.sorted.bam \
            -I $BAM/Patient_{0}/Non-Tumor_{1}/reorder.sorted.bam \
            -normal $name_normal \
            -tumor $name_tumor \
            -L $DIR_OUT/rvboost/{0}_s1s2/genes_pat{0}_s1s2_ensembl.bed \
            -O $DIR_OUT/output/pat{0}/s{1}/output_transcr_new.vcf".format(patient,sample), 
            shell = True, stderr=subprocess.PIPE)

    except IndexError:
        result=subprocess.run("BAM=/media/emir/Seagate/data/patients; \
        DIR_FILES=/media/emir/Storage/LINUX/gatk; \
        DIR_OUT=/media/emir/Storage/Cancer/mutect; \
        GATK=$DIR_OUT/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar; \
        java -jar $GATK \
        GetSampleName \
        -I $BAM/Patient_{0}/Tumor_{1}/reorder.sorted.bam \
        -O $BAM/Patient_{0}/Tumor_{1}/name; \
        name_tumor=$(less $BAM/Patient_{0}/Tumor_{1}/name); \
        rm $BAM/Patient_{0}/Tumor_{1}/name; \
        java -jar $GATK \
        GetSampleName \
        -I $BAM/Patient_{0}/Non-Tumor_{1}/reorder.sorted.bam \
        -O $BAM/Patient_{0}/Non-Tumor_{1}/name; \
        name_normal=$(less $BAM/Patient_{0}/Non-Tumor_{1}/name); \
        rm $BAM/Patient_{0}/Non-Tumor_{1}/name; \
        java -jar $GATK \
        Mutect2 \
        -R $DIR_FILES/ucsc_hg19.fasta \
        -I $BAM/Patient_{0}/Tumor_{1}/reorder.sorted.bam \
        -I $BAM/Patient_{0}/Non-Tumor_{1}/reorder.sorted.bam \
        -normal $name_normal \
        -tumor $name_tumor \
        -L $DIR_OUT/rvboost/{0}s{1}/genes_pat{0}s{1}_ensembl.bed \
        -O $DIR_OUT/output/pat{0}/s{1}/output_transcr_new.vcf".format(patient,sample), 
        shell = True, stderr=subprocess.PIPE)
    
    output_dir=f"/media/emir/Storage/Cancer/mutect/output/pat{patient}/s{sample}/"
    with open(os.path.join(output_dir,"search_somatic.log"), "w") as s:
        s.write(result.stderr.decode('utf-8'))
