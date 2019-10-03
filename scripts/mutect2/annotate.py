#!/usr/bin/env python

import subprocess
import os

def db():
    result = subprocess.run("DIR_OUT=/media/emir/Storage/Cancer/mutect; \
    SNPEFF=$DIR_OUT/annotation/snpEff/snpEff.jar; \
    CONFIG=$DIR_OUT/annotation/snpEff/snpEff.config; \
    java -jar $SNPEFF download \
    -c $CONFIG GRCh37.75", shell=True, stderr=subprocess.PIPE)
    return print(result.stderr.decode("utf-8"))


def annotate_variants(patient, sample):
     
    result = subprocess.run("DIR_OUT=/media/emir/Storage/Cancer/mutect; \
    MAIN_FOLDER=$DIR_OUT/annotation/snpEff; \
    CONFIG=$DIR_OUT/annotation/snpEff/snpEff.config; \
    cd $MAIN_FOLDER; \
    java -jar \"-Xmx4g\" snpEff.jar GRCh37.75 \
    $DIR_OUT/output/pat{0}/s{1}/output_transcr_new_filtered.vcf > \
    $DIR_OUT/output/pat{0}/s{1}/output_transcr_new_annotated.vcf".format(patient,sample),shell=True, stderr=subprocess.PIPE)
    output_dir=f"/media/emir/Storage/Cancer/mutect/output/pat{patient}/s{sample}/"
    with open(os.path.join(output_dir,"annotate_variants.log"), "w") as s:
        s.write(result.stderr.decode('utf-8'))


def func_pred(patient, sample):
     
    result = subprocess.run("DIR_OUT=/media/emir/Storage/Cancer/mutect; \
    SNPSIFT=$DIR_OUT/annotation/snpEff/SnpSift.jar; \
    CONFIG=$DIR_OUT/annotation/snpEff/snpEff.config; \
    DB_FOLDER=`dirname $SNPSIFT`; \
    java -jar $SNPSIFT dbnsfp \
    -c $CONFIG -v \
    -db $DB_FOLDER/dbNSFP2.9.txt.gz \
    -f CADD_phred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,SIFT_score,SIFT_pred,MutationTaster_score,MutationTaster_pred,MetaLR_score,MetaLR_pred,Uniprot_acc \
    $DIR_OUT/output/pat{0}/s{1}/output_transcr_new_annotated.vcf > \
    $DIR_OUT/output/pat{0}/s{1}/output_transcr_new_predicted.vcf".format(patient,sample),shell=True, stderr=subprocess.PIPE)
    output_dir=f"/media/emir/Storage/Cancer/mutect/output/pat{patient}/s{sample}/"
    with open(os.path.join(output_dir,"functional_prediction.log"), "w") as s:
        s.write(result.stderr.decode('utf-8'))
        
        
def make_dbsnp_file():
     
    result = subprocess.run("DATA=/media/emir/Seagate/data/; \
    gawk '{OFS=\"\t\"} \
    {if ($1 ~ /^#/) {print $0;} \
    else {$1=\"chr\"$1; print $0;}}' \
    $DATA/dbsnp_151_37.vcf > $DATA/dbsnp_151_37_w_chr_new.vcf",shell=True, stderr=subprocess.PIPE)
        
        
def dbsnp_annotate(patient, sample):
     
    result = subprocess.run("DIR_OUT=/media/emir/Storage/Cancer/mutect; \
    DIR_FILES=/media/emir/Storage/LINUX/gatk; \
    DBSNP=/media/emir/Seagate/data/dbsnp_151_37_w_chr_new.vcf.gz; \
    GATK=$DIR_OUT/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar; \
    java -jar $GATK VariantAnnotator \
    -R $DIR_FILES/ucsc_hg19.fasta \
    -V $DIR_OUT/output/pat{0}/s{1}/output_transcr_new_predicted.vcf \
    -D $DBSNP -O $DIR_OUT/output/pat{0}/s{1}/output_transcr_new_predicted_dbsnp.vcf".format(patient,sample),shell=True, stderr=subprocess.PIPE)
    output_dir=f"/media/emir/Storage/Cancer/mutect/output/pat{patient}/s{sample}/"
    with open(os.path.join(output_dir,"dbsnp_annotate.log"), "w") as s:
        s.write(result.stderr.decode('utf-8'))
        
