#!/usr/bin/env python

import subprocess
import os

def filter_somatic(patient, sample):
    
    result=subprocess.run("DIR_FILES=/media/emir/Storage/LINUX/gatk; \
    DIR_OUT=/media/emir/Storage/Cancer/mutect; \
    GATK=$DIR_OUT/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar; \
    java -jar $GATK \
    FilterMutectCalls \
    -R $DIR_FILES/ucsc_hg19.fasta \
    -V $DIR_OUT/output/pat{0}/s{1}/output_transcr_new.vcf \
    -O $DIR_OUT/output/pat{0}/s{1}/output_transcr_new_filtered.vcf".format(patient,sample), 
    shell = True, stderr=subprocess.PIPE)
    
    output_dir=f"/media/emir/Storage/Cancer/mutect/output/pat{patient}/s{sample}/"
    with open(os.path.join(output_dir,"filter_somatic.log"), "w") as s:
        s.write(result.stderr.decode('utf-8'))
