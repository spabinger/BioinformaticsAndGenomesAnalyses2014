# Practical 1

In this practical you will XXXX




## Required tools

* SAMtools - http://samtools.sourceforge.net/‎
* GATK - https://www.broadinstitute.org/gatk/download
* Picard - http://picard.sourceforge.net/command-line-overview.shtml
* Annovar
  * http://www.openbioinformatics.org/annovar
  * http://www.openbioinformatics.org/annovar/annovar_download.html
  * databases: (ljb23_all, 1000g2012apr, snp138)
* Qualimap - http://qualimap.bioinfo.cipf.es/
* Java (7) - http://www.oracle.com/technetwork/java/javase/downloads/index.html?ssSourceSiteId=otnjp
* IGV - https://www.broadinstitute.org/igv/home
* Freebayes - https://github.com/ekg/freebayes
* VCFtools - http://vcftools.sourceforge.net/
* Varscan 2 (2.3.6) - http://varscan.sourceforge.net/


## Information

* Original data from [1000 genomes](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/NA10847/exome_alignment/)



## Exercise

We assume that we have a properly mapped BAM file from quality checked reads.
For some variant callers we use a target file ADDLINK to shorten variant calling time.

__Important__ After each step inspect the generated output (cat, less, head, grep, ...).


#### SAMtools


__(*)__ How big is the BAM file

    ls -lah <file.bam>

__(*)__ View the BAM file

    samtools view <bam.file> | less
    
__(*)__ How many reads are in the BAM file?<br/>
Is there another way to count the reads (check the samtools view parameters - look for -v)
   
    samtools view <file.bam> | grep -v "^#" | wc -l
    
__(*)__ Answer the following questions by investigating the SAM file
* What version of the human assembly was used to perform the alignments?
* What version of bwa was used to align the reads?
* What is the name of the first read?
* At what position does the alignment of the read start?

    
__(*)__ Sort the BAM file

    samtools sort <file.bam> <sorted>
    
__(*)__ Index the bam file
    
    samtools index <sorted.bam>


#### Alignment stats
    samtools flagstat sorted.bam
    samtools idxstats sorted.bam



#### Qualimap
__(*)__ Inspect the BAM file in qualimap
    
    Open qualimap
    Load the BAM file (BAM QC) -> start analysis
    
    
#### Prepare reference genome
__(*)__ Prepare dict index
    
    java -jar CreateSequenceDictionary.jar R=hg19.fasta O=hg19.dict

__(*)__ Prepare fai index
    
    samtools faidx hg19.fasta 


#### BAM file preparations
__(*)__ Sort with Picard
    
    java -Xmx24g -Dsnappy.disable=true -jar SortSam.jar I=aln.bam O=sorted_picard.bam SORT_ORDER=coordinate


__(*)__ Mark duplicates
     
    java -Xmx24g -Dsnappy.disable=true -jar MarkDuplicates.jar I=sorted_picard.bam O=dedup.bam M=metrics.txt


__(*)__ Add ReadGroup
    
    java -Xmx24g -Dsnappy.disable=true -jar AddOrReplaceReadGroups.jar I= dedup.bam O=deduprg.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample1


__(*)__ Index with Picard
    
    java -Xmx24g -Dsnappy.disable=true -jar BuildBamIndex.jar I=deduprg.bam O=deduprg.bam.bai VALIDATION_STRINGENCY=SILENT


__(*)__ Questions
* How many reads were marked as duplicated?
* What are the other sorting possibilities for SortSam?


#### SAMtools variant calling

     samtools mpileup -uf hg19.fasta deduprg.bam | bcftools view -vcg - > samtools.vcf

#### FreeBayes variant calling

     ./freebayes -f hg19.fasta -t target.bed -v freebayes.vcf deduprg.bam



#### GATK variant calling

__(*)__ Known indel sites are here specified as variables - either copy the whole path or use variables as well

    KNOWN_INDELS_1="1000G_phase1.indels.hg19.vcf"
    KNOWN_INDELS_2="Mills_and_1000G_gold_standard.indels.hg19.vcf"


__(*)__ Realignment target creator

    java -Xmx8g -Dsnappy.disable=true -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R hg19.fasta -nt 8 -L target.bed -I deduprg.bam -known ${KNOWN_INDELS_1} -known ${KNOWN_INDELS_2} -o target_intervals.list

__(*)__ Perform realignment
    
    java -Xmx8g -Dsnappy.disable=true -jar GenomeAnalysisTK.jar -T IndelRealigner -R hg19.fasta -I deduprg.bam -targetIntervals target_intervals.list -known ${KNOWN_INDELS_1} -known ${KNOWN_INDELS_2} -o dedup_rg_real.bam


__(*)__ Base quality recalibration
    
    java -Xmx8g -Dsnappy.disable=true -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R hg19.fasta -I dedup_rg_real.bam -knownSites ${KNOWN_INDELS_1} -knownSites ${KNOWN_INDELS_2} -o recal_data_table.txt -L target.bed --maximum_cycle_value 800


__(*)__ Second pass of recalibration
     
     java -Xmx8g -Dsnappy.disable=true -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R hg19.fasta -I dedup_rg_real.bam -knownSites ${KNOWN_INDELS_1} -knownSites ${KNOWN_INDELS_2} -o post_recal_data_table.txt -L target.bed --maximum_cycle_value 800 -BQSR recal_data_table.txt 


__(*)__ Generate before after plots (requires R and ggplot2)
    
    java -Xmx24g -Dsnappy.disable=true -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R hg19.fasta -L target.bed -before recal_data_table.txt -after post_recal_data_table.txt -plots recalibration_plots.pdf



__(*)__ Print recalibrated reads
    
    java -Xmx24g -Dsnappy.disable=true -jar GenomeAnalysisTK.jar -T PrintReads -R hg19.fasta -L target.bed -I deduprgreal.bam -BQSR recal_data_table.txt -o dedup_rg_real_recal.bam


__(*)__ Now do variant calling
    
    java -Xmx24g -Dsnappy.disable=true -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R hg19.fasta -nct 8 -L target.bed -I dedup_rg_real_recal.bam --genotyping_mode DISCOVERY -o gatk.vcf

__(*)__ Questions
* Check out the before and after plots
* 



#### Merge VCFs and VCF stats

__(*)__ VCFlib - merge

    vcfcombine freebayes.vcf gatk.vcf samtools.vcf > vcf_lib_merged.vcf

__(*)__ VCFlib - stats - shown here for one VCF file - repeat for all 3

    vcfstats freeb_call.vcf > freeb_call_vcf_lib_stats.txt



__(*)__ VCFtools

    export PERL5LIB=<full-path>/vcftools_0.1.12a/perl
    export PATH=<full-path>/tabix/tabix-0.2.6:$PATH

    ## Index (tabix) and pack files
    cp gatk.vcf gatk_tab.vcf
    bgzip gatk_tab.vcf
    tabix -p vcf gatk_tab.vcf.gz

    ## repeat for the other two VCF files

    vcf-merge freebayes_tab.vcf.gz gatk_tab.vcf.gz samtools_tab.vcf.gz > vcf_tools_merged.vcf
    vcf-stats freebayes_tab.vcf.gz > freebayes_tab_stats.txt





#### Filter variants
     
     <full-path>/bcftools/vcfutils.pl varFilter -Q 20 -d 5 -D 200 samtools.vcf > samtools_filtered.vcf

__(*)__ Questions
* How many variants were filtered







insert size metrics



Transform SAM to BAM

transform to bam
samtools'view'–Sb'input.sam'>'tempfile.bam'





Filtering%SAM/BAM%ﬁles%

Required%ﬂag%(keep%if%matches)%
samtools'view'–f'

Filtering%ﬂag%(remove%if%matches)
samtools'view'–F'


Samtools variant calling

The BCF doesn't hold actual calls
∙ encodes likelihoods for all variants

GATK variant calling

freebayes variant calling

varscan variant calling



cnv

*) mendelscan



structural variations


Filtering of variants – VCF (one liners!)



Annovar

Web - SeattleAnnotation



IGV display files


ReadXplorer

Download the ReadXplorer from
http://www.uni-giessen.de/fbz/fb08/bioinformatik/software/ReadXplorer/access

Unzip and start it



vcf merge












Determine number of cores

cat /proc/cpuinfo


VCF Tools (vcf-annotate)
 Soft filter variants file for these biases
 Variants kept in the file – just annotated with potential bias affecting the
variant





How do the merged vcf files differ -> check the online documentation

*) Multiple cores
