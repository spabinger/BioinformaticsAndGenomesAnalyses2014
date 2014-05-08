# Practical 1

In this practical you will get to know basic tools for SAM/BAM manipulation and call variants using different programs.




## Required tools

* Annovar
  * http://www.openbioinformatics.org/annovar
  * http://www.openbioinformatics.org/annovar/annovar_download.html
  * databases: (ljb23_all, 1000g2012apr, snp138)
* Freebayes - https://github.com/ekg/freebayes
* GATK - https://www.broadinstitute.org/gatk/download
* IGV - https://www.broadinstitute.org/igv/home
* Java (7) - http://www.oracle.com/technetwork/java/javase/downloads/index.html?ssSourceSiteId=otnjp
* Picard - http://picard.sourceforge.net/command-line-overview.shtml
* Qualimap - http://qualimap.bioinfo.cipf.es/
* SAMtools - http://samtools.sourceforge.net/‎
* Varscan 2 (2.3.6) - http://varscan.sourceforge.net/
* VCFtools - http://vcftools.sourceforge.net/



## Information

* Original data from [1000 genomes](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/NA10847/exome_alignment/)



## Exercise

We assume that we have a properly mapped BAM file from quality checked reads.
For some variant callers we use a target file ADDLINK to shorten variant calling time.

__Important__: After each step inspect the generated output (cat, less, head, grep, ...).


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
    
    Change the heap size as described in http://qualimap.bioinfo.cipf.es/doc_html/faq.html#heapsize
    Open qualimap
    Load the BAM file (BAM QC) -> start analysis ---- Be sure to select the sorted BAM file
    Check out the different pages
    What does the "Coverage across Reference" tells you?
    
__(*)__ Hint:<br/>
If the BAM file is too big, try to view only a subset of it<br/>
(example with 1.000.000 rows [including header])

    samtools view -h sorted.bam | head -n 1000000 | samtools view -b -S - > sorted_small.bam
    
    
    
#### Prepare reference genome
__(*)__ Prepare dict index
    
    java -jar CreateSequenceDictionary.jar R=hg19.fasta O=hg19.dict

__(*)__ Prepare fai index
    
    samtools faidx hg19.fasta 


#### BAM file preparations
__(*)__ Sort with Picard
    
    java -Xmx8g -Dsnappy.disable=true -jar SortSam.jar I=aln.bam O=sorted_picard.bam SORT_ORDER=coordinate


__(*)__ Mark duplicates
     
    java -Xmx8g -Dsnappy.disable=true -jar MarkDuplicates.jar I=sorted_picard.bam O=dedup.bam M=metrics.txt


__(*)__ Add ReadGroup
    
    java -Xmx8g -Dsnappy.disable=true -jar AddOrReplaceReadGroups.jar I= dedup.bam O=deduprg.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample1


__(*)__ Index with Picard
    
    java -Xmx8g -Dsnappy.disable=true -jar BuildBamIndex.jar I=deduprg.bam O=deduprg.bam.bai VALIDATION_STRINGENCY=SILENT

__(*)__ Index with Picard
    
     java -Xmx8g -jar CollectInsertSizeMetrics.jar I=deduprg.bam O=insertSizeHistogram.txt H=insertSizeHistogram.pdf


__(*)__ Questions
* How many reads were marked as duplicated?
* What are the other sorting possibilities for SortSam?


#### SAMtools variant calling

     samtools mpileup -uf hg19.fasta deduprg.bam | bcftools view -vcg - > samtools.vcf

#### FreeBayes variant calling

     ./freebayes -f hg19.fasta -t target.bed -v freebayes.vcf deduprg.bam


#### Varscan variant calling

__(*)__ Convert to mpileup

    samtools mpileup -B -f ${GEN_REF} deduprg.bam > deduprg.pileup

__(*)__ Call SNPs

    java -jar <varscan.jar> mpileup2snp deduprg.pileup --output-vcf 1 > varscan_snp.vcf

__(*)__ Call Indels
    
    java -jar <varscan.jar> mpileup2indel deduprg.pileup --output-vcf 1 > varscan_indel.vcf




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
__(*)__ Using vcfutils
     
     <full-path>/bcftools/vcfutils.pl varFilter -Q 20 -d 5 -D 200 samtools.vcf > samtools_filtered.vcf

__(*)__ Questions
* How many variants were filtered




#### Display files in IGV

    (Download and open) IGV
    Load the BAM file and the VCF files into IGV
    Look at the mapping on Chr 11




#### Annovar
__(*)__ First convert vcf into Annovar format

    <annovar-path>/convert2annovar.pl -format vcf4 freebayes.vcf > freebayes.avinput

__(*)__ Annotate with Gene information TODODOD TEST!
    
    <annovar-path>/annotate_variation.pl -geneanno -buildver hg19 freebayes.avinput /home/stephan/bin/annovar/annovar/humandb/

__(*)__ Annotate with Region information - ljb23

     <annovar-path>/annotate_variation.pl -regionanno -dbtype ljb23_all -buildver hg19 freebayes.avinput /home/stephan/bin/annovar/annovar/humandb/

__(*)__ Annotate with Region information - snp138

     <annovar-path>/annotate_variation.pl -regionanno -dbtype snp138 -buildver hg19 freebayes.avinput /home/stephan/bin/annovar/annovar/humandb/


__(*)__ Questions
* Look at the annotated VCF files.
* What databases does "ljb23" include?



#### SeattleSeq Annotation

__(*)__ Access<br/>
http://snp.gs.washington.edu/SeattleSeqAnnotation138/

__(*)__ Annotate VCF file

    Upload VCF file
    Specify VCF as return type
    submit
    You should receive an annotated VCF file to the specified email address






#### Useful information

__(*)__ Determine number of cores

    cat /proc/cpuinfo

