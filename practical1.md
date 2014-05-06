# Practical 1

In this practical you will XXXX




## Required tools

* SAMtools - samtools.sourceforge.net/‎
* GATK - https://www.broadinstitute.org/gatk/download
* Picard - http://picard.sourceforge.net/command-line-overview.shtml
* Annovar
  * http://www.openbioinformatics.org/annovar
  * http://www.openbioinformatics.org/annovar/annovar_download.html
  * databases: (ljb23_all, 1000g2012apr, snp138)
* Qualimap - http://qualimap.bioinfo.cipf.es/
* Java (7)
* IGV - https://www.broadinstitute.org/igv/home
* Freebayes - https://github.com/ekg/freebayes
* VCFtools - http://vcftools.sourceforge.net/
* Varscan 2 (2.3.6) - http://varscan.sourceforge.net/


## Data

* Original data from [1000 genomes](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/NA10847/exome_alignment/)


## More information
XXXXX




## Exercise

We assume that we have a properly mapped BAM file from quality check reads. 


#### SAMtools


How big is the BAM file

    ls -lah <file.bam>

View BAM file

    samtools view <bam.file> | less
    
How many reads are in the BAM file?<br/>
Is there another way to count the reads (check the samtools view parameters - look for -v)
   
    samtools view <file.bam> | grep -v "^#" | wc -l
    
Sort BAM file

    samtools sort <file.bam> <sorted>
    
Index bam file
    
    samtools index <sorted.bam>



#### Qualimap




#### BAM file preparations




#### Mark Duplicates




how many reads removed



#### SAMtools variant calling






#### GATK variant calling




#### FreeBayes variant calling



#### Merge VCFs


#### Compare VCFs ??? HOW???









samtools/picard





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

vcf merge


*) Multiple cores
