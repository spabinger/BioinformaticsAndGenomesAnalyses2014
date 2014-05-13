# Practical 1

In this practical you will get to know basic tools for SAM/BAM manipulation and call variants using different programs. Please check out the [Useful information](#useful-information) section.

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

* Original data from [1000 genomes](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/NA10847/exome_alignment/).



## Exercise

We assume that we have a properly mapped BAM file from quality checked reads.
For some variant callers we use a [target file](target.bed) to shorten variant calling time.

#### Important

* After each step inspect the generated output (cat, less, head, grep, ...).
* Organize your data and results in folders.


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
__(*)__ Run qualimap on the command line
    
    ./qualimap bamqc -bam sorted.bam -nt <numberOfThreads> -outdir bamqc


__(*)__ Inspect the BAM file in qualimap

    Check out the generated report. What qc features does it include?

    
__(*)__ Alternatively start the graphical interface and inspect the file using the GUI
    
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
    
    java -Xmx8g -jar SortSam.jar I=aln.bam O=sorted_picard.bam SORT_ORDER=coordinate


__(*)__ Mark duplicates
     
    java -Xmx8g -jar MarkDuplicates.jar I=sorted_picard.bam O=dedup.bam M=metrics.txt


__(*)__ Add ReadGroup
    
    java -Xmx8g -jar AddOrReplaceReadGroups.jar I= dedup.bam O=deduprg.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample1


__(*)__ Index with Picard
    
    java -Xmx8g -jar BuildBamIndex.jar I=deduprg.bam O=deduprg.bam.bai VALIDATION_STRINGENCY=SILENT

__(*)__ Index with Picard
    
     java -Xmx8g -jar CollectInsertSizeMetrics.jar I=deduprg.bam O=insertSizeHistogram.txt H=insertSizeHistogram.pdf


__(*)__ Questions
* How many reads were marked as duplicated?
* What are the other sorting possibilities for SortSam?


#### SAMtools variant calling

__(*)__ Call

     samtools mpileup -uf hg19.fasta deduprg.bam | bcftools view -vcg - > samtools.vcf

__(*)__ Investigate result

    #How many variant were called
    grep -v "^#" samtools.vcf | wc -l
    #Print the variant that are between 1-300000 
    awk '!/^#/ && $2 < "300000"' samtools.vcf

#### FreeBayes variant calling

__(*)__ Call

     ./freebayes -f hg19.fasta -t target.bed -v freebayes.vcf deduprg.bam

__(*)__ Investigate result
  
    #Perform the same procedures as done for samtools
    #Do you notice differences?





#### Useful information

__(*)__ Determine number of cores

    cat /proc/cpuinfo  

__(*)__ Make file executable

    chmod +x <file.name>
    
__(*)__ Extract information from a file

    grep -i "info" <file.name>
    
__(*)__ Extract information from a file excluding the pattern and display the first 100 lines

    grep -v "^#" <file.name> | head -n 100


