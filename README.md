# RNA-Seq-differential-Expression
 we perform  data analysis for quality assessment relationship we between samples, perform differential gene expression analysis, and visually explore the results.

Download reference Genomes
Obtain a reference genome from Ensembl, iGenomes, NCBI or UCSC. In this example analysis we will use the human GRCh38 version of the genome from Ensembl.  

# getting the ref in gtf format 
  wget http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
 
 # Data retrieval for fastq files :
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_1.fastq.gz
  conda activate ngs1
 prefetch SRR1039508
 fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR1039508
 prefetch SRR1039509
 fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR1039509
 
 
# Setup enviornemnt (preparing R)

  conda activate ngs1
  
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/003/SRR1039523/SRR1039523_1.fastq.gz
 
  sudo apt-get install libopenblas-base r-base
  
  sudo apt-get install gdebi
  
 wget https://download1.rstudio.org/rstudio-xenial-1.1.419-amd64.deb
  
  sudo gdebi rstudio-xenial-1.1.379-amd64.deb
  
   R

{r install-tidyverse, eval = F}

install.packages("tidyverse", repos = 'https://cran.us.r-project.org')

##getting subsets of files as the other one is very large 
       tar xvzf airway_1.6.0.tar.gz 

# the files are in BAM already so convert it to Sam 
 
 for file in ./*.bam;
   do    
       echo $file ;     samtools view -h $file > ${file/.bam/.sam};
   done
 
  conda activate ngs1
  samtools

for file in ./*.bam;
     do  
       echo $file ; samtools view -h $file > ${file/.bam/.sam};
     done

 ##convert bam file to fastq
   for file in ./*.bam;
      do 
        echo $file ; samtools bam2fq *.bam > *.fastq
      done
     
##extracting reads ending with '/1' or '/2'
    for file in ./*.fastq
     do
       cat ./*.fastq | grep '^@.*/1$' -A 3 --no-group-separator > r1.fastq
       cat ./*.fastq | grep '^@.*/2$' -A 3 --no-group-separator > r2.fastq
     done
   
#initiating githup repo
   cd git_one/
   git init

# downloading ref genome , Trying to get the ref genome for the specifyed regions
  wget https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13
  wget Homo_sapiens.GRCh37.75_subset.fa.gz
 
   less SRR1039508_subset.sam 
  less SRR1039509_subset.sam 
   fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR1039508
  gunzip -k Homosapiens_GRCh37.fa.gz 
   less SRR1039508_subset.bam 
  wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz
   gunzip -k Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz 
 
  head Homo_sapiens.GRCh38.dna.chromosome.11.fa
   wc chr22_with_ERCC92.fa
  wc Homo_sapiens.GRCh38.dna.chromosome.11.fa
   head -n 425000 Homo_sapiens.GRCh38.dna.chromosome.11.fa | tail
  cat Homo_sapiens.GRCh38.dna.chromosome.11.fa |  grep -v ">" | perl -ne 'chomp $_; $bases{$_}++ for split //; if (eof){print "$_ $bases{$_}\n" for sort keys %bases}'

  
 # Step 1 Indexing
  REF_ERCC=./Homo_sapiens.GRCh38.dna.chromosome.11.fa 
  INDEX_ERCC=./Homo_sapiens.GRCh38.dna.chromosome.11

 hisat2-build $REF_ERCC $INDEX_ERCC
 
 INDEX=~/workdir/diff_exp/ref/ERCC92
 NDEX=~/Diff_proj/index/Homo_sapiens.GRCh38.dna.chromosome.11

RUNLOG=runlog.txt

READS_DIR=~/workdir/sample_data/
 READS_DIR=~/Diff_proj/sample_data/

 
# Step 2 (Alignment)
 for SAMPLE in SRR508 SRR509 SRR512 SRR513 SRR516 SRR517 SRR520 SRR521 ; 
    do    
 for REPLICATE in 1 2 3 4 5 6 7 8 ;    
    do         
 R1=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*read1.fastq.gz;       
 R2=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*read2.fastq.gz;       
 BAM=bam/${SAMPLE}_${REPLICATE}.bam;         
 hisat2 $INDEX -1 $R1 -2 $R2 | samtools sort > $BAM;      
 samtools index $BAM;  
    done; 
    done
 
 for SAMPLE in UNT TTT ; 
 do    
 for REPLICATE in 08 09 12 13 16 17 20 21;    
 do       
 R1=$READS_DIR/${SAMPLE}_SRR5${REPLICATE}*r1.fastq;     
 R2=$READS_DIR/${SAMPLE}_SRR5${REPLICATE}*r2.fastq;        
 BAM=bam/${SAMPLE}_${REPLICATE}.bam;        
 hisat2 $INDEX -1 $R1 -2 $R2 | samtools sort > $BAM
 samtools index $BAM 
 done; 
 done

# Step 3 (Quantification)
 
   #GTF=~/workdir/diff_exp/ref/ERCC92.gtf
   GTF=~/Diff_proj/Homo_sapiens.GRCh37.75_subset.gtf 

  conda install ngs1
 conda activate ngs1
cd Diff_proj/ 'bamfiles(copy)'/
  
 
 ##Generate the counts.
  featureCounts -a $GTF -g gene_name -o counts.txt  Bam_org/UNT*.bam  Bam_org/TTT*.bam
##Simplify the file to keep only the count columns.
  cat counts.txt | cut -f 1,7-12 > simple_counts.txt
  Head
   head results_deseq1.tsv 
  conda activate ngs1
  # Analyze the counts with DESeq1.
   cat simple_counts.txt | Rscript deseq1.r 3x3 > results_deseq1.tsv
   # View only rows with pval < 0.05
  cat results_deseq1.tsv | awk ' $8 < 0.05 { print $0 }' > filtered_results_deseq1.tsv
  cat filtered_results_deseq1.tsv | Rscript draw-heatmap.r > hisat_output.pdf
  head filtered_results_deseq1.tsv 

 

 
 
 
 

