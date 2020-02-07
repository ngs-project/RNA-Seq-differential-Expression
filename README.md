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
 1234  wget https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13
 1235  wget Homo_sapiens.GRCh37.75_subset.fa.gz
 1236  less SRR1039508
 1237  less SRR1039508_subset.sam 
 1238  less SRR1039509_subset.sam 
 1239  fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR1039508
 1240  gunzip -k Homosapiens_GRCh37.fa.gz 
 1241  less SRR1039508_subset.bam 
 1242  wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz
 1243  gunzip -k Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz 
 
 1245  head Homo_sapiens.GRCh38.dna.chromosome.11.fa
 1246  wc chr22_with_ERCC92.fa
 1247  wc Homo_sapiens.GRCh38.dna.chromosome.11.fa
 1248  head -n 425000 Homo_sapiens.GRCh38.dna.chromosome.11.fa | tail
 1249  cat Homo_sapiens.GRCh38.dna.chromosome.11.fa |  grep -v ">" | perl -ne 'chomp $_; $bases{$_}++ for split //; if (eof){print "$_ $bases{$_}\n" for sort keys %bases}'
 1250  cd ..
 1251  cd bamfiles
 1252  cd ..
 1253  ls
 1254  gunzip -k Homo_sapiens.GRCh38.86.gtf.gz 
 1255  conda activate ngs1
 1256  REF_ERCC=./Homo_sapiens.GRCh38.dna.chromosome.11.fa 
 1257  INDEX_ERCC=./Homo_sapiens.GRCh38.dna.chromosome.11
 1258  hisat2-build $REF_ERCC $INDEX_ERCC
 1259  INDEX=~/workdir/diff_exp/ref/ERCC92
 1260  INDEX=~/Diff_proj/index/Homo_sapiens.GRCh38.dna.chromosome.11
 1261  RUNLOG=runlog.txt
 1262  READS_DIR=~/workdir/sample_data/
 1263  READS_DIR=~/Diff_proj/sample_data/

 
 1267  mkdir bam
 1268  > txt.txt
"Step 2 (Alignment)"
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
 hisat2 $INDEX -1 $R1 -2 $R2 | samtools sort > $BAM;    
 samtools index $BAM;  
 done; 
 done

#featureCounts -a $GTF -g gene_name -o counts.txt  Bam_org/UNT*.bam  Bam_org/TTT*.bam
 1272  #GTF=~/workdir/diff_exp/ref/ERCC92.gtf
 1273  GTF=~/Diff_proj/Homo_sapiens.GRCh37.75_subset.gtf 
 1274  featureCounts -a $GTF -g gene_name -o counts.txt  Bam_org/UNT*.bam  Bam_org/TTT*.bam
 1275  conda install ngs1
 1276  conda activate ngs1
 1277  GTF=~/Diff_proj/Homo_sapiens.GRCh37.75_subset.gtf 
 1278  cd Diff_proj/
 1279  ls
 1280  cd 'bamfiles(copy)'/
 1281  ls
 1282  featureCounts -a $GTF -g gene_name -o counts.txt  Bam_org/UNT*.bam  Bam_org/TTT*.bam
 1283  less counts.txt
 1284  cat counts.txt | cut -f 1,7-12 > simple_counts.txt
 1285  less simple_counts.txt 
 1286  less counts.txt
 1287  c
 1288  hea results_deseq1.tsv 
 1289  head results_deseq1.tsv 
 1290  conda activate ngs1
 1291  cat simple_counts.txt | Rscript deseq1.r 3x3 > results_deseq1.tsv
 1292  cat results_deseq1.tsv | awk ' $8 < 0.05 { print $0 }' > filtered_results_deseq1.tsv
 1293  cat filtered_results_deseq1.tsv | Rscript draw-heatmap.r > hisat_output.pdf
 1294  head filtered_results_deseq1.tsv 

 

 
 
 
 

