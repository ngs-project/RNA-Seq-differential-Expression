# RNA-Seq-differential-Expression
 we perform  data analysis for quality assessment relationship we between samples, perform differential gene expression analysis, and visually explore the results.

#Download reference Genomes
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
  ls
  wget https://download1.rstudio.org/rstudio-xenial-1.1.419-amd64.deb
  sudo gdebi rstudio-xenial-1.1.379-amd64.deb
  conda activate ngs1
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/003/SRR1039523/SRR1039523_1.fastq.gz
 
  sudo apt-get install libopenblas-base r-base
  sudo apt-get install gdebi
  cd ~/Downloads
  wget https://download1.rstudio.org/rstudio-xenial-1.1.419-amd64.deb
  sudo gdebi rstudio-xenial-1.1.379-amd64.deb
  wget https://download1.rstudio.org/rstudio-xenial-1.1.419-amd64.deb
  sudo gdebi rstudio-xenial-1.1.379-amd64.deb
  ls
  sudo gdebi rstudio-xenial-1.1.419-amd64.deb
   R
  {r install-tidyverse, eval = F}
install.packages("tidyverse", repos = 'https://cran.us.r-project.org')

##getting subsets of files as the other one is very large 
       tar xvzf airway_1.6.0.tar.gz 

##the files are in BAM already so convert it to Sam 
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

 #convert bam file to fastq
   for file in ./*.bam;
      do 
        echo $file ; samtools bam2fq *.bam > *.fastq
      done
     
# extracting reads ending with '/1' or '/2'
    for file in ./*.fastq
     do
       cat ./*.fastq | grep '^@.*/1$' -A 3 --no-group-separator > r1.fastq
       cat ./*.fastq | grep '^@.*/2$' -A 3 --no-group-separator > r2.fastq
     done
