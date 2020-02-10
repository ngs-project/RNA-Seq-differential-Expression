
Aim :using the bioinformatic tools investigating the outbreak by assembling the genome of the deadly E. coli X strain. Specifically, we will provide you with Illumina reads from the TY2482 sample, which were generated at Beijing Genome Institute and deposited into the Short Read Archive (SRA) for public access.
 
 *[Data Retrieval](**url**)**
NCBI’s fastq-dump from sra-toolkit was used to download the short reads for NCBI short read archive (SRA).

#Using SRA-toolkit
```
 prefetch SRR292678
  fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR292678.sra

```
**[Prepare the referance data](url)**
```
mkdir -p ~/workdir/sample_data
   cd ~/workdir/sample_data

wget GCF_000221885.1_E.coli_0104_H4_Illumina_1.0_cds_from_genomic.fna.gz
 cp GCF_000221885.1_E.coli_0104_H4_Illumina_1.0_genomic.fna GCF_000221885.1_E.coli_0104_H4_Illumina_1.0_genomic.fa
 
 ` wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/221/885/GCF_000221885.1_E.coli_0104_H4_Illumina_1.0/GCF_000221885.1_E.coli_0104_H4_Illumina_1.0_genomic.gtf
```
[**Indexing**](url)
****```
  mkdir -p ~/workdir/hisat_align/hisatIndex 
 cd ~/workdir/hisat_align/hisatIndex
   
ln -s ~/workdir/sample_data/GCF_000221885.1_E.coli_0104_H4_Illumina_1.0/GCF_000221885.1_E.coli_0104_H4_Illumina_1.0_genomic.fa
  
 hisat2_extract_splice_sites.py ~/workdir/sample_data/GCF_000221885.1_E.coli_0104_H4_Illumina_1.0/GCF_000221885.1_E.coli_0104_H4_Illumina_1.0_genomic.gtf
  
  hisat2_extract_splice_sites.py ~/workdir/sample_data/GCF_000221885.1_E.coli_0104_H4_Illumina_1.0/GCF_000221885.1_E.coli_0104_H4_Illumina_1.0_genomic.gtf > splicesites.tsv
   
  less SGCF-000221885.1_E.coli_0104_H4_illumina_1.0_genomic.fna 
 
  ln -s ~/workdir/sample_data/ ln -s ~/workdir/sample_data/SGCF-000221885.1_E.coli_0104_H4_illumina_1.0_genomic.fna .
```
 [ install Hisat](url)
#try to install hisat2 in linux
```
sudo apt_get install hisat2
 sudo apt install hisat2

```
 ```
 mkdir -p ~/workdir/hisat_align/hisatIndex && cd ~/workdir/hisat_align/hisatIndex
   hisat2_extract_splice_sites.py ~/workdir/sample_data/GCF_000221885.1_E.coli_0104_H4_Illumina_1.0/GCF_000221885.1_E.coli_0104_H4_Illumina_1.0_genomic.gtf > splicesites.tsv
    hisat2_extract_splice_sites.py 
    hisat2_extract_splice
  hisat2 --help
 ##(Troubleshooting)
```

[
 install BWA on a Linux](url)
  bunzip2 bwa-0.5.9.tar.bz2
  bunzip2 bwa-0.7.17.tar.bz2 
   tar xvf bwa-0.5.9.tar
  tar xvf bwa-0.7.17.tar
  cd bwa-0.5.9
  cd bwa-0.7.17
    make
 export PATH=$PATH:/Downloads/to/bwa-0.7.17 
   source ~/.bashrc (troubleshooting)

[*install bwa](url)
 ```
 sudo apt-get install bwa
```
[**index my genome**](url)
```
 mkdir -p ~/workdir/bwa_align/bwaIndex 
 cd ~/workdir/bwa_align/bwaIndex

  ln -s ~/workdir/sample_data/SGCF-000221885.1_E.coli_0104_H4_illumina_1.0_genomic.fna .
   bwa index -a bwtsw GCF-000221885.1_E.coli_0104_H4_illumina_1.0_genomic.fna
 
```
[[**[sequence](url) alignment](url)**](url)
 cd ~/workdir/bwa_align
  R1="$/home/crizma/Desktop/ngs/project e.coli/NEW R1 & R2/fastq/SRR292678_pass_1.fastq.gz"
  R2="$/home/crizma/Desktop/ngs/project e.coli/NEW R1 & R2/fastq/SRR292678_pass_2.fastq.gz"
    /usr/bin/time -v bwa mem bwaIndex/GCF_000221885.1_E.coli_0104_H4_illumina_1.0_genomic.fna $R1 $R2 > /SRR292678.sam
   bwa mem bwaIndex/GCF_000221885.1_E.coli_0104_H4_illumina_1.0_genomic.fna $R1 $R2 > /SRR292678.sam
  bwa mem bwaIndex/GCF_000221885.1_E.coli_0104_H4_illumina_1.0_genomic.fna $R1 $R2 > SRR292678.sam
    ls
-the sam file was empty (troubleshooting)

  
  ```
 cd ~/workdir/bwa_align
  R1="$HOME/workdir/sample_data/SRR292678_pass_1.fastq.gz" 
  cat $R1
  ln -s ~/workdir/sample_data/GCF_000221885.1_E.coli_0104_H4_Illumina_1.0_genomic.fa .
  bwa index -a bwtsw GCF_000221885.1_E.coli_0104_H4_Illumina_1.0_genomic.fa
    cd ~/workdir/bwa_align
  R1="$HOME/workdir/sample_data/SRR292678_pass_1.fastq.gz" 
 R2="$HOME/workdir/sample_data/SRR292678_pass_2.fastq.gz" 

```
  bwa mem bwaIndex/GCF_000221885.1_E.coli_0104_H4_Illumina_1.0_genomic.fa $R1 $R2 > SRR292678.sam
  /usr/bin/time -v bwa mem bwaIndex/GCF_000221885.1_E.coli_0104_H4_Illumina_1.0_genomic.fa $R1 $R2 > SRR678.SAM

  head -n10 SRR678.SAM 
  head -n100 SRR678.SAM 
   head -n200 SRR678.SAM 
   
conda install -c conda-forge -y biopython
    conda install -y samtools
[   # install Samtools](url)
  ` sudo apt install samtools`
 # convert SAM file to BAM
 `samtools view -S -b SRR678.SAM -o SRR678.bam
 
```

[# Sorting the BAM file](url)
```
samtools sort SRR678.bam -o sorted_SRR678.bam
  head -n10 sorted_SRR678.bam
```
```
```

[# Indexing the BAM file](url)
```
`samtools index sorted_SRR678.bam`
 

```
in this point i stoped this project as i join to another group


project troubleshooting
  
 project1
TCRLAMC PCR Sezary Limiting dilution PB cDNA 1000ng Beta chain
Prepare the [data**](url)
Download the reference genome from Ensembl and use the human GRCh38 version of the genome.
`wget http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz`

[retrival of data ](url)
```
sudo apt-get install sra-toolkit
perfetch SRR5809585
fastq dump  SRR5809585
head -n5000000 SRR5809585.fastq  > SRR585.fastq
 fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR585.fastq 
```

**install fastqc**
 ```
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
  gunzip fastqc_v0.11.9.zip 
   chmod a+x *
  ./fastqc
  sudo apt-get install openjdk-8-jre-headless
    ./fastqc
   sudo apt-get install openjdk-8-jre-headless
   ./fastqc
    sudo apt-get install openjdk-8-jre-headless
    ./fastqc
   java -version
  sudo apt-get install openjdk-8-jre-headless
   ./fastqc
```
 using the fastqc tools this data is bad.

**project2- genome sequence of Rosa chinensis to elucidate ornamental [traits](**url**)**
 
```
perfetch SRR6449591
fast dump SRR6449591
 head -n5000000 SRR6449591.fastq > SRR91.fastq
 fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR585.fastq 

```
visualize the quality of data by fastqc
the data is  bad also.(trouble shooting)

[project 3](**url**)

```
Download the reference genome from Ensembl and use the human GRCh38 version of the genome.
prepare the reference data`
wget http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz`
[prefetch SRR2079547
 fastq-dump SRR2079547.sra 
head -n SRR2079547.fastq 
 head -n 50 SRR2079547.fastq 
  wc -l SRR2079547.fastq 
fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR2079547.fastq 
```
the data is too large
` 





 

  
 
 

