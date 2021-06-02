# Get STAR source
cd /home/Sundongxiao/WORKSPACE/lhc
wget https://github.com/alexdobin/STAR/archive/2.7.9a.tar.gz
tar -xzf 2.7.9a.tar.gz
cd STAR-2.7.9a

#Compile
cd STAR/source
make STAR

#Create environment variables
echo 'PATH=$PATH:/home/Sundongxiao/WORKSPACE/lhc/STAR-2.7.9a/bin/Linux_x86_64/' >> ~/.STAR
source ~/.STAR
STAR --help

#download genome FASTA file and GTF file of cattle
cd /home/Sundongxiao/WORKSPACE/lhc/
mkdir genome
cd genome
wget http://ftp.ensembl.org/pub/release-104/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz
gunzip Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz
sudo mv Bos_taurus.ARS-UCD1.2.dna.toplevel.fa Bos_taurus.ARS-UCD1.2.fa
wget http://ftp.ensembl.org/pub/release-104/gtf/bos_taurus/Bos_taurus.ARS-UCD1.2.104.gtf.gz
gunzip Bos_taurus.ARS-UCD1.2.104.gtf.gz

#Create index with cattle genome
cd /home/Sundongxiao/WORKSPACE/lhc/STAR-2.7.9a
mkdir index-cow
STAR \
--runMode genomeGenerate \
--runThreadN 10 \
--genomeDir /home/Sundongxiao/WORKSPACE/lhc/STAR-2.7.9a/index-cow \ #the dictionary of index
--genomeFastaFiles  /home/Sundongxiao/WORKSPACE/lhc/genome/Bos_taurus.ARS-UCD1.2.fa \ #input fasta file
--sjdbGTFfile /home/Sundongxiao/WORKSPACE/lhc/genome/Bos_taurus.ARS-UCD1.2.104.gtf \ #input gtf file
--sjdbOverhang 149 #length of reads - 1

#Align RNA-seq reads
#Jejunum
cd ~/data/Jejunum
for i in $(ls *_1.fq.gz)
do
i=${i/_1.fq.gz/}
STAR \
--runThreadN 10 \
--genomeDir /home/Sundongxiao/WORKSPACE/lhc/STAR-2.7.9a/index-cow \
--readFilesCommand zcat \ #if fastq files are conpressed, you need to add this
--outSAMtype BAM SortedByCoordinate \
--sjdbOverhang 149 \
--readFilesIn ${i}_1.fq.gz ${i}_2.fq.gz \
--outFileNamePrefix /home/Sundongxiao/WORKSPACE/lhc/STAR-2.7.9a/data/Jejunum/${i}.sorted \
--outBAMsortingThreadN 5
done

#PB1
cd ~/data/PB1
for i in $(ls *.fastq.gz)
do
i=${i/.fastq.gz/}
STAR \
--runThreadN 10 \
--genomeDir /home/Sundongxiao/WORKSPACE/lhc/STAR-2.7.9a/index-cow \
--readFilesCommand zcat \ 
--outSAMtype BAM SortedByCoordinate \
--sjdbOverhang 149 \
--readFilesIn ${i}.fastq.gz \
--outFileNamePrefix /home/Sundongxiao/WORKSPACE/lhc/STAR-2.7.9a/data/PB1/${i}.sorted \
--outBAMsortingThreadN 5
done

#ICV
cd ~/data/ICV
for i in $(ls *.fastq.gz)
do
i=${i/.fastq.gz/}
STAR \
--runThreadN 10 \
--genomeDir /home/Sundongxiao/WORKSPACE/lhc/STAR-2.7.9a/index-cow \
--readFilesCommand zcat \ 
--outSAMtype BAM SortedByCoordinate \
--sjdbOverhang 149 \
--readFilesIn ${i}.fastq.gz \
--outFileNamePrefix /home/Sundongxiao/WORKSPACE/lhc/STAR-2.7.9a/data/ICV/${i}.sorted \
--outBAMsortingThreadN 5
done

#PB2
cd ~/data/PB2
for i in $(ls *_1.fastq.gz)
do
i=${i/_1.fastq.gz/}
STAR \
--runThreadN 10 \
--genomeDir /home/Sundongxiao/WORKSPACE/lhc/STAR-2.7.9a/index-cow \
--readFilesCommand zcat \ 
--outSAMtype BAM SortedByCoordinate \
--sjdbOverhang 149 \
--readFilesIn ${i}_1.fq.gz ${i}_2.fq.gz \
--outFileNamePrefix /home/Sundongxiao/WORKSPACE/lhc/STAR-2.7.9a/data/PB2/${i}.sorted \
--outBAMsortingThreadN 5
done

#Gland
cd ~/data/Gland
for i in $(ls *_1.fastq.gz)
do
i=${i/_1.fastq.gz/}
STAR \
--runThreadN 10 \
--genomeDir /home/Sundongxiao/WORKSPACE/lhc/STAR-2.7.9a/index-cow \
--readFilesCommand zcat \ 
--outSAMtype BAM SortedByCoordinate \
--sjdbOverhang 149 \
--readFilesIn ${i}_1.fq.gz ${i}_2.fq.gz \
--outFileNamePrefix /home/Sundongxiao/WORKSPACE/lhc/STAR-2.7.9a/data/Gland/${i}.sorted \
--outBAMsortingThreadN 5
done

#get count value
#Jejunum
featurecount \
-T 5 \
-t exon \
-g gene_id \
-a /home/Sundongxiao/WORKSPACE/lhc/genome/Bos_taurus.ARS-UCD1.2.103.gtf \
-o counts_NP.txt \
NP1.sorted.bam NP2.sorted.bam NP3.sorted.bam

featurecount \
-T 5 \
-t exon \
-g gene_id \
-a /home/Sundongxiao/WORKSPACE/lhc/genome/Bos_taurus.ARS-UCD1.2.103.gtf \
-o counts_DP.txt \
DP1.sorted.bam DP2.sorted.bam DP3.sorted.bam
