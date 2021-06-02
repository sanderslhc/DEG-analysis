cd ~/data/Jejunum
mkdir fastQC
for i in $(ls *_1.fq.gz)
do
i=${i/_1.fq.gz/}
java -jar /home/WORKSPACE/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
-phred33 ${i}_1.fq.gz ${i}_2.fq.gz \
${i}_1.clean.fastq ${i}_1.unpaired.fastq ${i}_2.clean.fastq ${i}_2.unpaired.fastq \
ILLUMINACLIP:/home/WORKSPACE/Trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
fastqc -o ~/data/Jejunum/fastQC/ -t 6 {i}_1.clean.fastq
fastqc -o ~/data/Jejunum/fastQC/ -t 6 {i}_2.clean.fastq
done

cd ~/data/PB1
mkdir fastQC
for i in $(ls *.fastq.gz)
do
i=${i/.fastq.gz/}
java -jar /home/WORKSPACE/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
-phred33 ${i}.fastq.gz \
${i}.clean.fastq \
ILLUMINACLIP:/home/WORKSPACE/Trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
fastqc -o ~/data/PB1/fastQC/ -t 6 {i}.clean.fastq
done

cd ~/data/ICV
mkdir fastQC
for i in $(ls *.fastq.gz)
do
i=${i/.fastq.gz/}
java -jar /home/WORKSPACE/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
-phred33 ${i}.fastq.gz \
${i}.clean.fastq \
ILLUMINACLIP:/home/WORKSPACE/Trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
fastqc -o ~/data/ICV/fastQC/ -t 6 {i}.clean.fastq
done

cd ~/data/PB2
mkdir fastQC
for i in $(ls *_1.fastq.gz)
do
i=${i/_1.fastq.gz/}
java -jar /home/WORKSPACE/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
-phred33 ${i}_1.fastq.gz ${i}_2.fastq.gz \
${i}_1.clean.fastq ${i}_1.unpaired.fastq ${i}_2.clean.fastq ${i}_2.unpaired.fastq \
ILLUMINACLIP:/home/WORKSPACE/Trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
fastqc -o ~/data/PB2/fastQC/ -t 6 {i}_1.clean.fastq
fastqc -o ~/data/PB2/fastQC/ -t 6 {i}_2.clean.fastq
done

cd ~/data/Gland
mkdir fastQC
for i in $(ls *_1.fastq.gz)
do
i=${i/_1.fastq.gz/}
java -jar /home/WORKSPACE/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
-phred33 ${i}_1.fastq.gz ${i}_2.fastq.gz \
${i}_1.clean.fastq ${i}_1.unpaired.fastq ${i}_2.clean.fastq ${i}_2.unpaired.fastq \
ILLUMINACLIP:/home/WORKSPACE/Trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
fastqc -o ~/data/Gland/fastQC/ -t 6 {i}_1.clean.fastq
fastqc -o ~/data/Gland/fastQC/ -t 6 {i}_2.clean.fastq
done
