#1. Access to publihed RNA-seq data

#Create environment variables
echo 'PATH=$PATH:~/.aspera/connect/bin/' >> ~/.bashrc
source ~/.bashrc

#Copy the key file under the file to the main directory
cp ~/.aspera/connect/etc/asperaweb_id_dsa.openssh ~/ 

#Copy the certificate to the bin directory
cp ~/.aspera/connect/etc/aspera-license /usr/local/bin/ 

#Creat the data directories 
mkdir data
cd data
mkdir PB1 PB2 ICV Gland

#Download published data from EBI, using ascp
#(1) peripheral blood 1
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/018/SRR10113618/SRR10113618.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/019/SRR10113619/SRR10113619.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/020/SRR10113620/SRR10113620.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/021/SRR10113621/SRR10113621.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/022/SRR10113622/SRR10113622.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/023/SRR10113623/SRR10113623.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/024/SRR10113624/SRR10113624.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/025/SRR10113625/SRR10113625.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/026/SRR10113626/SRR10113626.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/027/SRR10113627/SRR10113627.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/028/SRR10113628/SRR10113628.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/029/SRR10113629/SRR10113629.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/030/SRR10113630/SRR10113630.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/031/SRR10113631/SRR10113631.fastq.gz

cd ~/data/PB1 

for i in {18..31}
do
a1='ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/0'$i
a2='/SRR101136'$i
a3='.fastq.gz .'
echo $a1$a2$a2$a3
done

#(2) Ileocecal valve
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/032/SRR10113632/SRR10113632.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/033/SRR10113633/SRR10113633.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/034/SRR10113634/SRR10113634.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/035/SRR10113635/SRR10113635.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/036/SRR10113636/SRR10113636.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/037/SRR10113637/SRR10113637.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/038/SRR10113638/SRR10113638.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/039/SRR10113639/SRR10113639.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/040/SRR10113640/SRR10113640.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/041/SRR10113641/SRR10113641.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/042/SRR10113642/SRR10113642.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/043/SRR10113643/SRR10113643.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/044/SRR10113644/SRR10113644.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/045/SRR10113645/SRR10113645.fastq.gz

cd ~/data/ICV

for i in {32..45}
do
a1='ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR101/0'$i
a2='/SRR101136'$i
a3='.fastq.gz .'
echo $a1$a2$a2$a3
done

#(3) peripheral blood 2
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR116/061/SRR11624561/SRR11624561_1.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR116/061/SRR11624561/SRR11624561_2.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR116/062/SRR11624562/SRR11624562_1.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR116/062/SRR11624562/SRR11624562_2.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR116/063/SRR11624563/SRR11624563_1.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR116/063/SRR11624563/SRR11624563_2.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR116/064/SRR11624564/SRR11624564_1.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR116/064/SRR11624564/SRR11624564_2.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR116/065/SRR11624565/SRR11624565_1.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR116/065/SRR11624565/SRR11624565_2.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR116/066/SRR11624566/SRR11624566_1.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR116/066/SRR11624566/SRR11624566_2.fastq.gz

cd ~/data/PB2
for i in {1..6}
do
a1='ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR116/06'$i
a2='/SRR1162456'$i
a3='_1.fastq.gz .'
a4='_2.fastq.gz .'
echo $a1$a2$a2$a3
echo $a1$a2$a2$a4
done

#(4) salivary gland
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR841/007/SRR8418397/SRR8418397_1.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR841/007/SRR8418397/SRR8418397_2.fastq.gz
                                                 ................................
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR841/005/SRR8418435/SRR8418435_1.fastq.gz
#era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR841/005/SRR8418435/SRR8418435_2.fastq.gz

cd ~/data/Gland
for i in {391..435}
do
a0='ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR841/00'$i
a1=$(($i % 10))
a2='/SRR8418'$i
a3='_1.fastq.gz .'
a4='_2.fastq.gz .'
echo $a0$a1$a2$a2$a3
echo $a0$a1$a2$a2$a4
done
