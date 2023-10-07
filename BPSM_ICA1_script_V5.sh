#!/bin/bash


#0.
#setting the work directory of this project,and set it as a variable to make users change their work directory easily.
#enter the work directory and make "fastq_file".
read -p "please enter your work directory:" Work_Address
echo ${Work_Address}
mkdir ./${Work_Address}
cd ./${Work_Address}
mkdir ./fastq_file

#copy the source fastq files to "fastq_file".
cp /localdisk/data/BPSM/ICA1/fastq/*.gz ./fastq_file
cp /localdisk/data/BPSM/ICA1/fastq/*.fqfiles ./fastq_file


#1.
#In ${Work_Address}.
#do fastqc.
#make "results_fastqc",do fastqc and ouput the results to "results_fastqc" directory.
mkdir ./results_fastqc 
fastqc  ./fastq_file/*.gz -o ./results_fastqc


#2.1 
#In ${Work_Address}.
#assess the numbers and quality of the raw sequence data based on the output(s) of fastqc.
#make "results_fastqc_unzip",and output the unzip files to this "results_fastqc_unzip" directory.
mkdir ./results_fastqc_unzip
unzip  "./results_fastqc/*.zip" -d ./results_fastqc_unzip 

#2.2 
#In "results_fastqc_unzip".
#Select "Per base sequence quality" as the judgment standard of sequence quality.
#output the number of PASS, FAIL, WARN at the same time, and display the content of which sequence "FAIL", "WARN".

cd ./results_fastqc_unzip
echo "If in this sample， the Per base sequence quality is PASS, the quality will be considered as good, and it will show you how many PASS, FAIL and WARN respectively."
for i in $(ls -d *); 
do 
sample_name=${i:0:10}
echo "In ${sample_name}:"
grep -wc "PASS" ./${i}/summary.txt | xargs echo "The number of PASS is"
grep -wc "FAIL" ./${i}/summary.txt | xargs echo "The number of FAIL is"
grep -wc "WARN" ./${i}/summary.txt | xargs echo "The number of WARN is"
awk 'BEGIN{FS="\t";}{if($1 == "FAIL"){print $2, "is FAIL";}}' ./ ${i}/summary.txt
awk 'BEGIN{FS="\t";}{if($1 == "WARN"){print $2, "is WARN";}}' ./${i}/summary.txt
awk 'BEGIN{FS="\t";}{if($1 == "PASS"&& $2=="Per base sequence quality"){print "The sequence quality of this sample is good";}}' ./${i}/summary.txt
done


#3.1
#change location to ${Work_Address}
#gunzip the fastq.gz files 
cd ../
gunzip ./fastq_file/*.gz 

#3.2 
#In ${Work_Address}
#obtain the genome sequence and gunzip them.
cp /localdisk/data/BPSM/ICA1/Tcongo_genome/*gz ./
ls
gunzip ./TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz 
ls

#3.3
#In ${Work_Address}
#build index, and the index's name is "Tco_index".
bowtie2-build ./TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta Tco_index

#3.4
#copy the Tco_index into the same directory of fastq files.
#change location to fastq_file
#cut the fields1, fields6, fields7 to create the variable for the for loop of align.

cp ./*bt2 ./fastq_file
cd ./fastq_file

cat ./Tco.fqfiles | cut -f1 ./Tco.fqfiles > ./Sample_name
sed -i '1d' ./Sample_name
cat ./Tco.fqfiles | cut -f6 ./Tco.fqfiles > ./Sample_name1
sed -i '1d' ./Sample_name1
cat ./Tco.fqfiles | cut -f7 ./Tco.fqfiles > ./Sample_name2
sed -i '1d' ./Sample_name2

#3.5
#In fastq_file
#align the read pairs by bowtie2 and trans the outputs of bowtie2 into bam format.
 
for i in {1..48}; 
do 
 sample1=$(sed -n $i'p' ./Sample_name1)
 name1="${sample1:0:13}" 
 sample2=$(sed -n $i'p' ./Sample_name2)                    
 name2="${sample2:0:13}" 
 sample=$(sed -n $i'p' ./Sample_name)
 bam_name=$sample"output.bam" 
 bowtie2 -p 10 -x ./Tco_index -1 ./$name1 -2 ./$name2 | samtools sort -O bam -@ 10 -o - > ./$bam_name; 
done


#4.1
#In fastq_file
#copy the "bedfile" to fastq_file.
#generate counts data by using bedtools to.Cov files.

cp  /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed ./
for i in {1..48}
 do
 sample=$(sed -n $i'p' ./Sample_name)
 bam_name=$sample"output.bam"
 Cov_name=$sample".Cov"
 bedtools coverage -a ./TriTrypDB-46_TcongolenseIL3000_2019.bed -b ./$bam_name > ./$Cov_name
done


#5.1
#In fastq_file.
#Grouped according to different experimental conditions.
#output the results to "groups" directory.

mkdir ./groups
while read SampleName SampleType Replicate Time Treatment rest
 do
 condition=${SampleType}_${Time}_${Treatment}
 echo ${SampleName} >> ./groups/${condition}.txt
 done < ./Tco.fqfiles
 #delete the invalid file
rm ./groups/SampleType*
 		
#5.2
#change location to groups.
#generate statistical mean (average) to "group_mean" files.
#'a="${a} ${s}_overlap" 'means that you could get a input after the for loop done.
#then you could use 'paste' to paste multiply files.
 		
cd ./groups
for f in $(ls *.txt)
 do 
 a=""
   for s in $(cat ${f})
    do
    a="${a} ${s}_overlap"
    Cov_name=${s}".Cov"
    cut -f6 ../$Cov_name > ./${s}"_overlap" 
   done
 group=${f%.txt*}
 group_mean=${group}"_mean"
 paste ${a} > ./${group}"_number"
 awk '{sum=0;for(i=1;i<=NF;i++) sum+=$i;printf("%.2f\n", sum/NF)}' ./${group}"_number" > ./$group_mean
done
 		                               
 		                               
#5.3.
#In groups.
#generate plain text tab-delimited output files to one total file and 15 separated files.
 		                               
z=""
for i in $(ls *_mean)
 do
 z="${z} $i"
 name=${i%_mean*}
 paste ../TriTrypDB-46_TcongolenseIL3000_2019.bed ./$i  > ./$name"_file"
done
paste ../TriTrypDB-46_TcongolenseIL3000_2019.bed ./$z  > ./group_mean_file.txt
 		                                    
#change location to {Work_Address}.
#move these results files to ${Work_Address}.
cd ../../
mkdir ./group_file
mv ./fastq_file/groups/*_file* ./group_file 
 		                                    
 		                                 
#6.1
#compare the fold change.
#"read -p" is a user interactive command.
#this for loop allows the user to choose which two groups to compare and whether to continue the comparison.
cd ./group_file	                                    
echo 
read -p "Do you want to compare two groups?  y/n: " choice
while [ "$choice" == "y" ]
 do
 read -p "enter the filename1 you want to compare(you could see the file names in 'group_file', which is in your 'Work Address' directory):" file1
 read -p "enter the filename2 you want to compare(you could see the file names in 'group_file', which is in your 'Work Address' directory):" file2
 		                                  
 cut -f 6 ./$file1 > ./file1_mean
 cut -f 6 ./$file2 > ./file2_mean
 compare=${file1}"vs"${file2}
 		                                    
 paste ./file1_mean ./file2_mean > ./$compare"_mean"
 awk '{a=log($1+0.01)/log(2);b=log($2+0.01)/log(2);fold=b-a;printf("%.2f\n",fold)}' ./$compare"_mean" > ./$compare"_fold_change"
 paste ../fastq_file/TriTrypDB-46_TcongolenseIL3000_2019.bed ./$compare"_fold_change" | sort -r -t$'\t' -k6,6 > ./$compare"_results"
 read -p "Do you want to keep comparing?  y/n: " choice
done
echo "The comparison is done, you could see the results in your 'group_file'"
echo "Thank you for using"
 		                                    
