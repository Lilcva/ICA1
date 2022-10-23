#!/usr/bin/env bash
#--------------------------direction----------------------
mkdir ~/ICA1
mkdir ~/ICA1/fastq
mkdir ~/ICA1/qc
mkdir ~/ICA1/index/
mkdir ~/ICA1/bam
mkdir ~/ICA1/bed/
mkdir ~/ICA1/count/
mkdir ~/ICA1/result/
mkdir ~/ICA1/group/
mkdir ~/ICA1/analysis

#------------------------address variables------------------
fastq_addresss=$HOME"/ICA1/fastq/"
sn_address=$fastq_addresss"samplename"
qc_address=$HOME"/ICA1/qc/" 
index_address=$HOME"/ICA1/index/"
genome="TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz"
bam_address=$HOME"/ICA1/bam/"
bed_address=$HOME"/ICA1/bed/"
bed_ref="/localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed"
count_address=$HOME"/ICA1/count/"
result_address=$HOME"/ICA1/result/"
group_address=$HOME"/ICA1/group/"
analysis_result=$HOME"/ICA1/analysis/"

touch $sn_address
touch $HOME"/ICA1/allQcReport"

#------------------------QC------------------

cp -r /localdisk/data/BPSM/ICA1/fastq/* $fastq_addresss
echo "Copy finished!"
# Read and store the sample name into a file called samplename
ls $fastq_addresss | grep ".gz$" | awk -F "_" '!a[$1]++{print $1}' > $sn_address
# qc
for i in {1..48}
do 
name=$(cat $sn_address | sed -n $i'p')
gunzip $fastq_addresss$name"_1.fq.gz"
gunzip $fastq_addresss$name"_2.fq.gz"
fastqc $fastq_addresss$name"_1.fq" -q -o $qc_address
fastqc $fastq_addresss$name"_2.fq" -q -o $qc_address
done
echo "QC finished!"
# write a overall QC report
echo -e "Basic Statistics\tPer base sequence quality\tPer sequence quality scores\tPer base sequence content\tPer base N content \n" > $HOME"/ICA1/allQcReport"
for i in {1..48}
do 
name=$(cat $sn_address | sed -n $i'p')
for j in {1..2}
do 
unzip -q $qc_address$name"_"$j"_fastqc.zip" -d $qc_address
summary=$qc_address$name"_"$j"_fastqc/summary.txt"
basic=$(cat $summary | sed -n 1p | cut -f 1)
baseQuality=$(cat $summary | sed -n 2p | cut -f 1)
sequenceScore=$(cat $summary | sed -n 3p | cut -f 1)
baseContent=$(cat $summary | sed -n 4p | cut -f 1)
nContent=$(cat $summary | sed -n 6p | cut -f 1)
echo -e $basic"\t"$baseQuality$"\t"$sequenceScore"\t"$baseContent"\t"$nContent"\n" >> $HOME"/ICA1/allQcReport"
done
done
echo "QC report generated!"

#------------------------Alignment------------------

#build the index of Trypanosoma congolense
cp "/localdisk/data/BPSM/ICA1/Tcongo_genome/"$genome $index_address
gunzip $index_address$genome
# build index
bowtie2-build -q $index_address"TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta" $index_address"TriTryp"
#alignment
for i in {1..48}
do 
name=$(cat $sn_address | sed -n $i'p')
name1=$name"_1.fq"
name2=$name"_2.fq"
bam_name=$bam_address$name".output.bam"
bowtie2 -q -p 10 -x $index_address"TriTryp" -1 $fastq_addresss$name1 -2 $fastq_addresss$name2 | samtools sort -O bam -@ 10 -o $bam_name
done
echo "alignment finished!"

#------------------------Gene count------------------

for i in {1..48}
do
name=$(cat $sn_address | sed -n $i'p')
bam_name=$bam_address$name".output.bam"
bedtools bamtobed -i $bam_name > $bed_address$name".bed"
bedtools coverage -a $bed_ref -b $bed_address$name".bed" -counts > $count_address$name".out"
done
echo "counts finished!"

#------------------------count mean------------------

clone=("Clone1" "Clone2" "WT")
time=("0" "24" "48")
treat=("Induced" "Uninduced")

# Group according to the treatment method, and store the sample name under each group
# totally, there are 15 groups: 3 X 3 X 2 -3 
# Since for Clone1, Clonee, WT with 0 hour, there is no Induced
for item in ${clone[@]}
do 
    for hour in ${time[@]}
    do
        for treatment in ${treat[@]}
        do
        all_sample=($(cat $fastq_addresss"Tco.fqfiles" | awk '{FS="\t";}{if($2=="'$item'"){print $0;}}' | awk 'BEGIN{FS="\t";}{if($4=="'$hour'"){print $0;}}' | awk 'BEGIN{FS="\t";}{if($5=="'$treatment'"){print $0;}}' | cut -f 1))
        echo ${all_sample[@]}>>$group_address$item"_"$hour"_"$treatment
        done
    done
done
# these three group do not exist 
rm -f $group_address"Clone1_0_Induced"
rm -f $group_address"WT_0_Induced"
rm -f $group_address"Clone2_0_Induced"
# This array is used to create a 0 array with the gene count length
ref_arr=($(cat $bed_ref | cut -f 1))
all_group_file=($(ls $group_address))

# calculate mean
for file in ${all_group_file[@]}
do 
all_sample=($(cat $group_address$file))
arr_len=${#all_sample[@]}
# Create an array with the same length and number of genes as 0, and add samples' count to this array
sum_arr=("${ref_arr[@]/*/0}")
    for sample in ${all_sample[@]}
    do
    sample=$(echo $sample | awk '{print substr($1, 4, 4)}')
    sample="Tco-"$sample
    sample_count=$count_address$sample".out"
    count_arr=($(cat $sample_count | cut -f 6))
        for i in "${!sum_arr[@]}"
        do
        sum_arr[i]=$[`expr ${sum_arr[i]}+${count_arr[i]}`]
        done
    done
    for i in "${!sum_arr[@]}"
    do
    mean=$[`expr ${sum_arr[i]}/$arr_len`]
    j=$[`expr $i+1`]
    gene_name=$(cat $bed_ref | sed -n $j'p' | cut -f 4)
    gene_exp=$(cat $bed_ref | sed -n $j'p' | cut -f 5)
    echo -e $mean"\t"$gene_name"\t"$gene_exp >>$result_address$file"_mean.out"
    done
done
echo "Group mean calculated!"

#------------------------analysis------------------

#! WT_0 vs Clone1_0 vs Clone2_0
#First column is WT, second is Clone1
length=$(cat $result_address"WT_0_Uninduced_mean.out" | wc -l)
for ((i=1;i<=$length;i++));
do
wt_count=$(cat $result_address"WT_0_Uninduced_mean.out" | sed -n $i'p' | cut -f 1)
clo1_count=$(cat $result_address"Clone1_0_Uninduced_mean.out" | sed -n $i'p' | cut -f 1)
gene_name=$(cat $bed_ref | sed -n $i'p' | cut -f 4)
gene_exp=$(cat $bed_ref | sed -n $i'p' | cut -f 5)
if [ $wt_count -ne 0 -a $clo1_count -ne 0 ]; then
    fold=$(printf "%.2f" `echo "scale=2;$clo1_count/$wt_count" | bc`)
    echo -e $wt_count"\t"$clo1_count"\t"$fold"\t"$gene_name"\t"$gene_exp>>$analysis_result"WT_vs_Clone1_0.out"
fi
done
#First column is WT, second is Clone2
for ((i=1;i<=$length;i++));
do
wt_count=$(cat $result_address"WT_0_Uninduced_mean.out" | sed -n $i'p' | cut -f 1)
clo2_count=$(cat $result_address"Clone2_0_Uninduced_mean.out" | sed -n $i'p' | cut -f 1)
gene_name=$(cat $bed_ref | sed -n $i'p' | cut -f 4)
gene_exp=$(cat $bed_ref | sed -n $i'p' | cut -f 5)
if [ $wt_count -ne 0 -a $clo2_count -ne 0 ]; then
    fold=$(printf "%.2f" `echo "scale=2;$clo2_count/$wt_count" | bc`)
    echo -e $wt_count"\t"$clo2_count"\t"$fold"\t"$gene_name"\t"$gene_exp>>$analysis_result"WT_vs_Clone2_0.out"
fi
done
#! sort
cat $analysis_result"WT_vs_Clone1_0.out" | sort -n -k 3r >  $analysis_result"WT_vs_Clone1_0_sorted.out"
cat $analysis_result"WT_vs_Clone2_0.out" | sort -n -k 3r >  $analysis_result"WT_vs_Clone2_0_sorted.out"

#! WT_24_in vs WT_24_un
#First column is induced, second is uninduced
for ((i=1;i<=$length;i++));
do
in_count=$(cat $result_address"WT_24_Induced_mean.out" | sed -n $i'p' | cut -f 1)
un_count=$(cat $result_address"WT_24_Uninduced_mean.out" | sed -n $i'p' | cut -f 1)
gene_name=$(cat $bed_ref | sed -n $i'p' | cut -f 4)
gene_exp=$(cat $bed_ref | sed -n $i'p' | cut -f 5)
if [ $in_count -ne 0 -a $un_count -ne 0 ]; then
    fold=$(printf "%.2f" `echo "scale=2;$in_count/$un_count" | bc`)
    echo -e $in_count"\t"$un_count"\t"$fold"\t"$gene_name"\t"$gene_exp>>$analysis_result"WT_24_induced_vs_WT_24_Uninduced.out"
fi
done
cat $analysis_result"WT_24_induced_vs_WT_24_Uninduced.out" | sort -n -k 3r >  $analysis_result"WT_24_induced_vs_WT_24_Uninduced_sorted.out"

#! WT_48_in vs WT_48_un
#First column is induced, second is uninduced
for ((i=1;i<=$length;i++));
do
in_count=$(cat $result_address"WT_48_Induced_mean.out" | sed -n $i'p' | cut -f 1)
un_count=$(cat $result_address"WT_48_Uninduced_mean.out" | sed -n $i'p' | cut -f 1)
gene_name=$(cat $bed_ref | sed -n $i'p' | cut -f 4)
gene_exp=$(cat $bed_ref | sed -n $i'p' | cut -f 5)
if [ $in_count -ne 0 -a $un_count -ne 0 ]; then
    fold=$(printf "%.2f" `echo "scale=2;$in_count/$un_count" | bc`)
    echo -e $in_count"\t"$un_count"\t"$fold"\t"$gene_name"\t"$gene_exp>>$analysis_result"WT_48_induced_vs_WT_48_Uninduced.out"
fi
done
cat $analysis_result"WT_48_induced_vs_WT_48_Uninduced.out" | sort -n -k 3r >  $analysis_result"WT_48_induced_vs_WT_48_Uninduced_sorted.out"
echo "All process finished!"












