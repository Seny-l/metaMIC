#!/bin/bash

contig=$1
bamfile=$2
output=$3
Samtools=$4

mkdir -p ${output}/temp/split
mkdir -p ${output}/temp/split/contigs
mkdir -p ${output}/temp/split/reads
mkdir -p ${output}/temp/split/sams
awk '/^>/{s=++num}{print > file"/split_"s".fa"}' file="${output}/temp/split/contigs" $contig


for file in `ls ${output}/temp/split/contigs`
do
    cat ${output}/temp/split/contigs/${file} | grep '>' | sed 's/>//g' | awk 'BEGIN{FS=" "}{print $1}'
done  > ${output}/temp/split/contig_name.txt


for file in `ls ${output}/temp/split/contigs`
do
echo ${file%.fa}
done > ${output}/temp/split/split_file_name.txt



for file in `cat ${output}/temp/split/contig_name.txt`
do
if [ ! -s "${output}/temp/split/reads/${file}.read.fa" ];
then
${Samtools} view -t 4  ${bamfile} ${file} > ${output}/temp/split/sams/${file}.sam
cat ${output}/temp/split/sams/${file}.sam | grep "=" | awk '{print ">"$1"\n"$10}' > ${output}/temp/split/reads/${file}.read.fa
fi
done
