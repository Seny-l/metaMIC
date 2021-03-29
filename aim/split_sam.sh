contig=$1
output=$2
Samtools=$3

mkdir -p ${output}/temp/split
mkdir -p ${output}/temp/split/contigs
mkdir -p ${output}/temp/split/reads
mkdir -p ${output}/temp/split/sams
cd ${output}/temp/split/contigs
awk '/^>/{s=++num}{print > "split_"s".fa"}' $contig 

if [ ! -s "${output}/temp/split/contig_name.txt" ];
then
for file in `ls`
do
    cat ${file} | grep '>' | sed 's/>//g' | awk 'BEGIN{FS=" "}{print $1}'  
done  > ${output}/temp/split/contig_name.txt
fi

if [ ! -s "${output}/temp/split/split_file_name.txt" ];
then
for file in `ls`
do
echo ${file%.fa}
done > ${output}/temp/split/split_file_name.txt
fi


for file in `cat ${output}/temp/split/contig_name.txt`
do
if [ ! -s "${output}/temp/split/reads/${file}.read.fa" ];
then
${Samtools} view -t 4  ${output}/temp/sam/contigs.filter.sort.bam ${file} > ${output}/temp/split/sams/${file}.sam
cat ${output}/temp/split/sams/${file}.sam | grep "=" | awk '{print ">"$1"\n"$10}' > ${output}/temp/split/reads/${file}.read.fa
fi
done    
