# metaMIC: Reference-free Misassembly Identification and Correction of metagenomic assemblies
metaMIC is a fully automated tool for identifying and correcting misassemblies of (meta)genomic assemblies with the following three steps. Firstly, metaMIC extracts various types of features from the alignment between paired-end sequencing reads and the assembled contigs.Secondly, the features extracted in the first step will be used as input of a random forest classifier for identifying misassembled metagenomic assemblies. Thirdly, metaMIC will localize misassembly breakpoints for each misassembled contig and then corrects misassemblies by splitting into parts at the breakpoints.


## Requirements and Installation
Make sure you have the dependencies below installed and accessible in your $PATH.

### Prepare dependencies

- [python 3.6-3.9](https://www.python.org/downloads/)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)
- [pysam](https://pypi.org/project/pysam/0.8.4/)
- [bwa 0.7.17](https://sourceforge.net/projects/bio-bwa/files/)
- [samtools 1.9](https://sourceforge.net/projects/samtools/files/samtools/)
- [jellyfish](http://www.cbcb.umd.edu/software/jellyfish/)

1. install python modules: pandas, numpy, pysam
```
conda install -c bioconda pandas numpy pysam 
```
or

```
pip install pysam pandas numpy
```

2. download and install samtools

```
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -jxvf samtools-1.9.tar.bz2
cd samtools-1.9
./configure
make
make install
export PATH=`pwd`:$PATH
```
3. download and install bwa

```
wget https://sourceforge.net/projects/bio-bwa/files/latest/download/bwa-0.7.17.tar.bz2
tar -jxvf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
export PATH=`pwd`:$PATH
```
4. download and install jellyfish


```
wget http://www.cbcb.umd.edu/software/jellyfish/jellyfish-1.1.10.tar.gz
tar zxvf jellyfish-1.1.10.tar.gz
cd jellyfish-1.1.10
./configure
make
make install
export PATH=`pwd`/bin:$PATH
```

### Installation

##### Install metaMIC via git

```
git clone https://github.com/Seny-l/metaMIC.git
```
download training models

```
de metaMIC
wget xx

```

## Quick Start
- First map paired-end reads to assembled contigs

```
bwa index $contig_file
bwa mem -a -t 8 $contig_file $read1 $read2 | samtools view -h -q 10 -m 50 -F 4 -b | samtools sort > $bam_file
```
- generate pileup file

```
samtools mpileup -C 50 -A -f $contig_file $bam_file |  awk '$3 != "N"' > $pileup_file
```
- run metaMIC

```
metaMIC/metaMIC.py --bam $bam_file -c $contig_file -o $output_dir --pileup $pileup_file -m meta 
```


