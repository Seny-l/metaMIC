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
- download [training models](https://zenodo.org/record/4717667#.YIQvu5MzZTY)

```
cd metaMIC/model
sh install_model.sh
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
- 
For metagenomics

```
metaMIC/metaMIC.py --bam $bam_file -c $contig_file -o $output_dir --pileup $pileup_file -m meta 
```
For isolate genomes

```
metaMIC/metaMIC.py --bam $bam_file -c $contig_file -o $output_dir --pileup $pileup_file -m single 
```


For more details about the usage of metaMIC, [read the docs](http:)

## Output
The output folder will contain
1. Misassembly score for each contig: **metaMIC_contig_score.txt** (for only metagenomics)
2. predicted misassembly breakpoints for misassembled contigs: **misassembly_breakpoint.txt**
3. Anomaly scores for each position in contigs: **anomaly_score.txt**
4. Fasta file of corrected contigs: **metaMIC_corrected_contigs.fa**
5. Some intermediate files

For more details about the output, [read the docs](http:)

## Complete option list
metaMIC:

```
Usage: metaMIC.py [options]

Options:
  -h, --help            show this help message and exit
  -1 READ1, --r1=READ1  paired-end #1 fasta/q files
  -2 READ2, --r2=READ2  paired-end #2 fasta/q files
  -p READ, --r=READ     smart pairing (ignoring #2 fasta/q)
  -c ASSEMBLIES, --contig=ASSEMBLIES
                        fasta file of assembled contigs
  --bam=BAMFILE         index bam file for alignment
  -a ASSEMBLER, --assembler=ASSEMBLER
                        The assembler-specific model or user-trained model
                        used for assembled fasta file [MEGAHIT/IDBA_UD/[new
                        training model specified by users]]
  -o OUTPUT, --output=OUTPUT
                        output directory for AIM results
  -m MODE, --mode=MODE  Applied to single genomic/metagenomic assemblies
                        [meta/single]
  -t THREADS, --threads=THREADS
                        Maximum number of threads [default: 8]
  -l MIN_LENGTH, --mlen=MIN_LENGTH
                        Minimum contig length [default: 5000bp]
  -s SPLIT_LENGTH, --slen=SPLIT_LENGTH
                        Minimum length of splitted fragments [default: 1000bp]
  --pileup=PILEUP       path to pileup file [samtools mpileup]
  --samtools=SAMTOOLS   path to samtools
  --bwa=BWA             path to bwa
  --jellyfish=JELLYFISH
                        path to jellyfish
  --train               Training on user-specific datasets
  --label=LABEL         Misassembly label of contigs for training assemblies
  --no-breakpoints      Do not locate possible breakpoints
  --no-correct          Do not break misassembled contigs at breakpoints
  --nb=BREAK_COUNT      Minimum number of read breakpoint counts for
                        correcting misassemblies in metagenomics
  --rb=BREAK_RATIO      Minimum read breakpoint ratio for correcting
                        misassemblies in metagenomics
  --at=ANOMALY_THRED    Minimum anomaly score for correcting misassemblies in
                        metagenomics
```
