# metaMIC
metaMIC is a fully automated tool for identifying and correcting misassemblies of (meta)genomic assemblies with the following three steps. Firstly, metaMIC extracts various types of features from the alignment between paired-end sequencing reads and the assembled contigs.Secondly, the features extracted in the first step will be used as input of a random forest classifier for identifying misassembled metagenomic assemblies. Thirdly, metaMIC will localize misassembly breakpoints for each misassembled contig and then corrects misassemblies by splitting into parts at the breakpoints.

![图片 1](https://user-images.githubusercontent.com/81548043/115959298-20b1a300-a53e-11eb-9bee-cd81265286ee.png)
