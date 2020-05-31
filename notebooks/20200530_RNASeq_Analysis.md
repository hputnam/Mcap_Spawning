# Required Programs

fastqc 
multiqc  
fastp  
samtools  
hisat2 
stringtie 



# Main Folder Spawn

`mkdir Spawn`   
`cd Spawn`

nohup scp -r -P 2292 /Users/hputnam/Downloads/Raw_Data/T4* hputnam@kitt.uri.edu:/home/hputnam/Spawn

 
### Genome Links
http://cyanophora.rutgers.edu/montipora/Mcap.genome_assembly.fa.gz
http://cyanophora.rutgers.edu/montipora/Mcap.GFFannotation.gff


# QC raw files

`mkdir fastqc_raw`

`cd fastqc_raw`

`fastqc ../Raw_Data/*.fastq.gz  -o /home/hputnam/Spawn/fastqc_raw`

`multiqc .`

`mv multiqc_report.html raw_multiqc_report.html`

# Trimming

`cd ../`  
`mkdir cleaned_reads`  
`cd cleaned_reads`

"SRR5311171" "SRR5311255" "SRR5311304" "SRR5311305" "SRR5311306" "SRR5311307" "SRR5311308" "SRR5311310" "SRR5311311" "SRR5311312" "SRR5311313" "SRR5311342" "SRR5311343" "SRR5311344" "SRR5311345" "SRR5311346" "SRR5311347" "SRR5311348"

```
sh -c 'for file in "SRR5311171" "SRR5311255" "SRR5311304" "SRR5311305" "SRR5311306" "SRR5311307" "SRR5311308" "SRR5311310" "SRR5311311" "SRR5311312" "SRR5311313" "SRR5311342" "SRR5311343" "SRR5311344" "SRR5311345" "SRR5311346" "SRR5311347" "SRR5311348"
do
fastp \
--in1 ../Raw_Data/${file}_1.fastq.gz \
--in2 ../Raw_Data/${file}_2.fastq.gz \
--out1 ${file}_1_clean.fastq.gz \
--out2 ${file}_2_clean.fastq.gz \
--failed_out ${file}_failed.txt \
--qualified_quality_phred 20 \
--unqualified_percent_limit 10 \
--length_required 50 detect_adapter_for_pe \
--cut_right cut_right_window_size 5 cut_right_mean_quality 20
done'
```

# QC trimmed files
mkdir fastqc_cleaned

`fastqc ../cleaned_reads/*.fastq.gz  -o /home/hputnam/Spawn/fastqc_cleaned`

`multiqc fastqc_cleaned`

`mv multiqc_report.html cleaned_multiqc_report.html`


# Download and view files

scp -r -P 2292 hputnam@kitt.uri.edu:/home/hputnam/Spawn/fastqc_raw/raw_multiqc_report.html /Users/hputnam/MyProjects/Mcap_Spawning_Timing/QC

scp -r -P 2292 hputnam@kitt.uri.edu:/home/hputnam/Spawn/fastqc_cleaned/cleaned_multiqc_report.html /Users/hputnam/MyProjects/Mcap_Spawning_Timing/QC

*HISAT2 is a fast and sensitive alignment program for mapping next-generation DNA and RNA sequencing reads to a reference genome.*

- Index the reference genome
- Align clean reads to the reference genome

#### Index the reference genome

Index the reference genome in the reference directory.

++HISAT2-build Alignment Arguments Used++:  
- <reference_in> - name of reference files  
- <gt2_base> -  basename of index files to write  
- -f -  reference file is a FASTA file

```
cd /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref

hisat2-build Mcap.genome_assembly.fa Mcap_ref
```

# Alignment of clean reads to the reference genome

### Aligning paired end reads
mkdir mapping
cd mapping


"SRR5311255" "SRR5311304" "SRR5311305" "SRR5311306" "SRR5311307" "SRR5311308" "SRR5311310" "SRR5311311" "SRR5311312" "SRR5311313" "SRR5311342" "SRR5311343" "SRR5311344" "SRR5311345" "SRR5311346" "SRR5311347" "SRR5311348"
## Mcap Mapping
`
sh -c 'for i in "SRR5311171" "SRR5311255" "SRR5311304" "SRR5311305" "SRR5311306" "SRR5311307" "SRR5311308" "SRR5311310" "SRR5311311" "SRR5311312" "SRR5311313" "SRR5311342" "SRR5311343" "SRR5311344" "SRR5311345" "SRR5311346" "SRR5311347" "SRR5311348"
do
hisat2 -p 8 --dta -q -x /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Mcap_ref \
-1 ../cleaned_reads/${i}_1_clean.fastq.gz \
-2 ../cleaned_reads/${i}_2_clean.fastq.gz -S ${i}.sam 
done'
`

# Convert and sort Sam to BAM

`
sh -c 'for i in "SRR5311171" "SRR5311255" "SRR5311304" "SRR5311305" "SRR5311306" "SRR5311307" "SRR5311308" "SRR5311310" "SRR5311311" "SRR5311312" "SRR5311313" "SRR5311342" "SRR5311343" "SRR5311344" "SRR5311345" "SRR5311346" "SRR5311347" "SRR5311348"
do
samtools sort -@ 8 -o ${i}.bam ${i}.sam
done'
`

# Remove Sam
`
rm *.sam

`

# Assemble aligned reads and quantify transcripts 

---

*StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts.*

- Reference-guided assembly with novel transcript discovery
- Merge output GTF files and assess the assembly performance
- Merged-GTF guided assembly without novel transcript discovery
- Compilation of GTF-files into gene and transcript count matrices

For our assembly, we will have to run StringTie twice. The first run will be a reference guided assembly that will allow for discovery of novel transcripts (by leaving out the -e option). Then we will merge the output GTF files and examine the sensitivity of our assembly. We will use the merged-GTF from our first assembly to guide our second StringTie run (including the ```-e``` option). This second run is necessary in order to compile our GTF files into gene and transcript count matrices that we will need for our differential expression analysis, because the StringTie script that compiles the GTF files ```prepDE.py``` only runs if the ```-e``` option is "on" during the previous assembly.

#### Reference-guided assembly with novel transcript discovery


Create the StringTie reference-guided assembly script, ```McapStringTie-assembly-to-ref.sh``` *inside of the StringTie program directory.*  

++StringTie Arguments Used++:  
- -p - Specify number of processers
- -G - Specify annotation file
- -o - Name of output file

#Mcap 
`
sh -c 'for i in "SRR5311171" "SRR5311255" "SRR5311304" "SRR5311305" "SRR5311306" "SRR5311307" "SRR5311308" "SRR5311310" "SRR5311311" "SRR5311312" "SRR5311313" "SRR5311342" "SRR5311343" "SRR5311344" "SRR5311345" "SRR5311346" "SRR5311347" "SRR5311348"
do
stringtie ${i}.bam -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Mcap.GFFannotation.gff -o ${i}.gtf 
done'
`

#### Assess the performance of the assembly

*Gffcompare is a tool that can compare, merge, annotate and estimate accuracy of GFF/GTF files when compared with a reference annotation*

Using the StringTie merge mode, merge the assembly-generated GTF files to assess how well the predicted transcripts track to the reference annotation file. This step requires the TXT file,  (mergelist.txt). This file lists all of the file names to be merged. *Make sure ```mergelist.txt``` is in the StringTie program directory*.

++StringTie Arguments Used++:  
- --merge - Distinct from the assembly usage mode used above, in the merge mode, StringTie takes as input a list of GTF/GFF files and merges/assembles these transcripts into a non-redundant set of transcripts.
- -p - Specify number of processers
- -G - Specify reference annotation file. With this option, StringTie assembles the transfrags from the input GTF files with the reference sequences
- -o - Name of output file
- <mergelist.txt> - File listing all filenames to be merged. Include full path.

#Mcap

`
nano Mcap_mergelist.txt

SRR5311171.gtf
SRR5311255.gtf
SRR5311304.gtf
SRR5311305.gtf
SRR5311306.gtf
SRR5311307.gtf
SRR5311308.gtf
SRR5311310.gtf
SRR5311311.gtf
SRR5311312.gtf
SRR5311313.gtf
SRR5311342.gtf
SRR5311343.gtf
SRR5311344.gtf
SRR5311345.gtf
SRR5311346.gtf
SRR5311347.gtf
SRR5311348.gtf
`

`
stringtie --merge -p 8 -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Mcap.GFFannotation.gff -o Mcap_stringtie_merged.gtf Mcap_mergelist.txt
`


Now we can use the program gffcompare to compare the merged GTF to our reference genome.

++Gffcompare Arguments Used++:  
- -r - Specify reference annotation file
- -G - Compare all the transcripts in our input file ```stringtie_merged.gtf```
- -o - Prefix of all output files

#Mcap
`
gffcompare -r /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Mcap.GFFannotation.gff -o Mcap_gffComp Mcap_stringtie_merged.gtf
`

Some of the output files you will see are... 
- merged.stats
- merged.tracking
- merged.annotated.gtf
- merged.stringtie_merged.gtf.refmap
- merged.loci
- merged.stringtie_merged.gtf.tmap

Move all of the gffcompare output files to the output directory. We are most interested in the files ```merged.annotation.gtf``` and ```merged.stats```. The file ```merged.annotated.gtf``` tells you how well the predicted transcripts track to the reference annotation file and the file ```merged.stats``` file shows the sensitivity and precision statistics and total number for different features (genes, exons, transcripts).  



#### Compilation of GTF-files into gene count matrices

The StringTie program includes a script, ```prepDE.py``` that compiles your assembly files into gene and transcript count matrices. This script requires as input the list of sample names and their full file paths, (sample_list.txt).

. Run ```prepDE.py``` to merge assembled files together into a DESeq2-friendly version.

++StringTie prepDE.py Arguments Used++:  
- -i - Specify that input is a TXT file
- -g - Require output gene count file, default name is ```gene_count_matrix.csv```
- -t - Require output transcript count gene count file, default name is ```transcript_count_matrix.csv```

# Mcap
nano Mcap_sample_list.txt

T4-1	SRR5311171.gtf
T4-10	SRR5311255.gtf
T4-16	SRR5311304.gtf
T4-17	SRR5311305.gtf
T4-6	SRR5311306.gtf
T4-8	SRR5311307.gtf
T5-1	SRR5311308.gtf
T5-10	SRR5311310.gtf
T5-16	SRR5311311.gtf
T5-17	SRR5311312.gtf
T5-6	SRR5311313.gtf
T5-8	SRR5311342.gtf
T7-1	SRR5311343.gtf
T7-10	SRR5311344.gtf
T7-16	SRR5311345.gtf
T7-17	SRR5311346.gtf
T7-6	SRR5311347.gtf
T7-8	SRR5311348.gtf

`
prepDE.py -i Mcap_sample_list.txt -g Mcap_gene_count_matrix.csv 
`

# Obtain the gene name from the MSTRG ID

scp -r -P 2292  hputnam@kitt.uri.edu:/home/hputnam/Spawn/mapping/Mcap_stringtie_merged.gtf /Users/hputnam/MyProjects/Mcap_Spawning/RAnalysis/Data


scp -r -P 2292  hputnam@kitt.uri.edu:/home/hputnam/Spawn/mapping/Mcap_gene_count_matrix.csv /Users/hputnam/MyProjects/Mcap_Spawning/RAnalysis/Data






