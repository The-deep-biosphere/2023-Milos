# Sequence processing pipeline
Here is a summary of the pipeline with the main points and explanations. Note that the real script is differently written, mainly because each step is being iterated over all samples. See the `pipeline.sh` file. The pipeline is based on the "alternative vsearch pipeline" proposed by [Torbjørn Rognes](https://github.com/torognes/vsearch/wiki/Alternative-VSEARCH-pipeline).

### Removal of the forward primers.
When we receive our sequences back, they are already demultiplexed, meaning that each sample has its own fastq file, and the adaptors and barcodes are already removed. However, the forward primer is still present. It is accepted that this part should be removed. Indeed, it does not contain any phylogenic information as all sequences are same for this part. We use the cutadapt tool. 
```bash
cutadapt -j 3 -g CAGCMGCCGCGGTAA --no-indels --discard-untrimmed -o $s.noprimers.fastq $f
```
-  `-j`: Run the command over 3 CPU.
- `-g `: The primer sequence to be removed.
- `--no-indels`:  We do not allow for insertion or deletion in the primer sequence, aka we want to find exactly the good sequence.
- `--discard-untrimmed`: Throw away the whole sequence if the primer is not found in it.
- `-o`: Name of the output.
- `$f`:  This is the input file.

### Trim at 220bp
The sequences that we get back from sequencing are decreasing in quality towards the end. This means the certainty of a given base pair becomes lower, and there are more errors. To use only the certain part of each sequence, we discard the end of them. We choose to trim at 220bp here. This is arbitrary, but decided after observing the quality of our sequences in a tool like FastQC. To trim our sequences, we use a command within the `vsearch` software. At the same time this will remove our sequences shorter than 220bp.
```bash
vsearch --fastq_filter $f \
	--fastq_maxns 0 \
	--fastq_trunclen 220 \
	--fastqout $s.trimmed.fastq
```
- `--fastq_filter`: Name of the input fastq file to be filtered.
- `--fastq_maxns`: Discard sequences with more than zero N's. N's are unknown base pairs.
- `--fastq_trunclen`: Trim the sequences at 220bp. Throw away the shorter ones.
- `--fastqout`: The name of the output file.

### Quality filtering
Fastq files contains both the sequences and the quality information related to each base pair. There are many ways to filter your sequences based on the quality. Here we use the max expected error principle, setting it to 2. That means the quality of our sequences allow us to expect no more than 2 errors in each sequence. Considering that we will later cluster our sequence at 97% similarity, 2 errors over 220bp should not make a major difference. Again, we use `vsearch`.
```bash
vsearch --fastq_filter $f \
	--relabel $s. \
	--fastq_maxee 2 \
	--fastaout $s.fasta \
	--fasta_width 0
```
- `--fastq_filter`: Name of the input fastq file to be filtered.
- `--relabel`: Relabel the name of the sequences within each file.
- `--fastq_maxee`:  The maximum amount of expected errors allowed.
- `--fastaout:` The name of the output fasta file.
- `--fasta_width`: The layout of the fasta file. Using 0, the sequence is written fully on one line.

### Dereplicate
We now have fasta files, meaning that only the sequences are present in the file, no more the quality information. We anyway do not need it anymore... Now we gonna dereplicate each files. This means that we gonna fuse together sequences that are 100% same. The aim is to decrease the size of the file, and therefore alleviate the computer power needed to process the data.
```bash 
vsearch --derep_fulllength $f \
--strand plus \
--sizeout \
--relabel $s. \
--fasta_width 0 \
--output $s.derep.dfa
```
- `--derep_fulllength`: What we want to do, dereplicate our fasta files.
- `--strand`: Do we want to compare only the same plus strand, or also the reverse complement? 
- `--sizeout`: Add the amount of time this sequence is present in the name of the sequence.
- `--relabel`: Relabel the name of the sequences within each file.
- `--fasta_width`: The layout of the fasta file. Using 0, the sequence is written fully on one line.
- `--output`: Name of the output file.

### Concatenate files
At this point we are done with the analysis of each sample separately. We can therefore concatenate all files together. This uses the `cat` command from Unix. It takes all the dfa files in the directory and pool them together into a new file.
```bash
cat *.dfa > AllSeqs.fa
```

### Plan from here
The "alternative vsearch pipeline" that we use here is based on the approach that usearch uses. In short, OTUs at 97% are produced from all sequences except the singleton sequences. Then, these OTUs are checked for chimeras. Finally, all original sequences (including singletons) are mapped to these created OTUs using the threshold of  97% similarity.

### Dereplicate again
Now that we have pooled all sequences together, there are likely some replicates again. Therefore we can dereplicate once more.
```bash
vsearch --derep_fulllength AllSeqs.fa \
	--minuniquesize 2 \
	--sizein \
	--sizeout \
	--fasta_width 0 \
	--output AllSeqs.derep.fasta
```
- `--minuniquesize`: The minimum of amount a sequence need to be present to be kept. Here we do not want to keep the singletons, so we use 2.
- `--sizein`: To say that there is quantity information already in the name of each sequence.

### Cluster at 97%
Now we can pool together sequences that are 97% similar into OTUs. While using raw sequences is a possibility, it adds some noise and computer power needs. 97% is a commonly used threshold.
```bash
vsearch --cluster_size AllSeqs.derep.fasta \
	--id 0.97 \
	--strand plus \
	--sizein \
	--sizeout \
	--fasta_width 0 \
	--centroids Centroids.fasta
```
- `--cluster_size`: What we want to do, cluster sequences.
- `--id`: The similarity that we want to use as a threshold, 97%.
- `--centroids`: The name of the output, a fasta files containing a fake centroid sequence for each OTU.

### Sort centroids by size
Before looking for chimera, we need to organise our centroids by size (quantity).
```bash
vsearch --sortbysize Centroids.fasta \
--sizein \
--sizeout \
--fasta_width 0 \
--output Centroids.sorted.fasta
```
- `--sortbysize`: The input file to be sorted.

### Look for chimeras.
Chimeras are artificial products of PCR amplification. We can look for them within the dataset (the algorithm looks for sequences made of several other sequences), or against a dataset (the algorithm looks for sequences made of several sequences from this dataset). Here we run both.
``` bash
vsearch --uchime_denovo Centroids.sorted.fasta \
	--sizein \
	--sizeout \
	--fasta_width 0 \
	--qmask none \
	--nonchimeras Centroids.denovo.nonchimeras.fasta \

vsearch --uchime_ref Centroids.denovo.nonchimeras.fasta \
	--db SILVA_138.1_SSURef_NR99_tax_silva.fasta \
	--sizein \
	--sizeout \
	--fasta_width 0 \
	--qmask none \
	--dbmask none \
	--nonchimeras Centroids.nonchimeras.fasta
```
- `--uchime_denovo`: The file that needs to be run through *denovo* chimera detection.
- `--qmask`: Masking is a method that detect and replaces areas with low complexity in the sequence (ie, that have mainly the same base pair many times, or repeats) and replace it by N's. The rational is that it may create issues in alignments and chimera detection. We do not want masking here.
- `--nonchimeras`: The output file.
- `--db`: The database to be used. Here we use SILVA138.
- `--dbmmask`: Like `--qmask`, but for the database

### Relabel OTUs as OTU_xxx
That is just for simplicity...
``` bash
vsearch --fastx_filter Centroids.nonchimeras.fasta \
--sizein \
--fasta_width 0 \
--relabel OTU_ \
--fastaout Otus.fasta
```

### Map OTUs
Now that we have a clean set of OTUs without chimeras, we can map our original sequences (just after the `cat` command) to them, and get our OTU table.
``` bash
vsearch --usearch_global AllSeqs.fa \
--db Otus.fasta \
--id 0.97 \
--strand plus \
--sizein \
--sizeout \
--fasta_width 0 \
--qmask none \
--dbmask none \
--otutabout Otutab.txt
```
- `--otutabout`: The name of the output.

### Clean with Lulu
We do one extra step of cleaning here, using the [Lulu](https://github.com/tobiasgf/lulu) algorithm. The software pools together OTUs that it believes are the results of sequencing errors.
First one need to produce a database using `makebalstdb`, and then blast our fasta file against itself using `blastn`.
``` bash
makeblastdb -in Otus.fasta -dbtype nucl

blastn -db Otus.fasta -outfmt '6 qseqid sseqid pident' \
	-out LULU_match_list.txt -qcov_hsp_perc 95 -perc_identity 95 \
	-query Otus.fasta -num_threads 12
```
- `-dbtype`: type of sequence, nucleotides here.
- `-db`: The database used. The one we just created in our case.
- `-outfmt`: The format of the output.
- `-out`: The output file.
- `-qcove_hsp_perc`: The coverage percentage needed per alignment.
- `-perc_identity`: The similarity percentage needed.
- `-query`: Our fasta file.
- `-num_threads`: Amount of CPUs to use.
``` R
setwd(".")
require(lulu)
require(methods)

matchlist = read.delim("LULU_match_list.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
otus.all = read.delim("Otutab.sorted.txt",row.names=1,header=T,sep="\t")
curated_result <- lulu(otus.all,matchlist, minimum_match = 97)
lulus = curated_result$curated_table
write.table(data.frame("OTU"=rownames(lulus),lulus),"Otutab_curated.tsv", row.names=FALSE, quote=F, sep="\t")
```
As Lulu only makes modifications to the OTU table, the following python3 script can be run to filter the fasta file too.
``` python
#!/usr/bin/env python

import sys

ot = open(sys.argv[2],"r")
fa = open(sys.argv[1],"r")
otuIDs = set()

firstLine = True
for line in ot:
    l = line.replace("\n","")
    if firstLine:
        firstLine = False
    else:
        otuID = l.split("\t")[0]
        otuIDs.add(otuID)

ot.close()
        
for line in fa:
    l = line.replace("\n","") 
    if l.startswith(">"):
        seqId = l[1:].split(" ")[0]
        if ";size" in seqId:
            seqId = seqId[:seqId.find(";size")]
        if seqId in otuIDs:
            otu_inc=True
            print(">%s" % seqId)
        else:
            otu_inc=False
    elif otu_inc:
        print(l)
    
fa.close()
```

### Taxonomic assignments
Finally, we can give a taxonomic assignment to the sequences that we have. For this we use the [Crest4](https://github.com/xapple/crest4) software that uses a last common ancestor algorithm. Assignments are given based on a modified SILVA138 database.
``` bash
crest4 -f Otus_curated.fasta
```

