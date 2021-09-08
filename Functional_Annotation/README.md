
### Functional annotation of *Porites astreoides* reference genome

Functional annotation tags putative genes in a reference genome or transcriptome with the known functions of homologous genes in other organisms. Homologous sequences are first found using the program BLAST, which searches the reference sequence against a database of reviewed protein sequences. Once homologous sequences are found, genes are tagged with known Gene Ontology or Kegg Pathway terms for the homologous sequences. Annotation allows us to better understand the biological processes that are linked to genes of interest.

---

#### Create and activate a conda environment
```
conda create -n functional_annotation
conda activate functional_annotation
```

#### Install all necessary programs within your conda environment

- [DIAMOND](http://www.diamondsearch.org/) Search Program v2.0.0
```
wget http://github.com/bbuchfink/diamond/releases/download/v2.0.11/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
```

- [Trinotate](https://informatics.fas.harvard.edu/trinotate-workflow-example-on-odyssey.html)
```
wget https://github.com/Trinotate/Trinotate/archive/v2.0.2.tar.gz -O Trinotate-2.0.2.tar.gz
tar xvf Trinotate-2.0.2.tar.gz

wget https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Trinotate_v2.0_RESOURCES/uniprot_sprot.trinotate_v2.0.pep.gz
wget https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Trinotate_v2.0_RESOURCES//Pfam-A.hmm.gz
wget "https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Trinotate_v2.0_RESOURCES/Trinotate.sprot_uniref90.20150131.boilerplate.sqlite.gz" -O Trinotate.sqlite.gz
```

- [InterProScan](https://github.com/ebi-pf-team/interproscan) v.5.46-81.0
    - Requires Java v11.0.2
<followed the instructions [here](https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html#obtaining-a-copy-of-interproscan)>    

- [KofamScan](https://github.com/takaram/kofam_scan) v1.3.0
```
conda install -c bioconda kofamscan
```
- KofamScan requires:
        
        - Ruby >= 2.4 (https://www.ruby-lang.org/en/documentation/installation/#building-from-source)
        ```
           ./configure --prefix=/data/home/mass/fscucchia/programs/ruby_installed
           make
           make install
        ```   

        - GNU parallel
        <compiled from [source](https://medium.com/analytics-vidhya/simple-tutorial-to-install-use-gnu-parallel-79251120d618)>

        - util-linux v2.34
        
        - HMMER >= 3.1
        ```
           conda install -c bioconda hmmer
        ```   

### Step 1: Find homologous sequences

#### i) Download/Update nr database

The nr, or non-redundant, database is a comprehensive collection of protein sequences that is compiled by the National Center for Biotechnology Information (NCBI). It contains non-identical sequences from GenBank CDS translations, PDB, Swiss-Prot, PIR, and PRF, and is updated on a daily basis.

Run the script [`nr_NCBI_download.sh`]() to download the nr database from [NCBI](ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz). 
Then, use Diamond's ```makedb``` command to format the database in a Diamond-friendly format. You can also use the command ```dbinfo``` to find version information for the database.

```
diamond makedb --in nr.gz -d nr
diamond dbinfo -d nr.dmnd
```

#### ii) Run DIAMOND Search

With your updated nr database, you can run DIAMOND. DIAMOND, like NCBI's BLAST tool, is a sequence aligner for protein and translated DNA searches. The tool is optimized for large datasets (>1 million sequences), and is 100x-20,000x faster BLAST without losing sensitivity. 

As input, DIAMOND requires your reference sequences (either protein or CDS nucleotides), and a path to your nr database. Below, I output the results to DIAMOND format and then convert to XML. I highly suggest the DIAMOND export format, as it can easily be converted into any other format that you may need for your analysis. 

*It may take a few days depending on the number of sequences you have the M. cap genome (~63,000 genes) took 4.5 days*


**Blastx: Align translated DNA query sequences against a protein reference database**

*Options:*
- **-d** - Path to nr database  
- **-q** - Path to reference fasta file  
- **-o** - Base output name  
- **-f** - Output format. **100**=DIAMOND output.     
- **-b** - Block size in billions of sequence letters to be processed at a time. Larger block sizes increase the use of memory and temporary disk space, but also improve performance. Set at **20**. 20 is the highest recommended value. CPU is about 6x this number (in GB).  
- **--more-sensitive** - slightly more sensitive than the --sensitive mode.  
- **-e** - Maximum expected value to report an alignment. **1E-05** is the cut-off typically used for sequence alignments.  
- **-k** - Maximum top sequences. Set at **1** because I only wanted the top sequence reported for each gene.  
- **--unal** - Report unaligned queries (yes=**1**, no=0). 

**View: Generate formatted output from DAA files**

*Options:*  
- **-a** - Path to input file  
- **-o** - Base output name  
- **-f** - Output format. **5**=XML output. 

```
#Run sequence alignment against the nr database
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q ../data/ref/Mcap.mRNA.fa -o Mcap.annot.200806 -f 100  -b 20 --more-sensitive -e 0.00001 -k1 --unal=1

#Converting format to XML format for BLAST2GO
diamond view -a Mcap.annot.200806.daa -o Mcap.annot.200806.xml -f 5
```
---

### Step 2: Map Gene ontology terms to genome  
*Can be done concurrently with Steps 1 and 3*

#### i) InterProScan

InterProScan searches the database InterPro database that compiles information about proteins' function multiple other resources. I used it to map Kegg and GO terms to my Mcap reference protein sequences.

The commands that I used are below, however, the script I used to execute this on bluewaves is available on my project [repository](https://github.com/echille/Montipora_OA_Development_Timeseries/blob/master/Scripts/IPS.sh). As input, InterProScan requires reference protein sequences. Below, I output the results to XML format, as it is the most data-rich output file and can be used as input into Blast2GO.

*Note: Many fasta files willl use an asterisk to denote a STOP codon. InterProScan does not accept special characters within the sequences, so I removed them prior to running the program using the code below:*

```
cp Mcap.protein.fa ./Mcap.IPSprotein.fa
sed -i 's/*//g' Mcap.IPSprotein.fa
```

**Interproscan.sh: Executes the InterProScan Program**

*Options:*

- **-version** - displays version number
- **-f** - output format
- **-i** - the input data
- **-b** - the output file base
- **-iprlookup** - enables mapping
- **-goterms** - map GO Terms
- **-pa** - enables Kegg term mapping

```
interproscan.sh -version
interproscan.sh -f XML -i ../data/ref/Mcap.IPSprotein.fa -b ./Mcap.interpro.200824  -iprlookup -goterms -pa 
```

#### ii) Trinotate

Place the Trinotate.sqlite database in the working directory






### Step 3: Map Kegg terms to genome  
*Uses KofamScan. Can be done concurrently with Steps 1 and 2. Currently Troubleshooting*



### Step 5: Compilation of the output of different methods

Done in RStudio. See RMarkdown [page](https://github.com/echille/Montipora_OA_Development_Timeseries/blob/master/RNAseq_Analyses/annot/Mcap_annot_compile.html).
