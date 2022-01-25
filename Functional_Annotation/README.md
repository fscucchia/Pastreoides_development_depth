
### Functional annotation of *Porites astreoides* reference genome

Homologous sequences are first found using the program BLAST, which searches the reference sequence against a database of reviewed protein sequences. Once homologous sequences are found, genes are tagged with known Gene Ontology terms for the homologous sequences. 

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

### Step 1: Find homologous sequences

#### Download/Update nr database

Run the script [`nr_NCBI_download.sh`]() to download the nr database from [NCBI](ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz). 
Then, use Diamond's ```makedb``` command to format the database in a Diamond-friendly format. You can also use the command ```dbinfo``` to find version information for the database.

```
diamond makedb --in nr.gz -d nr
diamond dbinfo -d nr.dmnd
```

#### Run DIAMOND Search

As input, DIAMOND requires your reference sequences (either protein or CDS nucleotides), and a path to your nr database. 

*It may take a few days depending on the number of sequences you have the M. cap genome (~63,000 genes) took 4.5 days*


**Blastx: Align translated DNA query sequences against a protein reference database**

```
#Run sequence alignment against the nr database
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q ../data/ref/Mcap.mRNA.fa -o Mcap.annot.200806 -f 100  -b 20 --more-sensitive -e 0.00001 -k1 --unal=1
```
---

### Step 2: Map Gene ontology terms to genome  

#### Uniprot

#### EggNog

Place the Trinotate.sqlite database in the working directory

### Step 3: Compilation of the output of different methods

Done in RStudio.
