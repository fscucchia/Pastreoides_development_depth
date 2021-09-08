
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
        
        - [Ruby](https://www.ruby-lang.org/en/documentation/installation/#building-from-source) >= 2.4 
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

#### InterProScan

As input, InterProScan requires reference protein sequences. 

*Note: Many fasta files willl use an asterisk to denote a STOP codon. InterProScan does not accept special characters within the sequences, so I removed them prior to running the program using the code below:*

```
cp Mcap.protein.fa ./Mcap.IPSprotein.fa
sed -i 's/*//g' Mcap.IPSprotein.fa
```

**Interproscan.sh: Executes the InterProScan Program**

```
interproscan.sh -version
interproscan.sh -f XML -i ../data/ref/Mcap.IPSprotein.fa -b ./Mcap.interpro.200824  -iprlookup -goterms -pa 
```

#### Trinotate

Place the Trinotate.sqlite database in the working directory

### Step 3: Map Kegg terms to genome  
*Uses KofamScan*

### Step 4: Compilation of the output of different methods

Done in RStudio.
