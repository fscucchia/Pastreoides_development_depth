
###### Concatenate proteome databases ######

/data/home/mass/tmass/stylophora_yotam/data/dbs/prot5/merge #old source dirs
/data/home/mass/tmass/stylophora_yotam/data/dbs

# Databases list 
/data/home/mass/tmass/stylophora_yotam/data/dbs/prot5/prot5.symbB.v1.2.augustus.prot.fa  #head SymbiodiniumMinutum_symbB.v1.2.000001.t1|scaffold5.1|size591573|1
/data/home/mass/tmass/stylophora_yotam/data/dbs/prot5/prot5.GCA_001939145.1_ASM193914v1_protein.faa  #head SymbiodiniumMicroadriaticumNCBI_OLP20498.1
/data/home/mass/tmass/stylophora_yotam/data/dbs/prot5/prot5.Smic.genome.annotation.pep.longest.fa  #head SymbiodiniumMicroadriaticum_Smic1
/data/home/mass/tmass/stylophora_yotam/data/dbs/prot5/prot5.SymbF.Gene_Models.PEP.fasta  #head SymbiodiniumKawa_SymbF.scaffold100.2
/data/home/mass/tmass/stylophora_yotam/data/dbs/prot5/prot5.SymbC1.Gene_Models.PEP.fasta  #head SymbiodiniumGor_SymbC1.scaffold1.10
/data/home/mass/tmass/stylophora_yotam/data/dbs/prot5/prot5.symC_aug_40.aa.fa  #head SymbiodiniumCladeC_s4526_g1.t1
/data/home/mass/tmass/stylophora_yotam/data/dbs/prot5/prot5.Cladocopium-C15-Porites-lutea-holobiont_SymbC15_plutea_v2.1.fna.evm.prot.final.faa  #head CladocopiumC15_evm.model.NODE_100098_length_1590_cov_40.8455_ID_76299152.1
/data/home/mass/tmass/stylophora_yotam/data/dbs/prot5/prot5.syma_aug_37.aa.fasta  #head SymbiodiniumCladeA3_s3212_g1.t1
/data/home/mass/tmass/stylophora_yotam/data/dbs/prot5/Breviolum_minutum_combined.fasta  #head BreviolumMinutum_10000|brevM.v1.2.021059.t1
/data/home/mass/tmass/stylophora_yotam/data/dbs/prot5/symbd_genemodels_prot.fa  #head

#<new symbB database, to switch with the one of marinegenomics. This comes from genome.jgi
/data/home/mass/tmass/stylophora_yotam/data/dbs/Breviolum_minutum_fromGenome.jgi
#<concatenate the two files of this one
$ cd /data/home/mass/tmass/stylophora_yotam/data/dbs/Breviolum_minutum_fromGenome.jgi
$ cat Bremi1_GeneModels_AllModels_20200813_aa.fasta Bremi1_GeneModels_ExternalModels_aa.fasta > /data/home/mass/tmass/stylophora_yotam/data/dbs/Breviolum_minutum_fromGenome.jgi/Breviolum_minutum_combined.fasta
#I created a new prot5.Breviolum_minutum_combined.fasta file to replace the names "jgi|Bremi1|10000|symbB.v1.2.021059.t1" with "BreviolumMinutum_" and "symbB" with "brevM", to match the other prot5 files that Assaf made.

#<more symbiont databases from González-Pech et al. 2021 https://espace.library.uq.edu.au/view/UQ:f1b3a11
/data/home/mass/tmass/stylophora_yotam/data/dbs/S_microadriaticum_CassKB8_tar
/data/home/mass/tmass/stylophora_yotam/data/dbs/S_tridacnidorum_CCMP2592_tar.gz
/data/home/mass/tmass/stylophora_yotam/data/dbs/S_natans_CCMP2548_tar.gz
/data/home/mass/tmass/stylophora_yotam/data/dbs/S_microadriaticum_04_503SCI_03_tar.gz
/data/home/mass/tmass/stylophora_yotam/data/dbs/S_necroappetens_CCMP2469_tar.gz
/data/home/mass/tmass/stylophora_yotam/data/dbs/S_linucheae_CCMP2456_tar.gz
/data/home/mass/tmass/stylophora_yotam/data/dbs/S_pilosum_CCMP2461_tar.gz #corrupted file

#<clade D from marinegenomics
/data/home/mass/tmass/stylophora_yotam/data/dbs/Durusdinium sp_marinegenomics.oist
#I created a new prot5.symbd_genemodels_prot.fa file to replace the names ">g1.t1" with ">SymbiodiniumDurusdinium_g1.t1", to match the other prot5 files that Assaf made.

# Commands
```
$ cd /data/home/mass/tmass/stylophora_yotam/data/dbs/prot5
$ cat prot5.Breviolum_minutum_combined.fasta prot5.symbd_genemodels_prot.fa S_linucheae_CCMP2456.proteins.faa S_microadriaticum_04_503SCI_03.proteins.faa S_microadriaticum_CassKB8.proteins.faa S_natans_CCMP2548.proteins.faa S_necroappetens_CCMP2469.proteins.faa S_tridacnidorum_CCMP2592.proteins.faa prot5.GCA_001939145.1_ASM193914v1_protein.faa prot5.Smic.genome.annotation.pep.longest.fa prot5.SymbF.Gene_Models.PEP.fasta prot5.SymbC1.Gene_Models.PEP.fasta prot5.symC_aug_40.aa.fa prot5.Cladocopium-C15-Porites-lutea-holobiont_SymbC15_plutea_v2.1.fna.evm.prot.final.faa prot5.syma_aug_37.aa.fasta > /data/home/mass/tmass/stylophora_yotam/data/dbs/prot5/merge/symbiont_combined_changedB_addedMore.fa
```

# Commands take 2 (excluded the ones from González-Pech et al. 2021)
```
$ cd /data/home/mass/tmass/stylophora_yotam/data/dbs/prot5
$ cat prot5.symbB.v1.2.augustus.prot.fa prot5.symbd_genemodels_prot.fa prot5.GCA_001939145.1_ASM193914v1_protein.faa prot5.Smic.genome.annotation.pep.longest.fa prot5.SymbF.Gene_Models.PEP.fasta prot5.SymbC1.Gene_Models.PEP.fasta prot5.symC_aug_40.aa.fa prot5.Cladocopium-C15-Porites-lutea-holobiont_SymbC15_plutea_v2.1.fna.evm.prot.final.faa prot5.syma_aug_37.aa.fasta prot5.GCF_002571385.1_Stylophora_pistillata_v1_protein.faa prot5.Pocillopora-amicornis_GCF_003704095.1_ASM370409v1_genomic.prot.fasta prot5.Galaxea-fascicularis_gfas_1.0.proteins.prot.fasta prot5.newCoral_stylo_NCBI_hints_predictions3.aa.100aa.fasta > /data/home/mass/tmass/stylophora_yotam/data/dbs/prot5/merge/symbiont_combined_changedB_addedMore.fa
```

<Use Diamond `makedb` command to format the database in a Diamond-friendly format> used the script diamond_makedb.sh
```
$ /data/home/mass/tmass/stylophora_yotam/data/dbs/prot5/merge
$ diamond makedb --in nr.gz -d nr
$ diamond dbinfo -d nr.dmnd
```

###### Run diamond and heatmap ######
```
$ bash diamond2_RUN2.sh 1
```

#for the heatmap run the R script
```
$ conda activate R_setting6
$ srun --pty -p hiveunlim,hive7d,hive1d R
```
