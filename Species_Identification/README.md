The following document contains the bioinformatic pipeline used for identifying the endosymbiotic algae species based on _P. astreoides_ raw RNA sequences.

---

**Tools used**  

Blastx search: [Diamond](https://github.com/bbuchfink/diamond)                                                                                                     
Final heatmap: [R](https://cran.r-project.org/)
---

### Concatenate proteome databases 

#### Databases list 
~/prot5.symbB.v1.2.augustus.prot.fa  #head SymbiodiniumMinutum_symbB.v1.2.000001.t1|scaffold5.1|size591573|1
~/prot5.GCA_001939145.1_ASM193914v1_protein.faa  #head SymbiodiniumMicroadriaticumNCBI_OLP20498.1
~/prot5.Smic.genome.annotation.pep.longest.fa  #head SymbiodiniumMicroadriaticum_Smic1
~/prot5.SymbF.Gene_Models.PEP.fasta  #head SymbiodiniumKawa_SymbF.scaffold100.2
~/prot5.SymbC1.Gene_Models.PEP.fasta  #head SymbiodiniumGor_SymbC1.scaffold1.10
~/prot5.symC_aug_40.aa.fa  #head SymbiodiniumCladeC_s4526_g1.t1
~/prot5.Cladocopium-C15-Porites-lutea-holobiont_SymbC15_plutea_v2.1.fna.evm.prot.final.faa  #head CladocopiumC15_evm.model.NODE_100098_length_1590_cov_40.8455_ID_76299152.1
~/prot5.syma_aug_37.aa.fasta  #head SymbiodiniumCladeA3_s3212_g1.t1
~/prot5.symbd_genemodels_prot.fa  #head SymbiodiniumDurusdinium_g1.t1

```
cd ~
cat prot5.symbB.v1.2.augustus.prot.fa prot5.symbd_genemodels_prot.fa prot5.GCA_001939145.1_ASM193914v1_protein.faa prot5.Smic.genome.annotation.pep.longest.fa prot5.SymbF.Gene_Models.PEP.fasta prot5.SymbC1.Gene_Models.PEP.fasta prot5.symC_aug_40.aa.fa prot5.Cladocopium-C15-Porites-lutea-holobiont_SymbC15_plutea_v2.1.fna.evm.prot.final.faa prot5.syma_aug_37.aa.fasta prot5.GCF_002571385.1_Stylophora_pistillata_v1_protein.faa prot5.Pocillopora-amicornis_GCF_003704095.1_ASM370409v1_genomic.prot.fasta prot5.Galaxea-fascicularis_gfas_1.0.proteins.prot.fasta prot5.newCoral_stylo_NCBI_hints_predictions3.aa.100aa.fasta prot5.Orbicella-faveolata_GCF_002042975.1_ofav_dov_v1_genomic.prot.fasta prot5.Porites-lutea_plut2v1.1.proteins.fasta > /data/home/mass/fscucchia/databases/symbiont_combined_prot/corals_symbionts_combined_prot.fa
```

### Format the protein database for Diamond

Run [`diamond_makedb.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/Species_Identification/diamond_makedb.sh) that formats the above created protein database in a Diamond-friendly format.

```
cd /data/home/mass/fscucchia/databases/symbiont_combined_prot
diamond makedb --in corals_symbionts_combined_prot.fa -d corals_symbionts_combined_prot.fa_makedb
```

### Run Diamond blastx

Run [`diamond2_RUN2.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/Species_Identification/diamond2_RUN2.sh) that uses the design [`design_diamond.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/Metadata/design_diamond.sh) to perform a blastx search by translating the nucleotide sequences and searching the above concatenated protein sequences. 
Here I used an evalue of 0.01 as cutoff.

### Visualize Diamond output

The R script [`diamond_heatmap.R`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/Species_Identification/diamond_heatmap.R)(run on Linux terminal within a conda environment) arranges the diamond blastx output into heatmaps of i) count of mapped reads and ii) proportion of mapped reads for both coral and endosymbionts species. It uses the design table [`metadata_diamond_Past.csv`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/Metadata/metadata_diamond_Past.csv), which contains information about the sample name, depth and developmental stage.
