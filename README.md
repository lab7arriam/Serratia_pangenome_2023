# The man, the plant, and the insect: shooting host specificity determinants in <i>Serratia marcescens</i> pangenome
The repository with working scripts, figures, data, and supplementary materials for pangenome analysis of 73 <i>Serratia marcescens</i> genomes.

<img src="https://github.com/lab7arriam/Serratia_pangenome_2023/blob/main/pics/fig3.svg?sanitize=true" width="700" height="700">


## Contents 

This repository contains scripts used for statistical analysis of pangenome, phylogeny, and phenotypical associations. Please consult the Methods section in the paper for extra details:

Shikov A.E., Merkushova A.V., Nizhnikov A.A., Antonets K.S. (2023) The man, the plant, and the insect: shooting host specificity determinants in <i>Serratia marcescens</i> pangenome


## Figures
Figures are available in the `pics/` directory. For the description, please consult the Results section of the article.

## Supplemenatry
All supplementary material is located in the `supplementary/` directory. 
`Supplementary_figures.docx` is a Microsoft Word document with a description of all supplementary figures.
`Supplementary_tables.xlsx` is an Exel table with all supplementary tables.

* `supplementary/pics/` folder contains individual supplementary figures.
* For convenience, in the `supplemenatry/tables/` folder individual supporting tables in CSV format are provided. For detailed description ogf the tables, please consult the `Supplementary_tables.xlsx` file.

## Data
Analyzed data are included in the `data/` directory.

* `data/pangenomes/` folder contains gene presence/absence tables obtained using Panaroo, PEPPAN, and Roary.
* `data/trees/` folder contains phylogenetic inferences in Newick format based on core gene alignments (aligned with MAFFT of Prank), either partitioned or unpartitioned, and the results of clusterizations based on gene presence/absence patterns.
* `data/scoary/` folder includes raw Scoary results attributed to particular hosts, namely, humans, plants, and insects.

## Scripts
The `scripts/` directory includes all code used for pangenome analysis.
* `IO_lib.py` is an ancillary script with functions for processing CSV files;
* `download_serratia_assemblies.sh` was used to download  <i>Serratia marcescens</i> genomes in fna and gff formats from the RefSeq database;
* `search_entrez.py` was applied to extract metadata of studied assemblies;
* `process_serratia_table.py` was used to calculate mean GC content, the number of CDS, and genome length and summarize assemblies' metadata to a single table;
* `add_numbers_to_presence.py` was utilized to add the number of genomes in which gene clusters are presented to gene presence/absence tables;
* `add_assembly_codes_to_peppan.py`- script for renaming labels in Newick tree generated by PEPPAN on the basis of gene presence/absence patterns;
* `calculate_mash_dist.py` - script for creating matrix based on ANI values using mash utility;
* `make_accessory_alignment.py` was used to generate pseudoalignment representing the presence/absence of accessory genes;
* `create_paritions.py` was applied to obtain partitioned alignment as well as the partitioning scheme for all core genes in the pangenome;
* `calculate_mean_support_new.py` - script for calculating mean bootstrap support of the phylogenetic tree;
* `compare_tree_topologies.py` was utilized for comparing phylogenetic inferences using quartet distance metrics;
* `make_table_for_scoary.py` - script for generating a table with phenotypic traits used for Scoary analysis;
* `filter_scoary_hits.py` was applied to extract positively associated gene clusters from Scoary-generated results;
* `calculate_number_of_host_genes_for_assemblies.py`- script for summarizing the number of host specificity factors per assembly;
* `add_genome_nums_to_VFDB.py` was applied to add the number of genomes in which the cluster with virulence factors was found;
* `filter_virDB_hits.py` - script for filtering hits matching the VFDB database using identity and coverage;
* `split_fasta_to_chunks.py` - script for re-naming Panaroo-attributed gene codes to real protein accession numbers;
* `assert_groups_to_genes_pos_scoary_hits.py` was used for grouping protein sequences into 8 categories, namely, pangenomic (core, accessory, and unique), virulence (core and accessory), and specificity genes attributed to a particular host (human, insect, plant);
* `prepare_ontology_data.py` was utilized to generate a universe of functional annotation terms for the over-representation test;
* `get_groups_genes_extracted_for_GO.py` - script for grouping GO terms according to eight aforementioned groups;
* `perform_enrichment_anaysis.py` was applied to perform an over-representation test using hypergeometric distribution;
* `create_matrix_for_enrichments.py` was used for generating a dissimilarity matrix based on significant functional terms;
* `Serratia_pangenome_plots.R` - main R script for generating all figures, clusterizations, assessing pangenome openness, graph-wise testing GO terms using topGO, summarizing data on virulences/specificity factors and mobile genetic elements.

