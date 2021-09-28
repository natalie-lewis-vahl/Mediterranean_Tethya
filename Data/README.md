### Dataset used for the analysis of *L. chondrodes* microbiome: 

Below a description of the datasets used and the analyses for which these datasets were used is provided

`cbas_tempVSctrl.otutab.csv`

This dataset can be considered the base dataset generated and used. It is an **OTU by sample** table providing **counts**. This dataset is not filtered in anyway, meaning is the derived from the `vsearch` output, and is the dataset that was for **all** analyses in the study.

Three files are related to this "base" dataset:

The file `cbas_tempVSctrl.otutab.info` is necessary to analyze the dataset `cbas_tempVSctrl.otutab.csv` with the DESeq2 and provides the information about the samples as required by DESeq2. We used DESeq2 in our enrichment/depletion.

The file `cbas_otu_taxonomy.csv` gives the taxonomic assignment for the OTUs found using `vsearch`. This file was used to combine the counts per OTU with the taxonomy of each OTU and derive richness and abundance per phylum chart, etc. Also, using this file we derived the core microbiome of *L. chondrodes*

The file `cbas_core_otus_final_names.fasta` includes the representative V4 16S rRNA sequence for all **core OTUs** detected in *L. chondrodes* with their taxonomy and a "final name". This "final name" comes from the fact that `vsearch` assigns names to the OTUs based on their abundance. Thus, OTU1 can be different if different data are used to generate the OTU table. We are "fixing" the names of the core community here found to allow future work on this species to refer to these OTUs in a standard way.

The following dataset are also provided and were used to assess whether normalizing the count data using bacterial load data had an effect on the conclusions reached:

`cbas_tempVSctrl.otutab.bactloadCorrected_valuesOnly.csv`

This dataset "corrects" the count data using the bacterial load data obtained via RT-qPCR for each sample and includes all the OTUs detected by `vsearch`.  

`cbas_tempVSctrl.otutab.bactloadCorrected_valuesOnly_wFreqs_more1PctCounts.csv`

This dataset "corrects" the count data using the bacterial load data obtained via RT-qPCR for each sample but includes only OTUs with a frequency of at least one percent across **all** samples.

`cbas_tempVSctrl.otutab.bactloadCorrected_valuesOnly_wFreqs_more50Counts.csv`

This dataset "corrects" the count data using the bacterial load data obtained via RT-qPCR for each sample but includes only OTUs with a frequency of at least 50 corrected counts **all** samples.

These datasets are provided for the sake of completeness and were used **only** to assess whether correcting by bacterial load had an effect on the results obtained using the raw uncorrected counts. The R scripts provided in this repository have commented lines of code pointing to the "corrected" data files used.

### A note on sample names:

Sample names are alphanumeric codes including information on **color morph** (G=Green or P=Purple), **color morph replicate** (1,2,...,n), **treatment** (C=Control, T=High temperature), **Treatment tank** (1,2,...,n)

code **G1C1** represents Green sponge number 1 in Control tank number 1.


Sergio