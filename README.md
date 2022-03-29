# The Car Tank Lid Bacteriome
This is a repository for the code used in "The car tank lid bacteriome: a reservoir of bacteria with potential in bioremediation of fuel", by Vidal-Verdú *et al.* (2022).

**Abstract**
Bioprospecting of microorganisms suitable for bioremediation of fuel or oil spills is often carried out in contaminated environments such as gas stations or polluted coastal areas. Using Next Generation Sequencing (NGS) we analyzed the microbiota thriving below the lids of the fuel deposits of diesel and gasoline cars. The microbiome colonizing the tank lids differed from the diversity found in other hydrocarbon-polluted environments, with *Proteobacteria* being the dominant phylum and without clear differences between gasoline or diesel fueled vehicles. We observed differential growth when samples were inoculated in cultures with gasoline or diesel as a main carbon source, as well as an increase in relative abundance of the genus Pseudomonas in diesel. A collection of culturable strains was established, mostly *Pseudomonas, Stenotrophomonas, Staphylococcus* and *Bacillus* genera. Strains belonging to *Bacillus, Pseudomonas, Achromobacter* and *Isoptericola* genera showed a clear diesel degradation pattern when analyzed by GC-MS, suggesting their potential use for bioremediation and a new species of Isoptericola was further characterized as hydrocarbon degrader.

# Qiime2 code
## Code for the main analysis

**Import demultipliex sequences**
```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path Manifest \
  --output-path ./demux-paired-end.qza \
  --input-format PairedEndFastqManifestPhred33
```

**Check parameters of the reads**
```bash
qiime demux summarize \
  --i-data ./demux-paired-end.qza \
  --o-visualization ./demux-paired-end.qzv
```
  
**Sequence quality control and feature table construction**
```bash
 qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./demux-paired-end.qza \
  --verbose \
  --p-trim-left-f 21 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 290
  --p-trunc-len-r 226
  --p-n-threads 8 \
  --o-table ./table.qza \
  --o-representative-sequences ./rep-seqs.qza \
  --o-denoising-stats ./denoising-stats.qza 

qiime metadata tabulate \
  --m-input-file ./denoising-stats.qza \
  --o-visualization ./stats-dada2.qzv
```

**Taxonomic classification. Database SILVA 132**
```bash
qiime feature-classifier classify-sklearn \
  --i-classifier /path/to/SILVA/132 \
  --i-reads ./rep-seqs.qza \
  --verbose \
  --p-n-jobs 8 \
  --o-classification ./taxonomy.qza
```  

**Export the qiime2 BiomTables (table.qza, taxonomy.qza) in such a way that they can be loaded into the R package phyloseq**

```bash
# Export table.qza and taxonomy.qza
qiime tools export \
 --input-path ./table-all.qza \
  --output-path  ./table

# Export taxonomy.qza and taxonomy.qza
qiime tools export \
 --input-path ./taxonomy.qza \
  --output-path  ./taxonomy


# Modify the biom-taxonomy tsv headers: change header "Feature ID" to "#OTUID"; "Taxon" to "taxonomy"; and "Confidence" to "confidence"
sed -i -e 's/Feature ID/#OTUID/g' ./taxonomy/taxonomy.tsv
sed -i -e 's/Taxon/taxonomy/g' ./taxonomy/taxonomy.tsv
sed -i -e 's/Confidence/confidence/g' ./taxonomy/taxonomy.tsv

# Add taxonomy data to .biom file
biom add-metadata -i ./table/feature-table.biom -o ./table-with-taxonomy.biom --observation-metadata-fp ./taxonomy/taxonomy.tsv --sc-separated taxonomy

# Convert to json format
biom convert -i ./table-with-taxonomy.biom -o ./table-with-taxonomy-json2.biom --table-type="OTU table" --to-json
```

-------------------------------------------------------------------------

## Code for the metaanalysis

Important: DADA2 should be run independently for each batch, as its error connection is batch-specific.

**Import demultipliex sequences: do it for each batch of sequences**
```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path Manifest \
  --output-path ./demux-paired-end.qza \
  --input-format PairedEndFastqManifestPhred33
```

**Check parameters of the reads**
```bash
qiime demux summarize \
  --i-data ./demux-paired-end.qza \
  --o-visualization ./demux-paired-end.qzv
```

--------------------------------
**Batch 1 - Solar panels**
 ```bash
 qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./demux-paired-end.qza \
  --verbose \
  --p-trim-left-f 21 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 289 \
  --p-trunc-len-r 222 \
  --p-n-threads 8 \
  --o-table ./table.qza \
  --o-representative-sequences ./rep-seqs.qza \
  --o-denoising-stats ./denoising-stats.qza 

qiime metadata tabulate \
  --m-input-file ./denoising-stats.qza \
  --o-visualization ./stats-dada2.qzv
```

--------------------------------

**Batch 2 - Tank lid**
```bash
 qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./demux-paired-end.qza \
  --verbose \
  --p-trim-left-f 21 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 289 \
  --p-trunc-len-r 250 \
  --p-n-threads 8 \
  --o-table ./table.qza \
  --o-representative-sequences ./rep-seqs.qza \
  --o-denoising-stats ./denoising-stats.qza 

qiime metadata tabulate \
  --m-input-file ./denoising-stats.qza \
  --o-visualization ./stats-dada2.qzv
```

--------------------------------

**Batch 3 - CosmBisphPhen**
```bash
 qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./demux-paired-end.qza \
  --verbose \
  --p-trim-left-f 21 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 290 \
  --p-trunc-len-r 227 \
  --p-n-threads 8 \
  --o-table ./table.qza \
  --o-representative-sequences ./rep-seqs.qza \
  --o-denoising-stats ./denoising-stats.qza 

qiime metadata tabulate \
  --m-input-file ./denoising-stats.qza \
  --o-visualization ./stats-dada2.qzv
```

--------------------------------

**Batch 4 - DieselSpillCompost**
```bash
 qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./demux-paired-end.qza \
  --verbose \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 240 \
  --p-n-threads 8 \
  --o-table ./table.qza \
  --o-representative-sequences ./rep-seqs.qza \
  --o-denoising-stats ./denoising-stats.qza 

qiime metadata tabulate \
  --m-input-file ./denoising-stats.qza \
  --o-visualization ./stats-dada2.qzv
```

--------------------------------

**Batch 5 - SRR5457-**
```bash
 qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./demux-paired-end.qza \
  --verbose \
  --p-trim-left-f 21 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 295 \
  --p-trunc-len-r 225 \
  --p-n-threads 8 \
  --o-table ./table.qza \
  --o-representative-sequences ./rep-seqs.qza \
  --o-denoising-stats ./denoising-stats.qza 


qiime metadata tabulate \
  --m-input-file ./denoising-stats.qza \
  --o-visualization ./stats-dada2.qzv
```

--------------------------------

**Batch 6 - SRR10261600**
```bash
 qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./demux-paired-end.qza \
  --verbose \
  --p-trim-left-f 21 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 245 \
  --p-trunc-len-r 245 \
  --p-n-threads 8 \
  --o-table ./table.qza \
  --o-representative-sequences ./rep-seqs.qza \
  --o-denoising-stats ./denoising-stats.qza 

qiime metadata tabulate \
  --m-input-file ./denoising-stats.qza \
  --o-visualization ./stats-dada2.qzv
```

--------------------------------

**Batch 7 - SRR125276-**
```bash
 qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./demux-paired-end.qza \
  --verbose \
  --p-trim-left-f 21 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 290 \
  --p-trunc-len-r 240 \
  --p-n-threads 8 \
  --o-table ./table.qza \
  --o-representative-sequences ./rep-seqs.qza \
  --o-denoising-stats ./denoising-stats.qza 

qiime metadata tabulate \
  --m-input-file ./denoising-stats.qza \
  --o-visualization ./stats-dada2.qzv
```

--------------------------------

**Batch 8 - SRR99641-**
```bash
 qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./demux-paired-end.qza \
  --verbose \
  --p-trim-left-f 21 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 290 \
  --p-trunc-len-r 225 \
  --p-n-threads 8 \
  --o-table ./table.qza \
  --o-representative-sequences ./rep-seqs.qza \
  --o-denoising-stats ./denoising-stats.qza 

qiime metadata tabulate \
  --m-input-file ./denoising-stats.qza \
  --o-visualization ./stats-dada2.qzv
```

--------------------------------

**Batch 9 - SRR82590-**
```bash
 qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./demux-paired-end.qza \
  --verbose \
  --p-trim-left-f 21 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 277 \
  --p-trunc-len-r 219 \
  --p-n-threads 8 \
  --o-table ./table.qza \
  --o-representative-sequences ./rep-seqs.qza \
  --o-denoising-stats ./denoising-stats.qza 

qiime metadata tabulate \
  --m-input-file ./denoising-stats.qza \
  --o-visualization ./stats-dada2.qzv
```

--------------------------------

**Batch 10 - SRR15809-**
```bash
 qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./demux-paired-end.qza \
  --verbose \
  --p-trim-left-f 21 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-n-threads 8 \
  --o-table ./table.qza \
  --o-representative-sequences ./rep-seqs.qza \
  --o-denoising-stats ./denoising-stats.qza 

qiime metadata tabulate \
  --m-input-file ./denoising-stats.qza \
  --o-visualization ./stats-dada2.qzv
```


**Join all the objects following this tutorial: https://docs.qiime2.org/2020.8/tutorials/fmt/
1st: merge tables**
```bash
qiime feature-table merge \
  --i-tables ConsBiphPhen/table.qza \
  --i-tables DieselSpillCompost/table.qza \
  --i-tables Solar-Panels/table.qza \
  --i-tables SRR10261600/table.qza \
  --i-tables SRR125276-/table.qza \
  --i-tables SRR15809-/table.qza \
  --i-tables SRR5457-/table.qza \
  --i-tables SRR82590-/table.qza \
  --i-tables SRR99641-/table.qza \
  --i-tables Tank-lid/table.qza \
  --o-merged-table table-all.qza
```

**2nd: merge rep-seqs**

```bash
qiime feature-table merge-seqs \
  --i-data ConsBiphPhen/rep-seqs.qza \
  --i-data DieselSpillCompost/rep-seqs.qza \
  --i-data Solar-Panels/rep-seqs.qza \
  --i-data SRR10261600/rep-seqs.qza \
  --i-data SRR125276-/rep-seqs.qza \
  --i-data SRR15809-/rep-seqs.qza \
  --i-data SRR5457-/rep-seqs.qza \
  --i-data SRR82590-/rep-seqs.qza \
  --i-data SRR99641-/rep-seqs.qza \
  --i-data Tank-lid/rep-seqs.qza \
  --o-merged-data rep-seqs-all.qza
```

**3rd: statistics after merging**
```bash
qiime feature-table summarize \
  --i-table table-all.qza \
  --o-visualization table-all.qzv
```

**Taxonomic classification. Database SILVA 138**
```bash
qiime feature-classifier classify-sklearn \
  --i-classifier /media/darwin/c416f3d9-3986-47fd-b7ea-5f73aae3a518/Databases/silva-138-99-nb-classifier-2020-08.qza \
  --i-reads ./rep-seqs-all.qza \
  --verbose \
  --p-n-jobs 8 \
  --o-classification ./taxonomy.qza
```

```bash
# Export table.qza and taxonomy.qza
qiime tools export \
 --input-path ./table-all.qza \
  --output-path  ./table

# Export taxonomy.qza and taxonomy.qza
qiime tools export \
 --input-path ./taxonomy.qza \
  --output-path  ./taxonomy

# Modify the biom-taxonomy tsv headers: change header "Feature ID" to "#OTUID"; "Taxon" to "taxonomy"; and "Confidence" to "confidence"
sed -i -e 's/Feature ID/#OTUID/g' ./taxonomy/taxonomy.tsv
sed -i -e 's/Taxon/taxonomy/g' ./taxonomy/taxonomy.tsv
sed -i -e 's/Confidence/confidence/g' ./taxonomy/taxonomy.tsv

# Add taxonomy data to .biom file
biom add-metadata -i ./table/feature-table.biom -o ./table-with-taxonomy.biom --observation-metadata-fp ./taxonomy/taxonomy.tsv --sc-separated taxonomy

# Convert to json format
biom convert -i ./table-with-taxonomy.biom -o ./table-with-taxonomy-json2.biom --table-type="OTU table" --to-json
```

# R code
### Data import

Import the BIOM table

```r
library(phyloseq)
physeqObjest <- import_biom("path/to/table-with-taxonomy-json2.biom", treefilename = "path/to/tree.nwk",
                 refseqArgs=NULL, parseFunction = parse_taxonomy_default, 
                 parallel=FALSE)
physeqObjest
```

And then, import (and arrange) the metadata file.

```r
# METADATA
mapfile = read.csv("path/to/metadata.csv", header = TRUE, sep=",")
samples = mapfile[,1]
rownames(mapfile) = samples
mapfile = mapfile[,2:length(colnames(mapfile))]
# Hay que cambiar el orden de las filas para adecuarlo al orden de las muestras en la OTU table.
ord = match(colnames(otu_table(physeqObjest)), rownames(mapfile))
mapfile = mapfile[ord,]
sampledata = sample_data(mapfile)

physeq1 = merge_phyloseq(physeqObjest, sampledata)
physeq1
```

Remove sample "17G" as it has not enough sequences for a proper analysis

```{r}
physeq1 = subset_samples(physeq1, SampleName != "17G")
physeq1
sample_sums(physeq1)
```

