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

# R code - Main Analysis
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

```r
physeq1 = subset_samples(physeq1, SampleName != "17G")
physeq1
sample_sums(physeq1)
```

### Rarefaction Curves

**Diesel**

```r
library("iNEXT")
Rarefaction <- iNEXT(as.data.frame(otu_table(physeq1)[,c("PS1", "PS2", "PS3", "PS4", "PS5", "PS6", "PS7", "PS8", "PS9", "PS10")]), q=0, datatype="abundance")
Rareplot1 <- ggiNEXT(x=Rarefaction, type=1) + theme_bw() +
  geom_line(size=0.1) +
  geom_point(size=0, na.rm = TRUE) +
  scale_shape_manual(values=c(rep(20,55))) + xlab("Number of sequences") + ylab("Richness")
Rareplot1
```

**Gasoline**

```r
library("iNEXT")
Rarefaction <- iNEXT(as.data.frame(otu_table(physeq1)[,c("PS11", "PS12", "PS13", "PS14", "PS15", "PS16", "PS18", "PS19", "PS20")]), q=0, datatype="abundance")
Rareplot1 <- ggiNEXT(x=Rarefaction, type=1) + theme_bw() +
  geom_line(size=0.1) +
  geom_point(size=0, na.rm = TRUE) +
  scale_shape_manual(values=c(rep(20,55))) + xlab("Number of sequences") + ylab("Richness")
Rareplot1
```

**Evolved**

```r
library("iNEXT")
Rarefaction <- iNEXT(as.data.frame(otu_table(physeq1)[,c("P3-1", "P3-9", "P3-10")]), q=0, datatype="abundance")
Rareplot1 <- ggiNEXT(x=Rarefaction, type=1) + theme_bw() +
  geom_line(size=0.1) +
  geom_point(size=0, na.rm = TRUE) +
  scale_shape_manual(values=c(rep(20,55))) + xlab("Number of sequences") + ylab("Richness")
Rareplot1
```

### Beta diversity

PCoA with Bray-Curtis dissimilarities

```r
# Calculate distance matrix
brayDist <- phyloseq::distance(physeqfull_relabund, method="bray")
# Calculate ordination
iMDS  <- ordinate(physeqfull_relabund, distance=brayDist, method = "PCoA")
## Make plot
# Create plot, store as temp variable, p
p <- plot_ordination(physeqfull_relabund, iMDS, color ="Type", shape = "Evolved", label = "SampleName")
# Add title to each plot
p <- p + ggtitle("PCoA using Bray-Curtis dissimilarities") + guides(color=guide_legend(ncol=3))
p
```

Adonis test

```r
library(vegan)
physeq33 <- as(sample_data(physeqfull_relabund), "data.frame")
adonis2(brayDist ~ Evolved, data = physeq33)
```

**Now I'll exclude the evolved samples**

```r
no.evol = subset_samples(physeq_R6, Evolved == "No")
no.evol.rel  = transform_sample_counts(no.evol, function(x) x / sum(x)*100 )
sample_sums(no.evol)
sample_sums(no.evol.rel)
```

```r
library(ggrepel)
# Calculate distance matrix
brayDist <- phyloseq::distance(no.evol.rel, method="bray")
# Calculate ordination
iMDS  <- ordinate(no.evol.rel, distance=brayDist, method = "PCoA")
## Make plot
# Create plot, store as temp variable, p
p <- plot_ordination(no.evol.rel, iMDS, color ="Type")
# Add title to each plot
p <- p + ggtitle("PCoA using Bray-Curtis dissimilarities") + guides(color=guide_legend(ncol=3)) + geom_text_repel(aes(label = SampleName))
p
```

Adonis test

```r
library(vegan)
physeq33 <- as(sample_data(no.evol.rel), "data.frame")
adonis2(brayDist ~ Type, data = physeq33)
```

# Alpha diversity

**Excluding the evolved samples**

```r
library(ggplot2)
no.evol.otu = subset_samples(physeq1, Evolved == "No" & SampleName != "17G")
no.evol.otu = rarefy_even_depth(no.evol.otu, sample.size = min(sample_sums(no.evol.otu)), rngseed = 711)
p = plot_richness(no.evol.otu, x="Type", color="Type", measures=c("Observed","Simpson")) + xlab("") + ylab("")
p
```

**Including the evolved samples**

```{r}
library(ggplot2)
no.evol.otu = subset_samples(physeq1, Type == "Diesel")
no.evol.otu = rarefy_even_depth(no.evol.otu, sample.size = min(sample_sums(no.evol.otu)), rngseed = 711)
p = plot_richness(no.evol.otu, x="Evolved", color="Evolved", measures=c("Observed","Simpson")) + xlab("") + ylab("") + geom_boxplot()
p
```

### Taxonomic composition

**Top 10: taxa in evolved samples**

```r
selected = names(sort(taxa_sums(evol.rel), decreasing = TRUE)[c(1:10)])
selected = otu_table(evol.rel)[selected, ]
other = c(100, 100, 100) - sample_sums(selected)
other
selected.df = rbind(as.data.frame(selected), other)
selected.df
colSums(selected.df)
```

Name the taxa

```r
rownames(selected.df) = c(tax_table(evol.rel)[rownames(selected),"Rank6"], "Other")
selected.df
```

Conver into long format

```r
arb1 = data.frame(selected.df$`P3-1`, rownames(selected.df), rep("P3-1", length(rownames(selected.df))))
colnames(arb1) = c("Relative_Abundance", "Genus", "Sample")
arb2 = data.frame(selected.df$`P3-9`, rownames(selected.df), rep("P3-9", length(rownames(selected.df))))
colnames(arb2) = c("Relative_Abundance", "Genus", "Sample")
arb3 = data.frame(selected.df$`P3-10`, rownames(selected.df), rep("P3-10", length(rownames(selected.df))))
colnames(arb3) = c("Relative_Abundance", "Genus", "Sample")
long.selected.df = rbind(arb1, arb2, arb3)
long.selected.df
```

Do the plot

```r
library("RColorBrewer")
library(ggplot2)
palette_Dark2 <- colorRampPalette(brewer.pal(11, "Dark2"))
p = ggplot(data=long.selected.df, aes(x=Sample, y=Relative_Abundance, fill=Genus)) +
  geom_bar(stat="identity", color = "black") + scale_fill_manual(values = palette_Dark2(14)) + xlab("Muestras") + ylab("Abundancia relativa (%)")
p
```

**Top 20 genera excluding the evolved samples**

```r
library(phyloseq)
selected = names(sort(taxa_sums(no.evol.rel), decreasing = TRUE)[c(1:20)])
selected = otu_table(no.evol.rel)[selected, ]
other = rep(100, length(colnames(selected))) - sample_sums(selected)
other
selected.df = rbind(as.data.frame(selected), other)
selected.df
colSums(selected.df)
```

```r
rownames(selected.df) = c(tax_table(no.evol.rel)[rownames(selected),"Rank6"], "Other")
selected.df
```

```r
arb1 = data.frame(selected.df$`PS1`, rownames(selected.df), rep("PS1", length(rownames(selected.df))))
colnames(arb1) = c("Relative_Abundance", "Genus", "Sample")
arb2 = data.frame(selected.df$`PS2`, rownames(selected.df), rep("PS2", length(rownames(selected.df))))
colnames(arb2) = c("Relative_Abundance", "Genus", "Sample")
arb3 = data.frame(selected.df$`PS3`, rownames(selected.df), rep("PS3", length(rownames(selected.df))))
colnames(arb3) = c("Relative_Abundance", "Genus", "Sample")
arb4 = data.frame(selected.df$`PS4`, rownames(selected.df), rep("PS4", length(rownames(selected.df))))
colnames(arb4) = c("Relative_Abundance", "Genus", "Sample")
arb5 = data.frame(selected.df$`PS5`, rownames(selected.df), rep("PS5", length(rownames(selected.df))))
colnames(arb5) = c("Relative_Abundance", "Genus", "Sample")
arb6 = data.frame(selected.df$`PS6`, rownames(selected.df), rep("PS6", length(rownames(selected.df))))
colnames(arb6) = c("Relative_Abundance", "Genus", "Sample")
arb7 = data.frame(selected.df$`PS7`, rownames(selected.df), rep("PS7", length(rownames(selected.df))))
colnames(arb7) = c("Relative_Abundance", "Genus", "Sample")
arb9 = data.frame(selected.df$`PS9`, rownames(selected.df), rep("PS9", length(rownames(selected.df))))
colnames(arb9) = c("Relative_Abundance", "Genus", "Sample")
arb10 = data.frame(selected.df$`PS10`, rownames(selected.df), rep("PS10", length(rownames(selected.df))))
colnames(arb10) = c("Relative_Abundance", "Genus", "Sample")
arb11 = data.frame(selected.df$`PS11`, rownames(selected.df), rep("PS11", length(rownames(selected.df))))
colnames(arb11) = c("Relative_Abundance", "Genus", "Sample")
arb12 = data.frame(selected.df$`PS12`, rownames(selected.df), rep("PS12", length(rownames(selected.df))))
colnames(arb12) = c("Relative_Abundance", "Genus", "Sample")
arb13 = data.frame(selected.df$`PS13`, rownames(selected.df), rep("PS13", length(rownames(selected.df))))
colnames(arb13) = c("Relative_Abundance", "Genus", "Sample")
arb14 = data.frame(selected.df$`PS14`, rownames(selected.df), rep("PS14", length(rownames(selected.df))))
colnames(arb14) = c("Relative_Abundance", "Genus", "Sample")
arb15 = data.frame(selected.df$`PS15`, rownames(selected.df), rep("PS15", length(rownames(selected.df))))
colnames(arb15) = c("Relative_Abundance", "Genus", "Sample")
arb16 = data.frame(selected.df$`PS16`, rownames(selected.df), rep("PS16", length(rownames(selected.df))))
colnames(arb16) = c("Relative_Abundance", "Genus", "Sample")
arb18 = data.frame(selected.df$`PS18`, rownames(selected.df), rep("PS18", length(rownames(selected.df))))
colnames(arb18) = c("Relative_Abundance", "Genus", "Sample")
arb19 = data.frame(selected.df$`PS19`, rownames(selected.df), rep("PS19", length(rownames(selected.df))))
colnames(arb19) = c("Relative_Abundance", "Genus", "Sample")
arb20 = data.frame(selected.df$`PS20`, rownames(selected.df), rep("PS20", length(rownames(selected.df))))
colnames(arb20) = c("Relative_Abundance", "Genus", "Sample")
arb8 = data.frame(selected.df$`PS8`, rownames(selected.df), rep("PS8", length(rownames(selected.df))))
colnames(arb8) = c("Relative_Abundance", "Genus", "Sample")
long.selected.df = rbind(arb1, arb2, arb3, arb4, arb5, arb6, arb7, arb8, arb9, arb10, arb11, arb12, arb13, arb14, arb15, arb16, arb18, arb19, arb20)
long.selected.df
```

```r
library("RColorBrewer")
library(ggplot2)
colourCount = 21
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
p = ggplot(data=long.selected.df, aes(x=Sample, y=Relative_Abundance, fill=Genus)) +
  geom_bar(stat="identity", color = "black") + scale_fill_manual(values = getPalette(colourCount)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Muestras") + ylab("Abundancia relativa (%)")
p
```

**Top 20 phyla excluding the evolved samples**

```r
physeq_R2 <- tax_glom(physeq1, taxrank = rank_names(physeq1)[2])
physeqfull_relabund  = transform_sample_counts(physeq_R2, function(x) x / sum(x)*100 )
```

```r
library(phyloseq)
selected = names(sort(taxa_sums(physeqfull_relabund), decreasing = TRUE))
selected = otu_table(physeqfull_relabund)[selected, ]
other = rep(100, length(colnames(selected))) - sample_sums(selected)
other
selected.df = rbind(as.data.frame(selected), other)
selected.df
colSums(selected.df)
```

```r
rownames(selected.df) = c(tax_table(physeqfull_relabund)[rownames(selected),"Rank2"], "Other")
selected.df
```

```r
arbP1 = data.frame(selected.df$`P3-1`, rownames(selected.df), rep("P3-1", length(rownames(selected.df))))
colnames(arbP1) = c("Relative_Abundance", "Genus", "Sample")
arbP9 = data.frame(selected.df$`P3-9`, rownames(selected.df), rep("P3-9", length(rownames(selected.df))))
colnames(arbP9) = c("Relative_Abundance", "Genus", "Sample")
arbP10 = data.frame(selected.df$`P3-10`, rownames(selected.df), rep("P3-10", length(rownames(selected.df))))
colnames(arbP10) = c("Relative_Abundance", "Genus", "Sample")
arb1 = data.frame(selected.df$`PS1`, rownames(selected.df), rep("PS1", length(rownames(selected.df))))
colnames(arb1) = c("Relative_Abundance", "Genus", "Sample")
arb2 = data.frame(selected.df$`PS2`, rownames(selected.df), rep("PS2", length(rownames(selected.df))))
colnames(arb2) = c("Relative_Abundance", "Genus", "Sample")
arb3 = data.frame(selected.df$`PS3`, rownames(selected.df), rep("PS3", length(rownames(selected.df))))
colnames(arb3) = c("Relative_Abundance", "Genus", "Sample")
arb4 = data.frame(selected.df$`PS4`, rownames(selected.df), rep("PS4", length(rownames(selected.df))))
colnames(arb4) = c("Relative_Abundance", "Genus", "Sample")
arb5 = data.frame(selected.df$`PS5`, rownames(selected.df), rep("PS5", length(rownames(selected.df))))
colnames(arb5) = c("Relative_Abundance", "Genus", "Sample")
arb6 = data.frame(selected.df$`PS6`, rownames(selected.df), rep("PS6", length(rownames(selected.df))))
colnames(arb6) = c("Relative_Abundance", "Genus", "Sample")
arb7 = data.frame(selected.df$`PS7`, rownames(selected.df), rep("PS7", length(rownames(selected.df))))
colnames(arb7) = c("Relative_Abundance", "Genus", "Sample")
arb9 = data.frame(selected.df$`PS9`, rownames(selected.df), rep("PS9", length(rownames(selected.df))))
colnames(arb9) = c("Relative_Abundance", "Genus", "Sample")
arb10 = data.frame(selected.df$`PS10`, rownames(selected.df), rep("PS10", length(rownames(selected.df))))
colnames(arb10) = c("Relative_Abundance", "Genus", "Sample")
arb11 = data.frame(selected.df$`PS11`, rownames(selected.df), rep("PS11", length(rownames(selected.df))))
colnames(arb11) = c("Relative_Abundance", "Genus", "Sample")
arb12 = data.frame(selected.df$`PS12`, rownames(selected.df), rep("PS12", length(rownames(selected.df))))
colnames(arb12) = c("Relative_Abundance", "Genus", "Sample")
arb13 = data.frame(selected.df$`PS13`, rownames(selected.df), rep("PS13", length(rownames(selected.df))))
colnames(arb13) = c("Relative_Abundance", "Genus", "Sample")
arb14 = data.frame(selected.df$`PS14`, rownames(selected.df), rep("PS14", length(rownames(selected.df))))
colnames(arb14) = c("Relative_Abundance", "Genus", "Sample")
arb15 = data.frame(selected.df$`PS15`, rownames(selected.df), rep("PS15", length(rownames(selected.df))))
colnames(arb15) = c("Relative_Abundance", "Genus", "Sample")
arb16 = data.frame(selected.df$`PS16`, rownames(selected.df), rep("PS16", length(rownames(selected.df))))
colnames(arb16) = c("Relative_Abundance", "Genus", "Sample")
arb18 = data.frame(selected.df$`PS18`, rownames(selected.df), rep("PS18", length(rownames(selected.df))))
colnames(arb18) = c("Relative_Abundance", "Genus", "Sample")
arb19 = data.frame(selected.df$`PS19`, rownames(selected.df), rep("PS19", length(rownames(selected.df))))
colnames(arb19) = c("Relative_Abundance", "Genus", "Sample")
arb20 = data.frame(selected.df$`PS20`, rownames(selected.df), rep("PS20", length(rownames(selected.df))))
colnames(arb20) = c("Relative_Abundance", "Genus", "Sample")
arb8 = data.frame(selected.df$`PS8`, rownames(selected.df), rep("PS8", length(rownames(selected.df))))
colnames(arb8) = c("Relative_Abundance", "Genus", "Sample")
long.selected.df = rbind(arbP1, arbP9, arbP10, arb1, arb2, arb3, arb4, arb5, arb6, arb7, arb8, arb9, arb10, arb11, arb12, arb13, arb14, arb15, arb16, arb18, arb19, arb20)
long.selected.df
```

```r
library("RColorBrewer")
library(ggplot2)
colourCount = 23
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
p = ggplot(data=long.selected.df, aes(x=Sample, y=Relative_Abundance, fill=Genus)) +
  geom_bar(stat="identity", color = "black") + scale_fill_manual(values = getPalette(colourCount)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Muestras") + ylab("Abundancia relativa (%)") + guides(fill=guide_legend(title="Phylum"))
p
```

**Top 20 taxa (class level) excluding the evolved samples**

Collapse the object

```{r}
physeq_R3 <- tax_glom(physeq1, taxrank = rank_names(physeq1)[3], NArm = FALSE)
physeqfull_relabund_R3  = transform_sample_counts(physeq_R3, function(x) x / sum(x)*100)
sample_sums(physeq_R3)
```

Write the table into a CSV file

```{r}
otu.table = data.frame(otu_table(physeqfull_relabund_R3))
tax.table = data.frame(tax_table(physeqfull_relabund_R3))
all = cbind(tax.table, otu.table)
write.table(all, row.names = FALSE, file = "./rel_otu_table_L3.csv", quote = FALSE, sep = '\t', dec = '.')
```

Modify the table in Excel and select only the top20 taxa. Add a column with the order (1st to 20th).

Load the modified table

```{r}
top20 = read.csv2("./top_L3.csv", sep = "\t", dec = ",")
colnames(top20)[c(3:length(colnames(top20)))] = c("1D","2D","3D","4D", "5D", "6D", "7D","8D", "9D","10D"," ", "11G","12G","13G","14G","15G","16G","18G","19G","20G")
top20
```

[Change to long format](http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/).

```{r results=FALSE}
library(tidyr)
# The arguments to gather():
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)
long.top20 <- gather(top20, Condition, Relative_abundance, "1D":"20G", factor_key=TRUE)
long.top20
```

Now do the plot

```{r}
library(ggplot2)
library(forcats)
library(RColorBrewer)
colourCount = 20

colors = c(colorRampPalette(brewer.pal(12, "Paired"))(colourCount))
colors[17] = "#d4bbdc"

p = ggplot(long.top20, aes(fill = fct_reorder(Class, Order, .desc = TRUE), y=Relative_abundance, x=Condition)) + 
    geom_bar( stat="identity", color = "black")
p = p + ylab("Relative abundance (%)") + xlab("") + guides(fill=guide_legend(ncol=2, title = "Taxa"))
p = p + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_rect(fill = "white", colour = "grey"))
p = p + scale_fill_manual(values = rev(colors))
p
```

# R code - Metaanalysis

Import the BIOM table

```r
library(phyloseq)
physeqObjest <- import_biom("path/to/table-with-taxonomy-json2.biom",
                            refseqArgs=NULL, parseFunction = parse_taxonomy_default, 
                            parallel=FALSE)
physeqObjest
```

Import (and format) the metadata file

```r
mapfile = read.csv("./metadata.csv", header = TRUE, sep=",")
samples = mapfile[,1]
rownames(mapfile) = samples
ord = match(colnames(otu_table(physeqObjest)), samples)
mapfile = mapfile[ord,]
sampledata = sample_data(mapfile)
physeq1 = merge_phyloseq(physeqObjest, sampledata)
physeq1
```

Remove chloroplasts, eukarya and mitochondria

```r
physeq_R6 <- tax_glom(physeq1, taxrank = rank_names(physeq1)[6], NArm = FALSE)
physeq_R6_rel  = transform_sample_counts(physeq_R6, function(x) x / sum(x)*100 )
sample_sums(physeq_R6_rel)
physeq_R6_rel
```

# Alfa diversity

Rarefy to the 2nd sample with less sequences (32,070).

**ASV level**

```r
library(ggplot2)
library(RColorBrewer)
p = plot_richness(rarefy_even_depth(physeq1, sample.size=32070, rngseed = 711), x="SampleName",measures=c("Observed", "Shannon","Simpson")) + xlab("") + ylab("") + geom_boxplot(alpha = .4, position=position_dodge(1)) + theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_brewer(palette = "Set2")
p
```

**Genus level**

```r
library(ggplot2)
library(RColorBrewer)
p = plot_richness(rarefy_even_depth(physeq_R6,sample.size = 32070,  rngseed = 711), x="SampleName",measures=c("Observed", "Shannon","Simpson")) + xlab("") + ylab("") + geom_boxplot(alpha = .4, position=position_dodge(1)) + theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_brewer(palette = "Set2")
p
```

**Heatmap. Top 20 Phyla**

```r
library(ampvis2)
library(ggplot2)
if(!require("devtools"))
  install.packages("devtools")
#source the phyloseq_to_ampvis2() function from the gist
devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")

# The colnames have to be changed
colnames(tax_table(physeq_R2_rel)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#go converty
ampvis2_obj <- phyloseq_to_ampvis2(physeq_R2_rel)
ampvis2_obj
# Change the metadata file just to include the desired names
rownames(ampvis2_obj$metadata) = ampvis2_obj$metadata[,"SampleName"]
colnames(ampvis2_obj$abund) = ampvis2_obj$metadata[,"SampleName"]
```

```r
p = amp_heatmap(
      data = ampvis2_obj,
      facet_by = "SampleName",
      normalise = FALSE,
      tax_show = 20,
      tax_aggregate = "Phylum",
      plot_values_size = 3,
      min_abundance = 0.01,
      color_vector = c("royalblue3",
                   "whitesmoke",
                   "lightcoral"),
      round = 2
    )

p
```

**Heatmap. Top 20 Genera**

```r
library(ampvis2)
library(ggplot2)
if(!require("devtools"))
  install.packages("devtools")
#source the phyloseq_to_ampvis2() function from the gist
devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")

# Necesitamos cambiar los colnames de la taxtable:
colnames(tax_table(physeq_R6_rel)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#go converty
ampvis2_obj <- phyloseq_to_ampvis2(physeq_R6_rel)
ampvis2_obj
# Cambiamos el metadata para que la figura salga con el orden y los nombres que queremos
rownames(ampvis2_obj$metadata) = ampvis2_obj$metadata[,"SampleName"]
colnames(ampvis2_obj$abund) = ampvis2_obj$metadata[,"SampleName"]
```

```r
p = amp_heatmap(
      data = ampvis2_obj,
      facet_by = "SampleName",
      normalise = FALSE,
      tax_show = 20,
      tax_aggregate = "Genus",
      tax_add = "Class",
      plot_values_size = 3.4,
      min_abundance = 0.01,
      color_vector = c("royalblue3",
                   "whitesmoke",
                   "lightcoral"),
      round = 2
    )
p
```


**IMPORTANT**: I've used the [phyloseq_to_df](https://rdrr.io/github/vmikk/metagMisc/man/phyloseq_to_df.html) function by vmikk

