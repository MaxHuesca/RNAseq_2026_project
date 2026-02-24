# Differential expression analysis in mice cortex tissue with upregulation of TREM2 
> De Los Santos Huesca Ismael Maximiliano
>
> A DE analysis using a dataset from recount3

## Summary 
[recount3](https://rna.recount.bio/) is an online resource consisting of RNA-seq gene, exon, and exon-exon junction counts as well as coverage bigWig files for 8,679 and 10,088 different studies for human and mouse respectively. It is the third generation of the ReCount project and part of recount.bio. With this database of uniformly processed data from RNA-seq, we can account for the differential expression analysis. 

This analysis focused on the use of [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment) objects to account for the preprocessing of the data until we can use the information stored in the 3 principal features of this type of object: the RowData table, the assay or assays table, and the colData table.

## Contents 
This repository has 4 main directories to store all the source and results of the analysis: 
* [`Results`](https://github.com/MaxHuesca/RNAseq_2026_project/tree/devel/Results): Directory that stores all the results that are produced by the entire DE analysis
  * [`objects_R`](https://github.com/MaxHuesca/RNAseq_2026_project/tree/devel/Results/objects_R): Here the code stores the SummarizedExperiment objects that the source uses across the different steps in the analysis
  * [`plots`](https://github.com/MaxHuesca/RNAseq_2026_project/tree/devel/Results/plots): Here are stored the visual results such as plots (boxplots, volcano plots, histograms, etc...)
   
* [`Code`](https://github.com/MaxHuesca/RNAseq_2026_project/tree/devel/code): This directory has all the code sources necessary for the DE analysis. It has two principal types of scripts: an .Rmd explaining the code and interpretation of the ongoing steps across the analysis. The entire pipeline could go ahead with the `Entire_analysis.Rmd` script but also each step has its own directory.
  * [`Explore_charge_data`](https://github.com/MaxHuesca/RNAseq_2026_project/tree/devel/code/Explore_charge_data): This is the data loading step. Here the SummarizedExperiment is created and some features like the *sra_attributes* are explored.
  * [`Vizualization_construction`](https://github.com/MaxHuesca/RNAseq_2026_project/tree/devel/code/Vizualization_construction): Here the preprocessing of the data and exploration of the variables we previously expanded are made, as well as the filtering and normalization process.
  * [`DE_analysis`](https://github.com/MaxHuesca/RNAseq_2026_project/tree/devel/code/DE_analysis): This is the final step where the statistical model is made and also the modeling of the data variability using the `limma` package.
* [`docs`](https://github.com/MaxHuesca/RNAseq_2026_project/tree/devel/docs): Here we store the information and binary files that are related to results and previous work with the data.
* [`raw-data`](https://github.com/MaxHuesca/RNAseq_2026_project/tree/devel/raw-data): This directory is intended for storing raw counts files like SRR of the experiment procedure. 

```bash
📦 Proyecto
├── 📁 Results
│   ├── 📁 objects_R
│   │   └── README.md
│   └── 📁 plots
│       ├── 📁 DE_analysis
│       └── 📁 proccesed_data  
│
├── 📁 code
│   ├── 📁 DE_analysis
│   │   ├── DE_anal.R
│   │   └── DE_anal.Rmd
│   ├── Entire_analysis.Rmd
│   ├── 📁 Explore_charge_data
│   │   ├── Explore_&_read.R
│   │   ├── Explore_&_read.Rmd
│   │   └── 📁 logs
│   ├── 📁 Vizualization_construction
│   │   ├── DataPreprocess.R
│   │   ├── PrePros_explore_contruc.Rmd
│   │   └── 📁 logs
│
├── 📁 docs
│   └── Art_from_data.pdf
│
└── 📁 raw-data
    ├── 📁 FASTQ
    │   └── README.md
    └── 📁 sample_info
        └── README.md
```

### Libraries used 
- [recount3](https://bioconductor.org/packages/recount3)
- [edgeR](https://bioconductor.org/packages/edgeR)
- [ggplot2](https://ggplot2.tidyverse.org)
- [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment)
- [iSEE](https://bioconductor.org/packages/iSEE)
- [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment)
- [scater](https://bioconductor.org/packages/scater)
- [dplyr](https://dplyr.tidyverse.org)
- [variancePartition](https://bioconductor.org/packages/variancePartition)
- [limma](https://bioconductor.org/packages/limma)
- [ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap)
- [ggrepel](https://ggrepel.slowkow.com)

## Results  

### SRP119675 project 
>This step relies on the script [Explore_&_read.Rmd](https://github.com/MaxHuesca/RNAseq_2026_project/blob/devel/code/Explore_charge_data/Explore_%26_read.Rmd)

The `SRP119675` is associated with the study of AD with the title "TREM2 Gene Dosage Increase Reprograms Microglia Responsivity and Ameliorates Pathological Phenotypes in Alzheimer's Disease Models", in which the authors explore the role of TREM2 upregulation in a typical AD model mouse 5XFAD. The object has `55500` genes within `144` samples. It also stores two different assay tables: one with the raw counts and another one with these normalized by recount3. The samples from this project have a total of 6 attributes, all are categorical so they can be managed as factors in R: 

- sra_attribute.age: 2m, 4m, 7m, are the age in months of the mice from which the samples were extracted. The 3 categories have 48 samples each, accounting for the total 144.

- sra_attribute.brain_region: In all the samples analyzed, the brain region where it was taken is the cortex.

- sra_attribute.gender: A two-level variable with the sex of the mouse, 72 samples for male and the other 72 for female.

- sra_attribute.genotype: A categorical variable that specifies the genetic background of the mouse. In this study there are 4 different genotypes with the number of samples representing them: 
  * 5xFAD: 36 
  * 5xFAD/BAC-TREM2: 36
  * BAC-TREM2: 36
  * WT: 36 

- sra_attribute.source_name: Character variable that stores the number and name of the mouse.

- sra_attribute.strain: Strain of the mouse; all mice belong to (C57BL/6J x FvB/NJ) F1.

- sra_attribute.tissue: The tissue where the sample was extracted: brain.

### Filtering the data
>This step relies on the script [PrePros_explore_contruc.Rmd](https://github.com/MaxHuesca/RNAseq_2026_project/blob/devel/code/Vizualization_construction/PrePros_explore_contruc.Rmd)

Once we have the data loaded we can start preprocessing it. For this task we have two major tasks to do. Here we describe the filter by gene proportion (the filter by gene expression is modeled by the edgeR function `FilterbyExp()`; you can review it in the script that belongs to this section), in which we explore the representation of the genes across the samples in order to find those that have an irregular behavior. To account for that we made two different plots: a histogram and a dot plot of this proportion across samples. 

| **Histogram of gene proportion** | **Dot plot for gene proportion** |
|:---:|:---:|
| ![Histogram](https://github.com/MaxHuesca/RNAseq_2026_project/blob/devel/Results/plots/proccesed_data/Histogram_distribution_geneProp.png) | ![Dot plot](https://github.com/MaxHuesca/RNAseq_2026_project/blob/devel/Results/plots/proccesed_data/Ponit_distribution_geneProp.png) |

### Dimensionality reduction analysis 
>This step relies on the script [PrePros_explore_contruc.Rmd](https://github.com/MaxHuesca/RNAseq_2026_project/blob/devel/code/Vizualization_construction/PrePros_explore_contruc.Rmd)

Once we have filtered the data and normalized it, we can start exploring the distribution of our samples according to the attributes they have. For this task a widely used technique is dimensionality reduction, so we use two principal techniques: Multidimensional Scaling (MDS) and Principal Component Analysis (PCA). 

#### PCA 
We can do a Principal Component Analysis (PCA) for viewing how our samples group based on the variation that the principal components (PCs) can explain. The objective of doing this is to relate a variable (qualitative, discrete variable) to the clustering of the samples. 

![PCA](https://github.com/MaxHuesca/RNAseq_2026_project/blob/devel/Results/plots/proccesed_data/preprocesed_PCA.png) 

#### MDS
As we did with PCA, Multidimensional Scaling is another technique widely used in large datasets for accounting similarity/dissimilarity between objects in a low-dimensional space, using distances in a dissimilarity matrix that contrasts with PCA, which uses statistical measures. Instead, it uses distances like Euclidean distance. 

![MDS](https://github.com/MaxHuesca/RNAseq_2026_project/blob/devel/Results/plots/proccesed_data/preprocessed_MDS.png) 

As we can note, neither PCA nor MDS could explain the variance in the expression levels by the genotype variable prior to the DE analysis, so we can do a secondary analysis with `variancePartition` to address the task of identifying the most relevant variables in our data; in other words, those variables that explain a high level of the variance in our data.  

### variancePartition 
>This step relies on the script [PrePros_explore_contruc.Rmd](https://github.com/MaxHuesca/RNAseq_2026_project/blob/devel/code/Vizualization_construction/PrePros_explore_contruc.Rmd) 

The [variancePartition](https://bioconductor.org/packages//release/bioc/vignettes/variancePartition/inst/doc/variancePartition.html) package implements a statistical method to quantify the contribution of multiple sources of variation and decouple between our samples. So we can adjust the mixed linear model that this library uses to fit our data in order to explain our variance. This can be visualized by a violin plot. 

![Variance_plot_preprocessed](https://github.com/MaxHuesca/RNAseq_2026_project/blob/devel/Results/plots/proccesed_data/Variance_plot_preprocessed.png) 

We can see that the most variance in our dataset is explained by the individual variance across the samples, followed by the variance of the age and genotype. Something important to see in these plots is that some genes' variance is explained by the age as well as the genotype. Those genes could not represent the total variance across samples, and that was because we could not separate the clusters in the dimensionality reduction analysis like PCA or MDS.

With this data preprocessing coupled with the variance analyses of the genes across the available variables, we can have an idea of what we should consider in our model. In this case we can consider the genotype, age, and gender as those that explain the most variance and are not correlated with each other. Particularly, we can note that there are some genes whose variance is more explained by the genotype and age. We are focusing on these variables to search for those genes that could be the difference among our samples. 

### DE analysis 

With the previous steps of a DE analysis, now we can narrow down the formal DE analysis. For this we are using the [limma](https://bioconductor.org/packages//release/bioc/html/limma.html) package. This package adjusts linear models to the genes to assess linear regressions about their expression and then adjusts those values to model the differentially expressed genes. Because it accounts for large amounts of testing, it retrieves a padj value associated with the log2 fold changes.  

### Creating a model  
As we could see in the preprocessing part of the analysis, the variables that are relevant for explaining variance in some genes that would be of interest are: 
- sra_attribute.genotype --> as the principal variable 
- sra_attribute.age 
- sra_attribute.gender

In this case we are using an interaction model in which we consider the gender variable for the variance in our model, but also we are assuming that there could be differences in the expression of the genotypes between different ages.

$$
Y=\beta_0​+\beta_G​+\beta_T​+\beta_G:T​+\beta_S​+\epsilon
$$

With this defined model we have these different coefficients: 
- "sra_attribute.genotype5xFAD/BAC-TREM2"
- "sra_attribute.genotypeBAC-TREM2"
- "sra_attribute.genotypeWT"
- "sra_attribute.age4m"
- "sra_attribute.age7m"
- "sra_attribute.gendermale"
- "sra_attribute.genotype5xFAD/BAC-TREM2:sra_attribute.age4m"
- "sra_attribute.genotypeBAC-TREM2:sra_attribute.age4m"
- "sra_attribute.genotypeWT:sra_attribute.age4m"
- "sra_attribute.genotype5xFAD/BAC-TREM2:sra_attribute.age7m"
- "sra_attribute.genotypeBAC-TREM2:sra_attribute.age7m"
- "sra_attribute.genotypeWT:sra_attribute.age7m"

The reference levels or baseline we are using are:

For the genotype, the reference or baseline is:
- 5xFAD

For the age:
- 2m

For the gender:
- female 

After fitting the linear model that `limma` uses, we can obtain the topTables with the solved coefficients for the model we defined. In this analysis, because of how we defined it, we are going to focus on the differences between the genotypes 5xFAD - 5xFAD/BAC-TREM2 and 5xFAD - BAC-TREM2 across time ("sra_attribute.genotype5xFAD/BAC-TREM2:sra_attribute.age4m" - "sra_attribute.genotype5xFAD/BAC-TREM2:sra_attribute.age7m", "sra_attribute.genotypeBAC-TREM2:sra_attribute.age4m" - "sra_attribute.genotypeBAC-TREM2:sra_attribute.age7m" coefficients). So to visualize the results we can use volcano plots and heatmaps. 

#### Volcano Plots 

| **5xFAD - 5xFAD/BAC-TREM2 : 4m** | **5xFAD - 5xFAD/BAC-TREM2 : 7m** |
|:---:|:---:|
| ![5xFAD - 5xFAD/BAC-TREM2 : 4m](https://github.com/MaxHuesca/RNAseq_2026_project/blob/devel/Results/plots/DE_analysis/genotype5xFAD.BAC.TREM2.sra_attribute.age4m_vs_5XFAD_vulcano.png) | ![5xFAD - 5xFAD/BAC-TREM2 : 7m](https://github.com/MaxHuesca/RNAseq_2026_project/blob/devel/Results/plots/DE_analysis/genotype5xFAD.BAC.TREM2.sra_attribute.age7m_vs_5XFAD_vulcano.png) |

| **5xFAD - BAC-TREM2 : 4m** | **5xFAD - BAC-TREM2 : 7m** |
|:---:|:---:|
| ![5xFAD - BAC-TREM2 : 4m](https://github.com/MaxHuesca/RNAseq_2026_project/blob/devel/Results/plots/DE_analysis/genotypeBAC.TREM2.sra_attribute.age4m_vs_5XFAD_vulcano.png) | ![5xFAD - BAC-TREM2 : 7m](https://github.com/MaxHuesca/RNAseq_2026_project/blob/devel/Results/plots/DE_analysis/genotypeBAC.TREM2.sra_attribute.age7m_vs_5XFAD_vulcano.png) |

We can see that in both comparison tehere is a time dependent changue in the expression of diferent genes, and this genes that are DE in the both phenotypes with alzheimer could have a rol in the preservation of cognitive functions that the autors report in the main article. 

#### Heat Maps 

| **5xFAD - BAC-TREM2** | **5xFAD - 5xFAD/BAC-TREM2** |
|:---:|:---:|
| ![5xFAD - BAC-TREM2 heatmap](https://github.com/MaxHuesca/RNAseq_2026_project/blob/devel/Results/plots/DE_analysis/5XFAD_vs_TREM2%20heat.png) | ![5xFAD - 5xFAD/BAC-TREM2 heatmap](https://github.com/MaxHuesca/RNAseq_2026_project/blob/devel/Results/plots/DE_analysis/5XFAD_vs_5XFAD-TREM2heat.png) | 

With this data we can know that there is a real effect over time in expression across both comparisons between genotypes. We can see that the ones that have both Alzheimer's genes have the most dramatic log fold-change. Enrichment analysis would be good for drawing out the biological information about the role of TREM2 in AD, particularly in the 5XFAD model.

## Discussion 
`recount3` is a very useful tool for accounting for the first steps when you are narrowing down a DE analysis. It is very useful to have that easy way to load the data in a `SummarizedExperiment` object with all the information. With that, you can start visualizing your data very fast after normalization and some intermediate steps to get an idea of how your data is distributed and its variance across the samples. 

Something very important is to make sure of what to consider in the model for differential expression, and what not to. For these, tools like `iSEE` or even `variancePartition` are very useful. 

Once your model is finally defined, you can proceed in a near straightforward way to model your variance using DE tools like `limma`, `edgeR`, `DeSeq2`, or even `NOIseq`. 

Last but not least, plotting the results is another important way to make real what you found in your analysis. This part is important to make sense of all the statistics the model retrieves for you. 

In this particular analysis, we can note that the variance across the genes was explained by the high variance across the samples, probably the reason why our samples do not have clusters among our phenotypes. So in that, we can search for the variables that contribute more to our data with tools like `variancePartition` to give us an idea of what to consider in our model. 

With the model well defined, we can narrow down the DE analysis to get into the final step related to the results interpretation. In this part we can use volcano and heatmap plots to classify the genes related to our model and its coefficient comparisons that are upregulated or downregulated with a metric for the confidence, represented as the adjusted p-value (adjusted because of the multiple testing effect).

### References  

> This project follows the workflow developed by Dr. Leonardo Collado-Torres and uses the recount3 resource for gene-level RNA-seq summaries.

#### Software and Resources

>Collado-Torres L (2023). Explore and download data from the recount3 project.Bioconductor package recount3, version 1.8.0. https://doi.org/10.18129/B9.bioc.recount3 https://github.com/LieberInstitute/recount3

>A Bioconductor-style differential expression analysis powered by SPEAQeasy. (s. f.). https://research.libd.org/SPEAQeasyWorkshop2023/articles/SPEAQeasyWorkshop2023.html#differential-expression-analysis-1

#### Original Study

> Lee, C. D., Daggett, A., Gu, X., Jiang, L., Langfelder, P., Li, X., Wang, N., Zhao, Y., Park, C. S., Cooper, Y., Ferando, I., Mody, I., Coppola, G., Xu, H., & Yang, X. W. (2018). Elevated TREM2 Gene Dosage Reprograms Microglia Responsivity and Ameliorates Pathological Phenotypes in Alzheimer’s Disease Models. Neuron, 97(5), 1032-1048.e5. https://doi.org/10.1016/j.neuron.2018.02.002
