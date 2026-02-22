## Load the requiered R packages and the previous proccesed objects
library(recount3)
library(edgeR)
library(ggplot2)
library(SummarizedExperiment)
library(iSEE)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(variancePartition)

#===============================================================================
#Pre loading and functions
#===============================================================================
load("Results/objects_R/rse_gene_SRP119675_parsed.Rdata")

outdir = "Results"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

out_dir_plots <- ("Results/plots/proccesed_data")

if (!dir.exists(out_dir_plots)) {
  dir.create(out_dir_plots, recursive = TRUE, showWarnings = FALSE)
}


#function for saving plots
save_plot <- function(
  plot_obj,
  filename,
  width = 8,
  height = 6,
  dpi = 300,
  path = "Results/plots/proccesed_data"
) {
  ggplot2::ggsave(
    filename = file.path(path, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
}

#===============================================================================
#Filtering the data by gene proportion and levels of expression
#===============================================================================

#first we filter the reads with poor rerpresemation acroos samples calculating the gene proportions
rse_gene_SRP119675$assigned_gene_prop <- rse_gene_SRP119675$recount_qc.gene_fc_count_all.assigned /
  rse_gene_SRP119675$recount_qc.gene_fc_count_all.total

rse_gene_SRP119675_filtered <- rse_gene_SRP119675[,
  rse_gene_SRP119675$assigned_gene_prop > 0.810
]
# we can visualized that distribution of outliers with this plots
df <- as.data.frame(colData(rse_gene_SRP119675))

scater_plot <- ggplot(
  df,
  aes(x = reorder(rownames(df), assigned_gene_prop), y = assigned_gene_prop)
) +
  geom_point(alpha = 0.6) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(
    title = "Assigned Gene Proportion per Sample",
    x = "Samples (ordered)",
    y = "Assigned Gene Proportion"
  )

histogram_plot <- ggplot(df, aes(x = assigned_gene_prop)) +
  geom_histogram(
    bins = 30,
    fill = "steelblue",
    color = "white",
    alpha = 0.8
  ) +
  labs(
    title = "Distribution of Assigned Gene Proportion",
    x = "Assigned Gene Proportion",
    y = "Frequency"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

#also we can filter the genes that do not contributes to our variable distribution of interest
keep <- filterByExpr(
  assay(rse_gene_SRP119675_filtered),
  group = colData(rse_gene_SRP119675_filtered)$sra_attribute.genotype
)
rse_gene_SRP119675_filtered <- rse_gene_SRP119675_filtered[keep, ]

#===============================================================================
#Normalization
#===============================================================================

#Also is important to normalize the data in order to have comparable features acroos samples
dge <- DGEList(
  counts = assay(rse_gene_SRP119675_filtered, "counts"),
  genes = rowData(rse_gene_SRP119675_filtered)
)
dge <- calcNormFactors(dge)

# normalization for dimentional and variance analysis
logCPM <- edgeR::cpm(dge, log = TRUE)
assay(rse_gene_SRP119675_filtered, "logCPM") <- logCPM

#===============================================================================
#Multi dimentional analysis
#===============================================================================

#make the SingleCellExperiment object to make reduction of dimensionality analysis

sce_gene_SRP119675_filtered <- as(
  rse_gene_SRP119675_filtered,
  "SingleCellExperiment"
)

#make the PCA analysis with log-normalized counts
sce_gene_SRP119675_filtered <- scater::runPCA(
  sce_gene_SRP119675_filtered,
  exprs_values = "logCPM",
  ncomponents = 10
)

#make the MDS
mat <- edgeR::cpm(dge, log = TRUE)
dist_mat <- dist(t(mat))
mds <- cmdscale(dist_mat, k = 2)
reducedDims(sce_gene_SRP119675_filtered)$MDS <- mds
red.dim <- reducedDim(sce_gene_SRP119675_filtered, "MDS")


#===============================================================================
#Reduction dimentionality Visualization
#===============================================================================

# now we initialized some variables for ploting the previous analysis
se <- sce_gene_SRP119675_filtered
colormap <- ExperimentColorMap()
se <- iSEE::cleanDataset(se)
colormap <- synchronizeAssays(colormap, se)

#adjust the data of the plot with the PCA data
red.dim <- reducedDim(se, "PCA")
plot.data <- data.frame(
  X = red.dim[, 1],
  Y = red.dim[, 2],
  row.names = colnames(se)
)

plot.data$ColorBy <- colData(se)[, "sra_attribute.genotype"]
set.seed(142)
plot.data <- plot.data[sample(nrow(plot.data)), , drop = FALSE]


PCA_plot <- ggplot() +
  geom_point(
    aes(x = X, y = Y, color = ColorBy),
    alpha = 1,
    plot.data,
    size = 1
  ) +
  labs(
    x = "Dimension 1",
    y = "Dimension 2",
    color = "sra_attribute.genotype",
    title = "PCA"
  ) +
  coord_cartesian(
    xlim = range(plot.data$X, na.rm = TRUE),
    ylim = range(plot.data$Y, na.rm = TRUE),
    expand = TRUE
  ) +
  scale_color_manual(
    values = colDataColorMap(
      colormap,
      "sra_attribute.genotype",
      discrete = TRUE
    )(4),
    na.value = 'grey50',
    drop = FALSE
  ) +
  scale_fill_manual(
    values = colDataColorMap(
      colormap,
      "sra_attribute.genotype",
      discrete = TRUE
    )(4),
    na.value = 'grey50',
    drop = FALSE
  ) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    legend.box = 'vertical',
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    title = element_text(size = 12)
  )


# Adjust with the MDS data
red.dim <- reducedDim(sce_gene_SRP119675_filtered, "MDS")
plot.data <- data.frame(
  X = red.dim[, 1],
  Y = red.dim[, 2],
  row.names = colnames(se)
)

plot.data$ColorBy <- colData(se)[, "sra_attribute.genotype"]
set.seed(142)
plot.data <- plot.data[sample(nrow(plot.data)), , drop = FALSE]


MDS_plot <- ggplot() +
  geom_point(
    aes(x = X, y = Y, color = ColorBy),
    alpha = 1,
    plot.data,
    size = 1
  ) +
  labs(
    x = "Dimension 1",
    y = "Dimension 2",
    color = "sra_attribute.genotype",
    title = "MDS"
  ) +
  coord_cartesian(
    xlim = range(plot.data$X, na.rm = TRUE),
    ylim = range(plot.data$Y, na.rm = TRUE),
    expand = TRUE
  ) +
  scale_color_manual(
    values = colDataColorMap(
      colormap,
      "sra_attribute.genotype",
      discrete = TRUE
    )(4),
    na.value = 'grey50',
    drop = FALSE
  ) +
  scale_fill_manual(
    values = colDataColorMap(
      colormap,
      "sra_attribute.genotype",
      discrete = TRUE
    )(4),
    na.value = 'grey50',
    drop = FALSE
  ) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    legend.box = 'vertical',
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    title = element_text(size = 12)
  )

#===============================================================================
#Variance Visualization
#===============================================================================

# Now we run the anaylisis for the variancePrtion visualization

form_fixed <- form <- ~ sra_attribute.genotype +
  sra_attribute.age +
  sra_attribute.gender

varPart <- variancePartition::fitExtractVarPartModel(
  logCPM,
  form_fixed,
  colData(rse_gene_SRP119675_filtered)
)
vp <- sortCols(varPart)


form_random <- ~ (1 | sra_attribute.genotype) +
  (1 | sra_attribute.age) +
  (1 | sra_attribute.gender) +
  (1 | sra_attribute.source_name)

C <- canCorPairs(form_random, colData(rse_gene_SRP119675_filtered))

#save the plots
variance_var_plot <- plotVarPart(vp)


#saving correlaition plot
out_file <- file.path(
  outdir,
  "plots",
  "proccesed_data",
  "Correlation_plot_preprocessed.png"
)

png(out_file, width = 2000, height = 2000, res = 300)
plotCorrMatrix(C)
dev.off()


#===============================================================================
#Box plot Visualization
#===============================================================================

#Finaly we can do the box plot of the genes proportions vs our variables of interest

box_plot_genProp <- ggplot(
  as.data.frame(colData(rse_gene_SRP119675_filtered)),
  aes(
    y = assigned_gene_prop,
    x = sra_attribute.genotype,
    fill = sra_attribute.genotype
  )
) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Age Group") +
  scale_fill_discrete(name = "Genotype")


#===============================================================================
#Saving objects
#===============================================================================
save(
  rse_gene_SRP119675_filtered,
  file = file.path(
    outdir,
    "objects_R",
    "rse_gene_SRP119675_preprocessed.Rdata"
  )
)

save(
  dge,
  file = file.path(
    outdir,
    "objects_R",
    "rse_gene_SRP119675_dge_normalized.Rdata"
  )
)

save_plot(scater_plot, "Ponit_distribution_geneProp.png")
save_plot(histogram_plot, "Histogram_distribution_geneProp.png")
save_plot(PCA_plot, "preprocesed_PCA.png")
save_plot(MDS_plot, "preprocessed_MDS.png")
save_plot(variance_var_plot, "Variance_plot_preprocessed.png")
save_plot(box_plot_genProp, "Boxplot_geneProp_genotype.png")
