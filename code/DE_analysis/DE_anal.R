## Load the requiered R packages and the previous proccesed objects
library(limma)
library(SummarizedExperiment)
library(ggplot2)
library(ComplexHeatmap)
library(ggrepel)

#===============================================================================
#Pre loading and functions
#===============================================================================
#load the previous generated objects
load("Results/objects_R/rse_gene_SRP119675_preprocessed.Rdata")
load("Results/objects_R/rse_gene_SRP119675_dge_normalized.Rdata")


outdir = "Results"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

out_dir_plots <- ("Results/plots/DE_analysis")

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
  path = out_dir_plots
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
#Define the model and fixed DE analysis based on LM model
#===============================================================================

#First we define a interaction model for related the more variance expliative variables genotype and age
mod <- model.matrix(
  ~ sra_attribute.genotype * sra_attribute.age + sra_attribute.gender,
  data = colData(rse_gene_SRP119675_filtered)
)

#now we can fixed a completed linear model to our data based on the model matrix que define
vGene <- voom(dge, mod, plot = TRUE)

eb_results <- eBayes(lmFit(vGene))

#===============================================================================
#Volcano plots
#===============================================================================

#we store the names of the
coef_names <- colnames(mod)[-1]

#apply anonymus function lambda to make a top table for all the coeficents solves
tt_list <- lapply(coef_names, function(coef_name) {
  topTable(
    eb_results,
    coef = coef_name,
    number = nrow(rse_gene_SRP119675_filtered),
    sort.by = "none"
  )
})

#save the names of our comparisons
names(tt_list) <- coef_names

#We set the coutes for significant diferential expressed genes
lfc_cut <- 1
fdr_cut <- 0.05

for (name in names(tt_list)) {
  #extract the top table we are ploting
  res <- tt_list[[name]]
  plot_title <- paste(sub("sra_attribute.", "", name), "_vs_5XFAD", sep = "")
  #Clasify its value for them that are significant DE expressed (for coloring them)
  res$group <- "NS"
  res$group[res$adj.P.Val < fdr_cut & res$logFC > lfc_cut] <- "Up"
  res$group[res$adj.P.Val < fdr_cut & res$logFC < -lfc_cut] <- "Down"

  #also mark all those genes that arae significant for ploting its name
  label_genes <- res[
    res$adj.P.Val < 10 & abs(res$logFC) > lfc_cut,
  ]

  plot <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = group)) +
    geom_point(alpha = 0.7, size = 1.8) +
    # Set lines for the tresholds
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed") +
    scale_color_manual(
      values = c(
        "Down" = "#2C7BB6", # azul
        "NS" = "grey80", # gris tenue
        "Up" = "#D7191C" # rojo
      )
    ) +
    theme_minimal(base_size = 14) +
    labs(
      title = plot_title,
      x = "log2 Fold Change",
      y = "-log10 Adjusted P-value",
      color = ""
    )
  plot_fileN <- paste(make.names(plot_title), "_vulcano.png", sep = "")
  save_plot(plot, plot_fileN)
}


#===============================================================================
#Heat maps
#===============================================================================

#first we have to changue the names of the coefficents in order to be parsed by R
colnames(mod) <- make.names(colnames(mod))

#Now we can changue the contrast of the solves for having the comparison betwen both times with a constant genotype
contrast.matrix <- makeContrasts(
  C5XFAD_TREM2_vs_5XFAD_4m = sra_attribute.genotype5xFAD.BAC.TREM2.sra_attribute.age4m,
  C5XFAD_TREM2_vs_5XFAD_7m = sra_attribute.genotype5xFAD.BAC.TREM2.sra_attribute.age7m,
  levels = mod
)

eb_contrats_5X.TREM2 <- contrasts.fit(eb_results, contrast.matrix)
eb_contrats_5X.TREM2 <- eBayes(eb_contrats_5X.TREM2)


contrast.matrix <- makeContrasts(
  TREM2_vs_5XFAD_4m = sra_attribute.genotypeBAC.TREM2.sra_attribute.age4m,
  TREM2_vs_5XFAD_7m = sra_attribute.genotypeBAC.TREM2.sra_attribute.age7m,
  levels = mod
)

eb_contrats.TREM2 <- contrasts.fit(eb_results, contrast.matrix)
eb_contrats.TREM2 <- eBayes(eb_contrats.TREM2)


#store the objects in a list
contrast_list <- list(eb_contrats.TREM2, eb_contrats_5X.TREM2)
names(contrast_list) <- c("5XFAD_vs_TREM2 ", "5XFAD_vs_5XFAD-TREM2")

#now we can plot the heat maps

annotation <- contrast_list[[1]]$genes
annotation$gene_id_clean <- sub("\\..*$", "", annotation$gene_id)
gene_map <- setNames(annotation$gene_name, annotation$gene_id_clean)


for (nm in names(contrast_list)) {
  contrast <- contrast_list[[nm]]
  plot_title <- nm

  res_4m <- topTable(contrast, coef = 1, number = Inf)
  res_7m <- topTable(contrast, coef = 2, number = Inf)

  sig_genes <- unique(c(
    rownames(res_4m[res_4m$adj.P.Val < 0.05, ]),
    rownames(res_7m[res_7m$adj.P.Val < 0.05, ])
  ))

  if (length(sig_genes) == 0) {
    next
  }

  sig_genes <- head(sig_genes, 50)

  mat <- cbind(
    "4m" = contrast$coefficients[sig_genes, 1],
    "7m" = contrast$coefficients[sig_genes, 2]
  )

  clean_ids <- sub("\\..*$", "", rownames(mat))

  rownames(mat) <- ifelse(
    is.na(gene_map[clean_ids]),
    clean_ids,
    gene_map[clean_ids]
  )

  plot_ht <- Heatmap(
    mat,
    name = "log2 Fold Change",
    column_title = plot_title,
    cluster_columns = FALSE,
    col = circlize::colorRamp2(
      c(-2, 0, 2),
      c("blue", "white", "red")
    )
  )
  title_plot_path <- paste(out_dir_plots, "/", plot_title, "heat.png", sep = "")
  png(title_plot_path, width = 2000, height = 3000, res = 300)
  draw(plot_ht)
  dev.off()
}
