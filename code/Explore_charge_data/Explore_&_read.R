## Load required packages first
library("recount3") ## for obtain the data of RNAseq

outdir = "Results"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

#Now weare going to initialize the object of summarized experiment to narrow down the Data loading a preporccesing
rse_gene_SRP119675_raw <- recount3::create_rse_manual(
  project = "SRP119675",
  project_home = "data_sources/sra",
  organism = "mouse",
  annotation = "gencode_v23",
  type = "gene"
)

#Now we compute the read counts for the project and we store it in the assay of the summarized experiment object
assay(rse_gene_SRP119675_raw, "counts") <- compute_read_counts(
  rse_gene_SRP119675_raw
)

#we filter the sra attributes to obtain only those samples that contain al the 96features
valid <- stringr::str_count(
  rse_gene_SRP119675_raw$sra.sample_attributes,
  "\\|"
) ==
  6

rse_gene_SRP119675 <- rse_gene_SRP119675_raw[, valid]

#with the data charged we can expnand the atrributes of the summarized experiment object to obtain more information about the project and the samples
rse_gene_SRP119675 <- expand_sra_attributes(rse_gene_SRP119675)

#now with the data filtered and sra attributes expanded we can changue the data types of those atributes for the stadistical anaysis
rse_gene_SRP119675$sra_attribute.gender <- as.factor(
  rse_gene_SRP119675$sra_attribute.gender
)
rse_gene_SRP119675$sra_attribute.genotype <- as.factor(
  rse_gene_SRP119675$sra_attribute.genotype
)
rse_gene_SRP119675$sra_attribute.tissue <- as.factor(
  rse_gene_SRP119675$sra_attribute.tissue
)
rse_gene_SRP119675$sra_attribute.source_name <- as.factor(
  rse_gene_SRP119675$sra_attribute.source_name
)
rse_gene_SRP119675$sra_attribute.age <- as.factor(
  rse_gene_SRP119675$sra_attribute.age
)
rse_gene_SRP119675$sra_attribute.brain_region <- as.factor(
  rse_gene_SRP119675$sra_attribute.brain_region
)

dir.create(
  file.path(outdir, "objects_R"),
  recursive = TRUE,
  showWarnings = FALSE
)

#Finaly we store the object that we generated for future scripts
save(
  rse_gene_SRP119675,
  file = file.path(
    outdir,
    "objects_R",
    "rse_gene_SRP119675_parsed.Rdata"
  )
)
