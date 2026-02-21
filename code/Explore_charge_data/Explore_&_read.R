## Load required packages first
library("recount3") ## for obtain the data of RNAseq

#Now weare going to initialize the object of summarized experiment to narrow down the Data loading a preporccesing
rse_gene_SRP187821_raw <- recount3::create_rse_manual(
  project = "SRP187821",
  project_home = "data_sources/sra",
  organism = "mouse",
  annotation = "gencode_v23",
  type = "gene"
)

#Now we compute the read counts for the project and we store it in the assay of the summarized experiment object
assay(rse_gene_SRP187821_raw, "counts") <- compute_read_counts(
  rse_gene_SRP187821_raw
)

#we filter the sra attributes to obtain only those samples that contain al the 9 features
valid <- stringr::str_count(
  rse_gene_SRP187821_raw$sra.sample_attributes,
  "\\|"
) ==
  9

rse_gene_SRP187821 <- rse_gene_SRP187821_raw[, valid]

#with the data charged we can expnand the atrributes of the summarized experiment object to obtain more information about the project and the samples
rse_gene_SRP187821 <- expand_sra_attributes(rse_gene_SRP187821)

#now with the data filtered and sra attributes expanded we can changue the data types of those atributes for the stadistical anaysis
rse_gene_SRP187821$sra_attribute.plate <- as.factor(
  rse_gene_SRP187821$sra_attribute.plate
)
rse_gene_SRP187821$sra_attribute.Sex <- as.factor(
  rse_gene_SRP187821$sra_attribute.Sex
)
rse_gene_SRP187821$sra_attribute.age <- as.numeric(
  rse_gene_SRP187821$sra_attribute.age
)
rse_gene_SRP187821$sra_attribute.genotype <- as.factor(
  rse_gene_SRP187821$sra_attribute.genotype
)
rse_gene_SRP187821$sra_attribute.lane <- as.factor(
  rse_gene_SRP187821$sra_attribute.lane
)

#also we create a new attribute for age for manage it as a factor on the 3 stages the experiment was made

rse_gene_SRP187821$sra_attribute.age_group <- factor(
  rse_gene_SRP187821$sra_attribute.age,
  levels = c(3, 6, 12, 21),
  labels = c("M3", "M6", "M12", "M21")
)

#Finaly we store the object that we generated for future scripts
save(
  rse_gene_SRP187821,
  file = "processed-data/objects_R/rse_gene_SRP187821_parsed.Rdata"
)
