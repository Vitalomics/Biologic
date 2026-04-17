#CastellaniLab Epigenetic Age Pipeline
############################################################################
############################################################################
#Install the following packages which act as dependencies for the pipeline
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Full Install
packages_to_install <- c(
  "ggplot2",
  "glmmTMB",
  "ggpubr",
  "devtools",
  "magick",
  "reshape2",
  "minfi",
  "FlowSorted.CordBlood.450k",
  "FlowSorted.CordBloodCombined.450k",
  "FlowSorted.Blood.EPIC",
  "BeadSorted.Saliva.EPIC",
  "FlowSorted.DLPFC.450k",
  "IlluminaHumanMethylation27kanno.ilmn12.hg19",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
  "IlluminaHumanMethylationMSAanno.ilm10a1.hg38",
  "IlluminaHumanMethylationEPICv2manifest",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylation27kmanifest",
  "IlluminaHumanMethylationMSAmanifest",
  "planet",
  "sesame",
  "methylclock",
  "ExperimentHub",
  "EpiDISH"
)

for (package in packages_to_install) {
  if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package)
  }
}

#instalacion libreria de Castellani Lab para analisis a partir de IDATs
remotes::install_github('CastellaniLab/EpigeneticAgePipeline')

#instalacion EpiSmokEr
library(devtools)
install_github("sailalithabollepalli/EpiSmokEr")

###############################################################################
###############################################################################
# instalacion methylCIPHER

devtools::install_github("MorganLevineLab/prcPhenoAge")
devtools::install_github("danbelsky/DunedinPoAm38")
devtools::install_github("MorganLevineLab/methylCIPHER")