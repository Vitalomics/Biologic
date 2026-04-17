#algoritmo de analisis Vitalomics (Betina, febrero 2026)

##############################################################################
##############################################################################

# A) EDAD EPIGENETICA

library(EpigeneticAgePipeline)
library(meffil)
library(methylclock)
library(methylCIPHER)
library(mgcv)

###############################################################################
##############################################################################

#procesamiento de IDATs y calculo de relojes con DunedinPace

EpigeneticAgePipeline::main(
  inputDirectory = file.path(getwd(), "IDAT files and sample sheet"),
  outputDirectory = file.path(getwd(), "output"),
  normalize = TRUE,
  useBeta = FALSE,
  arrayType = "EPIC",
  useSampleSheet = TRUE,
  doParallel = TRUE,
  writeBeta = TRUE
)

#############################################
#Calculo de clocks y Age acceleration (residuals) con methylclock

sample_sheet <- read.csv2(file.path(getwd(),
                          "IDAT files and sample sheet/Sample_Sheet.csv"))

betas <- read.csv(file.path(getwd(),"output/extractedBetaValues.csv"))
rownames(betas) <- betas$X

sex <- sample_sheet$Sex
sex[sex %in% c("m", "male", "1")] <- "male"
sex[sex %in% c("f", "female", "2")] <- "female"

clocks <- methylclock::DNAmAge(x=betas,
                              method = c("horvath", "hannum", "levine"),
                              age = sample_sheet$Age,
                              sex = sex,
                              cell.count=TRUE)


#calculo de otros relojes con methylCIPHER
otherClocks <- c("calcDNAmTL","calcEpiTOC2","calcHypoClock","calcMiAge",
                "calcZhang2019","calcHRSInChPhenoAge",
                "calcSmokingMcCartney")

clocks2 <- methylCIPHER::calcUserClocks(otherClocks, 
                                        t(betas[,-1]), 
                                        sample_sheet, 
                                        imputation = F)

#cargo medias y sds de cohorte de referencia para calcular z
ref_means_clocks2 <- readRDS(file.path(getwd(),"inputs/ref means clocks2.rds"))

ref_sds_clocks2 <- readRDS(file.path(getwd(),"inputs/ref sds clocks2.rds"))

# variables a estandarizar
vars <- c("DNAmTL","epiTOC2","hypoClock","MiAge","Smoking_McCartney")

# Estándarizar con medias y sd de referencia
for (v in vars) {
  # nombre de la nueva columna z
  z_col <- paste0(v, "_z")
  
  # crear la columna z score
  clocks2[[z_col]] <- (clocks2[[v]] - ref_means_clocks2[v]) / ref_sds_clocks2[v]
}

#calculo conteo de celulas como se entreno en cohorte de referencia
cell_counts <- meffil::meffil.estimate.cell.counts.from.betas(
  as.matrix(betas[,-1]),
  "blood gse35069 complete"
) 

#creo df con clocks2 estandarizados , cell counts, age y female para meter al modelo
#que corrige por edad, sexo y composicion celular

clocks2_z <- data.frame(
  clocks2[,paste0(vars,"_z")],
  cell_counts,
  age=clocks2$Age,
  sex=sample_sheet$sex
)

#cargar modelos
fit_TL <- readRDS(file.path(getwd(),"modelos/TL_model.rds"))
fit_TOC2 <- readRDS(file.path(getwd(),"modelos/TOC2_model.rds"))
fit_hypoClock <- readRDS(file.path(getwd(),"modelos/hypoClock_model.rds"))
fit_MiAge <- readRDS(file.path(getwd(),"modelos/MiAge_model.rds"))
fit_Smoking <- readRDS(file.path(getwd(),"modelos/Smoking_model.rds"))

clocks2$expected_TL <- predict(fit_TL, newdata = clocks2_z)
clocks2$expected_TOC2 <- predict(fit_TOC2, newdata = clocks2_z)
clocks2$expected_hypoClock <- predict(fit_hypoClock, newdata = clocks2_z)
clocks2$expected_MiAge <- predict(fit_MiAge, newdata = clocks2_z)
clocks2$expected_Smoking <- predict(fit_Smoking, newdata = clocks2_z)
clocks2$age <- sample_sheet$Age

#funcion GAM termina de eliminar la dependencia del modelo con la edad cronologica
clocks2$Accel_TL <- resid(
  mgcv::gam(expected_TL ~ s(age), data = clocks2, method = "REML")
)

clocks2$Accel_TOC2 <- resid(
  mgcv::gam(expected_TOC2 ~ s(age), data = clocks2, method = "REML")
)

clocks2$Accel_hypoClock <- resid(
  mgcv::gam(expected_hypoClock ~ s(age), data = clocks2, method = "REML")
)

clocks2$Accel_MiAge <- resid(
  mgcv::gam(expected_MiAge ~ s(age), data = clocks2, method = "REML")
)

clocks2$Accel_Smoking <- resid(
  mgcv::gam(expected_Smoking ~ s(age), data = clocks2, method = "REML")
)

DunedinPACE <- read.csv(file.path(getwd(),"output/DunedinPACESampleData.csv"))

epigenetic_age_out <- data.frame(
  SampleID=clocks$id,
  Sex = sample_sheet$Female,
  Age=clocks$age,
  Horvath = clocks$Horvath,
  AgeAccel.Horvath = clocks$ageAcc2.Horvath,
  Hannum=clocks$Hannum,
  AgeAccel.Hannum=clocks$ageAcc2.Hannum,
  Levine = clocks$Levine,
  AgeAccel.Levine = clocks$ageAcc2.Levine,
  PhenoAge = clocks2$HRSInChPhenoAge,
  Zhang = clocks2$Zhang2019,
  DunedinPace = DunedinPACE$DunedinPACE,
  TLength=clocks2$expected_TL,
  TLength_Accel=clocks2$Accel_TL,
  TOC2=clocks2$expected_TOC2,
  TOC2_Accel=clocks2$Accel_TOC2,
  hypoClock=clocks2$expected_hypoClock,
  hypoClock_Accel=clocks2$Accel_hypoClock,
  MiAge=clocks2$expected_MiAge,
  MiAge_Accel=clocks2$Accel_MiAge,
  SmokingScore=clocks2$expected_Smoking,
  SmokingScore_Accel=clocks2$Accel_Smoking
)

write.csv2(epigenetic_age_out,file.path(getwd(),"output final/Epigenetic_Age_out.csv"))

###############################################################################
###############################################################################

# B) FITNESS

library(glmnet)
library(glmnetUtils)
library(broom)
library(dplyr)

###############################################################################
###############################################################################

#fitness metrics FUNCTIONS: adaptadas para funcionar con DunedinPACE
### SOURCE FUNCTIONS FOR ESTIMATING DNAM FITNESS MODELS (10/24/2022 Kristen McGreevy)

DNAmFitnessModels <- readRDS(file.path(getwd(),"modelos/DNAmFitnessModelsandFitAge_Oct2022.rds"))

`%!in%` = Negate(`%in%`)

DNAmFitness_Xvars <- c("DNAmGait_noAge", 
                       "DNAmGrip_noAge", 
                       "DNAmVO2max", 
                       "DNAmGait_wAge", 
                       "DNAmGrip_wAge", 
                       "DNAmFEV1_wAge", 
                       "DunedinPACE")

## DATAFRAME NEEDS TO HAVE CPGS IN COLUMNS WITH A COLUMN FOR SAMPLEID, AGE, AND FEMALE ## 
data_prep2 <- function(dataset, idvariable){
  
  # keep only the CpG loci you need to calculate the DNAm fitness biomarkers
  extract_these <- colnames(dataset)[which(colnames(dataset) %in% DNAmFitnessModels$AllCpGs)]
  output_data <- dataset[, c(idvariable, "Female", "Age", extract_these)]
  
  # update output_data with medians if we are missing any CpGs # 
  if(length(extract_these) != length(DNAmFitnessModels$AllCpGs)){
    # separate by sex to impute medians by sex. 
    data_fem <- output_data[output_data$Female == 1, ]
    data_male <- output_data[output_data$Female == 0, ]
    
    cpgs_toadd <- colnames(DNAmFitnessModels$Female_Medians_All) %!in% extract_these
    
    # set missing CpG values to medians from our training data. 
    data_fem <- data.frame(data_fem, DNAmFitnessModels$Female_Medians_All[, cpgs_toadd])
    data_male <- data.frame(data_male, DNAmFitnessModels$Male_Medians_All[, cpgs_toadd])
    
    output_data <- rbind(data_fem, data_male)
    
    print(paste0("Total ", sum(cpgs_toadd), 
                 " Missing CpGs that are assigned median values from training data"))
  }
  return(output_data)
}

# Function to calculate all DNAm fitness estimates #
DNAmFitnessEstimators <- function(data, IDvar){
  
  data_fem <- data[data$Female ==1,]
  data_male <- data[data$Female ==0,]
  
  fem_est1 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Gait_noAge_Females, IDvar) # gait without age
  fem_est2 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Grip_noAge_Females, IDvar) # grip
  fem_est3 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$VO2maxModel, IDvar) # vo2max
  fem_est4 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Gait_wAge_Females, IDvar) # gait w age
  fem_est5 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Grip_wAge_Females, IDvar) # grip w age
  fem_est6 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$FEV1_wAge_Females, IDvar) # fev1 w age
  
  male_est1 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Gait_noAge_Males, IDvar) # gait
  male_est2 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Grip_noAge_Males, IDvar) # grip
  male_est3 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$VO2maxModel, IDvar) # vo2max
  male_est4 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Gait_wAge_Males, IDvar) # gait
  male_est5 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Grip_wAge_Males, IDvar) # grip
  male_est6 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$FEV1_wAge_Males, IDvar) # fev1
  
  # just remerging male and female estimates into one dataset. 
  fem_est1 <- rbind(fem_est1, male_est1)
  fem_est2 <- rbind(fem_est2, male_est2)
  fem_est3 <- rbind(fem_est3, male_est3)
  fem_est4 <- rbind(fem_est4, male_est4)
  fem_est5 <- rbind(fem_est5, male_est5)
  fem_est6 <- rbind(fem_est6, male_est6)
  
  # reassigning values to proper names 
  fem_est1[, DNAmFitness_Xvars[1]] <- fem_est1$DNAmEstimate
  fem_est2[, DNAmFitness_Xvars[2]] <- fem_est2$DNAmEstimate
  fem_est3[, DNAmFitness_Xvars[3]] <- fem_est3$DNAmEstimate
  fem_est4[, DNAmFitness_Xvars[4]] <- fem_est4$DNAmEstimate
  fem_est5[, DNAmFitness_Xvars[5]] <- fem_est5$DNAmEstimate
  fem_est6[, DNAmFitness_Xvars[6]] <- fem_est6$DNAmEstimate
  
  all_ests <- Reduce(function(x,y) merge(x = x, y = y, by = "ID", all.x = TRUE, all.y = TRUE), 
                     list(fem_est1[!is.na(fem_est1$ID), 2:3],
                          fem_est2[!is.na(fem_est2$ID), 2:3],
                          fem_est3[!is.na(fem_est3$ID), 2:3],
                          fem_est4[!is.na(fem_est4$ID), 2:3],
                          fem_est5[!is.na(fem_est5$ID), 2:3],
                          fem_est6[!is.na(fem_est6$ID), 2:3]))
  
  # combine with original dataframe
  match1 <- match(data[, IDvar], all_ests$ID)
  data_and_est <- data.frame(data, all_ests[match1, ])
  
  return(data_and_est)
}

# Function to provide estimates for any 1 DNAm fitness models 
# Use this if you want to calculate individual DNAmFitness estimates and not all 6
# TidyModel is a specific model in DNAmFitnessModels list
DNAmEstimatorAnyModel <- function(dataset, TidyModel, IDvar){
  
  int_length <- nrow(dataset)
  mod_length <- length(TidyModel$term)
  
  Xdat <- dataset[, colnames(dataset) %in% TidyModel$term]
  Xdat <- data.frame("Intercept" = rep(1, int_length), Xdat)
  Xdatnew <- as.matrix(Xdat[, c("Intercept", TidyModel$term[2:mod_length])])
  if(sum(colnames(Xdatnew)[2:mod_length] == TidyModel$term[2:mod_length]) == mod_length-1){
    estimate <- Xdatnew %*% TidyModel$estimate
  }
  if(sum(colnames(Xdatnew)[2:mod_length] == TidyModel$term[2:mod_length]) != mod_length-1){
    print("Not All Columns in New Dataframe")
    estimate <- rep(NA, int_length)
  }
  
  # put Estimates into one dataset;
  est_data <- data.frame(DNAmEstimate = estimate, ID = dataset[, {{IDvar}}])
  
  return(est_data)
}

### Nueva funcion Betina para calcular VitalFitScore (Dunedin en lugar de Grimage)

VitalFitEstimator_referenceZ <- function(
    patient_data,
    reference_stats,
    IDvar = "SampleID",
    sex_var = "Female"
){
  
  req_vars <- c(IDvar,
                sex_var,
                "DNAmVO2max",
                "DNAmGrip_noAge",
                "DNAmGait_noAge",
                "DunedinPACE")
  
  stopifnot(all(req_vars %in% colnames(patient_data)))
  
  # unir stats de referencia por sexo
  data <- patient_data %>%
    left_join(reference_stats, by = sex_var)
  
  # chequeo de merge
  if (any(is.na(data$vo2_mean)))
    stop("Hay pacientes sin stats de referencia para su sexo")
  
  # z-scores definidos por la referencia
  data <- data %>%
    mutate(
      z_vo2  = (DNAmVO2max      - vo2_mean)     / vo2_sd,
      z_grip = (DNAmGrip_noAge  - grip_mean)    / grip_sd,
      z_gait = (DNAmGait_noAge  - gait_mean)    / gait_sd,
      z_pace = (DunedinPACE         - dunedin_mean) / dunedin_sd
    )
  
  # índice VitalFitScore
  data$DNAm_VitalFitScore <-
    0.30 * data$z_vo2 +
    0.30 * data$z_grip +
    0.20 * data$z_gait -
    0.20 * data$z_pace
  
  return(
    data[, c(IDvar,
             "DNAm_VitalFitScore",
             "z_vo2",
             "z_grip",
             "z_gait",
             "z_pace")]
  )
}

##########################################
##########################################

# DATA preparation

beta_fit <- as.data.frame(t(betas[,-1]))
beta_fit$SampleID <- rownames(beta_fit)
dunedin <- read.csv(file.path(getwd(), "output/DunedinPACESampleData.csv"))

fit_input <- cbind(beta_fit,
                   DunedinPace=dunedin$DunedinPACE,
                   Age=sample_sheet$Age,
                   Female=sample_sheet$Female)

fit_prep <- data_prep2(
  dataset = fit_input,
  idvariable = "SampleID"
)

fit_fitness <- DNAmFitnessEstimators(
  data = fit_prep,
  IDvar = "SampleID"
)

fit_fitness$DunedinPACE <- dunedin$DunedinPACE[match(fit_fitness$SampleID,
                                                             dunedin$X)]

#calculo de VitalFit score

# stats de referencia
VitalFit_ref <- readRDS(file.path(getwd(),"inputs/VitalFit_reference_stats.rds"))

########################################################
########################################################
# cálculo de VitalFit score en pacientes

VitalFit_score <- VitalFitEstimator_referenceZ(
  patient_data   = fit_fitness,
  reference_stats = VitalFit_ref,
  IDvar = "SampleID",
  sex_var = "Female"
)
######################################################
######################################################

#creo df con VitalFit score, cell counts, age y female para meter al modelo
#que corrige el score por edad y composicion celular

VitalFit_df <- data.frame(
  VitalFit_score,
  cell_counts,
  age=sample_sheet$Age,
  Female=sample_sheet$Female
)

#cargar modelo
fit_VitalFitScore <- readRDS(file.path(getwd(),"modelos/VitalFitScore_model.rds"))
VitalFit_df$expected_VitalFit <- predict(fit_VitalFitScore, newdata = VitalFit_df)

#funcion GAM termina de eliminar la dependencia del modelo con la edad cronologica
VitalFit_df$VitalFitAccel <- resid(
  mgcv::gam(expected_VitalFit ~ s(age), data = VitalFit_df, method = "REML")
)

fit_out <- data.frame(SampleID=VitalFit_df$SampleID,
                      Age=sample_sheet$Age,
                      Female=sample_sheet$Female,
                      VitalFitScore=VitalFit_df$expected_VitalFit,
                      VitalFitAccel= VitalFit_df$VitalFitAccel,
                      Gait=fit_fitness$DNAmGait_noAge,
                      Grip=fit_fitness$DNAmGrip_noAge,
                      VO2max=fit_fitness$DNAmVO2max)

write.csv2(fit_out, file.path(getwd(),"output final/VitalFit_out.csv"))

###############################################################################

# C) EPITIPOS

library(readxl)

###############################################################################

# Funciones para calcular score por firma

# Funcion que arma la matriz de z scores z_matrix, estandarizando los betaValues
# del paciente por la media y sd de poblacion de referencia por sexo


build_cpg_zmatrix_by_sex <- function(
    betas,           # matrix CpG x SampleID (beta values)
    sample_sheet,    # data.frame con SampleID y Female (0/1)
    ref_stats        # data.frame con CpG, Female, mean_ref, sd_ref
) {
  
  ## -------------------------------
  ## Checks básicos
  ## -------------------------------
  stopifnot(
    is.matrix(betas),
    all(c("SampleID", "Female") %in% colnames(sample_sheet)),
    all(c("CpG", "Female", "mean_ref", "sd_ref") %in% colnames(ref_stats))
  )
  
  ## -------------------------------
  ## CpGs comunes entre cohortes
  ## -------------------------------
  common_cpgs <- intersect(
    rownames(betas),
    ref_stats$CpG
  )
  
  if (length(common_cpgs) == 0) {
    stop("No hay CpGs en común entre betas y ref_stats")
  }
  
  ## -------------------------------
  ## Inicializar matriz de z-scores
  ## -------------------------------
  z_matrix <- matrix(
    NA_real_,
    nrow = length(common_cpgs),
    ncol = nrow(sample_sheet),
    dimnames = list(common_cpgs, sample_sheet$SampleID)
  )
  
  ## -------------------------------
  ## Calcular z-scores por muestra
  ## -------------------------------
  for (i in seq_len(nrow(sample_sheet))) {
    
    pid   <- sample_sheet$SampleID[i]
    sex_i <- sample_sheet$Female[i]
    
    # stats de referencia para ese sexo
    ref_i <- ref_stats[
      ref_stats$Female == sex_i &
        ref_stats$CpG %in% common_cpgs,
    ]
    
    # asegurar mismo orden CpG
    ref_i <- ref_i[match(common_cpgs, ref_i$CpG), ]
    
    stopifnot(nrow(ref_i) == length(common_cpgs))
    
    z_matrix[, pid] <-
      (betas[common_cpgs, pid] - ref_i$mean_ref) /
      ref_i$sd_ref
  }
  
  ## -------------------------------
  ## Output
  ## -------------------------------
  return(z_matrix)
}


## Funcion para obtener matriz de z scores que pasan treshold por signature, y
# los multiplica por el signo de la firma, de modo que valores positivos de 
# z aligned indican concordancia, y negativos discordancia

score_sig_matrix <- function(signature, z_matrix, z_thr) {
  
  cpgs <- intersect(names(signature), rownames(z_matrix))
  #if (length(cpgs) < 5) return(NULL)
  
  z_sub <- z_matrix[cpgs, , drop = FALSE]
  dir   <- signature[cpgs]
  
  # alineación direccional
  z_aligned <- sweep(z_sub, 1, dir, "*")
  
  # aplicar umbral
  active <- abs(z_aligned) >= z_thr
  z_active <- z_aligned
  z_active[!active] <- NA
  
  list(
    z_aligned = z_active,
    n_active  = colSums(active, na.rm = TRUE),
    n_total   = length(cpgs)
  )
}

## Funcion que calcula el z aligned medio por signature, y le asigna un peso por
# tamanio de la signature

sig_score <- function(score_obj) {
  
  z_mat <- score_obj$z_aligned
  
  raw_score <- colMeans(z_mat, na.rm = TRUE)
  weight    <- score_obj$n_active / score_obj$n_total
  
  raw_score * weight
}

## Funcion que realiza un test estadistico de permutacion, para evaluar si el 
# z_aligned obtenido para la signature es distinto de uno obtenido con una CGs 
# seleccionadas al azar (random signature)

perm_test_sig <- function(
    signature,
    z_matrix,
    n_perm,
    z_thr
){
  
  ## score real
  real_mat   <- score_sig_matrix(signature, z_matrix, z_thr)
  real_score <- sig_score(real_mat)
  
  n_cpg    <- length(signature)
  all_cpgs <- rownames(z_matrix)
  
  ## permutaciones
  perm_scores <- replicate(n_perm, {
    
    rand_cpgs <- sample(all_cpgs, n_cpg)
    rand_dir  <- sample(c(-1, 1), n_cpg, replace = TRUE)
    rand_sig  <- setNames(rand_dir, rand_cpgs)
    
    rand_mat <- score_sig_matrix(rand_sig, z_matrix, z_thr)
    sig_score(rand_mat)
  })
  
  perm_scores <- t(perm_scores)  # pacientes × permutaciones
  
  # FORZAMOS dimensiones correctas
  perm_scores <- matrix(
    perm_scores,
    nrow = length(real_score),
    ncol = n_perm,
    dimnames = list(names(real_score), NULL)
  )
  
  ## p-values empíricos (one-sided, corrected)
  pvals <- sapply(seq_along(real_score), function(i){
    (1 + sum(perm_scores[i, ] >= real_score[i], na.rm = TRUE)) /
      (1 + n_perm)
  })
  
  names(pvals) <- names(real_score)
  
  ## z-score empírico
  z_empirical <- (real_score - rowMeans(perm_scores, na.rm = TRUE)) /
    apply(perm_scores, 1, sd, na.rm = TRUE)
  
  ## salida
  data.frame(
    real_score   = real_score,
    pval         = pvals,
    z_empirical  = z_empirical,
    mean_perm_scores  = rowMeans(perm_scores, na.rm = TRUE)
  )
}


## Funcion que calcula el sig_score y el resultado estadistico para todas las 
# signatures

run_signature_enrichment <- function(
    signatures,
    z_matrix,
    n_perm,
    z_thr,
    min_cpgs,
    verbose = TRUE
) {
  
  results <- list()
  
  for (sig_name in names(signatures)) {
    
    if (verbose)
      message("Procesando signature: ", sig_name)
    
    signature <- signatures[[sig_name]]
    
    # CpGs que realmente entran
    cpgs_use <- intersect(names(signature), rownames(z_matrix))
    
    if (length(cpgs_use) < min_cpgs) {
      if (verbose)
        message("  ↪ saltada (", length(cpgs_use), " CpGs)")
      next
    }
    
    stat <- perm_test_sig(
      signature = signature,
      z_matrix  = z_matrix,
      n_perm    = n_perm,
      z_thr     = z_thr
    )
    
    stat$SampleID  <- rownames(stat)
    stat$Signature <- sig_name
    stat$n_CpGs    <- length(cpgs_use)
    
    results[[sig_name]] <- stat
  }
  
  # combinar resultados
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  
  ## Corrección FDR por paciente
  out <- out %>%
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(
      FDR = p.adjust(pval, method = "BH")
    ) %>%
    dplyr::ungroup()
  
  # columnas finales
  out[, c(
    "SampleID",
    "Signature",
    "real_score",
    "z_empirical",
    "pval",
    "FDR",
    "mean_perm_scores",
    "n_CpGs"
  )]
}

# Funcion que Filtra firmas significativas (FDR + z empírico), anota firmas con Epitype 
# y Epi.DR, calcula un EpiScore ponderado por: dirección (D vs R), significancia (−log10 FDR)
# tamaño del efecto (z), resume por SampleID × Epitype, clasifica el efecto 
# (strong / moderate / neutral / mixed)

compute_EpiScores <- function(
    sig_enrichment,
    epitypes_df,
    FDR_cutoff = 0.05,
    z_cutoff = 1.96,
    lambda = 0.5
) {
  
  ## 1. Filtrar firmas significativas
  sig_results <- sig_enrichment %>%
    dplyr::filter(
      FDR < FDR_cutoff,
      abs(z_empirical) > z_cutoff
    )
  
  ## 2. Total de firmas por epitipo (denominador)
  total_sigs_epitype <- epitypes_df %>%
    dplyr::filter(!is.na(Epitype)) %>%
    dplyr::distinct(Epitype, Signature) %>%
    dplyr::count(Epitype, name = "total_signatures_epitype")
  
  ## 3. Anotación por epitipo
  sig_results_epitype <- sig_results %>%
    dplyr::mutate(
      Signature = as.character(Signature),
      SampleID  = as.character(SampleID)
    ) %>%
    dplyr::left_join(
      epitypes_df %>%
        dplyr::mutate(Signature = as.character(Signature)),
      by = "Signature"
    ) %>%
    dplyr::filter(!is.na(Epitype)) %>%
    dplyr::mutate(
      Epi.DR       = factor(Epi.DR, levels = c("D", "R")),
      FDR          = as.numeric(FDR),
      z_empirical  = as.numeric(z_empirical)
    )
  
  ## 4. Cálculo del EpiScore
  epi_scores <- sig_results_epitype %>%
    dplyr::mutate(
      w_DR     = ifelse(Epi.DR == "D", 1, lambda),
      w_FDR    = -log10(FDR),
      contrib  = w_DR * w_FDR * z_empirical,
      weight   = w_DR * w_FDR
    ) %>%
    dplyr::group_by(SampleID, Epitype) %>%
    dplyr::summarise(
      EpiScore = ifelse(sum(weight) > 0,
                        sum(contrib) / sum(weight),
                        NA_real_),
      n_signatures = sum(weight > 0),
      .groups = "drop"
    ) %>%
    dplyr::left_join(total_sigs_epitype, by = "Epitype") %>%
    dplyr::mutate(
      signature_ratio = n_signatures / total_signatures_epitype
    )
  
  ## 5. Clasificación del efecto
  epi_scores <- epi_scores %>%
    dplyr::mutate(
      effect = dplyr::case_when(
        abs(EpiScore) > 3   & signature_ratio > 0.4 ~ "strong",
        abs(EpiScore) > 1.5 & signature_ratio > 0.2 ~ "moderate",
        abs(EpiScore) < 1   & signature_ratio < 0.2 ~ "neutral",
        TRUE ~ "mixed"
      )
    )
  
  return(epi_scores)
}

################################################
################################################
#Analisis de firmas en muestras de pacientes

# cargo las firmas o signatures 
# Firmas EWAS: CpG y dirección 1/-1
signature_dir <- file.path(getwd(),"inputs/CG signatures")
signature_files <- list.files(signature_dir, pattern = "\\.xlsx$", full.names = TRUE)

# Leer todas las firmas en una lista
signatures <- lapply(signature_files, function(file) {
  df <- read_excel(file)
  setNames(as.numeric(df$trend), df$`Probe ID`)
})

names(signatures) <- tools::file_path_sans_ext(basename(signature_files))

#CGs stats (mean and sd) in reference cohorts
CG_ref_stats <- read.csv2(file.path(getwd(),"inputs/ref cohort CG stats.csv"))

#preparo matriz de betaValues
betas2 <- as.matrix(betas[,-1])
colnames(betas2) <- sample_sheet$SampleID #tienen que coincidir los samplesID

#####################################################

#construccion de matriz de z scores por stats segun sexo en ref
z_matrix <- build_cpg_zmatrix_by_sex(
              betas2,           # matrix CpG x SampleID (beta values)
              sample_sheet,    # data.frame con SampleID y Female (0/1)
              CG_ref_stats        # data.frame con CpG, Female, mean_ref, sd_ref
)


#######################################################

#firmas enriquecidas por paciente

sig_enrichment <- run_signature_enrichment(
  signatures = signatures,
  z_matrix   = z_matrix,
  n_perm     = 1000,
  z_thr      = 1.5,
  min_cpgs   = 1
)

#######################################################

# Calculo de Epitype Scores

epitypes_df <- read.csv2(file.path(getwd(),"Inputs/epitypes.csv"))

epi_scores <- compute_EpiScores(
  sig_enrichment = sig_enrichment,
  epitypes_df    = epitypes_df,
  FDR_cutoff     = 0.05,
  z_cutoff       = 1.96,
  lambda         = 0.5
)

write.csv2(epi_scores,file.path("output final/EpiScore.csv"))
