Algoritmo Biológico Vitalomics

Autor: Betina González

Resumen general
El sistema procesa datos epigenéticos crudos (.IDAT) de muestras de sangre humana y genera métricas biológicas agregadas por paciente.
El pipeline es batch-oriented, determinístico, y consta de tres módulos independientes (Edad epigenética, Fitness y Epitipos), cada uno de los cuales genera archivos .csv finales.
No requiere interacción manual una vez que los archivos de entrada están correctamente ubicados.

El algoritmo consta de 3 módulos:

A) Edad epigenética

B) Fitness

C) Epitipos

Configuración de archivos

 

Los archivos del laboratorio se colocan en la carpeta “IDAT files and sample sheet”: 2 .IDAT por paciente y 1 Sample_Sheet.csv con descripción de las muestras. 
Inputs tiene archivos inputs para los distintos modulos
Modelos tiene los modelos entrenados en población de referencia (sanos)
Cada módulo exporta un .csv a una carpeta “output final”.


Módulo A 

Utiliza funciones de librerías (librería::funcion) especializadas para:
•	preprocesamiento de datos crudos (.IDAT)
•	cálculo de métricas de edad epigenética publicadas
•	ajuste estadístico por covariables biológicas

Diagrama de flujo

.IDAT + Sample_Sheet.csv
        ↓
BetaValues (CpGs × pacientes)
        ↓
Cálculo de relojes epigenéticos
        ↓
Estandarización y ajustes (edad, sexo, células)
        ↓
epigenetic_age_out.csv

Los modelos de ajuste ya están entrenados previamente y se cargan desde disco; no se reentrenan durante la ejecución.

Descripción de funciones:

EpigeneticAgePipeline::main: procesa los archivos crudos .IDAT y los convierte en una matriz de datos betaValues con SampleID del paciente en columnas y el id de sondas de detección en filas. Además exporta otros archivos con métricas y graficos, de los cuales solo un .csv se usa de input en modulo B (DunedinPACESampleData.csv)

Input: 
Archivos en carpeta “IDAT files and sample sheet”: 2 .IDAT por paciente y 1 Sample_Sheet.csv con SampleID, Sexo, Edad (obligatorias) y otras covariables de interés.

Output: 
ExtractedBetaValues.csv <- matriz de datos de pacientes x sondas 
XXXXSample.Data.csv <- archivos de métricas
MatrixPloTXXXX <- gráficos en .pdf, 
Archivos de texto <- métricas en 2 .txt y 1 .md.

methylclock::DNAmAge: calcula métricas de edad epigenetica.

Inputs: 
matriz de betaValues y sample_sheet.csv

output:
df llamado clocks con varias métricas de edad epigenética y aceleración del envejecimiento.

methylCIPHER::calcUserClocks: calcula otras métricas de edad epigenetica

inputs:
matriz de betaValues y sample_sheet.csv

output:
df llamado clocks2 con varias métricas de edad epigenética. 

De cloks2 se toman 5 variables que se estandarizan con media y sd de cohorte de referencia (archivos en carpeta “inputs”), 

Los valores z obtenidos se ajustan por edad, sexo y composición celular de la muestra de sangre, usando modelos lineales que se encuentran en la carpeta “modelos”.

meffil::meffil.estimate.cell.counts.from.betas: función que calcula composición celular de la muestra de sangre.

Inputs:
matriz de betaValues y sample_sheet.csv

output:
df llamado cell_counts.

mgcv::gam: función que aplica un modelo aditivo generalizado que termina de eliminar la dependencia del score con la edad cronológica del paciente. Sobre esta función se calculan los residuos (distancia del paciente al valor poblacional predicho por el modelo)

Inputs:
variables en clocks2

output:
nuevas variables en clocks2

Este módulo exporta el archivo epigenetic_age_out con la selección de métricas finales
Módulo B

Utiliza funciones y modelos definidos en https://github.com/kristenmcgreevy/DNAmFitAge.git, (McGreevy et al. 2023 doi: 10.18632/aging.204538) y una función propia que calcula un VitalFitScore. 
El VitalFitScore es un índice compuesto calculado como combinación lineal ponderada de métricas de fitness estandarizadas, los pesos son fijos.

Descripción de funciones:

dataprep2: prepara los datos para el cálculo de métricas de fitness

Inputs: 
matriz de betaValues, DunedinPACESampleData.csv, sample_sheet, modelo “DNAmFitnessModelsandFitAge_Oct2022”

Output: 
df “fit_prep”

DNAmFitnessEstimators: calcula parámetros de fitness

Inputs: 
fit_prep,

Output: 
df “fit_fitness”

VitalFitEstimator_referenceZ: función que estandariza variables de fitness según media y sd de población de referencia, y calcula el VitalFitScore aplicando una combinación lineal ponderada de las variables estandarizadas

Input: 
fit_fitness, en carpeta inputs: VitalFit_ref (rds con mean y sd de cohorte de referencia)

Output: 
df “VitalFit_score”

El VitalFit_score se ajusta por edad, sexo y composición celular con un modelo lineal entrenado en cohorte de referencia.

Este módulo exporta el archivo fit_out

Módulo C

Descripción conceptual 
Para cada paciente, se evalúa si su perfil epigenético es compatible con firmas previamente definidas de estados fisiopatológicos. Esto se hace comparando su señal contra firmas reales y contra firmas aleatorias, generando una medida estadística de enriquecimiento.
Utiliza funciones inéditas diseñadas para cotejar el perfil de betaValues del paciente con firmas epigenéticas (signatures) de estados fisiopatológicos validadas en muestras de sangre humanas. Las signatures están agrupadas en epitipos, que definen estados generales como “riesgo cardiovascular”, “inflamación”, etc. 

Descripción de funciones:

build_cpg_zmatrix_by_sex: construye una matriz con valores z a partir de betaValues, estandarizados por media y sd por sexo en cohorte de referencia

Inputs:
matriz de betaValues, sample_sheet, en carpeta Inputs: CG_ref_stats.rmd

output:
z_matrix

run_signature_enrichment: función que utiliza tres funciones, para calcular el sig_score y el resultado estadístico para todas las signatures.

score_sig_matrix: calcula matriz de z scores que pasan treshold (z > 1.5) por signature, y los multiplica por el signo de la firma, de modo que valores positivos de z_aligned indican concordancia, y negativos discordancia.

sig_score: calcula el z aligned medio por signature, y le asigna un peso por tamaño de la signature (cuantas sondas CG la definen). Valores altos y positivos de este score indican alta representación de esa firma en el paciente.

perm_test_sig: función que realiza un test estadístico de permutación, para evaluar si el z_aligned obtenido para la signature es distinto de uno obtenido con el mismo número de CGs seleccionadas al azar (random signature). Se realizan 1000 permutaciones por firma, y se estima la probabilidad de que el z_aligned obtenido sea distinto al obtenido en la población al azar. Por ultimo corrige el pval con Benjamini Hotchberg pata obtener el False Discovery rate o FDR.

Inputs:
z_matrix, en carpeta inputs carpeta “CG signatures” (96 archivos .xlsx)

Output:
df llamado “sig_enrichment” 

compute_EpiScores: función que filtra firmas significativas (FDR + z empírico), anota firmas con Epitype y Epi.DR, calcula un EpiScore ponderado por: dirección (D vs R), significancia (−log10 FDR), tamaño del efecto (z), resume por SampleID × Epitype, clasifica el efecto (strong / moderate / neutral / mixed)

Inputs:
sig_enrichment, en carpeta inputs: epitypes.csv

output:
df llamado epi_scores

Este módulo exporta el archivo EpiScore.csv
