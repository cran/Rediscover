## ----LoadFunctions, echo=FALSE, message=FALSE, warning=FALSE, results='hide'----
library(knitr)
opts_chunk$set(error = FALSE)
library(Rediscover)
library(dplyr)
library(kableExtra)
library(maftools)
# library(TCGAbiolinks) no  available in current version of bioconductor
library(parallel)
library(qvalue)

## ----style, echo = FALSE, results = 'asis'------------------------------------
##BiocStyle::markdown()

## ---- eval=FALSE--------------------------------------------------------------
#  
#  install.packages("Rediscover")
#  

## ----getPM_matrix, eval=TRUE--------------------------------------------------

data("A_example")

PMA <- getPM(A_example)
PMA[1:4,1:4]


## ----getPM_, eval=TRUE--------------------------------------------------------

data("A_Matrix")
class(A_Matrix)
PMA <- getPM(A_Matrix)
PMA[1:4,1:4]


## ----getPM_COAD, eval=TRUE----------------------------------------------------

data("TCGA_COAD")
PM_COAD <- getPM(TCGA_COAD)


## ----getMutex_matrix, eval=TRUE-----------------------------------------------

data("A_example")

PMA <- getPM(A_example)

mymutex <- getMutex(A=A_example,PM=PMA)


## ----getMutex_Matrix_2, eval=TRUE---------------------------------------------

data("A_Matrix")

PMA_Matrix <- getPM(A_Matrix)

mymutex <- getMutex(A=A_Matrix,PM=PMA_Matrix)


## ----getMutex_COAD, eval=TRUE-------------------------------------------------

data("TCGA_COAD")
data("PM_COAD")

COAD_mutex <- getMutex(TCGA_COAD[1:100,], PM_COAD[1:100,])


## ----getMutex_COAD_exact, eval=TRUE-------------------------------------------
data("TCGA_COAD")
data("PM_COAD")
COAD_mutex_exact <- getMutex(TCGA_COAD[1:100,], PM_COAD[1:100,],mixed = TRUE,th = 0.001)

## ----getMutex_COAD_exact_2, eval=TRUE-----------------------------------------
data("TCGA_COAD")
data("PM_COAD")
COAD_mutex_exact <- getMutex(TCGA_COAD[1:100,], PM_COAD[1:100,],mixed = FALSE,method = "Binomial")

## ----getMutexAB_matrix, eval=TRUE---------------------------------------------

data("A_example")
data("B_example")

PMA <- getPM(A_example)
PMB <- getPM(B_example)

mismutex <- getMutexAB(A=A_example, PM=PMA, B=B_example, PMB = PMB)


## ----getMutexAB_Matrix_2, eval=TRUE-------------------------------------------

data("A_Matrix")
data("B_Matrix")

PMA <- getPM(A_Matrix)
PMB <- getPM(B_Matrix)

mismutex <- getMutexAB(A=A_Matrix, PM=PMA, B=B_Matrix, PMB = PMB)


## ----getMutexAB_COAD, eval=TRUE-----------------------------------------------

data("TCGA_COAD")
data("PM_COAD")
data("AMP_COAD")
data("PM_AMP_COAD")

common <- intersect(colnames(TCGA_COAD), colnames(AMP_COAD))

keep <- match(common,colnames(TCGA_COAD))
TCGA_COAD_100 <- TCGA_COAD[1:100,keep]
PM_TCGA_COAD_100 <- PM_COAD[1:100,keep]

keep <- match(common,colnames(AMP_COAD))
AMP_COAD_100 <- AMP_COAD[1:100,keep]
PM_AMP_COAD_100 <- PM_AMP_COAD[1:100,keep]

mismutex <- getMutexAB(A=TCGA_COAD_100, PMA=PM_TCGA_COAD_100, 
                       B=AMP_COAD_100, PMB = PM_AMP_COAD_100)



## ----getMutexAB_COAD_2, eval=TRUE---------------------------------------------

data("TCGA_COAD")
data("PM_COAD")
data("AMP_COAD")
data("PM_AMP_COAD")

common <- intersect(colnames(TCGA_COAD), colnames(AMP_COAD))

keep <- match(common,colnames(TCGA_COAD))
TCGA_COAD_100 <- TCGA_COAD[1:100,keep]
PM_TCGA_COAD_100 <- PM_COAD[1:100,keep]

keep <- match(common,colnames(AMP_COAD))
AMP_COAD_150 <- AMP_COAD[1:150,keep]
PM_AMP_COAD_150 <- PM_AMP_COAD[1:150,keep]

mismutex <- getMutexAB(A=TCGA_COAD_100, PMA=PM_TCGA_COAD_100, 
                       B=AMP_COAD_150, PMB = PM_AMP_COAD_150)



## ----getMutexGroup_example, eval=TRUE-----------------------------------------
data("A_example")

A2 <- A_example[,1:40]
A2[1,1:10] <- 1
A2[2,1:10] <- 0
A2[3,1:10] <- 0
A2[1,11:20] <- 0
A2[2,11:20] <- 1
A2[3,11:20] <- 0
A2[1,21:30] <- 0
A2[2,21:30] <- 0
A2[3,21:30] <- 1

PM2 <- getPM(A2)
A <- A2[1:3,]
PM <- PM2[1:3,]

## ----getMutexGroup_image, eval=TRUE, echo=FALSE-------------------------------
image(Matrix(A)) # These two genes are mutually exclusive (to a certain extent)

## ----getMutexGroup_studies, eval=TRUE-----------------------------------------
getMutexGroup(A, PM, "Impurity")
getMutexGroup(A, PM, "Coverage")
getMutexGroup(A, PM, "Exclusivity")


## ----maftools_COAD, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE, results='hide'----
#  coad.maf <- GDCquery_Maf("COAD", pipelines = "muse") %>% read.maf

## ----oncoplot_COAD, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE-------
#  oncoplot(coad.maf, top = 35)

## ----somInter_COAD, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE, results='hide'----
#  somaticInteractions(maf = coad.maf, top = 35, pvalue = c(1e-2, 2e-3))

## ----discsomInter_COAD, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE, results='hide'----
#  discoversomaticInteractions(maf = coad.maf, top = 35, pvalue = c(1e-2, 2e-3),getMutexMixed=FALSE)

## ----qvalue, eval=TRUE--------------------------------------------------------

data("TCGA_COAD")
data("PM_COAD")
COAD_mutex <- getMutex(TCGA_COAD[1:100,], PM_COAD[1:100,])
COAD_mutex_qvalue <- COAD_mutex
COAD_mutex_qvalue@x <- qvalue::qvalue(COAD_mutex_qvalue@x)$qvalue


## -----------------------------------------------------------------------------
sessionInfo()

