require(maSigPro)
require(dplyr)

#-------------------------------------------------------------------#
data("data.abiotic")   #expression matrix
data("edesign.abiotic")   #experimental design matrix

matrix_used <- data.abiotic
raw_design_used <- edesign.abiotic

class(matrix_used); dim(matrix_used)
class(raw_design_used); dim(raw_design_used)

#-------------------------------------------------------------------#
#Defining the regression model
design_used <- make.design.matrix(raw_design_used, degree = 2)

class(design_used)
length(design_used)

design_used$dis
design_used$groups.vector
design_used$edesign

#-------------------------------------------------------------------#
#Finding significant genes
fit_used <- p.vector(matrix_used, design_used, Q=0.05, MT.adjust="BH")   #Q=, MT.adjust=

class(fit_used)
length(fit_used)

#view the details of the output
#the output: a list with 14 elements
fit_used$dis %>% class(); fit_used$dis %>% dim()
fit_used$SELEC %>% class(); fit_used$SELEC %>% dim()
fit_used$p.vector %>% class(); fit_used$p.vector %>% dim()
fit_used$p.adjusted %>% class(); fit_used$p.adjusted %>% dim()
fit_used$G %>% class(); fit_used$G %>% dim()
fit_used$g %>% class(); fit_used$g %>% dim()
fit_used$FDR %>% class(); fit_used$FDR %>% dim()
fit_used$i %>% class(); fit_used$i %>% dim()
fit_used$dat %>% class(); fit_used$dat %>% dim()
fit_used$min.obs %>% class(); fit_used$min.obs %>% dim()
fit_used$Q %>% class(); fit_used$Q %>% dim()
fit_used$groups.vector %>% class(); fit_used$groups.vector %>% dim()
fit_used$edesign %>% class(); fit_used$edesign %>% dim()
fit_used$family %>% class(); fit_used$family %>% dim()

#-------------------------------------------------------------------#
#Finding significant differences
tstep <- T.fit(fit_used, step.method="backward", alfa=0.05)   #step.method=, alfa=

class(tstep)
length(tstep)

#view the details of the output
#the output: a list with 14 elements
tstep$sol %>% class(); tstep$sol %>% dim()
tstep$sig.profiles %>% class(); tstep$sig.profiles %>% dim()
tstep$coefficients %>% class(); tstep$coefficients %>% dim()
tstep$group.coeffs %>% class(); tstep$group.coeffs %>% dim()
tstep$t.score %>% class(); tstep$t.score %>% dim()
tstep$dis %>% class(); tstep$dis %>% dim()
tstep$variables %>% class(); tstep$variables %>% dim()
tstep$G %>% class(); tstep$G %>% dim()
tstep$g %>% class(); tstep$g %>% dim()
tstep$dat %>% class(); tstep$dat %>% dim()
tstep$step.method %>% class(); tstep$step.method %>% dim()
tstep$groups.vector %>% class(); tstep$groups.vector %>% dim()
tstep$edesign %>% class(); tstep$edesign %>% dim()
tstep$influ.info %>% class(); tstep$influ.info %>% dim()

#-------------------------------------------------------------------#
#Obtaining lists of significant genes
sigs_final <- get.siggenes(tstep, rsq=0.7, vars="groups")   #rsq=, vars=

class(sigs_final)
length(sigs_final)

#view the details of the output
sigs_final$summary %>% class(); sigs_final$summary %>% dim()   #a dataframe show the DEG for every variable (time or dummy variable)
sigs_final$sig.genes$Control %>% class(); sigs_final$sig.genes$Control %>% dim()

#-------------------------------------------------------------------#
#visualization
suma2Venn(sigs_final$summary[, c(2:4)])   #Venn plot

see.genes(sigs_final$sig.genes$ColdvsControl, show.fit=T, dis=design_used$dis, cluster.method="hclust", cluster.data=1, k=4)
