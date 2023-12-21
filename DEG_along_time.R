require(maSigPro)
require(dplyr)

#------------------------------------------------------------------------------#
#create a simulate data
Time <- rep(c(1,5,10,24), each = 3)
Replicates <- rep(c(1:4), each = 3)
Group <- rep(1,12)

ss.edesign <- cbind(Time,Replicates,Group)
rownames(ss.edesign) <- paste("Array", c(1:12), sep = "")

ss.GENE <- function(n, r, var11 = 0.01, var12 = 0.02, var13 = 0.02,
                    var14 = 0.02, a1 = 0, a2 = 0, a3 = 0, a4 = 0) {
  tc.dat <- NULL
  for (i in 1:n) {
    gene <- c(rnorm(r, a1, var11), rnorm(r, a1, var12),
              rnorm(r, a3, var13), rnorm(r, a4, var14))
    tc.dat <- rbind(tc.dat, gene)
  }
  tc.dat }

flat <-ss.GENE(n = 85, r = 3)
induc <- ss.GENE(n = 5, r = 3, a1 = 0, a2 = 0.2, a3 = 0.6, a4 = 1)
sat <- ss.GENE(n = 5, r = 3, a1 = 0, a2 = 1, a3 = 1.1, a4 = 1.2)
ord <- ss.GENE(n = 5, r = 3, a1 = -0.8, a2 = -1, a3 = -1.3, a4 =-0.9)
ss.DATA <- rbind(flat, induc,sat,ord)
rownames(ss.DATA) <- paste("feature", c(1:100), sep = "")
colnames(ss.DATA) <- paste("Array", c(1:12), sep = "")

#expression matrix: ss.DATA
#experimental design: ss.design
ss.edesign
head(ss.DATA)

#------------------------------------------------------------------------------#
design_used <- make.design.matrix(ss.edesign, degree = 2)
fit_used <- p.vector(ss.DATA, design_used, Q = 0.05, MT.adjust = "BH")
tstep <- T.fit(fit_used, step.method = "backward", alfa = 0.05)
sigs_final <- get.siggenes(tstep, rsq = 0.7, vars = "groups")
see.genes(sigs_final$sig.genes$Group, show.fit=T, dis=design_used$dis, cluster.method="hclust", cluster.data=1, k=4)


