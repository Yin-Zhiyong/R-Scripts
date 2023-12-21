require(Biobase)
require(splineTimeR)

#-------------------------------------------------------------------------------#
data(TCsimData)

#-------------------------------------------------------------------------------#
#example of expression matrix
matrix_used <- exprs(TCsimData)
class(matrix_used)
dim(matrix_used)
head(matrix_used)

#---------------------------#
#example of sample metadata
metadata_used <- pData(TCsimData)
rownames(metadata_used) <- metadata_used[,1]

class(metadata_used)
dim(metadata_used)
metadata_used

#-------------------------------------------------------------------------------#
phenoData_used <- new("AnnotatedDataFrame", data=metadata_used)   #construct metadata object
ExpressionSet_object <- ExpressionSet(assayData=matrix_used, phenoData=phenoData_used)   #construct ExpressionSet object
