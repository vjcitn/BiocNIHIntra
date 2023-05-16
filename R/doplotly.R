# this code will process 4k PBMCs from TENx with SingleR and scater
# compute some approximate PCs and produce an interactive scatterplot
# in PC space with a biplot

#' get HPCA reference
#' @param \dots passed to celldex::HumanPrimaryCellAtlasData
#' @return SummarizedExperiment
#' @export
getHPCAreference = function(...)
    celldex::HumanPrimaryCellAtlasData(...)

#' get 4000+ PBMCs from TENx genomics
#' @param doNorm logical(1) if TRUE, apply scater::logNormCounts to PBMC SingleCellExperiment
#' @param useSyms logical(1) if TRUE, remove features lacking Symbol in rowData, and use Symbol 
#' as rownames
#' @return SingleCellExperiment
#' @export
get4kPBMC = function(doNorm=TRUE, useSyms=TRUE) {
    vps = TENxPBMCData::TENxPBMCData("pbmc4k")
    if (useSyms) {
      hassym = which(!is.na(rowData(vps)$Symbol))
      vps = vps[hassym,]
      rownames(vps) = rowData(vps)$Symbol
      }
    if (!doNorm) {
      return(vps)
      }
    scater::logNormCounts(vps)
}
#library(SingleR)
#vsing = SingleR(vps, hd, hd$label.main)
#vps$label.main = vsing$labels
#vpssds = rowSds(assay(vps))
#kp = which(vpssds > quantile(vpssds, .8))
#vpslim = vps[kp,]
#mat = t(as.matrix(assay(vpslim,2)))
#library(irlba)
#apca = prcomp_irlba(mat, 4)
#bb = "ba60f07ac0d972de80d549f8659a19707a263edc"
#fbcode = devtools::source_url("https://gist.githubusercontent.com/vjcitn/ebcb221cdb8b913d3f41d81be43c1294/raw/3f88df12c169ef0afdd138c578879217545f3045/filtered_biplot.R", sha1=bb)
#filtered_biplot = fbcode$value
#library(scales)
#library(ggplot2)
#npl = filtered_biplot(apca, vpslim, nvar=10)
#library(plotly)
#ggplotly(npl)
