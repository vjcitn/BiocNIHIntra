#' produce a biplot based on sharing of information between a `prcomp`-like output and a SingleCellExperiment
#' @param prcomp_output instance of prcomp, as produced by, e.g., irlba::prcomp_irlba
#' @param sce a SingleCellExperiment instance
#' @param sampvar character(1) name of a variable in `colData(sce)` that will be used to label samples
#' @param colorvar character(1) names of a variable in `colData(sce)` that will be used to color points
#' @param which numeric(2), gives dimensions of `prcomp_output$x` to use for point locations
#' @param nvar numeric(1) number of highly ranked features (by sum of squares of loadings) to use in biplot
#' @param shr numeric(1) a fudge factor to 'shrink' projection span relative to boundaries of display
#' @param \dots not used
#' @return a ggplot instance
#' @export
filtered_biplot = function (prcomp_output, sce, sampvar = "Barcode", colorvar = "label.main", 
    which = c(1, 2), nvar = 5, shr = 0.6, ...) 
{
    stopifnot(sampvar %in% names(colData(sce)))
    rownames(prcomp_output$x) = sce[[sampvar]]
    rownames(prcomp_output$rotation) = rownames(sce)
    proj = prcomp_output$x[, which]
    rot = prcomp_output$rot[, which]
    sss = function(x) sum(x^2)
    lens = apply(rot, 1, sss)
    kprot = rot[order(lens, decreasing = TRUE)[1:nvar], ]
    sca = max(abs(proj))
    fac = sca/max(abs(kprot[, 1]))
    pcs = paste0("PC", which)
    df1 = data.frame(proj)
    df1$ctype = factor(sce[[colorvar]])
    print(class(df1$ctype))
    print(head(df1))
    bas = ggplot(df1, aes(x = PC1, y = PC2, colour = ctype)) + 
        geom_point(alpha = 0.4)
    df2 = data.frame(v1 = shr * fac * kprot[, 1], v2 = shr * 
        fac * kprot[, 2])
    df2$tag = rownames(kprot)
    df2$ctype = NA
    seg = geom_segment(data = df2, aes(x = 0, y = 0, xend = v1, 
        yend = v2), colour = muted("red"), alpha = 0.4, arrow = arrow(length = unit(0.1, 
        "inches")))
    bas + geom_text(data = df2, aes(x = v1, y = v2, label = tag)) + 
        seg
}
