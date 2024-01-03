psea <- function(maxZ, deltaZ, threshold, input)
{
    library(clusterProfiler)
    input <- input[order(input$gene), , drop=FALSE]
    gene_list <- sort(
        deltaZ[intersect(which(maxZ > threshold), which(deltaZ != 0))],
        decreasing=TRUE
    )

    term_to_gene <- input[, c("term", "gene")]
    out=GSEA(
        geneList=gene_list,
        TERM2GENE=term_to_gene,
        pvalueCutoff=1,
        minGSSize=3,
        maxGSSize=5000,
        eps=1e-30,
        verbose=FALSE,
        nPermSimple=10000,
        exponent=1
    )
    outtable_pre <- attributes(out)[[1]][,c(
        "ID", "enrichmentScore", "NES", "p.adjust",
        "core_enrichment", "pvalue", "qvalue"
    )]
    allpeptides=unlist(lapply(
        lapply(
            attributes(out)$geneSets,
            function(X) intersect(X, names(which(maxZ > threshold)))
        ),
        function(X) paste(X,collapse="/")
    ))
    alltestedpeptides <- allpeptides[match(
        row.names(outtable_pre),names(allpeptides)
    )]
    outtable <- cbind(outtable_pre, alltestedpeptides)
    speciesName <- S[match(outtable[, "ID"], S[, 2]), 1]
    outtable <- cbind(outtable, speciesName)
    return(outtable)
}
