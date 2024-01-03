max_z <- function(data, timepoint1, timepoint2)
{
    maxZ <- apply(data[, c(timepoint1, timepoint2)], 1, max)

    names(maxZ) <- rownames(maxZ)

    return(maxZ)
}
