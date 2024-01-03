delta_by_spline <- function(data, timepoint1, timepoint2)
{
    SS <- smooth.spline(data[, timepoint1], data[, timepoint2])
    deltaZ <- data[, timepoint2] - predict(SS, data[, timepoint1])$y

    names(deltaZ) <- rownames(deltaZ)

    # return(list(maxZ,deltaZ,SS$x,SS$y))
    return(deltaZ)
}
