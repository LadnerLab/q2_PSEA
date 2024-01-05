delta_by_spline <- function(timepoint1, timepoint2)
{
    SS <- smooth.spline(timepoint1, timepoint2)
    deltaZ <- timepoint2 - predict(SS, timepoint1)$y

    names(deltaZ) <- rownames(deltaZ)

    # return(list(maxZ,deltaZ,SS$x,SS$y))
    return(deltaZ)
}
