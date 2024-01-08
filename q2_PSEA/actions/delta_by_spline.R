delta_by_spline <- function(timepoint1, timepoint2)
{
    SS <- smooth.spline(timepoint1, timepoint2)
    deltaZ <- timepoint2 - predict(SS, timepoint1)$y

    names(deltaZ) <- rownames(deltaZ)

    # uncomment if smooth spline values needed in future
    # return(list(deltaZ, SS$x, SS$y))
    return(deltaZ)
}
