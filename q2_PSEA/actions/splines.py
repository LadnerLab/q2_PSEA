# Provides a place to collect spline operations.
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage


r_splines = """
cubic_spline <- function(x, y)
{
    library(splines2)

    knots <- summary(x)[c(2, 3, 5)]
    newx <- c(0, 0.27, 3)  # TODO: shoud the user have a chance to change this?
    sorted.x = sort(x)

    bsMat <- bSpline(
        x,
        knots = knots,
        degree = 3,  # TODO: should the user have a chance to change this?
        df = 7,  # TODO: should the user have a chance to change this?
        intercept = TRUE
    )

    cubic_spline_obj <- lm(y ~ bsMat)
    cubic_spline_preds <- predict(
        cubic_spline_obj,
        newdata = list(sorted.x),
        se=TRUE
    )
    cubic_spline_se_bands <- cbind(
        cubic_spline_preds$fit + 2 * cubic_spline_preds$se.fit,
        cubic_spline_preds$fit - 2 * cubic_spline_preds$se.fit
    )
    # cubic_spline_residuals <- y - cubic_spline_preds$fit

    return(cubic_spline_preds$fit)
}
"""

R_SPLINES = SignatureTranslatedAnonymousPackage(r_splines, "internal")
