#' @title Classic Dialektometrical ad-hoc classification
#'
#' Applies certain ad-hoc measures used in Dialektometrics. 
#'
#' The available classification methods are Med MinMwMax and MedMw.
#'
#' @param measure Object coercible to data.frame containing the data.
#' @param method Classification method, see details.
#' @param n Number of groups to be made.
#' @param ids vector of labels indicating which linguistical use as center. All will be used by default.
#'
#' @export
#'
#' @examples
#' 
#' data(dsample)
#' measure <- diaMeasure(dsample, location ~ question, 'answer', 'ird')
#' diaClass <- diaClassification(measure, 'Med', 6)
#'
diaClassification <- function(measure, method = c('Med', 'MinMwMax', 'MedMw'), n, ids)
    UseMethod('diaClassification')
    
#' @export    
diaClassification.diaMeasure <- function(measure, method = c('Med', 'MinMwMax', 'MedMw'), n, ids){
    method <- match.arg(method)
    if (method %in% c('MinMwMax', 'MedMw') & (n %% 2L != 0) )
        stop ('For method ', method, ' the parameter n must be multiple of 2')
    measureMat <- as.matrix(measure)
    
    if (missing(ids)){
        ids <- attr(measure, 'Labels')
        d <- ncol(measureMat) - 1L
        measureMat <- matrix(measureMat[row(measureMat) != col(measureMat)], d, d + 1, FALSE)
    } else {
        ids <- unique(ids)
        outlawIds <- ids[!(ids %in% attr(measure, 'Labels'))]
        if (length(outlawIds) != 0L){
            stop ('ids ', paste(outlawIds, collapse = ', '), ' are not in the "Labels" attribute of measure')
        }
        nrows <- ncol(measureMat) - 1L
        ncols <- length(ids)
        matched <- match(ids, attr(measure, 'Labels'))
        measureMat <- matrix(measureMat[row(measureMat) != col(measureMat) & col(measureMat) %in% matched], nrows, ncols, TRUE)
    }
    allIds <- attr(measure, 'Labels')
    result <- setNames(vector('list', length(ids)), ids)
    if (method == 'Med'){
        for (i in 1L:length(ids)){
            id <- allIds[i]
            breaks <- quantile(measureMat[, i], seq(0L, 1L, by = 1L / n))
            if (anyDuplicated(breaks)){
                warning ('Too few observations for significant classification in Med')
                breaks <- unique(breaks)
            }
            result[[id]] <- cut(measureMat[, i], breaks, labels = FALSE, include.lowest = TRUE)
            names(result[[id]]) <- setdiff(allIds, id)
        }
    } else if (method == 'MinMwMax'){
        for (i in 1L:length(ids)){
            id <- allIds[i]
            x <- measureMat[, i]
            meanx <- mean(x); minx <- min(x); maxx <- max(x)
            lw <- (meanx - minx) / (n / 2L)
            up <- (maxx - meanx) / (n / 2L)
            breaks <- c(1L:(n / 2L) * lw,  meanx + 1L:(n / 2L) * up)
            result[[id]] <- cut(x, breaks, labels = FALSE, include.lowest = TRUE)
            names(result[[id]]) <- setdiff(allIds, id)
        }
    } else { # MedMw
        for (i in 1L:length(ids)){
            id <- allIds[i]
            x <- measureMat[, i]
            splitPoint <- x[which.min(abs(x - mean(x)))]
            lowerBreaks <- quantile(x[x < splitPoint], seq(0, 1L, 2L / n))
            higherBreaks <- quantile(x[x >= splitPoint], seq(0, 1L, 2L / n))
            breaks <- c(lowerBreaks, higherBreaks)
            if (anyDuplicated(breaks)){
                warning ('Too few observations for significant classification in MedMw.')
                breaks <- unique(breaks)
            }
            result[[id]] <- cut(x, breaks, labels = FALSE, include.lowest = TRUE)
            names(result[[id]]) <- setdiff(allIds, id)
        }
    }
    attr(result, 'class') <- 'diaClassification'
    attr(result, 'range') <- 1L:n
    attr(result, 'Labels') <- allIds
    return (result)
}




