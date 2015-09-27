

#' @title Computes a Dialectometrical measure
#'
#' @param data object coercible to data.frame containing the data.
#' @param formula formula indicating which variables represent ids and variables, see details.
#' @param value.var character of length 1 indicating the variable containing.
#' @param measure character of length 1 indiciating the selected.
#' @param subset index indicating which subset of data to take.
#' 
#' @export
#'
#' @useDynLib diaMeasures diaMeasure_C
#' 
#' @examples
#' 
#' data(dsample)
#' measure <- diaMeasure(simpleSample, gender + location ~ question, 'answer', 'ird')
#' print(measure)
#'
diaMeasure <- function(data, formula, value.var, measure = c('ird', 'ipi'),
                       binaryIndex = c('dice', 'jaccard'), percentage = TRUE,
                       subset){

    mf <- match.call()
    m <- match(c('data', 'formula', 'value.var', 'subset'), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(data.table::dcast.data.table)
    mf$fun.agg <- quote(function(x) list(as.character(x)))
    if (!data.table::is.data.table(data)){
        mf$data <- substitute(data.table::as.data.table(data), list(data = mf$data))
    }
    mf <- as.list(eval(mf, parent.frame()))
    print(as.data.table(do.call(rbind, mf)))
    availableMethods <- c('ird', 'ipi')
    idnbr <- length(all.vars(formula[[2]]))
    measure <- match(match.arg(measure), availableMethods)
    binaryIndex <- match(match.arg(binaryIndex), c('dice', 'jaccard'))
    attrs <- list(Size = length(mf[[1]]), Labels = do.call(paste, mf[1:idnbr]),
                  class = 'diaMeasure', idVars = mf[1:idnbr], percentage = percentage)
    mf <- mf[-(1:idnbr)]
    result <- .Call(diaMeasure_C, mf, measure, binaryIndex, attrs)
    
    return (result)
}

#' Coertion to matrix
#'
#' Coerce the diaMeasure object into a  matrix class. Uses the as.matrix method from the package dist.
#' 
as.matrix.diaMeasure <- function(x){
    class(x) <- 'dist'
    y <- as.matrix(x)
    diag(y) <- attr(x, "diagv")
    return (y)
}

#' Prints diaMeasure
#'
#' Print lower or upper part of the coerced matrix. 
#' 
print.diaMeasure <- function(x, lower = TRUE, justify = 'none', digits = getOption('digits')){
    y <- as.matrix(x)
    z <- format(y, justify = 'none', digits = digits)
    if (lower){
        z[row(z) < col(z)] <- ''
    } else {
        z[row(z) > col(z)] <- ''
    }
    print(noquote(z))
    return (invisible(x))
}


