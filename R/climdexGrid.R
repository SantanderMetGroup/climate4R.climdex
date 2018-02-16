#     climdexGrid.R ETCCDI Climate Changes Indices in Climate4R
#
#     Copyright (C) 2018 Santander Meteorology Group (http://www.meteo.unican.es)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program. If not, see <http://www.gnu.org/licenses/>.

#' @title ETCCDI Climate Changes Indices in Climate4R
#' @description Calculation of the 27 core indices of the Expert Team on Climate Change Detection and Indices (ETCCDI).
#' The function is a wrapper of the package \pkg{climdex.pcic} for its seamless integration with climate4R objects.
#' @param tn A climate4R grid of daily minimum temperature
#' @param tx A climate4R grid of daily maximum temperature
#' @param pr A climate4R grid of daily precipitation
#' @param cal A calendar definition. Default to 365-day calendar. This argument is passed to \code{\link[PCICt]{as.PCICt}},
#' whose help documentation contains further details.
#' @param index.code Character string, indicating the specific code of the index according to the ETCCDI
#' definitions (see Details).
#' @param input.arg.list Optional. A list of arguments internally passed to \code{\link[climdex.pcic]{climdexInput.raw}}
#'  from package \pkg{climdex.pcic}
#' @param index.arg.list Optional (but depending on the specific index might be necessary). A list of specific arguments
#'  for the target index. See Details.
#' @template templateParallelParams
#' @importFrom transformeR getTimeResolution redim checkDim getShape getRefDates mat2Dto3Darray
#'  array3Dto2Dmat subsetGrid parallelCheck selectPar.pplyFun getYearsAsINDEX getCoordinates aggregateGrid
#' @importFrom parallel stopCluster
#' @importFrom magrittr %>% %<>%
#' @importFrom PCICt as.PCICt
#' @importFrom utils head
#' @import climdex.pcic
#' @details \code{\link{climdexShow}} will display on screen a full list of ETCCDI Core indices and their codes. The names of the
#' internal functions calculating each index is also displayed, whose help files can aid in the definition of index-specific arguments.
#'
#' \strong{Parallel computing}
#'
#' @template templateParallel
#'
#' @author J. Bedia
#' @export

climdexGrid <- function(index.code,
                        tn = NULL,
                        tx = NULL,
                        pr = NULL,
                        input.arg.list = list(),
                        index.arg.list = list(),
                        cal = "365_day",
                        parallel = FALSE,
                        max.ncores = 16,
                        ncores = NULL) {
    if (!is.null(tn)) {
        if (getTimeResolution(tn) != "DD") stop("Daily data is required as input", call. = FALSE)
    }
    if (!is.null(tx)) {
        if (getTimeResolution(tx) != "DD") stop("Daily data is required as input", call. = FALSE)
    }
    if (!is.null(pr)) {
        if (getTimeResolution(pr) != "DD") stop("Daily data is required as input", call. = FALSE)
    }
    index.code <- match.arg(index.code,
                            choices = c("FD", "SU", "ID", "TR", "GSL", "TXx",
                                        "TNx", "TXn", "TNn", "TN10p", "TX10p",
                                        "TN90p", "TX90p", "WSDI", "CSDI", "DTR",
                                        "Rx1day", "Rx5day", "SDII", "R10mm",
                                        "R20mm", "Rnnmm", "CDD", "CWD",
                                        "R95pTOT", "R99pTOT", "PRCPTOT"))
    aux <- read.master()
    metadata <- aux[grep(index.code, aux$code, fixed = TRUE), ]
    a <- c(!is.null(tn), !is.null(tx), !is.null(pr)) %>% as.numeric()
    b <- metadata[ , 4:6] %>% as.numeric()
    if (any(b - a > 0)) {
        stop("The required input variable(s) for ", index.code, " index calculation are missing\nType \'?", metadata$indexfun, "\' for help",call. = FALSE)
    }
    # Remove any possible uneeded input grid
    if (any(a - b > 0)) {
        ind <- which((a - b) > 0)
        rem <- c("tn", "tx", "pr")[ind]
        sapply(rem, function(x) assign(x, NULL)) %>% invisible()
        message("NOTE: some input grids provided are not required for ", index.code, " calculation and were removed")
    }
    # Ensure member is present data structures
    if (!is.null(tx)) tx %<>% redim(member = FALSE, var = FALSE)
    if (!is.null(tn)) tn %<>% redim(member = FALSE, var = FALSE)
    if (!is.null(pr)) pr %<>% redim(member = FALSE, var = FALSE)
    # Check structural consistency of data arrays when multiple
    if (index.code == "DTR") sapply(list(tn, tx), "checkDim", dimensions = c("member", "time", "lat", "lon")) %>% invisible()
    # name of the reference grid
    refGridName <- c("tn","tx","pr")[which(c(!is.null(tn), !is.null(tx), !is.null(pr)) %>% as.numeric() != 0)] %>% head(1)
    assign("refGrid", get(refGridName))
    # Number of members
    # n.mem <- getShape(refGrid, "member")
    coords <- expand.grid(refGrid$xyCoords$y, refGrid$xyCoords$x)[2:1]
    # coercion to PCICt
    refDates <- getRefDates(refGrid) %>% as.PCICt(cal = cal)
    refDates.tx <- refDates.tn <- refDates.pr <- NULL
    if (!is.null(tx)) refDates.tx <- refDates
    if (!is.null(tn)) refDates.tn <- refDates
    if (!is.null(pr)) refDates.pr <- refDates
    refDates <- NULL
    message("[", Sys.time(), "] Calculating ", index.code, " ...")
    # skip masked grid points
    rm.ind.tx <- rm.ind.tn <- rm.ind.pr <- c()
    aux.tx <- aux.tn <- aux.pr <- NULL
    if (!is.null(tx)) {
        aux.tx <- tx[["Data"]] %>% array3Dto2Dmat()
        rm.ind.tx <- which(apply(aux.tx, MARGIN = 2, FUN = function(x) all(is.na(x))))
    }
    if (!is.null(tn)) {
        aux.tn <- tn[["Data"]] %>% array3Dto2Dmat()
        rm.ind.tn <- which(apply(aux.tn, MARGIN = 2, FUN = function(x) all(is.na(x))))
    }
    if (!is.null(pr)) {
        aux.pr <- pr[["Data"]] %>% array3Dto2Dmat()
        rm.ind.pr <- which(apply(aux.pr, MARGIN = 2, FUN = function(x) all(is.na(x))))
    }
    # Remove missing values
    rm.ind <- Reduce(union, list(rm.ind.tx, rm.ind.tn, rm.ind.pr))
    if (length(rm.ind) > 0) {
        if (!is.null(tx)) {
            aux.tx <- aux.tx[, -rm.ind, drop = FALSE]
        }
        if (!is.null(tn)) {
            aux.tn <- aux.tn[, -rm.ind, drop = FALSE]
        }
        if (!is.null(pr)) {
            aux.pr <- aux.pr[, -rm.ind, drop = FALSE]
        }
    }
    nvalid.points <- nrow(coords) - length(rm.ind)
    valid.coords <- coords[-rm.ind, ]
    parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
    apply_fun <- selectPar.pplyFun(parallel.pars, .pplyFUN = "sapply")
    if (parallel.pars$hasparallel) on.exit(parallel::stopCluster(parallel.pars$cl))
    out <- apply_fun(1:nvalid.points, function(i) {
        input.arg.list[["northern.hemisphere"]] <- ifelse(valid.coords[i,2] >= 0, TRUE, FALSE)
        input.arg.list[["tmax"]] <- aux.tx[ , i]
        input.arg.list[["tmin"]] <- aux.tn[ , i]
        input.arg.list[["prec"]] <- aux.pr[ , i]
        input.arg.list[["tmax.dates"]] <- refDates.tx
        input.arg.list[["tmin.dates"]] <- refDates.tn
        input.arg.list[["prec.dates"]] <- refDates.pr
        # Define available temporal range as baseline by default
        if (!"base.range" %in% names(input.arg.list)) {
            input.arg.list[["base.range"]] <- getYearsAsINDEX(refGrid) %>% range() %>% as.integer()
        }
        ci <- do.call("climdexInput.raw", input.arg.list)
        index.arg.list[["ci"]] <- ci
        do.call(metadata$indexfun, index.arg.list)
    })
    # Recover original matrix with masked points
    aux <- matrix(NA, nrow = nrow(out), ncol = nrow(coords))
    aux[ ,setdiff(1:ncol(aux), rm.ind)] <- out
    out <- NULL
    # Transform to climate4R grid
    refGrid <- suppressMessages(aggregateGrid(refGrid, aggr.m = list(FUN = "mean")))
    attr(refGrid[["Variable"]], "monthly_agg_cellfun") <- metadata$indexfun
    if (nrow(aux) == getYearsAsINDEX(refGrid) %>% unique() %>% length()) {
        refGrid <- suppressMessages(aggregateGrid(refGrid, aggr.y = list(FUN = "mean")))
        attr(refGrid[["Variable"]], "annual_agg_cellfun") <- metadata$indexfun
    }
    refGrid[["Data"]] <- mat2Dto3Darray(aux, x = getCoordinates(refGrid)$x, y = getCoordinates(refGrid)$y)
    # Attributes
    refGrid[["Variable"]][["varName"]] <- metadata$code
    refGrid[["Variable"]][["varName"]] <- NULL
    attr(refGrid[["Variable"]], "longname") <- metadata$longname
    if (!is.na(metadata$units)) attr(refGrid[["Variable"]], "units") <- metadata$units
    message("[", Sys.time(), "] Done")
    invisible(refGrid)
}




#' @title List all the 27 ETCCDI Core Indices
#' @description Print a table with a summary of the 27 ETCCDI Core indices
#' @return Print a table on the screen with the following columns:
#' \itemize{
#' \item \strong{code}: Code of the index. This is the character string used as input value
#' for the argument \code{index.code} in \code{\link{climdexGrid}}
#' \item \strong{longname}: Long description of the index
#' \item \strong{index.fun}: The name of the internal function from package \pkg{\link{climdex.pcic}} used to calculate it
#' \item \strong{tn,tx,pr}: A logical value (0/1) indicating the input variables required for index calculation
#' \item \strong{units}: The units of the index (when different from those of the input variable)
#' }
#' @references The ETCDDI web page giving the definition of the 27 core indices:
#' \url{http://etccdi.pacificclimate.org/list_27_indices.shtml}
#' @author J. Bedia
#' @export
#' @importFrom magrittr %>%

climdexShow <- function() {
     read.master() %>% print()
}



#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom utils read.table

read.master <- function() {
    system.file("master", package = "climate4R.climdex") %>% read.table(header = TRUE,
                                                                        sep = ";",
                                                                        stringsAsFactors = FALSE,
                                                                        na.strings = "")
}



