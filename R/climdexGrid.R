#     climdexGrid.R ETCCDI Core Indices in Climate4R
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
#' @param tn A climate4R dataset of daily minimum temperature (degrees C)
#' @param tx A climate4R dataset of daily maximum temperature (degrees C)
#' @param pr A climate4R dataset of daily precipitation (mm)
#' @param cal A calendar definition. Default to 365-day calendar. This argument is passed to \code{\link[PCICt]{as.PCICt}},
#' whose help documentation contains further details.
#' @param index.code Character string, indicating the specific code of the index according to the ETCCDI
#' definitions (see Details).
#' @param input.arg.list Optional. A list of arguments internally passed to \code{\link[climdex.pcic]{climdexInput.raw}}
#'  from package \pkg{climdex.pcic}
#' @param index.arg.list Optional (but depending on the specific index might be necessary). A list of specific arguments
#'  for the target index. See Details.
#' @template templateParallelParams
#' @import transformeR
#' @importFrom parallel stopCluster
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom PCICt as.PCICt
#' @importFrom utils head
#' @import climdex.pcic
#' @details \code{\link{climdexShow}} will display on screen a full list of ETCCDI Core indices and their codes. The names of the
#' internal functions calculating each index is also displayed, whose help files can aid in the definition of index-specific arguments.
#'
#' \strong{Baseline period}
#' By default, the function will use as baseline period the full period of years encompassed by the input grid(s). This
#' is used, for instance, for calculating the relevant percentiles in some indices etc. To use a specific baseline period
#' use the \code{"base.range"} argument in the \code{input.arg.list} internally passed to
#'  \code{\link[climdex.pcic]{climdexInput.raw}}.
#'
#' @template templateParallel
#'
#' @examples \dontrun{
#' require(climate4R.climdex)
#' require(visualizeR)
#' data("tasmin.eobs")
#' ## FROST DAYS (Annual count of days when TN < 0 degC)
#' fd.grid <- climdexGrid(tn = tasmin.eobs, index.code = "FD")
#' spatialPlot(climatology(fd.grid), at = seq(0,165,10),
#'             main = "Mean annual number of frost days (1991-2010)")
#'
#'
#' # The following example will compute the CWD index
#' # (maximum number of consecutive days with RR ≥ 1mm)
#' # for the homogeneized VALUE dataset of stations over Europe
#' # (See Gutiérrez et al. 2018 DOI:10.1002/joc.5462 for details on this dataset)
#' library(loadeR)
#' library(transformeR)
#' library(visualizeR)
#' destfile = "/tmp/VALUE_ECA_86_v2.tar.gz"
#' # (~8Mb download, change destfile at your convenience)
#' download.file("http://meteo.unican.es/work/loadeR/data/VALUE_ECA_86_v2.tar.gz", destfile = destfile)
#' untar(destfile, exdir = "/tmp")
#' station.data <- loadStationData(dataset = "/tmp/VALUE_ECA_86_v2",
#'                                 var = "precip",
#'                                 years = 1981:2000)
#' cwd <- climdexGrid(index.code = "CWD", pr = station.data)
#' spatialPlot(climatology(cwd), backdrop.theme = "countries",
#'             main = "Mean number of consecutive annual wet days (1981-2000)")
#' }
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
        stop("The required input variable(s) for ", index.code,
             " index calculation are missing\nType \'?",
             metadata$indexfun, "\' for help", call. = FALSE)
    }
    # Remove any possible uneeded input grid
    if (any(a - b > 0)) {
        ind <- which((a - b) > 0)
        rem <- c("tn", "tx", "pr")[ind]
        sapply(rem, function(x) assign(x, NULL)) %>% invisible()
        message("NOTE: some input grids provided for ", index.code,
                " index calculation are not required and were removed")
    }
    # Ensure member is present in data structures / handle station <--> grid
    if (!is.null(tx)) {
        station <- ifelse(typeofGrid(tx) == "station", TRUE, FALSE)
        tx %<>% redim(member = TRUE, var = FALSE)
    }
    if (!is.null(tn)) {
        station <- ifelse(typeofGrid(tn) == "station", TRUE, FALSE)
        tn %<>% redim(member = TRUE, var = FALSE)
    }
    if (!is.null(pr)) {
        station <- ifelse(typeofGrid(pr) == "station", TRUE, FALSE)
        pr %<>% redim(member = TRUE, var = FALSE)
    }
    # Check structural consistency of data arrays when multiple
    if (index.code == "DTR" | index.code == "GSL") {
        sapply(list(tn, tx), "checkDim", dimensions = c("member", "time", "lat", "lon")) %>% invisible()
    }
    # name of the reference grid
    refGridName <- c("tn","tx","pr")[which(c(!is.null(tn),
                                             !is.null(tx),
                                             !is.null(pr)) %>% as.numeric() != 0)] %>% head(1)
    assign("refGrid", get(refGridName))
    # Number of members
    n.mem <- getShape(refGrid, "member")
    coords <- if (station) {
        getCoordinates(refGrid)
    } else {
        expand.grid(refGrid$xyCoords$y, refGrid$xyCoords$x)[2:1]
    }
    # coercion to PCICt
    refDates <- getRefDates(refGrid) %>% as.PCICt(cal = cal)
    refDates.tx <- refDates.tn <- refDates.pr <- NULL
    if (!is.null(tx)) refDates.tx <- refDates
    if (!is.null(tn)) refDates.tn <- refDates
    if (!is.null(pr)) refDates.pr <- refDates
    refDates <- NULL
    message("[", Sys.time(), "] Calculating ", index.code, " ...")
    if (n.mem > 1) {
        parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
        apply_fun <- selectPar.pplyFun(parallel.pars, .pplyFUN = "lapply")
        if (parallel.pars$hasparallel) on.exit(parallel::stopCluster(parallel.pars$cl))
    } else {
        if (isTRUE(parallel)) message("NOTE: Parallel processing was skipped (unable to parallelize one single member)")
        apply_fun <- lapply
    }
    out.list <- apply_fun(1:n.mem, function(x) {
        # if (n.mem > 1) message("[", Sys.time(), "] Calculating ",
                               # index.code, " for member ", x, " ...")
        rm.ind.tx <- rm.ind.tn <- rm.ind.pr <- c()
        aux.tx <- aux.tn <- aux.pr <- NULL
        if (!is.null(tx)) {
            tmp <- subsetGrid(tx, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            aux.tx <- tmp[["Data"]] %>% array3Dto2Dmat()
            rm.ind.tx <- which(apply(aux.tx, MARGIN = 2, FUN = function(x) all(is.na(x))))
        }
        if (!is.null(tn)) {
            tmp <- subsetGrid(tn, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            aux.tn <- tmp[["Data"]] %>% array3Dto2Dmat()
            rm.ind.tn <- which(apply(aux.tn, MARGIN = 2, FUN = function(x) all(is.na(x))))
        }
        if (!is.null(pr)) {
            tmp <- subsetGrid(pr, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            aux.pr <- tmp[["Data"]] %>% array3Dto2Dmat()
            rm.ind.pr <- which(apply(aux.pr, MARGIN = 2, FUN = function(x) all(is.na(x))))
        }
        # Remove missing values
        rm.ind <- Reduce(union, list(rm.ind.tx, rm.ind.tn, rm.ind.pr))
        valid.coords <- coords
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
            valid.coords <- coords[-rm.ind, ]
        }
        nvalid.points <- nrow(valid.coords)
        out <- sapply(1:nvalid.points, function(i) {
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
        tx <- tn <- pr <- NULL
        # Recover original matrix with masked points
        aux <- matrix(NA, nrow = nrow(out), ncol = nrow(coords))
        aux[ , setdiff(1:ncol(aux), rm.ind)] <- out
        out <- NULL
        # Transform to climate4R grid
        refGrid <- suppressMessages(aggregateGrid(tmp, aggr.m = list(FUN = "mean")))
        tmp <- NULL
        attr(refGrid[["Variable"]], "monthly_agg_cellfun") <- metadata$indexfun
        if (nrow(aux) == getYearsAsINDEX(refGrid) %>% unique() %>% length()) {
            refGrid <- suppressMessages(aggregateGrid(refGrid, aggr.y = list(FUN = "mean")))
            attr(refGrid[["Variable"]], "annual_agg_cellfun") <- metadata$indexfun
        }
        if (station) {
            refGrid[["Data"]] <- aux
            attr(refGrid[["Data"]], "dimensions") <- c("time", "loc")
        } else {
            refGrid[["Data"]] <- mat2Dto3Darray(aux, x = getCoordinates(refGrid)$x, y = getCoordinates(refGrid)$y)
        }
        # Add attributes
        refGrid[["Variable"]][["varName"]] <- metadata$code
        attr(refGrid[["Variable"]], "longname") <- metadata$longname
        attr(refGrid[["Variable"]], "wasDefinedBy") <- "ETCCDI"
        attr(refGrid[["Variable"]], "hasMainURL") <- "http://etccdi.pacificclimate.org/list_27_indices.shtml"
        attr(refGrid[["Variable"]], "description") <- metadata$description
        if (!is.na(metadata$units)) attr(refGrid[["Variable"]], "units") <- metadata$units
        return(refGrid)
    })
    message("[", Sys.time(), "] Done")
    out <- if (length(out.list) == 1) {
        out.list %>% extract2(1) %>% redim(drop = TRUE)
    } else {
        do.call("bindGrid", c(out.list, dimension = "member"))
    }
    if (station) out %<>% redim(drop = FALSE, loc = TRUE, member = FALSE)
    invisible(out)
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
#' @references The ETCCDI web page giving the definition of the 27 core indices:
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



