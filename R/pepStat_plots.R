# data.table data.table ':=' '.SD'

#' @import methods
#' @import data.table
NULL

.inter_group <- function(data){
  calls <- data.table(data[, colnames(data)[!colnames(data) %in%
                                              c("names", "space", "start", 
                                                "end", "width", "clade")]])
  calls <- calls[, lapply(.SD, mean), by = "position"]
  calls <- calls[order(position)]
  return(calls)
}

##
#' Plot frequency of response for each group
#' 
#' Plot an axis and the frequency of response of each group, averaged by peptides
#' at each position.
#' 
#' @param restab A \code{data.frame}. The result of a peptide microarray analysis,
#'   as returned by \code{pepStat}'s \code{restab} function.
#' @param from A \code{numeric}, the start coordinate of the plot.
#' @param to A \code{numeric}, the end coordinate of the plot.
#' @param ... Aditional arguments to be passed to \code{plotTracks}.
#' 
#' @seealso \code{restab}, \code{plot_clade}, \code{\link{plotTracks}}
#' @author Renan Sauteraud
#' @export
#' 
plot_inter <- function(restab, from = 0, to = max(restab$position), ...){
  if(from > to){
    stop(paste0("'from' (", from, ") is bigger than 'to' (", to, ")"))
  }
  data <- .inter_group(restab)
  pat <- ProteinAxisTrack(addNC = TRUE)
  dt <- DTrack(start=data$position, end=data$position, data=data[, 2:ncol(data), with=FALSE],
               groups =  colnames(data)[2:ncol(data)], name="Freq", legend=TRUE, type="l", background.title="darkgray")
  plotTracks(c(pat, dt), from = from, to = to, ...)
}

#' Plot frequency of response for a single clade.
#' 
#' Plot an axis and the frequency of response of a single selected clade.
#' 
#' @param restab A \code{data.frame}. The result of a peptide microarray analysis,
#'   as returned by \code{pepStat}'s \code{restab} function.
#' @param clade A \code{character}. The clade to plot.
#' @param from A \code{numeric}, the start coordinate of the plot.
#' @param to A \code{numeric}, the end coordinate of the plot.
#' @param ... Aditional arguments to be passed to \code{plotTracks}.
#' 
#' @seealso \code{restab}, \code{plot_inter}, \code{\link{plotTracks}}
#' @author Renan Sauteraud
#' @export
plot_clade <- function(restab, clade, from = 0, to = max(restab$position), ...){
  if(from > to){
    stop(paste0("'from' (", from, ") is bigger than 'to' (", to, ")"))
  }
  if(length(grep(",", unique(restab$clade))) > 0){
    warning("It looks like some peptides belong to multiple clades. Run restab 
            with 'long=TRUE' to have a single row for each peptide/clade.")
  }
  data <- restab[ restab$clade == clade, ]
  cn <- c("start", "end", "width", "position", "names", "space", "clade")
  groups <- colnames(data)[! colnames(data) %in% cn]
  pat <- ProteinAxisTrack(addNC = TRUE)
  dt <- DTrack(start = data$position, end = data$position,
               data = data[, groups], groups = groups, name = clade, 
               legend=TRUE, type = "l", background.title = "darkgray")
  plotTracks(c(pat, dt), from = from, to = to, ...)
}



#' CladeTrack
#' 
#' This track can be used to display the result of pepStat analysis for a single
#' clade. It contains \code{DTrack}.
#' 
#' @slot clade A \code{character}. The clade to display.
#' 
#' @seealso \code{DTrack}
#'
#' @name CladeTrack
#' @aliases CladeTrack-class
#' @author Renan Sauteraud
#' @export
setClass("CladeTrack", contains = "DTrack",
         representation = representation(clade = "character"))

#' @rdname CladeTrack
#' @param restab A \code{data.frame}. The result of a peptide microarray analysis,
#'   as returned by \code{pepStat}'s \code{restab} function.
#' @param clade A \code{character}. The clade to plot.
#' @param ... Additional argument to be passed to \code{DataTrack}. They will be
#'   treated as display parameters.
#' @export
CladeTrack <- function(restab, clade, ...){
  data <- restab[ restab$clade == clade, ]
  cn <- c("start", "end", "width", "position", "names", "space", "clade")
  groups <- colnames(data)[! colnames(data) %in% cn]
  new("CladeTrack", start = restab$positio, end = restab$position,
      data = restab, data = data[, groups], ...)
}

