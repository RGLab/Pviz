# data.table data.table ':=' '.SD'

#' @import methods
#' @import data.table
NULL

.inter_group <- function(data){
  calls <- data.table(data[, colnames(data)[!colnames(data) %in%
                                              c("peptide", "space", "start", 
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
#' @importFrom RColorBrewer brewer.pal
#' @export
#' 
plot_inter <- function(restab, from = 0, to = max(restab$position), ...){
  annotations <- data.frame(start = c(131L, 157L, 296L, 385L, 459L, 661L, 683L),
                            end = c(156L, 195L, 330L, 417L, 469L, 682L, 705L),
                            width = c(26L, 39L, 35L, 33L, 11L, 22L, 23L),
                            feature = c("V1", "V2", "V3", "V4", "V5", "MPER", "TM"))
  if(from > to){
    stop(paste0("'from' (", from, ") is bigger than 'to' (", to, ")"))
  }
  data <- .inter_group(restab)
  pat <- ProteinAxisTrack(addNC = TRUE)
  at <- ATrack(start = annotations$start, end = annotations$end, id = annotations$feature,
               fill=brewer.pal(n=nrow(annotations), "Paired"), fontcolor.feature = "black",
               background.title="darkgray", name = "Features")
  dt <- DTrack(start=data$position, end=data$position, data=data[, 2:ncol(data), with=FALSE],
               groups =  colnames(data)[2:ncol(data)], name="Freq", legend=TRUE, type="l", 
               background.title="darkgray")
  plotTracks(c(pat, at, dt), from = from, to = to, showFeatureId = TRUE, ...)
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
  uclade <- unique(restab$clade)
  if(length(grep(",", uclade)) > 0){
    warning("It looks like some peptides belong to multiple clades. Run 
            slidingMean with 'split.by.clade=TRUE' to have a single row for each
            peptide/clade.")
  }
  if(!all(clade %in% uclade)){
    stop(paste("All clades must be available in the result table.
         Clade(s)", paste(clade[!clade %in% uclade], collapse = ","), "invalid."))
  }
  pat <- ProteinAxisTrack(addNC = TRUE)
  annotations <- data.frame(start = c(131L, 157L, 296L, 385L, 459L, 661L, 683L),
                            end = c(156L, 195L, 330L, 417L, 469L, 682L, 705L),
                            width = c(26L, 39L, 35L, 33L, 11L, 22L, 23L),
                            feature = c("V1", "V2", "V3", "V4", "V5", "MPER", "TM"))
  at <- ATrack(start = annotations$start, end = annotations$end, id = annotations$feature,
               fill=brewer.pal(n=nrow(annotations), "Paired"), fontcolor.feature = "black",
               background.title="darkgray", name = "Features")
  tracks <- c(pat, at)
  #Clades
  cn <- c("start", "end", "width", "position", "peptide", "space", "clade")
  for(cur_clade in clade){
      data <- restab[ restab$clade == cur_clade, ]
      groups <- colnames(data)[! colnames(data) %in% cn]
      if(cur_clade == tail(clade, 1)){
        legend <- TRUE
      } else{
        legend <- FALSE
      }
      tracks <- c(tracks, DTrack(start = data$position, end = data$position,
                   data = data[, groups], groups = groups, name = cur_clade, 
                   legend = legend, type = "l", background.title = "darkgray", 
                   lwd = 2))
  }
  plotTracks(tracks, from = from, to = to, showFeatureId = TRUE, ...)
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

setMethod("initialize", "CladeTrack", function(.Object, clade, ...){
  .Object@clade <- clade
  .Object <- callNextMethod() #Call the initialize function of Pviz::DTrack
  return(.Object)
})

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
  mat <- as.matrix(data[, groups])
  DT <- DTrack(range = NULL, start = data$position, end = data$position,
                  groups = groups, data = t(mat), ...)
#  new("CladeTrack", clade = clade, range = DT@range, 
#      groups = DT@dp@pars$groups, data = DT@data,
#      chromosome = DT@chromosome, genome = DT@genome, ...)
}
