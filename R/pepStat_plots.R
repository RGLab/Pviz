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
#' @param restab A \code{data.frame}. The result of a peptide microarray
#'   analysis, as returned by \code{pepStat}'s \code{restab} function.
#' @param sequence A \code{character} or an \code{AAString}. If not NULL, the
#'   sequence of the \code{ProteinSequenceTrack} to plot. It should be the
#'   sequence of the reference genome used in the \code{peptideSet} that
#'   generated the \code{restab}.
#' @param from A \code{numeric}, the start coordinate of the plot.
#' @param to A \code{numeric}, the end coordinate of the plot.
#' @param ... Aditional arguments to be passed to \code{plotTracks}.
#'
#' @examples
#' if(require(pepDat)){
#'   data(restab_aggregate)
#'   plot_inter(restab_aggregate)
#' }
#'
#' @seealso \code{restab}, \code{plot_clade}, \code{\link{plotTracks}}
#' @author Renan Sauteraud
#' @export
#'
plot_inter <- function(restab, sequence = NULL, from = 0, to = max(restab$position), ...){
  if(from > to){
    stop(paste0("'from' (", from, ") is bigger than 'to' (", to, ")"))
  }
  data <- .inter_group(restab)
  tracks <- ProteinAxisTrack(addNC = TRUE)
  if(!is.null(sequence)){
    tracks <- c(tracks, ProteinSequenceTrack(sequence, name = "",
                                             background.title = "transparent"))
  }
  annotations <- data.frame(start = c(131L, 157L, 296L, 385L, 459L, 661L, 683L),
                            end = c(156L, 195L, 330L, 417L, 469L, 682L, 705L),
                            width = c(26L, 39L, 35L, 33L, 11L, 22L, 23L),
                            feature = c("V1", "V2", "V3", "V4", "V5", "MPER", "TM"))
  tracks <- c(tracks, ATrack(start = annotations$start, end = annotations$end,
               id = annotations$feature, fontcolor.feature = "black",
               fill = rep_len(c("darkgray", "lightgray"), length.out=nrow(annotations)),
               background.title="darkgray",
               name = "Features", rotation.item = 45))
  tracks <- c(tracks, DTrack(start=data$position, end=data$position, data=data[, 2:ncol(data), with=FALSE],
               groups =  colnames(data)[2:ncol(data)], name="Freq", legend=TRUE, type="l",
               background.title="darkgray"))
  plotTracks(tracks, from = from, to = to, showFeatureId = TRUE, ...)
}

#' Plot frequency of response for a single clade.
#'
#' Plot an axis and the frequency of response of a single selected clade.
#'
#' @param restab A \code{data.frame}. The result of a peptide microarray analysis,
#'   as returned by \code{pepStat}'s \code{restab} function.
#' @param clade A \code{character}. The clade to plot.
#' @param sequence An optional \code{character} or \code{AAString}. The
#'   sequence of the \code{ProteinSequenceTrack} to plot. It should be the
#'   sequence of the reference genome used in the \code{peptideSet} that
#'   generated the \code{restab}.
#' @param from A \code{numeric}, the start coordinate of the plot.
#' @param to A \code{numeric}, the end coordinate of the plot.
#' @param ... Aditional arguments to be passed to \code{plotTracks}.
#'
#' @examples
#' if(require(pepDat)){
#'   data(restab)
#'   plot_clade(restab, clade = c("A", "M"))
#' }
#'
#' @seealso \code{restab}, \code{plot_inter}, \code{\link{plotTracks}}
#' @author Renan Sauteraud
#' @export
plot_clade <- function(restab, clade, sequence = NULL, from = 0, to = max(restab$position), ...){
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
  tracks <- ProteinAxisTrack(addNC = TRUE)
  if(!is.null(sequence)){
    tracks <- c(tracks, ProteinSequenceTrack(sequence, name = "",
                                             background.title = "transparent"))
  }
  annotations <- data.frame(start = c(131L, 157L, 296L, 385L, 459L, 661L, 683L),
                            end = c(156L, 195L, 330L, 417L, 469L, 682L, 705L),
                            width = c(26L, 39L, 35L, 33L, 11L, 22L, 23L),
                            feature = c("V1", "V2", "V3", "V4", "V5", "MPER", "TM"))
  tracks <- c(tracks, ATrack(start = annotations$start, end = annotations$end,
                      id = annotations$feature, fontcolor.feature = "black",
                      fill = rep_len(c("darkgray", "lightgray"), length.out=nrow(annotations)),
                      background.title="darkgray", name = "Features",  rotation.item = 45))
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
#' @param name A \code{character}. The name of the track, used in the title
#'   panel when plotting. By default, the \code{clade} name.
#' @param ... Additional argument to be passed to \code{DataTrack}. They will be
#'   treated as display parameters.
#'
#' @examples
#' if(require(pepDat)){
#'   data(restab)
#'   ct <- CladeTrack(restab, clade = "M", type = "l", legend = TRUE)
#'   plotTracks(ct)
#' }

#'
#' @export
CladeTrack <- function(restab, clade, name = clade, ...){
  data <- restab[ restab$clade == clade, ]
  cn <- c("peptide", "start", "end", "width", "position", "names", "space", "clade")
  groups <- colnames(data)[! colnames(data) %in% cn]
  mat <- as.matrix(data[, groups])
  DT <- DTrack(range = NULL, start = data$position, end = data$position,
                  groups = groups, data = t(mat), name = name, ...)
  new("CladeTrack", clade = clade, range = DT@range,
      groups = DT@dp@pars$groups, data = DT@data,
      chromosome = DT@chromosome, genome = DT@genome,
      name = DT@name, ...)
}
