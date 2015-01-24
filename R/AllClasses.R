#' DTrack class
#'
#' This class contains Gviz's DataTrack and adds default values to the genome and
#' chromosome slot
#'
#' Refer to \code{DataTrack} for details regarding the constructor.
#'
#' @seealso \code{\link{DataTrack}}, \code{\link{GdObject}}
#'
#' @examples
#' dTrack <- DTrack(start=seq(1,1000, len=100), width=10, data=matrix(runif(400),
#'  nrow=4), name="random data")
#'
#' @name DTrack
#' @aliases DTrack-class
## @importClassesFrom Gviz DataTrack NumericTrack RangeTrack GdObject
#' @author Renan Sauteraud
#' @export
setClass("DTrack", contains = "DataTrack", representation = representation())

setMethod("initialize", "DTrack", function(.Object, ...){
    .Object <- callNextMethod() #Call the initialize function of Gviz::DataTrack
	return(.Object)
})

#' @rdname DTrack
#' @param range,start,end,width,data,name,... Arguments to be passed to
#'   \code{DataTrack}.
#' @export
DTrack <- function(range=NULL, start=NULL, end=NULL, width=NULL, data,
                   name="DTrack", ...){
	GvizDataTrack<-DataTrack(range=range, start=start, end=end, width=width,
			name=name, chromosome="chr0", strand="*",  genome="all",
			data=data, ...)
	#Create DTrack object using the slots of the DataTrack
	new("DTrack", range=GvizDataTrack@range,
			data=GvizDataTrack@data,
			chromosome=GvizDataTrack@chromosome,
			genome=GvizDataTrack@genome,
			strand=GvizDataTrack@strand,
			name=GvizDataTrack@name,
			...)
}


#' ATrack class
#'
#' This class contains Gviz's AnnotationTrack and adds default values to the genome and
#' chromosome slot
#'
#' @seealso \code{\link{AnnotationTrack}}, \code{\link{GdObject}}
#'
#' @examples
#' # Object construction
#' aTrack <- ATrack(start = c(20, 60), end = c(40, 100), name = "random.anno",
#' id=c("small","big"))
#' #Stacking example
#' a2Track <- ATrack(start = c(20, 30), end = c(40, 100), name = "stacking=dense",
#' id = c("small", "big"), stacking = "dense", fill=c("black", "orange"))
#' a3Track <- ATrack(start = c(20, 30), end = c(40, 100), name = "no stacking",
#' id = c("small", "big"), fill = c("black", "orange"))
#' #Plotting
#' plotTracks(trackList = c(aTrack, a2Track, a3Track), showFeatureId = TRUE)

#'
#' @name ATrack
#' @aliases ATrack-class
## @importClassesFrom Gviz AnnotationTrack StackedTrack RangeTrack GdObject
#' @author Renan Sauteraud
#' @export
setClass("ATrack", contains="AnnotationTrack", representation=representation())

setMethod("initialize", "ATrack", function(.Object, ...){
			.Object <- callNextMethod() #Call the initialize function of Gviz::AnnotationTrack
			return(.Object)
		})

#' @rdname ATrack
#' @param range,start,end,width,group,id,stacking,name,fun,selectFun,...
#'   Arguments to be passed to \code{AnnotationTrack}.
#' @export
ATrack <- function(range=NULL, start=NULL, end=NULL, width=NULL, group,
		id, stacking="squish", name="ATrack", fun, selectFun, ...){
	#Create an AnnotationTrack object
	GvizAnnotationTrack<-AnnotationTrack(range=range, start=start, end=end, width=width,
			id=id, name=name, chromosome="chr0", genome="all", stacking=stacking,
			fun=fun, selectFun=selectFun, ...)
	#Create ATrack using the AnnotationTrack
	new("ATrack", range=GvizAnnotationTrack@range,
			chromosome=GvizAnnotationTrack@chromosome,
			genome=GvizAnnotationTrack@genome,
			name=GvizAnnotationTrack@name,
			stacking=GvizAnnotationTrack@stacking,
			...)
}


#' ProbeTrack
#'
#' This track can be used to display the frequency of antibody binding for each
#' probe on an arrayas predicted by \code{pepStat}'s function \code{makeCalls}.
#'
#' @slot sequence A \code{character} vector. The probes sequence.
#' @slot probeStart A \code{numeric} vector. The start position of the probes.
#' @slot intensity A \code{numeric} vector. The frequency of response of each probe.
#'  Or the baseline corrected intensity of the signal.
#'
#' @seealso \code{\link{GdObject}}
#'
#' @inheritParams DTrack
#'
#' @examples
#' if(require(pepDat)){
#'   data(restab)
#'   pt <- ProbeTrack(sequence = restab$peptide,
#'                    intensity = restab$group2,
#'                    probeStart = restab$start)
#'   plotTracks(pt)
#'   plotTracks(pt, from = 460, to = 560, legend=TRUE)
#' }
#'
#' @name ProbeTrack
#' @aliases ProbeTrack-class
#' @author Renan Sauteraud
## @importClassesFrom Gviz GdObject
#' @export
setClass("ProbeTrack", contains = "GdObject",
		representation(sequence="list", intensity = "list",	probeStart = "list"),
		prototype(
				sequence=list(),
				intensity = list(),
				probeStart = list(),
				dp = DisplayPars(
          color=rev(heat.colors(8, alpha = 1)),
					size = 2, #it displays better when slightly bigger than annotation tracks
          cex=1
          )
        )
    )

setMethod("initialize", "ProbeTrack", function(.Object, sequence, intensity, probeStart, name, ...){
			.makeParMapping()
			.Object@name <- name
			.Object@sequence <- sequence
			.Object@intensity <- intensity
			.Object@probeStart <- probeStart
			.Object <- callNextMethod(.Object, ...) #Call the initialize function of Gviz::GdObject
			return(.Object)
      })

#' @rdname ProbeTrack
#' @param sequence A \code{character} vector. The sequence of peptides to
#'   display.
#' @param intensity A \code{numeric} vector. The frequency of binding or the
#'   baseline corrected intensity for the peptides.
#' @param probeStart A \code{numeric} vector. The start position of the peptides.
#' @param name A \code{character}. The name of the track used in the title panel
#'   when plotting
#' @param restab A \code{data.frame} containing all the above parameters, as
#'   outputted by \code{pepStat}'s \code{restab} function.
#' @param group A \code{character}. The group to display on the \code{ProbeTrak}.
#'   This is only required when \code{restab} is not NULL. See details section
#'   for more information.
#'
#' @details
#' The vectors for the arguments \code{sequence}, \code{freq} and
#' \code{probeStart} should be of the same length. If \code{restab} is provided,
#' the three previous arguments will be ignored and \code{group} must be
#' specified. \code{group} must be a valid column name in \code{restab},
#' \code{data.frame}.
#'
#' @seealso \code{restab}
#'
#' @export
ProbeTrack <- function(sequence, intensity, probeStart, restab = NULL,
                       group = NULL, name="ProbeTrack", ...){
  if(!is.null(restab) && class(restab) == "data.frame"){
    if(is.null(group) || !(group %in% colnames(restab))){
      stop("'group' has to be a valid column name in 'restab'")
    }
    sequence <- list(restab$names)
    probeStart <- list(restab$start)
    intensity <- list(restab[, group])
  } else{
    if(missing(probeStart)){
      stop("Need probeStart argument to know where to plot the data on the genome")
    }
    if(missing(sequence)){
      stop("Need sequence argument to know what to plot")
    }
    if(class(sequence)!="list"){
		  sequence<-list(sequence)
		  intensity<-list(intensity)
		  probeStart<-list(probeStart)
    }
  }
	##check the consistancy of number of entires
	if(!(identical(lapply(intensity,length),lapply(probeStart,length)) && identical(lapply(sequence,length),lapply(probeStart,length))))
		stop("sequence, intensity and probeStart need have identifcal structure!")
	##check the types
	if(any(unlist(lapply(sequence,class))!="character"))
		stop("sequence has to be character vector!")
	if(any(unlist(lapply(intensity,class))!="numeric")||any(!unlist(lapply(probeStart,class))%in%c("numeric","integer")))
		stop("intensity and probeStart have to be numeric vectors!")
	new("ProbeTrack", sequence=sequence, intensity=intensity, probeStart=probeStart, name=name, ...)
}


#' ProteinSequenceTrack
#'
#' A track to display peptides and protein sequences.
#'
#' @seealso \code{\link{SequenceTrack}}, \code{\link{DisplayPars}}
#'
#' @author Renan Sauteraud
#'
#' @name ProteinSequenceTrack
#' @aliases ProteinSequenceTrack-class
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom IRanges IRanges
#' @importFrom biovizBase getBioColor
#' @importFrom Biostrings AA_ALPHABET
## @importClassesFrom Gviz GdObject
#' @export
#'
setClass("ProteinSequenceTrack",
		contains = "GdObject",
		representation = representation(sequence = "character",
		                                chromosome = "character", genome = "character"),
    prototype = prototype(name = "Sequence", chromosome = "chr0", genome = "all")
)

#' @rdname ProteinSequenceTrack
#' @param sequence A \code{character} or \code{AAString} of length one. The
#'   sequence to display.
#' @param name A \code{character}. The name of the track used in the title panel
#'   when plotting
#' @param ... Additional items which will all be interpreted as display
#'   parameters.
#'
#' @examples
#' if(require(pepDat)){
#'   data(pep_hxb2)
#'   hxb2_seq <- metadata(pep_hxb2)$sequence
#'   st<-ProteinSequenceTrack(sequence=hxb2_seq, name="env")
#'
#'   # Plotting amino acids
#'   plotTracks(st, to = 20)
#'
#'   # When the range becomes wider, only coloured squares are displayed
#'   plotTracks(st, to = 100)
#'
#'   # When overplotting, a single line will mark the ProteinSequenceTrack
#'   plotTracks(st)
#' }
#'
#' @export
ProteinSequenceTrack <- function(sequence=NULL, name = "Sequence", ...){
  if(is.null(sequence)){
    stop("A 'character' or  'AAString' object should be supplied")
  }
  sequence<-as.character(sequence)
  pst <- new("ProteinSequenceTrack", chromosome = "chr0", genome = "all", name = name, ...)
  displayPars(pst)$fontcolor <- getBioColor(type="AA_ALPHABET")
  displayPars(pst)$size <- 0.25
  pst@sequence <- sequence
  return(pst)
}

#' ProteinAxisTrack
#'
#' A track to display an axis for protein or peptide sequences
#'
#' @seealso \code{\link{GenomeAxisTrack}}
#'
#' @author Renan Sauteraud
#'
#' @name ProteinAxisTrack
#' @aliases ProteinAxisTrack-class
#' @export
#'
setClass("ProteinAxisTrack",
		contains = "GenomeAxisTrack",
		representation(addNC = "logical"),
		prototype(addNC = FALSE, dp = DisplayPars(lwd = 2))
)

setMethod("initialize", "ProteinAxisTrack", function(.Object, addNC, name, ...){
	.makeParMapping()
	.Object@addNC<-addNC
	.Object@name<-name
	.Object<-callNextMethod(.Object, ...)
})

#' @rdname ProteinAxisTrack
#' @param range,name,id,... Arguments to be passed to \code{GenomeAxisTrack}.
#' @param addNC A \code{logical}. If TRUE, display the Amino-terminal and
#'   Carboxyl-terminal ends on the axis.
#'
#' @examples
#' # Object construction
#' paxTrack <- ProteinAxisTrack()
#' pax2 <- ProteinAxisTrack(addNC=TRUE)
#' pax3 <- ProteinAxisTrack(littleTicks=TRUE)
#' # Plotting
#' plotTracks(c(paxTrack,pax2,pax3), from=1, to=100)
#'
#' @export
ProteinAxisTrack<-function(range=NULL, name="Axis", addNC=FALSE, id=NULL, ...){
	GvizGAT<-new("GenomeAxisTrack", name=name, range=range, id=id, ...)
	new("ProteinAxisTrack",range=GvizGAT@range,
			name=GvizGAT@name,
			id=GvizGAT@id,
			addNC=addNC, ...)
}
