\name{ATrack-class}
\Rdversion{1.1}
\docType{class}
\alias{ATrack-class}
\alias{ATrack}
\alias{initialize,ATrack-method}

\title{Class \code{"ATrack"} and methods}
\description{
	A Track to display range annotations.
}
\section{Objects from the Class}{
	Objects can be created using the constructor function \code{ATrack}.
}
\usage{
	ATrack(range=NULL, start=NULL, end=NULL, width=NULL, group,
		id, strand="*", chromosome="chrR", genome="all", stacking="squish", name="ATrack", protein="prot",
		fun, selectFun, ...)
}

\arguments{
  
	\item{range}{An object of class \code{\linkS4class{GRanges}},
    \code{\linkS4class{IRanges}} or a \code{data.frame} which will be
    coerced into one. In the former case, there should be at least the
    following two columns present as \code{elementMetadata}:
		
    \describe{
      \item{}{\code{start}, \code{end}: the start and end coordinates
		for the track items.}
    }
    
    Alternatively, the \code{range} argument may be missing, in which
    case the relevant information has to be provided as individual
    function arguments (see below).
    }
    
    \item{start, end, width}{Integer vectors, giving the start and the end
    end coordinates for the individual track items, or their width. Two
    of the three need to be specified, and have to be of equal length or
    of length one, in which case this value will be recycled. Otherwise,
    the usual R recycling rules for vectors do not apply.}
    
    \item{id}{Character vector of track item identifiers. Needs to be of
    equal length as the \code{start, end} or \code{width} vectors, or of
    length 1.}
    
    \item{stacking}{The stacking type for overlapping items of the
    track. One in \code{c(hide, dense, squish)}.}
    
    \item{name}{Character scalar of the track's name used in the title
    panel when plotting.}
    
    \item{protein}{Name of the protein.}
    
    \item{\dots}{Additional items which will all be interpreted as further
    display parameters.}
    
    \item{chromosome}{Unused argument of class \code{"character"},
  	inherited from class \code{\linkS4class{RangeTrack}},
  	defaulted to "chrR".}

  	\item{strand}{Unused argument of class \code{"character"},
  	inherited from class \code{\linkS4class{AnnotationTrack}},
  	defaulted to "chrR".}

  	\item{genome}{Unused argument of class \code{"character"},
  	inherited from class \code{\linkS4class{RangeTrack}}, 
  	defaulted to "all".}
        
    \item{group}{Factor (or other vector that can be coerced into one),
    giving the group memberships for the individual track items.  Needs
    to be of equal length as the \code{start, end} or \code{width}
    vectors, or of length 1.
    Inherited from class \code{\linkS4class{AnnotationTrack}}.}

  	\item{fun}{Object inherited from\code{\linkS4class{AnnotationTrack}}.}

  	\item{selectFun}{Object inherited from\code{\linkS4class{AnnotationTrack}}.}

}

\value{
  The return value of the constructor function is a new object of class \code{ATrack}.  
}

\section{Slots}{
  \describe{
    \item{\code{protein}:}{Object of class \code{"character"},
      a name for the protein.}
    \item{\code{stacking}:}{Object of class \code{"character"},
      inherited from class \code{\linkS4class{StackedTrack}}}
    \item{\code{stacks}:}{Object of class \code{"numeric"},
      inherited from class \code{\linkS4class{StackedTrack}}}
    \item{\code{range}:}{Object of class \code{"GRangesOrIRanges"},
      inherited from class \code{\linkS4class{RangeTrack}}}
    \item{\code{chromosome}:}{Object of class \code{"character"},
      inherited from class \code{\linkS4class{RangeTrack}} }
    \item{\code{genome}:}{Object of class \code{"character"}, 
      inherited from class \code{\linkS4class{RangeTrack}}}
    \item{\code{dp}:}{Object of class \code{\linkS4class{DisplayPars}},
      inherited from class \code{\linkS4class{GdObject}}}
    \item{\code{name}:}{Object of class \code{"character"}, 
      inherited from class \code{\linkS4class{GdObject}}}
    \item{\code{imageMap}:}{Object of class \code{\linkS4class{ImageMap}}, 
      inherited from class \code{\linkS4class{GdObject}}}
  }
}

\section{Extends}{
Class \code{"\linkS4class{AnnotationTrack}"}, directly.
Class \code{"\linkS4class{StackedTrack}"}, by class "AnnotationTrack", distance 2.
Class \code{"\linkS4class{RangeTrack}"}, by class "AnnotationTrack", distance 3.
Class \code{"\linkS4class{GdObject}"}, by class "AnnotationTrack", distance 4.
}

\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "ATrack")}:
		Initialize the Object}
	 }
}

\author{
Renan Sauteraud
}

\seealso{
  \code{\linkS4class{AnnotationTrack}}
  
  \code{\linkS4class{DisplayPars}}

  \code{\linkS4class{GdObject}}

  \code{\link{plotTracks}}
}

\examples{
# Object construction
aTrack<-ATrack(start=c(20,60),end=c(40,100),name="random.anno",id=c("small","big"))
# Stacking example
a2Track<-ATrack(start=c(20,30),end=c(40,100),name="stacking=dense",id=c("small","big"),stacking="dense",fill=c("black","orange"))
a3Track<-ATrack(start=c(20,30),end=c(40,100),name="no stacking",id=c("small","big"),fill=c("black","orange"))
# Plotting
plotTracks(trackList=c(aTrack,a2Track,a3Track),showFeatureId=TRUE)
}

\keyword{classes}
