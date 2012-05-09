\name{SequenceTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{SequenceTrack-class}
\alias{SequenceTrack}
\alias{drawGD,SequenceTrack-method}
\alias{initialize,SequenceTrack-method}

\title{Class \code{"SequenceTrack"}}
\description{
	A track to display a reference sequence.
}

\section{Objects from the Class}{
	Objects can be created using the constructor function \code{SequenceTrack}.
}

\usage{
	SequenceTrack(anno = NULL, sequence = NULL, name = "Sequence", ...)
}

\arguments{
	\item{anno}{Object of class \code{\link[HIV.db:HivFeature-class]{"HivFeature"}}.
		The Sequence to be displayed needs to be the first element of this \code{\link[HIV.db:HivFeature-class]{"HivFeature"}}.
		It is possible to directl supply a character sequence instead.}
		
	\item{sequence}{Object of class \code{"character"}.
		The sequence to be displayed.}
		
	\item{name}{Character scalar of the track's name used in the title
   		panel when plotting.}
   		    
    \item{\dots}{Additional items which will all be interpreted as further
        display parameters.}
}

\value{
  The return value of the constructor function is a new object of class \code{SequenceTrack}.  
}

\section{Slots}{
  \describe{
    \item{\code{sequence}:}{Object of class \code{"character"},
    	The sequence itself}
    \item{\code{dp}:}{Object of class \code{\linkS4class{DisplayPars}},
      inherited from class \code{\linkS4class{GdObject}}}
    \item{\code{name}:}{Object of class \code{"character"}, 
      inherited from class \code{\linkS4class{GdObject}}}
    \item{\code{imageMap}:}{Object of class \code{\linkS4class{ImageMap}}, 
      inherited from class \code{\linkS4class{GdObject}}}
  }
}

\section{Extends}{
Class \code{"\linkS4class{GdObject}"}, directly.
}

\section{Methods}{
  \describe{
    \item{drawGD}{\code{signature(GdObject = "SequenceTrack")}: plot the
      object to a graphics device. The return value of this method is
      the input object, potentially updated during the plotting
      operation. Internally, there are two modes in which the method can
      be called. Either in 'prepare' mode, in which case no plotting is
      done but the object is preprocessed based on the
      available space, or in 'plotting' mode, in which case the actual
      graphical output is created.
      }
    \item{initialize}{\code{signature(.Object = "SequenceTrack")}:
    	Initialize the object
	  	}
	 }
}

\author{
Renan Sauteraud
}



\seealso{
  \code{\linkS4class{GdObject}}
  
  \code{\linkS4class{DisplayPars}}

  \code{\link{plotTracks}}
}

\examples{
showClass("SequenceTrack")
}
\keyword{classes}
