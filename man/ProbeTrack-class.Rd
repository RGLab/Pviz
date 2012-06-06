\name{ProbeTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{ProbeTrack-class}
\alias{ProbeTrack}
\alias{drawGD,ProbeTrack-method}
\alias{initialize,ProbeTrack-method}
\alias{ProbeIntensity<-,ProbeTrack-method}
\alias{ProbeIntensity,ProbeTrack-method}
\alias{ProbeSequence<-,ProbeTrack-method}
\alias{ProbeSequence,ProbeTrack-method}
\alias{ProbeStart<-,ProbeTrack-method}
\alias{ProbeStart,ProbeTrack-method}
\alias{start,ProbeTrack-method}
\alias{end,ProbeTrack-method}

\title{Class \code{"ProbeTrack"} and methods}
\description{
	A track to display microarray sequences intensity.
}
\section{Objects from the Class}{
Objects can be created using the constructor function \code{ProbeTrack}.
}

\usage{
ProbeTrack(sequence, intensity, probeStart, protein, name = "ProbeTrack", ...)
}

\arguments{
  \item{sequence}{\code{"list"} of \code{"character"},
  	list of sequences.}

  \item{intensity}{\code{"list"} of \code{"character"},
  	list of intensities for the sequences in the slot \code{sequences}.}

  \item{probeStart}{\code{"list"} of \code{"character"},
  	list of start for the sequences in the slot \code{sequences}.}
  	
  \item{protein}{Name of the protein.}

  \item{name}{Character scalar of the track's name used in the title
    panel when plotting.}
    
  \item{\dots}{Additional items which will all be interpreted as further
    display parameters.}
}

\value{
  The return value of the constructor function is a new object of class \code{ProbeTrack}.
}

\section{Slots}{
  \describe{
    \item{\code{sequence}:}{Object of class \code{"list"},
    	a list of \code{"character"}. Peptides microarray sequences.
    	The length of the list needs to be identical to the length of the lists in the \code{intensity} and \code{probeStart} slots.}
    \item{\code{intensity}:}{Object of class \code{"list"},
    	a list of \code{"numeric"}. Peptides microarray intensities.
    	The length of the list needs to be identical to the length of the lists in the \code{sequence} and \code{probeStart} slots.}
    \item{\code{probeStart}:}{Object of class \code{"list"},
    	a list of \code{"numeric"}. Peptides microarray start.
    	The length of the list needs to be identical to the length of the lists in the \code{sequence} and \code{intensity} slots.}
    \item{\code{protein}:}{Object of class \code{"character"},
    	a name for the protein.}
    \item{\code{dp}:}{Object of class\code{\linkS4class{DisplayPars}}, 
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
    \item{drawGD}{\code{signature(GdObject = "ProbeTrack")}: plot the
      object to a graphics device. The return value of this method is
      the input object, potentially updated during the plotting
      operation. Internally, there are two modes in which the method can
      be called. Either in 'prepare' mode, in which case no plotting is
      done but the object is preprocessed based on the
      available space, or in 'plotting' mode, in which case the actual
      graphical output is created.
      }
    \item{initialize}{\code{signature(.Object = "ProbeTrack")}:
		Initialize the Object
		}
    \item{ProbeIntensity<-}{\code{signature(obj = "ProbeTrack",value="list")}:
    	Set the value of the \code{intensity} slot of a \code{ProbeTrack} object. 
		\emph{Usage:}
		\code{ProbeIntensity<-(obj, value)}
		}
    \item{ProbeIntensity}{\code{signature(obj = "ProbeTrack")}:
    	Returns the \code{intensity} slot of a \code{ProbeTrack} object. 
		\emph{Usage:}
		\code{ProbeIntensity(obj)}
		}
    \item{ProbeSequence<-}{\code{signature(obj = "ProbeTrack",value="list")}:
    	Set the value of the \code{sequence} slot of a \code{ProbeTrack} object. 
		\emph{Usage:}
		\code{ProbeSequence<-(obj, value)}
		}
    \item{ProbeSequence}{\code{signature(obj = "ProbeTrack")}:
    	Returns the \code{sequence} slot of a \code{ProbeTrack} object. 
		\emph{Usage:}
		\code{ProbeSequence(obj)}
		}
    \item{ProbeStart<-}{\code{signature(obj = "ProbeTrack",value="list")}:
    	Set the value of the \code{probeStart} slot of a \code{ProbeTrack} object. 
		\emph{Usage:}
		\code{ProbeStart<-(obj, value)}
		}
    \item{ProbeStart}{\code{signature(obj = "ProbeTrack")}:
    	Returns the \code{probeStart} slot of a \code{ProbeTrack} object. 
		\emph{Usage:}
		\code{ProbeStart(obj)}
		}
	\item{start, end}{\code{signature(obj = "ProbeTrack")}:
		These functions are aliases for the \code{ProbeStart} method.
		They are implemented to be consistent with Gviz.}
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
showClass("ProbeTrack")
}
\keyword{classes}
