\name{ProteinAxisTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{ProteinAxisTrack-class}
\alias{ProteinAxisTrack}
\alias{drawGD,ProteinAxisTrack-method}
\alias{initialize,ProteinAxisTrack-method}

\title{Class \code{"ProteinAxisTrack"}}
\description{
	A track to display an axis and coordinates for the sequence.
}
\section{Objects from the Class}{
	Objects can be created using the constructor function \code{ProteinAxisTrack}.
}

\usage{
	ProteinAxisTrack(range=NULL, name="Axis", addNC=FALSE, id=NULL, ...)
}


\arguments{
  	\item{range}{Optional \code{\linkS4class{GRanges}} or
   		\code{\linkS4class{IRanges}} object to highlight certain regions on the axis.}
		
	\item{addNC}{Object of class \code{"logical"}.
		Displays the Amino-terminal and Carboxyl-terminal ends on the SequenceTrack.}
		
	\item{name}{Character scalar of the track's name used in the title
   		panel when plotting.}
   		
   	\item{id}{A character vector of the same length as \code{range}
   		containing identifiers for the ranges. If missing, the constructor
   		will try to extract the ids from \code{names(range)}.}
   		    
    \item{\dots}{Additional items which will all be interpreted as further
        display parameters.}
}

\value{
  The return value of the constructor function is a new object of class \code{ProteinAxisTrack}.  
}


\section{Slots}{
  \describe{
    \item{\code{addNC}:}{Object of class \code{"logical"},
    	Display the N/C ends of the sequence}
    \item{\code{range}:}{Object of class \code{\linkS4class{GRanges}},
    	highlighted on the axis}
    \item{\code{dp}:}{Object of class\code{\linkS4class{DisplayPars}}, 
    	inherited from class \code{\linkS4class{GdObject}}}
    \item{\code{name}:}{Object of class \code{"character"}, 
      	inherited from class \code{\linkS4class{GdObject}}}
    \item{\code{imageMap}:}{Object of class \code{\linkS4class{ImageMap}}, 
     	inherited from class \code{\linkS4class{GdObject}}}
  }
}

\section{Extends}{
Class \code{"\linkS4class{GenomeAxisTrack}"}, directly.
Class \code{"\linkS4class{GdObject}"}, by class "GenomeAxisTrack", distance 2.
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
    \item{initialize}{\code{signature(.Object = "ProteinAxisTrack")}:
    	Initialize the object }
	 }
}

\author{
Renan Sauteraud
}

\seealso{
  \code{\linkS4class{GenomeAxisTrack}}
  
  \code{\linkS4class{DisplayPars}}

  \code{\linkS4class{GdObject}}

  \code{\link{plotTracks}}
}

\examples{
# Object construction
paxTrack <- ProteinAxisTrack()
pax2 <- ProteinAxisTrack(addNC=TRUE)
pax3 <- ProteinAxisTrack(littleTicks=TRUE)

# Plotting
plotTracks(c(paxTrack,pax2,pax3), from=1, to=100)
}

\keyword{classes}
