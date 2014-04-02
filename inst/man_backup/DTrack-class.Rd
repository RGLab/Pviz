\name{DTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{DTrack-class}
\alias{DTrack}
\alias{initialize,DTrack-method}

\title{Class \code{"DTrack"} and methods}
\description{
	A track to display numeric data values along a sequence.
}
\section{Objects from the Class}{
	Objects can be created using the constructor function \code{DTrack}.
}

\usage{
DTrack(range = NULL, start = NULL, end = NULL, width = NULL, data, chromosome = "chrR", 
strand = "*", genome = "all", name = "DTrack", ...)
}

\arguments{
  
  \item{range}{An object of class \code{\linkS4class{IRanges}} or
    \code{\linkS4class{GRanges}} containing the data coordinates, or an
    object of class \code{data.frame} with the two mandatory columns
    'start' and 'end' to be coerced into an \code{IRanges} object. The
    coordinates (i.e., the length of the \code{IRanges} object) have to
    match the columns of the data matrix (see below).
    
    Alternatively, the \code{range} argument may be missing, in which
    case the relevant information has to be provided as individual
    function arguments (see below).

  }

  \item{start, end, width}{Integer vectors, giving the start and the end
    end coordinates for the individual track items, or their width. Two
    of the three need to be specified, and have to be of equal length or
    of length one, in which case this value will be recycled. Otherwise,
    the usual R recycling rules for vectors do not apply.}
  
  \item{data}{A numeric matrix of data points with number of columns
    equal to the number of coordinates in \code{range}, or a numeric
    vector of appropriate length that will be coerced in such a one-row
    matrix. Each individual row is supposed to contain data of a given
    sample, where the coordinates for each single observation are
    constant across samples. Depending on the plotting type of the data
    (see 'Details' and 'Display Parameters' sections), sample grouping
    or data aggregation may be available. Alternatively, this can be a
    character vector of column names that point into the element
    metadata of the \code{range} object. Naturally, this is only
    supported if \code{range} is of class \code{\linkS4class{GRanges}}.}

  \item{name}{Character scalar of the track's name used in the title
    panel when plotting.}

  \item{\dots}{Additional items which will all be interpreted as further
    display parameters.}
    
  \item{chromosome}{Unused argument of class \code{"character"},
  	inherited from class \code{\linkS4class{RangeTrack}},
  	defaulted to "chrR".}

  \item{strand}{Unused argument of class \code{"character"},
  	inherited from class \code{\linkS4class{DataTrack}},
  	defaulted to "chrR".}

  \item{genome}{Unused argument of class \code{"character"},
  	inherited from class \code{\linkS4class{RangeTrack}}, 
  	defaulted to "all".}
  	
}

\value{
  The return value of the constructor function is a new object of class \code{DTrack}.  
}

\section{Slots}{
  
  \describe{
  
    \item{\code{data}:}{Object of class \code{"matrix"}, containing the
      data values to be plotted. Individual rows of the matrix
      correspond to individual samples, and the number of columns has to
      be identical to the feature number of the \code{GRanges} object in
      the \code{range} slot.}
    \item{\code{range}:}{Object of class \code{\linkS4class{IRanges}},
      inherited from class \code{\linkS4class{RangeTrack}}. The
      coordinates for the data values. The length of the object needs to
      be identical to the number of columns of the data matrix in the
      \code{data} slot.}
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
Class \code{"\linkS4class{DataTrack}"}, directly.
Class \code{"\linkS4class{NumericTrack}"}, by class "DataTrack", distance 2.
Class \code{"\linkS4class{RangeTrack}"}, by class "DataTrack", distance 3.
Class \code{"\linkS4class{GdObject}"}, by class "DataTrack", distance 4.
}

\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "DTrack")}:
    Initialize the object }
	 }
}

\author{
Renan Sauteraud
}


\seealso{
  \code{\linkS4class{DataTrack}}
  
  \code{\linkS4class{DisplayPars}}

  \code{\linkS4class{GdObject}}

  \code{\link{plotTracks}}
}
\examples{
# Object construction
dTrack <- DTrack(start=seq(1,1000, len=100), width=10, data=matrix(runif(400), nrow=4), name="random data")

# Plotting
plotTracks(dTrack)
}
\keyword{classes}
