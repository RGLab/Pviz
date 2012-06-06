\name{readAlign}
\alias{readAlign}

\title{
	Function to get the reference sequence and scale from an alignment file.
}

\description{
	\code{readAlign} reads a fasta alignment file and returns the reference sequence with gaps 
	and the scale to convert the sequence coordinates into the extended system.
}

\usage{
	readAlign(filename=NULL)
}

\arguments{
	\item{filename}{The name of the alignment file used to extract the new information.}
}

\value{
	A list with the reference scale as the first element and the reference sequence in second position.
}

\author{
Renan Sauteraud
}

\keyword{functions}