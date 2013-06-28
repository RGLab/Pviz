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
	readAlign(filename=NULL, seqLine=2)
}

\arguments{
	\item{filename}{The name of the alignment fasta file used to extract the new information.
		The first sequence of the fasta file must be the reference.}
	\item{seqLine}{A \code{numeric}. The index of the line on which the target sequence is located.}
}

\value{
	A list with the reference scale as the first element and the reference sequence in second position.
}

\author{
Renan Sauteraud
}

\examples{
#Extract the refScale and gapped sequence from an alignment.
res<-readAlign(system.file("extdata/alignments/Musclehxb2_vaccineInsert.fasta", package="Pviz"))
res[[1]] #The refScale
res[[2]] #The gapped reference sequence
}

\keyword{functions}
