\name{convertDB}
\alias{convertDB}

\title{
	Function to convert a database based on an alignment file.
}

\description{
	\code{convertDB} modify the position, aligned and trimmed column of a database, based on an alignment file.
}

\usage{
	convertDB(db=HIV.db:::pep_hxb2,filename=NULL,refScale=NULL)
}

\arguments{
	\item{db}{An object of class \code{RangedData}, on the same format as pep_hxb2 from HIV.db package.
		It needs at least the columns ranges, aligned and trimmed.}
	\item{filename}{The name of the alignment file used to extract the new information.
		The first sequence of the fasta file must be the reference.}
	\item{refScale}{A \code{"numeric"} vector of the same length as the alignment.
		The refScale matches the position of the reference sequence with the alignment coordinates.}
}

\value{
	A new database of class \code{\linkS4class{RangedData}} with updated positions, aligned and trimmed columns with respect of the alignment.
}

\author{
Renan Sauteraud
}

\keyword{functions}