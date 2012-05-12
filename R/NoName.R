###
# Read alignment from fasta file to get scales and sequences
###

ReadAlign<-function(fileName)
{
	## EASY MODE:
	## Assume first line is >HXB2, second is hxb2 seq, third is >newSeqID,  fourth is newSeq seq
	path<-"/home/rsautera/Desktop/peptide_microarray/simpleAlign.fasta"
    alFile<-file(path,open="r")
	lineList<-list()
	lineCnt<-1
	while(length(oneLine <- readLines(alFile, n = 4, warn = FALSE))) #n=4 means read four lines
	{
		lineList<-c(lineList,oneLine)
	}
	close(alFile) #
	refSeq<-lineList[[2]]
	newSeq<-lineList[[4]]
	
	len<-nchar(refSeq)
	gapCnt<-0
	refScale<-numeric(len)
	for(i in 1:len)
	{
		if(substr(refSeq,i,i)=="-")
		{
			gapCnt<-gapCnt+1
		}
		refScale[i]<-i-gapCnt
	}
	refObj<-list()
	refObj[[1]]<-refScale #the scale
	refObj[[2]]<-refSeq   #the sequence with gaps
	return(refObj)
}


## Convert the coordinates of an object into a given scale
setGeneric("coord2ext", def=function(obj, refScale) standardGeneric("coord2ext"))

setMethod("coord2ext", signature=(obj="numeric"), function(obj, refScale)
{
	#works for both vectors and integers
	return(unlist(sapply(obj, function(x){min(which(refScale==x))})))
})
			
setMethod("coord2ext", signature=(obj="HivFeature"), function(obj, refScale)
{
	if(start(obj)[[1]]==0) { start(obj)[[1]]=1 } #To avoid Inf values
	
	extStart<-unlist(sapply(start(obj), function(x){min(which(refScale==x))}))
	extEnd<-unlist(sapply(end(obj), function(x){min(which(refScale==x))}))
	#assign new start coordinates after end to avoid width<0 issues
	end(obj)<-extEnd	
	start(obj)<-extStart
	return(obj)
})
