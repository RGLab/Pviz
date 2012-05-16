###
# Read alignment from fasta file to get scales and sequences
###
readAlign<-function(fileName)
{
	## EASY MODE:
	## Assume first line is >HXB2, second is hxb2 seq, third is >newSeqID,  fourth is newSeq seq
	path<-"/home/rsautera/Desktop/peptide_microarray/simpleAlign.fasta"
    alFile<-file(path,open="r")
	lineList<-list()
	lineCnt<-1
	while(length(oneLine <- readLines(alFile, n = 4, warn = FALSE))) #n=4 means read four lines (negative value for whole file)
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
	refObj[[3]]<-newSeq
	return(refObj)
}


## Convert the coordinates of an object into a given scale
setGeneric("coord2ext", def=function(obj, refScale) standardGeneric("coord2ext"))

setMethod("coord2ext", signature=(obj="numeric"), function(obj, refScale)
{
	#works for both vectors and integers
	return(unlist(sapply(obj, function(x){min(which(refScale==x))})))
})
			
setMethod("coord2ext", signature=(obj="RangedData"), function(obj, refScale)
{
	if(start(obj)[[1]]==0) { start(obj)[[1]]=1 } #To avoid Inf values
	
	extStart<-coord2ext(start(obj),refScale)
	extEnd<-coord2ext(end(obj),refScale)
	#assign new start coordinates after end to avoid width<0 issues
	end(obj)<-extEnd	
	start(obj)<-extStart
	return(obj)
})


## Insert gaps in peptides list (& change coord syst into extended?)
#  Input: peptideList
#         refScale
#         probeStartList
insertGaps<-function(probeList, probeStartList, refScale)
{
	#convert the probeStartList into extended system
	probeStartList<-coord2ext(probeStartList, refScale)
	newProbeList<-c() #list of peptide sequences with added gaps 
	for(i in 1:length(probeList)) #for each probe
	{
		#browser()
		j<-probeStartList[i]
		newProbe<-numeric()
		cnt<-1
		begin<-FALSE #get rid of the gaps before the peptide

		while(cnt<=nchar(probeList[i])) #loop till we reach the end of the peptide
		{
			if(j>1){prev<-refScale[j-1]}
			else{prev<-0}
			actual<-refScale[j]
			if(actual==prev) #gap
			{
				if(begin) #if gap is not the first element of the probe
				{
					newProbe<-c(newProbe,"-")
				}
			}
			else #no gap
			{

				newProbe<-c(newProbe,substr(probeList[i],cnt,cnt))
				begin<-TRUE
				cnt<-cnt+1
			}
			j<-j+1

		}
		newProbe<-paste(newProbe, collapse="")
		newProbeList<-c(newProbeList, newProbe)
	}
	return(newProbeList)
}


