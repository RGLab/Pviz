###
# Read alignment from fasta file to get scales and sequences
###
readAlign<-function(filename=NULL)
{
	## EASY MODE:
	## Assume first line is >HXB2, second is hxb2 seq, third is >newSeqID,  fourth is newSeq seq
	if(is.null(filename)){filename<-"/home/rsautera/Desktop/peptide_microarray/newMuscleMultipleAlignment.fasta"}
    alFile<-file(filename,open="r")
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


###
# convertDB
#  Changes the position, aligned,, trimmed and peptide columns
###
convertDB<-function(db=pep_hxb2,filename=NULL,refScale=NULL)
{
	## Read the file
	if(is.null(filename)){filename<-"/home/rsautera/Desktop/peptide_microarray/newMuscleMultipleAlignment.fasta"}
	alFile<-file(filename,open="r")
	lineList<-readLines(alFile, n=16) #16 lines, i.e ref+7 subtypes
	close(alFile) #
	
	refSeq<-lineList[2]
	len<-nchar(refSeq)
	gapCnt<-0
	
	##Get the refScale if needed
	if(is.null(refScale))
	{
		refScale<-numeric(len)
		for(i in 1:len)
		{
			if(substr(refSeq,i,i)=="-")
			{
				gapCnt<-gapCnt+1
			}
			refScale[i]<-i-gapCnt
		}
	}

	##Create the 7 subtype lists
	lCnt<-3 #no need to loop on the reference
	sTypeIDList<-character()
	sTypeSeq<-numeric()
	while(lCnt < length(lineList))
	{
		sType<-gsub("> ","",lineList[lCnt])
		sTypeIDList<-c(sTypeIDList,sType)
		sTypeSeq<-c(sTypeSeq,lineList[lCnt+1])
		lCnt<-lCnt+2
	}
	sTypeList<-list(sTypeIDList,sTypeSeq)
	print(sTypeList)

	##Convert the positions of the db
	ndb<-coord2ext(db,refScale)

	t1<-system.time({

	##Get the aligned sequence
	aligned<-reference<-character(length(rownames(ndb)))
	for(ID in sTypeIDList)
	{
		TFvec<-ndb[[ID]]
		idx<-which(TFvec==TRUE)
		newAligns<-sapply(idx, function(x){
					substr(sTypeSeq[which(sTypeIDList==ID)],start(ndb)[x],end(ndb)[x])
				})
		aligned[idx]<-newAligns
		
		newRef<-sapply(idx, function(x){
					substr(refSeq,start(ndb)[x],end(ndb)[x])
				})
		reference[idx]<-newRef
	}

	##Get the trimmed sequence
	trimmed<-c()
	for(seqCnt in 1:length(aligned))
	{
		trimmedS<-c()
		for(AACnt in 1:nchar(aligned[seqCnt]))
		{
			if(substr(reference[seqCnt],AACnt,AACnt)!="-")
			{
				trimmedS<-c(trimmedS,substr(aligned[seqCnt],AACnt,AACnt))
			}
		}
		trimmed<-c(trimmed,paste(trimmedS,collapse=""))
	}
						
	})#end t1

	##Set the values for the new db
	ndb$aligned<-aligned
	ndb$trimmed<-trimmed
	print(t1)
	return(ndb)
}



## Convert the coordinates of an object into the extended coordinate system
setGeneric("coord2ext", def=function(obj, refScale) standardGeneric("coord2ext"))

setMethod("coord2ext", signature=(obj="numeric"), function(obj, refScale)
{
			extVec<-sapply(obj, function(x){
								min(
										if(length(which(refScale==x)))
											which(refScale==x)
										else
											NaN
									)})#sapply#function#min
			extVec<-extVec[!is.na(extVec)]
			return(extVec)
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


### Insert gaps in peptides list (& change coord syst into extended?)
##  Input: peptideList
##         refScale
##         probeStartList
#insertGaps<-function(probeList, probeStartList, refScale=NULL)
#{
#	#convert the probeStartList into extended system
#	probeStartList<-coord2ext(probeStartList, refScale)
#	newProbeList<-c() #list of peptide sequences with added gaps 
#	for(i in 1:length(probeList)) #for each probe
#	{
#		#browser()
#		j<-probeStartList[i]
#		newProbe<-numeric()
#		cnt<-1
#		begin<-FALSE #get rid of the gaps before the peptide
#
#		while(cnt<=nchar(probeList[i])) #loop till we reach the end of the peptide
#		{
#			if(j>1){prev<-refScale[j-1]}
#			else{prev<-0}
#			actual<-refScale[j]
#			if(actual==prev) #gap
#			{
#				if(begin) #if gap is not the first element of the probe
#				{
#					newProbe<-c(newProbe,"-")
#				}
#			}
#			else #no gap
#			{
#
#				newProbe<-c(newProbe,substr(probeList[i],cnt,cnt))
#				begin<-TRUE
#				cnt<-cnt+1
#			}
#			j<-j+1
#
#		}
#		newProbe<-paste(newProbe, collapse="")
#		newProbeList<-c(newProbeList, newProbe)
#	}
#	return(newProbeList)
#}

#####
### Record the display parameters for each class once
#####
.makeParMapping <- function()
{
	classes <-  c("ATrack", "DTrack", "ProbeTrack","SequenceTrack","ProteinAxisTrack")
	defs <-  sapply(classes, function(x) as(getClassDef(x)@prototype@dp, "list"), simplify=FALSE)
	if(is.null(.parMappings))
		assignInNamespace(x=".parMappings", value=defs, ns="Pviz")
}
.parMappings <- NULL
