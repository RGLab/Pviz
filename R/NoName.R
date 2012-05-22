###
# Read alignment from fasta file to get scales and sequences
###
readAlign<-function(filename=NULL)
{
	## EASY MODE:
	## Assume first line is >HXB2, second is hxb2 seq, third is >newSeqID,  fourth is newSeq seq
	if(is.null(filename)){filename<-"/home/rsautera/Desktop/peptide_microarray/simpleAlign.fasta"}
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
# Create a database object from a multiple alignment fasta file
# Input : Fasta file with the reference
#           The first sequence should be the reference HXB2
#           The following 7 sequences should be the subtypes A,B,C,D,M,CRF01,CRF02 (in no specific order)
###
makeDB<-function(filename=NULL)
{
	tt1<-system.time(0)
	tt2<-system.time(0)

	mainDF<-data.frame()
	if(is.null(filename)){filename<-"/home/rsautera/Desktop/peptide_microarray/newMuscleMultipleAlignment.fasta"}
	alFile<-file(filename,open="r")
	lineList<-readLines(alFile, n=16) #16 lines, i.e ref+7 subtypes
	close(alFile) #
	
	refSeq<-lineList[2] #We assume that the reference sequence is the first on the alignment
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
	
	lCnt<-3 #no need to loop on the ref
				
	print(length(lineList))
	while(lCnt < length(lineList))
	{
		print(lCnt)
		seqTmp<-lineList[lCnt+1]#seq with gaps
			peptide<-c()
			aligned<-c()
			trimmed<-c()
			reference<-c()
			
			len<-nchar(seqTmp)
			gapCnt<-0
			gapTrimCnt<-0
			seqNoGap<-c()  #seq with no gaps
			seqTrimmed<-c()#seq with no insertions
			tmpScale<-numeric(len)
			for(i in 1:len)
			{
				if(substr(seqTmp,i,i)=="-")
				{
					gapCnt<-gapCnt+1
				}
				else
				{
					seqNoGap<-c(seqNoGap,substr(seqTmp,i,i))
				}
				tmpScale[i]<-i-gapCnt
			}
			#Aligned sequence
			seqNoGap<-paste(seqNoGap,collapse="")
			startNGList<-seq(1,nchar(seqNoGap)-15,3)
			endNGList<-startNGList+14
			#Peptide
			startList<-coord2ext(startNGList,tmpScale)
			endList<-coord2ext(endNGList,tmpScale)
			for(i in 1:length(startList))
			{
				peptide<-c(peptide, substr(seqNoGap,startNGList[i],endNGList[i]))
				aligned<-c(aligned, substr(seqTmp,startList[i],endList[i]))
				reference<-c(reference, substr(refSeq, startList[i],endList[i]))
			}

			#for each seq in aligned
			for(seqCnt in 1:length(aligned))
			{
				trimmedS<-c()
				#for each AA in the seq
				for(AACnt in 1:nchar(aligned[seqCnt]))
				{
					if(substr(reference[seqCnt],AACnt,AACnt)!="-")
					{
						trimmedS<-c(trimmedS,substr(aligned[seqCnt],AACnt,AACnt))
					}
				}
				trimmed<-c(trimmed,paste(trimmedS,collapse=""))
			}
#			browser()
			#create a data.frame object to store the infos for this sequence
			dfTmp<-data.frame(start=startList,end=endList,aligned,trimmed,reference,peptide)
			mainDF<-rbind(mainDF,dfTmp)#append the data.frames
			
		lCnt<-lCnt+2
	}
	#{---------->
	mainDF<-aggregate(mainDF,by=list(mainDF$peptide),FUN=function(x) x[1]) #remove duplicates
	mainDF<-mainDF[with(mainDF, order(start)),] #sort the big data.frame
	ir<-IRanges(start=mainDF[["start"]],end=mainDF[["end"]])
	alignObject<-RangedData(ir,aligned=mainDF[["aligned"]],trimmed=mainDF[["trimmed"]])
	rownames(alignObject)<-mainDF[["peptide"]]
	#----------->} ~10s
	return(alignObject)
	
}#end makeDB


#convertDB? Change pos/aligned/trimmed
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


