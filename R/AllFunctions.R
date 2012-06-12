###
# Read alignment from fasta file to get scales and sequences
###
readAlign<-function(filename=NULL)
{
	## EASY MODE:
	## Assume first line is >HXB2, second is hxb2 seq, third is >newSeqID,  fourth is newSeq seq
#	if(is.null(filename)){filename<-"/home/rsautera/Desktop/peptide_microarray/newMuscleMultipleAlignment.fasta"}
	if(is.null(filename)){filename<-system.file("alignment.fasta", package="Pviz")}
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

#####
# drawHighlight
#####
##FIXME l'overlay deborde sur la gauche si < minBase
.drawHighlight<-function(GdObject, minBase, maxBase)
{
	#gpar for grid.rect
	HRanges<- getPar(GdObject, "ranges.highlight")

	if(length(HRanges))
	{
		#DisplayPars
		alpha.highlight<- Gviz:::.dpOrDefault(GdObject, "alpha.highlight", 1)
		color.highlight<- Gviz:::.dpOrDefault(GdObject, "color.highlight", "black")
		fill.highlight<- Gviz:::.dpOrDefault(GdObject, "fill.highlight", "grey")
		lwd.highlight<- Gviz:::.dpOrDefault(GdObject, "lwd.highlight", 1)
		#Overlay on top of the legend or not
		legend<-Gviz:::.dpOrDefault(GdObject, "legend", FALSE)
		hLegend<-Gviz:::.dpOrDefault(GdObject, "legend.highlight", FALSE)
		if(!hLegend && legend)
		{
			lSpace <- getPar(GdObject, ".__verticalSpace")
			pushViewport(viewport(y=1, height=unit(1, "npc") - unit(lSpace, "inches"),
								just=c(0.5, 1)))
		}
		else
			pushViewport(viewport(xscale=c(minBase,maxBase),yscale=c(0,1)))	
		for(i in 1:length(HRanges))
		{
			grid.rect(x=1/(maxBase-minBase)*(start(HRanges[i])-minBase-0.5),
					width=1/(maxBase-minBase)*width(HRanges[i]), just="left",
					gp=gpar(col=color.highlight,fill=fill.highlight,alpha=alpha.highlight,
							lwd=lwd.highlight,linejoin="mitre"))
		}
		popViewport(1)
	}
}

#####
# .getGapPos
# Returns a vector of the gaps starts and end the sequence associated with the refScale
#####
.getGapPos<-function(refScale, minBase, maxBase)
{
	gapStartList<-c()
	gapEndList<-c()
	i<-minBase
	while(i<=maxBase)
	{
		#gap
		if(i>1 && refScale[i]==refScale[i-1])
		{
			gapStartList<-c(gapStartList,i)
			while(refScale[i]==refScale[i-1])
			{
				i=i+1
			}
			gapEndList<-c(gapEndList,i-1)
		}
		i=i+1
	}
	gapCoords<-list()
	gapCoords<-IRanges(start=gapStartList, end=gapEndList)
	return(gapCoords)
}

#substractRanges<-function()
#{
#	min<-1
#	nSL<-c()
#	nEL<-c()
#	for(i in 1:length(x))
#	{
#		
#	}
#	
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
