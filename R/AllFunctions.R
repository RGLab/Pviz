###
# Read alignment from fasta file to get scales and sequences
###
readAlign<-function(filename=NULL)
{
	## EASY MODE:
	## Assume first line is >HXB2, second is hxb2 seq, third is >newSeqID,  fourth is newSeq seq
	if(is.null(filename)){filename<-system.file("extdata/alignment.fasta", package="Pviz")}
        alFile<-file(filename,open="r")
	lineList<-list()
	while(length(oneLine <- readLines(alFile, n = 2, warn = FALSE))) #n=4 means read four lines (negative value for whole file)
	{
		lineList<-c(lineList,oneLine)
	}
	close(alFile) #
	refSeq<-lineList[[2]]
	
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


###
# convertPep
#  Changes the position, aligned,, trimmed and peptide columns
###
convertPep<-function(rd=HIV.db:::pep_hxb2,filename=NULL,refScale=NULL)
{
	## Read the file
	if(is.null(filename)){filename<-system.file("extdata/newMuscleMultipleAlignment.fasta", package="Pviz")}
	alFile<-file(filename,open="r")
	lineList<-readLines(alFile, n=-1L)#16) #16 lines, i.e ref+7 subtypes
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
	while(lCnt < 16) #length(lineList))
	{
		sType<-gsub("> ","",lineList[lCnt])
		sType<-gsub(">","",sType)
		sTypeIDList<-c(sTypeIDList,sType)
		sTypeSeq<-c(sTypeSeq,lineList[lCnt+1])
		lCnt<-lCnt+2
	}

	##Convert the positions of the rd
	nrd<-coord2ext(rd,refScale)

	#t1<-system.time({

	##Get the aligned sequence
	aligned<-reference<-character(length(rownames(nrd)))
	for(ID in sTypeIDList)
	{
		TFvec<-nrd[[ID]]
		idx<-which(TFvec==TRUE)
		newAligns<-sapply(idx, function(x){
					substr(sTypeSeq[which(sTypeIDList==ID)],start(nrd)[x],end(nrd)[x])
				})
		aligned[idx]<-newAligns
		
		newRef<-sapply(idx, function(x){
					substr(refSeq,start(nrd)[x],end(nrd)[x])
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
						
	#})#end t1

	##Set the values for the new rd
	nrd$aligned<-aligned
	nrd$trimmed<-trimmed
	#print(t1)
	return(nrd)
}

#####
# drawHighlight
#####
.drawHighlight<-function(GdObject, minBase, maxBase)
{
	#gpar for grid.rect
	HRanges<- getPar(GdObject, "ranges.highlight")
	if(length(HRanges))
		HRanges<-restrict(HRanges, start=minBase, end=maxBase)
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

#####
### SIVrd
###   return a RangedData like pep_hxb2
#####
.SIVrd<-function(fastaFile)
{
  zTable<-.zTable()
  
  #Parse fasta file
  IN<-file(fastaFile,open="r")
  lineList<-list()
  while(length(oneLine <- readLines(IN, n = 4, warn = FALSE))) #n=4 means read four lines (negative value for whole file)
  {
    lineList<-c(lineList,oneLine)
  }
  close(IN)

  #Lists to create RangedData
  peptideList<-c()
  startList<-c()
  alignedList<-c()
  mac239List<-E660List<-c()
  z1sumList<-z2sumList<-z3sumList<-z4sumList<-z5sumList<-c()
  
  for(seqCnt in c(2,4)) #mac239 E660
  {
    for(pepCnt in seq(1,nchar(lineList[seqCnt])-12,3))
    {
      peptide<-substr(lineList[seqCnt],pepCnt,pepCnt+14) #peptides of 15 AA
      aligned<-substr(lineList[2],pepCnt,pepCnt+14) #aligned on mac239

	  if(seqCnt!=2 & peptide==aligned) #if the peptide is the same in E660 & mac239 do not duplicate the peptide
	  {
          E660List[which(peptideList==peptide)]<-TRUE
		  break
	  }
	  z1sum<-z2sum<-z3sum<-z4sum<-z5sum<-0
	  for(AACnt in 1:15)
	  {
		  AA<-substr(peptide,AACnt,AACnt)
		  if(AA=="")
			  break
		  z1sum<-z1sum+zTable[AA,1]
		  z2sum<-z2sum+zTable[AA,2]
		  z3sum<-z3sum+zTable[AA,3]
		  z4sum<-z4sum+zTable[AA,4]
		  z5sum<-z5sum+zTable[AA,5]
	  }
	  if(seqCnt==2)
	  {
	    mac239List<-c(mac239List,TRUE)
	    E660List<-c(E660List,FALSE)
	  }
	  else if(seqCnt==4)
	  {
		E660List<-c(E660List,TRUE)
	    mac239List<-c(mac239List,FALSE)
	  }
	  peptideList<-c(peptideList,peptide)
	  startList<-c(startList,pepCnt)
      alignedList<-c(alignedList,aligned)
	  z1sumList<-c(z1sumList,z1sum)
	  z2sumList<-c(z2sumList,z2sum)
	  z3sumList<-c(z3sumList,z3sum)
	  z4sumList<-c(z4sumList,z4sum)
	  z5sumList<-c(z5sumList,z5sum)
    }
  }
#  browser()
  nrd<-RangedData(ranges=IRanges(start=startList,width=15), aligned=alignedList,
		  mac239=mac239List, E660=E660List,
		  z1sum=z1sumList,z2sum=z2sumList,z3sum=z3sumList,z4sum=z4sumList,z5sum=z5sumList)
  rownames(nrd)<-peptideList
  nrd<-nrd[order(start(nrd)),]
  return(nrd)

}



####
## SIV_2
####
.SIV_2<-function(fastaFile,orderFile)
{ 
  #Parse fasta file
  IN<-file(fastaFile,open="r")
  lineList<-list()
  while(length(oneLine <- readLines(IN, n = 4, warn = FALSE))) #n=4 means read four lines (negative value for whole file)
  {
	  lineList<-c(lineList,oneLine)
  }
  close(IN)
  #Store sequences
  seq1<-lineList[2] #239
  seq2<-lineList[4] #660
  
  #Parse order file
  IN<-file(orderFile,open="r")
  lineList<-list()
  while(length(oneLine <- readLines(IN, n = -1, warn = FALSE)))
  {
	  lineList<-c(lineList,oneLine)
  }  
  close(IN)
  lineList<-lineList[which(lineList!="" & lineList!=" ")]
  
  #Lists to create RangedData
  peptideList<-c()
  startList<-c()
  alignedList<-c()
  mac239List<-E660List<-c()
  z1sumList<-z2sumList<-z3sumList<-z4sumList<-z5sumList<-c()
  
  idx<-1
  while(idx<=length(lineList))
  {
    if(substr(strsplit(lineList[[idx]],"\\.")[[1]][2],1,1)=="1") #if the peptide is specific to its sequence
	{
	  #Get the sequences
	  #Couples are given with E660 before mac239
	  pep1<-strsplit(lineList[[idx]]," ")[[1]][2]
	  pep2<-strsplit(lineList[[idx+1]]," ")[[1]][2]
	  
	  peptideList<-c(peptideList,pep2,pep1)
	  alignedList<-c(alignedList,pep2,pep2) #mac239 is the ref
	  #Get the positions
	  start<-regexpr(pep2, seq1)[[1]]
	  startList<-c(startList,start,start) #the pos are given relative to mac239 coordinates
      #ID
	  mac239List<-c(mac239List,TRUE,FALSE)
	  E660List<-c(E660List,FALSE,TRUE)
	  #zScores
	  z1sumList<-c(z1sumList,.zSum(pep2,"z1"),.zSum(pep1,"z1"))  
	  z2sumList<-c(z2sumList,.zSum(pep2,"z2"),.zSum(pep1,"z2"))
	  z3sumList<-c(z3sumList,.zSum(pep2,"z3"),.zSum(pep1,"z3"))
	  z4sumList<-c(z4sumList,.zSum(pep2,"z4"),.zSum(pep1,"z4"))
	  z5sumList<-c(z5sumList,.zSum(pep2,"z5"),.zSum(pep1,"z5"))
	  
	  idx<-idx+2
    }
	else
	{
		#Sequence
		pep<-strsplit(lineList[[idx]]," ")[[1]][2]
		align<-pep
		peptideList<-c(peptideList,pep)
		alignedList<-c(alignedList,pep)
		#Position
		start<-regexpr(pep, seq1)[[1]]
		startList<-c(startList,start)
		#ID
		mac239List<-c(mac239List,TRUE)
		E660List<-c(E660List,TRUE)
		#zSum
		z1sumList<-c(z1sumList,.zSum(pep,"z1"))  
		z2sumList<-c(z2sumList,.zSum(pep,"z2"))
		z3sumList<-c(z3sumList,.zSum(pep,"z3"))
		z4sumList<-c(z4sumList,.zSum(pep,"z4"))
		z5sumList<-c(z5sumList,.zSum(pep,"z5"))
		
		idx<-idx+1	
	}

	
	
  }
  nrd<-RangedData(ranges=IRanges(start=startList,width=15), aligned=alignedList,
		  mac239=mac239List, E660=E660List,
		  z1sum=z1sumList,z2sum=z2sumList,z3sum=z3sumList,z4sum=z4sumList,z5sum=z5sumList)
  rownames(nrd)<-peptideList
#  nrd<-nrd[order(start(nrd)),]
  return(nrd)
}

#
# Returns the zZsum for the AAString
#
.zSum<-function(AAString, Z)
{
  rownames<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  colnames<-c("z1","z2","z3","z4","z5")
  zTable<-c(
0.24,	-2.32,	0.6,	-0.14,	1.3,
3.52,	2.5,	-3.5,	1.99,	-1.7,
3.05,	1.62,	1.04,	-1.15,	1.61,
3.98,	0.93,	1.93,	-2.46,	0.75,
0.84,	-1.67,	3.71,	0.18,	-2.65,
1.75,	0.5,	-1.44,	-1.34,	0.66,
3.11,	0.26,	-0.11,	-3.04,	-0.25,
2.05,	-4.06,	0.36,	-0.82,	-0.38,
2.47,	1.95,	0.26,	3.9,	0.09,
-3.89,	-1.73,	-1.71,	-0.84,	0.26,
-4.28,	-1.3,	-1.49,	-0.72,	0.84,
2.29,	0.89,	-2.49,	1.49,	0.31,
-2.85,	-0.22,	0.47,	1.94,	-0.98,
-4.22,	1.94,	1.06,	0.54,	-0.62,
-1.66,	0.27,	1.84,	0.7,	2,
2.39,	-1.07,	1.15,	-1.39,	0.67,
0.75,	-2.18,	-1.12,	-1.46,	-0.4,
-4.36,	3.94,	0.59,	3.44,	-1.59,
-2.54,	2.44,	0.43,	0.04,	-1.47,
-2.59,	-2.64,	-1.54,	-0.85,	-0.02)
  dim(zTable)<-c(20,5)
  dimnames(zTable)<-list(rownames,colnames)

  total<-0
  for(AA in 1:nchar(AAString))
  {
    total<-total+zTable[substr(AAString,AA,AA),Z]
  }
  
  return(total)

}
