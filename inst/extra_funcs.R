## Sequence alignment and outdated functions

# This function has been implemented in Gviz with 'HighlightTrack-class'
.drawHighlight<-function(GdObject, minBase, maxBase)
{
  minBase<-as.numeric(minBase)
  maxBase<-as.numeric(maxBase)
  #gpar for grid.rect
  HRanges<- getPar(GdObject, "ranges.highlight")
  if(length(HRanges))
    #if not GRanges, assume it is PE
    HRanges<-restrict(HRanges, start=minBase, end=maxBase)
  if(length(HRanges))
  {
    
    #DisplayPars
    alpha.highlight<- .dpOrDefault(GdObject, "alpha.highlight", 1)
    color.highlight<- .dpOrDefault(GdObject, "col.highlight", "transparent")
    fill.highlight<- .dpOrDefault(GdObject, "fill.highlight", "grey")
    lwd.highlight<- .dpOrDefault(GdObject, "lwd.highlight", 1)
    #Overlay on top of the legend or not
    legend<-.dpOrDefault(GdObject, "legend", FALSE)
    hLegend<-.dpOrDefault(GdObject, "legend.highlight", FALSE)
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
  classes <-  c("ATrack", "DTrack", "ProbeTrack","ProteinSequenceTrack","ProteinAxisTrack")
  defs <-  sapply(classes, function(x) as(getClassDef(x)@prototype@dp, "list"), simplify=FALSE)
  if(is.null(.parMappings))
    assignInNamespace(x=".parMappings", value=defs, ns="Pviz")
}
.parMappings <- NULL


###
# Read alignment from fasta file to get scales and sequences
###

# Here is a 4x faster way of getting the refSeq:
# function(x){cumsum(sapply(unlist(strsplit(seq, "")), function(x){!grepl("-", x)}, USE.NAMES=FALSE ))}

readAlign<-function(filename=NULL){
  if(is.null(filename)){filename<-system.file("extdata/alignment.fasta", package="PEP.db")}
  lineList <- readLines(filename)
  refSeq <- lineList[[2]]
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



setGeneric("coord2ext", def=function(obj, refScale){ standardGeneric("coord2ext") })

# Convert the coordinates of an object into the extended coordinate given a scale
# Input: An object and a reference scale
# Output: an object of the same type with coordinates in extended system
setMethod("coord2ext", signature=(obj="numeric"), function(obj, refScale){   
  extVec<-sapply(obj, function(x){
    min(
      if(x<0){
        x
      }else if(x==0){
        which(refScale==1)
      }else if(x>refScale[length(refScale)]){
        which(refScale==refScale[length(refScale)])
      }else if(length(which(refScale==x))){
        which(refScale==x)
      }else{        
        NaN
      })})
  extVec<-extVec[!is.na(extVec)]
  return(extVec)
})          

setMethod("coord2ext", signature=(obj="RangedData"), function(obj, refScale){
  if(start(obj)[[1]]==0) { start(obj)[[1]]=1 } #To avoid Inf values
  
  extStart<-coord2ext(start(obj),refScale)
  extEnd<-coord2ext(end(obj),refScale)
  #assign new start coordinates after end to avoid width<0 issues
  end(obj)<-extEnd
  start(obj)<-extStart
  return(obj)
})

setMethod("coord2ext", signature=(obj="IRanges"), function(obj, refScale){
  if(start(obj)[[1]]==0) { start(obj)[[1]]=1 } #To avoid Inf values
  
  extStart<-coord2ext(start(obj),refScale)
  extEnd<-coord2ext(end(obj),refScale)
  #assign new start coordinates after end to avoid width<0 issues
  end(obj)<-extEnd
  start(obj)<-extStart
  return(obj)
})
