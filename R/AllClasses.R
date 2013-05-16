####
## DataTrack class with defaulted chromosome/genome
####
setClass("DTrack",
    contains="DataTrack",
	representation=representation(),
	prototype=prototype(
			dp=DisplayPars(ranges.highlight=GRanges())
	)
)

setMethod("initialize", "DTrack", function(.Object, ...){
    .Object <- callNextMethod() #Call the initialize function of Gviz::DataTrack
	return(.Object)
})

## DTrack constructor
DTrack <- function(range=NULL, start=NULL, end=NULL, width=NULL, data, chromosome="chrR", strand="*",
		genome="all", name="DTrack", ...)
{
	#Create a DataTrack object
	GvizDataTrack<-DataTrack(range=range, start=start, end=end, width=width,
			name=name, chromosome=chromosome, strand=strand, genome=genome,
			data=data, ...)
	#Create DTrack object using the slots of the DataTrack
	new("DTrack", range=GvizDataTrack@range, 
			data=GvizDataTrack@data, 
			chromosome=GvizDataTrack@chromosome, 
			genome=GvizDataTrack@genome,
			strand=GvizDataTrack@strand, 
			name=GvizDataTrack@name, 
			...)
}


###
# AnnotationTrack
###
setClass("ATrack",
		contains="AnnotationTrack",
		representation=representation(),
		prototype()
)

setMethod("initialize", "ATrack", function(.Object, ...){
			.Object <- callNextMethod() #Call the initialize function of Gviz::AnnotationTrack
			return(.Object)
		})

## ATrack constructor
ATrack <- function(range=NULL, start=NULL, end=NULL, width=NULL, group,
		id, strand="*", chromosome="chrR", genome="all", stacking="squish", 
		name="ATrack", fun, selectFun, ...)
{
	#Create an AnnotationTrack object
	GvizAnnotationTrack<-AnnotationTrack(range=range, start=start, end=end, width=width,
			id=id, name=name, chromosome=chromosome, genome=genome, stacking=stacking, 
			fun=fun, selectFun=selectFun, ...)
	#Create ATrack using the AnnotationTrack
	new("ATrack", range=GvizAnnotationTrack@range,
			chromosome=GvizAnnotationTrack@chromosome,
			genome=GvizAnnotationTrack@genome,
			name=GvizAnnotationTrack@name,
			stacking=GvizAnnotationTrack@stacking,
			...)
}


###
# ProbeTrack
# Track to display intensities of peptide microarray datas
###
#Intensity is not supported by R basic display tool (use Cairo)
setClass("ProbeTrack", contains = "GdObject", 
		representation(sequence="list",
				intensity = "list",
				probeStart = "list"
		),
		prototype(
				sequence=list(),
				intensity = list(),
				probeStart = list(),
				dp = DisplayPars(
						color = "black",
						color=rev(heat.colors(8, alpha = 1)),
						size = 2, #it displays better when slightly bigger than annotation tracks
						cex=1
				)
		)
)

setMethod("initialize", "ProbeTrack", function(.Object, sequence, intensity, probeStart, name, ...)
		{
			.makeParMapping()
			.Object@name<-name
			.Object@sequence<-sequence
			.Object@intensity<-intensity
			.Object@probeStart<-probeStart
			.Object <- callNextMethod() #Call the initialize function of Gviz::GdObject
			return(.Object)
		})

## ProbeTrack constructor
ProbeTrack <- function(sequence, intensity, probeStart, name="ProbeTrack", ...)
{
	if(missing(probeStart)) stop("Need probeStart argument to know where to plot the data on the genome")
	if(missing(sequence)) stop("Need sequence argument to know what to plot")
	if(class(sequence)!="list")
	{
		sequence<-list(sequence)
		intensity<-list(intensity)
		probeStart<-list(probeStart)
	}
	##check the consistancy of number of entires 
	if(!(identical(lapply(intensity,length),lapply(probeStart,length))&&identical(lapply(sequence,length),lapply(probeStart,length))))
		stop("sequence ,intensity and probeStart need have identifcal structure!")
	##check the types
	if(any(unlist(lapply(sequence,class))!="character"))
		stop("sequence has to be character vector!")
	if(any(unlist(lapply(intensity,class))!="numeric")||any(!unlist(lapply(probeStart,class))%in%c("numeric","integer")))
		stop("intensity and probeStart have to be numeric vectors!")
	new("ProbeTrack", sequence=sequence,intensity=intensity,probeStart=probeStart, name=name, ...)
}

###
# ProteinSequenceTrack
# Track for the amino acid sequence of the reference
###
setClass("ProteinSequenceTrack",
		contains="RangeTrack",
		representation(
				sequence="character"),
		prototype(
		          dp = DisplayPars(size=0.25,
                          fontcolor=getBioColor(type="AA_ALPHABET"))
		)
)

setMethod("initialize", "ProteinSequenceTrack", function(.Object, sequence, ...)
{
	.makeParMapping()
	#.Object <- Gviz:::.updatePars(.Object, "ProteinSequenceTrack")
	.Object@sequence<-sequence
	.Object <- callNextMethod()
  return(.Object)
})

## ProteinSequenceTrack constructor
ProteinSequenceTrack<-function(sequence=NULL, anno=NULL, range=NULL, start=NULL, end=NULL, name="Sequence", chromosome="chrR", genome="all", ...)
{
	sequence<-as.character(sequence) #In case the user provide AA/DNAstring
	#chr and genome are needded to use methods inherited from RangeTrack
#	if(!is.null(anno))
#	{	
#		#Get the sequence from HivFeature object
#		sequence<-getHivFeatureSeq(anno)
#	}
	if(is.null(sequence))
	{
		stop("A sequence or HivFeature object should be supplied")
	}
  if(is.null(range) || !is(range, "GRAnges"))
  {
    if(is.null(start)){start<-1}
    if(is.null(end)){end<-nchar(sequence)}
    range<-GRanges(IRanges(start=start, end=end), seqnames=chromosome)
  }
	new("ProteinSequenceTrack", sequence=sequence, name=name, range=range, chromosome=chromosome, genome=genome, ...)
}



###
# ProteinAxisTrack
# Track to display the sequence axis
###
setClass("ProteinAxisTrack",
		contains="GenomeAxisTrack",
		representation(addNC="logical"),
		prototype(addNC=FALSE,
				dp=DisplayPars(
						lwd=2
					)
		)
)

setMethod("initialize", "ProteinAxisTrack", function(.Object, addNC, name, ...)
{
	.makeParMapping()
	.Object@addNC<-addNC
	.Object@name<-name
	.Object<-callNextMethod()
			
})
		

## ProteinAxisTrack constructor
ProteinAxisTrack<-function(range=NULL, name="Axis", addNC=FALSE, id=NULL, ...)
{
	GvizGAT<-new("GenomeAxisTrack", name=name, range=range, id=id, ...)
	new("ProteinAxisTrack",range=GvizGAT@range,
			name=GvizGAT@name,
			id=GvizGAT@id,
			addNC=addNC, ...)
}

