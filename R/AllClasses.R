##FIXME Gviz classes produce warnings when R CMD INSTALL 

###
# DataTrack class with defaulted chromosome/genome
# @protein: slot for the protein name (label/lengend + subsetting purposes)
###
setClass("DTrack",
    contains="DataTrack",
	representation=representation(protein="character"),
	prototype=prototype()
)

setMethod("initialize", "DTrack", function(.Object, protein, ...){
    .Object@protein<-protein
    .Object <- callNextMethod() #Call the initialize function of Gviz::DataTrack
	return(.Object)
})

## DTrack constructor
DTrack <- function(range=NULL, start=NULL, end=NULL, width=NULL, data, chromosome="chrR", strand="*",
		genome="all", name="DTrack", protein="prot", ...)
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
			protein=protein, ...)
}


###
# AnnotationTrack
# @protein: slot for the protein name
###
setClass("ATrack",
		contains="AnnotationTrack",
		representation=representation(protein="character"),
		prototype()
)

setMethod("initialize", "ATrack", function(.Object, protein, ...){
			.Object@protein<-protein
			.Object <- callNextMethod() #Call the initialize function of Gviz::AnnotationTrack
			return(.Object)
		})

## ATrack constructor
ATrack <- function(range=NULL, start=NULL, end=NULL, width=NULL, group,
		id, strand="*", chromosome="chrR", genome="all", stacking="squish", 
		name="ATrack", protein="prot", fun, selectFun, ...)
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
			protein=protein, ...)
}


###
# ProbeTrack
# Track to display intensities of peptide microarray datas
# @protein: slot for the protein name
###
#Intensity is not supported by R basic display tool (use Cairo)
setClass("ProbeTrack", contains = "GdObject", 
		representation(sequence="list",
				intensity = "list",
				probeStart = "list",
				protein="character"
		),
		prototype(
				sequence=list(),
				intensity = list(),
				probeStart = list(),
				dp = DisplayPars(
						color = "black",
						color.probe=heat.colors(8, alpha = 1),
						size = 2, #it displays better when slightly bigger than annotation tracks
						cex=5
				)
		)
)

setMethod("initialize", "ProbeTrack", function(.Object, sequence, intensity, probeStart, protein, name, ...)
		{
			.makeParMapping()
			.Object@name<-name
			.Object@protein<-protein
			.Object@sequence<-sequence
			.Object@intensity<-intensity
			.Object@probeStart<-probeStart
			.Object <- callNextMethod() #Call the initialize function of Gviz::GdObject
			return(.Object)
		})

## ProbeTrack constructor
ProbeTrack <- function(sequence, intensity, probeStart, protein, name="ProbeTrack", ...)
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
	new("ProbeTrack", sequence=sequence,intensity=intensity,probeStart=probeStart, protein=protein, name=name, ...)
}

###
# SequenceTrack
# Track for the amino acid sequence of the reference
###
setClass("SequenceTrack",
		contains="GdObject",
		representation(
				sequence="character"),
		prototype(
				dp = DisplayPars(size=0.25 #to keep it closer to a potential Axis
						)
		)
)

setMethod("initialize", "SequenceTrack", function(.Object, sequence, name, ...)
{
	.makeParMapping()
	.Object@sequence<-sequence
	.Object@name<-name
	.Object <- callNextMethod()
})

## SequenceTrack constructor
SequenceTrack<-function(anno=NULL, sequence=NULL, name="Sequence", ...)
{
	if(!is.null(anno))
	{	
		#Get the sequence from HivFeature object
		sequence<-getHivFeatureSeq(anno)
	}
	if(is.null(sequence))
	{
		stop("A sequence or HivFeature object should be supplied")
	}
	new("SequenceTrack", sequence=sequence, name=name, ...)
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

setMethod("initialize", "ProteinAxisTrack", function(.Object, addNC, ...)
{
	.makeParMapping()
	.Object@addNC<-addNC
	.Object<-callNextMethod()
			
})
		
		
## ProteinAxisTrack constructor
ProteinAxisTrack<-function(range=NULL, name="Axis", addNC=FALSE, ...)
{
	new("ProteinAxisTrack", range=range, name=name, addNC=addNC, ...)
}




