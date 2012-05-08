#Generics & Method

#ProbeTrack Accessors
setGeneric("ProbeStart",def=function(obj,...) standardGeneric("ProbeStart"))
setGeneric("ProbeIntensity",def=function(obj,...) standardGeneric("ProbeIntensity"))
setGeneric("ProbeSequence",def=function(obj,...) standardGeneric("ProbeSequence"))

setGeneric("ProbeStart<-",def=function(obj, value) standardGeneric("ProbeStart<-"))
setGeneric("ProbeIntensity<-",def=function(obj, value) standardGeneric("ProbeIntensity<-"))
setGeneric("ProbeSequence<-",def=function(obj, value) standardGeneric("ProbeSequence<-"))

setMethod("ProbeStart","ProbeTrack",function(obj) obj@probeStart)
setMethod("ProbeIntensity","ProbeTrack",function(obj) obj@intensity)
setMethod("ProbeSequence","ProbeTrack",function(obj) obj@sequence)

setReplaceMethod("ProbeStart", "ProbeTrack", function(obj, value){
			obj@probeStart<-value
			return(obj)
		})
setReplaceMethod("ProbeIntensity", "ProbeTrack", function(obj, value){
			obj@intensity<-value
			return(obj)
		})
setReplaceMethod("ProbeSequence", "ProbeTrack", function(obj, value){
			obj@sequence<-value
			return(obj)
		})

#Sequence Track Accessors
setGeneric("getHivFeatureSeq", def=function(HIVF) standardGeneric("getHivFeatureSeq"))
setMethod("getHivFeatureSeq", def=function(HIVF) as.character(HIVF@HIV_db$hxb2AA[1,]))

setGeneric("getSequenceSeq", def=function(GdObject) standardGeneric("getSequenceSeq"))
setMethod("getSequenceSeq", def=function(GdObject) GdObject@sequence)

####
## drawGD for ProbeTrack
## legend is the scale of intensities
####
setMethod(Gviz:::"drawGD", signature("ProbeTrack"), function(GdObject, minBase, maxBase, vpPosition, prepare=FALSE, subset=TRUE) 
#setMethod("drawGD", signature("ProbeTrack"), function(GdObject, minBase, maxBase, vpPosition, prepare=FALSE, subset=TRUE) 
{
		#Select only what is in range
		if(subset)
		{
			GdObject <- subset(GdObject, from=minBase, to=maxBase)
		}
		
		##FIXME make it possible to supply a list of length > 1
		indexVec<-which(ProbeStart(GdObject)[[1]]+15>minBase & ProbeStart(GdObject)[[1]]<maxBase)
		
		ProbeStart(GdObject)[[1]]<-ProbeStart(GdObject)[[1]][indexVec]
		ProbeIntensity(GdObject)[[1]]<-ProbeIntensity(GdObject)[[1]][indexVec]
		ProbeSequence(GdObject)[[1]]<-ProbeSequence(GdObject)[[1]][indexVec]
	
		probeStart<-ProbeStart(GdObject)
		sequence<-ProbeSequence(GdObject)
		intensity<-ProbeIntensity(GdObject)
		
		defIntensityRange<-c(min(GdObject@intensity[[1]]),max(GdObject@intensity[[1]]))
		range.legend<-Gviz:::.dpOrDefault(GdObject, "ylim", defIntensityRange)
		
		
		totalRows<-length(probeStart)
		
		### PREPARE MODE
		if(prepare)
		{
				##viewport for multipe probe subtracks
				pushViewport(viewport(xscale=c(minBase,maxBase),yscale=c(0,1)))#,clip=TRUE))
				popViewport(1)
				
				###LEGEND
				if(Gviz:::.dpOrDefault(GdObject, "legend", FALSE)){	
					
					color.probe<-getPar(GdObject,"color.probe")
					cex <- Gviz:::.dpOrDefault(GdObject, "cex.legend", 1)					
					fontsize <- Gviz:::.dpOrDefault(GdObject, "fontsize.legend", 12)
					fontface <- Gviz:::.dpOrDefault(GdObject, "fontface.legend", 1)
					lineheight <- Gviz:::.dpOrDefault(GdObject, "lineheight.legend", 1)
					fontfamily <- Gviz:::.dpOrDefault(GdObject, "fontfamily.legend", 1)
					pushViewport(viewport(width=unit(1, "npc")-unit(0.2,"inches"),
									gp=gpar(cex=cex, fontsize=fontsize, fontface=fontface,
											lineheight=lineheight)))
					boxSize <- 0.3
					spacing <- 0
					hspacing <- 0.02
					lengths<-as.numeric(convertWidth(stringWidth(as.character(range.legend)), "inches"))
					heights<-as.numeric(convertWidth(stringHeight(as.character(range.legend)), "inches"))

					colWidth <- max(lengths + boxSize + spacing*2)
					availSpace <- Gviz:::vpLocation()$isize
					colNum=length(color.probe)+3
					rowNum <- ceiling(length(color.probe)/colNum)
					rowHeight <- max(c(heights, 0.1))			
										
					vertSpace <- ((rowHeight * rowNum) + (hspacing * (rowNum-1)) + 0.2)#*1.5
							
					#Add the legend to the list of display parameters
					displayPars(GdObject) <- list(".__verticalSpace"=vertSpace, ".__layoutDims"=c(rowNum, colNum),
							".__boxSize"=boxSize, ".__spacing"=spacing, ".__groupLevels"=color.probe)
					popViewport(1)
				}
				
				return(invisible(GdObject))
		}
		
		### DRAWING MODE

		###LEGEND
		# it will display the scale of intensities with max and min
		intensityLevels <- Gviz:::.dpOrDefault(GdObject, ".__groupLevels")
		if(Gviz:::.dpOrDefault(GdObject, "legend", FALSE)){
			lSpace <- getPar(GdObject, ".__verticalSpace")
			pushViewport(viewport(y=1, height=unit(1, "npc") - unit(lSpace, "inches"), just=c(0.5, 1)))
			on.exit({popViewport(1)
						cex <- Gviz:::.dpOrDefault(GdObject, "cex.legend", 0.8)
						fontsize <- Gviz:::.dpOrDefault(GdObject, "fontsize.legend", 12)
						fontface <- Gviz:::.dpOrDefault(GdObject, "fontface.legend", 1)
						#lineheight <- Gviz:::.dpOrDefault(GdObject, "lineheight.legend", 1)
						fontfamily <- Gviz:::.dpOrDefault(GdObject, "fontfamily.legend", 1)
						fontcolor <- Gviz:::.dpOrDefault(GdObject, "fontcolor.legend", Gviz:::.DEFAULT_SHADED_COL)
						pushViewport(viewport(y=0, height=unit(lSpace, "inches"), just=c(0.5, 0),
										gp=gpar(cex=cex, fontsize=fontsize, fontface=fontface, fontcolor=fontcolor)))#,
												#lineheight=lineheight)))
						#grid.rect(gp=gpar(col="green"))
						pushViewport(viewport(width=unit(1, "npc") - unit(0.1, "inches"), height=unit(1, "npc") - unit(0.1, "inches")))

						boxSize <- getPar(GdObject, ".__boxSize")
						spacing <- getPar(GdObject, ".__spacing")
						dims <- getPar(GdObject, ".__layoutDims")
						
						cex<-prop<-Gviz:::.dpOrDefault(GdObject, "size.legend", 1) #vertical size of the legend and its text
						
						height<-lSpace*prop #height of the subspace allowed
						# Draw range text		
						pushViewport(viewport(width=(1/dims[2])/2,
											height=height,
											y=0.5,
											x=0
							))
						grid.text(label = formatC(range.legend[1], format="f",digits=2),
								gp=gpar(cex=cex,col=fontcolor),
								default.units="native",
								x=0.5,
								y=0.5
						)
						popViewport(1)#txtVP
				
						# Draw rectangles
						for(i in seq_along(intensityLevels))
						{
							row <- (((i)-1) %/% dims[2])+1
							col <- (((i)-1) %% dims[2])+1
							pushViewport(viewport(width=(1/dims[2])/2,
											height=height,
											y=0.5,
											x=((1/dims[2])*(col))/2
							))
							grid.rect(gp=gpar(col ="transparent", fill=intensityLevels[i]),)
							popViewport(1)#rectVP
						}
						# Draw range text
						pushViewport(viewport(width=(1/dims[2])/2, 
											height=height,
											y=0.5,
											x=((1/dims[2])*(length(intensityLevels)+1))/2
							))
						grid.text(label = formatC(range.legend[2], format="f",digits=2),
								x=0.5,
								y=0.5,
								gp=gpar(cex=cex,col=fontcolor),
								default.units="native")
						popViewport(1)#txtVP
						popViewport(2)
					})
			}

		### TRACK
		pushViewport(viewport(xscale=c(minBase,maxBase),yscale=c(0,1)))#,clip=TRUE))
		for(curR in 1:totalRows)
		{
			curProbName<-names(probeStart)[curR]
			probeArray<-data.frame(probeStart = probeStart[[curR]]
					,sequence = sequence[[curR]]
					,intensity = intensity[[curR]])
			probeArray<-subset(probeArray,probeStart>=minBase&probeStart<=maxBase)
			
			
			color<-getPar(GdObject,"fontcolor")
			color.probe<-getPar(GdObject,"color.probe")
			alpha<-getPar(GdObject,"alpha")
			width<-maxBase-minBase
			nSeq<-width/15
			
			#group the 15mers by positions
			probeArray$posInterval <- cut(probeArray$probeStart,seq(minBase,maxBase,20),right=F,include.lowest=T)
			
			#decide the maximum tiers where sequences to be displayed
			nRows<-max(table(probeArray$posInterval))
			##viewport for the current probe subtrack
			pushViewport(dataViewport(xData=c(minBase, maxBase), yscale=c(0, nRows+1), extension = 0,
							layout.pos.col = 1, layout.pos.row = curR))
			
			showSeq<-ifelse(width<=200,T,F)
			
			##include the range value for cutting function
			intVec<-c(probeArray$intensity,range.legend)
			#get color by intensity
			probeArray$intensityInterval <- cut(intVec,breaks=length(color.probe),right=F,include.lowest=T)[1:nrow(probeArray)]
			
			idCex<-getPar(GdObject,"cex")
			seqCex<-idCex/nSeq
			setPar(GdObject,"seqCex",seqCex)
			#calculate x,y coordinates of each character
			by(probeArray,probeArray$posInterval,function(curData){
						for(i in 1:nrow(curData))
						{
							AAsequence<-as.character(curData$sequence[i])
							xstart<-curData[i,"probeStart"]-0.5
							xright<-min(xstart+nchar(AAsequence),maxBase+0.5)
							grid.rect(x=xstart,y=i
									,width = xright-xstart,
									,height = 0.9#idCex/nSeq
									,gp =gpar(fill=color.probe[curData$intensityInterval[i]]
											,col="transparent"
											,alpha=alpha
									)
									,default.units = "native", just = c("left", "center"))
							if (showSeq)
							{
								for (j in 1:nchar(AAsequence))
								{
									xPos<-curData[i,"probeStart"]+j-1
									if(xPos<=maxBase)
									{
										grid.text(label=substr(AAsequence,j,j)
												,x=xPos
												,y=i,
												,gp = gpar(col=color
														,cex =seqCex
												)
												,default.units = "native")	 
									}
								}
							}
						}
						
					})
			popViewport(1)
		}
		popViewport(1)
	
		
})

####
## drawGD for SequenceTrack
####
setMethod(Gviz:::"drawGD", signature("SequenceTrack"), function(GdObject, minBase, maxBase, vpPosition, prepare=FALSE, subset=TRUE) 
{
	seq<-getSequenceSeq(GdObject)
	seq<-substr(seq,minBase,maxBase)
	pushViewport(dataViewport(xData=c(minBase, maxBase), yscale=c(0, 1), extension=0))
	
	### PREPARE MODE
	if(prepare)
	{
		popViewport(1)
		return(invisible(GdObject))
	}
	
	
	### DRAWING MODE
	
	cex<-getPar(GdObject,"cex")
	color<-getPar(GdObject,"fontcolor")
	len <- (maxBase-minBase + 1)
	
	#print each char of the sequence
	for(cnt in 1:len)
	{
		char<-substr(seq,cnt,cnt)
		vpLetter<-viewport(x=(1/(len-1))*(cnt-1),
				width=1/(len-1))
		pushViewport(vpLetter)
		grid.text(label=char,gp=gpar(cex=cex,col=color))
		popViewport(1)
	}

	popViewport(1)
	
})	






#####
#####
### Record the display parameters for each class once
.makeParMapping <- function()
{
    classes <-  c("ATrack", "DTrack", "ProbeTrack","SequenceTrack")
    defs <-  sapply(classes, function(x) as(getClassDef(x)@prototype@dp, "list"), simplify=FALSE)
    if(is.null(.parMappings))
        assignInNamespace(x=".parMappings", value=defs, ns="Rviz")
}
.parMappings <- NULL
