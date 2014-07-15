#Generics & Method

#show methods
setMethod("show", "DTrack", def = function(object){
  msg <- sprintf(paste("DTrack '%s'\n| positions: %s\n| samples:%s", sep = ""),
                 names(object), length(object),
                 nrow(values(object)))
  cat(paste(msg, collapse = "\n"), "\n")
})

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

setMethod("start", "ProbeTrack", function(x) if(length(x)) as.integer(ProbeStart(x)[[1]]) else NULL)
setMethod("end", "ProbeTrack", function(x) start(x))

setGeneric("getSequenceSeq", def=function(GdObject) standardGeneric("getSequenceSeq"))
setMethod("getSequenceSeq", def=function(GdObject) GdObject@sequence)

#ProteinAxisTrack Accessors
setGeneric("getNC", def=function(obj) standardGeneric("getNC"))
setMethod("getNC", def=function(obj) obj@addNC)


##################### drawGD redefinition
setMethod("drawGD", signature("DTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, subset=TRUE, ...){
	GdObject <- callNextMethod(GdObject,minBase=minBase, maxBase=maxBase, prepare=prepare, subset=subset) #Call the drawGD method of Gviz::DataTrack
	return(invisible(GdObject))
})

setMethod("drawGD", signature("ATrack"), function(GdObject, minBase, maxBase, prepare=FALSE, subset=TRUE, ...){
	GdObject <- callNextMethod(GdObject,minBase=minBase, maxBase=maxBase, prepare=prepare, subset=subset) #Call the drawGD method of Gviz::AnnotationTrack
	return(invisible(GdObject))
})

setMethod("drawGD", signature("ProbeTrack"), function(GdObject, minBase, maxBase, vpPosition, prepare=FALSE, subset=TRUE) {
  #Select only what is in range
  if(subset){
    GdObject <- subset(GdObject, from=minBase, to=maxBase)
  }

  if(minBase==0){minBase<-1}
  indexVec<-which(ProbeStart(GdObject)[[1]]+15>minBase & ProbeStart(GdObject)[[1]]<maxBase)

  ProbeStart(GdObject)[[1]]<-ProbeStart(GdObject)[[1]][indexVec]
  ProbeIntensity(GdObject)[[1]]<-ProbeIntensity(GdObject)[[1]][indexVec]
  ProbeSequence(GdObject)[[1]]<-ProbeSequence(GdObject)[[1]][indexVec]

  probeStart<-ProbeStart(GdObject)
  sequence<-ProbeSequence(GdObject)
  intensity<-ProbeIntensity(GdObject)

  defIntensityRange<-c(min(GdObject@intensity[[1]]),max(GdObject@intensity[[1]],1))
  range.legend<-.dpOrDefault(GdObject, "ylim", defIntensityRange)

  totalRows<-length(probeStart)

  ### PREPARE MODE
  if(prepare){
    ##viewport for multipe probe subtracks
    pushViewport(viewport(xscale=c(minBase,maxBase),yscale=c(0,1)))#,clip=TRUE))
    popViewport(1)

    ###LEGEND
    if(.dpOrDefault(GdObject, "legend", FALSE)){
      color<-getPar(GdObject,"color")
      cex <- .dpOrDefault(GdObject, "cex.legend", 1)
      fontsize <- .dpOrDefault(GdObject, "fontsize.legend", 12)
      fontface <- .dpOrDefault(GdObject, "fontface.legend", 1)
      lineheight <- .dpOrDefault(GdObject, "lineheight.legend", 1)
      fontfamily <- .dpOrDefault(GdObject, "fontfamily.legend", 1)
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
      colNum=length(color)+3
      rowNum <- ceiling(length(color)/colNum)
      rowHeight <- max(c(heights, 0.1))

      vertSpace <- ((rowHeight * rowNum) + (hspacing * (rowNum-1)) + 0.2)#*1.5

      #Add the legend to the list of display parameters
      displayPars(GdObject) <- list(".__verticalSpace"=vertSpace, ".__layoutDims"=c(rowNum, colNum),
      		".__boxSize"=boxSize, ".__spacing"=spacing, ".__groupLevels"=color)
      popViewport(1)
    }

    return(invisible(GdObject))
  }

  ### DRAWING MODE

  ###LEGEND
  # it will display the scale of intensities with max and min
  intensityLevels <- .dpOrDefault(GdObject, ".__groupLevels")
  if(.dpOrDefault(GdObject, "legend", FALSE)){
    lSpace <- getPar(GdObject, ".__verticalSpace")
    pushViewport(viewport(y=1, height=unit(1, "npc") - unit(lSpace, "inches"), just=c(0.5, 1)))
    on.exit({popViewport(1)
      cex <- .dpOrDefault(GdObject, "cex.legend", 0.8)
      fontsize <- .dpOrDefault(GdObject, "fontsize.legend", 12)
      fontface <- .dpOrDefault(GdObject, "fontface.legend", 1)
      #lineheight <- .dpOrDefault(GdObject, "lineheight.legend", 1)
      fontfamily <- .dpOrDefault(GdObject, "fontfamily.legend", 1)
      fontcolor <- .dpOrDefault(GdObject, "fontcolor.legend", Gviz:::.DEFAULT_SHADED_COL)
      pushViewport(viewport(y=0, height=unit(lSpace, "inches"), just=c(0.5, 0),
      				gp=gpar(cex=cex, fontsize=fontsize, fontface=fontface, fontcolor=fontcolor)))#,
      pushViewport(viewport(width=unit(1, "npc") - unit(0.1, "inches"), height=unit(1, "npc") - unit(0.1, "inches")))

      boxSize <- getPar(GdObject, ".__boxSize")
      spacing <- getPar(GdObject, ".__spacing")
      dims <- getPar(GdObject, ".__layoutDims")

      cex<-prop<-.dpOrDefault(GdObject, "size.legend", 1) #vertical size of the legend and its text

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
      					x=((1/dims[2])*(col))/2))
      	grid.rect(gp=gpar(col ="transparent", fill=intensityLevels[i]),)
      	popViewport(1)#rectVP
      }
      # Draw range text
      pushViewport(viewport(width=(1/dims[2])/2,
      					height=height,
      					y=0.5,
      					x=((1/dims[2])*(length(intensityLevels)+1))/2))
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
  pushViewport(viewport(xscale=c(minBase,maxBase),yscale=c(0,1)))
  for(curR in 1:totalRows)
  {
    curProbName<-names(probeStart)[curR]
    probeArray<-data.frame(probeStart = probeStart[[curR]]
    		,sequence = sequence[[curR]]
    		,intensity = intensity[[curR]])
    probeArray<-subset(probeArray,probeStart>=minBase&probeStart<=maxBase)


    fontcolor<-getPar(GdObject,"fontcolor")
    color<-getPar(GdObject,"color")
    alpha<-getPar(GdObject,"alpha")
    width<-maxBase-minBase
    nSeq<-width/15

    #group the 15mers by positions
    probeArray$posInterval <- cut(probeArray$probeStart,seq(minBase,maxBase,20),
                                  right=FALSE, include.lowest=TRUE)

    #decide the maximum tiers where sequences to be displayed
    nRows<-max(table(probeArray$posInterval))
    ##viewport for the current probe subtrack
    pushViewport(dataViewport(xData=c(minBase, maxBase), yscale=c(0, nRows+1), extension = 0,
    				layout.pos.col = 1, layout.pos.row = curR))

    showSeq<-ifelse(width<=200, TRUE, FALSE)

    ##include the range value for cutting function
    intVec<-c(probeArray$intensity,range.legend)
    #get color by intensity
    probeArray$intensityInterval <- cut(intVec, breaks=length(color),
                                        right=FALSE, include.lowest=TRUE)[1:nrow(probeArray)]

    idCex<-getPar(GdObject,"cex")
    seqCex<-idCex*5/nSeq

    #calculate x,y coordinates of each character
    by(probeArray,probeArray$posInterval,function(curData){
      for(i in 1:nrow(curData)){
      	AAsequence<-as.character(curData$sequence[i])
      	xstart<-curData[i,"probeStart"]-0.5
      	xright<-min(xstart+nchar(AAsequence),maxBase+0.5)
      	grid.rect(x=xstart,y=i
      			,width = xright-xstart,
      			,height = 0.9#idCex/nSeq
      			,gp =gpar(fill=color[curData$intensityInterval[i]]
      					,col="transparent"
      					,alpha=alpha
      		        )
      			,default.units = "native", just = c("left", "center"))
      	if (showSeq){
      		for (j in 1:nchar(AAsequence))
      		{
      			xPos<-curData[i,"probeStart"]+j-1
      			if(xPos<=maxBase)
      			{
      				grid.text(label=substr(AAsequence,j,j)
      						,x=xPos
      						,y=i,
      						,gp = gpar(col=fontcolor
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
  return(invisible(GdObject))

})

# ####
## drawGD for ProteinSequenceTrack
# ####
setMethod("drawGD", signature("ProteinSequenceTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, ...) {
  seq<-getSequenceSeq(GdObject)
  lenSeq<-nchar(seq)
  if(maxBase-minBase<1){ maxBase<-lenSeq}
  seq<-substr(seq,minBase,maxBase)
#   if(maxBase>lenSeq){
#     endSpace<-rep(" ", maxBase-lenSeq)
#     endSpace<-paste(endSpace,collapse="")
#     seq<-paste(seq,endSpace,sep="")
#   }
  cex<-getPar(GdObject,"cex")
  fontsize<-getPar(GdObject,"fontsize")
  fontfamily<-getPar(GdObject,"fontfamily")
  fontface<-getPar(GdObject,"fontface")
  lineheight<-getPar(GdObject,"lineheight")
  pushViewport(dataViewport(xData=c(minBase, maxBase), yscale=c(0, 1), extension=0,
                            gp=gpar(cex=cex, fontzise=fontsize, fontfamily=fontfamily,
                                    fontface=fontface, lineheight=lineheight )
  ))

  ### PREPARE MODE
  if(prepare){
    popViewport(1)
    return(invisible(GdObject))
  }
  len<-nchar(seq)
  AA_ALPHABET<-c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
                 "M", "F", "P", "S", "T", "W", "Y", "V", "U", "B", "Z", "X",
                 "*", "-", "+")
  fcol <- .dpOrDefault(GdObject, "fontcolor", getBioColor(type="AA_ALPHABET"))
  delta <- maxBase-minBase
  if(delta==0)
    return(invisible(GdObject))
  lwidth <- max(as.numeric(convertUnit(stringWidth(AA_ALPHABET),"inches")))
  perLetter <- Gviz:::vpLocation()$isize["width"]/(maxBase-minBase+1)
  diff <- Gviz:::.pxResolution(.dpOrDefault(GdObject, "min.width", 2), coord="x")
  if(diff>1 || (maxBase-minBase+1)>=10e6){
    grid.lines(x=unit(c(minBase, nchar(seq)), "native"), y=0.5,
               gp=gpar(col=getPar(GdObject, "col"),
                       lwd=getPar(GdObject, "lwd")))
  }else{
    sequence <- unlist(strsplit(seq, ""))
    at <- seq((minBase+0.5), minBase + length(sequence) - 1 + 0.5, by=1)
    #     sequence[sequence=="-"] <- ""
    if(perLetter<0.5){
      sequence[c(1, length(sequence))] <- ""
    }
    col <- fcol[toupper(sequence)]
    if(lwidth<perLetter && !.dpOrDefault(GdObject, "noLetters", FALSE)){
      grid.text(x=unit(at, "native"), y=0.5, label=sequence,
                gp=gpar(col=col))
    } else {
      grid.rect(x=unit(at, "native"), y=0.05, width=unit(1, "native"), height=0.9,
                gp=gpar(fill=col, col="white"), just=c(0.5, 0))
    }
  }
  popViewport(1)
  return(invisible(GdObject))
})


####
## drawGD for ProteinAxisTrack
####
setMethod("drawGD", signature("ProteinAxisTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, subset=TRUE, ...) {
	if(minBase==0) minBase<-1
	if((maxBase-minBase)<=0)
		return(invisible(GdObject))
	if(subset)
		GdObject <- subset(GdObject, from=minBase, to=maxBase)
	pushViewport(dataViewport(xData=c(minBase, maxBase), yscale=c(0, 1), extension=0))
	cex <- .dpOrDefault(GdObject, "cex", 0.8)
	labelPos <- .dpOrDefault(GdObject, "labelPos", "alternating")
	lwd<-getPar(GdObject,"lwd")
	pres <- Gviz:::.pxResolution()
	textYOff <-  pres["y"]*3
	lwdAdd <- (lwd-1)/2
	dfact <- max(1, .dpOrDefault(GdObject, "distFromAxis", 1))
	littleTicks <- .dpOrDefault(GdObject, "littleTicks", FALSE)
	tickHeight <- (ifelse(littleTicks, 2, 1) * 3 * dfact + lwdAdd) * pres["y"]
	pyOff <- pres["y"]*1.5
	pxOff <- pres["x"]*5

	### PREPARE MODE
	#Figure out the optimal vertical size
	if(prepare){
		nsp <- if(is.null(.dpOrDefault(GdObject, "scale", NULL))){
					(sum(tickHeight, pyOff*2, textYOff*2 + (as.numeric(convertHeight(stringHeight("1"),"native"))/2)*cex)*2*1.3)/pres["y"]
				} else {
					labelPos <- match.arg(labelPos, c("alternating", "revAlternating", "above", "below", "beside"))
					if(labelPos %in% c("above", "below")){
						(sum(tickHeight, pyOff*2 + (as.numeric(convertHeight(stringHeight("1"),"native"))/2)*cex)*2)/pres["y"]
					} else {
						(sum(tickHeight, pyOff*2 + (as.numeric(convertHeight(stringHeight("1"),"native"))/2)*cex))/pres["y"]
					}
				}
		displayPars(GdObject) <- list("neededVerticalSpace"=nsp)
		popViewport(1)
		return(invisible(GdObject))
	}

	### DRAWING MODE
	#Get DisplayParameters
	color <- .dpOrDefault(GdObject, "col", "darkgray")[1]
	alpha<-getPar(GdObject,"alpha")
	fontface <- .dpOrDefault(GdObject, "fontface", 1)

		refScale<-getPar(GdObject, "refScale")
		#For reference coordinates system
		if(!is.null(refScale))
		{
			minB<-refScale[as.numeric(minBase)]
			maxB<-min(refScale[as.numeric(maxBase)],refScale[length(refScale)], na.rm=TRUE)
			axRange<-c(minB,maxB)
			if(length(axRange)==1){axRange<-c(1, axRange)}
			col.gap <- .dpOrDefault(GdObject,"col.gap","gray")
			col.range <- .dpOrDefault(GdObject,"col.range","black")
			len <- (maxBase-minBase)

			tckTmp<-Gviz:::.ticks(axRange)
			tckTmp <- tckTmp[tckTmp<axRange[2]-pxOff*2 & tckTmp>axRange[1]+pxOff*2]
			label<-as.character(tckTmp)
			tck<-numeric(0)

		maxIR<-min(maxBase, length(refScale))
		gapCoords<-.getGapPos(refScale, minBase, maxIR)
		newIRanges<-gaps(gapCoords, minBase, maxIR)
		HR<-range(GdObject)
		#Draw axis without gaps
		for(i in 1:length(newIRanges))
		{
			ifelse(start(newIRanges[i])<=minBase,extS<-0.5,extS<-0)
			ifelse(end(newIRanges[i])>=maxBase,extE<-0.5,extE<-0)
			vpAxisPos<-viewport(x=(1/(len))*(start(newIRanges[i])-minBase-0.5+extS),width=(width(newIRanges[i])-extS-extE)/(len),just=c("left"))
			pushViewport(vpAxisPos)
			grid.segments(x0=0, y0=0.5, x1=1,  y1=0.5,
						default.units="native",
						gp=gpar(col=color, alpha=alpha, lwd=lwd*2, lineend="butt"))
			popViewport(1)

		}
		#Draw highlighted ranges
		if(length(HR))
		{
		for(i in 1:length(HR))
		{
			ifelse(start(HR[i])<=minBase,extS<-minBase-start(HR[i])+0.5,extS<-0)
			ifelse(end(HR[i])>=maxBase,extE<-end(HR[i])-maxBase+0.5,extE<-0)
			vpAxisPos<-viewport(x=(1/(len))*(start(HR[i])-minBase-0.5+extS),width=(width(HR[i])-extS-extE)/(len),just=c("left"))
			pushViewport(vpAxisPos)
			grid.segments(x0=0, y0=0.5, x1=1,  y1=0.5,
					default.units="native",
					gp=gpar(col=col.range, alpha=alpha, lwd=lwd*2, lineend="butt"))
			popViewport(1)
		}
		}
		#Draw gaps
		if(length((gapCoords)))
		{
		for(i in 1:length(gapCoords))
		{

			for(gPos in 1:width(gapCoords[i]))
			{
				ifelse(start(gapCoords[i])<=minBase && gPos==1,x0pos<-0.5,x0pos<-0.2)
				ifelse(end(gapCoords[i])>=maxBase && gPos==width(gapCoords[i]),x1pos<-0.5,x1pos<-0.8)
				vpAxisPos<-viewport(x=(1/(len))*(start(gapCoords[i])-minBase+gPos-1-0.5),width=1/(len),just=c("left"))
				pushViewport(vpAxisPos)
				grid.segments(x0=x0pos, y0=0.5, x1=x1pos,  y1=0.5,
						default.units="native",
						gp=gpar(col=col.gap, alpha=alpha, lwd=lwd))
				popViewport(1)
			}
		}
		}
			tck<-coord2ext(tckTmp,refScale)
		}




		else
		{
			axRange<-c(as.numeric(minBase),as.numeric(maxBase))
			#draw the axis
			grid.segments(x0=minBase, y0=0.5, x1=maxBase,  y1=0.5,
			default.units="native",
			gp=gpar(col=color, alpha=alpha, lwd=lwd*2, lineend="butt"))

			#vertical ticks

			#Get label position (relative to the axis)
			labelPos <- match.arg(labelPos, c("alternating", "revAlternating", "above", "below", "beside"))
			#Get the optimal ticks number and coordinates
			tck <- Gviz:::.ticks(axRange)
			tck <- tck[tck<axRange[2]-pxOff*2 & tck>axRange[1]+pxOff*2]
			tckText <- tck
			label<-as.character(tckText)
		}
		if(length(GdObject) && is.null(refScale))
		{
			rcolor <- .dpOrDefault(GdObject, "col.range", "black")
			diff <- Gviz:::.pxResolution(coord="x")
			#GdObject <- collapseTrack(GdObject, diff=diff, xrange=c(minBase, maxBase))
			start(GdObject) <- pmax(axRange[1], start(GdObject))
			end(GdObject) <- pmin(axRange[2], end(GdObject))
			coords <- cbind(start(GdObject), -0.1, end(GdObject), 0.1)
			y0t<-y1t<-rep(0.5,length(start(GdObject)))
			grid.segments(x0=start(GdObject),x1=end(GdObject),y0=y0t,y1=y1t,
					default.units="native",
					gp=gpar(col=rcolor,alpha=alpha,lwd=lwd*2, lineend="square"))
			vals <- values(GdObject)
#			if(showIds)
#				grid.text(ids, x=start(GdObject) + width(GdObject)/2, y=0,
#						gp=gpar(col=rcol, cex=rcex, fontface=fontface),
#						default.units="native", just=c("center", "center"))
			## Calculate the coordinates for the image map
			map <- as.matrix(Gviz:::.getImageMap(coords))
			rownames(map) <- paste("region", seq_len(nrow(map)), sep="_")
			tags <- lapply(list(title=rownames(map), start=as.character(start(GdObject)), end=as.character(end(GdObject))),
					function(x){ names(x) <- rownames(map); x})
			#imageMap(GdObject) <- Gviz:::ImageMap(coords=map, tags=tags)
		}
		y0t <- rep(c(0.5), length(tck))[1:length(tck)]
		y1t <- y0t + rep(c(tickHeight, -tickHeight), length(tck))[1:length(tck)]
		y0t <- switch(labelPos, "alternating"=y0t, "revAlternating"=-y0t, "above"=abs(y0t), "below"=-abs(y0t), "beside"=y0t)
		y1t <- switch(labelPos, "alternating"=y1t, "revAlternating"=-y1t, "above"=abs(y1t), "below"=-abs(y1t), "beside"=y1t)
		grid.segments(x0=tck, x1=tck, y0=y0t, y1=y1t,  default.units="native", gp=gpar(col=color, alpha=alpha, lwd=lwd, lineend="square"))

		ylabs<-(y1t-0.5)*3+0.5
		ylabs <- y1t + (ifelse(y1t>0.5, 1, -1) * (textYOff + (as.numeric(convertHeight(stringHeight("1"),"native"))/2)*cex))
		grid.text(label=label, x=tck, y=ylabs, just=c("centre", "centre"),
				gp=gpar(cex=cex, fontface=fontface), default.units="native")


		## The scecond level ticks and labels if necessary
		lcex <- cex*0.75 #slightly smaller labels
		if (.dpOrDefault(GdObject, "littleTicks", FALSE) && length(tck)>1)
		{
			if(!is.null(refScale))
			{
				tck<-as.numeric(label)
			}
			avSpace <- min(diff(tck))
			spaceFac <- 1.8
			spaceNeeded <- min(as.numeric(convertWidth(stringWidth(if(is.character(label)) label else "000000000"),"native"))/2)*lcex*spaceFac
			nTcks <- (avSpace %/% spaceNeeded)
			if(nTcks%%2 == 0)
				nTcks <- nTcks-1
			btck <- tck
			if (!(minBase %in% btck))
				btck <- c(minBase, btck)
			if (!(maxBase %in% btck))
				btck <- c(btck, maxBase)
			y0lt <- y1lt <- ltck <- NULL
			for(i in seq_len(length(btck)-1))
			{
				toFill <- btck[i:(i+1)]
				ttck <- if(i==1) rev(toFill[2]-(avSpace/nTcks)*seq_len(nTcks-1)) else toFill[1]+(avSpace/nTcks)*seq_len(nTcks-1)
				ltck <- c(ltck, ttck)
				ord <- if(i==1){ if(y1t[1]>0.5) c(1,-1) else c(-1,1) } else if(y1t[i-1]<0.5) c(1,-1) else c(-1,1)
				y0 <- rep(0.5, length(ttck))[1:length(ttck)]
				y1 <- y0 + rep(ord*tickHeight/2, length(ttck))[1:length(ttck)]
				y0lt <- c(y0lt, switch(labelPos, "alternating"=y0, "revAlternating"=y0, "above"=abs(y0), "below"=-abs(y0)))
				y1lt <- c(y1lt, switch(labelPos, "alternating"=y1, "revAlternating"=y1, "above"=abs(y1), "below"=-abs(y1)))
			}
			endPadding <- pres["x"]*15
			sel <- ltck > min(tck, axRange+endPadding) & ltck < max(tck, axRange-endPadding)
			labelRef<-NULL
			if(!is.null(refScale))
			{
				labelRef<-as.character(as.integer(ltck))
				ltck<-coord2ext(as.integer(ltck),refScale)
			}
			if(length(ltck[sel]) && min(diff(tck))>nTcks)
			{
				grid.segments(x0=ltck[sel], x1=ltck[sel], y0=y0lt[sel], y1=y1lt[sel],  default.units="native",
						gp=gpar(col=color, alpha=alpha, lwd=lwd, lineend="square"))
				ltckText <- ltck[sel]
				if(is.null(labelRef))
					llabel<-as.character(as.integer(ltckText))
				else
					llabel<-labelRef[sel]
				ytlabs <- y1lt + (ifelse(y1lt>0.5, 1, -1) * (textYOff + (as.numeric(convertHeight(stringHeight("1"),"native"))/2)*lcex))
				if(is.character(label))
					grid.text(label=llabel, x=ltck[sel], y=ytlabs[sel], just=c("centre", "centre"),
							gp=gpar(cex=lcex, fontface=fontface), default.units="native")
				else
					for(i in seq_along(llabel))
						grid.text(label=llabel[[i]], x=ltck[sel][i], y=ytlabs[sel][i], just=c("centre", "centre"),
								gp=gpar(cex=lcex, fontface=fontface), default.units="native")
			}
		}


		#Draw NC ends if needed
		if(getNC(GdObject))
		{
			grid.text(label="NH", x=0, y=1, gp=gpar(cex=cex, fontface=fontface))
			grid.text(label="COOH         ", x=1, y=1, gp=gpar(cex=cex, fontface=fontface))
		}

	popViewport(1)
	return(invisible(GdObject))
})
######################################## End drawGD
