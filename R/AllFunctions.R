# Record the display parameters for Pviz extra classes
.parMappings <- NULL
.makeParMapping <- function(){
	classes <-  c("ATrack", "DTrack", "ProbeTrack","ProteinSequenceTrack","ProteinAxisTrack")
	defs <-  sapply(classes, function(x) as(getClassDef(x)@prototype@dp, "list"), simplify=FALSE)
	if(is.null(.parMappings))
		assignInNamespace(x=".parMappings", value=defs, ns="Pviz")
}

.dpOrDefault<-function(GdObject, par, default=NULL){
  val<-getPar(GdObject, par)
  if(is.null(val)){
    val<-default
  }
  return(val)
}