### R code from vignette source 'pepViz.Rnw'
###################################################
### code chunk number 1: load datas
###################################################
library(pepStat)
library(HIV.db)

###################################################
### code chunk number 2: loadPackage
###################################################
library(Pviz)


##################Peptide set construction
mapFile<-"/home/rsautera/workspace/ex_files/smapping.csv"
pSet1<-makePeptideSet(files=NULL,path="/home/rsautera/workspace/ex_files/RV2/",mapping.file=mapFile,use.flags=FALSE,rm.control.list=c("empty","none","JPT-","Ig","Cy","landmark"),log=TRUE,verbose=TRUE)
#
#data(pSet1.rda)
#
data(pep_hxb2)
ndb<-convertDB(db=pep_hxb2) #Any object created using ndb will be in extended coordinates.
psSet1<-summarizePeptides(pSet1,summary="median",position=ndb)
Zpep = read.csv("Programs/git/pepStat/misc/peptides_zpep.csv") 
feat = peptide(psSet1)
ind = match(feat, as.character(Zpep$peptide))
Zpep = Zpep[ind,]
pnSet1 = NormalizeArray(psSet1, robust=TRUE, standard=FALSE, method="Zpep",centered=TRUE
		, verbose=TRUE)
pnmSet1<-slidingMean(pnSet1, width=9, verbose=TRUE)

load("/home/rsautera/workspace/ex_files/treatment.rda")
#treatment<-treatment[1:198]
treatment<-treatment[1:16]
RV144.V<-pnSet1[,!grepl("P",treatment)]
RV144s.V<-psSet1[,!grepl("P",treatment)]
RV144.V.smooth<-slidingMean(pnSet1[,!grepl("P",treatment)],width=9)
cutoff<-1.1
RV144.V.calls<-makeCalls(RV144.V.smooth,cutoff=cutoff,method="absolute")
RVVcalls<-RV144.V.calls
dim(RVVcalls)<-c(1423,1)
RV144.V.freq<-apply(RVVcalls,1,mean)


#seqByClade<-apply(clade(RV144.V.smooth),2,function(x,y){y[x]},y=peptide(RV144.V.smooth))
seqByClade<-apply(clade(RV144.V.smooth),2,function(x,y){y[x]},y=RV144.V.smooth@featureRange$aligned)
posByClade<-apply(clade(RV144.V.smooth),2,function(x,y){y[x]},y=start(RV144.V.smooth))
intByClade<-apply(clade(RV144.V.smooth),2,function(x,y){y[x]},y=RV144.V.freq)
#
#data(pepMicroarrayEx)
#seqByClade<-pepMicroarrayEx$probeSeq
#posByClade<-pepMicroarrayEx$probePos
#intByClade<-pepMicroarrayEx$probeFreq
#

##################Scale
alignObj <- readAlign(filename = system.file("extdata/alignment.fasta",package = "Pviz"))
refScale<-alignObj[[1]]
refSeq<-alignObj[[2]]

##################Database
HIV_db<-loadFeatures(ref="env", refScale=refScale)  #gp160
envBase<-getFeature(HIV_db)#, name=c("env"))
envStart=getHXB2Coordinates(envBase)[1,][1] 
envEnd=getHXB2Coordinates(envBase)[1,][2] 

##################Getting Features from database for annotation
proteins<-getFeature(HIV_db,category="protein",start=envStart,end=envEnd,frame=getFrame(envBase))
antis<-getEpitope(envBase,name=c("VRC01"))
helix<-getChildren(envBase,category=c("helix"))

##################Track construction
dT6<-t(exprs(pnmSet1)[,8:9,drop=FALSE])
#dPos<-coord2ext(pepStat::position(pnmSet1),refScale)
dPos<-pepStat::position(pnmSet1)
rd6<-DTrack(range=IRanges(start=dPos,end=dPos),groups=rownames(dT6),data=dT6,genome='hxb2',protein="gp120",alpha=0.5,col=c("blue","red"),cex=1,fontsize.title=20)
ra2<-ATrack(id=proteins@values@unlistData@listData$name,start=start(proteins),end=end(proteins),chromosome="chrgp160",genome='hxb2',name="Protein",protein="gp120",fill="red",size=1)
ra3<-ATrack(id=helix@values@unlistData@listData$name,start=start(helix),end=end(helix),genome='hxb2',name="Helix",fill="gray",protein="gp120")
ra6<-ATrack(id=antis@values@unlistData@listData$name,start=start(antis),end=end(antis),genome='hxb2',name="Epitopes",fill="gray",protein="gp120")


rpext<-ProteinAxisTrack( col.range="black",littleTicks=TRUE)
rpref<-ProteinAxisTrack(refScale=refScale)
rpr<-ProteinAxisTrack(range=IRanges(start=c(40,110), end=c(50,220)), littleTicks=TRUE,
		refScale=refScale, col.range="black",col.gap="green") 
rs<-SequenceTrack(refSeq)#,ranges.highlight=IRanges(start=c(5,24), end=c(10,45)),alpha.highlight=0.5,lwd.highlight=0)
library(RColorBrewer)
mypalette<-brewer.pal(9,"YlOrRd")
rp<-ProbeTrack(seqByClade$B, intByClade$B, posByClade$B, protein="gp120", name="sequence(B)",
		color.probe=mypalette
)
g1<-GenomeAxisTrack(littleTicks=TRUE)

ra2<-ATrack(id=proteins@values@unlistData@listData$name,start=start(proteins),end=end(proteins),chromosome="chrgp160",
		genome='hxb2',name="Protein",protein="gp120",fill="red",size=1,
		ranges.highlight=range(ra3),alpha.highlight=0.5,lwd.highlight=0)

rd6<-DTrack(range=IRanges(start=dPos,end=dPos),
		groups=rownames(dT6),data=dT6,genome='hxb2',protein="gp120",
		col=c("black","orange"),cex=1,fontsize.title=20,
		ranges.highlight=IRanges(start=c(5,24), end=c(10,45)), lwd.highlight=0,
		alpha.highlight=0.5,legend.highlight=TRUE)
rd7<-DTrack(range=IRanges(start=dPos,width=0),#end=dPos),
		groups=rownames(dT6),data=dT6,genome='hxb2',protein="gp120",
		col=c("black","orange"),cex=1,fontsize.title=20,
		ranges.highlight=range(ra3), lwd.highlight=0,
		alpha.highlight=0.5,legend=TRUE)

plotTracks(trackList=c(rpext,rpref,rs,rpr,rs,ra2,rd6,ra3),from=120, to=150)#, type=c("p", "smooth"),showFeatureId=TRUE)
plotTracks(c(rpr), from=1, to=500)
#####
save(pepMicroarrayEx,file="/home/rsautera/workspace/BioC/devel/Pviz/data/pepMicroarrayEx.rda")
save(pepExprEx,file="/home/rsautera/workspace/BioC/devel/Pviz/data/pepExprEx.rda")
#####


gp<-.getGapPos(refScale,1,890)
nr<-.substractRanges(gp)
