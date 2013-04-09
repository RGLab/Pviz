### R code from vignette source 'Pviz.Rnw'

###################################################
### code chunk number 1: loading-package
###################################################
library(Pviz)


###################################################
### code chunk number 2: ATrack-example
###################################################
at<-ATrack(start=c(250,480), end=c(320,520), name="Annotations")
plotTracks(at, from=1, to=600)


###################################################
### code chunk number 3: DTrack-example
###################################################
data(PvizVignette)
dt<-DTrack(data=freqEx, start=posEx, width=15, name="Freq")
plotTracks(dt, from=1, to=850, type="l")


###################################################
### code chunk number 4: ProteinAxisTrack-basic
###################################################
pat<-ProteinAxisTrack()
plotTracks(pat, from=1, to=850)	


###################################################
### code chunk number 5: ProteinAxisTrack-options
###################################################
pat<-ProteinAxisTrack(addNC=TRUE, littleTicks=TRUE)
plotTracks(pat, from=1, to=850)


###################################################
### code chunk number 6: ProteinSequenceTrack-basic
###################################################
seqEx
st<-ProteinSequenceTrack(sequence=seqEx, name="env")
plotTracks(trackList=c(pat,st), from=1, to=50)


###################################################
### code chunk number 7: ProteinSequenceTrack-unreadable
###################################################
plotTracks(trackList=c(pat,st), from=1, to=850, cex=0.5)


###################################################
### code chunk number 8: ProbeTrack-basic
###################################################
pt<-ProbeTrack(pepEx, freqEx, posEx) 
plotTracks(pt, from=460, to=530)


###################################################
### code chunk number 9: ProbeTrack-wide-ranges
###################################################
plotTracks(pt)


###################################################
### code chunk number 10: ProbeTrack-legend
###################################################
plotTracks(pt, from=460, to=530, legend=TRUE)	


###################################################
### code chunk number 11: complex-plot
###################################################
pt<-ProbeTrack(pepEx, freqEx, posEx, cex=0.8) 
plotTracks(trackList=c(pat, st, at, pt, dt), from=460, to=530, type="l", legend=TRUE)


