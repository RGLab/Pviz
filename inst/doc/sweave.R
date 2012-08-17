#R CMD BATCH '--args rnw="Pviz.Rnw"' sweave.R
args=(commandArgs(TRUE))
if (length(args)>0) for(i in 1:length(args)) eval(parse(text=args[[i]]))
Sweave(rnw)
q(save="no")
