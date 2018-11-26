
glassoParallel=function(S,RhoS,approx=-1,shrink=-1,threads,cutoffs=0.0001){
dyn.load("src/GlassoFastParallel.so")

ptm <- proc.time()
out=.C("main",cov=as.numeric(S),L=as.numeric(RhoS),size=as.integer(nrow(S)),approximation=as.integer(approx),
     shrink=as.integer(shrink),thread=as.integer(1),cutoff=as.numeric(cutoffs))
(runTimeOurs=proc.time()-ptm)

}
