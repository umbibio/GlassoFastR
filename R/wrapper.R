
glassoParallel=function(S,
                        RhoS,
                        approx=-1,
                        shrink=-1,
                        threads,
                        cutoffs=0.0001){
dyn.load("src/GlassoFastParallel.so")

ptm <- proc.time()
out=a=.C("main",
         cov=as.numeric(S),
         L=as.numeric(RhoS),
         size=as.integer(nrow(S)),
         approximation=as.integer(approx),
         shrink=as.integer(shrink),
         thread=as.integer(threads),
         cutoff=as.numeric(cutoffs),
         numberofIter=as.integer(0),
         wTmp=as.vector(as.numeric(S)),
         wiTmp=as.vector(as.numeric(S)))
runTimeOurs=proc.time()-ptm

results <- list(w=out[[9]], wi=out[[10]], niter=LASSO[[8]],runTime=runTimeOurs)
return(results)

}
