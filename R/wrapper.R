#library(Rcpp)
#library(RcppArmadillo)
#library(RcppEigen)

#sourceCpp("include/matMulRcpp.cpp")




#call GlassoFastR
GlassoFastR=function(S,
                        rho,
                        approx=-1,
                        shrink=-1,
                        threads=1,
                        cutoff=0.0001){



if(!is.matrix(rho)){
  tmp=rho[1]
  rho=S
  rho[!is.na(rho)]=tmp
}
cat("checking symmetric for rho: ",isSymmetric(rho),"\n dimension of rho is: ",dim(rho),"\n dimension of S is: ",dim(S))  
  rho=t(rho)
S=t(S)
ptm = proc.time()
out=a=.C("main2",
         cov=as.numeric(S),
         L=as.numeric(rho),
         size=as.integer(nrow(S)),
         approximation=as.integer(approx),
         shrink=as.integer(shrink),
         thread=as.integer(threads),
         cutoff=as.numeric(cutoff),
         numberofIter=as.integer(0),
         wTmp=as.vector(as.numeric(S)),
         wiTmp=as.vector(as.numeric(S)))
runTimeOurs=proc.time()-ptm

results = list(w=matrix(out[[9]],nrow=nrow(S),ncol=ncol(S)), wi=matrix(out[[10]],nrow=nrow(S),ncol=ncol(S)), niter=out[[8]],runTime=runTimeOurs)
  return(results)

}
