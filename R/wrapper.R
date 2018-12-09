library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)

sourceCpp("include/matMulRcpp.cpp")


Distance1= function(m0, m1i){
  mult=eigenMatMult(A=as.matrix(m1i),B=as.matrix(m0))
  trace.m1inv_m0=sum(diag(mult))
  trace.m1inv_m0
}

Distance2= function(m0,  m1i){
  trace.m1inv_m0=sum(diag(m1i)*diag(m0))
  trace.m1inv_m0
}

#call best graph function

bestGraph= function(S,confidence=0.99,threads=1,cutoff=1e-4){

  confidenceCutoff=ifelse(1-confidence<0,0.1,1-confidence)

  rho.Max = max(abs(S))
  rho.Min=ifelse(quantile(abs(S),0.1)<0.01,quantile(abs(S),0.1),0.01)
  rhoS=S
  rhoS[!is.na(rhoS)]=rho.Max
  rhoList = seq(rho.Max,rho.Min , length=101)
  DistanceList = rep(0, 100)
  nonZeroList=rep(0,100)

  glasso.result = GlassoFastR(S, rho=rhoS,threads)
  m0 = glasso.result$wi
  Flag=F
   for (i in 1:100) {

    rhoS[!is.na(rhoS)]=rhoList[i+1]
    glasso.result = GlassoFastR(S, rho=rhoS,threads)
    adjacency = abs(glasso.result$wi) >cutoff
    adjacency=matrix(sapply(adjacency,as.numeric),nrow = nrow(S))
    diag(adjacency)=0
    nonZeroList[i]=length(which(adjacency>0))/(nrow(S)^2)
    m1i = glasso.result$w
    DistanceList[i] = abs(Distance1(m0,m1i)-Distance2(m0,m1i))

    if( abs(DistanceList[i])*nrow(S)*confidence>confidenceCutoff){
        if(Flag){
            if( nonZeroList[i]/nonZeroList[i-1]>1){
              i=i-1
              break()
            }
        }else{
          Flag=T
        }
    }

    m0 =  glasso.result$wi
  }
  rhoPath = data.frame(rho = rhoList[1:i], DistanceList = DistanceList[1:i],EdgeRatio=nonZeroList[1:i])
  plot(rhoPath$rho, rhoPath$DistanceList,xlab = "rho",ylab="trace difference",main="rho path")
  abline(h = 0,col="blue")
  result=list(rhoPath=rhoPath,bestRho=rhoList[i],confidence=confidence)
  return (result)
}


#call GlassoFastR
GlassoFastR=function(S,
                        rho,
                        approx=-1,
                        shrink=-1,
                        threads=1,
                        cutoff=0.0001){



if(is.numeric(rho)){
  tmp=rho
  rho=S
  rho[!is.na(rho)]=tmp
}
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
