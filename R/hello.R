# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}

R CMD SHLIB psicov.c
dyn.load("src/psicov.so")

dyn.load("src/test.so")
x=as.vector( c(13,14,1))
y=as.vector( c(15,16,1))
d=as.vector( c(17,18,1))
z=rbind(x,y,d)

z2=z+10

dyn.load("src/test.so")
.C("main",n=as.integer(z2),m=as.integer(z),k=as.integer(z2))


S=matrix(runif(16,0,1),ncol=4,nrow = 4)
S=cov(S)
rhoS=matrix(seq(1,16,1),ncol=4)
dyn.load("fastGlasso/psicov.so")
rhoPath=seq(0.001,0.1,by=0.005)
resultDF=data.frame(user=-1 , system=-1, elapsed=-1,m=0.0,stringsAsFactors = F)
for(i in rhoPath){

  rhoS=S
  rhoS[!is.na(rhoS)]=i
ptm <- proc.time()
a=.C("main",cov=as.numeric(S),L=as.numeric(rhoS),size=as.integer(nrow(S)),approximation=as.integer(-1),
     shrink=as.integer(1),thread=as.integer(200),cutoff=as.numeric(0.0001))
(runTimeOurs=proc.time()-ptm)

ptm <- proc.time()
a=.C("main",cov=as.numeric(S),L=as.numeric(rhoS),size=as.integer(nrow(S)),approximation=as.integer(-1),
     shrink=as.integer(-1),thread=as.integer(200),cutoff=as.numeric(0.0001))
(runTimeOursNoShrinkage=proc.time()-ptm)

ptm <- proc.time()
out=glassoFast(S,rho=i)
(runTimeGlassoFast=proc.time()-ptm)

ptm <- proc.time()
out=glasso(S,rho=i)
(runTimeGlasso=proc.time()-ptm)
resultDF=rbind(resultDF,c(runTimeOurs[1:3],m=0.1),c(runTimeOursNoShrinkage[1:3],m=0.4),c(runTimeGlassoFast[1:3],m=0.2),c(runTimeGlasso[1:3],m=0.3))
}

precision=read.table("Precision.Tab",stringsAsFactors = F)

error=round(abs(precision-out$wi),5)
sum(error)

length(which(out$wi!=0))
length(which(precision!=0))
