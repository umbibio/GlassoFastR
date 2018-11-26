\name{GlassoFast2}
\alias{GlassoFast2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parralel implementation of graphical lasso
%%  ~~function to do ... ~~
}
\description{
This function is a faster alternative to the "glasso" function in the glasso package (Friedman et al. 2008). This package is based on the algorithm (Fort subroutine) of Sustik and Calderhead (2012) and Jones, David T., et al. "PSICOV" (2011)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GlassoFast2(S, rho, thr = 1e-04, thread=1, maxIt = 10000, start = c("cold", "warm"),
w.init = NULL, wi.init = NULL, shrinkage = F,approximation=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
%%     ~~Describe \code{S} here~~
Covariance matrix (a p by p symmetric matrix)
}
  \item{rho}{
The regularization parameter for lasso. (a non-negative value or a p by p matrix of regularization parameters)
%%     ~~Describe \code{rho} here~~
}
  \item{thread}{
Optional parameter for  maximum number of threads in order to run the program in parrallel utilizing open mp. Default 1.
%%     ~~Describe \code{trace} here~~
}
  \item{thr}{
Optional threshold for convergence. Default is 1.e-4.
%%     ~~Describe \code{thr} here~~
}
  \item{maxIt}{
Optional maximum number of iterations of outer loop. Default is 10,000.
%%     ~~Describe \code{maxIt} here~~
}
  \item{start}{
Optional type of start. Cold start is default. Using warm start, can provide starting values for w and wi.
%%     ~~Describe \code{start} here~~
}
  \item{w.init}{
Optional starting values for estimated covariance matrix (p by p). Only needed when start="warm" is specified

%%     ~~Describe \code{w.init} here~~
}
  \item{wi.init}{
Optional starting values for estimated inverse covariance matrix (p by p). Only needed when start="warm" is specified

%%     ~~Describe \code{wi.init} here~~
}
  \item{shrinkage}{
Optional flag for perturbing covariance matrix S and testing by Cholesky decomposition to avoid any negative eignevalues in order to speed up convergence. Default FALSE.
%%     ~~Describe \code{trace} here~~
}
  \item{approximation}{
Optional flag for useing approximate Lasso algorithm. Default FALSE.
%%     ~~Describe \code{trace} here~~
}

}
\details{
Estimate a sparse inverse covariance matrix using a lasso (L1) penalty, following the Friedman et al. (2008) approach.
The function is a wrapper of the faster and corrected (for non-termination convergence issues) which described by Sustik and Calderhead (2012) along with parallel implementation using open mp and shrinkaging sample covariance matrix to speed up the convergence time.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{w }{Estimated covariance matrix}
\item{wi }{Estimated inverse covariance matrix}
\item{errflag }{Memory allocation error flag: 0 means no error; !=0 means memory allocation error}
\item{niter }{Number of iterations of outer loop}
}
\references{
Friedman J., Hastie T., Tibshirani R. 2008. Sparse inverse covariance estimation with the graphical lasso. Biostatistics. 9:432-441.

Sustik M.A., Calderhead B. 2012. GLASSOFAST: An efficient GLASSO implementation. UTCS Technical Report TR-12-29:1-3.

Jones, David T., et al. "PSICOV: precise structural contact prediction using sparse inverse covariance estimation on large multiple sequence alignments." Bioinformatics 28.2 (2011): 184-190.

%% ~put references to the literature/web site here ~
}
\author{
Amir Vajdi
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{glasso}
\code{glassoFast}
}
\examples{

# Make a random covariance matrix
S=matrix(runif(100,0,1),ncol=10,nrow = 10)
S=cov(S)


# Compute the LASSO estimates
GlassoFast2(S, rho=.001)

}
