#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define FALSE 0
#define TRUE 1

#define SQR(x) ((x)*(x))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

#define MAXSEQLEN 50000

/* Dump a rude message to standard error and exit */
void fail(char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt) ;
  fprintf(stderr, "*** ");
  vfprintf(stderr, fmt, ap);
  fputc('\n', stderr);

  exit(-1);
}



/* Allocate matrix */
void           *allocmat(int rows, int columns, int size)
{
  int             i;
  void          **p, *rp;

  rp = malloc(rows * sizeof(void *) + sizeof(int));

  if (rp == NULL)
    fail("allocmat: malloc [] failed!");

  *((int *)rp) = rows;

  p = rp + sizeof(int);

  for (i = 0; i < rows; i++)
    if ((p[i] = calloc(columns, size)) == NULL)
      fail("allocmat: malloc [][] failed!");

    return p;
}

/* Allocate vector */
void           *allocvec(int columns, int size)
{
  void          *p;

  p = calloc(columns, size);

  if (p == NULL)
    fail("allocvec: calloc failed!");

  return p;
}


#define EPS (1.1e-15)
#define BIG (1e9)

int glassofast(const int n, double **S, double **L, const double thr, const int maxit, int approxflg, int warm, double **X, double **W)
{
  int i, j, ii, iter, jj, tid;
  double a, b, c, delta, dlx, dw, shr, sum, thrlasso, tmp, wd[MAXSEQLEN], wxj[MAXSEQLEN];

  for (shr=ii=0; ii<n; ii++)
    for (jj=0; jj<n; jj++)
      shr += fabs(S[ii][jj]);

  for (i=0; i<n; i++)
    shr -= fabs(S[i][i]);

  if (shr == 0.0)
  {
    /* S is diagonal. */
    for (ii=0; ii<n; ii++)
      for (jj=0; jj<n; jj++)
        W[ii][jj] = X[ii][jj] = 0.0;

    for (i=0; i<n; i++)
      W[i][i] = W[i][i] + L[i][i];

    for (ii=0; ii<n; ii++)
      for (jj=0; jj<n; jj++)
        X[ii][jj] = 0.0;

    for (i=0; i<n; i++)
      X[i][i] = 1.0 / MAX(W[i][i], EPS);

    return 0;
  }

  shr *= thr/(n-1);
  thrlasso = shr/n;
  if (thrlasso < 2*EPS)
    thrlasso = 2*EPS;

  if (!warm)
  {
    for (ii=0; ii<n; ii++)
      for (jj=0; jj<n; jj++)
      {
        W[ii][jj] = S[ii][jj];
        X[ii][jj] = 0.0;
      }
  }
  else
  {
    for (i=0; i<n; i++)
    {
      for (ii=0; ii<n; ii++)
        X[i][ii] = -X[i][ii]/X[i][i];
      X[i][i] = 0.0;
    }
  }

  for (i=0; i<n; i++)
  {
    wd[i] = S[i][i] + L[i][i];
    W[i][i] = wd[i];
  }

  for (iter = 1; iter<=maxit; iter++)
  {
    dw = 0.0;

#pragma omp parallel for private(i,j,ii,wxj,a,b,c,dlx,delta,sum)
    for (j=0; j<n; j++)
    {
   //   tid = omp_get_thread_num();
    //  Rprintf("%d\t",tid);

      for (ii=0; ii<n; ii++)
        wxj[ii] = 0.0;

      for (i=0; i<n; i++)
        if (X[j][i] != 0.0)
          for (ii=0; ii<n; ii++)
            wxj[ii] += W[i][ii] * X[j][i];

        for (;;)
        {
          dlx = 0.0;

          for (i=0; i<n; i++)
          {
            if (i != j && L[j][i] < BIG)
            {
              a = S[j][i] - wxj[i] + wd[i] * X[j][i];
              b = fabs(a) - L[j][i];
              if (b <= 0.0)
                c = 0.0;
              else if (a >= 0.0)
                c = b / wd[i];
              else
                c = -b / wd[i];

              delta = c - X[j][i];
              if (delta != 0.0 && (!approxflg || fabs(delta) > 1e-6))
              {
                X[j][i] = c;

                for (ii=0; ii<n; ii++)
                  wxj[ii] += W[i][ii] * delta;

                if (fabs(delta) > dlx)
                  dlx = fabs(delta);
              }
            }
          }

          if (dlx < thrlasso)
            break;
        }

        wxj[j] = wd[j];

        for (sum=ii=0; ii<n; ii++)
          sum += fabs(wxj[ii] - W[j][ii]);

#pragma omp critical
        if (sum > dw)
          dw = sum;

        for (ii=0; ii<n; ii++)
          W[j][ii] = wxj[ii];
        for (ii=0; ii<n; ii++)
          W[ii][j] = wxj[ii];
    }

    if (dw <= shr)
      break;
  }

  for (i=0; i<n; i++)
  {
    for (sum=ii=0; ii<n; ii++)
      sum += X[i][ii] * W[i][ii];

    tmp = 1.0 / (wd[i] - sum);

    for (ii=0; ii<n; ii++)
      X[i][ii] = -tmp * X[i][ii];
    X[i][i] = tmp;
  }

  for (i=0; i<n-1; i++)
  {
    for (ii=i+1; ii<n; ii++)
    {
      X[i][ii] = 0.5 * (X[i][ii] + X[ii][i]);
      X[ii][i] = X[i][ii];
    }
  }

  return iter;
}


/* Test Cholesky decomposition on matrix */
int test_cholesky(double **a, const int n)
{
  int i, j, k;
  double sum, diag;

  for (i=0; i<n; i++)
  {
    sum = a[i][i];

    for (k=i-1; k >= 0; k--)
      sum -= a[i][k]*a[i][k];

    if (sum <= 0.0)
      return TRUE;

    diag = sqrt(sum);

#pragma omp parallel for
    for (j=i+1; j<n; j++)
    {
      double sum = a[i][j];

      for (k=0; k<i; k++)
        sum -= a[i][k]*a[j][k];

      a[j][i] = sum / diag;
    }
  }

  return FALSE;
}




int main2(double *cov,double *L,int *size,int *approximation, int *shrink,int *threads, double *cutoff, int *numberofIter, double *wTmp, double *wiTmp)
{
  int             a, b, i, j, k,  s, dimension,  filtflg=0, approxflg=0, initflg=0, maxit=10000, shrinkflg=1,   minseqsep = 5, overrideflg=0, ntries;
  double thresh=1e-4, del,  trialrho, rhodefault = -1.0;
  double sum, score, **pa,  lambda, low_lambda, high_lambda, smean,  r2, mean, sd;
  double  besttd = 1.0, bestrho = 0.001;


  if(*approximation==1)
    approxflg=1;
  if(*cutoff>0)
    thresh=*cutoff;

#ifdef _OPENMP
  omp_set_num_threads(*threads);
#endif

  dimension=*size;

  double **cmat, **rho, **ww, **wwi, **tempmat;


  cmat = allocmat(dimension, dimension, sizeof(double));
  tempmat = allocmat(dimension, dimension, sizeof(double));
  rho = allocmat(dimension, dimension, sizeof(double));
  ww = allocmat(dimension, dimension, sizeof(double));
  wwi = allocmat(dimension, dimension, sizeof(double));

  /* Form the covariance matrix */
  double x;

#pragma omp parallel for default(shared) private(j)
  for (i=0; i<dimension; i++)
  {
    for (j=0; j<dimension; j++)
    {
      cmat[i][j] =cov[i*dimension+j];
      rho[i][j] = L[i*dimension+j];

    }
  }


  /* Shrink sample covariance matrix towards shrinkage target F = Diag(1,1,1,...,1) * smean */

  if (*shrink==1)
  {
    for (smean=i=0; i<dimension; i++)
      smean += cmat[i][i];

    smean /= (double)dimension;

    high_lambda = 1.0;
    low_lambda = 0.0;
    lambda = 0.5;
    int counter=0;
    for (;;)
    {
#pragma omp parallel for default(shared) private(j)
      for (i=0; i<dimension; i++)
        for (j=0; j<dimension; j++)
          if (i != j)
            tempmat[i][j] = cmat[i][j] * (1.0 - lambda);
          else
            tempmat[i][j]=smean * lambda + (1.0 - lambda) * cmat[i][j];
          /* Test if positive definite using Cholesky decomposition */
          if (!test_cholesky(tempmat, dimension))
          {
            if (high_lambda - low_lambda < 0.01)
              break;

            high_lambda = lambda;
            lambda = 0.5 * (lambda + low_lambda);
          }
          else
          {
            low_lambda = lambda;
            lambda = 0.5 * (lambda + high_lambda);
          }
          counter++;
    }
    for (i=0; i<dimension; i++)
      for (j=0; j<dimension; j++)
        if (i != j)
          cmat[i][j] *= (1.0 - lambda);
        else
          cmat[i][j] = smean * lambda + (1.0 - lambda) * cmat[i][j];
  }



  numberofIter[0]=glassofast(dimension, cmat, rho, thresh, maxit, approxflg, initflg, wwi, ww);

#pragma omp parallel for default(shared) private(j)
  for (i=0;i<dimension;i++){
    for(j=0;j<dimension;j++){
      wTmp[i*dimension+j]=ww[i][j];
      wiTmp[i*dimension+j]=wwi[i][j];
    }
  }

         
 FILE *foutwi = fopen("Penlaty.Tab", "w");
 FILE *foutw = fopen("CovarianceEstimated.Tab", "w");

 for (i=0;i<dimension;i++){
    for(j=0;j<dimension;j++){
      fprintf(foutwi,"%lf\t",rho[i][j]);  
       fprintf(foutw,"%lf\t",ww[i][j]);  
   
    }
      fprintf(foutwi,"\n");  
      fprintf(foutw,"\n");  
  } 
    fclose(foutwi);
    fclose(foutw);  


  return 0;
}
