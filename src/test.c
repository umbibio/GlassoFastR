#include <R.h>
#include <stdio.h>

int** main(int *n,int **m,int ***k){

for(int i=0;i<3;i++){
Rprintf("%d",*(m+i));
//*(n+i)=*(m+i);
}

return m;
}
