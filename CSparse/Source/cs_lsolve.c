#include "cs.h"
/* solve Lx=b where x and b are dense.  x=b on input, solution on output. */
csi cs_lsolve (const cs *L, double *x)
{
    csi p, j, n, *Lp, *Li ;
    double *Lx ;
    if (!CS_CSC (L) || !x) return (0) ;                     /* check inputs */
    n = L->n ; Lp = L->p ; Li = L->i ; Lx = L->x ;
    for (j = 0 ; j < n ; j++)
    {
        //x_j=x_j/l_{jj}
        //默认了第一个元素就是对角元
        x [j] /= Lx [Lp [j]] ;
        // for each i > j and l_{ij} !=0
        // x_i=x_i-l_{ij}*x_j
        for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
        {
            x [Li [p]] -= Lx [p] * x [j] ;
        }
    }
    return (1) ;
}
