#include "cs.h"
/* y = A*x+y */
CS_INT cs_gaxpy (const cs *A, const CS_ENTRY *x, CS_ENTRY *y)
{
    CS_INT p, j, n, *Ap, *Ai ;
    CS_ENTRY *Ax ;
    //判断输入是否为空，且判断A矩阵是否为compressed column form,如果不是则退出
    if (!CS_CSC (A) || !x || !y) return (0) ;       /* check inputs */
    n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    //对列进行循环，同一列放存连续
    for (j = 0 ; j < n ; j++)
    {   
        //对列中的每行非零元进行访问。对于triplet form,首先获得该列行指标的地址范围并对该范围进行循环
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            //进行乘法，对Ax中为零的矩阵元，进行计算因为相乘后仍为零
            //y=Ax可以理解为利用x中的元素对Ax的列做一个线性组合然后相加
            y [Ai [p]] += Ax [p] * x [j] ;
        }
    }
    return (1) ;
}
