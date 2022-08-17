#include "cs.h"
//C_*j=\sum_k A_*k B_kj
//x代表的是C_*j
//beta是B_kj
/* x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse */
csi cs_scatter (const cs *A, csi j, double beta, csi *w, double *x, csi mark,
    cs *C, csi nz)
{
    csi i, p, *Ap, *Ai, *Ci ;
    double *Ax ;
    //检查矩阵A的格式
    if (!CS_CSC (A) || !w || !CS_CSC (C)) return (-1) ;     /* check inputs */
    //Ap 表示A矩阵的列
    //Ai表示A矩阵的行
    //Ax表示矩阵元
    //Ci表示C矩阵的行指标
    Ap = A->p ; Ai = A->i ; Ax = A->x ; Ci = C->i ;
    //对A矩阵的各列进行循环
    for (p = Ap [j] ; p < Ap [j+1] ; p++)
    {
        //得到行指标,A(i,j)非零
        i = Ai [p] ;                            /* A(i,j) is nonzero */
        //w中存储的是什么东西？ 
        //w[i]被用来判断j列中是否有i指标(行指标)
        //mark=j+1 
        //w is a dense matrix
        //在哪里判定的B_*j不等于零？应该是在外面,cs_scatter会在B的每个非零元处被调用
        //第一次循环时，给w[i]赋值，并告诉C(i,j)不等于零，后面的循环则无需对Ci赋值因为前面已经进行了标记
        //w[i]的作用是防止给Ci进行重复赋值
        if (w [i] < mark)
        {
            w [i] = mark ;                      /* i is new entry in column j */
            Ci [nz++] = i ;                     /* add i to pattern of C(:,j) */
            if (x) x [i] = beta * Ax [p] ;      /* x(i) = beta*A(i,j) */
        }
        else if (x) x [i] += beta * Ax [p] ;    /* i exists in C(:,j) already */
    }
    return (nz) ;
}
