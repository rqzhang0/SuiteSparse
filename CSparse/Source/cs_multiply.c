#include "cs.h"
/* C = A*B */
cs *cs_multiply (const cs *A, const cs *B)
{
    //声明变量
    csi p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values, *Bi ;
    double *x, *Bx, *Cx ;
    cs *C ;
    //检查输入矩阵A,B的类型，是否为compressed column 格式
    if (!CS_CSC (A) || !CS_CSC (B)) return (NULL) ;      /* check inputs */
    //检查A,B矩阵的维度，A矩阵的列数要与B矩阵的行数相等
    if (A->n != B->m) return (NULL) ;
    m = A->m ; anz = A->p [A->n] ;
    n = B->n ; Bp = B->p ; Bi = B->i ; Bx = B->x ; bnz = Bp [n] ;
    //w向量，长度为A向量的列数，B向量的行数
    w = cs_calloc (m, sizeof (csi)) ;                    /* get workspace */
    values = (A->x != NULL) && (Bx != NULL) ;
    //若A,B矩阵的矩阵元都不为空，则分配内存，长度与A的列向量数相同
    x = values ? cs_malloc (m, sizeof (double)) : NULL ; /* get workspace */
    //分配C矩阵的内存, anz是A矩阵非零元素个数，bnz是B矩阵的非零元个数,给C矩阵分配内存为A,B矩阵中非零元个数之和
    //为什么要分配为A,B矩阵非零元个数之和？
    //对于矩阵乘法，C矩阵的矩阵元比A,B矩阵的矩阵元多是非常常见的情况(比如A:9*3; B:3*16 => C:9*16)
    C = cs_spalloc (m, n, anz + bnz, values, 0) ;        /* allocate result */
    //如果分配失败或者A,B矩阵元为空以及x为空，则退出
    if (!C || !w || (values && !x)) return (cs_done (C, w, x, 0)) ;
    Cp = C->p ;
    //对B的列进行循环
    for (j = 0 ; j < n ; j++)
    {
        //nz的初值为零，根据语境，nz表示的应该是矩阵C当前非零元个数
        //如果下一列遍历时会出现内存溢出则非配多余内存
        //将realloc放在条件判断中，因为是&&符号，所以在每次判断时都会先判断第一个条件再判断第二个条件
        //这样避免了在每次判断时都非配内存，而是只在需要时分配内存
        if (nz + m > C->nzmax && !cs_sprealloc (C, 2*(C->nzmax)+m))
        {
            return (cs_done (C, w, x, 0)) ;             /* out of memory */
        } 
        //如果重新分配了内存，指针或许会改变所以每次都重新进行赋值
        Ci = C->i ; Cx = C->x ;         /* C->i and C->x may be reallocated */
        //得到第j列的起始指针
        Cp [j] = nz ;                   /* column j of C starts here */
        for (p = Bp [j] ; p < Bp [j+1] ; p++)
        {
            nz = cs_scatter (A, Bi [p], Bx ? Bx [p] : 1, w, x, j+1, C, nz) ;
        }
        if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    cs_sprealloc (C, 0) ;               /* remove extra space from C */
    return (cs_done (C, w, x, 1)) ;     /* success; free workspace, return C */
}
