#include "cs.h"
/* C = alpha*A + beta*B */
cs *cs_add (const cs *A, const cs *B, double alpha, double beta)
{
    csi p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values ;
    double *x, *Bx, *Cx ;
    cs *C ;
    //判断矩阵存储方式是否正确以及是否为空
    if (!CS_CSC (A) || !CS_CSC (B)) return (NULL) ;         /* check inputs */
    //判断矩阵存储维度是否相同,行列
    if (A->m != B->m || A->n != B->n) return (NULL) ;
    //m是行数，n是列数，anz是A中非零元个数
    //Bp是B的列指标，Bx是B的矩阵元数组，bnz是B中的非零元个数
    m = A->m ; anz = A->p [A->n] ;
    n = B->n ; Bp = B->p ; Bx = B->x ; bnz = Bp [n] ;
    //w的长度为AB矩阵的行数
    w = cs_calloc (m, sizeof (csi)) ;                       /* get workspace */
    //A，B矩阵元素是否为空
    values = (A->x != NULL) && (Bx != NULL) ;
    //若A B矩阵元素不同时为空，则分配m长度，若同时为空不分配空间
    x = values ? cs_malloc (m, sizeof (double)) : NULL ;    /* get workspace */
    //给C矩阵分配存储空间
    C = cs_spalloc (m, n, anz + bnz, values, 0) ;           /* allocate result*/
    //如果存在空矩阵，则异常退出
    if (!C || !w || (values && !x)) return (cs_done (C, w, x, 0)) ;
    //Cp是C矩阵的列指标向量，Ci是C的行指标向量，Cx是C的矩阵元数组
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    //对列进行循环
    for (j = 0 ; j < n ; j++)
    {
        //用nz来表示矩阵中现有的非零元个数
        Cp [j] = nz ;                   /* column j of C starts here */
        //计算C_*j+=alpha*A_*j
        nz = cs_scatter (A, j, alpha, w, x, j+1, C, nz) ;   /* alpha*A(:,j)*/
        //计算C_*j+=beta*B_*j
        nz = cs_scatter (B, j, beta, w, x, j+1, C, nz) ;    /* beta*B(:,j) */
        //将x中的数据拷贝到C中的元素矩阵
        //对于每列，无需对x进行初始化，cs_scatter会自动进行初始化
        if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    cs_sprealloc (C, 0) ;               /* remove extra space from C */
    return (cs_done (C, w, x, 1)) ;     /* success; free workspace, return C */
}
