#include "cs.h"
/* solve Gx=b(:,k), where G is either upper (lo=0) or lower (lo=1) triangular */
csi cs_spsolve (cs *G, const cs *B, csi k, csi *xi, double *x, const csi *pinv,
    csi lo)
{
    csi j, J, p, q, px, top, n, *Gp, *Gi, *Bp, *Bi ;
    double *Gx, *Bx ;
    if (!CS_CSC (G) || !CS_CSC (B) || !xi || !x) return (-1) ;
    Gp = G->p ; Gi = G->i ; Gx = G->x ; n = G->n ;
    Bp = B->p ; Bi = B->i ; Bx = B->x ;
    top = cs_reach (G, B, k, xi, pinv) ;        /* xi[top..n-1]=Reach(B(:,k)) */
    for (p = top ; p < n ; p++) x [xi [p]] = 0 ;    /* clear x */
    for (p = Bp [k] ; p < Bp [k+1] ; p++) x [Bi [p]] = Bx [p] ; /* scatter B, 将向量B存储到x向量的相应位置,相当于x=b */
    for (px = top ; px < n ; px++) //计算非零的x元素，为零的跳过
    {
        j = xi [px] ;                               /* xi中存储的是节点,此处拿到非零节点 */


                                                    //如果pinv不为空则使用pinv中的第j个元素，否则使用j
                                                    //pinv是第六章中所要使用的参数，目前假设pinv为空（也就是LU分解）
        J = pinv ? (pinv [j]) : j ;                 /* j maps to col J of G */

        //如果第j列为空则跳过此循环
        if (J < 0) continue ;                       /* column J is empty */

        //lo为0,下三角矩阵，则取第一个，lo为1，上三角矩阵，取最后一个
        x [j] /= Gx [lo ? (Gp [J]) : (Gp [J+1]-1)] ;/* x(j) /= G(j,j) */

        //拿到循环的上下限，即第J列的非零元序号
        //若为上三角矩阵，去掉第一个元素，因为第一个元素为L_jj
        //若为下三角矩阵，去掉最后一个元素，因为最后一个元素为U_jj
        p = lo ? (Gp [J]+1) : (Gp [J]) ;            /* lo: L(j,j) 1st entry */
        q = lo ? (Gp [J+1]) : (Gp [J+1]-1) ;        /* up: U(j,j) last entry  */
        //对第j列进行循环
        for ( ; p < q ; p++)
        {
            x [Gi [p]] -= Gx [p] * x [j] ;          /* x(i) -= G(i,j) * x(j) */
        }
    }
    return (top) ;                                  /* return top of stack */
}
/*
这个函数返回解x的非零模式，其被存储在xi中，为xi[top:n]。
xi有2n个元素，前n个元素为输出栈，后n个元素为cs_dfs中会用到的栈。
数值被存放在稠密向量x中，它在输入时无需初始化,cs_spsolve会自动将其值清空
为了解Lx=b, pinv必须设置为null, lo必须为非零，pinv在LU分解时有用。
*/