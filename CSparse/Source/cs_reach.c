#include "cs.h"
/* xi [top...n-1] = nodes reachable from graph of G*P' via nodes in B(:,k).
 * xi [n...2n-1] used as workspace */
csi cs_reach (cs *G, const cs *B, csi k, csi *xi, const csi *pinv)
{
    csi p, n, top, *Bp, *Bi, *Gp ;
    if (!CS_CSC (G) || !CS_CSC (B) || !xi) return (-1) ;    /* check inputs */
    n = G->n ; Bp = B->p ; Bi = B->i ; Gp = G->p ;
    top = n ;//top是啥？第一个非零节点，初始值设置为n，如果没有任何节点进行搜索则不会标记任何节点
    for (p = Bp [k] ; p < Bp [k+1] ; p++)//B的第k列
    {
        if (!CS_MARKED (Gp, Bi [p]))    /* start a dfs at unmarked node i */
        {//如果没有被标记，则对该节点进行深度优先搜索
            top = cs_dfs (Bi [p], G, top, xi, xi+n, pinv) ;//将top赋值为第一个节点
        }
    }
    for (p = top ; p < n ; p++) CS_MARK (Gp, xi [p]) ;  /* restore G */
    return (top) ;
}
