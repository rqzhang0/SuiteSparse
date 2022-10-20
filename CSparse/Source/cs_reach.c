#include "cs.h"
/* xi [top...n-1] = nodes reachable from graph of G*P' via nodes in B(:,k).
 * xi [n...2n-1] used as workspace */

/* 此处介绍xi向量的作用。
 * xi [top...n-1]被用存储B(:,k)中的可达到的节点，也就是x中的非零元素
 * 如果xi[top...n-1]是被用来存储矩阵边的信息,他决定了节点j搜索的起始边,用于搜索相邻节点*/

csi cs_reach (cs *G, const cs *B, csi k, csi *xi, const csi *pinv)
{
    csi p, n, top, *Bp, *Bi, *Gp ;
    if (!CS_CSC (G) || !CS_CSC (B) || !xi) return (-1) ;    /* check inputs */
    n = G->n ; Bp = B->p ; Bi = B->i ; Gp = G->p ;
    top = n ;//top是啥？使用来记录可用节点的
             //初值设置为n，如果没有任何节点进行搜索则不会标记任何节点
    for (p = Bp [k] ; p < Bp [k+1] ; p++)//B的第k列
    {
        if (!CS_MARKED (Gp, Bi [p]))    /* start a dfs at unmarked node i */
        {//如果没有被标记，则对该节点进行深度优先搜索
            top = cs_dfs (Bi [p], G, top, xi, xi+n, pinv) ;//深度优先搜索之后，
                                    //top的值会被更新以保证能够正确的在向量xi中拿到存储的已被搜索的节点
        }
    }
    for (p = top ; p < n ; p++) CS_MARK (Gp, xi [p]) ;  /* restore G */ /*恢复G的内容，将标记过的对角元去除标记*/
    return (top) ;
}

/*
 * cs_reach, 找到所有Lx=b中b的非零元。
 */