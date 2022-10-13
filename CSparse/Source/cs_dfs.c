#include "cs.h"

//cs_dfs与书中的dfsr有什么区别？
//dfsr会有递归的问题，他会调用dfsr
//但是cs_dfs没有递归调用的问题，不需要递归调用cs_dfs
/* depth-first-search of the graph of a matrix, starting at node j */
csi cs_dfs (csi j, cs *G, csi top, csi *xi, csi *pstack, const csi *pinv)
{
    csi i, p, p2, done, jnew, head = 0, *Gp, *Gi ;
    if (!CS_CSC (G) || !xi || !pstack) return (-1) ;    /* check inputs */
    Gp = G->p ; Gi = G->i ;
                //第一个节点为j
    xi [0] = j ;                /* initialize the recursion stack */
    while (head >= 0)//将回调改成循环
    {
                                //拿到要搜索的node
        j = xi [head] ;         /* get j from the top of the recursion stack */
                                //pinv是什么？ jnew是下一个搜索节点？
        jnew = pinv ? (pinv [j]) : j ;
        if (!CS_MARKED (Gp, j)) //判断是否被标记,如果没有标记则进行后面的计算
        {
                                    //首先标记当前节点
            CS_MARK (Gp, j) ;       /* mark node j as visited */
            pstack [head] = (jnew < 0) ? 0 : CS_UNFLIP (Gp [jnew]) ; //如果jnew小于零，则该节点被标记过，
                                                                    //否则对该节点的数值进行翻转
        }
                                    //用done来判定是否访问了所有相邻节点，初值设置为1，默认没有相邻节点
        done = 1 ;                  /* node j done if no unvisited neighbors */
        p2 = (jnew < 0) ? 0 : CS_UNFLIP (Gp [jnew+1]) ;//对jnew进行判定，与之前不同此处是jnew+1

                                                //对所有的相邻节点进行遍历，下限为pstack[head]，上限为p2
                                                //pstack[head]代表
                                                //p2代表
        for (p = pstack [head] ; p < p2 ; p++)  /* examine all neighbors of j */
        {                           //拿到行指标i
            i = Gi [p] ;            /* consider neighbor node i */
                                    //如果该节点已经被标记过，则跳过此次循环
            if (CS_MARKED (Gp, i)) continue ;   /* skip visited node i */
                                    //暂停第j个节点的遍历, 用pstack数组来存放？
            pstack [head] = p ;     /* pause depth-first search of node j */
                                    //将第i个节点也放入搜索序列中, 将xi中的下一个元素中存入第i个节点
            xi [++head] = i ;       /* start dfs at node i */
                                    //将第j个节点标记为未完结
            done = 0 ;              /* node j is not done */
                                    //break直接退出循环,开始对节点i的搜索
                                    //从其相邻节点中寻找，知道找到第一个未被标记的节点再跳出循环,否则继续循环
            break ;                 /* break, to start dfs (i) */
        }
        if (done)               /* depth-first search at node j is done */
        {
            head-- ;            /* remove j from the recursion stack */
            xi [--top] = j ;    /* and place in the output stack */
        }
    }
    return (top) ;
}
