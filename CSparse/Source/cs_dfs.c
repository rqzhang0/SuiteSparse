#include "cs.h"


//cs_dfs与书中的dfsr有什么区别？
//dfsr会有递归的问题，他会调用dfsr
//但是cs_dfs没有递归调用的问题，不需要递归调用cs_dfs
/* depth-first-search of the graph of a matrix, starting at node j */

/*对G中的第j个节点开始进行dfs搜索*/
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

                                //pinv是什么？pinv在Lx=b中没有用，后面在LU分解中会用到。 jnew是下一个搜索节点？jnew是要进行搜索的列
        jnew = pinv ? (pinv [j]) : j ; //对于Lx=b,默认为jnew=j
        if (!CS_MARKED (Gp, j)) //判断是否被标记,如果没有标记则进行后面的计算
        {
                                    //首先标记当前节点, Gp[j]实际上标记的是L_{jj}，也就是L的对角元
            CS_MARK (Gp, j) ;       /* mark node j as visited */
                                    //pstack的第一次赋值
            pstack [head] = (jnew < 0) ? 0 : CS_UNFLIP (Gp [jnew]) ; //如果jnew小于零，则该节点被标记过，
                                                                    //否则对该节点的数值进行翻转, 在此处不对Gp[jnew]做任何事，只是拿到jnew的第一条边
        }
                                    //用done来判定是否访问了所有相邻节点，初值设置为1，默认没有相邻节点
        done = 1 ;                  /* node j done if no unvisited neighbors */

        p2 = (jnew < 0) ? 0 : CS_UNFLIP (Gp [jnew+1]) ;//对jnew进行判定，与之前不同此处是jnew+1, 更新循环的上限
                                                       //如果jnew是小于零的，p2等于0,后面的循环不会继续
                                                       //如果jnew大于零，p2做为下一个节点的初值，也就是下一列的初值,以p2作为循环的上限，则可以遍历所有的相邻边

                                                //对所有的相邻节点进行遍历，下限为pstack[head]，上限为p2
                                                //pstack[head]代表节点j的指向外的第一条边,循环会从第一条边开始一直寻找相邻节点，知道找到未被标记的节点为止
                                                //p2代表下一列的初始指标，代表了循环的上限
        for (p = pstack [head] ; p < p2 ; p++)  /* examine all neighbors of j */
        {                           //拿到行指标i
            i = Gi [p] ;            /* consider neighbor node i */
                                    //如果该节点已经被标记过，则跳过此次循环
            if (CS_MARKED (Gp, i)) continue ;   /* skip visited node i */

                                    //暂停第j个节点的遍历, 用pstack数组来存放？
                                    //pstack的第二次赋值
            pstack [head] = p ;     /* pause depth-first search of node j */

                                    //将第i个节点也放入搜索序列中, 将xi中的下一个元素中存入第i个节点
            xi [++head] = i ;       /* start dfs at node i */
                                    //将第j个节点标记为未完结
            done = 0 ;              /* node j is not done */
                                    //break直接退出循环,开始对节点i的搜索
                                    //从其相邻节点中寻找，知道找到第一个未被标记的节点再跳出循环,否则继续循环

                                    //break，while 会跳转到对i节点的搜索
            break ;                 /* break, to start dfs (i) */
        }                        /*当完成该节点的所有子节点的搜索后，done设置为1*/
        if (done)               /* depth-first search at node j is done */
        {
                                //如果一个层级搜索完成了，则返回上一级进行搜索,所以head要减一
            head-- ;            /* remove j from the recursion stack */

                                //xi[--top], top是被用来当作output stack的起始值的
                                //top的初值是n，被搜索的节点从后往前存放,所以实际上我们去的时候会从最外层节点开始，依次往次级节点分布   
            xi [--top] = j ;    /* and place in the output stack */
        }
    }
    return (top) ;
}

/*
cs_dfs函数开始的时候会将j放置在递归栈的第0个位置xi[0]。
每次循环中，cs_dfs会开始或继续对第j个节点的(相邻节点)搜索。
如果j在递归栈xi中且未被标记，这说明这是第一次去访问它，
在此情况下，首先标记节点，然后pstack[head]被设置为从节点j中出来的第一条边。
如果一个没有被标记的相邻节点i被发现, 则将其存储到递归栈xi中，对第j个节点的搜索会被暂停。
当对第j个节点的递归栈的搜索最后完成时，他会被从递归栈中移除，且被放置到输出栈中。
*/

/*xi存储dfs的节点信息
 *pstack用来存储dfs边的信息
**/