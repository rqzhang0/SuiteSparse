
#ifndef NDEBUG
#define NDEBUG
#endif

#include <assert.h>
#include "colamd.h"
#include <limits.h>
#include <math.h>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "matrix.h"
#endif 

#if !defined (NPRINT) || !defined (NDEBUG)
#include <stdio.h>
#endif

#ifndef NULL
#define NULL ((void *) 0)
#endif


#ifdef DLONG

#define Int SuiteSparse_long
#define ID  SuiteSparse_long_id
#define Int_MAX SuiteSparse_long_max

#define COLAMD_recommended colamd_l_recommended
#define COLAMD_set_defaults colamd_l_set_defaults
#define COLAMD_MAIN colamd_l
#define SYMAMD_MAIN symamd_l

#else

#define Int int
#define ID "%d"
#define Int_MAX INT_MAX

#define COLAMD_recommended colamd_recommended
#define COLAMD_set_defaults colamd_set_defaults
#define COLAMD_MAIN colamd
#define SYMAMD_MAIN symamd
#define COLAMD_report colamd_report
#define SYMAMD_report symamd_report

#endif

typedef struct Colamd_Col_struct
{
    Int start ;		
			
    Int length ;	
    union
    {
	Int thickness ;	
	Int parent ;	
    } shared1 ;
    union
    {
	Int score ;	
	Int order ;	
    } shared2 ;
    union
    {
	Int headhash ;	
	Int hash ;	
	Int prev ;	
    } shared3 ;
    union
    {
	Int degree_next ;	
	Int hash_next ;		
    } shared4 ;

} Colamd_Col ;

typedef struct Colamd_Row_struct
{
    Int start ;		
    Int length ;	
    union
    {
	Int degree ;	
	Int p ;		
    } shared1 ;
    union
    {
	Int mark ;	
	Int first_column ;
    } shared2 ;

} Colamd_Row ;


#define PUBLIC
#define PRIVATE static

#define DENSE_DEGREE(alpha,n) \
    ((Int) MAX (16.0, (alpha) * sqrt ((double) (n))))

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define ONES_COMPLEMENT(r) (-(r)-1)


#ifndef TRUE
#define TRUE (1)
#endif

#ifndef FALSE
#define FALSE (0)
#endif


#define EMPTY	(-1)

#define ALIVE	(0)
#define DEAD	(-1)

#define DEAD_PRINCIPAL		(-1)
#define DEAD_NON_PRINCIPAL	(-2)

#define ROW_IS_DEAD(r)			ROW_IS_MARKED_DEAD (Row[r].shared2.mark)
#define ROW_IS_MARKED_DEAD(row_mark)	(row_mark < ALIVE)
#define ROW_IS_ALIVE(r)			(Row [r].shared2.mark >= ALIVE)
#define COL_IS_DEAD(c)			(Col [c].start < ALIVE)
#define COL_IS_ALIVE(c)			(Col [c].start >= ALIVE)
#define COL_IS_DEAD_PRINCIPAL(c)	(Col [c].start == DEAD_PRINCIPAL)
#define KILL_ROW(r)			{ Row [r].shared2.mark = DEAD ; }
#define KILL_PRINCIPAL_COL(c)		{ Col [c].start = DEAD_PRINCIPAL ; }
#define KILL_NON_PRINCIPAL_COL(c)	{ Col [c].start = DEAD_NON_PRINCIPAL ; }


#if defined (MATLAB_MEX_FILE) || defined (MATHWORKS)
#define INDEX(i) ((i)+1)
#else
#define INDEX(i) (i)
#endif


PRIVATE Int init_rows_cols
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int p [],
    Int stats [COLAMD_STATS]
) ;

PRIVATE void init_scoring
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int head [],
    double knobs [COLAMD_KNOBS],
    Int *p_n_row2,
    Int *p_n_col2,
    Int *p_max_deg
) ;

PRIVATE Int find_ordering
(
    Int n_row,
    Int n_col,
    Int Alen,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int head [],
    Int n_col2,
    Int max_deg,
    Int pfree,
    Int aggressive
) ;

PRIVATE void order_children
(
    Int n_col,
    Colamd_Col Col [],
    Int p []
) ;

PRIVATE void detect_super_cols
(


    Colamd_Col Col [],
    Int A [],
    Int head [],
    Int row_start,
    Int row_length
) ;

PRIVATE Int garbage_collection
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int *pfree
) ;

PRIVATE Int clear_mark
(
    Int tag_mark,
    Int max_mark,
    Int n_row,
    Colamd_Row Row []
) ;





#define ASSERT(expression) (assert (expression))

PRIVATE void colamd_get_debug	
(
    char *method
) ;

PRIVATE void debug_deg_lists
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int head [],
    Int min_score,
    Int should,
    Int max_deg
) ;

PRIVATE void debug_mark
(
    Int n_row,
    Colamd_Row Row [],
    Int tag_mark,
    Int max_mark
) ;

PRIVATE void debug_matrix
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A []
) ;

PRIVATE void debug_structures
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int n_col2
) ;

#define DEBUG0(params) ;
#define DEBUG1(params) ;
#define DEBUG2(params) ;
#define DEBUG3(params) ;
#define DEBUG4(params) ;

// #define ASSERT(expression)



static size_t t_add (size_t a, size_t b, int *ok)
{
    (*ok) = (*ok) && ((a + b) >= MAX (a,b)) ;
    return ((*ok) ? (a + b) : 0) ;
}

static size_t t_mult (size_t a, size_t k, int *ok)
{
    size_t i, s = 0 ;
    for (i = 0 ; i < k ; i++)
    {
	s = t_add (s, a, ok) ;
    }
    return (s) ;
}

#define COLAMD_C(n_col,ok) \
    ((t_mult (t_add (n_col, 1, ok), sizeof (Colamd_Col), ok) / sizeof (Int)))

#define COLAMD_R(n_row,ok) \
    ((t_mult (t_add (n_row, 1, ok), sizeof (Colamd_Row), ok) / sizeof (Int)))


PUBLIC size_t COLAMD_recommended	
(
    Int nnz,			
    Int n_row,			
    Int n_col			
)
{
    size_t s, c, r ;
    int ok = TRUE ;
    if (nnz < 0 || n_row < 0 || n_col < 0)
    {
	return (0) ;
    }
    s = t_mult (nnz, 2, &ok) ;	    
    c = COLAMD_C (n_col, &ok) ;	    
    r = COLAMD_R (n_row, &ok) ;	    
    s = t_add (s, c, &ok) ;
    s = t_add (s, r, &ok) ;
    s = t_add (s, n_col, &ok) ;	    
    s = t_add (s, nnz/5, &ok) ;	    
    ok = ok && (s < Int_MAX) ;
    return (ok ? s : 0) ;
}



PUBLIC void COLAMD_set_defaults
(
    double knobs [COLAMD_KNOBS]		
)
{
    Int i ;

    if (!knobs)
    {
	return ;			
    }
    for (i = 0 ; i < COLAMD_KNOBS ; i++)
    {
	knobs [i] = 0 ;
    }
    knobs [COLAMD_DENSE_ROW] = 10 ;
    knobs [COLAMD_DENSE_COL] = 10 ;
    knobs [COLAMD_AGGRESSIVE] = TRUE ;	
}



PUBLIC Int SYMAMD_MAIN			
(
    Int n,				
    Int A [],				
    Int p [],				
    Int perm [],			
    double knobs [COLAMD_KNOBS],	
    Int stats [COLAMD_STATS],		
    void * (*allocate) (size_t, size_t),
    void (*release) (void *)
)
{
    Int *count ;		
    Int *mark ;			
    Int *M ;			
    size_t Mlen ;		
    Int n_row ;			
    Int nnz ;			
    Int i ;			
    Int j ;			
    Int k ;			
    Int mnz ;			
    Int pp ;			
    Int last_row ;		
    Int length ;		

    double cknobs [COLAMD_KNOBS] ;		
    double default_knobs [COLAMD_KNOBS] ;	


    if (!stats)
    {
	DEBUG0 (("symamd: stats not present\n")) ;
	return (FALSE) ;
    }
    for (i = 0 ; i < COLAMD_STATS ; i++)
    {
	stats [i] = 0 ;
    }
    stats [COLAMD_STATUS] = COLAMD_OK ;
    stats [COLAMD_INFO1] = -1 ;
    stats [COLAMD_INFO2] = -1 ;

    if (!A)
    {
    	stats [COLAMD_STATUS] = COLAMD_ERROR_A_not_present ;
	DEBUG0 (("symamd: A not present\n")) ;
	return (FALSE) ;
    }

    if (!p)		
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_p_not_present ;
	DEBUG0 (("symamd: p not present\n")) ;
    	return (FALSE) ;
    }

    if (n < 0)		
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_ncol_negative ;
	stats [COLAMD_INFO1] = n ;
	DEBUG0 (("symamd: n negative %d\n", n)) ;
    	return (FALSE) ;
    }

    nnz = p [n] ;
    if (nnz < 0)	
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_nnz_negative ;
	stats [COLAMD_INFO1] = nnz ;
	DEBUG0 (("symamd: number of entries negative %d\n", nnz)) ;
	return (FALSE) ;
    }

    if (p [0] != 0)
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_p0_nonzero ;
	stats [COLAMD_INFO1] = p [0] ;
	DEBUG0 (("symamd: p[0] not zero %d\n", p [0])) ;
	return (FALSE) ;
    }

    if (!knobs)
    {
	COLAMD_set_defaults (default_knobs) ;
	knobs = default_knobs ;
    }

    count = (Int *) ((*allocate) (n+1, sizeof (Int))) ;
    if (!count)
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
	DEBUG0 (("symamd: allocate count (size %d) failed\n", n+1)) ;
	return (FALSE) ;
    }

    mark = (Int *) ((*allocate) (n+1, sizeof (Int))) ;
    if (!mark)
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
	(*release) ((void *) count) ;
	DEBUG0 (("symamd: allocate mark (size %d) failed\n", n+1)) ;
	return (FALSE) ;
    }


    stats [COLAMD_INFO3] = 0 ;  

    for (i = 0 ; i < n ; i++)
    {
    	mark [i] = -1 ;
    }

    for (j = 0 ; j < n ; j++)
    {
	last_row = -1 ;

	length = p [j+1] - p [j] ;
	if (length < 0)
	{
	    stats [COLAMD_STATUS] = COLAMD_ERROR_col_length_negative ;
	    stats [COLAMD_INFO1] = j ;
	    stats [COLAMD_INFO2] = length ;
	    (*release) ((void *) count) ;
	    (*release) ((void *) mark) ;
	    DEBUG0 (("symamd: col %d negative length %d\n", j, length)) ;
	    return (FALSE) ;
	}

	for (pp = p [j] ; pp < p [j+1] ; pp++)
	{
	    i = A [pp] ;
	    if (i < 0 || i >= n)
	    {
		stats [COLAMD_STATUS] = COLAMD_ERROR_row_index_out_of_bounds ;
		stats [COLAMD_INFO1] = j ;
		stats [COLAMD_INFO2] = i ;
		stats [COLAMD_INFO3] = n ;
		(*release) ((void *) count) ;
		(*release) ((void *) mark) ;
		DEBUG0 (("symamd: row %d col %d out of bounds\n", i, j)) ;
		return (FALSE) ;
	    }

	    if (i <= last_row || mark [i] == j)
	    {
		stats [COLAMD_STATUS] = COLAMD_OK_BUT_JUMBLED ;
		stats [COLAMD_INFO1] = j ;
		stats [COLAMD_INFO2] = i ;
		(stats [COLAMD_INFO3]) ++ ;
		DEBUG1 (("symamd: row %d col %d unsorted/duplicate\n", i, j)) ;
	    }

	    if (i > j && mark [i] != j)
	    {
		count [i]++ ;
		count [j]++ ;
	    }

	    mark [i] = j ;

	    last_row = i ;
	}
    }

    perm [0] = 0 ;
    for (j = 1 ; j <= n ; j++)
    {
	perm [j] = perm [j-1] + count [j-1] ;
    }
    for (j = 0 ; j < n ; j++)
    {
	count [j] = perm [j] ;
    }


    mnz = perm [n] ;
    n_row = mnz / 2 ;
    Mlen = COLAMD_recommended (mnz, n_row, n) ;
    M = (Int *) ((*allocate) (Mlen, sizeof (Int))) ;
    DEBUG0 (("symamd: M is %d-by-%d with %d entries, Mlen = %g\n",
    	n_row, n, mnz, (double) Mlen)) ;

    if (!M)
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
	(*release) ((void *) count) ;
	(*release) ((void *) mark) ;
	DEBUG0 (("symamd: allocate M (size %g) failed\n", (double) Mlen)) ;
	return (FALSE) ;
    }

    k = 0 ;

    if (stats [COLAMD_STATUS] == COLAMD_OK)
    {
	for (j = 0 ; j < n ; j++)
	{
	    ASSERT (p [j+1] - p [j] >= 0) ;
	    for (pp = p [j] ; pp < p [j+1] ; pp++)
	    {
		i = A [pp] ;
		ASSERT (i >= 0 && i < n) ;
		if (i > j)
		{
		    M [count [i]++] = k ;
		    M [count [j]++] = k ;
		    k++ ;
		}
	    }
	}
    }
    else
    {
	DEBUG0 (("symamd: Duplicates in A.\n")) ;
	for (i = 0 ; i < n ; i++)
	{
	    mark [i] = -1 ;
	}
	for (j = 0 ; j < n ; j++)
	{
	    ASSERT (p [j+1] - p [j] >= 0) ;
	    for (pp = p [j] ; pp < p [j+1] ; pp++)
	    {
		i = A [pp] ;
		ASSERT (i >= 0 && i < n) ;
		if (i > j && mark [i] != j)
		{
		    M [count [i]++] = k ;
		    M [count [j]++] = k ;
		    k++ ;
		    mark [i] = j ;
		}
	    }
	}
    }

    (*release) ((void *) count) ;
    (*release) ((void *) mark) ;	
    ASSERT (k == n_row) ;


    for (i = 0 ; i < COLAMD_KNOBS ; i++)
    {
	cknobs [i] = knobs [i] ;
    }

    cknobs [COLAMD_DENSE_ROW] = -1 ;
    cknobs [COLAMD_DENSE_COL] = knobs [COLAMD_DENSE_ROW] ;

    (void) COLAMD_MAIN (n_row, n, (Int) Mlen, M, perm, cknobs, stats) ;

    stats [COLAMD_DENSE_ROW] = stats [COLAMD_DENSE_COL] ;

    (*release) ((void *) M) ;
    DEBUG0 (("symamd: done.\n")) ;
    return (TRUE) ;

}


Int COLAMD_MAIN		
(
    Int n_row,			
    Int n_col,			
    Int Alen,			
    Int A [],			
    Int p [],			
    double knobs [COLAMD_KNOBS],
    Int stats [COLAMD_STATS]	
)
{
    Int i ;			
    Int nnz ;			
    size_t Row_size ;		
    size_t Col_size ;		
    size_t need ;		
    Colamd_Row *Row ;		
    Colamd_Col *Col ;		
    Int n_col2 ;		
    Int n_row2 ;		
    Int ngarbage ;		
    Int max_deg ;		
    double default_knobs [COLAMD_KNOBS] ;	
    Int aggressive ;		
    int ok ;


    if (!stats)
    {
	DEBUG0 (("colamd: stats not present\n")) ;
	return (FALSE) ;
    }
    for (i = 0 ; i < COLAMD_STATS ; i++)
    {
	stats [i] = 0 ;
    }
    stats [COLAMD_STATUS] = COLAMD_OK ;
    stats [COLAMD_INFO1] = -1 ;
    stats [COLAMD_INFO2] = -1 ;

    // if (!A)		
    // {
	// stats [COLAMD_STATUS] = COLAMD_ERROR_A_not_present ;
	// DEBUG0 (("colamd: A not present\n")) ;
	// return (FALSE) ;
    // }

    // if (!p)		
    // {
	// stats [COLAMD_STATUS] = COLAMD_ERROR_p_not_present ;
	// DEBUG0 (("colamd: p not present\n")) ;
    // 	return (FALSE) ;
    // }

    // if (n_row < 0)	
    // {
	// stats [COLAMD_STATUS] = COLAMD_ERROR_nrow_negative ;
	// stats [COLAMD_INFO1] = n_row ;
	// DEBUG0 (("colamd: nrow negative %d\n", n_row)) ;
    // 	return (FALSE) ;
    // }

    // if (n_col < 0)	
    // {
	// stats [COLAMD_STATUS] = COLAMD_ERROR_ncol_negative ;
	// stats [COLAMD_INFO1] = n_col ;
	// DEBUG0 (("colamd: ncol negative %d\n", n_col)) ;
    // 	return (FALSE) ;
    // }

    nnz = p [n_col] ;
    if (nnz < 0)	
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_nnz_negative ;
	stats [COLAMD_INFO1] = nnz ;
	// DEBUG0 (("colamd: number of entries negative %d\n", nnz)) ;
	return (FALSE) ;
    }

    if (p [0] != 0)
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_p0_nonzero	;
	stats [COLAMD_INFO1] = p [0] ;
	// DEBUG0 (("colamd: p[0] not zero %d\n", p [0])) ;
	return (FALSE) ;
    }

    if (!knobs)
    {
	COLAMD_set_defaults (default_knobs) ;
	knobs = default_knobs ;
    }

    aggressive = (knobs [COLAMD_AGGRESSIVE] != FALSE) ;


    ok = TRUE ;
    Col_size = COLAMD_C (n_col, &ok) ;	    
    Row_size = COLAMD_R (n_row, &ok) ;	    

    need = t_mult (nnz, 2, &ok) ;
    need = t_add (need, n_col, &ok) ;
    need = t_add (need, Col_size, &ok) ;
    need = t_add (need, Row_size, &ok) ;

    if (!ok || need > (size_t) Alen || need > Int_MAX)
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_A_too_small ;
	stats [COLAMD_INFO1] = need ;
	stats [COLAMD_INFO2] = Alen ;
	// DEBUG0 (("colamd: Need Alen >= %d, given only Alen = %d\n", need,Alen));
	return (FALSE) ;
    }

    Alen -= Col_size + Row_size ;
    Col = (Colamd_Col *) &A [Alen] ;
    Row = (Colamd_Row *) &A [Alen + Col_size] ;

    if (!init_rows_cols (n_row, n_col, Row, Col, A, p, stats))
    {
	// DEBUG0 (("colamd: Matrix invalid\n")) ;
	return (FALSE) ;
    }

    init_scoring (n_row, n_col, Row, Col, A, p, knobs,
	&n_row2, &n_col2, &max_deg) ;

    ngarbage = find_ordering (n_row, n_col, Alen, Row, Col, A, p,
	n_col2, max_deg, 2*nnz, aggressive) ;


    order_children (n_col, Col, p) ;

    stats [COLAMD_DENSE_ROW] = n_row - n_row2 ;
    stats [COLAMD_DENSE_COL] = n_col - n_col2 ;
    stats [COLAMD_DEFRAG_COUNT] = ngarbage ;
    // DEBUG0 (("colamd: done.\n")) ; 
    return (TRUE) ;
}





PRIVATE Int init_rows_cols	
(

    Int n_row,			
    Int n_col,			
    Colamd_Row Row [],		
    Colamd_Col Col [],		
    Int A [],			
    Int p [],			
    Int stats [COLAMD_STATS]	
)
{

    Int col ;			
    Int row ;			
    Int *cp ;			
    Int *cp_end ;		
    Int *rp ;			
    Int *rp_end ;		
    Int last_row ;		

    for (col = 0 ; col < n_col ; col++)
    {
	Col [col].start = p [col] ;
	Col [col].length = p [col+1] - p [col] ;

	if (Col [col].length < 0)
	{
	    stats [COLAMD_STATUS] = COLAMD_ERROR_col_length_negative ;
	    stats [COLAMD_INFO1] = col ;
	    stats [COLAMD_INFO2] = Col [col].length ;
	    DEBUG0 (("colamd: col %d length %d < 0\n", col, Col [col].length)) ;
	    return (FALSE) ;
	}

	Col [col].shared1.thickness = 1 ;
	Col [col].shared2.score = 0 ;
	Col [col].shared3.prev = EMPTY ;
	Col [col].shared4.degree_next = EMPTY ;
    }



    stats [COLAMD_INFO3] = 0 ;	

    for (row = 0 ; row < n_row ; row++)
    {
	Row [row].length = 0 ;
	Row [row].shared2.mark = -1 ;
    }

    for (col = 0 ; col < n_col ; col++)
    {
	last_row = -1 ;

	cp = &A [p [col]] ;
	cp_end = &A [p [col+1]] ;

	while (cp < cp_end)
	{
	    row = *cp++ ;

	    if (row < 0 || row >= n_row)
	    {
		stats [COLAMD_STATUS] = COLAMD_ERROR_row_index_out_of_bounds ;
		stats [COLAMD_INFO1] = col ;
		stats [COLAMD_INFO2] = row ;
		stats [COLAMD_INFO3] = n_row ;
		DEBUG0 (("colamd: row %d col %d out of bounds\n", row, col)) ;
		return (FALSE) ;
	    }

	    if (row <= last_row || Row [row].shared2.mark == col)
	    {
		stats [COLAMD_STATUS] = COLAMD_OK_BUT_JUMBLED ;
		stats [COLAMD_INFO1] = col ;
		stats [COLAMD_INFO2] = row ;
		(stats [COLAMD_INFO3]) ++ ;
		DEBUG1 (("colamd: row %d col %d unsorted/duplicate\n",row,col));
	    }

	    if (Row [row].shared2.mark != col)
	    {
		Row [row].length++ ;
	    }
	    else
	    {
		Col [col].length-- ;
	    }

	    Row [row].shared2.mark = col ;

	    last_row = row ;
	}
    }


    Row [0].start = p [n_col] ;
    Row [0].shared1.p = Row [0].start ;
    Row [0].shared2.mark = -1 ;
    for (row = 1 ; row < n_row ; row++)
    {
	Row [row].start = Row [row-1].start + Row [row-1].length ;
	Row [row].shared1.p = Row [row].start ;
	Row [row].shared2.mark = -1 ;
    }

    if (stats [COLAMD_STATUS] == COLAMD_OK_BUT_JUMBLED)
    {
	for (col = 0 ; col < n_col ; col++)
	{
	    cp = &A [p [col]] ;
	    cp_end = &A [p [col+1]] ;
	    while (cp < cp_end)
	    {
		row = *cp++ ;
		if (Row [row].shared2.mark != col)
		{
		    A [(Row [row].shared1.p)++] = col ;
		    Row [row].shared2.mark = col ;
		}
	    }
	}
    }
    else
    {
	for (col = 0 ; col < n_col ; col++)
	{
	    cp = &A [p [col]] ;
	    cp_end = &A [p [col+1]] ;
	    while (cp < cp_end)
	    {
		A [(Row [*cp++].shared1.p)++] = col ;
	    }
	}
    }
    for (row = 0 ; row < n_row ; row++)
    {
	Row [row].shared2.mark = 0 ;
	Row [row].shared1.degree = Row [row].length ;
    }

    if (stats [COLAMD_STATUS] == COLAMD_OK_BUT_JUMBLED)
    {
    	DEBUG0 (("colamd: reconstructing column form, matrix jumbled\n")) ;


	Col [0].start = 0 ;
	p [0] = Col [0].start ;
	for (col = 1 ; col < n_col ; col++)
	{
	    Col [col].start = Col [col-1].start + Col [col-1].length ;
	    p [col] = Col [col].start ;
	}

	for (row = 0 ; row < n_row ; row++)
	{
	    rp = &A [Row [row].start] ;
	    rp_end = rp + Row [row].length ;
	    while (rp < rp_end)
	    {
		A [(p [*rp++])++] = row ;
	    }
	}
    }

    return (TRUE) ;
}



PRIVATE void init_scoring
(
    Int n_row,			
    Int n_col,			
    Colamd_Row Row [],		
    Colamd_Col Col [],		
    Int A [],			
    Int head [],		
    double knobs [COLAMD_KNOBS],
    Int *p_n_row2,		
    Int *p_n_col2,		
    Int *p_max_deg		
)
{
    Int c ;			
    Int r, row ;		
    Int *cp ;			
    Int deg ;			
    Int *cp_end ;		
    Int *new_cp ;		
    Int col_length ;		
    Int score ;			
    Int n_col2 ;		
    Int n_row2 ;		
    Int dense_row_count ;	
    Int dense_col_count ;	
    Int min_score ;		
    Int max_deg ;		
    Int next_col ;		


    if (knobs [COLAMD_DENSE_ROW] < 0)
    {
	dense_row_count = n_col-1 ;
    }
    else
    {
	dense_row_count = DENSE_DEGREE (knobs [COLAMD_DENSE_ROW], n_col) ;
    }
    if (knobs [COLAMD_DENSE_COL] < 0)
    {
	dense_col_count = n_row-1 ;
    }
    else
    {
	dense_col_count =
	    DENSE_DEGREE (knobs [COLAMD_DENSE_COL], MIN (n_row, n_col)) ;
    }

    DEBUG1 (("colamd: densecount: %d %d\n", dense_row_count, dense_col_count)) ;
    max_deg = 0 ;
    n_col2 = n_col ;
    n_row2 = n_row ;

    for (c = n_col-1 ; c >= 0 ; c--)
    {
	deg = Col [c].length ;
	if (deg == 0)
	{
	    Col [c].shared2.order = --n_col2 ;
	    KILL_PRINCIPAL_COL (c) ;
	}
    }
    DEBUG1 (("colamd: null columns killed: %d\n", n_col - n_col2)) ;

    for (c = n_col-1 ; c >= 0 ; c--)
    {
	if (COL_IS_DEAD (c))
	{
	    continue ;
	}
	deg = Col [c].length ;
	if (deg > dense_col_count)
	{
	    Col [c].shared2.order = --n_col2 ;
	    cp = &A [Col [c].start] ;
	    cp_end = cp + Col [c].length ;
	    while (cp < cp_end)
	    {
		Row [*cp++].shared1.degree-- ;
	    }
	    KILL_PRINCIPAL_COL (c) ;
	}
    }
    DEBUG1 (("colamd: Dense and null columns killed: %d\n", n_col - n_col2)) ;

    for (r = 0 ; r < n_row ; r++)
    {
	deg = Row [r].shared1.degree ;
	ASSERT (deg >= 0 && deg <= n_col) ;
	if (deg > dense_row_count || deg == 0)
	{
	    KILL_ROW (r) ;
	    --n_row2 ;
	}
	else
	{
	    max_deg = MAX (max_deg, deg) ;
	}
    }
    DEBUG1 (("colamd: Dense and null rows killed: %d\n", n_row - n_row2)) ;

    for (c = n_col-1 ; c >= 0 ; c--)
    {
	if (COL_IS_DEAD (c))
	{
	    continue ;
	}
	score = 0 ;
	cp = &A [Col [c].start] ;
	new_cp = cp ;
	cp_end = cp + Col [c].length ;
	while (cp < cp_end)
	{
	    row = *cp++ ;
	    if (ROW_IS_DEAD (row))
	    {
		continue ;
	    }
	    *new_cp++ = row ;
	    score += Row [row].shared1.degree - 1 ;
	    score = MIN (score, n_col) ;
	}
	col_length = (Int) (new_cp - &A [Col [c].start]) ;
	if (col_length == 0)
	{
	    DEBUG2 (("Newly null killed: %d\n", c)) ;
	    Col [c].shared2.order = --n_col2 ;
	    KILL_PRINCIPAL_COL (c) ;
	}
	else
	{
	    ASSERT (score >= 0) ;
	    ASSERT (score <= n_col) ;
	    Col [c].length = col_length ;
	    Col [c].shared2.score = score ;
	}
    }
    DEBUG1 (("colamd: Dense, null, and newly-null columns killed: %d\n",
    	n_col-n_col2)) ;




    for (c = 0 ; c <= n_col ; c++)
    {
	head [c] = EMPTY ;
    }
    min_score = n_col ;
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	if (COL_IS_ALIVE (c))
	{
	    DEBUG4 (("place %d score %d minscore %d ncol %d\n",
		c, Col [c].shared2.score, min_score, n_col)) ;

	    score = Col [c].shared2.score ;

	    ASSERT (min_score >= 0) ;
	    ASSERT (min_score <= n_col) ;
	    ASSERT (score >= 0) ;
	    ASSERT (score <= n_col) ;
	    ASSERT (head [score] >= EMPTY) ;

	    next_col = head [score] ;
	    Col [c].shared3.prev = EMPTY ;
	    Col [c].shared4.degree_next = next_col ;

	    if (next_col != EMPTY)
	    {
		Col [next_col].shared3.prev = c ;
	    }
	    head [score] = c ;

	    min_score = MIN (min_score, score) ;

	}
    }


    *p_n_col2 = n_col2 ;
    *p_n_row2 = n_row2 ;
    *p_max_deg = max_deg ;
}



PRIVATE Int find_ordering	
(
    Int n_row,			
    Int n_col,			
    Int Alen,			
    Colamd_Row Row [],		
    Colamd_Col Col [],		
    Int A [],			
    Int head [],		
    Int n_col2,			
    Int max_deg,		
    Int pfree,			 
    Int aggressive
)
{

    Int k ;			
    Int pivot_col ;		
    Int *cp ;			
    Int *rp ;			
    Int pivot_row ;		
    Int *new_cp ;		
    Int *new_rp ;		
    Int pivot_row_start ;	
    Int pivot_row_degree ;	
    Int pivot_row_length ;	
    Int pivot_col_score ;	
    Int needed_memory ;		
    Int *cp_end ;		
    Int *rp_end ;		
    Int row ;			
    Int col ;			
    Int max_score ;		
    Int cur_score ;		
    unsigned Int hash ;		
    Int head_column ;		
    Int first_col ;		
    Int tag_mark ;		
    Int row_mark ;		
    Int set_difference ;	
    Int min_score ;		
    Int col_thickness ;		
    Int max_mark ;		
    Int pivot_col_thickness ;	
    Int prev_col ;		
    Int next_col ;		
    Int ngarbage ;		

    max_mark = INT_MAX - n_col ;	
    tag_mark = clear_mark (0, max_mark, n_row, Row) ;
    min_score = 0 ;
    ngarbage = 0 ;
    DEBUG1 (("colamd: Ordering, n_col2=%d\n", n_col2)) ;

    for (k = 0 ; k < n_col2 ; )
    {


	ASSERT (min_score >= 0) ;
	ASSERT (min_score <= n_col) ;
	ASSERT (head [min_score] >= EMPTY) ;

	while (head [min_score] == EMPTY && min_score < n_col)
	{
	    min_score++ ;
	}
	pivot_col = head [min_score] ;
	ASSERT (pivot_col >= 0 && pivot_col <= n_col) ;
	next_col = Col [pivot_col].shared4.degree_next ;
	head [min_score] = next_col ;
	if (next_col != EMPTY)
	{
	    Col [next_col].shared3.prev = EMPTY ;
	}

	ASSERT (COL_IS_ALIVE (pivot_col)) ;

	pivot_col_score = Col [pivot_col].shared2.score ;

	Col [pivot_col].shared2.order = k ;

	pivot_col_thickness = Col [pivot_col].shared1.thickness ;
	k += pivot_col_thickness ;
	ASSERT (pivot_col_thickness > 0) ;
	DEBUG3 (("Pivot col: %d thick %d\n", pivot_col, pivot_col_thickness)) ;

	needed_memory = MIN (pivot_col_score, n_col - k) ;
	if (pfree + needed_memory >= Alen)
	{
	    pfree = garbage_collection (n_row, n_col, Row, Col, A, &A [pfree]) ;
	    ngarbage++ ;
	    ASSERT (pfree + needed_memory < Alen) ;
	    tag_mark = clear_mark (0, max_mark, n_row, Row) ;

	}

	pivot_row_start = pfree ;

	pivot_row_degree = 0 ;

	Col [pivot_col].shared1.thickness = -pivot_col_thickness ;

	cp = &A [Col [pivot_col].start] ;
	cp_end = cp + Col [pivot_col].length ;
	while (cp < cp_end)
	{
	    row = *cp++ ;
	    DEBUG4 (("Pivot col pattern %d %d\n", ROW_IS_ALIVE (row), row)) ;
	    if (ROW_IS_ALIVE (row))
	    {
		rp = &A [Row [row].start] ;
		rp_end = rp + Row [row].length ;
		while (rp < rp_end)
		{
		    col = *rp++ ;
		    col_thickness = Col [col].shared1.thickness ;
		    if (col_thickness > 0 && COL_IS_ALIVE (col))
		    {
			Col [col].shared1.thickness = -col_thickness ;
			ASSERT (pfree < Alen) ;
			A [pfree++] = col ;
			pivot_row_degree += col_thickness ;
		    }
		}
	    }
	}

	Col [pivot_col].shared1.thickness = pivot_col_thickness ;
	max_deg = MAX (max_deg, pivot_row_degree) ;

	cp = &A [Col [pivot_col].start] ;
	cp_end = cp + Col [pivot_col].length ;
	while (cp < cp_end)
	{
	    row = *cp++ ;
	    DEBUG3 (("Kill row in pivot col: %d\n", row)) ;
	    KILL_ROW (row) ;
	}

	pivot_row_length = pfree - pivot_row_start ;
	if (pivot_row_length > 0)
	{
	    pivot_row = A [Col [pivot_col].start] ;
	    DEBUG3 (("Pivotal row is %d\n", pivot_row)) ;
	}
	else
	{
	    pivot_row = EMPTY ;
	    ASSERT (pivot_row_length == 0) ;
	}
	ASSERT (Col [pivot_col].length > 0 || pivot_row_length == 0) ;


	DEBUG3 (("** Computing set differences phase. **\n")) ;


	DEBUG3 (("Pivot row: ")) ;
	rp = &A [pivot_row_start] ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    col = *rp++ ;
	    ASSERT (COL_IS_ALIVE (col) && col != pivot_col) ;
	    DEBUG3 (("Col: %d\n", col)) ;

	    col_thickness = -Col [col].shared1.thickness ;
	    ASSERT (col_thickness > 0) ;
	    Col [col].shared1.thickness = col_thickness ;

	    cur_score = Col [col].shared2.score ;
	    prev_col = Col [col].shared3.prev ;
	    next_col = Col [col].shared4.degree_next ;
	    ASSERT (cur_score >= 0) ;
	    ASSERT (cur_score <= n_col) ;
	    ASSERT (cur_score >= EMPTY) ;
	    if (prev_col == EMPTY)
	    {
		head [cur_score] = next_col ;
	    }
	    else
	    {
		Col [prev_col].shared4.degree_next = next_col ;
	    }
	    if (next_col != EMPTY)
	    {
		Col [next_col].shared3.prev = prev_col ;
	    }

	    cp = &A [Col [col].start] ;
	    cp_end = cp + Col [col].length ;
	    while (cp < cp_end)
	    {
		row = *cp++ ;
		row_mark = Row [row].shared2.mark ;
		if (ROW_IS_MARKED_DEAD (row_mark))
		{
		    continue ;
		}
		ASSERT (row != pivot_row) ;
		set_difference = row_mark - tag_mark ;
		if (set_difference < 0)
		{
		    ASSERT (Row [row].shared1.degree <= max_deg) ;
		    set_difference = Row [row].shared1.degree ;
		}
		set_difference -= col_thickness ;
		ASSERT (set_difference >= 0) ;
		if (set_difference == 0 && aggressive)
		{
		    DEBUG3 (("aggressive absorption. Row: %d\n", row)) ;
		    KILL_ROW (row) ;
		}
		else
		{
		    Row [row].shared2.mark = set_difference + tag_mark ;
		}
	    }
	}

	DEBUG3 (("** Adding set differences phase. **\n")) ;
	rp = &A [pivot_row_start] ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    col = *rp++ ;
	    ASSERT (COL_IS_ALIVE (col) && col != pivot_col) ;
	    hash = 0 ;
	    cur_score = 0 ;
	    cp = &A [Col [col].start] ;
	    new_cp = cp ;
	    cp_end = cp + Col [col].length ;

	    DEBUG4 (("Adding set diffs for Col: %d.\n", col)) ;

	    while (cp < cp_end)
	    {
		row = *cp++ ;
		ASSERT(row >= 0 && row < n_row) ;
		row_mark = Row [row].shared2.mark ;
		if (ROW_IS_MARKED_DEAD (row_mark))
		{
		    DEBUG4 ((" Row %d, dead\n", row)) ;
		    continue ;
		}
		DEBUG4 ((" Row %d, set diff %d\n", row, row_mark-tag_mark));
		ASSERT (row_mark >= tag_mark) ;
		*new_cp++ = row ;
		hash += row ;
		cur_score += row_mark - tag_mark ;
		cur_score = MIN (cur_score, n_col) ;
	    }

	    Col [col].length = (Int) (new_cp - &A [Col [col].start]) ;

	    if (Col [col].length == 0)
	    {
		DEBUG4 (("further mass elimination. Col: %d\n", col)) ;
		KILL_PRINCIPAL_COL (col) ;
		pivot_row_degree -= Col [col].shared1.thickness ;
		ASSERT (pivot_row_degree >= 0) ;
		Col [col].shared2.order = k ;
		k += Col [col].shared1.thickness ;
	    }
	    else
	    {

		DEBUG4 (("Preparing supercol detection for Col: %d.\n", col)) ;

		Col [col].shared2.score = cur_score ;

		hash %= n_col + 1 ;

		DEBUG4 ((" Hash = %d, n_col = %d.\n", hash, n_col)) ;
		ASSERT (((Int) hash) <= n_col) ;

		head_column = head [hash] ;
		if (head_column > EMPTY)
		{
		    first_col = Col [head_column].shared3.headhash ;
		    Col [head_column].shared3.headhash = col ;
		}
		else
		{
		    first_col = - (head_column + 2) ;
		    head [hash] = - (col + 2) ;
		}
		Col [col].shared4.hash_next = first_col ;

		Col [col].shared3.hash = (Int) hash ;
		ASSERT (COL_IS_ALIVE (col)) ;
	    }
	}

	DEBUG3 (("** Supercolumn detection phase. **\n")) ;

	detect_super_cols (


		Col, A, head, pivot_row_start, pivot_row_length) ;

	KILL_PRINCIPAL_COL (pivot_col) ;

	tag_mark = clear_mark (tag_mark+max_deg+1, max_mark, n_row, Row) ;


	DEBUG3 (("** Finalize scores phase. **\n")) ;

	rp = &A [pivot_row_start] ;
	new_rp = rp ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    col = *rp++ ;
	    if (COL_IS_DEAD (col))
	    {
		continue ;
	    }
	    *new_rp++ = col ;
	    A [Col [col].start + (Col [col].length++)] = pivot_row ;

	    cur_score = Col [col].shared2.score + pivot_row_degree ;

	    max_score = n_col - k - Col [col].shared1.thickness ;

	    cur_score -= Col [col].shared1.thickness ;

	    cur_score = MIN (cur_score, max_score) ;
	    ASSERT (cur_score >= 0) ;

	    Col [col].shared2.score = cur_score ;

	    ASSERT (min_score >= 0) ;
	    ASSERT (min_score <= n_col) ;
	    ASSERT (cur_score >= 0) ;
	    ASSERT (cur_score <= n_col) ;
	    ASSERT (head [cur_score] >= EMPTY) ;
	    next_col = head [cur_score] ;
	    Col [col].shared4.degree_next = next_col ;
	    Col [col].shared3.prev = EMPTY ;
	    if (next_col != EMPTY)
	    {
		Col [next_col].shared3.prev = col ;
	    }
	    head [cur_score] = col ;

	    min_score = MIN (min_score, cur_score) ;
	}

	if (pivot_row_degree > 0)
	{
	    Row [pivot_row].start  = pivot_row_start ;
	    Row [pivot_row].length = (Int) (new_rp - &A[pivot_row_start]) ;
	    ASSERT (Row [pivot_row].length > 0) ;
	    Row [pivot_row].shared1.degree = pivot_row_degree ;
	    Row [pivot_row].shared2.mark = 0 ;

	    DEBUG1 (("Resurrect Pivot_row %d deg: %d\n",
			pivot_row, pivot_row_degree)) ;
	}
    }

    return (ngarbage) ;
}



PRIVATE void order_children
(

    Int n_col,			
    Colamd_Col Col [],		
    Int p []			
)
{

    Int i ;			
    Int c ;			
    Int parent ;		
    Int order ;			


    for (i = 0 ; i < n_col ; i++)
    {
	ASSERT (COL_IS_DEAD (i)) ;
	if (!COL_IS_DEAD_PRINCIPAL (i) && Col [i].shared2.order == EMPTY)
	{
	    parent = i ;
	    do
	    {
		parent = Col [parent].shared1.parent ;
	    } while (!COL_IS_DEAD_PRINCIPAL (parent)) ;

	    c = i ;
	    order = Col [parent].shared2.order ;

	    do
	    {
		ASSERT (Col [c].shared2.order == EMPTY) ;

		Col [c].shared2.order = order++ ;
		Col [c].shared1.parent = parent ;
		c = Col [c].shared1.parent ;
	    } while (Col [c].shared2.order == EMPTY) ;

	    Col [parent].shared2.order = order ;
	}
    }

    for (c = 0 ; c < n_col ; c++)
    {
	p [Col [c].shared2.order] = c ;
    }
}



PRIVATE void detect_super_cols
(

    Colamd_Col Col [],		
    Int A [],			
    Int head [],		
    Int row_start,		
    Int row_length		
)
{

    Int hash ;			
    Int *rp ;			
    Int c ;			
    Int super_c ;		
    Int *cp1 ;			
    Int *cp2 ;			
    Int length ;		
    Int prev_c ;		
    Int i ;			
    Int *rp_end ;		
    Int col ;			
    Int head_column ;		
    Int first_col ;		

    rp = &A [row_start] ;
    rp_end = rp + row_length ;
    while (rp < rp_end)
    {
	col = *rp++ ;
	if (COL_IS_DEAD (col))
	{
	    continue ;
	}

	hash = Col [col].shared3.hash ;
	ASSERT (hash <= n_col) ;

	head_column = head [hash] ;
	if (head_column > EMPTY)
	{
	    first_col = Col [head_column].shared3.headhash ;
	}
	else
	{
	    first_col = - (head_column + 2) ;
	}


	for (super_c = first_col ; super_c != EMPTY ;
	    super_c = Col [super_c].shared4.hash_next)
	{
	    ASSERT (COL_IS_ALIVE (super_c)) ;
	    ASSERT (Col [super_c].shared3.hash == hash) ;
	    length = Col [super_c].length ;

	    prev_c = super_c ;

	    for (c = Col [super_c].shared4.hash_next ;
		 c != EMPTY ; c = Col [c].shared4.hash_next)
	    {
		ASSERT (c != super_c) ;
		ASSERT (COL_IS_ALIVE (c)) ;
		ASSERT (Col [c].shared3.hash == hash) ;

		if (Col [c].length != length ||
		    Col [c].shared2.score != Col [super_c].shared2.score)
		{
		    prev_c = c ;
		    continue ;
		}

		cp1 = &A [Col [super_c].start] ;
		cp2 = &A [Col [c].start] ;

		for (i = 0 ; i < length ; i++)
		{
		    ASSERT (ROW_IS_ALIVE (*cp1))  ;
		    ASSERT (ROW_IS_ALIVE (*cp2))  ;
		    if (*cp1++ != *cp2++)
		    {
			break ;
		    }
		}

		if (i != length)
		{
		    prev_c = c ;
		    continue ;
		}

		ASSERT (Col [c].shared2.score == Col [super_c].shared2.score) ;

		Col [super_c].shared1.thickness += Col [c].shared1.thickness ;
		Col [c].shared1.parent = super_c ;
		KILL_NON_PRINCIPAL_COL (c) ;
		Col [c].shared2.order = EMPTY ;
		Col [prev_c].shared4.hash_next = Col [c].shared4.hash_next ;
	    }
	}


	if (head_column > EMPTY)
	{
	    Col [head_column].shared3.headhash = EMPTY ;
	}
	else
	{
	    head [hash] = EMPTY ;
	}
    }
}



PRIVATE Int garbage_collection  
(
    Int n_row,			
    Int n_col,			
    Colamd_Row Row [],		
    Colamd_Col Col [],		
    Int A [],			
    Int *pfree			
)
{
    Int *psrc ;			
    Int *pdest ;		
    Int j ;			
    Int r ;			
    Int c ;			
    Int length ;		

    pdest = &A[0] ;
    for (c = 0 ; c < n_col ; c++)
    {
	if (COL_IS_ALIVE (c))
	{
	    psrc = &A [Col [c].start] ;

	    ASSERT (pdest <= psrc) ;
	    Col [c].start = (Int) (pdest - &A [0]) ;
	    length = Col [c].length ;
	    for (j = 0 ; j < length ; j++)
	    {
		r = *psrc++ ;
		if (ROW_IS_ALIVE (r))
		{
		    *pdest++ = r ;
		}
	    }
	    Col [c].length = (Int) (pdest - &A [Col [c].start]) ;
	}
    }

    for (r = 0 ; r < n_row ; r++)
    {
	if (ROW_IS_DEAD (r) || (Row [r].length == 0))
	{
	    KILL_ROW (r) ;
	}
	else
	{
	    psrc = &A [Row [r].start] ;
	    Row [r].shared2.first_column = *psrc ;
	    ASSERT (ROW_IS_ALIVE (r)) ;
	    *psrc = ONES_COMPLEMENT (r) ;
	}
    }

    psrc = pdest ;
    while (psrc < pfree)
    {
	if (*psrc++ < 0)
	{
	    psrc-- ;
	    r = ONES_COMPLEMENT (*psrc) ;
	    ASSERT (r >= 0 && r < n_row) ;
	    *psrc = Row [r].shared2.first_column ;
	    ASSERT (ROW_IS_ALIVE (r)) ;
	    ASSERT (Row [r].length > 0) ;
	    ASSERT (pdest <= psrc) ;
	    Row [r].start = (Int) (pdest - &A [0]) ;
	    length = Row [r].length ;
	    for (j = 0 ; j < length ; j++)
	    {
		c = *psrc++ ;
		if (COL_IS_ALIVE (c))
		{
		    *pdest++ = c ;
		}
	    }
	    Row [r].length = (Int) (pdest - &A [Row [r].start]) ;
	    ASSERT (Row [r].length > 0) ;
	}
    }
    ASSERT (debug_rows == 0) ;

    return ((Int) (pdest - &A [0])) ;
}

PRIVATE Int clear_mark	
(

    Int tag_mark,	
    Int max_mark,	

    Int n_row,		
    Colamd_Row Row []	
)
{

    Int r ;

    if (tag_mark <= 0 || tag_mark >= max_mark)
    {
	for (r = 0 ; r < n_row ; r++)
	{
	    if (ROW_IS_ALIVE (r))
	    {
		Row [r].shared2.mark = 0 ;
	    }
	}
	tag_mark = 1 ;
    }

    return (tag_mark) ;
}
