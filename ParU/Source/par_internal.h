typedef struct /* Tuple */
{
    /* The (e,f) tuples for element lists */
    int e,   /*  element number */
        f;  /*   offest */
} Tuple;

typedef struct /* Element */
{
    int 
        nrows,
        ncols;
    /* followed by 
     * int col[0..ncols-1], //column indices of the element
     * int row[0..nrows-1]; row indices of the element
     * double C[0..nrows*ncols-1] * Numerical values
     */
} Element
