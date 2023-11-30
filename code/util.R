# Make a rectangular matrix with n_row rows and n_col columns, where the upper
# triangular elements are ordered along their rows. That is, for row 1, the
# first element is distributed as N(0, 1) (unconstrained), but the value in the
# second column is guaranteed to be larger that in the first column, the third
# larger than in the second column, etc. The diagonals are unconstrained. Note
# that this is constructed using an improper probability distribution over the
# ordered variables, and so cannot be sampled a priori in greta.
triangular_ordered_matrix <- function(n_row, n_col) {
    
    # matrix to populate
    x <- zeros(n_row, n_col)
    
    # indices
    upper <- upper.tri(x, diag = FALSE)
    diag <- row(x) == col(x)
    
    # pull out the diagonal elements for ordered rows to add onto these orders
    bases <- x[diag]
    
    # create a matrix with ordered variables in the upper triangular
    ordered <- zeros(n_row, n_col)
    n_ordered_rows <- n_col - 1
    for (row in seq_len(n_ordered_rows)) {
        cols <- (row + 1):n_col
        n_col_fill <- length(cols)
        # make an ordered (or positive-constrained if only one) variable to add
        if (n_col_fill > 1) {
            pos <- ordered_variable(n_col_fill)
        } else {
            pos <- variable(lower = 0)
        }
        ordered[row, cols] <- bases[row] + pos
    }
    
    # put free parameters on the diagonal and lower triangular parts
    x[!upper] <- normal(0, 1, dim = sum(!upper))
    # and ordered on the upper triangular part
    x[upper] <- ordered[upper]
    x
    
}
