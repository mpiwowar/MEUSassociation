% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ZCoefficientFunctions.R
\name{z_coefficient}
\alias{z_coefficient}
\title{Z coefficient evaluation generic function}
\usage{
z_coefficient(m, col_groups = NULL, row_groups = NULL)
}
\arguments{
\item{m}{A matrix. It can be a contingency table with counts or a table with 
probabilities. It is expected that there are no negative entries in a 
matrix and that mariginal sums are positive for all columns and rows.}

\item{col_groups, row_groups}{A vector or factor giving the grouping, with one
element per column (row) of \code{m}. Missing values will be treated as 
another group and a warning will be given. When a row (column) grouping is 
specified it is assumed that mariginal sums for each group is positive.}
}
\value{
The single Z coefficient, calculated for  matrix \code{m}.
}
\description{
\code{z_coefficient} returns the Z coefficient calculated for a given matrix
\code{M} (contingency table).
}
\details{
Z coefficient \code{Z[A_i,B_j]} for a matrix \code{m} is defined by the 
following formula: \deqn{Z(A_i,B_j)^2 = 1 - 
[P(B_j)*P(A_i|B_j)/P(A_i)*P(A_i^c|B_j)/P(A_i^c) + 
P(B_j^c)*P(A_i|B_j^c)/P(A_i)*P(A_i^c|B_j^c)/P(A_i^c)],} where \code{A_i} 
corresponds to the row \code{i} and \code{B_j} corresponds to the column 
\code{j}. It is further generalized for more rows (row groups) and columns 
(column groups). For more information see \url{link to the articile to be
added once published}

This function calculates a generic version of Z coefficient for a matrix 
\code{m} calculated for all rows (row groups) and columns (column groups).
}
\examples{
#load gene mutations in cancer dataset together with gene groups
data("cancer_mutations")
data("cancer_mutations_gene_groups")

#calculating Z coefficient for the whole matrix, also using gene groups
z_coefficient(cancer_mutations)
z_coefficient(cancer_mutations, col_groups = cancer_mutations_gene_groups)
}
\seealso{
Other Z coefficient: \code{\link{z_coefficient_matrix}},
  \code{\link{z_coefficient_ranks}}
}
\concept{
Z,coefficient,association
}
