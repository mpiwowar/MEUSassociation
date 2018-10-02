#' Package MEUSAssociation implements Z association coefficient.
#' 
#' Z coefficient \code{Z[A,B]} for two event \code{A} and \code{B} proposed by 
#' Z. Meus is defined by the following formula: \deqn{Z(A,B)^2 = 1 - 
#' [P(B)*(1-P(A|B))/(1-P(A))*(1-P(A^c|B))/(1-P(A^c)) + 
#' P(B^c)*(1-P(A|B^c))/(1-P(A))*(1-P(A^c|B^c))/(1-P(A^c))].} It is further 
#' generalized for a more generic case of contingency tables . More information
#' are available in the paper that is prepared for publication and the citation
#' will be added in the new version of the package as soon as it boceoms
#' available.
#' 
#' Packages provides the following set of functions for calculating the Z 
#' association coefficient. \itemize{ \item{\code{\link{z_coefficient}} Generic 
#' version of Z coefficient } \item{\code{\link{z_coefficient_matrix}} Matrix 
#' version of Z coefficient} \item{\code{\link{z_coefficient_ranks}} Ordered Z 
#' coefficients} }
#' 
#' The package have also the following sample datasets: \itemize{ 
#' \item{\code{\link{cancer_mutations}} Gene mutations in different types of 
#' cancer} \item{\code{\link{cancer_mutations_gene_groups}} Gene groups defined 
#' for gene mutations in different types of cancer} }
"_PACKAGE"
