#' Z coefficient evaluation function in a matrix version
#' 
#' The function \code{z_coefficient_matrix} returns a matrix of Z coefficients.
#' 
#' Z coefficient \code{Z[A_i,B_j]} for a matrix \code{m} is defined by the
#' following formula: \deqn{Z(A_i,B_j)^2 = 1 -
#' [P(B_j)*P(A_i|B_j)/P(A_i)*P(A_i^c|B_j)/P(A_i^c) + 
#' P(B_j^c)*P(A_i|B_j^c)/P(A_i)*P(A_i^c|B_j^c)/P(A_i^c)],} where \code{A_i}
#' corresponds to the row \code{i} and \code{B_j} corresponds to the column
#' \code{j}.
#' 
#' The function \code{z_coefficient_matrix} calculates a matrix, where each entry
#' \code{m[i,j]} of a matrix \code{m} is the Z coefficient \code{Z(A_i,B_j)}.
#' 
#' @param m A matrix. It can be a contingency table with counts or a table with
#'   probabilities. It is expected that there are no negative entries in a
#'   matrix and that mariginal sums are positive for all columns and rows.
#' @param col_groups,row_groups A vector or factor giving the grouping, with one
#'   element per column (row) of \code{m}. Missing values will be treated as
#'   another group and a warning will be given. When a row (column) grouping is
#'   specified it is assumed that mariginal sums for each group is positive.
#' @return A matrix of Z coefficients.
#' @concept Z,coefficient,association
#' @family Z coefficient
#' @examples
#' #load gene mutations in cancer dataset together with gene groups
#' data("cancer_mutations")
#' data("cancer_mutations_gene_groups")
#' 
#' #calculating Z coefficients in a matrix form, also using gene groups
#' #head used to show only top 6 rows
#' head(z_coefficient_matrix(cancer_mutations))
#' head(z_coefficient_matrix(cancer_mutations, col_groups = cancer_mutations_gene_groups))
#' @export
z_coefficient_matrix <-function(m, col_groups = NULL, row_groups = NULL){
  m = as.matrix(m)
  
  if(length(which(m < 0)) > 0) warning("Expecting non negative values in the matrix m")
  
  if(!is.null(col_groups)) m = t(rowsum(t(m), col_groups))
  if(!is.null(row_groups)) m = rowsum(m, row_groups)
  
  s <- sum(m)
  pa = apply(m, 1, sum)/s
  pb = apply(m, 2, sum)/s
  if(length(which(pa == 0))>0) warning("Expecting positive marginal sums for rows (row groups) in the matrix m")
  if(length(which(pb == 0))>0) warning("Expecting positive marginal sums for columns (column groups) in the matrix m")
  
  z <- matrix(ncol=ncol(m), nrow=nrow(m))
  colnames(z)<-colnames(m)
  rownames(z)<-rownames(m)
  
  #Z(A:B)^2 = 1 - [pb * PAB/pa * PA_B/PA_) + PB_ * PAB_/pa * PA_B_/PA_] = 1 - [pb * PAB/pa * (1-PAB)/(1-pa) + (1-pb) * PAB_/pa * (1-PAB_)/(1-pa)]
  #using simplified formula Z(A:B) = (PAB - pa * pb)/(pa * PA_ * pb * PB_)
  for (i in 1: nrow (m)) {
    for (j in 1: ncol(m)){
      z[i, j] = abs(m[i,j]/s - pa[i] * pb[j])/sqrt(pa[i] * (1 - pa[i]) * pb[j] * (1 - pb[j]))
    }
  }
  
  z
}

#' Ordered Z coefficient evaluation function
#' 
#' The function \code{z_coefficient_ranks} returns ordered Z coefficients 
#' calculated for each entry of a matrix \code{m}.
#' 
#' Z coefficient \code{Z[A_i,B_j]} for a matrix \code{m} is defined by the 
#' following formula: \deqn{Z(A_i,B_j)^2 = 1 - 
#' [P(B_j)*P(A_i|B_j)/P(A_i)*P(A_i^c|B_j)/P(A_i^c) + 
#' P(B_j^c)*P(A_i|B_j^c)/P(A_i)*P(A_i^c|B_j^c)/P(A_i^c)],} where \code{A_i} 
#' corresponds to the row \code{i} and \code{B_j} corresponds to the column 
#' \code{j}.
#' 
#' The function \code{z_coefficient_ranks} calculates the Z coefficient for each 
#' entry \code{m[i,i]} of a matrix \code{m} and returns an ordered list of Z 
#' coefficients.
#' 
#' @param m A matrix. It can be a contingency table with counts or a table with 
#'   probabilities. It is expected that there are no negative entries in a 
#'   matrix and that mariginal sums are positive for all columns and rows.
#' @param col_groups,row_groups A vector or factor giving the grouping, with one 
#'   element per column (row) of \code{m}. Missing values will be treated as 
#'   another group and a warning will be given. When a row (column) grouping is 
#'   specified it is assumed that mariginal sums for each group is positive.
#' @return A data frame containing an ordered list of Z coefficients. The 
#'   resulting data frame has three columns, specifying rows (or row groups, if 
#'   they are specified), columns (or column groups, it they are specified) and 
#'   Z coefficients. The data frame is sorted according to descending order of Z
#'   coefficients.
#' @concept Z,coefficient,association
#' @family Z coefficient
#' @examples
#' #load gene mutations in cancer dataset together with gene groups
#' data("cancer_mutations")
#' data("cancer_mutations_gene_groups")
#' 
#' #calculating ranked Z coefficients, also using gene groups
#' #head used to show only top 6 rows
#' head(z_coefficient_ranks(cancer_mutations))
#' head(z_coefficient_ranks(cancer_mutations, col_groups = cancer_mutations_gene_groups))
#' @export
z_coefficient_ranks <-function(m, col_groups = NULL, row_groups = NULL){
  m = as.matrix(m)
  
  if(length(which(m < 0)) > 0) warning("Expecting non negative values in the matrix m")
  
  if(!is.null(col_groups)) m = t(rowsum(t(m), col_groups))
  if(!is.null(row_groups)) m = rowsum(m, row_groups)
  
  nr = nrow(m)
  nc = ncol(m)

  s <- sum(m)
  pa = apply(m, 1, sum)/s
  pb = apply(m, 2, sum)/s
  if(length(which(pa == 0))>0) warning("Expecting positive marginal sums for rows (row groups) in the matrix m")
  if(length(which(pb == 0))>0) warning("Expecting positive marginal sums for columns (column groups) in the matrix m")
  
  if(is.null(rownames(m))) rows = 1:nr else rows = rownames(m)

  if(is.null(colnames(m))) cols = 1:nc else cols = colnames(m)

  z = data.frame(
    row = rep(rows, each = nc),
    col = rep(cols, nr),
    z_coeff = rep(NA, nr*nc)
  )
  
  #Z(A:B)^2 = 1 - [pb * PAB/pa * PA_B/PA_) + PB_ * PAB_/pa * PA_B_/PA_] = 1 - [pb * PAB/pa * (1-PAB)/(1-pa) + (1-pb) * PAB_/pa * (1-PAB_)/(1-pa)]
  #using simplified formula Z(A:B) = (PAB - pa * pb)/(pa * PA_ * pb * PB_)
  for (i in 1: nr) {
    for (j in 1: nc){
      z$z_coeff[(i-1)*nc + j] = abs(m[i,j]/s - pa[i] * pb[j])/sqrt(pa[i] * (1 - pa[i]) * pb[j] * (1 - pb[j]))
    }
  }
  
  z.sorted = z[order(-z$z_coeff),]
  
  rownames(z.sorted) = 1:(nr*nc)
  
  z.sorted
}

#' Gene mutations in cancers
#' 
#' The \code{cancer_mutations} is a matrix (contingency table) representing 
#' different gene mutations in different cancer types.
#' 
#' Sample dataset was taken from (Kandoth et al., (2013) Mutational
#' landscape and significance across 12 major cancer types. Nature, 502,
#' 333–339.) and consist of information related to mutated genes (with
#' point mutations and small insertions/deletions) from 3281 tumours across 12
#' tumour types, groupedd into 11 categories.
#' Rows correspond to different tumore types and columns to different genes.
#' 
#' @docType data
#' @usage data(cancer_mutations)
#' @format A data frame representing coningency table.
#' @family cancer mutations
#' @source The \code{cancer_mutations} is taken from the following article 
#'   \url{https://www.nature.com/articles/nature12634}.
"cancer_mutations"

#' Gene groups for gene mutations in cancers
#' 
#' The \code{cancer_mutations_gene_groups} is a vector of factors specifying gene
#' groups for different gene mutations in different cancer types stored in
#' \code{cancer_mutations}.
#' 
#' The factor vector \code{cancer_mutations_gene_groups} defines a biochemical
#' process in which particular gene is taking part.
#' 
#' @docType data
#' @usage data(cancer_mutations_gene_groups)
#' @format A factor vector
#' @family cancer mutations
#' @source Information regarding biochemical processes in which genes are taking
#'   part is taken from the following article/database/website: \url{to be supplemented}.
"cancer_mutations_gene_groups"

#' Z coefficient evaluation generic function
#' 
#' \code{z_coefficient} returns the Z coefficient calculated for a given matrix
#' \code{M} (contingency table).
#' 
#' Z coefficient \code{Z[A_i,B_j]} for a matrix \code{m} is defined by the 
#' following formula: \deqn{Z(A_i,B_j)^2 = 1 - 
#' [P(B_j)*P(A_i|B_j)/P(A_i)*P(A_i^c|B_j)/P(A_i^c) + 
#' P(B_j^c)*P(A_i|B_j^c)/P(A_i)*P(A_i^c|B_j^c)/P(A_i^c)],} where \code{A_i} 
#' corresponds to the row \code{i} and \code{B_j} corresponds to the column 
#' \code{j}. It is further generalized for more rows (row groups) and columns 
#' (column groups). For more information see \url{link to the articile to be
#' added once published}
#' 
#' This function calculates a generic version of Z coefficient for a matrix 
#' \code{m} calculated for all rows (row groups) and columns (column groups).
#' 
#' @param m A matrix. It can be a contingency table with counts or a table with 
#'   probabilities. It is expected that there are no negative entries in a 
#'   matrix and that mariginal sums are positive for all columns and rows.
#' @param col_groups,row_groups A vector or factor giving the grouping, with one
#'   element per column (row) of \code{m}. Missing values will be treated as 
#'   another group and a warning will be given. When a row (column) grouping is 
#'   specified it is assumed that mariginal sums for each group is positive.
#' @return The single Z coefficient, calculated for  matrix \code{m}.
#' @concept Z,coefficient,association
#' @family Z coefficient
#' @examples
#' #load gene mutations in cancer dataset together with gene groups
#' data("cancer_mutations")
#' data("cancer_mutations_gene_groups")
#' 
#' #calculating Z coefficient for the whole matrix, also using gene groups
#' z_coefficient(cancer_mutations)
#' z_coefficient(cancer_mutations, col_groups = cancer_mutations_gene_groups)
#' @export
z_coefficient = function(m, col_groups = NULL, row_groups = NULL){
  m = as.matrix(m)
  
  if(length(which(m < 0)) > 0) warning("Expecting non negative values in the matrix m")
  
  if(!is.null(col_groups)) m = t(rowsum(t(m), col_groups))
  if(!is.null(row_groups)) m = rowsum(m, row_groups)
  
  nr = nrow(m)
  nc = ncol(m)
  
  s = sum(m)
  P = m/s
  pa = apply(m, 1, sum)/s
  pb = apply(m, 2, sum)/s

  if(length(which(pa == 0))>0) warning("Expecting positive marginal sums for rows (row groups) in the matrix m")
  if(length(which(pb == 0))>0) warning("Expecting positive marginal sums for columns (column groups) in the matrix m")
  
  pa1 = matrix(rep(pa, nc), nrow = nr, byrow = FALSE)
  pb1 = matrix(rep(pb, nr), nrow = nr, byrow = TRUE)
  z1 = (1-P/pb1)/(1-pa1)
  z2 = apply(z1, 2, prod)
  z3 = sum(z2*pb)
  #correcting rounding errors
  if(z3>1) z3=1
  if(z3<0) z3=0
  z = sqrt(1 - z3)
  z
}
