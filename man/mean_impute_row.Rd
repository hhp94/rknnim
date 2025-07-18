% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rknnim.R
\name{mean_impute_row}
\alias{mean_impute_row}
\title{Row Mean Imputation}
\usage{
mean_impute_row(obj)
}
\arguments{
\item{obj}{A numeric matrix where rows represent features and columns represent samples.}
}
\value{
The matrix with NA values replaced by row means.
}
\description{
Imputes missing values (NA) in a matrix by replacing them with the mean of their
respective rows, computed excluding NA values.
}
\examples{
\dontrun{
mat <- matrix(c(1, NA, 3, 4, 5, NA), nrow = 2)
imputed <- mean_impute_row(mat)
# Result: row 1 mean ( (1+3)/2 = 2 ) for NA; row 2 mean ( (4+5)/2 = 4.5 ) for NA
}
}
