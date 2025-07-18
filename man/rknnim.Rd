% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rknnim.R
\name{rknnim}
\alias{rknnim}
\title{Rolling KNN Imputation}
\usage{
rknnim(
  obj,
  n.feat,
  n.overlap = 10,
  coord = NULL,
  k = 10,
  rowmax = 0.5,
  colmax = 0.8,
  rng.seed = 362436069,
  .progress = FALSE,
  .parallel = FALSE
)
}
\arguments{
\item{obj}{A numeric matrix where rows represent features and columns represent samples.}

\item{n.feat}{Integer specifying the number of features (rows) per imputation window.}

\item{n.overlap}{Integer specifying the number of features of overlap
between consecutive windows (default: 10).}

\item{coord}{Optional coordinate matrix (default: NULL). Currently unused in the function.}

\item{k}{Number of neighbors to be used in the
    imputation (default=10)}

\item{rowmax}{The maximum percent missing data allowed in any row
    (default 50\%). For any rows with more than \code{rowmax}\% missing
    are imputed using the overall mean per sample.}

\item{colmax}{The maximum percent missing data allowed in any column
    (default 80\%). If any column has more than \code{colmax}\% missing data,
    the program halts and reports an error.}

\item{rng.seed}{The seed used for the random number generator (default
    362436069) for reproducibility.}

\item{.progress}{Whether to show a progress bar. Use \code{TRUE} to turn on
a basic progress bar, use a string to give it a name, or see
\link[purrr]{progress_bars} for more details.}

\item{.parallel}{Logical indicating whether to use parallel processing for imputation
(default: FALSE). Requires the 'mirai' package if enabled.}
}
\value{
A numeric matrix of the same dimensions as \code{obj} with missing values imputed.
}
\description{
Performs rolling window KNN imputation on a matrix with missing values.
The function divides the input matrix into overlapping windows along the rows,
imputes each window using KNN imputation, and combines the results by averaging
values in overlapping regions. Any remaining NA values are imputed with row means.
}
\examples{
\dontrun{
# Example matrix with missing values
set.seed(123)
mat <- matrix(rnorm(100), nrow = 20, ncol = 5)
mat[sample(100, 10)] <- NA
imputed <- rknnim(mat, n.feat = 10, n.overlap = 2, k = 5)
}
}
