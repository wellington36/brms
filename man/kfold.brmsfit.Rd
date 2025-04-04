% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kfold.R
\name{kfold.brmsfit}
\alias{kfold.brmsfit}
\alias{kfold}
\title{K-Fold Cross-Validation}
\usage{
\method{kfold}{brmsfit}(
  x,
  ...,
  K = 10,
  Ksub = NULL,
  folds = NULL,
  group = NULL,
  joint = FALSE,
  compare = TRUE,
  resp = NULL,
  model_names = NULL,
  save_fits = FALSE,
  recompile = NULL,
  future_args = list()
)
}
\arguments{
\item{x}{A \code{brmsfit} object.}

\item{...}{Further arguments passed to \code{\link{brm}}.}

\item{K}{The number of subsets of equal (if possible) size
into which the data will be partitioned for performing
\eqn{K}-fold cross-validation. The model is refit \code{K} times, each time
leaving out one of the \code{K} subsets. If \code{K} is equal to the total
number of observations in the data then \eqn{K}-fold cross-validation is
equivalent to exact leave-one-out cross-validation.}

\item{Ksub}{Optional number of subsets (of those subsets defined by \code{K})
to be evaluated. If \code{NULL} (the default), \eqn{K}-fold cross-validation
will be performed on all subsets. If \code{Ksub} is a single integer,
\code{Ksub} subsets (out of all \code{K}) subsets will be randomly chosen.
If \code{Ksub} consists of multiple integers or a one-dimensional array
(created via \code{as.array}) potentially of length one, the corresponding
subsets will be used. This argument is primarily useful, if evaluation of
all subsets is infeasible for some reason.}

\item{folds}{Determines how the subsets are being constructed.
Possible values are \code{NULL} (the default), \code{"stratified"},
\code{"grouped"}, or \code{"loo"}. May also be a vector of length
equal to the number of observations in the data. Alters the way
\code{group} is handled. More information is provided in the 'Details'
section.}

\item{group}{Optional name of a grouping variable or factor in the model.
What exactly is done with this variable depends on argument \code{folds}.
More information is provided in the 'Details' section.}

\item{joint}{Indicates which observations' log likelihoods shall be
considered jointly in the ELPD computation. If \code{"obs"} or \code{FALSE}
(the default), each observation is considered separately. This enables
comparability of \code{kfold} with \code{loo}. If \code{"fold"} or
\code{TRUE}, the joint log likelihoods per fold are used. If
\code{"group"}, the joint log likelihoods per group within folds are used
(only available if argument \code{group} is specified).}

\item{compare}{A flag indicating if the information criteria
of the models should be compared to each other
via \code{\link{loo_compare}}.}

\item{resp}{Optional names of response variables. If specified, predictions
are performed only for the specified response variables.}

\item{model_names}{If \code{NULL} (the default) will use model names
derived from deparsing the call. Otherwise will use the passed
values as model names.}

\item{save_fits}{If \code{TRUE}, a component \code{fits} is added to
the returned object to store the cross-validated \code{brmsfit}
objects and the indices of the omitted observations for each fold.
Defaults to \code{FALSE}.}

\item{recompile}{Logical, indicating whether the Stan model should be
recompiled. This may be necessary if you are running \code{reloo} on
another machine than the one used to fit the model.}

\item{future_args}{A list of further arguments passed to
\code{\link[future:future]{future}} for additional control over parallel
execution if activated.}
}
\value{
\code{kfold} returns an object that has a similar structure as the
  objects returned by the \code{loo} and \code{waic} methods and
  can be used with the same post-processing functions.
}
\description{
Perform exact K-fold cross-validation by refitting the model \eqn{K}
times each leaving out one-\eqn{K}th of the original data.
Folds can be run in parallel using the \pkg{future} package.
}
\details{
The \code{kfold} function performs exact \eqn{K}-fold
  cross-validation. First the data are partitioned into \eqn{K} folds
  (i.e. subsets) of equal (or as close to equal as possible) size by default.
  Then the model is refit \eqn{K} times, each time leaving out one of the
  \code{K} subsets. If \eqn{K} is equal to the total number of observations
  in the data then \eqn{K}-fold cross-validation is equivalent to exact
  leave-one-out cross-validation (to which \code{loo} is an efficient
  approximation). The \code{compare_ic} function is also compatible with
  the objects returned by \code{kfold}.

  The subsets can be constructed in multiple different ways:
  \itemize{
  \item If both \code{folds} and \code{group} are \code{NULL}, the subsets
  are randomly chosen so that they have equal (or as close to equal as
  possible) size.
  \item If \code{folds} is \code{NULL} but \code{group} is specified, the
  data is split up into subsets, each time omitting all observations of one
  of the factor levels, while ignoring argument \code{K}.
  \item If \code{folds = "stratified"} the subsets are stratified after
  \code{group} using \code{\link[loo:kfold-helpers]{loo::kfold_split_stratified}}.
  \item If \code{folds = "grouped"} the subsets are split by
  \code{group} using \code{\link[loo:kfold-helpers]{loo::kfold_split_grouped}}.
  \item If \code{folds = "loo"} exact leave-one-out cross-validation
  will be performed and \code{K} will be ignored. Further, if \code{group}
  is specified, all observations corresponding to the factor level of the
  currently predicted single value are omitted. Thus, in this case, the
  predicted values are only a subset of the omitted ones.
  \item If \code{folds} is a numeric vector, it must contain one element per
  observation in the data. Each element of the vector is an integer in
  \code{1:K} indicating to which of the \code{K} folds the corresponding
  observation belongs. There are some convenience functions available in
  the \pkg{loo} package that create integer vectors to use for this purpose
  (see the Examples section below and also the
  \link[loo:kfold-helpers]{kfold-helpers} page).
  }

  When running \code{kfold} on a \code{brmsfit} created with the
  \pkg{cmdstanr} backend in a different \R session, several recompilations
  will be triggered because by default, \pkg{cmdstanr} writes the model
  executable to a temporary directory. To avoid that, set option
  \code{"cmdstanr_write_stan_file_dir"} to a nontemporary path of your choice
  before creating the original \code{brmsfit} (see section 'Examples' below).
}
\examples{
\dontrun{
fit1 <- brm(count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
           data = epilepsy, family = poisson())
# throws warning about some pareto k estimates being too high
(loo1 <- loo(fit1))
# perform 10-fold cross validation
(kfold1 <- kfold(fit1, chains = 1))

# use joint likelihoods per fold for ELPD evaluation
kfold(fit1, chains = 1, joint = "fold")

# use the future package for parallelization of models
# that is to fit models belonging to different folds in parallel
library(future)
plan(multisession, workers = 4)
kfold(fit1, chains = 1)
plan(sequential)

## to avoid recompilations when running kfold() on a 'cmdstanr'-backend fit
## in a fresh R session, set option 'cmdstanr_write_stan_file_dir' before
## creating the initial 'brmsfit'
## CAUTION: the following code creates some files in the current working
## directory: two 'model_<hash>.stan' files, one 'model_<hash>(.exe)'
## executable, and one 'fit_cmdstanr_<some_number>.rds' file
set.seed(7)
fname <- paste0("fit_cmdstanr_", sample.int(.Machine$integer.max, 1))
options(cmdstanr_write_stan_file_dir = getwd())
fit_cmdstanr <- brm(rate ~ conc + state, data = Puromycin,
                    backend = "cmdstanr", file = fname)

# now restart the R session and run the following (after attaching 'brms')
set.seed(7)
fname <- paste0("fit_cmdstanr_", sample.int(.Machine$integer.max, 1))
fit_cmdstanr <- brm(rate ~ conc + state,
                    data = Puromycin,
                    backend = "cmdstanr",
                    file = fname)
kfold_cmdstanr <- kfold(fit_cmdstanr, K = 2)
}

}
\seealso{
\code{\link{loo}}, \code{\link{reloo}}
}
