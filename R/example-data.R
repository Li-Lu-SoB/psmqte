#' Load an example data set for psmqte
#'
#' Loads an example data set suitable for demonstrating propensity-score
#' matching QTE/QTT estimation. The default `"demo"` source is a fixed synthetic
#' data set shipped inside this package. The optional `"nhefs"` source uses the
#' public NHEFS example data from the `causaldata` package without redistributing
#' those data inside `psmqte`.
#'
#' @param source One of `"demo"`, `"simulated"`, or `"nhefs"`.
#' @param n Sample size for the simulated example.
#' @param seed Optional random seed for the simulated example.
#'
#' @return A data frame with outcome `Y`, treatment `D`, and covariates.
#' @export
psmqte_example_data <- function(source = c("demo", "simulated", "nhefs"),
                                n = 1000L, seed = 1L) {
  source <- match.arg(source)
  if (source == "demo") {
    path <- system.file("extdata", "psmqte_demo.csv", package = "psmqte",
                        mustWork = TRUE)
    return(utils::read.csv(path, stringsAsFactors = FALSE))
  }

  if (source == "simulated") {
    return(simulate_psmqte_data(n = n, seed = seed))
  }

  if (!requireNamespace("causaldata", quietly = TRUE)) {
    stop("Install the causaldata package to use source = 'nhefs': ",
         "install.packages('causaldata')", call. = FALSE)
  }

  env <- new.env(parent = emptyenv())
  utils::data("nhefs_complete", package = "causaldata", envir = env)
  if (!exists("nhefs_complete", envir = env, inherits = FALSE)) {
    stop("Could not load causaldata::nhefs_complete.", call. = FALSE)
  }
  dat <- as.data.frame(env$nhefs_complete)
  keep <- c(
    "wt82_71", "qsmk", "age", "sex", "race", "school", "smokeintensity",
    "smokeyrs", "exercise", "active", "wt71"
  )
  missing <- setdiff(keep, names(dat))
  if (length(missing) > 0L) {
    stop("The installed causaldata::nhefs_complete data set is missing: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- dat[, keep]
  names(out)[names(out) == "wt82_71"] <- "Y"
  names(out)[names(out) == "qsmk"] <- "D"
  out
}
