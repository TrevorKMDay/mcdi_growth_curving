gomp2.fit <- function(data,
                      response_var = "inventory", t_var = "age",
                      min = 0, max = 681,
                      max.iter = 50) {

  # Assign max to age
  if (!is.na(max)) {
    A <- max
  } else {
    A <- max(data[, response_var])
    message(paste("Using max:", A))
  }

  # Assign min to W_0, but substitute extremely small value if 0 is supplied
  if (min == 0) {
    W_0 <- .Machine$double.eps
  } else {
    W_0 <- min
  }

  the_formula <- paste(response_var,
                       "~ W_0 * (A / W_0) ^ (1 - exp(-k_g *", t_var, "))") %>%
    as.formula()

  # This is the formula that works best to solve in R
  # 19 refers to its number in the Tjorve and Tjorve paper
  fit19 <- NULL
  try( fit19 <- nls(the_formula,
                    data = data,
                    start = list(k_g = .1),
                    control = list(maxiter = max.iter)) )

  if (is.null(fit19))
    fit19 <- NA

  return(fit19)

}

extract.kg <- function(fits, omit.NA = FALSE) {

  kgs <- sapply(fits, function(x)
      # If list element is logical it is NA/missing
      # Only extract coefs from real objects
      ifelse(typeof(x) == "logical", NA, summary(x)$coefficients[1, 1])
    )

  if (omit.NA) {
    kgs <- na.omit(kgs)
  }

  return(kgs)

}

calculate.curves <- function(list_of_kgs, max, alpha = 1, color = "black") {

  W_0 <- .Machine$double.eps
  A <- max

  curves <- lapply(list_of_kgs, function(x)
    stat_function(data = NULL,
                  fun = function(age) { W_0 * (A / W_0) ^ (1 - exp(-x * age)) },
                  color = color,
                  alpha = alpha))

  return(curves)

}

plot.curves <- function(list_of_curves, max_val,
                        x_limits, x_breaks = waiver(),
                        y_limits) {

  p <- ggplot(NULL) +
    scale_x_continuous(limits = x_limits,
                       breaks = x_breaks,
                       minor_breaks = NULL) +
    scale_y_continuous(limits = y_limits) +
    labs(x = "Age (mo.)", y = "Inventory") +
    geom_hline(yintercept = c(max_val, max_val / exp(1)),
               size = 1,
               linetype = "dashed",
               color = "darkred") +
    theme_bw()

  for (i in list_of_curves) {
    p <- p + i
  }

  return(p)
}

# Solve for a given X or Y, given kg
solve.gomp2 <- function(x = NA, y = NA, k_g, A = 681) {

  W_0 <- .Machine$double.eps

  if (is.na(y) & !is.na(x))
    result <-  W_0 * (A / W_0) ^ (1 - exp(-k_g * x))
  else if (is.na(x) & !is.na(y))
    result <- -log(1 - log(y / W_0) / log(A / W_0)) / k_g
  else
    message("Exactly one of x/y must be specified")

  return(result)

}
