
#' Template Multiple Regression
#'
#' This function performs a multiple regression on the provided source table with the specified
#' response and covariates using the chosen regression method.
#'
#' @param source_tab A data frame containing the source maps.
#' @param response A character string specifying the column name of the response variable.
#' @param covars A character vector specifying the column names of the covariates.
#' @param method A character vector of available regression methods. Default is c("lm", "rlm", "nnls", "logistic").
#'   The selected method will be used for the regression analysis.
#'   - "lm": Linear regression (default).
#'   - "rlm": Robust linear regression.
#'   - "nnls": Non-negative least squares.
#'   - "logistic": Logistic regression.
#' @param intercept A logical value indicating whether to include an intercept term in the model.
#'   Default is TRUE.
#'
#' @return A data frame with the source table augmented with a new column, multireg, containing
#'   the results of the multiple regression analysis.
#' @importFrom dplyr rowwise bind_cols
#' @importFrom broom tidy
#' @importFrom MASS rlm
#' @importFrom nnls nnls
#' @importFrom stats lm glm
#' @export
template_multireg <- function(source_tab, response, covars, method=c("lm", "rlm", "nnls", "logistic"), intercept=TRUE) {
  ret <- source_tab %>% rowwise() %>% do( {
    #browser()
    y <- as.vector(.[[response]]$z/sum(.[[response]]$z))
    xs <- lapply(covars, function(x) as.vector(.[[x]]$z/sum(.[[x]]$z)))
    names(xs) <- covars
    xs[[".response"]] <- y
    dfx <- try(as.data.frame(xs))
    #if (inherits(dfx, "try-error")) {
    #  browser()
    #}

    form <- if (intercept) {
      as.formula(paste(".response", "~", paste(covars, collapse=" + ")))
    } else {
      as.formula(paste(".response", "~", paste(covars, collapse=" + "), " -1"))
    }


   # browser()

    fit <- if (method== "lm") {
      tidy(stats::lm(form, data=dfx))
    } else if (method == "rlm") {
      tidy(MASS::rlm(form, data=dfx))
    } else if (method == "logistic") {
      tidy(glm(form, data=dfx, family="binomial"))
    } else if (method == "nnls") {
      A <- as.matrix(dfx[covars])
      b <- dfx[[".response"]]
      fit <- nnls(A,b)
      data.frame(term=covars, estimate=coef(fit))
    } else {
      stop()
    }

    tibble(multireg=list(fit))
  })

  #ret %>% mutate(multireg_result=ret)
  ret <- bind_cols(source_tab, ret)
  ret

}




#' Template Regression
#'
#' This function performs a template regression on the provided reference and source tables to estimate
#' the beta weights for the baseline and source maps.
#'
#' @param ref_tab A data frame containing the reference maps.
#' @param source_tab A data frame containing the source maps.
#' @param match_on A character string specifying the column name to be used for matching between the
#'   reference and source tables.
#' @param baseline_tab A data frame containing the baseline maps.
#' @param baseline_key A character string specifying the column name to be used for matching between the
#'   baseline table and the source table.
#' @param method A character vector of available regression methods. Default is c("lm", "rlm", "rank").
#'   The selected method will be used for the regression analysis.
#'   - "lm": Linear regression (default).
#'   - "rlm": Robust linear regression.
#'   - "rank": Rank-based correlation.
#'
#' @return A data frame with the source table augmented with two new columns, beta_baseline and beta_source,
#'   representing the estimated beta weights for the baseline and source maps, respectively.
#' @importFrom dplyr ungroup mutate filter rowwise
#' @importFrom quantreg rq
#' @importFrom ppcor pcor
#' @importFrom MASS rlm
#' @export
template_regression <- function(ref_tab, source_tab, match_on,
                                baseline_tab, baseline_key,
                                method=c("lm", "rlm", "rank")) {
  method <- match.arg(method)
  matchind <- match(source_tab[[match_on]], ref_tab[[match_on]])
  source_tab <- source_tab %>% ungroup() %>% mutate(matchind=matchind)

  if (any(is.na(matchind))) {
    warning("did not find matching template map for all source maps. Removing non-matching elements.")
    source_tab <- source_tab %>% filter(!is.na(matchind))
    matchind <- matchind[!is.na(matchind)]
  }


  ret <- source_tab %>% rowwise() %>% do( {
    id <- which(baseline_tab[[baseline_key]] == .[[baseline_key]][1])
    bdens <- baseline_tab$density[[id]]
    d1 <- ref_tab$density[[.$matchind]]
    d2 <- .$density

    df1 <- data.frame(y=as.vector(d2$z), baseline=as.vector(bdens$z), x2=as.vector(d1$z))

    est <- if (method == "lm") {
      res <- stats::lm(y ~ baseline + x2, data=df1)
      coef(res)[2:3]
    } else if (method == "rlm") {
      res <- MASS::rlm(y ~ baseline + x2, data=df1, maxit=100)
      coef(res)[2:3]
    } else if (method == "rank") {
      #browser()
      #res <- rfit(y ~ baseline + x2, data=df1)
      #coef(res)[2:3]
      res <- ppcor::pcor(df1, method="spearman")
      res$estimate[2:3,1]
    } else {
      stop()
    }

    data.frame(b0=est[1], b1=est[2])
  })


  source_tab %>% mutate(beta_baseline=ret$b0, beta_source=ret$b1)

}

#' Sample density maps with coordinates derived from fixation groups.
#'
#' this function extracts the density for any arbitrary time point in a trial. It simply extracts the value of the density map
#' for the fixation at time t. The fixations are taken from the `fixgroup` variable and may be associated, for example,
#' with an independent set of trials.
#'
#' @param source_tab the name of the table containing the density map and fixations
#' @param template the name of the template density variable
#' @param fixgroup the name of the fixation group supplying the spatiotemporal coordinates used to sample the template
#' @param time the time points used to extract coordinates from the the `fixation_group`
#' @param outcol the name of the output variable
#' @export
#' @importFrom purrr pmap
#' @importFrom tibble add_column
#' @importFrom rlang sym
template_sample <- function(source_tab, template, fixgroup="fixgroup", time=NULL, outcol="sample_out") {
  x1 <- rlang::sym(template)
  x2 <- rlang::sym(fixgroup)

  ## filters out NULLs
  ret <- source_tab %>% select(a=!!x1, b=!!x2) %>% filter(!(is.null(a) | is.null(b))) %>% pmap(function(a,b) {
    sample_density(a,b, time)
  })


  # res <- source_tab %>% rowwise() %>% do({
  #   temp <- .[[template]]
  #   fix <- .[[fixgroup]]
  #   cds <- as.matrix(cbind(fix$x, fix$y))
  #   s <- sample_density(temp, fix, time)
  #   browser()
  #   as_tibble(rbind(.)) %>% mutate(samples=list(s))
  # })

  oname <- rlang::sym(outcol)
  source_tab %>% add_column(!!oname := ret)
}
