
#' @importFrom broom tidy
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
      tidy(lm(form, data=dfx))
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



#' @export
#' @importFrom quantreg rq
#' @importFrom ppcor pcor
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
      res <- lm(y ~ baseline + x2, data=df1)
      coef(res)[2:3]
    } else if (method == "rlm") {
      res <- rlm(y ~ baseline + x2, data=df1, maxit=100)
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

template_sample <- function(source_tab, template, fixgroup="fixgroup", time=NULL, outcol="sample_out") {
  x1 <- rlang::sym(template)
  x2 <- rlang::sym(fixgroup)

  ret <- source_tab %>% select(a=!!x1, b=!!x2) %>% pmap(function(a,b) {
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
