
template_multireg <- function(source_tab, response, covars, method=c("lm", "rlm", "nnls", "rq"), intercept=TRUE) {
  ret <- source_tab %>% rowwise() %>% do( {
    y <- as.vector(.[[response]]$z)
    xs <- lapply(covars, function(x) as.vector(.[[x]]$z))
    names(xs) <- covars
    xs[[".response"]] <- y
    dfx <- as.data.frame(xs)

    form <- as.formula(paste(".response", "~", paste(covars, collapse=" + ")))
    lm.1 <- lm(form, data=dfx)
    tibble(multireg=list(tidy(lm.1)))
  })

  ret %>% mutate(multireg_result=ret)
  ret <- bind_cols(source_tab, ret)

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
