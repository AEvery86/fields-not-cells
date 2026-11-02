# Functions for the R re-analysis of A7 (mixed models and differential-abundance tools).
# Reference implementations: lme4::glmer, glmmTMB::glmmTMB (beta-binomial), MASS::glm.nb,
# speckle::getTransformedProps + speckle::propeller.ttest (limma).

suppressPackageStartupMessages({
  library(lme4)
  library(glmmTMB)
  library(MASS)
  library(limma)
})
HAVE_SPECKLE <- requireNamespace("speckle", quietly = TRUE)
if (HAVE_SPECKLE) {
  getTransformedProps <- speckle::getTransformedProps
  propeller.ttest <- speckle::propeller.ttest
} else {
  # Fallback used only if speckle cannot be installed: speckle's two-group path on the official limma functions.
  getTransformedProps <- function(clusters, sample, transform = "logit") {
    tab <- table(sample, clusters)
    props <- tab / rowSums(tab)
    if (transform == "asin") {
      prop.trans <- asin(sqrt(props))
    } else {
      props.pseudo <- (tab + 0.5) / rowSums(tab + 0.5)
      prop.trans <- log(props.pseudo / (1 - props.pseudo))
    }
    list(Counts = t(tab), TransformedProps = t(prop.trans), Proportions = t(props))
  }
  propeller.ttest <- function(prop.list, design, contrasts, robust = TRUE, trend = FALSE, sort = TRUE) {
    fit <- limma::lmFit(prop.list$TransformedProps, design)
    fit.cont <- limma::eBayes(limma::contrasts.fit(fit, contrasts = contrasts), robust = robust, trend = trend)
    data.frame(Tstatistic = fit.cont$t[, 1], P.Value = fit.cont$p.value[, 1],
               FDR = p.adjust(fit.cont$p.value[, 1], method = "BH"), row.names = rownames(prop.list$TransformedProps))
  }
}

# Evaluate an expression, returning NULL on error and recording whether any warning was raised.
quiet <- function(expr) {
  warned <- FALSE
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) NULL),
    warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") },
    message = function(m) invokeRestart("muffleMessage"))
  list(value = val, warned = warned)
}

lrt_p <- function(m1, m0) {
  if (is.null(m1) || is.null(m0)) return(NA_real_)
  l1 <- tryCatch(as.numeric(logLik(m1)), error = function(e) NA_real_)
  l0 <- tryCatch(as.numeric(logLik(m0)), error = function(e) NA_real_)
  if (!is.finite(l1) || !is.finite(l0)) return(NA_real_)
  pchisq(max(2 * (l1 - l0), 0), df = 1, lower.tail = FALSE)
}

wald_p <- function(m, term = "grp") {
  if (is.null(m)) return(NA_real_)
  cf <- tryCatch(coef(summary(m)), error = function(e) NULL)
  if (is.null(cf)) return(NA_real_)
  if (is.list(cf) && !is.matrix(cf)) cf <- cf$cond          # glmmTMB
  if (!term %in% rownames(cf)) return(NA_real_)
  as.numeric(cf[term, ncol(cf)])
}

# Binomial GLMMs on field counts.
#   glmm_ri  : individual random intercept, adaptive Gauss-Hermite quadrature with 7 nodes
#   glmm_olre: individual + field (observation-level) random intercepts, Laplace approximation
glmm_tests <- function(y, n, ind, grp) {
  d <- data.frame(y = y, f = n - y, ind = factor(ind), grp = grp, obs = factor(seq_along(y)))
  ctl <- glmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
  r1 <- quiet(glmer(cbind(y, f) ~ grp + (1 | ind), data = d, family = binomial, nAGQ = 7, control = ctl))
  r0 <- quiet(glmer(cbind(y, f) ~ 1 + (1 | ind), data = d, family = binomial, nAGQ = 7, control = ctl))
  o1 <- quiet(glmer(cbind(y, f) ~ grp + (1 | ind) + (1 | obs), data = d, family = binomial, nAGQ = 1, control = ctl))
  o0 <- quiet(glmer(cbind(y, f) ~ 1 + (1 | ind) + (1 | obs), data = d, family = binomial, nAGQ = 1, control = ctl))
  c(glmm_ri_lrt = lrt_p(r1$value, r0$value), glmm_ri_wald = wald_p(r1$value),
    glmm_olre_lrt = lrt_p(o1$value, o0$value), glmm_olre_wald = wald_p(o1$value),
    warn_glmm_ri = as.numeric(r1$warned || r0$warned), warn_glmm_olre = as.numeric(o1$warned || o0$warned))
}

# Beta-binomial regression on individual totals (glmmTMB), likelihood-ratio test.
bb_test <- function(yi, ni, gi) {
  d <- data.frame(y = yi, f = ni - yi, g = gi)
  m1 <- quiet(glmmTMB(cbind(y, f) ~ g, family = betabinomial(link = "logit"), data = d))
  m0 <- quiet(glmmTMB(cbind(y, f) ~ 1, family = betabinomial(link = "logit"), data = d))
  c(betabin = lrt_p(m1$value, m0$value), warn_betabin = as.numeric(m1$warned || m0$warned))
}

# Negative-binomial GLM on individual totals with log(total) offset (MASS::glm.nb), likelihood-ratio test.
nb_test <- function(yi, ni, gi) {
  d <- data.frame(y = yi, g = gi, off = log(ni))
  m1 <- quiet(glm.nb(y ~ g + offset(off), data = d))
  m0 <- quiet(glm.nb(y ~ 1 + offset(off), data = d))
  c(nb_glm = lrt_p(m1$value, m0$value), warn_nb_glm = as.numeric(m1$warned || m0$warned))
}

# propeller (speckle) on an individuals x clusters count matrix P; gi = 0/1 per row of P.
# Cell-level vectors are rebuilt from the counts so that the official getTransformedProps() is used.
propeller_speckle <- function(P, gi, transform, robust) {
  cnt <- as.vector(t(P))
  sample <- factor(rep(rep(rownames(P), each = ncol(P)), times = cnt), levels = rownames(P))
  clusters <- factor(rep(rep(colnames(P), times = nrow(P)), times = cnt), levels = colnames(P))
  pl <- suppressMessages(getTransformedProps(clusters = clusters, sample = sample, transform = transform))
  for (nm in intersect(c("Counts", "TransformedProps", "Proportions"), names(pl))) {
    m <- unclass(as.matrix(pl[[nm]]))
    if (!all(colnames(P) %in% rownames(m))) m <- t(m)          # clusters x samples
    pl[[nm]] <- m[colnames(P), rownames(P), drop = FALSE]
  }
  grp <- factor(gi, levels = c(0, 1))
  design <- model.matrix(~ 0 + grp)
  res <- tryCatch(suppressMessages(suppressWarnings(
           propeller.ttest(prop.list = pl, design = design, contrasts = c(-1, 1),
                           robust = robust, trend = FALSE, sort = FALSE))),
         error = function(e) NULL)
  if (is.null(res)) return(setNames(rep(NA_real_, ncol(P)), colnames(P)))
  setNames(res[colnames(P), "P.Value"], colnames(P))
}

propeller_all <- function(P, gi, lineage = NULL) {
  out <- list()
  for (tr in c("logit", "asin")) for (rb in c(FALSE, TRUE)) {
    key <- paste0("propeller_", tr, if (rb) "_robust" else "")
    pv <- propeller_speckle(P, gi, tr, rb)
    if (!is.null(lineage)) {
      L <- t(rowsum(t(P), lineage[colnames(P)]))
      pvl <- propeller_speckle(L, gi, tr, rb)
      names(pvl) <- paste("[lineage]", names(pvl))
      pv <- c(pv, pvl)
    }
    out[[key]] <- pv
  }
  out
}

# One comparison: F = field table (data.frame), X = field x label count matrix (colnames = labels),
# ind = individual per field, g1 = individuals in group 1, cats = tested categories.
one_contrast <- function(X, ind, g1, cats, lineage = NULL) {
  ind <- as.character(ind); grp <- as.integer(ind %in% g1)
  n <- rowSums(X)
  P <- rowsum(X, ind)
  gi <- as.integer(rownames(P) %in% g1)
  prop <- propeller_all(P, gi, lineage)
  rows <- lapply(cats, function(cc) {
    if (startsWith(cc, "[lineage] ")) {
      keep <- lineage[colnames(X)] == sub("^\\[lineage\\] ", "", cc)
      y <- rowSums(X[, keep, drop = FALSE])
    } else {
      y <- X[, cc]
    }
    yi <- as.numeric(rowsum(y, ind)[, 1]); ni <- as.numeric(rowsum(n, ind)[, 1])
    gii <- as.integer(rownames(rowsum(y, ind)) %in% g1)
    v <- c(glmm_tests(y, n, ind, grp), bb_test(yi, ni, gii), nb_test(yi, ni, gii),
           sapply(prop, function(p) unname(p[cc])))
    data.frame(category = cc, as.list(v), check.names = FALSE, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
