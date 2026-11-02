# Compare the R reference-package results with the Python implementation used in the manuscript.
# Called at the end of 01_validate_A7.R; can also be run alone: Rscript 02_compare_with_python.R [quick|full]

METHODS <- c("glmm_ri_lrt", "glmm_ri_wald", "glmm_olre_lrt", "glmm_olre_wald", "betabin", "nb_glm",
             "propeller_logit", "propeller_asin")
LABEL <- c(glmm_ri_lrt = "GLMM (individual) LRT", glmm_ri_wald = "GLMM (individual) Wald",
           glmm_olre_lrt = "GLMM (individual + field) LRT", glmm_olre_wald = "GLMM (individual + field) Wald",
           betabin = "Beta-binomial (individual)", nb_glm = "Negative binomial (individual)",
           propeller_logit = "propeller (logit)", propeller_asin = "propeller (arcsine)")
RPKG <- c(glmm_ri_lrt = "lme4::glmer, nAGQ = 7", glmm_ri_wald = "lme4::glmer, nAGQ = 7",
          glmm_olre_lrt = "lme4::glmer, Laplace", glmm_olre_wald = "lme4::glmer, Laplace",
          betabin = "glmmTMB betabinomial", nb_glm = "MASS::glm.nb",
          propeller_logit = "speckle::propeller.ttest", propeller_asin = "speckle::propeller.ttest")

compare_with_python <- function(mode = "quick") {
  rdz <- function(f) read.csv(gzfile(f), check.names = FALSE, stringsAsFactors = FALSE)
  pyz <- function(f) read.delim(gzfile(f), check.names = FALSE, stringsAsFactors = FALSE)
  Rk <- rdz(sprintf("output/R_kidney_%s.csv.gz", mode))
  Rh <- rdz(sprintf("output/R_intestine_%s.csv.gz", mode))
  Rs <- rdz(sprintf("output/R_simulation_%s.csv.gz", mode))
  Pk <- pyz("python_results/kidney.tsv.gz")
  Ph <- do.call(rbind, lapply(Sys.glob("python_results/hickey_*.tsv.gz"), pyz))
  Ps <- do.call(rbind, lapply(Sys.glob("python_results/sim_*.tsv.gz"), pyz))

  mk <- function(R, P, keys, dataset) {
    m <- merge(R, P[, c(keys, METHODS)], by = keys, suffixes = c(".R", ".py"))
    m$dataset <- dataset(m)
    m
  }
  K <- mk(Rk, Pk, c("source", "scenario", "split", "category"),
          function(m) ifelse(m$scenario == "real_contrast", paste(m$source, "real contrast"), paste0(m$source, " null, ", m$scenario)))
  H <- mk(Rh, Ph, c("site", "scenario", "split", "category"), function(m) paste0("Intestine donor null, ", m$scenario))
  S <- mk(Rs, Ps, c("delta", "n_per_group", "k", "rep"),
          function(m) paste0("Sim ", ifelse(m$delta == 0, "null", "power"), ", ", m$n_per_group, "/group, ", m$k, " fields"))
  cat(sprintf("\nmatched rows: kidney %d of %d, intestine %d of %d, simulation %d of %d\n",
              nrow(K), nrow(Rk), nrow(H), nrow(Rh), nrow(S), nrow(Rs)))

  agree <- function(d, by) {
    out <- list()
    for (g in split(d, d[[by]], drop = TRUE)) for (mth in METHODS) {
      pr <- g[[paste0(mth, ".R")]]; pp <- g[[paste0(mth, ".py")]]
      ok <- is.finite(pr) & is.finite(pp)
      dl <- abs(log10(pmax(pr[ok], 1e-300)) - log10(pmax(pp[ok], 1e-300)))
      out[[length(out) + 1]] <- data.frame(
        group = g[[by]][1], method = mth, test = LABEL[[mth]], R_function = RPKG[[mth]],
        n = sum(ok), R_missing = sum(!is.finite(pr)), python_missing = sum(!is.finite(pp)),
        python_reject_pct = round(100 * mean(pp[ok] < 0.05), 1), R_reject_pct = round(100 * mean(pr[ok] < 0.05), 1),
        decision_agreement_pct = round(100 * mean((pp[ok] < 0.05) == (pr[ok] < 0.05)), 1),
        spearman = round(tryCatch(suppressWarnings(cor(pp[ok], pr[ok], method = "spearman")), error = function(e) NA_real_), 4),
        median_abs_dlog10p = signif(median(dl), 3), p95_abs_dlog10p = signif(unname(quantile(dl, 0.95)), 3),
        stringsAsFactors = FALSE)
    }
    do.call(rbind, out)
  }
  null_rows <- rbind(K[K$scenario != "real_contrast", c("dataset", paste0(METHODS, ".R"), paste0(METHODS, ".py"))],
                     H[, c("dataset", paste0(METHODS, ".R"), paste0(METHODS, ".py"))],
                     S[, c("dataset", paste0(METHODS, ".R"), paste0(METHODS, ".py"))])
  by_set <- agree(null_rows, "dataset")
  null_rows$all <- "all comparisons"
  overall <- agree(null_rows, "all")
  write.csv(by_set, sprintf("output/R_vs_Python_by_dataset_%s.csv", mode), row.names = FALSE)
  write.csv(overall, sprintf("output/R_vs_Python_overall_%s.csv", mode), row.names = FALSE)

  # real contrasts: BH discoveries (q < 0.05), as in Table S9
  real <- K[K$scenario == "real_contrast", ]
  disc <- do.call(rbind, lapply(split(real, real$source), function(g) data.frame(
    source = g$source[1], method = METHODS, test = unname(LABEL[METHODS]),
    python_BH_discoveries = sapply(METHODS, function(m) sum(p.adjust(ifelse(is.finite(g[[paste0(m, ".py")]]), g[[paste0(m, ".py")]], 1), "BH") < 0.05)),
    R_BH_discoveries = sapply(METHODS, function(m) sum(p.adjust(ifelse(is.finite(g[[paste0(m, ".R")]]), g[[paste0(m, ".R")]], 1), "BH") < 0.05)),
    stringsAsFactors = FALSE)))
  write.csv(disc, sprintf("output/R_vs_Python_real_contrasts_%s.csv", mode), row.names = FALSE)

  # robust propeller (speckle default) and warning frequencies, R only
  allR <- rbind(cbind(dataset = "kidney", Rk[, grep("^(propeller|warn_)", names(Rk))]),
                cbind(dataset = "intestine", Rh[, grep("^(propeller|warn_)", names(Rh))]),
                cbind(dataset = "simulation", Rs[, grep("^(propeller|warn_)", names(Rs))]))
  extra <- aggregate(. ~ dataset, data = transform(allR,
             propeller_logit_robust = as.numeric(propeller_logit_robust < 0.05),
             propeller_asin_robust = as.numeric(propeller_asin_robust < 0.05),
             propeller_logit = as.numeric(propeller_logit < 0.05),
             propeller_asin = as.numeric(propeller_asin < 0.05)), FUN = function(x) round(100 * mean(x, na.rm = TRUE), 1), na.action = na.pass)
  write.csv(extra, sprintf("output/R_only_robust_propeller_and_warnings_%s.csv", mode), row.names = FALSE)

  sink(sprintf("output/R_vs_Python_report_%s.txt", mode), split = TRUE)
  cat("R re-analysis of the A7 comparison (mode:", mode, ")\n")
  cat("Rows matched: kidney", nrow(K), "intestine", nrow(H), "simulation", nrow(S), "\n\n")
  cat("Overall agreement across all null and power comparisons\n")
  print(overall[, c("test", "R_function", "n", "R_missing", "python_reject_pct", "R_reject_pct",
                    "decision_agreement_pct", "spearman", "median_abs_dlog10p", "p95_abs_dlog10p")], row.names = FALSE)
  cat("\nRejection at P < 0.05 by dataset (python / R, %)\n")
  wide <- reshape(by_set[, c("group", "test", "python_reject_pct", "R_reject_pct")], idvar = "test", timevar = "group", direction = "wide")
  print(wide, row.names = FALSE)
  cat("\nReal disease-versus-control contrasts: BH discoveries (python / R)\n")
  print(disc[, c("source", "test", "python_BH_discoveries", "R_BH_discoveries")], row.names = FALSE)
  cat("\nR only: rejection (%) of propeller with robust = TRUE (speckle default) versus robust = FALSE, and share of fits with warnings (%)\n")
  print(extra, row.names = FALSE)
  sink()
  invisible(list(overall = overall, by_dataset = by_set, real = disc))
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  compare_with_python(if (length(a)) a[1] else "quick")
}
