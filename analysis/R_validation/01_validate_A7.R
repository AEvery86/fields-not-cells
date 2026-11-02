# R re-analysis of A7 with reference packages.
# usage: Rscript 01_validate_A7.R [quick|full] [cores]
#   quick: 20 random splits per kidney source and scenario, 5 per intestinal site and scenario, simulation reps 0-9
#   full : every split and all 200 simulation replicates (as in the Python analysis)
# Results of each job are cached in output/cache, so an interrupted run can be restarted and a full run
# reuses the jobs of a quick run.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[1] else "quick"
stopifnot(mode %in% c("quick", "full"))
ncores <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores(logical = FALSE) - 1L)
if (is.na(ncores) || ncores < 1) ncores <- 1L

script_dir <- function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
setwd(script_dir())
dir.create(file.path("output", "cache"), recursive = TRUE, showWarnings = FALSE)
source("A7_R_functions.R")
D <- "data"
rd <- function(f) read.csv(file.path(D, f), check.names = FALSE, stringsAsFactors = FALSE)

cat(sprintf("[%s] mode = %s, cores = %d, R %s\n", format(Sys.time()), mode, ncores, getRversion()))
pk <- c("lme4", "glmmTMB", "MASS", "limma", "speckle")
cat("packages:", paste(sprintf("%s %s", pk, sapply(pk, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "not installed"))), collapse = "; "), "\n")
if (!HAVE_SPECKLE) cat("speckle not installed: propeller uses the limma-based fallback in A7_R_functions.R\n")

# ------------------------------------------------------------------ build job list
jobs <- list()
add <- function(j) jobs[[length(jobs) + 1]] <<- j

for (src in c("GSE325587", "Lupus")) {
  F <- rd(sprintf("kidney_%s_fields.csv", src))
  lab <- rd(sprintf("kidney_%s_labels.csv", src))
  X <- as.matrix(F[, lab$column]); colnames(X) <- lab$label
  lineage <- setNames(lab$lineage, lab$label)
  cats <- rd(sprintf("kidney_%s_categories.csv", src))$category
  sp <- rd(sprintf("kidney_%s_splits.csv", src))
  dis <- as.character(rd(sprintf("kidney_%s_disease_individuals.csv", src))$individual)
  ids <- sort(unique(sp$split))
  if (mode == "quick") { set.seed(1); ids <- sort(sample(ids, 20)) }
  sel <- list(all_fields = F$disease == 1, `3_fields` = F$disease == 1 & F$in_3_fields == 1)
  add(list(id = sprintf("kidney_%s_real", src), type = "kidney", source = src, scenario = "real_contrast", split = -1,
           X = X, ind = F$patient, g1 = dis, cats = cats, lineage = lineage))
  for (sc in names(sel)) for (s in ids) {
    r <- sel[[sc]]
    add(list(id = sprintf("kidney_%s_%s_%d", src, sc, s), type = "kidney", source = src, scenario = sc, split = s,
             X = X[r, , drop = FALSE], ind = F$patient[r], g1 = as.character(sp$member[sp$split == s]),
             cats = cats, lineage = lineage))
  }
}

sites <- rd("intestine_sites.csv")
for (k in seq_len(nrow(sites))) {
  si <- sites$site_index[k]
  F <- rd(sprintf("intestine_site%d_fields.csv", si))
  lab <- rd(sprintf("intestine_site%d_labels.csv", si))
  X <- as.matrix(F[, lab$column]); colnames(X) <- lab$label
  cats <- rd(sprintf("intestine_site%d_categories.csv", si))$category
  sp <- rd(sprintf("intestine_site%d_splits.csv", si))
  ids <- sort(unique(sp$split))
  if (mode == "quick") ids <- head(ids, 5)
  sel <- list(all_fields = rep(TRUE, nrow(F)), `10_fields` = F$in_10_fields == 1)
  for (sc in names(sel)) for (s in ids) {
    r <- sel[[sc]]
    add(list(id = sprintf("intestine_%d_%s_%d", si, sc, s), type = "intestine", site = sites$site[k], scenario = sc, split = s,
             X = X[r, , drop = FALSE], ind = F$donor[r], g1 = as.character(sp$member[sp$split == s]),
             cats = cats, lineage = NULL))
  }
}

S <- read.csv(file.path(D, "simulation_fields.csv.gz"))
reps <- if (mode == "quick") 0:9 else sort(unique(S$rep))
for (r in reps) add(list(id = sprintf("sim_%d", r), type = "sim", rep = r, S = S[S$rep == r, ]))
rm(S)

# ------------------------------------------------------------------ job runner
run_job <- function(j) {
  cf <- file.path("output", "cache", paste0(j$id, ".rds"))
  if (file.exists(cf)) return(readRDS(cf))
  res <- tryCatch(run_job_core(j), error = function(e) { message("job ", j$id, " failed: ", conditionMessage(e)); NULL })
  if (!is.null(res)) saveRDS(res, cf)
  res
}
run_job_core <- function(j) {
  if (j$type == "sim") {
    out <- list()
    for (delta in c(0, 0.4)) for (nper in c(3, 6, 12)) for (k in c(3, 10, 36)) {
      x <- j$S[abs(j$S$delta - delta) < 1e-9 & j$S$individual < nper & j$S$field_order <= k, ]
      C <- as.matrix(x[, paste0("T", 1:6)]); colnames(C) <- paste0("T", 1:6)
      ind <- paste(x$group, x$individual, sep = "_")
      res <- one_contrast(C, ind, unique(ind[x$group == 1]), "T1", NULL)
      res$category <- NULL
      out[[length(out) + 1]] <- cbind(data.frame(delta = delta, n_per_group = nper, k = k, rep = j$rep), res)
    }
    res <- do.call(rbind, out)
  } else {
    res <- one_contrast(j$X, j$ind, j$g1, j$cats, j$lineage)
    if (j$type == "kidney") res <- cbind(data.frame(source = j$source, scenario = j$scenario, split = j$split), res)
    else res <- cbind(data.frame(site = j$site, scenario = j$scenario, split = j$split), res)
  }
  res
}

# ------------------------------------------------------------------ run in batches with progress
todo <- which(!file.exists(file.path("output", "cache", paste0(sapply(jobs, `[[`, "id"), ".rds"))))
cat(sprintf("[%s] %d jobs, %d already cached, %d to run\n", format(Sys.time()), length(jobs), length(jobs) - length(todo), length(todo)))
t0 <- Sys.time()
if (length(todo)) {
  cl <- NULL
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    wd <- getwd()
    parallel::clusterCall(cl, function(w) { setwd(w); source("A7_R_functions.R"); NULL }, wd)
    parallel::clusterExport(cl, c("run_job", "run_job_core"))
  }
  bs <- max(ncores * 4, 8)
  for (b in split(todo, ceiling(seq_along(todo) / bs))) {
    if (is.null(cl)) lapply(jobs[b], run_job) else parallel::parLapplyLB(cl, jobs[b], run_job)
    done <- max(which(todo == max(b)))
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("[%s] %d/%d jobs done, %.1f min elapsed, about %.0f min left\n", format(Sys.time(), "%H:%M:%S"),
                done, length(todo), el, el / done * (length(todo) - done)))
    flush.console()
  }
  if (!is.null(cl)) parallel::stopCluster(cl)
}

# ------------------------------------------------------------------ collect
collect <- function(type) {
  ids <- sapply(Filter(function(j) j$type == type, jobs), `[[`, "id")
  f <- file.path("output", "cache", paste0(ids, ".rds"))
  if (any(!file.exists(f))) cat(sprintf("WARNING: %d %s jobs failed and are missing (rerun to retry)\n", sum(!file.exists(f)), type))
  do.call(rbind, lapply(f[file.exists(f)], readRDS))
}
k <- collect("kidney"); h <- collect("intestine"); s <- collect("sim")
write.csv(k, gzfile(sprintf("output/R_kidney_%s.csv.gz", mode)), row.names = FALSE)
write.csv(h, gzfile(sprintf("output/R_intestine_%s.csv.gz", mode)), row.names = FALSE)
write.csv(s, gzfile(sprintf("output/R_simulation_%s.csv.gz", mode)), row.names = FALSE)
writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
cat(sprintf("[%s] finished; results written to output/\n", format(Sys.time())))

source("02_compare_with_python.R")
compare_with_python(mode)
