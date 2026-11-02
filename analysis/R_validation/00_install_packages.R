# Install the reference packages for the A7 re-analysis (skips packages that are already installed).
# CRAN: lme4, glmmTMB, MASS, BiocManager; Bioconductor: limma, speckle.
# A mirror in China is tried first, then the official repositories.

lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(lib)) { dir.create(lib, recursive = TRUE, showWarnings = FALSE); .libPaths(c(lib, .libPaths())) }

cran_mirrors <- c("https://mirrors.tuna.tsinghua.edu.cn/CRAN/", "https://mirrors.ustc.edu.cn/CRAN/", "https://cloud.r-project.org/")
bioc_mirrors <- c("https://mirrors.tuna.tsinghua.edu.cn/bioconductor", "https://mirrors.ustc.edu.cn/bioc", "https://bioconductor.org")

have <- function(p) requireNamespace(p, quietly = TRUE)

install_cran <- function(pkgs) {
  for (m in cran_mirrors) {
    need <- pkgs[!sapply(pkgs, have)]
    if (!length(need)) return(TRUE)
    cat("Installing from", m, ":", paste(need, collapse = ", "), "\n")
    try(install.packages(need, repos = m, lib = .libPaths()[1], dependencies = c("Depends", "Imports", "LinkingTo")), silent = FALSE)
  }
  all(sapply(pkgs, have))
}

install_bioc <- function(pkgs) {
  for (m in bioc_mirrors) {
    need <- pkgs[!sapply(pkgs, have)]
    if (!length(need)) return(TRUE)
    options(BioC_mirror = m, repos = c(CRAN = cran_mirrors[1]))
    cat("Installing from Bioconductor via", m, ":", paste(need, collapse = ", "), "\n")
    try(BiocManager::install(need, update = FALSE, ask = FALSE, lib = .libPaths()[1]), silent = FALSE)
  }
  all(sapply(pkgs, have))
}

ok1 <- install_cran(c("BiocManager", "lme4", "glmmTMB", "MASS"))
ok2 <- install_bioc(c("limma", "speckle"))
if (!have("speckle") && have("limma"))
  cat("\nNOTE: speckle could not be installed; propeller will be run with a copy of speckle's two-group",
      "code path (proportion transform + limma lmFit/contrasts.fit/eBayes) on top of the official limma package.\n")

for (p in c("lme4", "glmmTMB", "MASS", "limma", "speckle"))
  cat(sprintf("%-8s %s\n", p, if (have(p)) as.character(packageVersion(p)) else "NOT INSTALLED"))
if (!ok1 || !have("limma")) { cat("\nRequired packages are missing; see the messages above.\n"); quit(status = 1) }
cat("\nAll required packages are available.\n")
