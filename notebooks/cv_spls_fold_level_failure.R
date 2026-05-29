#!/usr/bin/env Rscript
# Plain-R mirror of cv_spls_fold_level_failure.Rmd (same chunks, in order).
# Run in the BayGMST image (has R + `spls`):
#   docker run --rm -v "${PWD}/notebooks:/nb:ro" --entrypoint Rscript baygmst:local /nb/cv_spls_fold_level_failure.R

library(spls)

cat("## 1. The guard inside spls()\n")
writeLines(grep("zero variance", deparse(body(spls)), value = TRUE))

cat("\n## 2. A reducer-like calibration design (3 dense + 5 gappy proxies)\n")
set.seed(1)
n <- 30
y <- rnorm(n)
dense  <- matrix(rnorm(n * 3), n, 3)
sparse <- matrix(0, n, 5)
for (j in 1:5) sparse[sample(n, 1), j] <- rnorm(1)   # each non-zero in ONE year
X <- cbind(dense, sparse)
colnames(X) <- c(paste0("dense", 1:3), paste0("sparse", 1:5))
print(dim(X))

cat("\n## 3. Global zero-variance filter (current fix): every column passes\n")
gv <- apply(X, 2, var)
print(round(gv, 4))
cat("survive var>0:", sum(gv > 0), "of", ncol(X), "\n")

cat("\n## 4. One CV fold makes a 'good' column constant\n")
omit <- which(sparse[, 1] != 0)
cat("fold omits year:", omit, "\n")
Xtrain <- X[-omit, , drop = FALSE]
print(round(apply(Xtrain, 2, var), 4))

cat("\n## 5. spls() on that fold -> the exact reducer error\n")
print(tryCatch(spls(Xtrain, y[-omit], K = 2, eta = 0.5),
               error = function(e) conditionMessage(e)))

cat("\n## 6. cv.spls() on the globally-filtered matrix still fails\n")
print(tryCatch(cv.spls(X, y, fold = 5, K = 1:2, eta = c(0.3, 0.6), plot.it = FALSE),
               error = function(e) conditionMessage(e)))

cat("\n## 8.1 Pigeonhole boundary: m >= ceil(n/fold)+1 is safe\n")
n <- 30; fold <- 5
hmax <- ceiling(n / fold)
cat(sprintf("n=%d, fold=%d -> hmax = ceil(n/fold) = %d\n", n, fold, hmax))
for (m in c(hmax - 1, hmax, hmax + 1))
  cat(sprintf("  m=%d obs: a single fold can blank the column? %s\n", m, m <= hmax))

cat("\n## 8.4 Scale: density filter keeps dense cols -> cv.spls completes\n")
set.seed(42)
n <- 30; fold <- 5; hmax <- ceiling(n / fold)
p_dense <- 4; p_gappy <- 40
y  <- rnorm(n)
Xd <- matrix(rnorm(n * p_dense), n, p_dense)
Xg <- matrix(0, n, p_gappy)
for (j in 1:p_gappy) Xg[sample(n, sample(1:3, 1)), j] <- rnorm(1)
X <- cbind(Xd, Xg)
cat("raw cv.spls            -> ",
    tryCatch({ cv.spls(X, y, fold = fold, K = 1:3, eta = 0.5, plot.it = FALSE); "completed" },
             error = function(e) paste("ERROR:", conditionMessage(e))), "\n")
keep <- colSums(X != 0) >= hmax + 1
cat(sprintf("density filter (>= %d obs) keeps %d of %d cols\n", hmax + 1, sum(keep), ncol(X)))
cat("density-filtered cv.spls -> ",
    tryCatch({ cv.spls(X[, keep, drop = FALSE], y, fold = fold, K = 1:3, eta = 0.5, plot.it = FALSE); "completed" },
             error = function(e) paste("ERROR:", conditionMessage(e))), "\n")
