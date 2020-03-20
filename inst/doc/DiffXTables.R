## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=F, warning=F, results = "hide"----------------------------
require(FunChisq)
require(DiffXTables)

## ---- out.width='60%',fig.asp=0.6---------------------------------------------
tables <- list(
 matrix(c(
   14,  0,  4,
    0,  8,  0,
    4,  0, 12), byrow=TRUE, nrow=3),
 matrix(c(
    7,  0,  2,
    0,  4,  0,
    2,  0,  6), byrow=TRUE, nrow=3)
)
par(mfrow=c(1,2), cex=0.5, oma=c(0,0,2,0))
plot_table(tables[[1]], highlight="none", xlab=NA, ylab=NA)
plot_table(tables[[2]], highlight="none", xlab=NA, ylab=NA)
mtext("Conserved patterns", outer = TRUE)

cp.chisq.test(tables)

sharma.song.test(tables)

heterogeneity.test(tables)

## ---- out.width='60%',fig.asp=0.6---------------------------------------------
tables <- list(
  matrix(c(
    16, 4, 20,
     4, 1,  5,
    20, 5, 25), nrow = 3, byrow = TRUE),
  matrix(c(
     1, 1,  8,
     1, 1,  8,
     8, 8, 64), nrow = 3, byrow = TRUE)
  )
par(mfrow=c(1,2), cex=0.5, oma=c(0,0,2,0))
plot_table(tables[[1]], highlight="none", col="cornflowerblue", xlab=NA, ylab=NA)
plot_table(tables[[2]], highlight="none", col="cornflowerblue", xlab=NA, ylab=NA)
mtext("First-order differential patterns", outer = TRUE)

cp.chisq.test(tables)

sharma.song.test(tables)

heterogeneity.test(tables)

## ---- out.width='60%',fig.asp=0.6---------------------------------------------
tables <- list(
  matrix(c(
    8,  1, 1, 38, 4,
    5,  1, 1, 17, 1,
    2,  1, 1,  9, 1,
    2,  1, 1,  4, 1), nrow=4, byrow = TRUE),
  matrix(c(
    1,  2, 1,  1, 2,
    2,  9, 1,  1, 4,
    2, 13, 1,  1, 1,
    3, 45, 2,  1, 7), nrow=4, byrow = TRUE)
)
par(mfrow=c(1,2), cex=0.5, oma=c(0,0,2,0))
plot_table(tables[[1]], highlight="none", col="cornflowerblue", xlab=NA, ylab=NA)
plot_table(tables[[2]], highlight="none", col="cornflowerblue", xlab=NA, ylab=NA)
mtext("First-order differential patterns", outer = TRUE)

cp.chisq.test(tables)

sharma.song.test(tables)

heterogeneity.test(tables)

## ---- out.width='60%',fig.asp=0.6---------------------------------------------
tables <- list(
  matrix(c(
    4, 0, 0,
    0, 4, 0,
    0, 0, 4
  ), byrow=TRUE, nrow=3),
  matrix(c(
    0, 4, 4,
    4, 0, 4,
    4, 4, 0
  ), byrow=TRUE, nrow=3)
)
par(mfrow=c(1,2), cex=0.5, oma=c(0,0,2,0))
plot_table(tables[[1]], highlight="none", col="salmon", xlab=NA, ylab=NA)
plot_table(tables[[2]], highlight="none", col="salmon", xlab=NA, ylab=NA)
mtext("Second-order differential patterns", outer = TRUE)
cp.chisq.test(tables)

sharma.song.test(tables)

heterogeneity.test(tables)

## ---- out.width='60%',fig.asp=0.6---------------------------------------------
tables <- list(
  matrix(c(
    50,  0, 0,  0,
     0,  0, 1,  0,
     0, 50, 0,  0,
     1,  0, 0,  0,
     0,  0, 0, 50
  ), byrow=T, nrow = 5),
  matrix(c(
     1,  0,  0, 0,
     0,  0, 50, 0,
     0,  1,  0, 0,
    50,  0,  0, 0,
     0,  0,  0, 1
  ), byrow=T, nrow = 5)
)
par(mfrow=c(1,2), cex=0.5, oma=c(0,0,2,0))
plot_table(tables[[1]], highlight="none", col="orange", xlab=NA, ylab=NA)
plot_table(tables[[2]], highlight="none", col="orange", xlab=NA, ylab=NA)
mtext("Differential patterns", outer = TRUE)
cp.chisq.test(tables)

sharma.song.test(tables)

heterogeneity.test(tables)

