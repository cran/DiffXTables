# test_DiffXTables.R
#
# Updated:
#  June 12, 2019. Fixed a test case due to wrong code before.
#  October 29, 2019 Added a test case for 2X2 tables
#  Jan 2020. Some basic tests are added for the differential 
#    simulator function.
#  March 18, 2020. Improved readability of code.
#  March 19, 2020. More tests are added for the simulator function.

library(testthat)
library(DiffXTables)

context("Testing sharma.song.test()")

test_that("Testing sharma.song.test", {
  
  # 2X2 matrix test
  x <- list(
    matrix(c(4,0,
             0,4), nrow=2),
    matrix(c(0,4,
             4,0), nrow=2)
  )
  h <- sharma.song.test(x)
  expect_equivalent(signif(h$p.value, 8), 
                    pchisq(16, 1, lower.tail=FALSE))
  expect_equivalent(signif(h$statistic, 8), 16)
  expect_equivalent(h$parameter, 1)
  
  # 2X2 all zeros 
  x <- list(
    matrix(c(0,0,
             0,0), nrow=2),
    matrix(c(0,0,
             0,0), nrow=2)
  )
  h <- sharma.song.test(x)
  expect_equivalent(signif(h$p.value, 8), 1)
  expect_equivalent(signif(h$statistic, 8), 0)
  expect_equivalent(h$parameter, 0)
  

  # all zero matrices

  x0 <- matrix(c(0,0,0,
                 0,0,0,
                 0,0,0), nrow=3)

  x <- list(x0, x0, x0)

  h <- sharma.song.test(x)
  expect_equivalent(signif(h$p.value, 8), 1)
  expect_equivalent(signif(h$statistic, 8), 0)
  expect_equivalent(h$parameter, 0)

  # conserved patterns
  x0 <- matrix(c(4,0,0,
                 0,4,0,
                 0,0,4), nrow=3)

  x <- list(x0, x0, x0)

  h <- sharma.song.test(x)
  expect_equivalent(signif(h$p.value, 8), 1)
  expect_equivalent(signif(h$statistic, 8), 0)
  expect_equivalent(h$parameter, 8)

  # second order differential
  x <- list(
    matrix(c(4,0,0,
             0,4,0,
             0,0,4), nrow=3),
    matrix(c(0,4,4,
             4,0,4,
             4,4,0), nrow=3)
  )

  h <- sharma.song.test(x)
 
  expect_equivalent(signif(h$p.value, 8),
                      pchisq(36, 4, lower.tail=FALSE))
  expect_equivalent(signif(h$statistic, 8), 36)
  expect_equivalent(h$parameter, 4)
  
  # indepdendent differential patterns
  x <- list(
    matrix(c(16, 4, 20,
              4, 1,  5,
             20, 5, 25), nrow = 3, byrow = TRUE),
    matrix(c(1, 1,  8,
             1, 1,  8,
             8, 8, 64), nrow = 3, byrow = TRUE)
  )

  h <- sharma.song.test(x)
  expect_equivalent(signif(h$p.value, 8), 1)
  expect_equivalent(signif(h$statistic, 8), 0)
  expect_equivalent(h$parameter, 4)

  x0 <- matrix(c(2,  1, 4, 7,
                 5, 10, 8, 2,
                 0,  7, 5, 3), nrow = 3, byrow=T)

  x <- list(x0, x0*2, x0*5)

  h <- sharma.song.test(x)
  expect_equivalent(signif(h$p.value, 8), 1)
  expect_equivalent(signif(h$statistic, 8), 0)
  expect_equivalent(h$parameter, 12)

  x <- list(
    matrix(c(3, 0, 0,
             0, 3, 0,
             0, 0, 3), nrow = 3, byrow=T),
    matrix(c(1, 1, 1,
             1, 1, 1,
             1, 1, 1), nrow = 3, byrow=T)
  )

  h <- sharma.song.test(x)
  expect_equivalent(signif(h$p.value, 8),
                    pchisq(9, 4, lower.tail=FALSE))
  expect_equivalent(signif(h$statistic, 8), 9)
  expect_equivalent(h$parameter, 4)

  x <- list(
    matrix(c(3, 0, 0,
             0, 3, 0,
             0, 0, 3), nrow = 3, byrow=T),
    matrix(c(0, 0, 3,
             0, 3, 0,
             3, 0, 0), nrow = 3, byrow=T)
  )

  h <- sharma.song.test(x)
  expect_equivalent(signif(h$p.value, 8),
                    pchisq(18, 4, lower.tail=FALSE))
  expect_equivalent(signif(h$statistic, 8), 18)
  expect_equivalent(h$parameter, 4)

  x <- list(
    matrix(c(3, 1, 0,
             0, 3, 0,
             4, 0, 3), nrow = 3, byrow=T),
    matrix(c(7, 0, 3,
             0, 3, 5,
             3, 0, 0), nrow = 3, byrow=T)
  )

  h2 <- sharma.song.test(x)
  x.exp <- lapply(x, expected)
  h1 <- cp.chisq.test(x.exp)
  h <- cp.chisq.test(x)

})

context("Testing cp.chisq.test()")

# Testing comparative chi square test
test_that("Testing comparative chi-squared test", {
  
  # 2X2 matrix test
  x <- list(
    matrix(c(4,0,
             0,4), nrow=2),
    matrix(c(0,4,
             4,0), nrow=2)
  )
  h <- cp.chisq.test(x)
  expect_equivalent(signif(h$p.value, 8), 
                    pchisq(16, 1, lower.tail=FALSE))
  expect_equivalent(signif(h$statistic, 8), 16)
  expect_equivalent(h$parameter, 1)
  

  # all zero matrices
  x0 <- matrix(c(0,0,0,
                 0,0,0,
                 0,0,0), nrow=3)

  x <- list(x0, x0, x0)

  h <- cp.chisq.test(x)
  expect_equivalent(signif(h$p.value, 8), 1)
  expect_equivalent(signif(h$statistic, 8), 0)
  expect_equivalent(h$parameter, 0)

  h <- cp.chisq.test(x, method="nchisq")
  expect_equivalent(signif(h$p.value, 8), 1)

  # Conserved pattern
  x0 <- matrix(c(4,0,0,
                 0,4,0,
                 0,0,4), nrow=3)

  x <- list(x0, x0, x0)

  h <- cp.chisq.test(x, method="nchisq")
  expect_equivalent(signif(h$p.value, 8), 0.97724987)

  # Second order differential
  x <- list()

  x[[1]] <- matrix(c(4,0,0,
                     0,4,0,
                     0,0,4), nrow=3)

  x[[2]] <- matrix(c(0,4,4,
                     4,0,4,
                     4,4,0), nrow=3)

  h <- cp.chisq.test(x)
  expect_equivalent(signif(h$p.value, 8), 2.8936962e-07)
  expect_equivalent(signif(h$statistic, 8), 36)
  expect_equivalent(h$parameter, 4)

  h <- cp.chisq.test(x, method="nchisq")
  expect_equivalent(signif(h$p.value, 8), 0)

  # First order differential patterns
  x <- matrix(c(4,0,4,
                0,4,0,
                1,0,1), 3)
  y <- t(x)
  z <- matrix(c(1,0,1,
                4,0,4,
                0,4,0), 3)
  data <- list(x,y,z)
  h <- cp.chisq.test(data)
  expect_equivalent(signif(h$p.value, 8), 1.3542453e-06)
  expect_equivalent(signif(h$statistic, 8), 42)
  expect_equivalent(h$parameter, 8)

  h <- cp.chisq.test(data, method="nchisq")
  expect_equivalent(signif(h$p.value, 8), 9.4795348e-18)

})

context("Testing heterogeneity.test()")
# Testing heterogeneity test
test_that("Testing heterogeneity test", {
  
  # 2X2 matrix test
  x <- list(
    matrix(c(4,0,
             0,4), nrow=2),
    matrix(c(0,4,
             4,0), nrow=2)
  )
  h <- heterogeneity.test(x)
  expect_equivalent(signif(h$p.value, 8), 
                    pchisq(9, 1, lower.tail=FALSE))
  expect_equivalent(signif(h$statistic, 8), 9)
  expect_equivalent(h$parameter, 1)

  # all zero matrices
  x <- list()

  x[[1]] <- matrix(c(0,0,0,
                     0,0,0,
                     0,0,0), nrow=3)
  x[[2]] <- x[[1]]
  x[[3]] <- x[[1]]

  h <- heterogeneity.test(x)
  expect_equivalent(signif(h$p.value, 8), 1)
  expect_equivalent(signif(h$statistic, 8), 0)
  expect_equivalent(h$parameter, 0)


  # Conserved pattern
  x <- list()

  x[[1]] <- matrix(c(4,0,0,
                     0,4,0,
                     0,0,4), nrow=3)
  x[[2]] <- x[[1]]
  x[[3]] <- x[[1]]
  h <- heterogeneity.test(x)
  expect_equivalent(signif(h$p.value, 8), 1)
  expect_equivalent(signif(h$statistic, 8), 0)
  expect_equivalent(h$parameter, 8)


  # second order differential
  x <- list()

  x[[1]] <- matrix(c(4,0,0,
                     0,4,0,
                     0,0,4), nrow=3)

  x[[2]] <- matrix(c(0,4,4,
                     4,0,4,
                     4,4,0), nrow=3)

  h <- heterogeneity.test(x)
  expect_equivalent(signif(h$p.value, 8), 2.8936962e-07)
  expect_equivalent(signif(h$statistic, 8), 36)
  expect_equivalent(h$parameter, 4)

  x <- matrix(c(4,0,4,
                0,4,0,
                1,0,1), nrow=3)
  y <- t(x)
  z <- matrix(c(1,0,1,
                4,0,4,
                0,4,0), nrow=3)
  data <- list(x,y,z)
  h <- heterogeneity.test(data)
  expect_equivalent(signif(h$p.value, 8), 0.0001144345 )
  expect_equivalent(signif(h$statistic, 8), 31.5)
  expect_equivalent(h$parameter, 8)
  
})

context("Testing simulate_diff_tables()")

# Testing simulateDiffXTables 
test_that("Testing simulate_diff_tables method", {
  
  # First-order differential 3 x 4 tables for K = 3 conditions
  x <- simulate_diff_tables(3, 3, 4, 100, 50, type = "first-order")
  expect_equivalent(length(x$probability.tables), 3)
  expect_equivalent(length(x$contingency.tables), 3)
  
  for(i in seq_along(x$contingency.tables))
  {
    expect_equivalent(sum(x$contingency.tables[[i]]), 100)
    expect_equivalent(sum(x$probability.tables[[i]]), 1)
    
    if(i > 1) {
      expect_equivalent(
        sum(x$contingency.tables[[i]] != x$contingency.tables[[i-1]]) > 0,
        TRUE
      )
      expect_equivalent(
        sum(x$probability.tables[[i]] != x$probability.tables[[i-1]]) > 0,
        TRUE
      )
    }
  }
  
  # 3X5 second-order tables for K = 4 conditions
  x <- simulate_diff_tables(4, 3, 5, 140, 70, type = "second-order")
  expect_equivalent(length(x$contingency.tables), 4)
  expect_equivalent(length(x$probability.tables), 4)
  
  for(i in seq_along(x$contingency.tables))
  {
    expect_equivalent(sum(x$contingency.tables[[i]]), 140)
    expect_equivalent(sum(x$probability.tables[[i]]), 1)
    
    if(i > 1) {
      expect_equivalent(
        sum(x$contingency.tables[[i]] != x$contingency.tables[[i-1]]) > 0,
        TRUE
      )
      expect_equivalent(
        sum(x$probability.tables[[i]] != x$probability.tables[[i-1]]) > 0,
        TRUE
      )
    }
  }
  
  # 4X4 second-order tables for K = 2 conditions
  x <- simulate_diff_tables(2, 4, 4, 140, 80, type = "full-order")
  expect_equivalent(length(x$contingency.tables), 2)
  expect_equivalent(length(x$probability.tables), 2)
  
  for(i in seq_along(x$contingency.tables))
  {
    expect_equivalent(sum(x$contingency.tables[[i]]), 140)
    expect_equivalent(sum(x$probability.tables[[i]]), 1)
    
    if(i > 1) {
      expect_equivalent(
        sum(x$contingency.tables[[i]] != x$contingency.tables[[i-1]]) > 0,
        TRUE
      )
      expect_equivalent(
        sum(x$probability.tables[[i]] != x$probability.tables[[i-1]]) > 0,
        TRUE
      )
    }
  }
  
})


