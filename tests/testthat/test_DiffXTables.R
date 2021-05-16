# test_DiffXTables.R
#
# Updated:
#  June 12, 2019. Fixed a test case due to wrong code before.
#  October 29, 2019 Added a test case for 2X2 tables
#  Jan 2020. Some basic tests are added for the differential 
#    simulator function.
#  March 18, 2020. Improved readability of code.
#  March 19, 2020. More tests are added for the simulator function.
#  July 17, 2020. Added test cases for marginal.change.test
#  July 20, 2020. Added test cases for strength.test
#  July 21, 2020. Added test cases for type.analysis method.
#  December 17, 2020, Added test cases for uniform null table marginal.
#  April 20, 2021, Added test cases for compensated sharma.song test.

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

  # First row and first column are zero
  x <- list(
    matrix(c(0, 0, 0,
             0, 3, 0,
             4, 0, 3), nrow = 3, byrow=T),
    matrix(c(0, 0, 3,
             0, 3, 5,
             0, 0, 5), nrow = 3, byrow=T)
  )
  
  h2 <- sharma.song.test(x)
  x.exp <- lapply(x, expected)
  h1 <- cp.chisq.test(x.exp)
  h <- cp.chisq.test(x)
  
  ## Only one non-zero value in the column, using compensated 
  ## parameter
  
  x <- list(
    matrix(c(30,20,0,
             10,5,0,
             0,0,0),  ncol = 3),
    matrix(c(30,20,0,
             10,5,0,
             0,0,1), ncol = 3)
  )
  h2 <- sharma.song.test(x, compensated = TRUE)
  expect_equivalent(signif(h2$p.value, 8),
                    pchisq(8.320419, 4, lower.tail=FALSE))
  expect_equivalent(signif(h2$statistic, 8), 8.32, tolerance = 0.0362 )
  expect_equivalent(h2$parameter, 4)
  
  
  x <- list(
    matrix(data = c(410,39,0,
                    46, 5, 0,
                    0, 0, 0), nrow = 3, ncol = 3), 
    matrix(data = c(430,21,0,
                    35,11,0,
                    1, 0, 1), nrow = 3, ncol = 3)
  )
  h2 <- sharma.song.test(x, compensated = TRUE)
  expect_equivalent(signif(h2$p.value, 8),
                    pchisq(33.518, 4, lower.tail=FALSE))
  expect_equivalent(signif(h2$statistic, 8), 33.1, tolerance = 0.0362 )
  expect_equivalent(h2$parameter, 4)
  
  
  x <- list(
    matrix(data = c(434,24,1,
                    14,26, 0,
                    0, 1, 0), nrow = 3, ncol = 3), 
    matrix(data = c(434,24,0,
                    14,26,0,
                    0, 0, 0), nrow = 3, ncol = 3)
  )
  h2 <- sharma.song.test(x, compensated = TRUE)
  expect_equivalent(signif(h2$p.value, 8),
                    pchisq(18.6914, 4, lower.tail=FALSE))
  expect_equivalent(signif(h2$statistic, 8), 18.7, tolerance = 0.0496 )
  expect_equivalent(h2$parameter, 4)
  
  #### Testing sharma.song.test with uniform marginal in null population
  
  # 2X2 matrix test
  x <- list(
    matrix(c(4,0,
             0,4), nrow=2),
    matrix(c(0,4,
             4,0), nrow=2)
  )
  h <- sharma.song.test(x, null.table.marginal = "uniform")
  expect_equivalent(signif(h$p.value, 8), 
                    pchisq(16, 1, lower.tail=FALSE))
  expect_equivalent(signif(h$statistic, 8), 16)
  expect_equivalent(h$parameter, 1)
  
  x <- list(
    matrix(c(0,0,
             0,0), nrow=2),
    matrix(c(0,0,
             0,0), nrow=2)
  )
  h <- sharma.song.test(x, null.table.marginal = "uniform")
  expect_equivalent(signif(h$p.value, 8), 1)
  expect_equivalent(signif(h$statistic, 8), 0)
  expect_equivalent(h$parameter, 0)
  
  
  # second-order differential patterns
  x <- list(
    matrix(c(4,0,0,
             0,4,0,
             0,0,4), nrow=3),
    matrix(c(0,4,4,
             4,0,4,
             4,4,0), nrow=3)
  )
  
  h <- sharma.song.test(x, null.table.marginal = "uniform")
  
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
  
  h <- sharma.song.test(x, null.table.marginal = "uniform")
  expect_equivalent(signif(h$p.value, 8), 1)
  expect_equivalent(signif(h$statistic, 8), 0)
  expect_equivalent(h$parameter, 4)
  
  ## Only one non-zero value in the column, using uniform null distribution 
  ## to demote sparse patterns.
  x <- list(
    matrix(c(411,37,0,
             47,5,0,
             0,0,0), byrow = 3, ncol = 3),
    matrix(c(431,21,0,
             35,11,0,
             1,0,1), byrow = 3, ncol = 3)
  )
  h <- sharma.song.test(x, null.table.marginal = "uniform")
  expect_equivalent(signif(h$p.value, 8), 
                    pchisq(2.0874, 4, lower.tail=FALSE), tolerance = 0.0005)
  expect_equivalent(signif(h$statistic, 8), 2.0874, tolerance = 0.0496 )
  expect_equivalent(h$parameter, 4)
  
  ## Only one non-zero value in the row, using uniform null distribution 
  ## to demote these kind of patterns
  x <- list(
    matrix(data = c(434,24,1,
                    14,26, 0,
                    0, 1, 0), nrow = 3, ncol = 3), 
    matrix(data = c(434,24,0,
                    14,26,0,
                    0, 0, 0), nrow = 3, ncol = 3)
  )
  h2 <- sharma.song.test(x, null.table.marginal = "uniform")
  expect_equivalent(signif(h2$p.value, 8),
                    pchisq(0.0242, 4, lower.tail=FALSE), tolerance = 0.0005)
  expect_equivalent(signif(h2$statistic, 8), 0.0242, tolerance = 0.0496)
  expect_equivalent(h2$parameter, 4)
  
  
  x <- list(
    matrix(c(
      55,  1,  0,
      0,  110,  0,
      0,  0, 16), byrow=TRUE, nrow=3),
    matrix(c(
      13,  0,  0,
      0,  83,  0,
      1,  0,  84), byrow=TRUE, nrow=3)
  )
  h2 <- sharma.song.test(x, null.table.marginal = "uniform")
  expect_equivalent(signif(h2$p.value, 8),
                    pchisq(117.43, 4, lower.tail=FALSE))
  expect_equivalent(signif(h2$statistic, 8), 117.43, tolerance = 0.0496 )
  expect_equivalent(h2$parameter, 4)
  
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
                    pchisq(16, 1, lower.tail=FALSE))
  expect_equivalent(signif(h$statistic, 8), 16)
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

context("Testing marginal.change.test()")

# Testing marginal.change.test
test_that("Testing marginal.change.test method", {
  
  # Same marginal 
  x <- list(
    matrix(c(4,0,
             0,4), nrow=2),
    matrix(c(0,4,
             4,0), nrow=2)
  )
  h <- marginal.change.test(x)
  expect_equivalent(signif(h$p.value, 8), 1)
  expect_equivalent(signif(h$statistic, 8), 0)
  expect_equivalent(h$parameter, 2)
  
  # 2X2 all zeros 
  x <- list(
    matrix(c(0,0,
             0,0), nrow=2),
    matrix(c(0,0,
             0,0), nrow=2)
  )
  h <- marginal.change.test(x)
  expect_equivalent(signif(h$p.value, 8), 1)
  expect_equivalent(signif(h$statistic, 8), 0)
  expect_equivalent(h$parameter, 0)
  
  ## Conserved pattern
  x0 <- matrix(c(4,0,0,
                 0,4,0,
                 0,0,4), nrow=3)
  
  x <- list(x0, x0, x0)
  
  h <- marginal.change.test(x)
  expect_equivalent(signif(h$p.value, 8), 1)
  expect_equivalent(signif(h$statistic, 8), 0)
  expect_equivalent(h$parameter, 8)
  
  ## Marginal table example
  x0 <- matrix(c(8,0,0,
                 0,10,0,
                 0,0,20), nrow=3)
  
  x1 <- matrix(c(8,0,0,
                 0,20,0,
                 0,0,10), nrow=3)
  
  x <- list(x0, x1, x0)
  h <- marginal.change.test(x)
  expect_equivalent(signif(h$p.value, 8), 0.02122649  )
  expect_equivalent(signif(h$statistic, 8),  18 )
  expect_equivalent(h$parameter, 8)
  
  ## Marginal table example
  x0 <- matrix(c(0,30,0,
                 0,10,0,
                 0,0,0), nrow=3)
  
  x1 <- matrix(c(0,0,20,
                 0,50,0,
                 10,0,0), nrow=3)
  
  x <- list(x0, x1)
  h <- marginal.change.test(x)
  expect_equivalent(signif(h$p.value, 8), 7.424e-10)
  expect_equivalent(signif(h$statistic, 8), 48.5 )
  expect_equivalent(h$parameter, 4)
  
})

context("Testing strength.test()")

# Testing strength.test
test_that("Testing strength.test method", {

   x <- list(
     matrix(c(30,0,0,
              0,10,0,
              0,0,20), nrow=3),
     matrix(c(10,0,0,
               0,20,0,
               0,0,30), nrow=3)
  )
  h <- strength.test(x)
  expect_equivalent(signif(h$p.value, 8), 2.264417e-47)
  expect_equivalent(signif(h$statistic, 8), 240)
  expect_equivalent(h$parameter, 8)

  # One table has strong association:
  x <- list(
      matrix(c(4,0,0,
               0,4,0,
               0,0,4), nrow=3),
      matrix(c(4,0,4,
               8,4,8,
               4,0,4), nrow=3)
  )
  h <- strength.test(x)
  expect_equivalent(signif(h$p.value, 8), 0.0005566)
  expect_equivalent(signif(h$statistic, 8), 27.6)
  expect_equivalent(h$parameter, 8)

  # Both tables has no association:
  x <- list(
      matrix(c(4,0,4,
               8,4,8,
               4,0,4), nrow=3),
      matrix(c(4,0,4,
               8,4,8,
               4,0,4), nrow=3)
  )
  h <- strength.test(x)
  expect_equivalent(signif(h$p.value, 8), 0.5152161)
  expect_equivalent(signif(h$statistic, 8), 7.2)
  expect_equivalent(h$parameter, 8)
})




context("Testing type.analysis()")

# Testing association.order.ana method 
test_that("Testing type.analysis method", {
  
  # Type-2 example 
  x <- list(
    matrix(c(20,0,
             0,20), nrow=2),
    matrix(c(0,20,
             20,0), nrow=2)
  )
  invisible(h <- type.analysis(x, 0.05))
  expect_equivalent(h$Type, 2)
  
  # Type Null example 
  x <- list(
    matrix(c(0,0,
             0,0), nrow=2),
    matrix(c(0,0,
             0,0), nrow=2)
  )
  h <- type.analysis(x, 0.05)
  expect_equivalent(h$Type, NULL)
  
  # Type Null example 
  x <- list(
    matrix(c(2,4,
             4,8), nrow=2),
    matrix(c(2,4,
             4,8), nrow=2)
  )
  h <- type.analysis(x, 0.05)
  expect_equivalent(h$Type, NULL)
  
  ## Type-0 example
  x0 <- matrix(c(4,0,0,
                 0,4,0,
                 0,0,4), nrow=3)
  
  x <- list(x0, x0, x0)
  
  h <- type.analysis(x, 0.05)
  expect_equivalent(h$Type, 0)
  
  ## Type-1 example
  x0 <- matrix(c(8,0,0,
                 0,10,0,
                 0,0,20), nrow=3)
  
  x1 <- matrix(c(8,0,0,
                 0,20,0,
                 0,0,10), nrow=3)
  
  x <- list(x0, x1, x0)
  h <- type.analysis(x, 0.05)
  expect_equivalent(h$Type, 1)
  
})

