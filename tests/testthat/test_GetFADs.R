library(NullTurnover)
library(phytools)

context("Test GetFADs()")

myTree <- pbtree(1, 0.5, n=10)

test_that("error catching works", {
          expect_error(GetFADs("aba"),"phylo")
          expect_error(GetFADs(pbtree(n=2)), "tip.label")
      })

test_that("output is numeric", 
          expect_is(GetFADs(myTree),"numeric")
)