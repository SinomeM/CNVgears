
context(paste("Multiple methods CNVresults objects comparison and merge into",
              "a signle one. One of the two most complicated functions of the",
              "packege."))

# This is a very intricate task (and function) as CNVRs computation.
# I'll progressivewly try to sumplify it as the package mature

library(data.table)

sl <- data.table(sample_ID = "NA12878", sex = "2", role = "sibling",
                 fam_ID = "1463")

test_that("Overlapping calls are merged only if of comparable sizes", {
  pp <- data.table(chr = 22, start = c(20892562, 21060265),
                   end = c(20931293, 21090100), GT = 1, CN = 1,
                   sample_ID = "NA12878", meth_ID = "P")
  qq <- data.table(chr = 22, start = c(20909801, 21060265),
                   end = c(20931293, 21090100), GT = 1, CN = 1,
                   sample_ID = "NA12878", meth_ID = "Q")
  ii <- data.table(chr = 22, start = 20878863,
                   end = 21191994, GT = 1, CN = 1,
                   sample_ID = "NA12878", meth_ID = "I")
  class(pp) <- c("CNVresults", class(pp))
  class(qq) <- c("CNVresults", class(qq))
  class(ii) <- c("CNVresults", class(ii))

  invisible(capture.output(
    mm <- inter_res_merge(res_list = list(pp, qq, ii), sample_list = sl,
                          chr_arms = hg19_chr_arms)))
  rr <- data.table(chr = 22 , start = c(20892562, 21060265, 20878863),
                   end = c(20931293, 21090100, 21191994))
  expect_equal(nrow(mm), 3)
  # expect_equal(mm[,1:3], rr) # qui bisogna prima decidere ls strruttura delle colonne (num, int, chr, etc)
  expect_equal(mm$GT, c(1,1,1))
  expect_equal(mm$CN, c(1,1,1))
})

test_that("Overlapping calls are not merged if GT is not compatible ", {
  pp <- data.table(chr = 22, start = c(20892562, 21060265),
                   end = c(20931293, 21090100), GT = 2, CN = c(3,4),
                   sample_ID = "NA12878", meth_ID = "P")
  qq <- data.table(chr = 22, start = c(20909801, 21060265),
                   end = c(20931293, 21090100), GT = 1, CN = 1,
                   sample_ID = "NA12878", meth_ID = "Q")
  ii <- data.table(chr = 22, start = 20878863,
                   end = 21191994, GT = 1, CN = 1,
                   sample_ID = "NA12878", meth_ID = "I")
  class(pp) <- c("CNVresults", class(pp))
  class(qq) <- c("CNVresults", class(qq))
  class(ii) <- c("CNVresults", class(ii))

  invisible(capture.output(
    mm <- inter_res_merge(res_list = list(pp, qq, ii), sample_list = sl,
                          chr_arms = hg19_chr_arms)))
  rr <- data.table(chr = 22 ,
                   start = c(20892562, 20909801, 21060265, 21060265, 20878863),
                   end = c(20931293, 20931293, 21090100, 21090100, 21191994))

  expect_equal(nrow(mm), 5)
  # expect_equal(mm[,1:3], rr) # qui bisogna prima decidere ls strruttura delle colonne (num, int, chr, etc)
  expect_equal(mm$GT, c(2,1,2,1,1))
  expect_equal(mm$CN, c(3,1,4,1,1))
})

test_that("Test internal function get_region()", {
  pp <- data.table(chr = 22, start = c(20892562, 21060265),
                   end = c(20931293, 21090100))
  r1 <- CNVgears:::get_region(pp[1], 1)
  r2 <- CNVgears:::get_region(pp[2], 0.5)

  expect_equal(r1[1], 22)
  expect_equal(r1[2], 20892562)
  expect_equal(r1[3], 20931293)
  expect_equal(r1[4], 20931293-20892562+1)
  expect_equal(r2[4], (21090100-21060265+1)*0.5)
})

test_that("Test internal function compare_neighbors()", {

})
