# Invoke with Rscript EstimateSize.R input.file [ maxErrorRate(0.075 default)]

args <- commandArgs(TRUE)
file <- args[1]

maxErrorRate <- 0.075
if (length(args) == 2) {
  maxErrorRate <- args[2]
}

data <- read.table(file, header = TRUE, sep="\t")
rawKmers <- data[1]
uniqueKmers <- data[3]
d <- data.frame(rawKmers=rawKmers, uniqueKmers=uniqueKmers)

fun3 <- function(x, ax, bx) {
  val <- ax * x + bx - bx * (( ax + bx - 1) / bx ) ^ x
  cx <- (( bx + ax - 1) / bx)
  # derivatives from http://www.solvemymath.com/online_math_calculator/calculus/derivative_calculator/index.php
  dy.da <- ( ( -bx * ( cx ^ x ) + bx + ax - 1) * x ) / (bx + ax - 1)
  dy.db <- ( (ax - 1) * ( cx ^ x ) * x + ( - bx - ax + 1) * ( cx ^ x ) + bx + ax - 1) / ( bx + ax - 1)
  #dy.dx <- ax - bx * ( cx ^ x ) * log ( cx )

  attr(val, "gradient") <- cbind( ax = dy.da, bx = dy.db)
  val
}

cont <- nls.control(maxiter=500, tol=1e-5, minFactor=1/1024,warnOnly=TRUE)
res <- nls(uniqueKmers ~ fun3(rawKmers, errorRate, genomeSize), d , trace=FALSE, control=cont, start=list(errorRate=0.01, genomeSize=5000000), algorithm='port', lower=list(errorRate=0.0001, genomeSize=500), upper=list(errorRate=maxErrorRate, genomeSize=100000000000))

#print(summary(res))

#res

coef(res)
