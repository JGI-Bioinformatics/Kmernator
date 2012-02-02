args <- commandArgs(TRUE)

file <- args[1]

data <- read.table(file, header = TRUE, sep="\t")
rawKmers <- data[3]
uniqueKmers <- data[2]
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
res <- nls(uniqueKmers ~ fun3(rawKmers, errorRate, genomeSize), d , trace=TRUE, control=cont, start=list(errorRate=0.01, genomeSize=5000000), algorithm='port', lower=list(errorRate=0.0001, genomeSize=500), upper=list(errorRate=0.05, genomeSize=100000000000))

print(summary(res))

res

