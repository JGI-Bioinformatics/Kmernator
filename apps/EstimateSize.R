# Invoke with Rscript EstimateSize.R input.file [ maxErrorRate(0.25 default)]

args <- commandArgs(TRUE)
file <- args[1]

maxErrorRate <- 0.50
if (length(args) == 2) {
  maxErrorRate <- args[2]
}

inputData <- read.table(file, header = TRUE, sep="\t")

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

predictedAsymptote <- function(rawKmerCount, errorRate, genomeSize) {
  errorRate * rawKmerCount + genomeSize
}

fitModel <- function(rawKmers, uniqueKmers) { 

  d <- data.frame(rawKmers=rawKmers, uniqueKmers=uniqueKmers)
  w <- unlist(uniqueKmers)
  w <- w / sum(as.numeric(w))


  cont <- nls.control(maxiter=500, tol=1e-5, minFactor=1/1024,warnOnly=TRUE)
  res <- nls(uniqueKmers ~ fun3(rawKmers, errorRate, genomeSize), d , trace=TRUE, control=cont, start=list(errorRate=0.0001, genomeSize=5000000), algorithm='port', lower=list(errorRate=0.00005, genomeSize=5000), upper=list(errorRate=maxErrorRate, genomeSize=100000000000), weights=w)

  #print(summary(res))
  #res
  coef(res)

}

inputRawKmers <- inputData[1]
inputUniqueKmers <- inputData[3]
bestFit <- fitModel(inputRawKmers, inputUniqueKmers)
start <- 20
errorRates <- vector(length=start-1)
genomeSizes <- vector(length=start-1)
for(i in 1: start-1) {
  errorRates[i] <- NA
  genomeSizes[i] <- NA
}
for(i in start:dim(inputRawKmers)[1]) { 
  i2 <- i
  ss <- inputData[ seq(1,i), ]
  rk <- ss[1]
  uk <- ss[3]
  x <- fitModel(rk, uk)
  if (x['errorRate'] < maxErrorRate) {
    errorRates <- append(errorRates, x['errorRate'])
    genomeSizes <- append(genomeSizes, x['genomeSize'])
  } else {
    errorRates <- append(errorRates, NA)
    genomeSizes <- append(genomeSizes, NA)
  }  
}

pl <- data.frame(RawKmers=inputRawKmers, EstimatedErrorRate=errorRates)
pl2 <- data.frame(RawKmers=inputRawKmers, EstimatedGenomeSize=genomeSizes)
pl3 <- data.frame(RawKmers=inputRawKmers, UniqueKmers=inputUniqueKmers)

plotfile <- paste(file, sep='', '.png')
png(filename=plotfile, width=640, height=480)
plot(pl3, xlim=c(0, max(pl3[,1])*1.10), ylim=c(0, max(pl3[,2])*1.10), xaxs='i', yaxs='i', main='Estimated Genome Size', sub=file)

predictedUniqKmers <- predictedAsymptote(inputRawKmers, bestFit['errorRate'], bestFit['genomeSize'])
#pred <- data.frame(RawKmers=inputRawKmers, Predicted=predictedUniqKmers)
#lines(pred, col='red')
abline(bestFit['genomeSize'],bestFit['errorRate'], col='red')
l<- paste("y =",  round(bestFit['genomeSize']), "+", bestFit['errorRate'], "* x")
legend("topleft", col=c("black","red","green"),lty=1,legend=c("DataPoints",l,"Estimated Genomic Kmers"))
points(pl2, col='green')

bestFit

