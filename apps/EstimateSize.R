# Invoke with Rscript EstimateSize.R input.file [ maxErrorRate(0.5 default)]
#   do *not* trust any estimate with a predicted error rate that high
#   unless the data is *really* bad
#

# read the input file
args <- commandArgs(TRUE)
file <- args[1]

# set the maximum error rate for the model
maxErrorRate <- 0.50
if (length(args) == 2) {
  maxErrorRate <- args[2]
}

# read in the data
inputData <- read.table(file, header = TRUE, sep="\t")

# define the asymptotic function model and its gradient
# uniqueKmers = fun3(rawKmers, errorRate, genomeSize)
fun3 <- function(x, ax, bx) {
  cx  <- (( bx + ax - 1) / bx)
  val <- ax * x + bx - bx * ( cx ) ^ x
  # derivatives from http://www.solvemymath.com/online_math_calculator/calculus/derivative_calculator/index.php
  dy.da <- ( ( -bx * ( cx ^ x ) + bx + ax - 1) * x ) / (bx + ax - 1)
  dy.db <- ( (ax - 1) * ( cx ^ x ) * x + ( - bx - ax + 1) * ( cx ^ x ) + bx + ax - 1) / ( bx + ax - 1)
  #dy.dx <- ax - bx * ( cx ^ x ) * log ( cx )

  attr(val, "gradient") <- cbind( ax = dy.da, bx = dy.db)
  val
}

# function to compare results to model function
predictedAsymptote <- function(rawKmerCount, errorRate, genomeSize) {
  errorRate * rawKmerCount + genomeSize
}

# Run the NLS fit for a data set
fitModel <- function(rawKmers, uniqueKmers) { 

  d <- data.frame(rawKmers=rawKmers, uniqueKmers=uniqueKmers)
  w <- log(unlist(rawKmers))
  w <- w / sum(as.numeric(w))


  cont <- nls.control(maxiter=500, tol=1e-5, minFactor=1/1024,warnOnly=TRUE)
  res <- nls(uniqueKmers ~ fun3(rawKmers, errorRate, genomeSize), d , trace=TRUE, control=cont, start=list(errorRate=0.0001, genomeSize=5000000), algorithm='port', lower=list(errorRate=0.00005, genomeSize=5000), upper=list(errorRate=maxErrorRate, genomeSize=100000000000), weights=w)

  #print(summary(res))
  #res
  coef(res)

}


fun4 <- function(x, ax, bx, fx) {
  E <- fx * bx
  E_x <- ((E-1)/E) ^ x
  B_x <- ((bx-1)/bx) ^ x
  val <- ax * E * (1 - E_x) + bx * (1 - B_x)
  dy.da  <- E * (1 - E_x)
  dy.df1 <- ax * bx * ( 1 - E_x )
  dy.df2 <- ax * E * E * E_x * x * ( (1/fx) - (E - 1) / (E * fx) ) / (E - 1)
  dy.df  <- dy.df1 - dy.df2
  dy.db1 <- ax * E * E * E_x * x * ( (1/bx) - (E - 1) / (bx * E) ) / (E - 1)
  dy.db2 <- B_x * ( (1/bx) - (bx - 1) / (bx*bx) ) * bx * bx * x / (bx - 1)
  dy.db3 <- ax * fx * (1 - E_x)
  dy.db  <- 0 - dy.db1 - dy.db2 + dy.db3 - B_x + 1

  attr(val, "gradient") <- cbind(ax = dy.da, bx = dy.db, fx = dy.df)
  val
}

# Run the NLS fit for a data set
fitModel4 <- function(rawKmers, uniqueKmers) { 

  d <- data.frame(rawKmers=rawKmers, uniqueKmers=uniqueKmers)
  w <- log(unlist(rawKmers))
  w <- w / sum(as.numeric(w))


  cont <- nls.control(maxiter=500, tol=1e-5, minFactor=1/1024,warnOnly=TRUE)
  res <- nls(uniqueKmers ~ fun4(rawKmers, errorRate, genomeSize, errorFactor), d , trace=TRUE, control=cont, start=list(errorRate=0.001, genomeSize=5000000, errorFactor=50), algorithm='port', lower=list(errorRate=0.0005, genomeSize=5000, errorFactor=50), upper=list(errorRate=maxErrorRate, genomeSize=100000000000, errorFactor=1000000), weights=w)

  #print(summary(res))
  #res
  coef(res)

}


# get the best fit with all the data
inputRawKmers <- inputData[1] 
inputUniqueKmers <- inputData[3] 
bestFit1 <- fitModel(inputRawKmers, inputUniqueKmers)
bestFit2 <- fitModel4(inputRawKmers, inputUniqueKmers)

# track every estimate with incremental data sets
start <- 3
errorRates <- vector(length=start-1)
genomeSizes <- vector(length=start-1)
errorFactors <- vector(length=start-1)
for(i in 1: start-1) {
  errorRates[i] <- NA
  genomeSizes[i] <- NA
  errorFactors[i] <- NA
}

numDataPoints <- dim(inputRawKmers)[1]
delta <- vector(length=numDataPoints-start+1)
delta[start-1] <- (inputUniqueKmers[start-1,] - inputUniqueKmers[start-2,]) / (inputRawKmers[start-1,] - inputRawKmers[start-2,])
deltadelta <- vector(length=numDataPoints-start+1)

for(i in start:numDataPoints) { 
  i2 <- i
  ss <- inputData[ seq(1,i), ]
  delta[i] <- (inputUniqueKmers[i,] - inputUniqueKmers[i-1,]) / (inputRawKmers[i,] - inputRawKmers[i-1,])
  deltadelta[i] <- (delta[i] - delta[i-1]) / (inputRawKmers[i,] - inputRawKmers[i-1,])
  
  rk <- ss[1]
  uk <- ss[3]
  x <- fitModel(rk, uk)
  if (x['errorRate'] < maxErrorRate) {
    errorRates <- append(errorRates, x['errorRate'])
    genomeSizes <- append(genomeSizes, x['genomeSize'])
    errorFactors <- append(errorFactors, x['errorFactor'])
  } else {
    errorRates <- append(errorRates, NA)
    genomeSizes <- append(genomeSizes, NA)
    errorFactors <- append(errorFactors, NA)
  }  
}

# plot the data to a PNG file
pl <- data.frame(RawKmers=inputRawKmers, EstimatedErrorRate=errorRates)
pl2 <- data.frame(RawKmers=inputRawKmers, EstimatedGenomeSize=genomeSizes)
pl3 <- data.frame(RawKmers=inputRawKmers, UniqueKmers=inputUniqueKmers)


plotfile <- paste(file, sep='', '.png')
png(filename=plotfile, width=640, height=480)

plot(pl3, xlim=c(0, max(pl3[,1])*1.10), ylim=c(0, max(pl3[,2])*1.10), xaxs='i', yaxs='i', main='Estimated Genome Size', sub=file)

predictedUniqKmers1 <- predictedAsymptote(inputRawKmers, bestFit1['errorRate'], bestFit1['genomeSize'])
predictedUniqKmers2 <- predictedAsymptote(inputRawKmers, bestFit2['errorRate'], bestFit2['genomeSize'])

#pred <- data.frame(RawKmers=inputRawKmers, Predicted=predictedUniqKmers)
#lines(pred, col='red')
abline(bestFit1['genomeSize'],bestFit1['errorRate'], col='red')
abline(bestFit2['genomeSize'],bestFit2['errorRate'], col='blue')
abline(bestFit2['genomeSize'] * bestFit2['errorFactor'], 0, col='orange')
l1<- paste("y_simple =",  round(bestFit1['genomeSize']), "+", bestFit1['errorRate'], "* x")
l2<- paste("y_saturated =",  round(bestFit2['genomeSize']), "+", bestFit2['errorRate'], "* x")
legend("topleft", col=c("black","red","blue","green"),lty=1,legend=c("DataPoints",l1,l2,"Estimated Genomic Kmers"))
points(pl2, col='green')

for(i in 20:numDataPoints) {
  if (!is.na(genomeSizes[i])) {
     #abline(genomeSizes[i], errorRates[i], col='yellow')
  }
}

# return the parameters for the best fit last
bestFit1
bestFit2

