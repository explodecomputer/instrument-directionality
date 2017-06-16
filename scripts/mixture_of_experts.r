# Mixture of experts

library(tidyverse)

meth1 <- function(dat)
{
	mod <- summary(lm(y ~ x, dat))
	return(coefficients(mod)[2,])
}

meth2 <- function(dat)
{
	mod <- summary(lm(y ~ x + poly(x,2), dat))
	return(coefficients(mod)[3,])
}

meth3 <- function(dat)
{
	mod <- summary(lm(y ~ 0 + poly(x, 3), dat))
	return(coefficients(mod)[3,])
}

sims1 <- list()
for(i in 1:100)
{
	x <- rnorm(100)
	y <- x + rnorm(100) + runif(1)
	sims1[[i]] <- data.frame(x=x, y=y)
}

sims2 <- list()
for(i in 1:100)
{
	x <- rnorm(100)
	y <- x^2 + rnorm(100) + runif(1)
	sims2[[i]] <- data.frame(x=x, y=y)
}

sims3 <- list()
for(i in 1:100)
{
	x <- rnorm(100)
	y <- x^3 + rnorm(100)
	sims3[[i]] <- data.frame(x=x, y=y)
}


meth1(sims1[[1]])
meth2(sims1[[1]])
meth3(sims1[[1]])

meth1(sims2[[1]])
meth2(sims2[[1]])
meth3(sims2[[1]])

meth1(sims3[[1]])
meth2(sims3[[1]])
meth3(sims3[[1]])


sims <- c(sims1, sims2, sims3)

modsel <- rep(1:3, each=100)
x <- modsel + matrix(rnorm(15*300), 300, 15)

library(randomForest)

dat <- data.frame(y=as.factor(modsel), x)

rf <- randomForest(y ~ ., dat, ntrees=150)


mat <- matrix(0, 300, 3)
for(i in 1:300)
{
	mat[i, 1] <- -log10(meth1(sims[[i]])[4])
	mat[i, 2] <- -log10(meth2(sims[[i]])[4])
	mat[i, 3] <- -log10(meth3(sims[[i]])[4])
}

mod1 <- anova(lm(mat[,1] ~ x))
mod2 <- anova(lm(mat[,2] ~ x))

dat1 <- data.frame(y=mat[,1], x)
dat2 <- data.frame(y=mat[,2], x)
dat3 <- data.frame(y=mat[,3], x)

rf1 <- randomForest(y ~ ., dat1)
rf2 <- randomForest(y ~ ., dat2)
rf3 <- randomForest(y ~ ., dat3)

predict(rf1, dat1[1,-1])
predict(rf1, dat1[1,-1])
predict(rf1, dat1[1,-1])




sims1t <- list()
for(i in 1:100)
{
	x <- rnorm(100)
	y <- x + rnorm(100) + runif(1)
	sims1t[[i]] <- data.frame(x=x, y=y)
}

sims2t <- list()
for(i in 1:100)
{
	x <- rnorm(100)
	y <- x^2 + rnorm(100) + runif(1)
	sims2t[[i]] <- data.frame(x=x, y=y)
}

sims3t <- list()
for(i in 1:100)
{
	x <- rnorm(100)
	y <- x^3 + rnorm(100)
	sims3t[[i]] <- data.frame(x=x, y=y)
}

simst <- c(sims1t, sims2t, sims3t)

modsel <- rep(1:3, each=100)
xt <- modsel + matrix(rnorm(15*300), 300, 15)


p1 <- predict(rf1, data.frame(xt))
p2 <- predict(rf2, data.frame(xt))
p3 <- predict(rf3, data.frame(xt))



p <- cbind(p1, p2, p3)

psel <- apply(p, 1, which.max)


table(psel, modsel)

