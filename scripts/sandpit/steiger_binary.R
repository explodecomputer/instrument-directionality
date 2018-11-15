nid <- 500000
eff <- 0.015
g <- rbinom(nid, 2, 0.5)
y <- scale(g) * sqrt(eff) + rnorm(nid, sd=sqrt(1-eff))

cc1 <- y
cc1[y < quantile(y, 0.5)] <- 1
cc1[y >= quantile(y, 0.5)] <- 0

cc2 <- y
cc2[y < quantile(y, 0.1)] <- 1
cc2[y >= quantile(y, 0.1)] <- 0

cc3 <- y
cc3[y < quantile(y, 0.05)] <- 1
cc3[y >= quantile(y, 0.05)] <- 0

cor(g, y)^2
cor(g, cc1)^2
cor(g, cc2)^2
cor(g, cc3)^2

m1 <- glm(cc1 ~ g, family="binomial")
m2 <- glm(cc2 ~ g, family="binomial")
m3 <- glm(cc3 ~ g, family="binomial")

NagelkerkeR2(m1)
NagelkerkeR2(m2)
NagelkerkeR2(m3)


p <- 0.5
1 - (p^p * (1-p)^(1-p))^2
1 - exp(-m1$null/nid)

1 - (p^p * (1-p)^(1-p))^2 - 1
- exp(-m1$null/nid)

log(-(1 - (p^p * (1-p)^(1-p))^2 - 1))
- m1$null/nid


log(-(1 - (p^p * (1-p)^(1-p))^2 - 1)) * -nid
m1$null

or <- log(m1$coefficients[1])

or

f1 <- glm(cc1 ~ sample(g), family="binomial")
f1$dev
f1$null

log(f1$null)

m1
m2
m3
(m1^2 * 2 * 0.5 * 0.5 / var(cc1))
(m2^2 * 2 * 0.5 * 0.5 / var(cc2))
(m3^2 * 2 * 0.5 * 0.5 / var(cc3))

0.5/0.5

dnorm(0.1)^2


log(m1$null)


1-exp(-n1$null/nid)
NagelkerkeR2(m3)

(1 - NagelkerkeR2(m3)$R2) ^ (1/(2/nid))

cor()


index1 <- cc1 == 1
index1[sample(which(cc1 == 0), sum(index1), replace=FALSE)] <- TRUE
table(index1, cc1)

index2 <- cc2 == 1
index2[sample(which(cc2 == 0), sum(index2), replace=FALSE)] <- TRUE
table(index2, cc2)

index3 <- cc3 == 1
index3[sample(which(cc3 == 0), sum(index3), replace=FALSE)] <- TRUE
table(index3, cc3)



(m1c <- glm(cc1[index1] ~ g[index1], family="binomial"))
(m2c <- glm(cc2[index2] ~ g[index2], family="binomial")$coefficients[2])
(m3c <- glm(cc3[index3] ~ g[index3], family="binomial")$coefficients[2])

sum(g[index1]) / (2 * sum(index1))
sum(g[index2]) / (2 * sum(index2))
sum(g[index3]) / (2 * sum(index3))

(m1c^2 * 2 * 0.5 * 0.5 / var(cc1))^2
cor()


NagelkerkeR2()




set.seed(17)
i<-0
sim <- replicate(1e4, {
  while(TRUE) {
    x <- matrix(round(rexp(4), 2), 2, 2)
    if(all(rowSums(x) > 0) && all(colSums(x) > 0) && x[1,2]*x[2,1] > 0) break
  }
  x <- x / sum(x)
  beta <- rowSums(x)[1]
  gamma <- colSums(x)[1]
  rho <- x[1,1]*x[2,2] / (x[1,2]*x[2,1])

  y<- f(beta, gamma, rho)
  delta <- try(zapsmall(c(1, sqrt(crossprod(as.vector(x-y)))))[2])
  if ("try-error" %in% class(delta)) cat("Error processing ", x, "\n")
  delta
})
max(sim)

f(1.2, 0.2, 0.5)
f(log(1.2), 0.2, 0.5)

A <- 200
B <- 400
C <- 200
D <- 1000
OR <- A*D/(B*C)
N <- A+B+C+D
f((A+B)/N, (A+C)/N, OR)

LR <- (B/(A+B)) / (D/(D+C))

log(-(1 - (p^p * (1-p)^(1-p))^2 - 1)) * -nid

A/N
B/N
C/N
D/N


cont_from_or <- function(prevalence, allele_frequency, oddsratio, eps=1e-15)
{
	a <- oddsratio-1
	b <- (prevalence+allele_frequency)*(1-oddsratio)-1
	c_ <- oddsratio*prevalence*allele_frequency

	if (abs(a) < eps) {
	z <- -c_ / b
	} else {
	d <- b^2 - 4*a*c_
	if (d < eps*eps) s <- 0 else s <- c(-1,1)
	z <- (-b + s*sqrt(max(0, d))) / (2*a)
	}
	y <- vapply(z, function(a) zapsmall(matrix(c(a, allele_frequency-a, prevalence-a, 1+a-prevalence-allele_frequency), 2, 2)), matrix(0.0, 2, 2))
	i <- apply(y, 3, function(u) all(u >= 0))
	return(y[,,i])
}

cont_from_or(0.5, 0.5, exp(m1$coefficients[2]))
sum(cc1 & )/nid

(sum(cc1 & g==2)*2 + sum(cc1 & g==1)) / (2*nid)

