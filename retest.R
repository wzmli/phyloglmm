library(lme4)

x <- rnorm(3000, sd=1)
obs <- rep(c("a","b","c"),each=1000)
y <- rep(c(100,200,300),each=1000)
z <- x+y

dat <- data.frame(x,z,y,obs)

mod1 <- lm(z~1)
mod2 <- lmer(z~1+(1|obs))


print(sd(x))
print(sd(y))


print(summary(dat))


print(summary(mod1))

print(summary(mod2))
