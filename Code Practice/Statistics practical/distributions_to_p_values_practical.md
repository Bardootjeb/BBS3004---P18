% Simulating sampling to get p-values

In this practical you will use R's sampling functions to:

1. Determine the distributions of test statistics under the null hypothesis.
2. Use these distributions to infer p-values.
3. Perform *in silico* experiments by sampling from non-identical
distributions.
4. Compare the inferred p-values which those obtained from R's built in
functions.

Note that the above are objectives and that you will not do them strictly in order.

## Sample sets of numbers from different distributions.

To start with we will define the distributions from which we will sample, the
sample sizes (the replicate number, `k`) and the number of samples (`N`). To
keep things simple I will describe sampling from normal distributions using
`rnorm` as it is reasonably easy to understand how these are defined. You can
also try to see what happens when you use other distributions.

```R
## Normal distributions are defined by their central value (mean) and their
## standard deviation. We will sample from several different distributions
## and to facilitate this we first create a matrix with two columns defining
## the mean values and standard deviations:

## norm.par, for normal parameters
norm.par <- cbind("mean"=c(5, 7, 10), "sd"=c(2, 2, 5))

## As always, look at the data structure you have created before doing
## anything else. Make sure that it is what you expect from the above
## command.
## 'cbind', is short for column bind; you should know what it does.

## To sample from the distributions defined by norm.par we first want to define
## our sample sizes and numbers. I will use k=3 and N=100000.
## Please consider using other values.

k <- 3
N <- 100000

## We can then use `lapply` to obtain a set of samples for each row of
## norm.par

norm.samples <- lapply(1:nrow(norm.par), function(i){
    t(sapply(1:N, function(j){ rnorm(n=k, mean=norm.par[i,'mean'],
                                     sd=norm.par[i,'sd']) }))
})

## as usual have a look at what you obtained from this:
## you can do:
summary(norm.samples)  ## ??

## but I prefer:
class(norm.samples)

sapply(norm.samples, class)

sapply(norm.samples, dim)

## I usually like to have a look at a bit of the objects:
head(norm.samples[[1]])
head(norm.samples[[3]])

### NOTE
## I used `lapply` on the rows of the norm.par table.
## Generally to apply a function to the rows of a matrix
## we would use `apply`. However, here I wanted a list,
## but apply gave me a merged matrix. The rules of what
## the different apply functions give you can be difficult
## to remember. So always check what you obtain.

## Look at the statistical properties of the samples:

sapply(norm.samples, mean)

sapply(norm.samples, sd)

## Are these what you expect? If they are not, or you are not sure
## what to expect, then please ask. It should be fairly obvious.

## to look at the distributions you can try:
norm.samples.h <- lapply(norm.samples, hist, plot=FALSE, breaks=100)

## lets plot the density component of these:
## but before that inspect the return values so that you understand
## what you are doing:

class(norm.samples.h)

class(norm.samples.h[[1]])

names(norm.samples.h[[1]])

sapply(norm.samples.h[[1]], class)

sapply(norm.samples.h[[1]], length)


with(norm.samples.h[[1]], plot(mids, density, type='l'))

## the above is equivalent to writing:
plot( norm.samples.h[[1]]$mids, norm.samples.h[[1]]$density, type='l')

## we can add lines to this plot. In this case we can try the values from
## the other histograms and the lines() function:
with(norm.samples.h[[2]], lines(mids, density, col='red'))
with(norm.samples.h[[3]], lines(mids, density, col='blue'))

## hmm, we can't see very much of the blue samples. We can increase the xlimits
## by specifying these in the first plot() command that sets up the plotting
## surface:
with(norm.samples.h[[1]], plot(mids, density, type='l',
                               xlim=range(norm.samples.h[[3]]$mids)))
with(norm.samples.h[[2]], lines(mids, density, col='red'))
with(norm.samples.h[[3]], lines(mids, density, col='blue'))

## since we know how these were sampled, we actually know the true distributions
## we can use these to add on the theoretical distributions as above.
## This makes use of the dnorm function which gives you the density of
## the distribution at specifie points.
with(norm.samples.h[[1]], lines( mids,
                                dnorm(mids, mean=norm.par[1,'mean'],
                                      sd=norm.par[1,'sd']),
                                col='black', lty=2, lwd=3))
with(norm.samples.h[[2]], lines( mids,
                                dnorm(mids, mean=norm.par[2,'mean'],
                                      sd=norm.par[2,'sd']),
                                col='red', lty=2, lwd=3))
with(norm.samples.h[[3]], lines( mids,
                                dnorm(mids, mean=norm.par[3,'mean'],
                                      sd=norm.par[3,'sd']),
                                col='blue', lty=2, lwd=3))

## Note the use of lty (line type) and lwd (line width) to show more lines.
## To find out how to use these and other graphical parameters, look at
## the help for the 'par()' function (?par). There are many to choose form.

## The distributions should look very similar to the theoretical values.
## This is because we used a large value of N.
## Try the same, but with a much smaller value of N

## You should also note that the commands given above are very similar to each other,
## differing primarily by a number. That's usually a sign of bad code. It might
## be better to rewrite it as a loop or to specify a function or something similar.
## Especially as we are repeating the same number many times in the same line.

```

When I comment on code being bad, ugly and horrible that usually means two
things:

1. The code is repetitive and inefficient.
2. It is easy to make mistakes. In the code snippet above I have copied and
   pasted a line several times and then modified it by replacing '1' with '2'
   or '3' in several places. That's bad, because it's really easy to miss
   out one of the substitutions. It's also bad because it's not reusable.
   
Code that is difficult to understand is also bad. Even if it can be
demonstrated to be correct you will have a hard time understanding it in the
future; and an even harder time trying to modify it to do something different.

You should always think about whether the code you write is elegant or ugly,
and how you can improve it to be more understandable, reusable and
safer[^safe].

[^safe]: If you hang around computer scientists you may find them talking
    about *code safety*. What they are they referring to is how likely
    mistakes are in a given program or language. Mistakes are dangerous as
    they can lead to more than just error statements.

You should be aware that the way that I have written the code above is not
very efficient, nor very elegant. I have *tried* to make it simple by making
everything explicit. But certainly, the sampling could be done much more
efficiently. For example to get a matrix with three columns and N rows with
values sampled from a single distribution it would be simpler to write:

```R
## but note that you need to define mean and sd values first
## or replace the mean and sd with something else in
## the following command.
s <- matrix( rnorm(3 * N, mean=mean, sd=sd), ncol=3 )
```

Though obviously you would have to use variables that exist in your session.

To get a matrix with three columns with numbers sampled from three different
normal distributions you could also do:

```R
## use norm.par, and N as before
s <- matrix(rnorm(nrow(norm.par) * N,
                  mean=norm.par[,'mean'], sd=norm.par[,'sd']),
            ncol=nrow(norm.par), byrow=TRUE )

## check that that did the correct thing:
dim(s)
apply(s, 2, mean)
apply(s, 2, sd)
## compare to:
norm.par
```

In this statement I make a single call to `rnorm()`, with the length set to
`N` multiplied by the number of different distributions that I want. I then
specify the mean and standard deviations as vectors taken from `norm.par`.
These will be recycled, so every third value will have the first mean and
standard deviation, and so on. If I fill that into a matrix with the
appropriate number of columns by row, then I'll end up with what I want.

That is much more efficient, but it is more difficult to understand and so
I've avoided it. Feel free to replace the code I'm giving you with more
efficient one when you think you know how. Just make sure to check that it's
really doing what you expect.

## Obtain some distributions under the *null* hypothesis

The values you obtained by sampling in the previous section can be considered
as simulations of experimental sampling of different distributions. To compare
these to sampling under the *null* hypothesis we can sample two sets of
samples from each individual distribution. We can also consider sampling from
a merged distribution. I don't *actually* know how this should be done, but
I guess that a combined distribution will have mean of the means and a
standard deviation that is the square root of the sums of the individual variances:

$$
sd_m = \sqrt{ \frac{\sum_{i=1}^{n} \sigma_i^2}{n} }
$$

Where $sd_m$ is the merged variance and $\sigma_i$ is the standard deviation
for the $i^{th}$ distribution.

Lets call these the *null* samples or something similar:

```R
## first lets define the parameters of a merged distribution:
mean.m <- mean( norm.par[,'mean'] )
sd.m <- sqrt( sum( norm.par[,'sd']^2 ) / nrow(norm.par) )

## have a look and see how that looks..
## the with here is simply to obtain the mids component
## from the histogram.
with(norm.samples.h[[3]], plot(mids, dnorm(mids, mean=mean.m,
                                           sd=sd.m),
                               type='l', col='brown',
                               ylim=range(norm.samples.h[[1]]$density )))

## and add the other lines to see how they fit:
with(norm.samples.h[[1]], lines(mids, density, type='l'))
with(norm.samples.h[[2]], lines(mids, density, type='l', col='red'))
with(norm.samples.h[[3]], lines(mids, density, type='l', col='blue'))

## the merged distribution doesn't look like the distribution
## of the full set of values. This is because it doesn't take
## into account the differences in the mean values. 

## To get that you could do:
mean( unlist(norm.samples) ) ## 7.335299 (the same)
sd( unlist(norm.samples) ) ## 3.897953 (larger)

## Add the merged distribution to the norm.par
## DO NOT excute this more than once!
norm.par <- rbind(norm.par, c(mean.m, sd.m))
## and look at it!

## then lets get two samples from each individual distribution here.
## (this is a bit slow, and could be done much, much faster)
null.samples <- apply( norm.par, 1, function(p){
    list('x'=t(sapply(1:N, function(i){ rnorm(n=k, mean=p['mean'], sd=p['sd']) })),
         'y'=t(sapply(1:N, function(i){ rnorm(n=k, mean=p['mean'], sd=p['sd']) })))
})

## and look at the structure of what you got:
class(null.samples)

class(null.samples[[1]])

sapply(null.samples, names)

dim(null.samples[[1]]$x)
dim(null.samples[[1]]$y)

head(null.samples[[1]]$x)
## and so on..

## Hopefully you are getting what you expect.

```

## Calculate some summary statistics for the samples

Now that we have lots of samples we can calculate some summary statistics and
look at their statistical properties. First lets calculate sample means and
sample variances for all of the individual sets:

```R

## For the norm.samples it's rather simple:
norm.samples.mean <- lapply( norm.samples, rowMeans )

## unfortunately, there is nor rowVar function, so we'll need
## a nested apply statement:
norm.samples.var <- lapply( norm.samples, function(x){
    apply(x, 1, var) })

## That was very slow. There is a much faster way to do this calculation
## using vectorised functions. I will show this in the last section 

## For the null samples we have a list of two entries for each distribution
## So we will need deeper nesting. This can be a bit confusing, but do
## try to work out how the following works.

null.samples.mean <- lapply( null.samples, function(x){
    lapply(x, rowMeans)
})

null.samples.var <- lapply( null.samples, function(x){
    lapply(x, function(y){
        apply(y, 1, var)
    })
})

## Have a look at the resulting data structures. They are not actually
## that easy to visualise. Use class, length, dim, head to work out
## what they represent.

### And look at these statistics:
## First make sure that they are what you expect:
class( norm.samples.mean )
length( norm.samples.mean )
sapply( norm.samples.mean, class )
sapply( norm.samples.mean, length )

class( null.samples.mean )
length( null.samples.mean )
sapply( null.samples.mean, class )
sapply( null.samples.mean[[1]], class )
names( null.samples.mean[[1]] )

head(null.samples.mean[[1]]$x)
head(null.samples.var[[1]]$x)

## then have a look at how the values are distributed:
## to get a sufficiently large plotting area you might want to:
x11()
## and then:
par(mfrow=c(1,3)) ## this will give you a row of three plots
hist( norm.samples.mean[[1]] )
hist( norm.samples.mean[[2]] )
hist( norm.samples.mean[[3]] )

hist( norm.samples.var[[1]] )
hist( norm.samples.var[[2]] )
hist( norm.samples.var[[3]] )

## these may be a bit more familiar:
hist( sqrt(norm.samples.var[[1]]) )
hist( sqrt(norm.samples.var[[2]]) )
hist( sqrt(norm.samples.var[[3]]) )
## Why does the use of sqrt( ) make sense here?

### These distributions should look familiar to you. If it's
### not immediately obvious then consider how you obtained the
### data sets.

## You can also look at the relationship between mean and standard deviation
## for a given set. I'm not sure if you'll see much here, but that is
## a good reason for having a look:

par(mfrow=c(1,1))
plot( norm.samples.mean[[1]], sqrt(norm.samples.var[[1]]), cex=0.3)
## too many points to see anything. There does not seem to much relationship
## between the values. (The cex arguments sets the size of the points).

## This uses alpha blending to get a better idea of the number of overlapping
## points. The colour is specified by it's red, green, blue and alpha
## components.
plot( norm.samples.mean[[1]], sqrt(norm.samples.var[[1]]),
     cex=0.3, col=rgb(0,0,0,0.04) )

## If you prefer blue points:
plot( norm.samples.mean[[1]], sqrt(norm.samples.var[[1]]),
     cex=0.3, col=rgb(0,0,1,0.04) )

## So there doesn't look like there is any relationship between the
## variance and the mean values within a given set. That's pretty
## reasonable.

### DO ALSO look at the distributions of the null samples. They
### should look the same as here. However, you will need to
### change how you access the values since we have two sets
### of numbers for each distribution. If you can't do this,
### then you DO need to ask for help (but not until you've tried
### it yourself).

```

## Uncertainty in summary stats

The mean and variances calculated in the previous session are summary
statistics that are used to estimate the parameters of the distribution
from which they were sampled. We can consider how variable they are
as well:

```R
mean.sd <- sapply(norm.samples.mean, sd)
var.sd <- sapply(norm.samples.var, sd)

## we can check if the uncertainty itself forms a normal distribution
par(mfrow=c(1,2))
m.h1 <- hist(norm.samples.mean[[1]], probability=TRUE)
with(m.h1, lines(mids, dnorm(mids, mean=mean(norm.samples.mean[[1]]),
                             sd=mean.sd[1]), col='red', lwd=3))

m.h2 <- hist(norm.samples.mean[[2]], probability=TRUE)
with(m.h2, lines(mids, dnorm(mids, mean=mean(norm.samples.mean[[2]]),
                             sd=mean.sd[2]), col='red', lwd=3))

```

Not surprisingly the distribution of sample means also follows a normal
distribution. We can ask how this uncertainty changes with the number
of samples (i.e., replicates `k`).

```R
## for this we will create a new set of samples, called k.smp
kv <- 1:10
kv.smp <- lapply(kv, function(x){ t(sapply(1:N, function(i){
    rnorm(x, mean=norm.par[1,'mean'], sd=norm.par[1,'sd'])}))
})

## check what you got:
sapply(kv.smp, dim)

## sometimes R does strange things. Here I'll want to do:
kv.smp[[1]] <- t(kv.smp[[1]])
## in order to make things a bit more consistent
sapply(kv.smp, dim)

## we can then calculate the sample means and variances
kv.smp.means <- sapply(kv.smp, rowMeans)
kv.smp.var <- sapply(kv.smp, function(x){ apply(x, 1, var) })
## note that these are matrices (check class and dim)

## you can look at various aspects:
par(mfrow=c(1,3))
hist(kv.smp.means[,1])
hist(kv.smp.means[,5])
hist(kv.smp.means[,10])

hist(kv.smp.var[,2])
hist(kv.smp.var[,5])
hist(kv.smp.var[,10])

## you can also consider:
hist(sqrt(kv.smp.var[,2]))
hist(sqrt(kv.smp.var[,5]))
hist(sqrt(kv.smp.var[,10]))

## then get the variances of these estimates
kv.smp.means.var <- apply( kv.smp.means, 2, var )
kv.smp.var.var <- apply( kv.smp.var, 2, var )
kv.smp.sd.var <-  apply( sqrt(kv.smp.var), 2, var )

## and have a look:
par(mfrow=c(1,3))
plot(kv, kv.smp.means.var, type='b')
plot(kv, kv.smp.var.var, type='b')
plot(kv, kv.smp.sd.var, type='b')

plot(1/kv, kv.smp.means.var, type='b')
plot(1/kv, kv.smp.var.var, type='b')
plot(1/kv, kv.smp.sd.var, type='b')

## note also:
plot(sqrt(1/kv), sqrt(kv.smp.means.var), type='b')
## and
abline(0, 2, col='red')
## consider where the two came from.

```

From these we can see that there is a linear relationship between the standard
deviation of the sample means and the inverse of the sample number (k). When
the sample number is 1, then the variance of the estimate is the same as of
the underlying distribution. The standard error is the uncertainty in the
estimate of the mean expressed as a standard deviation (square root of the
variance). It is thus the (estimated) population standard deviation
divided by the square root of the sample number.

## Summary stats to differential stats

Having obtained the summary statistics we can now start to use these to
calculate statistics that relate to how different the sets of values are. For
this exercise we will consider the sample sets as being paired. That is that
we consider the rows of the sample matrices to somehow be linked to each
other; a bit as if each row represented one gene and the columns are replicate
measurements of each individual gene. For the time being lets not consider
these as test statistics per se, but just as descriptive terms.

The simplest statistic of difference is just that; the difference ($x - y$). 
But we can consider how this changes with the variance. Here that may not be
that interesting as all the value sets should have similar variances being
obtained from the same samples. You can think about how you could modify the
sampling procedure to have different variances for every row.

```R

## Lets calculate all against all differences for mean values. This is redundant
## (i.e. we compare set one against itself, set one against set two, and
## set two against set one, and so on), and wasteful. Later on you will
## see how you can remove this redundancy.

norm.samples.diff <- lapply( norm.samples.mean, function(x){
    lapply( norm.samples.mean, function(y){
        x - y
    })
})

## then inspect the object obtained to make sure it is what you think it is.
## use the normal means of looking at it.


## after you've worked out what the lists contain have a look at:
head(norm.samples.diff[[1]][[1]]) ## ?? 

head(norm.samples.diff[[1]][[2]])
head(norm.samples.diff[[2]][[1]])
## Note the redundancy

## what about:
plot(norm.samples.diff[[1]][[2]], norm.samples.diff[[2]][[1]])

## These values should make sense to you.. 

## Have a look at the distributions of some of these:
hist( norm.samples.diff[[1]][[2]])
hist( norm.samples.diff[[1]][[3]])

## consider the relationship to norm.par

## maybe consider their relationship to their standard deviations:
## but here we have two sets of standard deviations. Perhaps we should
## combine them. Try to work out what this does:

plot(sqrt(norm.samples.var[[1]] + norm.samples.var[[2]] ) * sqrt(1/k),
     norm.samples.diff[[1]][[2]], cex=0.3, col=rgb(0,0,0,0.04) )
## this draws a line. Guess which one
abline(a=0, b=1, col='red')

## no particular pattern here; which is probably what we want to see

## We could ask the question, do we see bigger differences between
## the different samples. In this case I'll use a quantiles approach:

plot(sort(norm.samples.diff[[1]][[2]]), sort(norm.samples.diff[[1]][[3]]), cex=0.3)
abline(a=0, b=1, col='red')
## That makes a lot of sense as the there is a bigger difference between distributions
## 1 and 3 than between 2 and 3.

## We can also calculate something like the t-statistic from these values. We can compare
## this later to the built in statistic.

norm.samples.t <- lapply(1:length(norm.samples), function(i){
    lapply(1:length(norm.samples), function(j){
        norm.samples.diff[[i]][[j]] /
            (sqrt(norm.samples.var[[i]] + norm.samples.var[[j]]) * sqrt(1/k))
    })
})

## again look at the distributions of these values:
hist(norm.samples.t[[1]][[1]]) ## ??
head(norm.samples.t[[1]][[1]]) ## ahah..

hist(norm.samples.t[[1]][[2]])
hist(norm.samples.t[[1]][[3]])
hist(norm.samples.t[[3]][[1]]) ## note the relationships

## Try to look at some quantiles examples as well (eg. just sort the values).

```

## Differential stats when the *null* hypothesis is true

In the previous section you will have obtained a load of numbers giving some
indication of difference and variance for the individual samples. Next we will
do the same but for samples obtained by sampling from within the same
distributions.

```R

null.samples.diff <- lapply(null.samples.mean, function(x){
    x$x - x$y })

null.samples.t <- lapply(1:length(null.samples.diff), function(i){
    null.samples.diff[[i]] /
        (sqrt(null.samples.var[[i]]$x + null.samples.var[[i]]$y)
            * sqrt(1/k))
})

## and have a quick look at these as usual..
hist(null.samples.diff[[1]])
hist(null.samples.diff[[2]])

hist(null.samples.t[[1]])
hist(null.samples.t[[2]])

plot( sort(null.samples.t[[1]]), cex=0.3)

## compare these to the t-distribution from R
null.samples.t.h <- lapply(null.samples.t, hist, breaks=100, plot=FALSE)

## then do:
par(mfrow=c(2,2))
invisible( lapply(null.samples.t.h, function(x){
    plot(x, freq=FALSE)
    with(x, lines(mids, dt(mids, df=4), col='red', lwd=2))
}))
## the invisible, stops the R from printing to screen anything returned
## by the lapply statement.

## compare these to what you get with statistics obtained from different
## distributions

plot( sort(norm.samples.t[[1]][[2]]), cex=0.3 )

## and a quantiles quantiles plot
## this compares the sorted values for each set
plot( sort(null.samples.t[[1]]), sort(norm.samples.t[[1]][[2]]), cex=0.3)
abline(a=0, b=1, col='red')
## that's a bit weird:

## absolute values may be more reasonable...
plot( sort(abs(null.samples.t[[1]])), sort(abs(norm.samples.t[[1]][[2]])), cex=0.3)
abline(a=0, b=1, col='red')
## As we might expect the t-values have lower magnitude when we sample
## from the same distributions (the null samples).

plot( sort(abs(null.samples.t[[2]])), sort(abs(norm.samples.t[[1]][[2]])), cex=0.3)
abline(a=0, b=1, col='red')

## We can also look at just the difference values:
## where we see the same expected pattern.
plot(sort(abs(null.samples.diff[[1]])), sort(abs(norm.samples.diff[[1]][[2]])), cex=0.3)
abline(a=0, b=1, col='red')

plot(sort(abs(null.samples.diff[[1]])), sort(abs(norm.samples.diff[[1]][[3]])), cex=0.3)
abline(a=0, b=1, col='red')

## You should also look at the following. This is actually important!
## And explains why the t-test is useful:
par(mfrow=c(1,2))
plot(sort(null.samples.diff[[1]]), cex=0.3, type='l' )
lines(1:N, sort(null.samples.diff[[2]]), cex=0.3, col='red' )
lines(1:N, sort(null.samples.diff[[3]]), cex=0.3, col='blue' )
lines(1:N, sort(null.samples.diff[[4]]), cex=0.3, col='green' )
## 
plot(sort(null.samples.t[[1]]), cex=0.3, type='l' )
lines(1:N, sort(null.samples.t[[2]]), cex=0.3, col='red' )
lines(1:N, sort(null.samples.t[[3]]), cex=0.3, col='blue' )
lines(1:N, sort(null.samples.t[[4]]), cex=0.3, col='green' )
## remember the variances associated with each distribution.
## check norm.par if you don't remember
##
## the t-statistic removes the effect of the variance
## and means that it doesn't matter here.
```

## Estimating p-values using the *null* distributions

We can use the null distributions to estimate the likelihood of observing
specific t-values when the samples are obtained from the same
distribution. However, it's not clear as to which distribution we should use
or if this even matters. We can consider other ways of defining our *null*
t-distributions later. Even if it doesn't matter here, it may be useful in
other circumstances.

Here I will first consider using the t-values. However, it might be possible
to also use the raw difference values if the *null* distribution can be
constrained in an appropriate manner.

```R
## The p-value is the likelihood of observing a given or more extreme statistic
## under the null hypothesis. That's basically the proportion of tests under the
## null hypothesis that have this statistic. Hence given the null distribution
## and a test statistic we can define the following function:

## null.d should be a set of statistics obtained under the null hypothesis
## stat.v is a statistic to be tested.
## note that here we assume a two sided test and that we are only interested
## in magnitudes.
stat.to.p <- function(stat.v, null.d){
    ## note that this function is monstrously inefficient
    ## but we don't care here as we only want to use it for
    ## demonstration purposes
    sum( abs(null.d) >= abs(stat.v) ) / length(null.d)
}

## lets try it out..
stat.to.p( norm.samples.t[[1]][[2]][1], null.samples.t[[1]] )
stat.to.p( norm.samples.t[[1]][[3]][3], null.samples.t[[1]] )

## the above function is really slow, so we can try
## to time it
system.time(
    tmp <- sapply(norm.samples.t[[1]][[2]][1:1000], stat.to.p,
                  null.d=null.samples.t[[1]] )
)
## on my machine that gave:
##  0.384   0.000   0.384
## so running all of the samples will take around 40 seconds.

## It might take something like 6 minutes
## to do the full sets. So lets do 10%
## and also avoid tests of samples against
## themselves. (Try to understand how this
## loop construct works).
norm.samples.p <- lapply(1:(length(norm.samples.t)-1), function(i){
    lapply((i+1):length(norm.samples.t), function(j){
        sapply(norm.samples.t[[i]][[j]][1:10000], stat.to.p,
               null.d = null.samples.t[[1]] )
    })
})

## and lets do the same for the null samples..
null.samples.p <- lapply(null.samples.t, function(x){
    sapply(x[1:10000], stat.to.p, null.d = x )
})

## have a look at distributions..
par(mfrow=c(1,3))
hist( norm.samples.p[[1]][[1]], breaks=20)
hist( norm.samples.p[[1]][[2]], breaks=20)
hist( norm.samples.p[[2]][[1]], breaks=20)

## Note that in this case, the null hypothesis is false in all
## cases. Yet there are not many p-values below 0.05.
## Consider what this tells you.

## compare to:
hist( null.samples.p[[1]], breaks=20 )
hist( null.samples.p[[2]], breaks=20 )
hist( null.samples.p[[3]], breaks=20 )
## and that is kind of what you expect to see.

## can also try:
hist( log10(norm.samples.p[[1]][[2]]), breaks=20 )
plot(sort(log10(norm.samples.p[[1]][[2]])), type='p', cex=0.3)

par(mfrow=c(2,2))
plot(sort(log10(null.samples.p[[1]])), type='p', cex=0.3)
plot(sort(log10(null.samples.p[[2]])), type='p', cex=0.3)
plot(sort(log10(null.samples.p[[3]])), type='p', cex=0.3)
plot(sort(log10(null.samples.p[[4]])), type='p', cex=0.3)
## note the range of p-values from the null p-values
## This is why multiple testing correction is necessary.

## you can also ask for each one of these, how many tests
## have p-values less than 0.05:

sapply(null.samples.p, function(x){ sum(x <= 0.05) })

## and the proportion of values that this represents:
sapply(null.samples.p, function(x){ sum(x <= 0.05) / length(x) })

## this defines the p-value;

```

Note that I used `null.samples.t[[1]]` as the *null* distribution for all of
the tests. Why is that OK?

## Compare to the built in t-test

R has a `t.test()` function built in. We can try to compare the values we get
from that to our own. 

```R
## we will restrict ourselves to the first 10,000 rows for each data set:

norm.samples.tt <- lapply(1:(length(norm.samples)-1), function(i){
    lapply((i+1):length(norm.samples), function(j){
        lapply(1:10000, function(k){
            t.test(norm.samples[[i]][k,], norm.samples[[j]][k,],
                   alternative="two.sided", var.equal=FALSE)
        })
    })
})

null.samples.tt <- lapply(null.samples, function(x){
    lapply(1:10000, function(i){
        t.test(x$x[i,], x$y[i,], alternative="two.sided", var.equal=TRUE)
    })
})

norm.samples.tt.p <- lapply(norm.samples.tt, function(x){
    lapply(x, function(y){
        sapply(y, function(tt){ tt$p.value })
    })
})

null.samples.tt.p <- lapply(null.samples.tt, function(x){
    sapply(x, function(tt){ tt$p.value })
})

## check some distributions
par(mfrow=c(1,2))
hist(null.samples.tt.p[[1]])
hist(norm.samples.tt.p[[1]][[2]])

## compare to our own estimates:
par(mfrow=c(1,2))
plot( norm.samples.p[[1]][[1]], norm.samples.tt.p[[1]][[1]],
     cex=0.3, col=rgb(0,0,0,0.04) )
abline(0, 1, col='red')
plot( norm.samples.p[[1]][[2]], norm.samples.tt.p[[1]][[2]],
     cex=0.3, col=rgb(0,0,0,0.04) )
abline(0, 1, col='red')

par(mfrow=c(1,2))
plot( null.samples.p[[1]], null.samples.tt.p[[1]],
     cex=0.3, col=rgb(0,0,0,0.1) )
abline(0, 1, col='red')
plot( null.samples.p[[2]], null.samples.tt.p[[2]],
     cex=0.3, col=rgb(0,0,0,0.1) )
abline(0, 1, col='red')

## The p-values we have estimated are well correlated with the
## ones from the built in test; but for the norm.samples they
## are not identical. Ours are slightly higher. If you look at:

norm.samples.tt[[1]][[1]][[1]]

## you may be able to work out why.

## We can look at the t-statistics
## obtained by the t-test:

norm.samples.tt.t <- lapply(norm.samples.tt, function(x){
    lapply(x, function(y){
        sapply(y, function(tt){ tt$statistic })
    })
})

## and ask how that compares to what we obtained:
plot( norm.samples.t[[1]][[2]][1:10000], norm.samples.tt.t[[1]][[1]] )
abline(0, 1, col='red')

## so the t-stats agree better than the p-values.
## what do you think is going one here?

```

When we compare our estimated p-values to the ones from the built in
t-test we find correlated but not identical values.
Sometimes the differences are quite large.

To work out why this happens you can have a look at the direct output
of the `t.test()` function as above:

```R
## you can look at any of them:
norm.samples.tt[[1]][[1]][[1]]
```

If you look at the details you will find that the test uses a weird looking
degrees of freedom value (it's non-integral). If you read the documentation
you will find that this is because we have set the `var.equal` option to
`FALSE`, and that `R` will use the Welch modification of degrees of freedom
to correct for the estimated difference in the variance. This highlights a
potential flaw in the simulate we use to derive the distribution of the
test statistic under the *null* hypothesis. Most models are flawed in some
way or other, and it is important to know that the p-value obtained is an
estimate that itself is subject to variance.

In this exercise I have used the distribution of t statistics to estimate
p-values. It would also be possible to use the simple difference statistic to
do the same. However, in that case you would need to be more careful about the
variances of the distributions that you sample from. The beauty of the
t-statistic is that it normalises by the amount of variance and so we can use
the same distribution (as long as the sample sizes are the same). 

We can also esimate the p-values directly from the `pt()` function. Please try
this yourself after reading the manual `?pt`.

## Bonus: a vectorised `rowVar` function

In this exercise we have used `apply(x, 1, var)` to calculate the variance for
each row of the matrix `x`. This is slow. There are many ways in which you can
make this faster in `R` by performing vectorised operations.

```R

## for simplicity we will make a matrix x here..
x <- norm.samples[[1]]

## A function to calculate the variance of each row:
rowVars <- function(x){
    m <- rowMeans( x )
    rowSums( (x - m)^2 ) / (ncol(x) - 1)
}

## test it and make sure it works as it should.
system.time(
    x.v1 <- rowVars(x)
)
##    user  system elapsed 
##   0.004   0.000   0.003 

## compare to the apply version..
system.time(
    x.v2 <- apply(x, 1, var)
)
##  user  system elapsed 
## 3.929   0.028   3.956

## so the vectorised function is:
3.956 / 0.003
## about 1300 times faster.

## do we get the same resuls:
plot(x.v1, x.v2, cex=0.1, col=rgb(0,0,0,0.04))

## or:
range(x.v1 - x.v2)
## [1] -7.105427e-15  7.105427e-15
## seems like within reasonable margins of error.
```

My `rowVars` function calculates pretty similar variance estimates, but here
it takes 0.003 seconds instead of almost 4 seconds. In other words it is more
than 1000 times faster. Sometimes doing it the straightforward way in `R` is
too slow and you do need to optimise calculations.


