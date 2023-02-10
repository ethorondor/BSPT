\name{BostonHomicide}
\alias{BostonHomicide}
\title{Youth Homicides in Boston}
\usage{data(BostonHomicide)}
\description{
Data about the number of youth homicides in Boston during the
`Boston Gun Project'---a policing initiative aiming at lowering
homicide victimization among young people in Boston.
}
\format{
A data frame containing 2 monthly time series
and a factor coding seasonality.
 \describe{
    \item{homicides}{time series. Number of youth homicides,}
    \item{population}{time series. Boston population (aged 25-44),
       linearly interpolated from annual data.}
    \item{season}{factor coding the month.}
}
}

\details{The `Boston Gun Project' is a policing initiative aiming at lowering
youth homicides in Boston. The project began in early 1995 and implemented the
so-called `Operation Ceasefire' intervention which began in the late spring of 1996.

More information is available at:

  \url{http://www.ksg.harvard.edu/criminaljustice/research/bgp.htm}
}

\source{Cooper et al. (2001), Figure 1 and Table I.}

\references{
Cooper S.J., Piehl A.M., Braga A.A., Kennedy D.M. (2001), Testing for Structural
Breaks in the Evaluation of Programs, Unpublished Paper, John F. Kennedy School
of Government.

Kennedy D.M., Piehl A.M., Braga A.A. (1996), Youth Violence in Boston: Gun Markets,
Serious Youth Offenders, and a Use-Reduction Strategy, \emph{Law and Contemporary Problems},
\bold{59}, 147-183.
}

\examples{
data(BostonHomicide)

fm <- glm(homicides ~ population + season, data = BostonHomicide,
          family = poisson)
anova(fm, test = "F")
}

\keyword{datasets}

\eof
\name{Fstats}
\alias{Fstats}
\alias{print.Fstats}
\title{F Statistics}
\description{Computes a series of F statistics for a specified data window.}
\usage{
Fstats(formula, from = 0.15, to = NULL, data,
    cov.type = c("const", "HC", "HC1"))}
\arguments{
 \item{formula}{a symbolic description for the model to be tested}
 \item{from, to}{numeric. If \code{from} is smaller than 1 they are
     interpreted as percentages of data and by default \code{to} is taken to be
     1 - \code{from}. F statistics will be calculated for the observations
     \code{(n*from):(n*to)}, when \code{n} is the number of observations in the
     model. If \code{from} is greater than 1 it is interpreted to be the index
     and \code{to} defaults to \code{n - from}. If \code{from} is a vector with
     two elements, then \code{from} and \code{to} are interpreted as time
     specifications like in \code{\link{ts}}, see also the examples.}
 \item{data}{an optional data frame containing the variables in the model. By
     default the variables are taken from the environment which \code{Fstats} is
     called from.}
 \item{cov.type}{a string indicating which type of covariance matrix
     estimator should be used. Constant homoskedastic variances are assumed
     if set to \code{"const"} and White's heteroskedasticity consistent
     estimator is used if set to \code{"HC"}. And \code{"HC1"} stands for a
     standardized estimator of \code{"HC"}, see also \code{\link{covHC}}.}
}
\details{For every potential change point in \code{from:to} a F statistic (Chow
test statistic) is computed. For this an OLS model is fitted for the
observations before and after the potential change point, i.e. \code{2k}
parameters have to be estimated, and the error sum of squares is computed (ESS).
Another OLS model for all obervations with a restricted sum of squares (RSS) is
computed, hence \code{k} parameters have to be estimated here. If \code{n} is 
the number of observations and \code{k} the number of regressors in the model,
the formula is:

\deqn{F = \frac{(RSS - ESS)}{ESS/(n - 2 k)}}{F = (RSS-ESS)/ESS * (n-2*k)}
}

\value{\code{Fstats} returns an object of class \code{"Fstats"}, which contains
mainly a time series of F statistics. The function \code{\link{plot}} has a
method to plot the F statistics or the
corresponding p values; with \code{sctest} a
supF-, aveF- or expF-test on structural change can be performed.}

\references{
Andrews D.W.K. (1993), Tests for parameter instability and structural
change with unknown change point, \emph{Econometrica}, \bold{61}, 821-856.

Hansen B. (1992), Tests for parameter instability in regressions with I(1)
processes, \emph{Journal of Business & Economic Statistics}, \bold{10}, 321-335.

Hansen B. (1997), Approximate asymptotic p values for structural-change
tests, \emph{Journal of Business & Economic Statistics}, \bold{15}, 60-67. }

\seealso{\code{\link{plot.Fstats}}, \code{\link{sctest.Fstats}},
\code{\link{boundary.Fstats}}}

\examples{
require(ts)

## Nile data with one breakpoint: the annual flows drop in 1898
## because the first Ashwan dam was built
data(Nile)
plot(Nile)

## test the null hypothesis that the annual flow remains constant
## over the years
fs.nile <- Fstats(Nile ~ 1)
plot(fs.nile)
sctest(fs.nile)
## visualize the breakpoint implied by the argmax of the F statistics
plot(Nile)
lines(breakpoints(fs.nile))

## UK Seatbelt data: a SARIMA(1,0,0)(1,0,0)_12 model
## (fitted by OLS) is used and reveals (at least) two
## breakpoints - one in 1973 associated with the oil crisis and
## one in 1983 due to the introduction of compulsory
## wearing of seatbelts in the UK.
data(UKDriverDeaths)
seatbelt <- log10(UKDriverDeaths)
seatbelt <- cbind(seatbelt, lag(seatbelt, k = -1), lag(seatbelt, k = -12))
colnames(seatbelt) <- c("y", "ylag1", "ylag12")
seatbelt <- window(seatbelt, start = c(1970, 1), end = c(1984,12))
plot(seatbelt[,"y"], ylab = expression(log[10](casualties)))

## compute F statistics for potential breakpoints between
## 1971(6) (corresponds to from = 0.1) and 1983(6) (corresponds to
## to = 0.9 = 1 - from, the default)
## compute F statistics
fs <- Fstats(y ~ ylag1 + ylag12, data = seatbelt, from = 0.1)
## this gives the same result
fs <- Fstats(y ~ ylag1 + ylag12, data = seatbelt, from = c(1971, 6),
             to = c(1983, 6))
## plot the F statistics
plot(fs, alpha = 0.01)
## plot F statistics with aveF boundary
plot(fs, aveF = TRUE)
## perform the expF test
sctest(fs, type = "expF")
}
\keyword{regression}

\eof
\name{GermanM1}
\alias{GermanM1}
\alias{historyM1}
\alias{monitorM1}
\non_function{}
\title{German M1 Money Demand}
\usage{data(GermanM1)}
\description{
German M1 money demand.
}
\format{
\code{GermanM1} is a data frame containing 12 quarterly time series
from 1961(1) to 1995(4) and two further variables. \code{historyM1}
is the subset of \code{GermanM1} up to 1990(2), i.e., the data before
the German monetary unification on 1990-06-01. \code{monitorM1}
is the complement of \code{historyM1}, i.e., the data after
the unification. All three data frames contain the variables
 \describe{
    \item{m}{time series. Logarithm of real M1 per capita,}
    \item{p}{time series. Logarithm of a price index,}
    \item{y}{time series. Logarithm of real per capita gross
             national product,}
    \item{R}{time series. Long-run interest rate,}
    \item{dm}{time series. First differences of \code{m},}
    \item{dy2}{time series. First differences of lag 2 of \code{y},}
    \item{dR}{time series. First differences of \code{R},}
    \item{dR1}{time series. First differences of lag 1 of \code{R},}
    \item{dp}{time series. First differences of \code{p},}
    \item{m1}{time series. Lag 1 of \code{m},}
    \item{y1}{time series. Lag 1 of \code{y},}
    \item{R1}{time series. Lag 1 of \code{R},}
    \item{season}{factor coding the seasonality,}
    \item{ecm.res}{vector containing the OLS residuals of
                   the Lütkepohl et al. (1999) model fitted in the history
		   period.}
}
}

\details{
Lütkepohl et al. (1999) investigate the linearity and
stability of German M1 money demand: they find a stable regression relation
for the time before the monetary union on 1990-06-01 but a clear structural
instability afterwards.

Zeileis et al. (2002) use a model with
\code{ecm.res} instead of \code{m1}, \code{y1} and \code{R1}, which
leads to equivalent results in the history period but slightly
different results in the monitoring period. The reason for the
replacement is that stationary regressors are needed for the
structural change tests. See references and the examples below for
more details.
}

\source{The data is provided by the German central bank and is
available online in the data archive of the Journal of Applied
Econometrics
\url{http://qed.econ.queensu.ca/jae/1999-v14.5/lutkepohl-terasvirta-wolters/}.}

\references{
Lütkepohl H., Teräsvirta T., Wolters J. (1999), Investigating
Stability and Linearity of a German M1 Money Demand Function,
\emph{Journal of Applied Econometrics}, \bold{14}, 511-525.

Zeileis A., Leisch F., Kleiber C., Hornik K. (2002), Monitoring
Structural Change in Dynamic Econometric Models,
Report 64, SFB "Adaptive Information Systems and Modelling in Economics
and Management Science", Vienna University of Economics,
\url{http://www.wu-wien.ac.at/am/reports.htm#78}.}

\examples{
data(GermanM1)
## Lütkepohl et al. (1999) use the following model
LTW.model <- dm ~ dy2 + dR + dR1 + dp + m1 + y1 + R1 + season
## Zeileis et al. (2002) use
M1.model <- dm ~ dy2 + dR + dR1 + dp + ecm.res + season


## historical tests
ols <- efp(LTW.model, data = GermanM1, type = "OLS-CUSUM")
plot(ols)
re <- efp(LTW.model, data = GermanM1, type = "fluctuation")
plot(re)
fs <- Fstats(LTW.model, data = GermanM1, from = 0.1)
plot(fs)

## monitoring
M1 <- historyM1
ols.efp <- efp(M1.model, type = "OLS-CUSUM", data = M1)
newborder <- function(k) 1.5778*k/118
ols.mefp <- mefp(ols.efp, period = 2)
ols.mefp2 <- mefp(ols.efp, border = newborder)
M1 <- GermanM1
ols.mon <- monitor(ols.mefp)
ols.mon2 <- monitor(ols.mefp2)
plot(ols.mon)
lines(boundary(ols.mon2), col = 2)

## dating
bp <- breakpoints(LTW.model, data = GermanM1)
summary(bp)
plot(bp)

plot(fs)
lines(confint(bp))
}

\keyword{datasets}

\eof
\name{Grossarl}
\alias{Grossarl}
\non_function{}
\title{Marriages and Births in Grossarl}
\usage{data(Grossarl)}
\description{
Data about the number of marriages and illegitimate and legitimate
births in the Austrian Alpine village Grossarl in the 18th and 19th
century.
}
\format{
\code{Grossarl} is a data frame containing 5 annual time series
(1700 - 1899) and 3 factors coding policy interventions.
 \describe{
    \item{marriages}{time series. Number of marriages,}
    \item{illegitimate}{time series. Number of illegitimate births,}
    \item{legitimate}{time series. Number of legitimate births,}
    \item{fraction}{time series. Fraction of illegitimate births,}
    \item{lag.marriages}{time series. Number of marriages in the previous year,}
    \item{politics}{factor coding 4 different political regimes,}
    \item{morals}{factor coding 5 different moral regulations,}
    \item{nuptiality}{factor coding 5 different marriage restrictions.}
}
}

\details{The data frame contains historical demographic data from
Grossarl, a village in the Alpine region of Salzburg, Austria,
during the 18th and 19th century.
During this period, the total population of Grossarl did not vary much on the whole,
with the very exception of the period of the protestant emigrations
in 1731/32.

Especially
during the archbishopric, moral interventions aimed at lowering
the proportion of illegitimate baptisms. For details see the references.}

\source{Parish registers provide the basic
demographic series of baptisms (which is almost equivalent to births
in the study area) and marriages. For
more information see Veichtlbauer et al. (2002).}

\references{
Veichtlbauer O., Hanser E., Zeileis A., Leisch F. (2002),
The Impact Of Policy Interventions on a Pre-Industrial
Population System in the Austrian Alps, forthcoming.

Zeileis A., Veichtlbauer O. (2002), Policy Interventions
Affecting Illegitimacy in Preindustrial Austria:
A Structural Change Analysis, In R. Dutter (ed.),
\emph{Festschrift 50 Jahre Österreichische Statistische Gesellschaft}, 133-146,
Österreichische Statistische Gesellschaft,
\url{http://www.statistik.tuwien.ac.at/oezstat/}.
}

\examples{
data(Grossarl)

## illegitimate births
######################
## lm + MOSUM
plot(Grossarl$fraction)
fm.min <- lm(fraction ~ politics, data = Grossarl)
fm.max <- lm(fraction ~ politics + morals + nuptiality + lag.marriages,
             data = Grossarl)
fm.final <- step(fm.max)
lines(ts(fitted(fm.min), start = 1700), col = 3)
lines(ts(fitted(fm.final), start = 1700), col = 4)
mos.min <- efp(fraction ~ politics, data = Grossarl, type = "OLS-MOSUM")
mos.final <- efp(fraction ~ politics + morals + nuptiality, data = Grossarl,
                 type = "OLS-MOSUM")
plot(mos.min)
lines(mos.final, lty = 2)

## dating
bp <- breakpoints(fraction ~ 1, data = Grossarl, h = 0.1)
summary(bp)
## RSS, BIC, AIC
plot(bp)
plot(0:8, AIC(bp), type = "b")

## probably use 5 (or maybe 6) breakpoints and compare with
## coding of the factors as used by us
##
## politics                   1803      1816 1850
## morals      1736 1753 1771 1803
## nuptiality                 1803 1810 1816      1883
##
## m = 5            1753 1785           1821 1856 1878
## m = 6       1734 1754 1785           1821 1856 1878
##              6    2    5              1    4    3

fm.bp <- lm(fraction ~ breakfactor(breakpoints(bp, breaks = 6)),
            data = Grossarl)

plot(Grossarl$fraction)
lines(fitted(fm.final), col = 3)
lines(fitted(fm.bp), col = 4)



## marriages
############
## lm + MOSUM
plot(Grossarl$marriages)
fm.min <- lm(marriages ~ politics, data = Grossarl)
fm.final <- lm(marriages ~ politics + morals + nuptiality, data = Grossarl)
lines(ts(fitted(fm.min), start = 1700), col = 3)
lines(ts(fitted(fm.final), start = 1700), col = 4)
mos.min <- efp(marriages ~ politics, data = Grossarl, type = "OLS-MOSUM")
mos.final <- efp(marriages ~ politics + morals + nuptiality, data = Grossarl,
                 type = "OLS-MOSUM")
plot(mos.min)
lines(mos.final, lty = 2)

## dating
bp <- breakpoints(marriages ~ 1, data = Grossarl, h = 0.1)
summary(bp)
## RSS, BIC, AIC
plot(bp)
plot(0:8, AIC(bp), type = "b")

## probably use 3 (or maybe 4) breakpoints and compare with
## coding of the factors as used by us
##
## politics                   1803      1816 1850
## morals      1736 1753 1771 1803
## nuptiality                 1803 1810 1816      1883
##
## m = 3       1738                     1813      1875
## m = 4       1738      1794           1814      1875
##              2         4              1         3
fm.bp <- lm(marriages ~ breakfactor(breakpoints(bp, breaks = 4)),
            data = Grossarl)

plot(Grossarl$marriages)
lines(fitted(fm.final), col = 3)
lines(fitted(fm.bp), col = 4)
}

\keyword{datasets}

\eof
\name{USIncExp}
\alias{USIncExp}
\title{Income and Expenditures in the US}
\description{Personal income and personal
consumption expenditures in the US between January 1959 and
February 2001 (seasonally adjusted at annual rates).}
\usage{data(USIncExp)}
\format{
 A multivariate monthly time series from 1959(1) to 2001(2) with variables
 \describe{
    \item{income}{monthly personal income (in billion US dollars),}
    \item{expenditure}{monthly personal consumption expenditures
                       (in billion US Dollars).}
 }
}
\source{\url{http://www.economagic.com/} }
\references{

A. Zeileis, F. Leisch, K. Hornik, C. Kleiber (2002),
strucchange: An R Package for Testing for Structural Change in Linear
Regression Models.
\emph{Journal of Statistical Software} \bold{7}(2), 1--38.
}
\examples{
## These example are presented in the vignette distributed with this
## package, the code was generated by Stangle("strucchange-intro.Rnw")

###################################################
### chunk number 1: data
###################################################
library(strucchange)
data(USIncExp)
plot(USIncExp, plot.type = "single", col = 1:2, ylab = "billion US$")
legend(1960, max(USIncExp), c("income", "expenditures"),
       lty = c(1,1), col = 1:2, bty = "n")


###################################################
### chunk number 2: subset
###################################################
library(strucchange)
data(USIncExp)
library(ts)
USIncExp2 <- window(USIncExp, start = c(1985,12))


###################################################
### chunk number 3: ecm-setup
###################################################
coint.res <- residuals(lm(expenditure ~ income, data = USIncExp2))
coint.res <- lag(ts(coint.res, start = c(1985,12), freq = 12), k = -1)
USIncExp2 <- cbind(USIncExp2, diff(USIncExp2), coint.res)
USIncExp2 <- window(USIncExp2, start = c(1986,1), end = c(2001,2))
colnames(USIncExp2) <- c("income", "expenditure", "diff.income",
                         "diff.expenditure", "coint.res")
ecm.model <- diff.expenditure ~ coint.res + diff.income


###################################################
### chunk number 4: ts-used
###################################################
plot(USIncExp2[,3:5], main = "")


###################################################
### chunk number 5: efp
###################################################
ocus <- efp(ecm.model, type="OLS-CUSUM", data=USIncExp2)
me <- efp(ecm.model, type="ME", data=USIncExp2, h=0.2)


###################################################
### chunk number 6: efp-boundary
###################################################
bound.ocus <- boundary(ocus, alpha=0.05)


###################################################
### chunk number 7: OLS-CUSUM
###################################################
plot(ocus)


###################################################
### chunk number 8: efp-boundary2
###################################################
plot(ocus, boundary = FALSE)
lines(bound.ocus, col = 4)
lines(-bound.ocus, col = 4)


###################################################
### chunk number 9: ME-null
###################################################
plot(me, functional = NULL)


###################################################
### chunk number 10: efp-sctest
###################################################
sctest(ocus)


###################################################
### chunk number 11: efp-sctest2
###################################################
sctest(ecm.model, type="OLS-CUSUM", data=USIncExp2)


###################################################
### chunk number 12: Fstats
###################################################
fs <- Fstats(ecm.model, from = c(1990, 1), to = c(1999,6), data = USIncExp2)


###################################################
### chunk number 13: Fstats-plot
###################################################
plot(fs)


###################################################
### chunk number 14: pval-plot
###################################################
plot(fs, pval=TRUE)


###################################################
### chunk number 15: aveF-plot
###################################################
plot(fs, aveF=TRUE)


###################################################
### chunk number 16: Fstats-sctest
###################################################
sctest(fs, type="expF")


###################################################
### chunk number 17: Fstats-sctest2
###################################################
sctest(ecm.model, type = "expF", from = 49, to = 162, data = USIncExp2)


###################################################
### chunk number 18: mefp
###################################################
USIncExp3 <- window(USIncExp2, start = c(1986, 1), end = c(1989,12))
me.mefp <- mefp(ecm.model, type = "ME", data = USIncExp3, alpha = 0.05)


###################################################
### chunk number 19: monitor1
###################################################
USIncExp3 <- window(USIncExp2, start = c(1986, 1), end = c(1990,12))
me.mefp <- monitor(me.mefp)


###################################################
### chunk number 20: monitor2
###################################################
USIncExp3 <- window(USIncExp2, start = c(1986, 1))
me.mefp <- monitor(me.mefp)
me.mefp


###################################################
### chunk number 21: monitor-plot
###################################################
plot(me.mefp)


###################################################
### chunk number 22: mefp2
###################################################
USIncExp3 <- window(USIncExp2, start = c(1986, 1), end = c(1989,12))
me.efp <- efp(ecm.model, type = "ME", data = USIncExp3, h = 0.5)
me.mefp <- mefp(me.efp, alpha=0.05)


###################################################
### chunk number 23: monitor3
###################################################
USIncExp3 <- window(USIncExp2, start = c(1986, 1))
me.mefp <- monitor(me.mefp)


###################################################
### chunk number 24: monitor-plot2
###################################################
plot(me.mefp)

}

\keyword{datasets}

\eof
\name{boundary.Fstats}
\alias{boundary.Fstats}
\title{Boundary for F Statistics}
\description{Computes boundary for an object of class \code{"Fstats"}}
\usage{
\method{boundary}{Fstats}(x, alpha = 0.05, pval = FALSE, aveF = FALSE,
    asymptotic = FALSE, ...)}
\arguments{
 \item{x}{an object of class \code{"Fstats"}.}
 \item{alpha}{numeric from interval (0,1) indicating the confidence level for
     which the boundary of the supF test will be computed.}
 \item{pval}{logical. If set to \code{TRUE} a boundary for the corresponding p
     values will be computed.}
 \item{aveF}{logical. If set to \code{TRUE} the boundary of the aveF (instead of
    the supF) test will be computed. The resulting boundary then is a boundary
    for the mean of the F statistics rather than for the F statistics
    themselves.}
 \item{asymptotic}{logical. If set to \code{TRUE} the asymptotic (chi-square)
     distribution instead of the exact (F) distribution will be used to compute
     the p values (only if \code{pval} is \code{TRUE}).}
 \item{...}{currently not used.}
}

\value{an object of class \code{"ts"} with the same time properties as
the time series in \code{x}}

\seealso{\code{\link{Fstats}}, \code{\link{plot.Fstats}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains
## constant over the years for potential break points between 1941
## (corresponds to from = 0.5) and 1962 (corresponds to to = 0.85)
## compute F statistics
fs <- Fstats(nhtemp ~ 1, from = 0.5, to = 0.85)
## plot the p values without boundary
plot(fs, pval = TRUE, alpha = 0.01)
## add the boundary in another colour
lines(boundary(fs, pval = TRUE, alpha = 0.01), col = 2)
}
\keyword{regression}

\eof
\name{boundary}
\alias{boundary}
\title{Boundary Function for Structural Change Tests}
\description{A generic function computing boundaries for structural change
tests}
\usage{
boundary(x, ...)}

\arguments{
  \item{x}{an object. Use \code{\link{methods}} to see which
    \code{\link{class}} has a method for boundary.}
  \item{...}{additional arguments affecting the boundary.}
}

\value{an object of class \code{"ts"} with the same time properties as
the time series in \code{x}}

\seealso{\code{\link{boundary.efp}}, \code{\link{boundary.mefp}},
  \code{\link{boundary.Fstats}}}

\keyword{regression}

\eof
\name{boundary.efp}
\alias{boundary.efp}
\title{Boundary for Empirical Fluctuation Processes}
\description{Computes boundary for an object of class \code{"efp"}}
\usage{
\method{boundary}{efp}(x, alpha = 0.05, alt.boundary = FALSE,
   functional = "max", ...)
}

\arguments{
  \item{x}{an object of class \code{"efp"}.}
  \item{alpha}{numeric from interval (0,1) indicating the confidence level for
     which the boundary of the corresponding test will be computed.}
  \item{alt.boundary}{logical. If set to \code{TRUE} alternative boundaries
     (instead of the standard linear boundaries) will be computed (for Brownian
     bridge type processes only).}
  \item{functional}{indicates which functional should be applied to the
     empirical fluctuation process. See also \code{\link{plot.efp}}.}
  \item{\dots}{currently not used.}
}

\value{an object of class \code{"ts"} with the same time properties as
the process in \code{x}}

\seealso{\code{\link{efp}}, \code{\link{plot.efp}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains constant
## over the years
## compute OLS-CUSUM fluctuation process
temp.cus <- efp(nhtemp ~ 1, type = "OLS-CUSUM")
## plot the process without boundaries
plot(temp.cus, alpha = 0.01, boundary = FALSE)
## add the boundaries in another colour
bound <- boundary(temp.cus, alpha = 0.01)
lines(bound, col=4)
lines(-bound, col=4)
}
\keyword{regression}

\eof
\name{boundary.mefp}
\alias{boundary.mefp}
\title{Boundary Function for Monitoring of Structural Changes}
\description{Computes boundary for an object of class \code{"mefp"}}
\usage{
\method{boundary}{mefp}(x, ...)}
\arguments{
  \item{x}{an object of class \code{"mefp"}.}
  \item{...}{currently not used.}
}
\value{an object of class \code{"ts"} with the same time properties as
the monitored process}
\seealso{\code{\link{mefp}}, \code{\link{plot.mefp}}}
\examples{
df1 <- data.frame(y=rnorm(300))
df1[150:300,"y"] <- df1[150:300,"y"]+1
me1 <- mefp(y~1, data=df1[1:50,,drop=FALSE], type="ME", h=1,
              alpha=0.05)
me2 <- monitor(me1, data=df1)

plot(me2, boundary=FALSE)
lines(boundary(me2), col="green", lty="44")
}
\keyword{regression}


\eof
\name{breakdates}
\alias{breakdates}
\alias{breakdates.breakpoints}
\alias{breakdates.confint.breakpoints}
\title{Breakdates Corresponding to Breakpoints}
\description{
 A generic function for computing the breakdates corresponding
 to breakpoints (and their confidence intervals).
}
\usage{
breakdates(obj, format.times = FALSE, ...)
}
\arguments{
  \item{obj}{An object of class \code{"breakpoints"}, \code{"breakpointsfull"} or their
             confidence intervals as returned by \code{\link{confint}}.}
  \item{format.times}{logical. If set to \code{TRUE} a vector of
    strings with the formatted breakdates. See details for more
    information.}
  \item{\dots}{currently not used.}
}
\details{
  Breakpoints are the number of observations that are the last in one
  segment and breakdates are the corresponding points on the underlying
  time scale. The breakdates can be formatted which enhances readability
  in particular for quarterly or monthly time series. For example the
  breakdate \code{2002.75} of a monthly time series will be formatted to
  \code{"2002(10)"}.
}

\value{
  A vector or matrix containing the breakdates.
}

\seealso{\code{\link{breakpoints}}, \code{\link{confint}}}

\examples{
require(ts)

## Nile data with one breakpoint: the annual flows drop in 1898
## because the first Ashwan dam was built
data(Nile)
plot(Nile)

bp.nile <- breakpoints(Nile ~ 1)
summary(bp.nile)
plot(bp.nile)

## compute breakdates corresponding to the
## breakpoints of minimum BIC segmentation
breakdates(bp.nile)

## confidence intervals
ci.nile <- confint(bp.nile)
breakdates(ci.nile)
ci.nile

plot(Nile)
lines(ci.nile)
}

\keyword{regression}

\eof
\name{breakfactor}
\alias{breakfactor}
\title{Factor Coding of Segmentations}
\description{
 Generates a factor encoding the segmentation given by
 a set of breakpoints.
}
\usage{
breakfactor(obj, breaks = NULL, labels = NULL, ...)
}
\arguments{
  \item{obj}{An object of class \code{"breakpoints"} or
    \code{"breakpointsfull"} respectively.}
  \item{breaks}{an integer specifying the number of breaks
    to extract (only if \code{obj} is of class \code{"breakpointsfull"}),
    by default the minimum BIC partition is used.}
  \item{labels}{a vector of labels for the returned factor,
    by default the segments are numbered starting from
    \code{"segment1"}.}
  \item{\dots}{further arguments passed to \code{factor}.}
}

\value{
  A factor encoding the segmentation.
}

\seealso{\code{\link{breakpoints}}}

\examples{
require(ts)

## Nile data with one breakpoint: the annual flows drop in 1898
## because the first Ashwan dam was built
data(Nile)
plot(Nile)

## compute breakpoints
bp.nile <- breakpoints(Nile ~ 1)

## fit and visualize segmented and unsegmented model
fm0 <- lm(Nile ~ 1)
fm1 <- lm(Nile ~ breakfactor(bp.nile, breaks = 1))

lines(fitted(fm0), col = 3)
lines(fitted(fm1), col = 4)
lines(bp.nile, breaks = 1)
}

\keyword{regression}

\eof
\name{breakpoints}
\alias{breakpoints}
\alias{breakpoints.formula}
\alias{breakpoints.breakpointsfull}
\alias{breakpoints.Fstats}
\alias{summary.breakpoints}
\alias{summary.breakpointsfull}
\alias{plot.breakpointsfull}
\alias{plot.summary.breakpointsfull}
\alias{print.breakpoints}
\alias{print.summary.breakpointsfull}
\alias{lines.breakpoints}
\title{Dating Breaks}

\description{
  Computation of breakpoints in regression relationships. Given a number
  of breaks the function computes the optimal breakpoints.
}

\usage{
\method{breakpoints}{formula}(formula, h = 0.15, breaks = NULL,
    data = list(), ...)
\method{breakpoints}{breakpointsfull}(obj, breaks = NULL, ...)
\method{summary}{breakpointsfull}(object, breaks = NULL, sort = TRUE,
    format.times = NULL, ...)
\method{lines}{breakpoints}(x, breaks = NULL, lty = 2, ...)
}

\arguments{
  \item{formula}{a symbolic description for the model in which breakpoints
    will be estimated.}
  \item{h}{minimal segment size either given as fraction relative to the
    sample size or as an integer giving the minimal number of observations
    in each segment.}
  \item{breaks}{integer specifying the maximal number of breaks to be calculated.
    By default the maximal number allowed by \code{h} is used.}
  \item{data}{an optional data frame containing the variables in the model. By
     default the variables are taken from the environment which \code{breakpoints} is
     called from.}
  \item{\dots}{currently not used.}
  \item{obj, object}{an object of class \code{"breakpointsfull"}.}
  \item{sort}{logical. If set to \code{TRUE} \code{summary} tries to match
    the breakpoints from partitions with different numbers of breaks.}
  \item{format.times}{logical. If set to \code{TRUE} a vector of
    strings with the formatted breakdates is printed. See \code{\link{breakdates}}
    for more information.}
  \item{x}{an object of class \code{"breakpoints"}.}
  \item{lty}{line type.}
}

\details{
  All procedures in this package are concerned with testing or assessing
  deviations from stability in the classical linear regression model

  \deqn{y_i = x_i^\top \beta + u_i}{y_i = x_i' b + u_i}

  In many applications it is reasonable to assume
  that there are \eqn{m} breakpoints, where the coefficients shift from
  one stable regression relationship to a different one. Thus,
  there are \eqn{m+1} segments in which the regression coefficients are
  constant, and the model can be rewritten as

  \deqn{y_i = x_i^\top \beta_j + u_i
  \qquad (i = i_{j-1} + 1, \dots, i_j, \quad j = 1, \dots, m+1)}{y_i =
    x_i' b_j + u_i   (i = i_{j-1} + 1, \dots, i_j,   j = 1, \dots, m+1)}

  where \eqn{j} denotes the segment index. In practice the breakpoints \eqn{i_j}
  are rarely given exogenously, but have to be estimated.
  \code{breakpoints} estimates these breakpoints by minimizing the residual sum of
  squares (RSS) of the equation above.

  The foundation for estimating breaks in time series regression models
  was given by Bai (1994) and was extended to multiple breaks by Bai (1997ab)
  and Bai & Perron (1998). \code{breakpoints} implements the algorithm
  described in Bai & Perron (2003) for simultanous estimation of
  multiple breakpoints. The distribution function used for the confidence
  intervals for the breakpoints is given in Bai (1997b). The ideas behind
  this implementation are described in Zeileis et al. (2003).

  The algorithm for computing the optimal breakpoints given the number
  of breaks is based on a dynamic programming approach. The underlying
  idea is that of the Bellman principle. The main computational effort
  is to compute a triangular RSS matrix, which gives the residual
  sum of squares for a segment starting at observation \eqn{i} and
  ending at \eqn{i'} with \eqn{i} < \eqn{i'}.

  Given a \code{formula} as the first argument, \code{breakpoints} computes
  an object of class \code{"breakpointsfull"} which inherits from \code{"breakpoints"}.
  This contains in particular the triangular RSS
  matrix and functions to extract an optimal segmentation. A \code{summary}
  of this object will give the breakpoints (and associated) breakdates
  for all segmentations up to the maximal number of breaks together
  with the associated RSS and BIC. These will be plotted if \code{plot}
  is applied and thus visualize the minimum BIC estimator of the number
  of breakpoints. From an object of class \code{"breakpointsfull"} an
  arbitrary number of \code{breaks} (admissable by the minimum segment
  size \code{h}) can be extracted by another application of
  \code{breakpoints}, returning an object of class \code{"breakpoints"}.
  This contains only the breakpoints for the specified number of breaks
  and some model properties (number of observations, regressors, time
  series properties and the associated RSS) but not the triangular RSS
  matrix and related extractor functions. The set of breakpoints which
  is associated by default with a \code{"breakpointsfull"} object is
  the minimum BIC partition.

  Breakpoints are the number of observations that are the last in one
  segment, it is also possible to compute the corresponding \code{breakdates}
  which are the breakpoints on the underlying time scale. The breakdates
  can be formatted which enhances readability in particular for quarterly
  or monthly time series. For example the breakdate \code{2002.75} of a monthly
  time series will be formatted to \code{"2002(10)"}. See \code{\link{breakdates}}
  for more details.

  From a \code{"breakpointsfull"} object confidence intervals for the breakpoints
  can be computed using the method of \code{\link{confint}}.
  The breakdates corresponding to the breakpoints can again be computed
  by \code{\link{breakdates}}. The breakpoints and their confidence
  intervals can be visualized by \code{lines}.

  The log likelihood as well as some information criteria can be computed
  using the methods for the \code{\link{logLik}} and \code{\link{AIC}}. As
  for linear models the log likelihood is computed on a normal model and
  the degrees of freedom are the number of regression coefficients multiplied
  by the number of segements plus the number of estimated breakpoints plus
  1 for the error variance. More details can be found on the help page of
  the method \code{\link{logLik.breakpoints}}.

  As the maximum of a sequence of F statistics is equivalent to the minimum
  OLS estimator of the breakpoint in a 2-segment partition it can be
  extracted by \code{breakpoints} from an object of class \code{"Fstats"}
  as computed by \code{\link{Fstats}}. However, this cannot be used to extract
  a larger number of breakpoints.

  For illustration see the commented examples below and Zeileis et al. (2003).
}
\section{value}{
  An object of class \code{"breakpoints"} is a list with the following
  elements:
  \describe{
   \item{breakpoints}{the breakpoints of the optimal partition with the
     number of breaks specified,}
   \item{RSS}{the associated RSS,}
   \item{nobs}{the number of observations,}
   \item{nreg}{the number of regressors,}
   \item{call}{the function call,}
   \item{datatsp}{the time series properties \code{tsp} of the data,
     if any, \code{c(0, 1, nobs)} otherwise.}
  }

  If applied to a \code{formula} as first argument, \code{breakpoints} returns an object of class
  \code{"breakpointsfull"} (which inherits from \code{"breakpoints"}), that
  contains some additional (or slightly different) elements such as:
  \describe{
   \item{breakpoints}{the breakpoints of the minimum BIC partition,}
   \item{RSS}{a function which takes two arguments \code{i,j} and computes
     the residual sum of squares for a segment starting at observation \code{i} and
     ending at \code{j} by looking up the corresponding element in the triangular
     RSS matrix \code{RSS.triang},}
   \item{RSS.triang}{a list encoding the triangular RSS matrix.}
  }
}

\references{
Bai J. (1994), Least Squares Estimation of a Shift in Linear Processes,
  \emph{Journal of Time Series Analysis}, \bold{15}, 453-472.

Bai J. (1997a), Estimating Multiple Breaks One at a Time,
  \emph{Econometric Theory}, \bold{13}, 315-352.

Bai J. (1997b), Estimation of a Change Point in Multiple Regression Models,
  \emph{Review of Economics and Statistics}, \bold{79}, 551-563.

Bai J., Perron P. (1998), Estimating and Testing Linear Models With Multiple Structural
   Changes, \emph{Econometrica}, \bold{66}, 47-78.

Bai J., Perron P. (2003), Computation and Analysis of Multiple Structural Change
  Models, \emph{Journal of Applied Econometrics}, \bold{18}, 1-22.

Zeileis A., Kleiber C., Krämer W., Hornik K. (2003), Testing and Dating of
  Structural Changes in Practice, \emph{Computational Statistics and Data Analysis},
  forthcoming.
}

\examples{
require(ts)

## Nile data with one breakpoint: the annual flows drop in 1898
## because the first Ashwan dam was built
data(Nile)
plot(Nile)

## F statistics indicate one breakpoint
fs.nile <- Fstats(Nile ~ 1)
plot(fs.nile)
breakpoints(fs.nile)
lines(breakpoints(fs.nile))

## or
bp.nile <- breakpoints(Nile ~ 1)
summary(bp.nile)

## the BIC also chooses one breakpoint
plot(bp.nile)
breakpoints(bp.nile)

## fit null hypothesis model and model with 1 breakpoint
fm0 <- lm(Nile ~ 1)
fm1 <- lm(Nile ~ breakfactor(bp.nile, breaks = 1))
plot(Nile)
lines(fitted(fm0), col = 3)
lines(fitted(fm1), col = 4)
lines(bp.nile)

## confidence interval
ci.nile <- confint(bp.nile)
ci.nile
lines(ci.nile)


## UK Seatbelt data: a SARIMA(1,0,0)(1,0,0)_12 model
## (fitted by OLS) is used and reveals (at least) two
## breakpoints - one in 1973 associated with the oil crisis and
## one in 1983 due to the introduction of compulsory
## wearing of seatbelts in the UK.
data(UKDriverDeaths)
seatbelt <- log10(UKDriverDeaths)
seatbelt <- cbind(seatbelt, lag(seatbelt, k = -1), lag(seatbelt, k = -12))
colnames(seatbelt) <- c("y", "ylag1", "ylag12")
seatbelt <- window(seatbelt, start = c(1970, 1), end = c(1984,12))
plot(seatbelt[,"y"], ylab = expression(log[10](casualties)))

## testing
re.seat <- efp(y ~ ylag1 + ylag12, data = seatbelt, type = "RE")
plot(re.seat)

## dating
bp.seat <- breakpoints(y ~ ylag1 + ylag12, data = seatbelt, h = 0.1)
summary(bp.seat)
lines(bp.seat, breaks = 2)

## minimum BIC partition
plot(bp.seat)
breakpoints(bp.seat)
## the BIC would choose 0 breakpoints although the RE and supF test
## clearly reject the hypothesis of structural stability. Bai &
## Perron (2003) report that the BIC has problems in dynamic regressions.
## due to the shape of the RE process of the F statistics choose two
## breakpoints and fit corresponding models
bp.seat2 <- breakpoints(bp.seat, breaks = 2)
fm0 <- lm(y ~ ylag1 + ylag12, data = seatbelt)
fm1 <- lm(y ~ breakfactor(bp.seat2)/(ylag1 + ylag12) - 1, data = seatbelt)

## plot
plot(seatbelt[,"y"], ylab = expression(log[10](casualties)))
time.seat <- as.vector(time(seatbelt))
lines(time.seat, fitted(fm0), col = 3)
lines(time.seat, fitted(fm1), col = 4)
lines(bp.seat2)

## confidence intervals
ci.seat2 <- confint(bp.seat, breaks = 2)
ci.seat2
lines(ci.seat2)
}

\keyword{regression}

\eof
\name{confint.breakpointsfull}
\alias{confint.breakpointsfull}
\alias{lines.confint.breakpoints}
\alias{print.confint.breakpoints}
\title{Confidence Intervals for Breakpoints}
\description{
  Computes confidence intervals for breakpoints.
}
\usage{
\method{confint}{breakpointsfull}(object, parm = NULL, level = 0.95,
    breaks = NULL, het.reg = TRUE, het.err = TRUE, ...)
\method{lines}{confint.breakpoints}(x, col = 2, angle = 90, length = 0.05,
    code = 3, at = NULL, breakpoints = TRUE, ...)
}
\arguments{
  \item{object}{an object of class \code{"breakpointsfull"} as computed by
    \code{\link{breakpoints}} from a \code{formula}.}
  \item{parm}{the same as \code{breaks}, only one of the two should be
    specified.}
  \item{level}{the confidence level required.}
  \item{breaks}{an integer specifying the number of breaks to be used.
    By default the breaks of the minimum BIC partition are used.}
  \item{het.reg}{logical. Should heterogenous regressors be assumed? If set
    to \code{FALSE} the distribution of the regressors is assumed to be
    homogenous over the segments.}
  \item{het.err}{logical. Should heterogenous errors be assumed? If set
    to \code{FALSE} the distribution of the errors is assumed to be
    homogenous over the segments.}
  \item{x}{an object of class \code{"confint.breakpoints"} as returned by
    \code{confint}.}
  \item{col, angle, length, code}{arguments passed to \code{\link{arrows}}.}
  \item{at}{position on the y axis, where the confidence arrows should be
    drawn. By default they are drawn at the bottom of the plot.}
  \item{breakpoints}{logical. If \code{TRUE} vertical lines for the breakpoints
    are drawn.}
  \item{\dots}{\emph{currently not used}.}
}

\details{
  As the breakpoints are integers (observation numbers) the corresponding
  confidence intervals are also rounded to integers.

  The distribution function used for the computation of confidence
  intervals of breakpoints is given in Bai (1997). The procedure, in
  particular the usage of heterogenous regressors and/or errors, is
  described in more detail in Bai & Perron (2003).

  The breakpoints should be computed from a formula with \code{breakpoints},
  then the confidence intervals for the breakpoints can be derived by
  \code{confint} and these can be visualized by \code{lines}. For an
  example see below.
}

\value{
  A matrix containing the breakpoints and their lower and upper
  confidence boundary for the given level.
}

\references{
Bai J. (1997), Estimation of a Change Point in Multiple Regression Models,
  \emph{Review of Economics and Statistics}, \bold{79}, 551-563.

Bai J., Perron P. (2003), Computation and Analysis of Multiple Structural Change
  Models, \emph{Journal of Applied Econometrics}, \bold{18}, 1-22.
}

\seealso{\code{\link{breakpoints}}}

\examples{
require(ts)

## Nile data with one breakpoint: the annual flows drop in 1898
## because the first Ashwan dam was built
data(Nile)
plot(Nile)

## dating breaks
bp.nile <- breakpoints(Nile ~ 1)
ci.nile <- confint(bp.nile, breaks = 1)
lines(ci.nile)
}

\keyword{regression}

\eof
\name{covHC}
\alias{covHC}
\title{Heteroskedasticity-Consistent Covariance Matrix Estimation}
\description{
 Heteroskedasticity-consistent estimation of the covariance matrix of the
 coefficient estimates in a linear regression model.
}
\usage{
covHC(formula, type = c("HC2", "const", "HC", "HC1", "HC3"),
   data=list())
}
\arguments{
 \item{formula}{a symbolic description for the model to be fitted.}
  \item{type}{a character string specifying the estimation type. For
     details see below.}
 \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment which \code{covHC}
    is called from.}
}
\details{When \code{type = "const"} constant variances are assumed and
  and \code{covHC} gives the usual estimate of the covariance matrix of
  the coefficient estimates:

  \deqn{\hat \sigma^2 (X^\top X)^{-1}}{sigma^2 (X'X)^{-1}}

  All other methods do not assume constant variances and are suitable in case of
  heteroskedasticity. \code{"HC"} gives White's estimator; for details see the
  references.
}
\value{A matrix containing the covariance matrix estimate.
}
\references{

MacKinnon J. G., White H. (1985),
Some heteroskedasticity-consistent
covariance matrix estimators with improved finite sample properties.
\emph{Journal of Econometrics} \bold{29}, 305-325

}

\seealso{\code{\link{lm}}}

\examples{
## generate linear regression relationship
## with homoskedastic variances
x <- sin(1:100)
y <- 1 + x + rnorm(100)
## compute usual covariance matrix of coefficient estimates
covHC(y~x, type="const")

sigma2 <- sum(residuals(lm(y~x))^2)/98
sigma2 * solve.crossprod(cbind(1,x))
}
\keyword{htest}

\eof
\name{durab}
\alias{durab}
\non_function{}
\title{US Labor Productivity}
\usage{data(durab)}
\description{
US labor productivity in the manufacturing/durables sector.
}
\format{
\code{durab} is a multivariate monthly time series from 1947(3)
to 2001(4) with variables
 \describe{
    \item{y}{growth rate of the Industrial Production Index to
             average weekly labor hours in the manufacturing/durables sector,}
    \item{lag}{lag 1 of the series \code{y},}
}
}

\source{The data set is available from Bruce Hansen's homepage
\url{http://www.ssc.wisc.edu/~bhansen/}. For more information
see Hansen (2001).}

\references{
Hansen B. (2001), The New Econometrics of Structural Change:
Dating Breaks in U.S. Labor Productivity,
\emph{Journal of Economic Perspectives}, \bold{15}, 117-128.

Zeileis A., Leisch F., Kleiber C., Hornik K. (2002), Monitoring
Structural Change in Dynamic Econometric Models,
Report 64, SFB "Adaptive Information Systems and Modelling in Economics
and Management Science", Vienna University of Economics,
\url{http://www.wu-wien.ac.at/am/reports.htm#78}.}

\examples{
data(durab)
## use AR(1) model as in Hansen (2001) and Zeileis et al. (2002)
durab.model <- y ~ lag

## historical tests
## OLS-based CUSUM process
ols <- efp(durab.model, data = durab, type = "OLS-CUSUM")
plot(ols)
## F statistics
fs <- Fstats(durab.model, data = durab, from = 0.1)
plot(fs)
## F statistics based on heteroskadisticy-consistent covariance matrix
fsHC <- Fstats(durab.model, data = durab, from = 0.1, cov.type = "HC")
plot(fsHC)

## monitoring
Durab <- window(durab, start=1964, end = c(1979, 12))
ols.efp <- efp(durab.model, type = "OLS-CUSUM", data = Durab)
newborder <- function(k) 1.5778*k/192
ols.mefp <- mefp(ols.efp, period=2)
ols.mefp2 <- mefp(ols.efp, border=newborder)
Durab <- window(durab, start=1964)
ols.mon <- monitor(ols.mefp)
ols.mon2 <- monitor(ols.mefp2)
plot(ols.mon)
lines(boundary(ols.mon2), col = 2)

## dating
bp <- breakpoints(durab.model, data = durab)
summary(bp)
plot(summary(bp))

plot(ols)
lines(breakpoints(bp, breaks = 1), col = 3)
lines(breakpoints(bp, breaks = 2), col = 4)
plot(fs)
lines(breakpoints(bp, breaks = 1), col = 3)
lines(breakpoints(bp, breaks = 2), col = 4)
}
\keyword{datasets}

\eof
\name{efp}
\alias{efp}
\alias{print.efp}
\title{Empirical Fluctuation Processes}
\description{Computes an empirical fluctuation process according
to a specified method from the generalized fluctuation test
framework, which includes CUSUM and MOSUM tests based on recursive
or OLS residuals, parameter estimates or ML scores (OLS first order
conditions).}
\usage{
efp(formula, data, type = <<see below>>, h = 0.15,
    dynamic = FALSE, rescale = TRUE)}
\arguments{
  \item{formula}{a symbolic description for the model to be tested.}
  \item{data}{an optional data frame containing the variables in the model. By
     default the variables are taken from the environment which \code{efp} is
     called from.}
  \item{type}{specifies which type of fluctuation process will be
     computed. For details see below.}
  \item{h}{a numeric from interval (0,1) sepcifying the bandwidth. determins the
     size of the data window relative to sample size (for MOSUM and ME processes
     only).}
  \item{dynamic}{logical. If \code{TRUE} the lagged observations are included as
     a regressor.}
  \item{rescale}{logical. If \code{TRUE} the estimates will be standardized by
    the regressor matrix of the corresponding subsample according to Kuan & Chen
    (1994); if \code{FALSE} the whole regressor matrix will be used.
    (only if \code{type} is either \code{"RE"} or \code{"ME"})}
}

\details{If \code{type} is one of \code{"Rec-CUSUM"}, \code{"OLS-CUSUM"},
\code{"Rec-MOSUM"} or \code{"OLS-MOSUM"} the function \code{efp} will return a
one-dimensional empiricial process of sums of residuals. Either it will be based
on recursive residuals or on OLS residuals and the process will contain
CUmulative SUMs or MOving SUMs of residuals in a certain data window.
For the MOSUM and ME processes all estimations are done for the
observations in a moving data window, whose size is determined by \code{h} and
which is shifted over the whole sample.

If \code{type} is either \code{"RE"} or \code{"ME"} a
\emph{k}-dimensional process will be returned, if \emph{k} is the number of
regressors in the model, as it is based on recursive OLS estimates of the
regression coefficients or moving OLS estimates respectively. The recursive
estimates test is also called fluctuation test, therefore setting \code{type}
to \code{"fluctuation"} was used to specify it in earlier versions of
strucchange. It still can be used now, but will be forced to \code{"RE"}.

If \code{type} is \code{"Score-CUSUM"} or \code{"Score-MOSUM"} a \emph{k+1}-dimensional
process will be returned, one for each score of the regression coefficients and one for
the scores of the variance. The process gives the decorrelated cumulative sums of the ML
scores (in a gaussian model) or first order conditions respectively (in an OLS framework).

If there is a single structural change point \eqn{t^*}, the recursive CUSUM path
starts to depart from its mean 0 at \eqn{t^*}. The Brownian bridge type paths
will have their respective peaks around \eqn{t^*}.
The Brownian bridge increments type paths should have a strong change at \eqn{t^*}.

The function \code{\link{plot}}
has a method to plot the empirical fluctuation process; with
\code{sctest} the corresponding test on structural change can be
performed.
}

\value{
  \code{efp} returns a list of class \code{"efp"} with components inlcuding
  \itemize
  \item{process}{the fitted empirical fluctuation process of class
    \code{"ts"} or \code{"mts"} respectively,}
  \item{type}{a string with the \code{type} of the process fitted,}
  \item{nreg}{the number of regressors,}
  \item{nobs}{the number of observations,}
  \item{par}{the bandwidth \code{h} used.}
}

\references{Brown R.L., Durbin J., Evans J.M. (1975), Techniques for
testing constancy of regression relationships over time, \emph{Journal of the
Royal Statistal Society}, B, \bold{37}, 149-163.

Chu C.-S., Hornik K., Kuan C.-M. (1995), MOSUM tests for parameter
constancy, \emph{Biometrika}, \bold{82}, 603-617.

Chu C.-S., Hornik K., Kuan C.-M. (1995), The moving-estimates test for
parameter stability, \emph{Econometric Theory}, \bold{11}, 669-720.

Hansen B. (1992), Testing for Parameter Instability in Linear Models,
\emph{Journal of Policy Modeling}, \bold{14}, 517-533.

Hjort N.L., Koning A. (2002), Tests for Constancy of Model Parameters
Over Time, \emph{Nonparametric Statistics}, \bold{14}, 113-132.

Krämer W., Ploberger W., Alt R. (1988), Testing for structural change in
dynamic models, \emph{Econometrica}, \bold{56}, 1355-1369.

Kuan C.-M., Hornik K. (1995), The generalized fluctuation test: A
unifying view, \emph{Econometric Reviews}, \bold{14}, 135 - 161.

Kuan C.-M., Chen (1994), Implementing the fluctuation and moving estimates
tests in dynamic econometric models, \emph{Economics Letters}, \bold{44},
235-239.

Ploberger W., Krämer W. (1992), The CUSUM test with OLS residuals,
\emph{Econometrica}, \bold{60}, 271-285.

Zeileis A., Leisch F., Hornik K., Kleiber C. (2002), \code{strucchange}:
An R Package for Testing for Structural Change in Linear Regression Models,
\emph{Journal of Statistical Software}, \bold{7}(2), 1-38.

Zeileis A., Hornik K. (2003), Generalized M-Fluctuation Tests for Parameter
Instability, Report 80, SFB "Adaptive Information Systems and Modelling in Economics
and Management Science", Vienna University of Economics,
\url{http://www.wu-wien.ac.at/am/reports.htm#80}.
}

\seealso{\code{\link{plot.efp}}, \code{\link{print.efp}},
\code{\link{sctest.efp}}, \code{\link{boundary.efp}}}

\examples{
require(ts)

## Nile data with one breakpoint: the annual flows drop in 1898
## because the first Ashwan dam was built
data(Nile)
plot(Nile)

## test the null hypothesis that the annual flow remains constant
## over the years
## compute OLS-based CUSUM process and plot
## with standard and alternative boundaries
ocus.nile <- efp(Nile ~ 1, type = "OLS-CUSUM")
plot(ocus.nile)
plot(ocus.nile, alpha = 0.01, alt.boundary = TRUE)
## calculate corresponding test statistic
sctest(ocus.nile)

## UK Seatbelt data: a SARIMA(1,0,0)(1,0,0)_12 model
## (fitted by OLS) is used and reveals (at least) two
## breakpoints - one in 1973 associated with the oil crisis and
## one in 1983 due to the introduction of compulsory
## wearing of seatbelts in the UK.
data(UKDriverDeaths)
seatbelt <- log10(UKDriverDeaths)
seatbelt <- cbind(seatbelt, lag(seatbelt, k = -1), lag(seatbelt, k = -12))
colnames(seatbelt) <- c("y", "ylag1", "ylag12")
seatbelt <- window(seatbelt, start = c(1970, 1), end = c(1984,12))
plot(seatbelt[,"y"], ylab = expression(log[10](casualties)))

## use RE process
re.seat <- efp(y ~ ylag1 + ylag12, data = seatbelt, type = "RE")
plot(re.seat)
plot(re.seat, functional = NULL)
sctest(re.seat)
}
\keyword{regression}

\eof
\name{logLik.breakpoints}
\alias{logLik.breakpoints}
\alias{logLik.breakpointsfull}
\alias{AIC.breakpointsfull}
\title{Log Likelihood and Information Criteria for Breakpoints}

\description{
  Computation of log likelihood and AIC type information criteria
  for partitions given by breakpoints.
}

\usage{
\method{logLik}{breakpointsfull}(object, breaks = NULL, ...)
\method{AIC}{breakpointsfull}(object, breaks = NULL, ..., k = 2)
}

\arguments{
  \item{object}{an object of class \code{"breakpoints"} or \code{"breakpointsfull"}.}
  \item{breaks}{if \code{object} is of class \code{"breakpointsfull"} the
    number of breaks can be specified.}
  \item{\dots}{\emph{currently not used}.}
  \item{k}{the penalty parameter to be used, the default \code{k = 2}
    is the classical AIC, \code{k = log(n)} gives the BIC, if \code{n}
    is the number of observations.}
}

\details{
  As for linear models the log likelihood is computed on a normal model and
  the degrees of freedom are the number of regression coefficients multiplied
  by the number of segements plus the number of estimated breakpoints plus
  1 for the error variance.

  If \code{AIC} is applied to an object of class \code{"breakpointsfull"}
  \code{breaks} can be a vector of integers and the AIC for each corresponding
  partition will be returned. By default the maximal number of breaks stored
  in the \code{object} is used. See below for an example.
}

\value{
  An object of class \code{"logLik"} or a simple vector containing
  the AIC respectively.
}

\seealso{\code{\link{breakpoints}}}

\examples{
require(ts)

## Nile data with one breakpoint: the annual flows drop in 1898
## because the first Ashwan dam was built
data(Nile)
plot(Nile)

bp.nile <- breakpoints(Nile ~ 1)
summary(bp.nile)
plot(bp.nile)

## BIC of partitions with0 to 5 breakpoints
plot(0:5, AIC(bp.nile, k = log(bp.nile$nobs)), type = "b")
## AIC
plot(0:5, AIC(bp.nile), type = "b")

## BIC, AIC, log likelihood of a single partition
bp.nile1 <- breakpoints(bp.nile, breaks = 1)
AIC(bp.nile1, k = log(bp.nile1$nobs))
AIC(bp.nile1)
logLik(bp.nile1)
}

\keyword{regression}

\eof
\name{mefp}
\alias{mefp}
\alias{mefp.formula}
\alias{mefp.efp}
\alias{print.mefp}
\alias{monitor}
\title{Monitoring of Empirical Fluctuation Processes}
\description{
  Online monitoring of structural breaks in a linear regression model. A
  sequential fluctuation test based on parameter estimates or OLS residualas
  signals structural breaks.
}
\usage{
mefp(obj, ...)

\method{mefp}{formula}(formula, type = c("OLS-CUSUM", "OLS-MOSUM", "RE", "ME",
    "fluctuation"), data, h = 1, alpha = 0.05,
    functional = c("max", "range"), period = 10,
    tolerance = .Machine$double.eps^0.5, CritvalTable = NULL,
    rescale = NULL, border = NULL, ...)

\method{mefp}{efp}(obj, alpha=0.05, functional = c("max", "range"),
    period = 10, tolerance = .Machine$double.eps^0.5,
    CritvalTable = NULL, rescale = NULL, border = NULL, ...)

monitor(obj, data = NULL, verbose = TRUE)
}
\arguments{
  \item{formula}{a symbolic description for the model to be tested.}
  \item{data}{an optional data frame containing the variables in the model. By
     default the variables are taken from the environment which \code{efp} is
     called from.}
  \item{type}{specifies which type of fluctuation process will be
     computed.}
  \item{h}{(only used for MOSUM/ME processes). A numeric scalar from interval
    (0,1) specifying the size of the data window relative to the sample
    size.}
  \item{obj}{Object of class \code{"efp"} (for \code{mefp}) or
    \code{"mefp"} (for \code{monitor}).}
  \item{alpha}{Significance level of the test, i.e., probability of
    type I error.}
  \item{functional}{Determines if maximum or range of parameter
    differences is used as statistic.}
  \item{period}{(only used for MOSUM/ME processes). Maximum time (relative
    to the history period) that will be monitored. Default is 10 times
    the history period.}
  \item{tolerance}{Tolerance for numeric \code{==} comparisons.}
  \item{CritvalTable}{Table of critical values, this table is
    interpolated to get critical values
    for arbitrary \code{alpha}s. The default depends on the \code{type}
    of fluctuation process (pre-computed tables are available for all
    types). \emph{This argument is under development.}}
  \item{rescale}{If \code{TRUE} the estimates will be standardized by
    the regressor matrix of the corresponding subsample similar to
    Kuan & Chen (1994); if \code{FALSE} the historic regressor matrix will
    be used. The default is to rescale the monitoring processes of type
    \code{"ME"} but not of \code{"RE"}.}
  \item{border}{An optional user-specified border function for the
      empirical process. \emph{This argument is under development.}}
  \item{verbose}{If \code{TRUE}, signal breaks by text output.}
  \item{...}{Currently not used.}
}
\details{
  \code{\link{mefp}} creates an object of class \code{"mefp"} either
  from a model formula or from an object of class \code{"efp"}. In
  addition to the arguments of \code{\link{efp}}, the type of statistic
  and a significance level for the monitoring must be specified. The
  monitoring itself is performed by \code{monitor}, which can be
  called arbitrarily often on objects of class \code{"mefp"}. If new
  data have arrived, then the empirical fluctuation process is computed
  for the new data. If the process crosses the boundaries corresponding
  to the significance level \code{alpha}, a structural break is detected
  (and signaled).

  The typical usage is to initialize the monitoring by creation of an
  object of class \code{"mefp"} either using a formula or an
  \code{"efp"} object. Data available at this stage are considered the
  \emph{history sample}, which is kept fixed during the complete
  monitoring process, and may not contain any structural changes.

  Subsequent calls to \code{monitor} perform a sequential test of the
  null hypothesis of no structural change in new data against the
  general alternative of changes in one or more of the coefficients of
  the regression model.

  The recursive
  estimates test is also called fluctuation test, therefore setting \code{type}
  to \code{"fluctuation"} was used to specify it in earlier versions of
  strucchange. It still can be used now, but will be forced to \code{"RE"}
}
\seealso{\code{\link{plot.mefp}}, \code{\link{boundary.mefp}}}
\references{
  Leisch F., Hornik K., Kuan C.-M. (2000), Monitoring
  Structural Changes with the Generalized Fluctuation Test,
  \emph{Econometric Theory}, \bold{16}, 835-854.

  Zeileis A., Leisch F., Kleiber C., Hornik K. (2003), Monitoring Structural Change
  in Dynamic Econometric Models, Report 64, SFB "Adaptive Information Systems and
  Modelling in Economics and Management Science", Vienna University of Economics,
  \url{http://www.wu-wien.ac.at/am/reports.htm#64}.
}
\examples{
df1 <- data.frame(y=rnorm(300))
df1[150:300,"y"] <- df1[150:300,"y"]+1

## use the first 50 observations as history period
e1 <- efp(y~1, data=df1[1:50,,drop=FALSE], type="ME", h=1)
me1 <- mefp(e1, alpha=0.05)

## the same in one function call
me1 <- mefp(y~1, data=df1[1:50,,drop=FALSE], type="ME", h=1,
              alpha=0.05)

## monitor the 50 next observations
me2 <- monitor(me1, data=df1[1:100,,drop=FALSE])
plot(me2)

# and now monitor on all data
me3 <- monitor(me2, data=df1)
plot(me3)


## Load dataset "USIncExp" with income and expenditure in the US
## and choose a suitable subset for the history period
data(USIncExp)
USIncExp3 <- window(USIncExp, start=c(1969,1), end=c(1971,12))
## initialize the monitoring with the formula interface
me.mefp <- mefp(expenditure~income, type="ME", rescale=TRUE,
                   data=USIncExp3, alpha=0.05)

## monitor the new observations for the year 1972
USIncExp3 <- window(USIncExp, start=c(1969,1), end=c(1972,12))
me.mefp <- monitor(me.mefp)

## monitor the new data for the years 1973-1976
USIncExp3 <- window(USIncExp, start=c(1969,1), end=c(1976,12))
me.mefp <- monitor(me.mefp)
plot(me.mefp, functional = NULL)
}
\keyword{regression}

\eof
\name{plot.Fstats}
\alias{plot.Fstats}
\alias{lines.Fstats}
\title{Plot F Statistics}
\description{Plotting method for objects of class \code{"Fstats"}}
\usage{
\method{plot}{Fstats}(x, pval = FALSE, asymptotic = FALSE, alpha = 0.05,
    boundary = TRUE, aveF = FALSE, xlab = "Time", ylab = NULL,
    ylim = NULL, ...)
}

\arguments{
 \item{x}{an object of class \code{"Fstats"}.}
 \item{pval}{logical. If set to \code{TRUE} the corresponding p values instead
     of the original F statistics will be plotted.}
 \item{asymptotic}{logical. If set to \code{TRUE} the asymptotic (chi-square)
     distribution instead of the exact (F) distribution will be used to compute
     the p values (only if \code{pval} is \code{TRUE}).}
 \item{alpha}{numeric from interval (0,1) indicating the confidence level for
     which the boundary of the supF test will be computed.}
 \item{boundary}{logical. If set to \code{FALSE} the boundary will be computed
     but not plotted.}
 \item{aveF}{logical. If set to \code{TRUE} the boundary of the aveF test will
     be plotted. As this is a boundary for the mean of the F statistics rather
     than for the F statistics themselves a dashed line for the mean of the F
     statistics will also be plotted.}
 \item{xlab, ylab, ylim, ...}{high-level \code{\link{plot}} function parameters.}}

\references{
Andrews D.W.K. (1993), Tests for parameter instability and structural
change with unknown change point, \emph{Econometrica}, \bold{61}, 821-856.

Hansen B. (1992), Tests for parameter instability in regressions with I(1)
processes, \emph{Journal of Business & Economic Statistics}, \bold{10}, 321-335.

Hansen B. (1997), Approximate asymptotic p values for structural-change
tests, \emph{Journal of Business & Economic Statistics}, \bold{15}, 60-67. }

\seealso{\code{\link{Fstats}}, \code{\link{boundary.Fstats}},
\code{\link{sctest.Fstats}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains
## constant over the years for potential break points between 1941
## (corresponds to from = 0.5) and 1962 (corresponds to to = 0.85)
## compute F statistics
fs <- Fstats(nhtemp ~ 1, from = 0.5, to = 0.85)
## plot the F statistics
plot(fs, alpha = 0.01)
## and the corresponding p values
plot(fs, pval = TRUE, alpha = 0.01)
## perform the aveF test
sctest(fs, type = "aveF")
}
\keyword{hplot}

\eof
\name{plot.efp}
\alias{plot.efp}
\alias{lines.efp}
\title{Plot Empirical Fluctuation Process}
\description{Plot and lines method for objects of class \code{"efp"}}
\usage{
\method{plot}{efp}(x, alpha = 0.05, alt.boundary = FALSE, boundary = TRUE,
    functional = "max", main = NULL,  ylim = NULL,
    ylab = "empirical fluctuation process", ...)
\method{lines}{efp}(x, functional = "max", ...)
}

\arguments{
  \item{x}{an object of class \code{"efp"}.}
  \item{alpha}{numeric from interval (0,1) indicating the confidence level for
     which the boundary of the corresponding test will be computed.}
  \item{alt.boundary}{logical. If set to \code{TRUE} alternative boundaries
     (instead of the standard linear boundaries) will be plotted (for CUSUM
     processes only).}
  \item{boundary}{logical. If set to \code{FALSE} the boundary will be computed
     but not plotted.}
  \item{functional}{indicates which functional should be applied to the
     process before plotting and which boundaries should be used. If set to \code{NULL}
     a multiple process with boundaries for the \code{"max"} functional is plotted.
     For more details see below.}
  \item{main, ylim, ylab, ...}{high-level \code{\link{plot}} function
 parameters.}
}

\details{Plots are available for the \code{"max"} functional for all process types.
For Brownian bridge type processes the maximum or mean squared Euclidian norm
(\code{"maxL2"} and \code{"meanL2"}) can be used for aggregating before plotting.
No plots are available for the \code{"range"} functional.

Alternative boundaries that are proportional to the standard deviation
of the corresponding limiting process are available for processes with Brownian
motion or Brownian bridge limiting processes.
}

\value{\code{\link{efp}} returns an object of class \code{"efp"} which inherits
from the class \code{"ts"} or \code{"mts"} respectively. The function
\code{\link{plot}} has a method to plot the
empirical fluctuation process; with \code{sctest} the corresponding test for
structural change can be performed.}

\references{Brown R.L., Durbin J., Evans J.M. (1975), Techniques for
testing constancy of regression relationships over time, \emph{Journal of the
Royal Statistal Society}, B, \bold{37}, 149-163.

Chu C.-S., Hornik K., Kuan C.-M. (1995), MOSUM tests for parameter
constancy, \emph{Biometrika}, \bold{82}, 603-617.

Chu C.-S., Hornik K., Kuan C.-M. (1995), The moving-estimates test for
parameter stability, \emph{Econometric Theory}, \bold{11}, 669-720.

Krämer W., Ploberger W., Alt R. (1988), Testing for structural change in
dynamic models, \emph{Econometrica}, \bold{56}, 1355-1369.

Kuan C.-M., Hornik K. (1995), The generalized fluctuation test: A
unifying view, \emph{Econometric Reviews}, \bold{14}, 135 - 161.

Kuan C.-M., Chen (1994), Implementing the fluctuation and moving estimates
tests in dynamic econometric models, \emph{Economics Letters}, \bold{44},
235-239.

Ploberger W., Krämer W. (1992), The CUSUM test with OLS residuals,
\emph{Econometrica}, \bold{60}, 271-285.

Zeileis A. (2000), p Values and Alternative Boundaries for CUSUM Tests,
Working Paper 78, SFB "Adaptive Information Systems and Modelling in Economics
and Management Science", Vienna University of Economics,
\url{http://www.wu-wien.ac.at/am/wp00.htm#78}.
}

\seealso{\code{\link{efp}}, \code{\link{boundary.efp}},
\code{\link{sctest.efp}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains
## constant over the years
## compute Rec-CUSUM fluctuation process
temp.cus <- efp(nhtemp ~ 1)
## plot the process
plot(temp.cus, alpha = 0.01)
## and calculate the test statistic
sctest(temp.cus)

## compute (recursive estimates) fluctuation process
## with an additional linear trend regressor
lin.trend <- 1:60
temp.me <- efp(nhtemp ~ lin.trend, type = "fluctuation")
## plot the bivariate process
plot(temp.me, functional = NULL)
## and perform the corresponding test
sctest(temp.me)
}
\keyword{hplot}

\eof
\name{plot.mefp}
\alias{plot.mefp}
\alias{lines.mefp}
\title{Plot Methods for mefp Objects}
\description{This is a method of the generic \code{\link{plot}} function for
  for \code{"mefp"} objects as returned by \code{\link{mefp}} or
  \code{\link{monitor}}. It plots the emprical fluctuation process (or a
  functional therof) as a time series plot, and includes boundaries
  corresponding to the significance level of the monitoring procedure.
}
\usage{
\method{plot}{mefp}(x, boundary = TRUE, functional = "max", main = NULL,
    ylab = "empirical fluctuation process", ylim = NULL, ...)
}
\arguments{
  \item{x}{an object of class \code{"mefp"}.}
  \item{boundary}{if \code{FALSE}, plotting of boundaries is suppressed.}
  \item{functional}{indicates which functional should be applied to a
    multivariate empirical process. If set to \code{NULL} all dimensions
    of the process (one process per coefficient in the linear model) are
    plotted. }
  \item{main, ylab, ylim, ...}{high-level \code{\link{plot}} function parameters.}
}
\seealso{\code{\link{mefp}}}
\examples{
df1 <- data.frame(y=rnorm(300))
df1[150:300,"y"] <- df1[150:300,"y"]+1
me1 <- mefp(y~1, data=df1[1:50,,drop=FALSE], type="ME", h=1,
              alpha=0.05)
me2 <- monitor(me1, data=df1)

plot(me2)
}
\keyword{hplot}

\eof
\name{recresid}
\alias{recresid}
\alias{recresid.default}
\alias{recresid.formula}
\alias{recresid.lm}
\title{Recursive Residuals}
\description{
  A generic function for computing the recursive residuals
  (standardized one step prediction errors) of a linear regression model.
}
\usage{
\method{recresid}{default}(x, y, \dots)
\method{recresid}{formula}(formula, data = list(), \dots)
\method{recresid}{lm}(x, data = list(), \dots)
}
\arguments{
  \item{x, y, formula}{specification of the linear regression model:
    either by a regressor matrix \code{x} and a response variable \code{y},
    or by a \code{formula} or by a fitted object \code{x} of class \code{"lm"}.}
  \item{data}{an optional data frame containing the variables in the model. By
     default the variables are taken from the environment which \code{recresid} is
     called from. Specifying \code{data} might also be necessary when applying
     \code{recresid} to a fitted model of class \code{"lm"} if this does not
     contain the regressor matrix and the response.}
  \item{\dots}{\emph{currently not used.}}
}
\details{
  Under the usual assumptions for the linear regression model the
  recdursive residuals are (asymptotically) normal and
  i.i.d. (see Brown, Durbin, Evans (1975) for details).
}
\value{
  A vector containing the recursive residuals.
}
\references{
Brown R.L., Durbin J., Evans J.M. (1975), Techniques for
testing constancy of regression relationships over time, \emph{Journal of the
Royal Statistal Society}, B, \bold{37}, 149-163.
}


\seealso{\code{\link{efp}}}

\examples{
x <- rnorm(100)
x[51:100] <- x[51:100] + 2
rr <- recresid(x ~ 1)
plot(cumsum(rr), type = "l")

plot(efp(x ~ 1, type = "Rec-CUSUM"))
}
\keyword{regression}

\eof
\name{root.matrix}
\alias{root.matrix}
\alias{solve.crossprod}
\alias{solve.pd}
\title{Root of a Matrix}
\description{Computes the root of a symmetric and positive semidefinite matrix.}
\usage{
root.matrix(X)
solve.pd(X)
}
\arguments{
  \item{X}{a symmetric and positive semidefinite matrix}
}

\value{a symmetric matrix of same dimensions as \code{X}}

\examples{
X <- matrix(c(1,2,2,8), ncol=2)
test <- root.matrix(X)
## control results
X
test \%*\% test
}
\keyword{algebra}

\eof
\name{sctest.Fstats}
\alias{sctest.Fstats}
\title{supF-, aveF- and expF-Test}
\description{Performs the supF-, aveF- or expF-test}
\usage{
\method{sctest}{Fstats}(x, type = c("supF", "aveF", "expF"),
    asymptotic = FALSE, ...)
}
\arguments{
 \item{x}{an object of class \code{"Fstats"}.}
 \item{type}{a character string specifying which test will be performed.}
 \item{asymptotic}{logical. Only necessary if \code{x} contains just a single F
     statistic and type is \code{"supF"} or \code{"aveF"}. If then set to
     \code{TRUE} the asymptotic (chi-square) distribution instead of the exact
     (F) distribution will be used to compute the p value.}
 \item{...}{currently not used.}
}

\details{If \code{x} contains just a single F statistic and type is
\code{"supF"} or \code{"aveF"} the Chow test will be performed.

The original GAUSS code for computing the p values of the supF-, aveF- and
expF-test was written by Bruce Hansen and is available from
\url{http://www.ssc.wisc.edu/~bhansen/}. R port by Achim Zeileis.
}

\value{an object of class \code{"htest"} containing:
\itemize
\item{statistic}{the test statistic}
\item{p.value}{the corresponding p value}
\item{method}{a character string with the method used}
\item{data.name}{a character string with the data name}}

\references{
Andrews D.W.K. (1993), Tests for parameter instability and structural
change with unknown change point, \emph{Econometrica}, \bold{61}, 821-856.

Andrews D.W.K., Ploberger W. (1994), Optimal tests when a nuisance parameter
is present only under the alternative, \emph{Econometrica}, \bold{62}, 1383-1414.

Hansen B. (1992), Tests for parameter instability in regressions with I(1)
processes, \emph{Journal of Business & Economic Statistics}, \bold{10}, 321-335.

Hansen B. (1997), Approximate asymptotic p values for structural-change
tests, \emph{Journal of Business & Economic Statistics}, \bold{15}, 60-67. }

\seealso{\code{\link{Fstats}}, \code{\link{plot.Fstats}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains
## constant over the years for potential break points between 1941
## (corresponds to from = 0.5) and 1962 (corresponds to to = 0.85)
## compute F statistics
fs <- Fstats(nhtemp ~ 1, from = 0.5, to = 0.85)
## plot the F statistics
plot(fs, alpha = 0.01)
## and the corresponding p values
plot(fs, pval = TRUE, alpha = 0.01)
## perform the aveF test
sctest(fs, type = "aveF")
}
\keyword{htest}

\eof
\name{sctest.formula}
\alias{sctest}
\alias{sctest.formula}
\title{Structural Change Tests}
\description{Performs tests for structural change.}
\usage{
\method{sctest}{formula}(formula, type = <<see below>>, h = 0.15,
    alt.boundary = FALSE, functional = c("max", "range",
    "maxL2", "meanL2"), from = 0.15, to = NULL, point = 0.5,
    asymptotic = FALSE, data, ...)
}

\arguments{
  \item{formula}{a formula describing the model to be tested.}
  \item{type}{a character string specifying the structural change test that ist
     to be performed. Besides the tests types described in \code{\link{efp}}
     and \code{\link{sctest.Fstats}}. The Chow test and the Nyblom-Hansen test
     can be performed by setting type to \code{"Chow"} or \code{"Nyblom-Hansen"},
     respectively.}
  \item{h}{numeric from interval (0,1) specifying the bandwidth. Determins the
     size of the data window relative to sample size (for MOSUM and ME tests
     only).}
  \item{alt.boundary}{logical. If set to \code{TRUE} alternative boundaries
     (instead of the standard linear boundaries) will be used (for CUSUM
     processes only).}
  \item{functional}{indicates which functional should be used to aggregate
     the empirical fluctuation processes to a test statistic.}
 \item{from, to}{numerics. If \code{from} is smaller than 1 they are
     interpreted as percentages of data and by default \code{to} is taken to be
     the 1 - \code{from}. F statistics will be calculated for the observations
     \code{(n*from):(n*to)}, when \code{n} is the number of observations in the
     model. If \code{from} is greater than 1 it is interpreted to be the index
     and \code{to} defaults to \code{n - from}. (for F tests only)}
  \item{point}{parameter of the Chow test for the potential change point.
    Interpreted analogous to the \code{from} parameter. By
     default taken to be \code{floor(n*0.5)} if \code{n} is the  number of
     observations in the model.}
  \item{asymptotic}{logical. If \code{TRUE} the asymptotic (chi-square)
     distribution instead of the exact (F) distribution will be used to compute
     the p value (for Chow test only).}
  \item{data}{an optional data frame containing the variables in the model. By
     default the variables are taken from the environment which
     \code{sctest} is called from.}
  \item{...}{further arguments passed to \code{\link{efp}} or
     \code{\link{Fstats}}.} }

\details{\code{sctest.formula} is mainly a wrapper for \code{\link{sctest.efp}}
and \code{\link{sctest.Fstats}} as it fits an empirical fluctuation process
first or computes the F statistics respectively and subsequently performs the
corresponding test. The Chow test and the Nyblom-Hansen test are available explicitely here.}

\value{an object of class \code{"htest"} containing:
\itemize
\item{statistic}{the test statistic}
\item{p.value}{the corresponding p value}
\item{method}{a character string with the method used}
\item{data.name}{a character string with the data name}}

\seealso{\code{\link{sctest.efp}}, \code{\link{sctest.Fstats}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains
## constant over the years with the Standard CUSUM test
sctest(nhtemp ~ 1)
## with the Chow test (under the alternative that there is a change 1941)
sctest(nhtemp ~ 1, type = "Chow", point = c(1941,1))
}
\keyword{htest}

\eof
\name{sctest.efp}
\alias{sctest.efp}
\title{Generalized Fluctuation Tests}
\description{Performs a generalized fluctuation test.}
\usage{
\method{sctest}{efp}(x, alt.boundary = FALSE,
    functional = c("max", "range", "maxL2", "meanL2"), ...)
}

\arguments{
  \item{x}{an object of class \code{"efp"}.}
  \item{alt.boundary}{logical. If set to \code{TRUE} alternative boundaries
 (instead of
     the standard linear boundaries) will be used (for CUSUM
     processes only).}
  \item{functional}{indicates which functional should be applied to the
     empirical fluctuation process.}
  \item{...}{currently not used.}
}

\details{The critical values for the MOSUM tests and the ME test are just
tabulated for confidence levels between 0.1 and 0.01, thus the p
value approximations will be poor for other p values. Similarly the
critical values for the maximum and mean squared Euclidian norm (\code{"maxL2"}
and \code{"meanL2"}) are tabulated for confidence levels between 0.2 and 0.005.}

\value{an object of class \code{"htest"} containing:
\itemize
\item{statistic}{the test statistic}
\item{p.value}{the corresponding p value}
\item{method}{a character string with the method used}
\item{data.name}{a character string with the data name}}

\references{Brown R.L., Durbin J., Evans J.M. (1975), Techniques for
testing constancy of regression relationships over time, \emph{Journal of the
Royal Statistal Society}, B, \bold{37}, 149-163.

Chu C.-S., Hornik K., Kuan C.-M. (1995), MOSUM tests for parameter
constancy, \emph{Biometrika}, \bold{82}, 603-617.

Chu C.-S., Hornik K., Kuan C.-M. (1995), The moving-estimates test for
parameter stability, \emph{Econometric Theory}, \bold{11}, 669-720.

Krämer W., Ploberger W., Alt R. (1988), Testing for structural change in
dynamic models, \emph{Econometrica}, \bold{56}, 1355-1369.

Kuan C.-M., Hornik K. (1995), The generalized fluctuation test: A
unifying view, \emph{Econometric Reviews}, \bold{14}, 135 - 161.

Kuan C.-M., Chen (1994), Implementing the fluctuation and moving estimates
tests in dynamic econometric models, \emph{Economics Letters}, \bold{44},
235-239.

Ploberger W., Krämer W. (1992), The CUSUM Test with OLS Residuals,
\emph{Econometrica}, \bold{60}, 271-285.

Zeileis A. (2000), p Values and Alternative Boundaries for CUSUM Tests,
Working Paper 78, SFB "Adaptive Information Systems and Modelling in Economics
and Management Science", Vienna University of Economics,
\url{http://www.wu-wien.ac.at/am/wp00.htm#78}.
}

\seealso{\code{\link{efp}}, \code{\link{plot.efp}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains
## constant over the years compute OLS-CUSUM fluctuation process
temp.cus <- efp(nhtemp ~ 1, type = "OLS-CUSUM")
## plot the process with alternative boundaries
plot(temp.cus, alpha = 0.01, alt.boundary = TRUE)
## and calculate the test statistic
sctest(temp.cus)

## compute moving estimates fluctuation process
temp.me <- efp(nhtemp ~ 1, type = "ME", h = 0.2)
## plot the process with functional = "max"
plot(temp.me)
## and perform the corresponding test
sctest(temp.me)
}
\keyword{htest}

\eof
\name{strucchange.internal}
\alias{sc.beta.sup}
\alias{sc.beta.ave}
\alias{sc.beta.exp}
\alias{sc.me}
\alias{sc.meanL2}
\alias{sc.maxL2}
\alias{pvalue.efp}
\alias{pvalue.Fstats}
\alias{monitorMECritval}
\alias{monitorMECritvalData}
\alias{monitorMECritvalTable}
\alias{monitorRECritval}
\alias{monitorRECritvalData}
\alias{monitorRECritvalTable}
\alias{pargmaxV}
\title{Internal strucchange objects}
\description{
  These are not to be called by the user.
}
\keyword{internal}

\eof
