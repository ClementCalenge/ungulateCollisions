% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/code.R
\name{prepareFit}
\alias{prepareFit}
\title{Bayesian Model To Identify Factors Affecting Wildlife-Vehicle Collisions}
\usage{
prepareFit(X, alphas, collisions, nYear, Area, departement)
}
\arguments{
\item{X}{a data.frame containing the numeric variables supposed to
have an effect on the wildlife-vehicle collisions (columns) for
spatial unit (rows).}

\item{alphas}{a character string vector with length equal to
\code{ncol(X)} giving the name, for each column of \code{X}, of the
group of variables to which each column of \code{X} belongs. This
vector is used to describe how to group the variables so that they
are either in the model or out of the model together.  In general,
\code{alphas[i]} is set equal to \code{names(X)[i]}, except if
several columns are to be considered as one variable. For example,
if \code{X} contains several dummy variables (i.e. coded 0/1)
corresponding to the same factor "habitat" (e.g. Urban,
agriculture, forest), then, we might set \code{names(X)[i] <-
"habitat"} for all i corresponding to this habitat.}

\item{collisions}{an integer vector with length equal to
\code{nrow(X)} defining the number of wildlife-ungulate collisions
in each studied spatial unit.}

\item{nYear}{an integer vector with length equal to \code{nrow(X)}
defining the number of Years during which wildlife-ungulate
collisions have been recorded in each studied spatial unit.}

\item{Area}{a numeric vector with length equal to \code{nrow(X)}
defining the area of each studied spatial unit.}

\item{departement}{a character vector with length equal to
\code{nrow(X)} defining the department of each studied spatial unit.}
}
\value{
a list with all elements required for the fit of the model
with JAGS, that is: (i) \code{data4jags}: the list of the data
required by the model, to be passed to the argument \code{data} of
the function \code{jags.model} of the package \code{rjags}, (ii)
\code{ini}: list of starting values for the parameters, to be
passed to the argument \code{init} of the function
\code{jags.model}, (iii) \code{modelstring}: a character string
containing the model fit by JAGS, (iv) \code{coefnames}: vector of
character strings containing the names of the coefficients of
interest in the model, to be passed to the argument
\code{variable.names} of the function \code{coda.samples} of the
package \code{rjags}.
}
\description{
Bayesian Model To Identify Factors Affecting Wildlife-Vehicle Collisions
}
\details{
\code{prepareFit} prepares the elements required for the
fit of the Kuo-Mallik model (a Bayesian model used to identify
factors affecting wildlife-vehicle collisions).
}
\examples{
## Load the data:
data("dataCollision")

## Prepare the data
X <- dataCollision$RedDeer$X

## Sets the elevation <600 to NULL (redundant with other elevations:
## if the elevation is neither >1500, nor between 600 and 1500,
## it is necessarily < 600
X$elev600 <- NULL

## Set the variable roads to NULL (not useful here, as noted in the
## paper)
X$roads <- NULL

## prepares the alpha:
## vector alphas: the name of alphas is the same as the columns of X...
alphas <- names(X)
## ... except for elevation and habitat,
## which correspond to several variables in X
alphas[7:8] <- "elev"
alphas[11:14] <- "hab"

## Note that the variables Agriculture, Open, Urban, and Forest
## are strongly correlated together. To reduce the correlation and
## improve mixing, we transform these variables with the help of a PCA
## (package ade4)
hab <- X[,11:14]
X <- X[,-c(11:14)]
pc <- ade4::dudi.pca(hab,scannf = FALSE, nf=4)

## We scale the variables to improve mixing
X <- as.data.frame(scale(X))

## and we add the transformed habitat variables from the PCA:
li1 <- pc$l1
X$AX1 <- li1[,1]
X$AX2 <- li1[,2]
X$AX3 <- li1[,3]
X$AX4 <- li1[,4]

## Use of the function preparefit
pf <- prepareFit(X, alphas, dataCollision$RedDeer$coll,
                 dataCollision$RedDeer$Y, dataCollision$RedDeer$Area,
                 dataCollision$RedDeer$departement)

\dontrun{
## WARNING: very long execution (about 1 hour)!!
## the results are stored in the dataset "modelRedDeer"
## But if you want to try it, the data are now ready for the fit:

mo <- jags.model(textConnection(pf$modelstring), ini=pf$ini, data=pf$data4jags)
update(mo, n.iter=1000)
modelRedDeer <- coda.samples(mo, variable.names = pf$coefnames,
                             n.iter = 500000, thin=100)
}
}
\author{
Clement Calenge, \email{clement.calenge@oncfs.gouv.fr}
}
\seealso{
\code{\link[rjags]{jags.model}}, \code{\link[rjags]{coda.samples}}
}

