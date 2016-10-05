#' ungulateCollisions: Identifying the Factors Affecting Wildlife Ungulate Collisions
#'
#' The package ungulateCollisions provides the dataset as
#' well as the R code used by Saint Andrieux et al. (in prep.) to
#' identify the environmental factors affecting the wildlife Ungulate
#' collisions. This is done using a Bayesian model of the type
#' described by Kuo and Mallik (1998).
#'
#' @references Kuo, L. and Mallik, B. 1998.Variable selection for
#' regression models. Sankhya: The Indian Journal of Statistics,
#' Series B, 60, 65-81.
#'
#' @references Saint Andrieux, C., Calenge, C. and Bonenfant, C. in
#' prep. Comparison of ecological, biological and anthropogenic causes
#' of vehicle-wildlife collisions among three large mammalian species.
#'
#' @docType package
#' @name ungulateCollisions
NULL


#' @title Ungulate-Vehicle Collisions in Nine French Departements
#'
#' @details Dataset describing the number of Ungulate-Vehicle
#' collisions in the management units of Nine French Departements.
#'
#'
#' @format A list with three elements named \code{RedDeer},
#' \code{RoeDeer} and \code{WildBoar}, each element being itself a
#' list containing the following elements:
#'
#' \describe{
#'
#' \item{X}{a data.frame containing the variables describing the
#' management units, and supposed to have an effect on the collisions:
#' (i) \code{forFrag}: forest fragmentation, (ii) \code{urbFrag}:
#' urban fragmentation, (iii) \code{locr}: density of local roads,
#' (iv) \code{regr}: density of regional roads, (v) \code{natr}:
#' density of national roads, (vi) \code{motr}: density of motorways,
#' (vii) \code{elev600} percentage of the management unit covered by
#' areas with elevation lower than 600 m asl, (viii)
#' \code{elev600_1500} percentage of the management unit covered by
#' areas with elevation comprised between 600 m and 1500 m asl, (ix)
#' \code{elev600} percentage of the management unit covered by areas
#' with elevation greater than 1500 m asl, (x) \code{sinus}: average
#' sinuosity of the roads in the management unit, (xi) \code{hunt}:
#' the hunting bag, (xii) \code{Agriculture,Open,Urban,Forest}:
#' proportion of the management unit covered by these habitat types.}
#'
#' \item{coll}{the number of collisions recorded in each management
#' unit}
#'
#' \item{Y}{the number of years during which collisions were recorded
#' in each management unit}
#'
#' \item{Area}{the Area of each management unit}
#'
#' \item{departement}{the departement corresponding to each management unit}
#' }
#'
#' @references Saint-Andrieux, C., Calenge, C. and Bonenfant, C. in
#' prep. Comparison of ecological, biological and anthropogenic causes
#' of vehicle-wildlife collisions among three large mammalian species
#'
"dataCollision"



#' @title Kuo and Mallik Model of the Factors Affecting Ungulate-Vehicle Collisions
#'
#' @details Dataset containing the results of the MCMC fit of the
#' model used to identify the variables affecting Red Deer-Vehicle
#' collisions (as described by Kuo and Mallik, 1998). The fit was
#' performed using the package rjags. See the examples of the function
#' \code{prepareFit} to see how this dataset was generated.
#'
#'
#' @format An object of class \code{mcmc.list}
#'
#' @seealso \code{\link{prepareFit}}
#'
#' @references Kuo, L. and Mallik, B. 1998.Variable selection for
#' regression models. Sankhya: The Indian Journal of Statistics,
#' Series B, 60, 65-81.
#'
"modelRedDeer"



#' @title Final Model of the Factors Affecting Ungulate-Vehicle Collisions
#'
#' @details Dataset containing the results of the MCMC fit of the
#' final model predicting the Red Deer-Vehicle collisions as a
#' function of the forest fragmentation and hunting bag. The fit was
#' performed using the package rjags. See the examples of the function
#' \code{finalModel} to see how this dataset was generated.
#'
#'
#' @format An object of class \code{mcmc.list}
#'
#' @seealso \code{\link{finalModel}}
#'
#' @references Kuo, L. and Mallik, B. 1998.Variable selection for
#' regression models. Sankhya: The Indian Journal of Statistics,
#' Series B, 60, 65-81.
#'
"finalModelRedDeer"






#' @title Bayesian Model To Identify Factors Affecting Wildlife-Vehicle Collisions
#' @export
#' @aliases prepareFit
#' @details \code{prepareFit} prepares the elements required for the
#' fit of the Kuo-Mallik model (a Bayesian model used to identify
#' factors affecting wildlife-vehicle collisions).
#'
#' @param X a data.frame containing the numeric variables supposed to
#' have an effect on the wildlife-vehicle collisions (columns) for
#' spatial unit (rows).
#' @param alphas a character string vector with length equal to
#' \code{ncol(X)} giving the name, for each column of \code{X}, of the
#' group of variables to which each column of \code{X} belongs. This
#' vector is used to describe how to group the variables so that they
#' are either in the model or out of the model together.  In general,
#' \code{alphas[i]} is set equal to \code{names(X)[i]}, except if
#' several columns are to be considered as one variable. For example,
#' if \code{X} contains several dummy variables (i.e. coded 0/1)
#' corresponding to the same factor "habitat" (e.g. Urban,
#' agriculture, forest), then, we might set \code{names(X)[i] <-
#' "habitat"} for all i corresponding to this habitat.
#' @param collisions an integer vector with length equal to
#' \code{nrow(X)} defining the number of wildlife-ungulate collisions
#' in each studied spatial unit.
#' @param nYear an integer vector with length equal to \code{nrow(X)}
#' defining the number of Years during which wildlife-ungulate
#' collisions have been recorded in each studied spatial unit.
#' @param Area a numeric vector with length equal to \code{nrow(X)}
#' defining the area of each studied spatial unit.
#' @param departement a character vector with length equal to
#' \code{nrow(X)} defining the department of each studied spatial unit.
#'
#'
#' @return a list with all elements required for the fit of the model
#' with JAGS, that is: (i) \code{data4jags}: the list of the data
#' required by the model, to be passed to the argument \code{data} of
#' the function \code{jags.model} of the package \code{rjags}, (ii)
#' \code{ini}: list of starting values for the parameters, to be
#' passed to the argument \code{init} of the function
#' \code{jags.model}, (iii) \code{modelstring}: a character string
#' containing the model fit by JAGS, (iv) \code{coefnames}: vector of
#' character strings containing the names of the coefficients of
#' interest in the model, to be passed to the argument
#' \code{variable.names} of the function \code{coda.samples} of the
#' package \code{rjags}.
#'
#' @seealso \code{\link[rjags]{jags.model}}, \code{\link[rjags]{coda.samples}}
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
#' @import stats
#' @import coda
#' @import rjags
#' @examples
#'
#' ## Load the data:
#' data("dataCollision")
#'
#' ## Prepare the data
#' X <- dataCollision$RedDeer$X
#'
#' ## Sets the elevation <600 to NULL (redundant with other elevations:
#' ## if the elevation is neither >1500, nor between 600 and 1500,
#' ## it is necessarily < 600
#' X$elev600 <- NULL
#'
#' ## prepares the alpha:
#' ## vector alphas: the name of alphas is the same as the columns of X...
#' alphas <- names(X)
#' ## ... except for elevation and habitat,
#' ## which correspond to several variables in X
#' alphas[7:8] <- "elev"
#' alphas[11:14] <- "hab"
#'
#' ## Note that the variables Agriculture, Open, Urban, and Forest
#' ## are strongly correlated together. To reduce the correlation and
#' ## improve mixing, we transform these variables with the help of a PCA
#' ## (package ade4)
#' hab <- X[,11:14]
#' X <- X[,-c(11:14)]
#' pc <- ade4::dudi.pca(hab,scannf = FALSE, nf=4)
#'
#' ## We scale the variables to improve mixing
#' X <- as.data.frame(scale(X))
#'
#' ## and we add the transformed habitat variables from the PCA:
#' li1 <- pc$l1
#' X$AX1 <- li1[,1]
#' X$AX2 <- li1[,2]
#' X$AX3 <- li1[,3]
#' X$AX4 <- li1[,4]
#'
#' ## Use of the function preparefit
#' pf <- prepareFit(X, alphas, dataCollision$RedDeer$coll,
#'                  dataCollision$RedDeer$Y, dataCollision$RedDeer$Area,
#'                  dataCollision$RedDeer$departement)
#'
#' \dontrun{
#' ## WARNING: very long execution (about 1 hour)!!
#' ## the results are stored in the dataset "modelRedDeer"
#' ## But if you want to try it, the data are now ready for the fit:
#'
#' mo <- jags.model(textConnection(pf$modelstring), ini=pf$ini, data=pf$data4jags)
#' update(mo, n.iter=1000)
#' modelRedDeer <- coda.samples(mo, variable.names = pf$coefnames,
#'                              n.iter = 500000, thin=100)
#' }
#'
prepareFit <- function(X, alphas, collisions, nYear, Area, departement)
{
    if (!inherits(X, "data.frame"))
        stop("X should be a data.frame")
    if (nrow(X)!=length(collisions))
        stop("The number of rows of X should be equal to length(collisions)")
    if (nrow(X)!=length(nYear))
        stop("The number of rows of X should be equal to length(nYear)")
    if (nrow(X)!=length(Area))
        stop("The number of rows of X should be equal to length(Area)")
    if (nrow(X)!=length(departement))
        stop("The number of rows of X should be equal to length(departement)")
    if (!is.character(departement))
        stop("departement should be character")
    if (length(alphas)!=ncol(X))
        stop("The length of alpha should be equal to ncol(X)")

    ## Data preparation
    data4jags <- lapply(X, function(x) x)

    names(data4jags) <- toupper(names(data4jags))

    data4jags$collisions <- collisions
    data4jags$departement <- as.numeric(factor(departement))
    data4jags$nYear <- as.numeric(factor(nYear))
    data4jags$Area <- as.numeric(factor(Area))
    data4jags$ndep <- length(unique(departement))
    data4jags$nUG <- nrow(X)

    ##
    ##
    ## Starting values

    ## starting values for the parameters alpha. At the beginning
    ## all the variables are in the model
    alph <- paste("alpha",unique(alphas),sep="_")
    init1 <- lapply(alph, function(x) 1)
    names(init1) <- alph

    ## starting values for the parameters beta. Random draw
    ## from the normal distribution
    beta <- paste("beta", names(X), sep="_")
    init2 <- lapply(beta, function(x) rnorm(1))
    names(init2) <- beta

    ## All starting values
    ini <- c(init1, init2)

    ## formula describing the core of the model
    spl1 <- split(beta, alphas)
    spl2 <- split(toupper(names(X)), alphas)
    na <- paste("alpha",names(spl1), sep="_")
    formcoeur <- paste(paste(na,paste0("(",sapply(1:length(spl1), function(i) {
                                         paste(paste0(spl1[[i]],"*",spl2[[i]], "[i]"),collapse = "+")
                                     }),")"), sep="*"), collapse="+")

    ##
    ##
    ## Model Code

    formulalik <- paste0("lambda[i] <- nYear[i]*Area[i]*exp(",
                         "efal[departement[i]]+",
                         formcoeur,
                         "+rs[i])")

    priors <- c("model {",
                "## Priors",
                "## Coefficients alpha (inclusion in the model)",
                paste0(alph,"~dbern(0.5)"),
                "\n",
                "## Coefficients beta (importance in the model)",
                paste0(beta, "~dnorm(0.0, 0.01)"),
                "\n",
                "## precision of random effects departement",
                "sigefa~dgamma(0.01, 0.01)",
                "muefa~dnorm(0, 0.01)",
                "\n",
                "## Random effect departement",
                "for (d in 1:ndep) {",
                "efal[d]~dnorm(muefa, sigefa)",
                "}",
                "\n",
                "## precision of overdispersion residuals",
                "sigeps~dgamma(0.01, 0.01)",
                "\n",
                "## Likelihood",
                "for (i in 1:nUG) {",
                "rs[i]~dnorm(0,sigeps)",
                formulalik,
                "collisions[i]~dpois(lambda[i])",
                "}",
                "\n",
                "}\n")
    modelstring <- (paste(priors, collapse = "\n"))

    ## Coefficient names
    coefnames <- c(alph, beta, "sigeps", "sigefa", "muefa")

    ## Resultats
    resu <- list(data4jags=data4jags,
                 ini=ini,
                 modelstring=modelstring,
                 coefnames=coefnames)

    ## output
    return(resu)

}






#' @title Probability of Variables and Models Calculated with a Kuo-Mallik Approach
#' @export
#' @aliases probabilityKM
#' @details \code{probabilityKM} calculates the probability of
#' inclusion of the variables, as well as the probability that each
#' possible model is true, using a Kuo-Mallik approach. The model used
#' should have been fitted with the package rjags, using the model
#' string returned by the function \code{prepareFit}.
#'
#' @param x an object of class \code{mcmc.list}, returned by the
#' function \code{coda.samples} from the package \code{rjags}.  The
#' model string used to fit the model should have been generated by
#' the function \code{prepareFit}.
#'
#'
#' @return a list with two elements containing: (i) \code{variables}:
#' a data.frame containing the probability of inclusion of each
#' variable in the model, (ii) \code{models}: a data.frame containing
#' the probability that each possible model (i.e. a given combination
#' of the variables used in the model) is true. In both cases, the
#' variables and models are sorted by decreasing value of probability.
#'
#' @seealso \code{\link{prepareFit}}, \code{\link[rjags]{coda.samples}}
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
#' @examples
#'
#' data("modelRedDeer")
#' probabilityKM(modelRedDeer)
#'
probabilityKM <- function(x)
{
    if (!inherits(x, "mcmc.list"))
        stop("x should be of class mcmc.list")
    if (length(grep("alpha_", colnames(x[[1]])))==0)
        stop("no alpha coefficients detected in x.\nThis model has not been fitted using models created\nby the function prepareFit")
    xb <- do.call(rbind,x)
    results <- list()
    probas <- sort(apply(xb[,grep("alpha", colnames(xb))],2,mean),decreasing=TRUE)
    nam <- gsub("alpha_","",names(probas))
    df <- data.frame(variable=nam,proba=probas)
    row.names(df) <- 1:nrow(df)
    results$variables <- df

    re <- xb[,grep("alpha", colnames(xb))]
    na <- colnames(xb[,grep("alpha", colnames(xb))])
    ta <- table(apply(re,1, function(x)
        paste(gsub("alpha_","",na[unlist(x)==1]), collapse="+")))
    ta <- ta/sum(ta)
    ta <- sort(ta, decreasing=TRUE)
    dfb <- data.frame(model=names(ta), proba=as.vector(ta))
    row.names(dfb) <- 1:nrow(dfb)
    results$models <- dfb

    return(results)
}




#' @title Bayesian Model To Identify Factors Affecting Wildlife-Vehicle Collisions
#' @export
#' @aliases finalModel
#' @details \code{finalModel} fits the final Bayesian model used to
#' predict the number of collisions between ungulates and vehicle as a
#' function of a linear combination of a set of environmental
#' variables.
#'
#' @param X a data.frame containing the numeric variables supposed to
#' have an effect on the wildlife-vehicle collisions (columns) for
#' spatial unit (rows).
#' @param vectorFinalVariables a vector of character strings
#' containing the names of the final variables used in the
#' combination.
#' @param collisions an integer vector with length equal to
#' \code{nrow(X)} defining the number of wildlife-ungulate collisions
#' in each studied spatial unit.
#' @param nYear an integer vector with length equal to \code{nrow(X)}
#' defining the number of Years during which wildlife-ungulate
#' collisions have been recorded in each studied spatial unit.
#' @param Area a numeric vector with length equal to \code{nrow(X)}
#' defining the area of each studied spatial unit.
#' @param departement a character vector with length equal to
#' \code{nrow(X)} defining the department of each studied spatial unit.
#'
#'
#' @return a list with all elements required for the fit of the model
#' with JAGS, that is: (i) \code{data4jags}: the list of the data
#' required by the model, to be passed to the argument \code{data} of
#' the function \code{jags.model} of the package \code{rjags}, (ii)
#' \code{ini}: list of starting values for the parameters, to be
#' passed to the argument \code{init} of the function
#' \code{jags.model}, (iii) \code{modelstring}: a character string
#' containing the model fit by JAGS, (iv) \code{coefnames}: vector of
#' character strings containing the names of the coefficients of
#' interest in the model, to be passed to the argument
#' \code{variable.names} of the function \code{coda.samples} of the
#' package \code{rjags}.
#'
#' @seealso \code{\link[rjags]{jags.model}},
#' \code{\link[rjags]{coda.samples}}, \code{\link{prepareFit}}
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
#' @import stats
#' @examples
#'
#' ## The data used for the fit
#' data("dataCollision")
#'
#' ## Consider the Kuo-Mallik model fit (see ?prepareFit to see how
#' ## this model was fitted)
#' data("modelRedDeer")
#'
#' ## we consider the probability associated to each model
#' probabilityKM(modelRedDeer)
#'
#' ## The best model seems to be the model including forFrag and hunt
#'
#' ## Prepare the data
#' X <- dataCollision$RedDeer$X
#'
#' ## We scale the variables to improve mixing
#' X <- as.data.frame(scale(X))
#'
#' ## Prepare the final model, including only the variables forFrag and hunt
#' pfm <- finalModel(X, c("forFrag","hunt"), dataCollision$RedDeer$coll,
#'                   dataCollision$RedDeer$Y, dataCollision$RedDeer$Area,
#'                   dataCollision$RedDeer$departement)
#'
#' \dontrun{
#' ## WARNING: long execution (about 10 min)!!
#' ## the results are stored in the dataset "finalModelRedDeer"
#' ## But if you want to try it, the data are now ready for the fit:
#'
#' mo <- jags.model(textConnection(pfm$modelstring), ini=pfm$ini, data=pfm$data4jags)
#' update(mo, n.iter=1000)
#' finalModelRedDeer <- coda.samples(mo, variable.names = pfm$coefnames,
#'                                   n.iter = 500000, thin=100)
#'
#' }
#'
#'
finalModel <- function(X, vectorFinalVariables, collisions, nYear, Area, departement)
{
    if (!inherits(X, "data.frame"))
        stop("X should be a data.frame")
    if (nrow(X)!=length(collisions))
        stop("The number of rows of X should be equal to length(collisions)")
    if (nrow(X)!=length(nYear))
        stop("The number of rows of X should be equal to length(nYear)")
    if (nrow(X)!=length(Area))
        stop("The number of rows of X should be equal to length(Area)")
    if (nrow(X)!=length(departement))
        stop("The number of rows of X should be equal to length(departement)")
    if (!is.character(departement))
        stop("departement should be character")
    if (!all(vectorFinalVariables%in%names(X)))
        stop("some variables in vectorFinalVariables not in X")

    ## Data preparation
    X <- X[vectorFinalVariables]
    data4jags <- lapply(X, function(x) x)

    names(data4jags) <- toupper(names(data4jags))

    data4jags$collisions <- collisions
    data4jags$departement <- as.numeric(factor(departement))
    data4jags$nYear <- as.numeric(factor(nYear))
    data4jags$Area <- as.numeric(factor(Area))
    data4jags$ndep <- length(unique(departement))
    data4jags$nUG <- nrow(X)

    ##
    ##
    ## Starting values

    ## starting values for the parameters beta. Random draw
    ## from the normal distribution
    beta <- paste("beta", names(X), sep="_")
    ini <- lapply(beta, function(x) rnorm(1))
    names(ini) <- beta

    ## formula describing the core of the model
    formcoeur <- paste(paste(paste("beta", names(X), sep="_"),
                             paste0(toupper(names(X)), "[i]"), sep="*"),
                       collapse="+")

    ##
    ##
    ## Model Code

    formulalik <- paste0("lambda[i] <- nYear[i]*Area[i]*exp(",
                         "efal[departement[i]]+",
                         formcoeur,
                         "+rs[i])")

    priors <- c("model {",
                "## Priors",
                "## Coefficients beta (importance in the model)",
                paste0(beta, "~dnorm(0.0, 0.01)"),
                "\n",
                "## precision of random effects departement",
                "sigefa~dgamma(0.01, 0.01)",
                "muefa~dnorm(0, 0.01)",
                "\n",
                "## Random effect departement",
                "for (d in 1:ndep) {",
                "efal[d]~dnorm(muefa, sigefa)",
                "}",
                "\n",
                "## precision of overdispersion residuals",
                "sigeps~dgamma(0.01, 0.01)",
                "\n",
                "## Likelihood",
                "for (i in 1:nUG) {",
                "rs[i]~dnorm(0,sigeps)",
                formulalik,
                "collisions[i]~dpois(lambda[i])",
                "}",
                "\n",
                "}\n")
    modelstring <- (paste(priors, collapse = "\n"))

    ## Coefficient names
    coefnames <- c(beta, "sigeps", "sigefa", "muefa")

    ## Resultats
    resu <- list(data4jags=data4jags,
                 ini=ini,
                 modelstring=modelstring,
                 coefnames=coefnames)

    ## output
    return(resu)

}
