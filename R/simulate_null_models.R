require(geoR)
require(RandomFields)
require(stars)

#' Simulate null models
#'
#' This is the main function, simulating a set of null models based on a given model and data, keeping the original
#' spatial autocorrelation structure unchanged, but relaxing any original structural relationship between the response
#' and predictors.
#'
#' The method 'shift' randomly move the data points (all in the same way) by first rotating the dataset around the data
#' centroid by a random angle, and then shifting them by in a random direction by a distance generated as a uniform random
#' variable on a (0, radius) distance. Then, new predictor values are extracted at these new locations from the corresponding
#' predictor raster passed through the argument pred_ras. Methods 'shift_only' and 'rotate_only' perform only the corresponding
#' part of the procedure (i.e., either random shift or random rotation).
#'
#' The method 'RFsim' generate the new datasets by unconditional simulation of random fields with given variogram. the simulation
#' is performed using the function RFsimulate from the package RandomFields (Schlather et al. 2022). The variogram is either
#' known and passed through the variog argument, or unknown and estimated from data using maximum likelihood method (function
#' likfit from the geoR (Ribeiro et al. 2020) package).
#'
#' The method 'kriging' is based on first random rotation and shift of the data points (same as for the method 'shift'), and then
#' estimation of the new values at these locations by performing kriging. Similarly to the method RFsim, if the argument variog is
#' passed (i.e., is not NULL), this variogram is used for the kriging. If variog is NULL, the variogram needed for the kriging
#' is estimated from data by maximum likelihood method.
#'
#' The method 'Viladomat' is based on random permutation of the predictor values and then smoothing the values so that their
#' variogram approximately matches the original variogram. the method is described in Viladomat et al. 2014.
#' @param model A fitted model object. Currently, objects of classes 'lm', 'glm', 'gam', and 'ranger' are supported.
#' @param data A data frame with the original dataset used to fit the model, including spatial coordinates of the data points.
#' @param preds A character vector of names of predictors, for which the null models should be generated. If NULL, all predictors
#' from the model are used (except coordinates themselves, in case they were used as predictors). Note that for some methods,
#' corresponding predictor rasters must be passed through the argument pred_ras (see Details).
#' @param pred_ras A named list of predictor rasters, with one raster for each predictor. (Only one predictor is currently
#' supported). The rasters in the list must be named the same as the original predictors. The argument is optional, only
#' required when method is one of 'shift', 'shift_only', and 'rotate only'.
#' @param variog A variogram object of the class 'variogram' from the package geoR (Ribeiro et al. 2020). An optional argument,
#' ignored unless the method is 'RFsim'. See Details.
#' @param coords A character vector of lenght two with names of the columns in data where spatial coordinates are stored.
#' Defaults are 'x' and 'y'.
#' @param method A method how the null models are generated. Possible values are 'shift','RFsim','Viladomat','kriging',
#' 'shift_only',and 'rotate_only'. See Details for their description. Default is 'shift'.
#' @param radius A maximum distance for random shift of the data points. Only used when method is 'shift' or 'shift_only'
#' (see Details).
#' @param nsim The number of simulated null models to return.
#' @return A list of models of the same class as the input model and the lentgh given by nsim.
#' @export
simulate_null_models <- function(model, data, preds=NULL, pred_ras=NULL, variog=NULL, coords=c('x','y'),
                                 method=c('shift','RFsim','Viladomat','kriging','shift_only','rotate_only'),
                                 radius=NULL,
                                 nsim=1000)
{
  if (method %in% c('shift','kriging','shift_only','rotate_only') & class(radius)!="numeric")
    stop("Error: radius must be specified for the selected method!")
  if (method %in% c('shift','shift_only','rotate_only') & !inherits(pred_ras, "stars"))
    stop("Error: pred_ras is missing or wrong but needed for the selected method!")

  if (class(coords)!="character") {stop("Error: coords must be a character vector!")}
  else if (min(coords %in% colnames(data))==0) {stop("Error: specified coordinate columns cannot be found in data!")}

  if (is.null(preds)) {
    if (inherits(model, "ranger")){
      preds <- strsplit(as.character(model$call)[2], split=" ~ ", fixed=T)[[1]][2]
      preds <- strsplit(preds, split=" + ", fixed=T)[[1]]
      preds <- preds[preds %in% colnames(data) & !(preds %in% coords)]
    } else {
      preds <- attr(model$terms , "term.labels")
      preds <- preds[preds %in% colnames(data) & !(preds %in% coords)]
    }
  }
  if (method %in% c('shift','shift_only','rotate_only')) preds <- preds[preds %in% names(pred_ras)]

  lapply(1:nsim, function(i){
    newdata <- simulate_data(data=data, preds=preds, coords=coords, pred_ras=pred_ras, variog=variog, method=method)
    update(model, data=newdata)
  })
}

# Simulate new data
simulate_data <- function(data, preds, coords=c('x','y'), pred_ras=NULL, variog=NULL,
                          method=c('shift','RFsim','Viladomat','kriging','shift_only','rotate_only'),
                          radius=NULL)
{
  if (method=='shift') {
    newdata <- shift_rotate(data, coords, radius)
    for (pred in preds){
      newdata[,pred] <- st_extract(pred_rasters[[pred]], cbind(newdata[,coords[1]],newdata[,coords[2]]))[[1]]
    }
  } else if (method=='shift_only') {
    newdata <- shift(data, coords, radius)
    for (pred in preds){
      newdata[,pred] <- st_extract(pred_rasters[[pred]], cbind(newdata[,coords[1]],newdata[,coords[2]]))[[1]]
    }
  } else if (method=='rotate_only') {
    newdata <- rotate(data, coords)
    for (pred in preds){
      newdata[,pred] <- st_extract(pred_rasters[[pred]], cbind(newdata[,coords[1]],newdata[,coords[2]]))[[1]]
    }
  } else if (method=='RFsim') {

  } else if (method=='Viladomat') {

  } else if (method=='kriging') {
    newdata <- shift_rotate(data, coords, radius)
  }
  newdata
}

# Data random shift
shift <- function(data, coords=c("x","y"), radius=50, angle=NULL, distance=NULL){
  if (is.null(angle)) angle <- runif(1, 0, 2*pi)
  if (is.null(distance)) distance <- runif(1, 0, radius)
  newdat <- data
  newdat[,coords[1]] <- data[,coords[1]] + distance*cos(angle)
  newdat[,coords[2]] <- data[,coords[2]] + distance*sin(angle)
  newdat
}

# Data random rotation
rotate <- function(data, coords=c("x","y"), angle=NULL){
  if (is.null(angle)) angle <- runif(1, 0, 2*pi)
  x <- data[,coords[1]]
  y <- data[,coords[2]]
  mx <- mean(data[,coords[1]])
  my <- mean(data[,coords[2]])
  newdat <- data
  newdat[,coords[1]] <- (x-mx)*cos(angle) - (y-my)*sin(angle) + mx
  newdat[,coords[2]] <- (x-mx)*sin(angle) + (y-my)*cos(angle) + my
  newdat
}

# Data random rotation and shift
shift_rotate <- function(data, coords=c("x","y"), radius=50, angle1=NULL, angle2=NULL, distance=NULL){
  if (is.null(angle1)) angle1 <- runif(1, 0, 2*pi)
  x <- data[,coords[1]]
  y <- data[,coords[2]]
  mx <- mean(data[,coords[1]])
  my <- mean(data[,coords[2]])
  x.rot <- (x-mx)*cos(angle1) - (y-my)*sin(angle1) + mx
  y.rot <- (x-mx)*sin(angle1) + (y-my)*cos(angle1) + my
  if (is.null(angle2)) angle2 <- runif(1, 0, 2*pi)
  if (is.null(distance)) distance <- runif(1, 0, radius)
  newdat <- data
  newdat[,coords[1]] <- x.rot + distance*cos(angle2)
  newdat[,coords[2]] <- y.rot + distance*sin(angle2)
  newdat
}

# Update model - fit model with simulated data
# update_model <- function(model, newdata){
#   if (inherits(model, "ranger")){
#
#   }
#   else if (inherits(model, "gam"))
# }
