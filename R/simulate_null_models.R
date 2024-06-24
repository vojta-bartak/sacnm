require(geoR)
require(RandomFields)

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
#' @param output There will be an explanation here...
#' @param formulas There will be an explanation here...
#' @param radius A maximum distance for random shift of the data points. Only used when method is 'shift' or 'shift_only'
#' (see Details).
#' @param nsim The number of simulated null models to return.
#' @return A list of summaries of simulated models. The first item correspond to the original model, the rest correpond to the
#' models fitted to simulated data. The length of the list is nsim + 1.
#'
#' Here, the description of the individual list items will be included.
#' @export
simulate_null_models <- function(model, data, preds = NULL, pred_ras = NULL, variog = NULL, coords = c('x','y'),
                                 method = c('shift','RFsim','Viladomat','kriging','shift_only','rotate_only',
                                            'rotate_warp','kriging_warp'),
                                 output = c('coef','p','dev','AIC','AUC','R2','MSE'),
                                 formulas = NULL,
                                 radius = NULL,
                                 nsim = 1000,
                                 var_model = "mat",
                                 kappa = 1.5,
                                 fixed_nugget = TRUE,
                                 nugget = 0,
                                 center_data = TRUE)
{
  method=method[1]
  if (method %in% c('shift','kriging','shift_only') & class(radius)!="numeric")
    stop("Error: radius must be specified for the selected method!")
  if (method %in% c('shift','shift_only','rotate_only')) {
    if (inherits(pred_ras, "list")) {
      if (!min(sapply(pred_ras, function(ras) inherits(ras, "stars"))))
        stop("Error: pred_ras is missing or wrong but needed for the selected method!")
    }
    else if (!inherits(pred_ras, "stars"))
      stop("Error: pred_ras is missing or wrong but needed for the selected method!")
  }
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

  prctile = NULL
  nuggets = NULL
  if (method %in% c('shift','shift_only','rotate_only','rotate_warp',"kriging_warp")) preds <- preds[preds %in% names(pred_ras)]
  else if (method == "kriging"){
    require(spaMM)
    variog <- lapply(preds, function(pred) fitme(as.formula(paste(pred,"~1+Matern(1|",coords[1],"+",coords[2],")", sep="")), data=data))
    names(variog) <- preds
  } else if (method == "Viladomat"){
    require(geoR)
    prctile <- quantile(dist(data[,coords]), probs = 0.25)
    variog <- lapply(preds, function(pred)
      variog(data = data[,pred], coords = data[,coords], max.dist = prctile, option = "bin", messages = FALSE)
      )
    names(variog) <- preds
  } else if (method == "RFsim"){
    require(geoR)
    nuggets <- numeric(0)
    variog <- list()
    for (pred in preds){
      vgm <- likfit(data=data[,pred],
                    coords = cbind(data[,coords[1]], data[,coords[2]]),
                    cov.model = var_model, kappa = kappa,
                    ini.cov.pars = c(var(data[,pred]), max(max(data[,coords[1]])-min(data[,coords[1]]),
                                                           max(data[,coords[2]])-min(data[,coords[2]]))),
                    fix.nugget = fixed_nugget, nugget = nugget,
                    messages = F)
      nuggets <- c(nuggets, if (!fixed_nugget) vgm$nugget else nugget)
      variog <- c(variog, switch(var_model,
                                 "mat" = RMwhittle(nu=kappa, var = vgm$cov.pars[1], scale = vgm$cov.pars[2])))
    }
    names(variog) <- preds
    names(nuggets) <- preds
  }

  output <- lapply(1:nsim, function(i){
    simdata <- simulate_data(data=data, preds=preds, coords=coords, pred_ras=pred_ras, variog=variog, method=method, radius=radius,
                             prctile=prctile, nuggets = nuggets, center_data = center_data)
    newdata <- simdata$newdata
    if (inherits(model, "ranger")) newdata <- na.omit(newdata)
    newmodel <- try(update(model, data=newdata), silent = T)
    if (inherits(newmodel, "try-error")) {NULL} else {
      list(model = newmodel, indices = simdata$indices)
    }
  })

  output <- append(output, list(list(model=model, indices = 1:nrow(data))), 0)
  output
}

#' Summarize null models
#'
#' Description will be here...
#'
#' @param null_models ...
#' @param data ...
#' @return ...
#'
#' @export
summarize_null_models <- function(null_models, data){
  lapply(null_models, function(model) summarize_model(model, data))
}


summarize_model <- function(model_where, data){
  model = model_where$model
  where = model_where$indices
  preds <- if (!inherits(model, "ranger")) all.vars(delete.response(terms(model)))
  else strsplit(strsplit(as.character(model$call)[2], " ~ ")[[1]][2], " + ", fixed=T)[[1]]
  output <- list(
    coefs = coef(model),
    deviance = deviance(model),
    mse = NULL,
    r2 = NULL,
    AIC = NULL,
    d2 =NULL,
    preds.p = NULL,
    preds.dev = NULL,
    preds.AIC = NULL,
    preds.r2 = NULL,
    preds.importanceRF = NULL
  )
  if (inherits(model, "ranger")){
    output$mse <- model$prediction.error
    output$r2 <- model$r.squared
    output$preds.importanceRF <- model$variable.importance
  } else if (inherits(model, "gam")) {
    s <- summary(model)
    output$r2 <- s$r.sq
    output$d2 <- s$dev.expl
    output$AIC <- AIC(model)
  } else if (inherits(model, "glm")) {
    s <- summary(model)
    output$r2 <- 1 - s$deviance/s$null.deviance
    output$d2 <- output$r2
    output$AIC <- AIC(model)


  } else if (inherits(model, "lm")){
    s <- summary(model)
    output$mse <- s$sigma
    output$r2 <- s$r.squared
    output$d2 <- output$r2
    output$AIC <- AIC(model)
  }

  # if (inherits(model, "lm") & !inherits(model, "gam")){
  #   terms <- attr(terms(model),"term.labels")
  #   updated.models <- lapply(preds, function(p) {
  #     formula_update <- paste("~.-", paste(terms[grepl(p, terms)], collapse = "-"), sep="")
  #     update(model, formula_update, data=data[where,])
  #   })
  #   anovas <- lapply(updated.models, function(updated.model) {
  #     anova(model, updated.model, test="LR")
  #   })
  #   output$preds.p <- sapply(anovas, function(a) a[2,5])
  #   names(output$preds.p) <- preds
  #   output$preds.dev <- sapply(anovas, function(a) a[2,4])
  #   names(output$preds.dev) <- preds
  #   output$preds.AIC <- sapply(updated.models, function(mod) AIC(model) - AIC(mod))
  #   names(output$preds.AIC) <- preds
  #   if (inherits(model, "glm")){
  #     output$preds.r2 <- sapply(updated.models, function(mod) {
  #       s <- summary(mod)
  #       output$r2 - 1 + s$deviance/s$null.deviance
  #     })
  #   } else {
  #     output$preds.r2 <- sapply(updated.models, function(mod) {
  #       s <- summary(mod)
  #       output$r2 - s$r.squared
  #     })
  #   }
  #
  #   names(output$preds.r2) <- preds
  # }

  output
}




# Simulate new data
#' @export
simulate_data <- function(data, preds, coords=c('x','y'), pred_ras=NULL, variog=NULL,
                          method=c('shift','RFsim','Viladomat','kriging','shift_only','rotate_only','rotate_warp','kriging_warp'),
                          radius=NULL, prctile=NULL, nuggets = NULL, center_data = TRUE)
{
  preds <- preds[preds %in% colnames(data)]
  if (method %in% c("shift", "shift_only", "rotate_only")) require(stars)

  if (method=='shift') {
    newdata <- shift_rotate(data, coords, radius)
    for (pred in preds){
      newdata[,pred] <- st_extract(pred_ras[[pred]], cbind(newdata[,coords[1]],newdata[,coords[2]]))[[1]]
    }
  } else if (method=='shift_only') {
    newdata <- shift(data, coords, radius)
    for (pred in preds){
      newdata[,pred] <- st_extract(pred_ras[[pred]], cbind(newdata[,coords[1]],newdata[,coords[2]]))[[1]]
    }
  } else if (method=='rotate_only') {
    newdata <- rotate(data, coords)
    for (pred in preds){
      newdata[,pred] <- st_extract(pred_ras[[pred]], cbind(newdata[,coords[1]],newdata[,coords[2]]))[[1]]
    }
  } else if (method=='rotate_warp') {
    newdata <- rotate_warp(data, coords, radius = radius)
    for (pred in preds){
      newdata[,pred] <- st_extract(pred_ras[[pred]], cbind(newdata[,coords[1]],newdata[,coords[2]]))[[1]]
    }
  } else if (method=='RFsim') {
    newdata <- data
    for (pred in preds){
      newdata[,pred] <- RFsimulate(variog[[pred]], x=data[,coords[1]], y=data[,coords[2]])$variable1 +
        (if (nuggets[pred] > 0) RFsimulate(RMnugget(var=nuggets[pred]))$variable1 else 0)
    }
  } else if (method=='Viladomat') {
    newdata <- data
    perm <- sample(1:nrow(data), size = nrow(data), replace = F)
    for (pred in preds) {
      newdata[,pred] <- variog.matching(data[perm, pred],
                                        coords=data[,coords],
                                        Delta = seq(0.1,0.9,0.1),
                                        target_variog = variog[[pred]],
                                        prctile = prctile)
    }

  } else if (method=='kriging') {
    newdata <- shift_rotate(data, coords, radius)
    for (pred in preds) newdata[,pred] <- predict(variog[[pred]], newdata=newdata)[,1]
  } else if (method=='kriging_warp') {
    newdata <- rotate_warp(data, coords, radius = radius)
    for (pred in preds) newdata[,pred] <- predict(variog[[pred]], newdata=newdata)[,1]
  }
  indices <- which(apply(newdata, 1, function(row) max(is.na(row))) == 0)
  if (center_data) for (pred in preds) newdata[,pred] <- newdata[,pred] + mean(data[,pred], na.rm=T) - mean(newdata[,pred], na.rm=T)
  for (pred in preds) newdata <- newdata[!is.na(newdata[,pred]),] # na.omit instead?
  list(newdata=newdata, indices=indices)
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

# Data random rotation and warp
rotate_warp <- function(data, coords=c("x","y"), radius=NULL, angle=NULL, distance=NULL){
  if (is.null(radius)) radius <- diff(range(data[,coords[1]]))
  if (is.null(angle)) angle <- runif(1, 0, 2*pi)
  x <- data[,coords[1]]
  y <- data[,coords[2]]
  mx <- mean(data[,coords[1]])
  my <- mean(data[,coords[2]])
  x.rot <- (x-mx)*cos(angle) - (y-my)*sin(angle) + mx
  y.rot <- (x-mx)*sin(angle) + (y-my)*cos(angle) + my
  if (is.null(distance)) distance <- runif(1, 0, radius)
  newdat <- data
  newdat[,coords[1]] <- x.rot
  newdat[,coords[2]] <- ifelse(y.rot + distance > max(newdat[,coords[2]]),
                               min(newdat[,coords[2]]) + y.rot + distance - max(newdat[,coords[2]]),
                               y.rot + distance)
  newdat
}

# Update model - fit model with simulated data
# update_model <- function(model, newdata){
#   if (inherits(model, "ranger")){
#
#   }
#   else if (inherits(model, "gam"))
# }

variog.matching <- function(x, coords, Delta, target_variog, prctile) {
  require(geoR)
  require(locfit)
  variog.X.delta <- vector(mode = "list", length = length(Delta))
  linear.fit <- variog.X.delta
  hat.X.delta <- variog.X.delta
  resid.sum.squares <- rep(0, length(Delta))
  for (k in 1:length(Delta)) {
    # smooth X.randomized using locfit:
    fit <- locfit(x ~ lp(coords[,1], coords[,2], nn = Delta[k], deg = 0), kern = "gauss", maxk = 300)
    X.delta <- fitted(fit)

    # variogram of X.delta:
    variog.X.delta[[k]] <- variog(data = X.delta, coords = coords, option = "bin", max.dist = prctile, messages = FALSE)

    # linear regression between the target and X.delta variograms:
    linear.fit[[k]] <- lm(target_variog$v ~ 1 + variog.X.delta[[k]]$v)

    # least square estimates:
    bet.hat <- as.numeric(linear.fit[[k]]$coefficients)

    # transformed X.delta:
    hat.X.delta[[k]] <- X.delta * sqrt(abs(bet.hat[2])) + rnorm(length(X.delta)) * sqrt(abs(bet.hat[1]))
    variog.hat.X.delta <- variog(data = hat.X.delta[[k]], coords = coords, option = "bin", max.dist = prctile, messages = FALSE)

    # hat.X.delta[[k]] <- X.delta
    # variog.hat.X.delta <- variog(data = X.delta, coords = coords, option = "bin", max.dist = prctile, messages = FALSE)

    # sum of squares of the residuals:
    resid.sum.squares[k] <- sum((variog.hat.X.delta$v-target_variog$v) ^ 2)
  }

  # delta that minimizes the residual sum of squares:
  delta.star.id <- which.min(resid.sum.squares)
  hat.X.delta.star <- hat.X.delta[[delta.star.id]]
  # print(resid.sum.squares[delta.star.id])

  return(hat.X.delta.star)
}
