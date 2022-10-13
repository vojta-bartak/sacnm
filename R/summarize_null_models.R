#' Model coefficients inference table based on a given set of null model objects
#'
#' This function summarizes the null distribution of model coefficients based on a given set of null model objects
#' (typically created by the function simulate_null_models from teh same package) and compares the observed value of
#' each coefficient with these null distributions.
#'
#' @param model A fitted model object. Currently, objects of classes 'lm', 'glm', and 'gam' are supported.
#' @param null_models A list of model objects of the same class as model representing simulations of a null model.
#' @param coefs An optional character vector specifying for which model coefficients the inference should be made.
#' If not specified,the inference is made for all coefficients.
#' @return A data frame summarizing the null distribution for each model coefficient, including a Monte Carlo p-value.
#' @export
coef_tab <- function(model, null_models, coefs=NULL)
{
  if (is.null(coefs)) coefs <- names(coef(model))
  obs.coef <- coef(model)[coefs]
  null.coef <- do.call(rbind, lapply(null_models, function(m) coef(m)[coefs]))
  sums <- as.data.frame(t(apply(null.coef, 2, quantile, probs=c(0,0.25,0.5,0.75,1), na.rm=T)))
  colnames(sums) <- c("Min.","1st Qu.", "Median", "3rd Qu.", "Max.")
  sums$Mean <- apply(null.coef, 2, mean, na.rm=T)
  sums$Observed <- obs.coef
  sums$P.value <- sapply(names(obs.coef), function(prm){
    p <- rank(c(obs.coef[prm],null.coef[,prm]))[1]/(nrow(null.coef)+1)
    2*min(p, 1-p)
  })
  sums$Significance <- sapply(sums$P.value, function(p){
    if(p<0.001){"***"}else if (p<0.01){"**"}else if (p<0.05){"*"}else if(p<0.1){"."}else{""}
  })
  sums
}

#' Model coefficients inference plots based on a given set of null model objects
#'
#' This function plots the observed value of each model coefficient and compares it with a null distribution of this
#' coefficient based on a set of null model objects. These null models are typically created by the function
#' simulate_null_models from the same package. The plots are based on (and thus require) ggplot2 package.
#'
#' @param model A fitted model object. Currently, objects of classes 'lm', 'glm', and 'gam' are supported.
#' @param null_models A list of model objects of the same class as model representing simulations of a null model.
#' @param coefs An optional character vector specifying which model coefficients should be inferred. If not specified,
#' the inference is made for all model coefficients.
#' @return A ggplot with a facet for each model coefficient, plotting its observed value (red line) on a density of
#' its null distribution as well as 0.025, 0.5, and 0.975 quantiles.
#' @export
coef_plot <- function(model, null_models, coefs=NULL)
{
  require(ggplot2)
  if (is.null(coefs)) coefs <- names(coef(model))
  obs.coef <- as.data.frame(coef(model)[coefs])
  colnames(obs.coef) <- "value"
  obs.coef$name <- rownames(obs.coef)
  obs.coef$quant <- "Observed"
  null.coef <- as.data.frame(do.call(rbind, lapply(null_models, function(m) coef(m)[coefs])))
  ylims <- apply(null.coef, 2, function(col) max(stats::density(col)$y))
  obs.coef$ylim <- ylims
  ci <- t(apply(null.coef, 2, function(col) quantile(col, c(0.025,0.5,0.975))))
  ci <- rbind(data.frame(value=ci[,1],name=rownames(ci),quant=colnames(ci)[1],ylim=ylims),
              data.frame(value=ci[,2],name=rownames(ci),quant=colnames(ci)[2],ylim=ylims),
              data.frame(value=ci[,3],name=rownames(ci),quant=colnames(ci)[3],ylim=ylims))
  ci <- rbind(ci, obs.coef)
  null.coef <- do.call(rbind, lapply(names(null.coef), function(col){
    data.frame(value=null.coef[,col], name=col)
  }))
  ggplot(null.coef, aes(x=value)) +
    geom_density(fill="grey", alpha=.3) +
    geom_segment(aes(x=value, xend=value, y=0, yend=ylim, lty=quant, color=quant), data=ci) +
    facet_wrap(~name, scales="free") +
    scale_linetype_manual(values = c(2,1,2,1)) +
    scale_color_manual(values=c("black","black","black","red")) +
    labs(y="Density", linetype="Quantiles", color="Quantiles")
}

#' Plotting the effects of individual continuous predictors
#'
#' This function plots the effects of individual continuous predictors, keeping the values of other predictors at their
#' mean values, and showing 95% envelope based on given null model simulations. It can be used for a visual check whether
#' the observed effect could be explained by the null model or not. Typically, the null model relaxes any structural
#' relationships between the response and the predictors while keeping the autocorrelation structure unchanged.
#'
#' @param model A fitted model object. Currently, objects of classes 'lm', 'glm', and 'gam' are supported.
#' @param null_models A list of model objects of the same class as model representing simulations of a null model.
#' @param data A data frame used to fit the model.
#' @param preds An optional character vector of names of predictors to be plotted. If not specified, all continuous
#' predictors used to fit the model are used.
#' @param lower The lower quantile of the simulation envelope to be plotted. Default is 0.025.
#' @param upper The upper quantile of the simulation envelope to be plotted. Default is 0.975.
#' @param nval The number of values (between minimum and maximum from its observed values) of each predictor to be used
#' for plotting the effect. Default is 100.
#' @param plotdata A logical indicating whether the data points should be plotted together with the effect. Default is TRUE.
#' @return A ggplot object with a facet for each model predictor, plotting the predicted values (red line) together with
#' simulation envelopes (grey shaded area) and median (black line).
#' @export
effect_plot <- function(model, null_models, data, preds=NULL, lower=c(0.025), upper=c(0.975), nval=100, plotdata=TRUE){
  require(ggplot2)
  null_models <- Filter(function(x) !is.null(x), null_models)
  if (inherits(model, "ranger")){
    if (is.null(preds))  {
      preds <- strsplit(as.character(model$call)[2], split=" ~ ", fixed=T)[[1]]
      preds <- preds[2:length(preds)]
      preds <- preds[class(data[,preds])=="numeric"]
    }
    df <- do.call(rbind, lapply(preds, function(pred){
      newdata <- prepare_data(data, pred, nval=nval)
      predictions <- do.call(cbind,lapply(null_models, function(m){
        predict(m, data=newdata)$predictions
      }))
      data.frame(
        x = newdata[,pred],
        y = predict(model, data=newdata)$predictions,
        mean = apply(predictions, 1, mean),
        lwr = apply(predictions, 1, quantile, probs=lower),
        upr = apply(predictions, 1, quantile, probs=upper),
        predictor = pred
      )
    }))
  } else {
    if (is.null(preds)) {
      preds <- attr(terms(model), "term.labels")
      preds <- preds[class(data[,preds])=="numeric"]
    }
    df <- do.call(rbind, lapply(preds, function(pred){
      newdata <- prepare_data(data, pred, nval=nval)
      predictions <- do.call(cbind,lapply(null_models, function(m){
        predict(m, newdata=newdata, type="response")
      }))
      data.frame(
        x = newdata[,pred],
        y = predict(model, newdata=newdata, type="response"),
        mean = apply(predictions, 1, mean),
        lwr = apply(predictions, 1, quantile, probs=lower),
        upr = apply(predictions, 1, quantile, probs=upper),
        predictor = pred
      )
    }))
  }
  resp <- strsplit(as.character(model$call)[2], split=" ~ ", fixed=TRUE)[[1]][1]
  p <- ggplot(df, aes(x=x, y=mean)) +
    geom_line() +
    geom_ribbon(aes(ymin=lwr, ymax=upr), fill="grey", alpha=0.3) +
    geom_line(aes(y=y), color="red") +
    facet_wrap(~predictor, scales = "free") +
    labs(x="predictor values", y=resp)
  if (plotdata) {
    datadf <- do.call(rbind, lapply(preds, function(pred){
      data.frame(
        x=data[,pred],
        mean=data[,resp]
      )
    }))
    p <- p + geom_point(data=datadf, alpha=.1)
  }
  p
}

#' Function to prepare data for the effect plot, preparing the sequence of values of the predictor pred
#' while keeping the values of other continuous predictors at their mean value.
prepare_data <- function(data, pred, nval){
  data[,pred] <- seq(min(data[,pred]), max(data[,pred]), l=nval)
  for (p in names(data)[names(data)!=pred]) if (class(data[,p])=="numeric") data[,p] <- mean(data[,p])
  data
}
