#' Significance of predictors based on simulated null models
#'
#' This function summarizes the null distribution of predictor p values based on a given set of null model objects
#' (typically created by the function simulate_null_models from teh same package) and compares the observed p value
#' with this null distribution.
#'
#' @param model A fitted model object. Currently, objects of classes 'lm', 'glm', and 'gam' are supported.
#' @param null_models A list of model objects of the same class as model representing simulations of a null model.
#' @param preds An optional character vector specifying for which model coefficients the inference should be made.
#' If not specified,the inference is made for all coefficients.
#' @return A data frame summarizing the null distribution for each model coefficient, including a Monte Carlo p-value.
#' @export
drop_one <- function(model, null_models, preds=NULL)
{
  if (is.null(coefs)) coefs <- names(coef(model))
  obs.coef <- coef(model)[coefs]
  null.coef <- do.call(rbind, lapply(null_models, function(m) coef(m)[coefs]))
  sums <- as.data.frame(t(apply(null.coef, 2, quantile, probs=c(0,0.25,0.5,0.75,1), na.rm=T)))
  colnames(sums) <- c("Min.","1st Qu.", "Median", "3rd Qu.", "Max.")
  sums$Mean <- apply(null.coef, 2, mean, na.rm=T)
  sums$Observed <- obs.coef
  sums$P.value <- sapply(names(obs.coef), function(prm){
    if (is.na(obs.coef[prm])) {NA} else {
      p <- rank(c(obs.coef[prm],null.coef[,prm]))[1]/(nrow(null.coef)+1)
      2*min(p, 1-p)
    }
  })
  sums$Significance <- sapply(sums$P.value, function(p){
    if(is.na(p)){NA}else if(p<0.001){"***"}else if (p<0.01){"**"}else if (p<0.05){"*"}else if(p<0.1){"."}else{""}
  })
  sums
}

compare_two <- function(model, null_models, new_formula, measures=c("p","dAIC"), plot=T){
  out <- list()
  if (plot) {
    require(ggplot2)
    require(gridExtra)
    plots <- list()
  }
  reduced_model <- update(model, formula=new_formula)
  a <- anova(model, reduced_model)
  if ("p" %in% measures){
    p <- a[2,6]
    ps <- sapply(null_models, function(nm){
      anova(nm, update(nm, formula=new_formula))[2,6]
    })
    out$p_p <- rank(c(p,ps))[1]/(length(ps)+1)
    if (plot) {
      plots <- append(plots,
                      ggplot(data.frame(x=ps), aes(x=x)) +
                        geom_density() +
                        annotate("segment", x=p, xend=p, y=0, yend=max(stats::density(na.omit(ps))$y)))
    }
  } else if ("dAIC" %in% measures){
    dAIC <- AIC(model) - AIC(reduced_model)
  }
  if (plot) grid.arrange(grobs=plots)
  return(out)
}
