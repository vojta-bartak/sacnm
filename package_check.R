library(RandomFields)
library(spatstat)
library(sf)
library(stars)
library(tidyverse)
library(tmap)
library(gstat)
library(car)
library(ranger)
library(mgcv)

# predictor parameters
scale.pred <- 50
var.pred <- 1
kappa.pred <- 1.5
mod.pred <- RMwhittle(nu=kappa.pred, var=var.pred, scale=scale.pred)

# response parameters
scale.resp <- 50
var.resp <- 1
kappa.resp <- 1.5
alpha <- 1
mod.resp <- RMwhittle(nu=kappa.resp, var=var.resp, scale=scale.resp)

# plot rasters
r.pred <- RFsimulate(mod.pred, x=0:100, y=0:100, grid=T) %>% st_as_stars
r2.pred <- RFsimulate(mod.pred, x=0:100, y=0:100, grid=T) %>% st_as_stars
r.resp <- r.pred + alpha*(RFsimulate(mod.resp, x=0:100, y=0:100, grid=T) %>% st_as_stars)
tmap_arrange(
  tm_shape(r.pred) + tm_raster(style="cont"),
  tm_shape(r.resp) + tm_raster(style="cont")
)

# data
x <- runif(100, 0, 100)
y <- runif(100, 0, 100)
p1 <- st_extract(r.pred, cbind(x,y))$variable1
p2 <- st_extract(r2.pred, cbind(x,y))$variable1
response <- p1 + alpha*st_extract(r.resp, cbind(x,y))$variable1
df <- data.frame(x=x, y=y, p1=p1, response=response, p2=p2)
plot(p1,response)
plot(p2,response)
cor.test(p1, response)
cor.test(p2, response)

source("R/simulate_null_models.R")
source("R/summarize_null_models.R")

plot(df %>% st_as_sf(coords=c("x","y")) %>% st_geometry)
df2 <- df %>% shift_rotate
df3 <- df %>% shift_rotate
plot(df2 %>% st_as_sf(coords=c("x","y")) %>% st_geometry, add=T, col="red")
plot(df3 %>% st_as_sf(coords=c("x","y")) %>% st_geometry, add=T, col="green")
df2$p1 <- predict(mm, newdata=df2)[,1]
df3$p1 <- predict(mm, newdata=df3)[,1]

plot(df2$predictor, df3$predictor)

coef(ms[[1]])
coef(ms[[2]])

# models ==============================================================================================================
df$resp <- as.numeric(df$response > -3)
table(df$resp)
m.lm <- lm(response~p1+p2, data=df)
m.glm <- glm(resp~p1 + p2, data=df, family=binomial)
m.gam <- gam(response~s(p1) + s(p2), data=df)
m.rf <- ranger(response~p1+p2, data=df)

m <- m.rf
summary(m)
ms <- simulate_null_models(m, df, radius = 100, pred_ras=list(predictor=r.pred), method="kriging")
sums <- summarize_null_models(ms, df)
coef_tab(sums)
pred_tab(sums)
coef_plot(m, sums)
effect_plot(m, ms, df)




ms.rf <- simulate_null_models(m.rf, df, radius = 100, pred_ras=list(predictor=r.pred), method="kriging")

table(sapply(ms, is.null))

coef_tab(sums)
pred_tab(sums)
coef_plot(m, sums)

effect_plot(m, ms, df)
effect_plot(m, ms, df, preds=c("p1"))

params <- expand.grid(method=c("kriging","shift"), radius=c(10,50,100), scale=c(1,10,50))
plots <- lapply(1:nrow(params), function(row){
  print(params[row,])
  method=params[row,"method"]
  scale=params[row,"scale"]
  radius=params[row,"radius"]
  mod.pred <- RMwhittle(nu=1.5, var=1, scale=scale)
  mod.resp <- RMwhittle(nu=1.5, var=1, scale=scale)
  r.pred <- RFsimulate(mod.pred, x=0:100, y=0:100, grid=T) %>% st_as_stars
  r.resp <- r.pred + (RFsimulate(mod.resp, x=0:100, y=0:100, grid=T) %>% st_as_stars)
  x <- runif(100, 0, 100)
  y <- runif(100, 0, 100)
  predictor <- st_extract(r.pred, cbind(x,y))$variable1
  response <- st_extract(r.resp, cbind(x,y))$variable1
  df <- data.frame(x=x, y=y, predictor=predictor, response=response)
  m <- lm(response~predictor, data=df)
  ms <- simulate_null_models(m, df, radius = radius, pred_ras=list(predictor=r.pred), method=method)
  effect_plot(m, ms, df) +
    labs(title=paste("Method = ",method,", Radius = ",radius,", Scale = ",scale,sep=""))
})
library(gridExtra)
grid.arrange(grobs=plots, ncol=6)
