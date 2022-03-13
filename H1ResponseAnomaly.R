# --------------- Code --------------
#  Bioclimatic context of species’ populations determines community stability'
# H.1 Negative covariance in response to temperature anomalies at thermal range edges

# The code runs on a slice of the data as population data for a project must be requested from the 
# Finnish, Catalonian, BMS schemes. Consequently, results will not be identical to those presented in the manuscript

library(tidyverse)
library(rethinking)


# ----------- load and prep data -----------

# load data 
butterflym1<-read.delim("/BioclimaticContextCommunityStability/DataFiles/H1Slice.txt", sep= " ")


# make a tidy list for analysis
makelistsp <- function(dt) {
  # need to set indexes for site and species out here
  dt$year_index <- dt$YEAR - 1999 #  INLA likes indexes 1 to n
  dt$transect_index <- as.integer(as.factor(dt$TRANSECT_ID))
  dt$species_index <- as.integer(as.factor(dt$SPECIES_NAME))
  dt$country_index <- as.integer(as.factor(dt$site))
  
  dt.list1 <- list(
    # integers
    ndata = nrow(dt),
    ntransect = length(unique(dt$TRANSECT_ID)),
    nspecies = length(unique(dt$SPECIES_NAME)),
    nyear = length(unique(dt$YEAR)),
    
    # outcome
    n_change =  dt$n_change,
    
    # indexes
    transect =  dt$transect_index,
    year = dt$year_index,
    species = dt$species_index,
    country = dt$country_index,
    
    # main coefficients
    lognminus1 = dt$lognminus1,
    
    span = dt$siteanomalyspring,
    span_2 = dt$siteanomalyspring_sq,
    span_3 = dt$siteanomalyspring_cu,
    
    rng = dt$rng,
    
    # interactions
    span_in = dt$rng_span,
    span_in_2 = dt$rng_span_sq,
    span_in_3 = dt$rng_span_cu,
    # space
    xcor = dt$transect_lon,
    ycor = dt$transect_lat
    
  )
  

}

# turn into lists 
butterfly.dt<-makelistsp(butterflym1)


# ---------- Fitting non spatial models -------------

# load in packages 
library(INLA)
library(brinla)

# code download https://www.r-inla.org/learnmore/books
source('/INLAbookfunc/R/spde-book-functions.R')

# simple model for baseline prior on random effects
lmod <- lm(n_change ~ 1 + lognminus1+ span +span_2 + span_3  + rng + span_in + span_in_2 + span_in_3, data=butterfly.dt)
summary(lmod)

# priors
sdres <- sqrt(var(residuals(lmod)))
rand_prior <- list(prec = list(prior="pc.prec", param = c(3*sdres,0.01))) 

# specify random slopes
butterfly.dt$species.0 = butterfly.dt$species
butterfly.dt$species.1 = butterfly.dt$species
butterfly.dt$species.2 = butterfly.dt$species
butterfly.dt$species.3 = butterfly.dt$species
butterfly.dt$species.4 = butterfly.dt$species
butterfly.dt$species.5 = butterfly.dt$species
butterfly.dt$species.6 = butterfly.dt$species

# fit model (Note: more complex models with increased numbers of random terms were fit during development, this is the model without any random effects estimated close to 0)
inla.m3 <- inla(n_change ~ 1 + lognminus1+ span +span_2 + span_3 + rng + span_in + span_in_2 + span_in_3 +
                  f(transect,model="iid",hyper=rand_prior) + f(year,model="iid",hyper=rand_prior) +  # transect + year iid effects
                  f(species.0,model="iid",hyper=rand_prior) +
                  f(species,lognminus1,model="iid",hyper=rand_prior) +  # density per species
                  f(species.1,span,model="iid",hyper=rand_prior) + f(species.2,span_2,model="iid",hyper=rand_prior)+
                  f(species.4,span_in,model="iid") + f(species.5,span_in_2,model="iid"), # interaction with range per species
                  data=butterfly.dt,
                  control.inla =list(int.strategy ="eb",tolerance=1e-6,reordering ="metis",strategy="gaussian"),
                  control.compute = list(config=TRUE,waic = TRUE,dic=TRUE),
                  control.fixed = list(
                  mean=0, prec= 2, # fixed effect priors
                  mean.intercept= 0, prec.intercept= 1 )# priors on the intercept
)


#  base model 
summary(inla.m3)
bri.hyperpar.plot(inla.m3,together = F)

# check obs vs pred - fits pretty well but observed responses possible more extreme at edges SEE Appendix 3
library(scales)
plot(inla.m3$summary.linear.predictor$mean,butterfly.dt$n_change,col = alpha("grey", 0.4),xlab="Pred",ylab="obs")
abline(0,1)

# compare to linear model for a sanity check
library(lme4)
library(lmerTest)

# see how the priors are influencing relative to MLE
modelcompsp2<- lmer(n_change ~lognminus1+ span+rng +span_2+span_3  + span_in + span_in_2 + span_in_3 +
                      (poly(span,2, raw=TRUE)+span_in+span_in_2+lognminus1|species)+(1|transect)+(1|year),
                    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), na.action=na.fail,data = butterfly.dt, REML = T)


summary(modelcompsp2)



# --------- Fitting Spatial model -------------

# add in the spatial component to account for spatial auto correlation 
#get points 
xye<-as.matrix(cbind(butterfly.dt$xcor,butterfly.dt$ycor))

# draw a non-convex hull mesh for Europe
domain <- inla.nonconvex.hull(xye,
                              convex = -0.03, concave = -0.05,
                              resolution = c(40, 40))
meshe <- inla.mesh.2d(boundary = domain, max.edge = c(15000,100000),cutoff = 10000)
plot(meshe)
points(xye)

# make the projector matrix 
A <- inla.spde.make.A(meshe,loc=xye)
# prior on spatial effect
spde <- inla.spde2.pcmatern(mesh = meshe,
                            prior.range = c(50000, 0.05), #  less than prob
                            prior.sigma = c(5, 0.05))   # more than prob


# make index
index <- inla.spde.make.index("w", n.spde=spde$n.spde)

stack <- inla.stack(
  data = list(y = butterfly.dt$n_change), # response
  A = list(A,1),
  effects =  
    list( c(index, list(Intercept = 1)),
          list(span=butterfly.dt$span, 
               span_2 =butterfly.dt$span_2,
               span_3 =butterfly.dt$span_3,
               span_in=butterfly.dt$span_in, 
               span_in_2 =butterfly.dt$span_in_2,
               span_in_3 =butterfly.dt$span_in_3,
               rng = butterfly.dt$rng,
               species = butterfly.dt$species,
               species.0 = butterfly.dt$species,
               species.1 = butterfly.dt$species,
               species.2 = butterfly.dt$species,
               species.4 = butterfly.dt$species,
               species.5 = butterfly.dt$species,
               lognminus1 =butterfly.dt$lognminus1,
               transect = butterfly.dt$transect,
               year=butterfly.dt$year)),
  tag= "data"
)

# adding the spatial component doesn't seem to make a huge difference to coefficients but have much
# better fitting scores (wAIC)
inla.space1e2<-inla(y~0+Intercept + lognminus1+ span +span_2 + span_3 + rng + span_in + span_in_2 + span_in_3 +
                      f(w,model=spde)+ # spatial field 
                      f(transect,model="iid",hyper=rand_prior) + f(year,model="iid",hyper=rand_prior) +  # transect + year iid effects
                      f(species.0,model="iid",hyper=rand_prior) +
                      f(species,lognminus1,model="iid",hyper=rand_prior) +  # density per species
                      f(species.1,span,model="iid",hyper=rand_prior) + f(species.2,span_2,model="iid",hyper=rand_prior)+ # anomaly response per species
                      f(species.4,span_in,model="iid") + f(species.5,span_in_2,model="iid"),# interaction with range per species
                      data = inla.stack.data(stack),family="Gaussian",
                      control.inla =list(int.strategy ="eb",tolerance=1e-6,reordering ="metis",strategy="gaussian"),
                      control.compute = list(config=TRUE,waic = TRUE,dic=TRUE),
                      control.predictor = list(A = inla.stack.A(stack),compute=TRUE),
                      control.fixed = list(
                      mean=0, prec= 2, # fixed effect priors 
                      prec.intercept=1 ),# priors on the intercept, #  mean.intercept= 0, prec.intercept= 1
                    #  verbose = T
)


summary(inla.space1e2)

# plot spatial field
book.plot.field(
  inla.space1e2$summary.random$w$mean , 
  mesh = meshe) # the north east south west lower intercept is mainta

book.plot.field(
  inla.space1e2$summary.random$w$mean + inla.space1e2$summary.fixed$mean[1] , 
  mesh = meshe) # the north east south west lower intercept is mainta

# plot out effects 
points<-meshe$loc
step <- 2000 # m grid

east.range <- diff(range(points[,1]))  # calculate the length of the Easting range
north.range <- diff(range(points[,2])) # calculate the length of the Northing range

nxy <- round(c(east.range, north.range)/step)  # Calculate the number of cells in the x and y ranges

# Project the spatial field on the mesh vertices using the inla.mesh.projector() function
projgrid <- inla.mesh.projector(meshe,
                                xlim = range(points[,1]),
                                ylim = range(points[,2]),
                                dims = nxy)

xmean <- inla.mesh.project(projgrid,
                           inla.space1e2$summary.random$w$mean)

xsd <- inla.mesh.project(projgrid,
                         inla.space1$summary.random$w$sd)

library(raster)

xmean2 <- t(xmean)
xmean3 <- xmean2[rev(1:length(xmean2[,1])),]
xmean_ras <- raster(xmean3,
                    xmn = range(projgrid$x)[1], xmx = range(projgrid$x)[2],
                    ymn = range(projgrid$y)[1], ymx = range(projgrid$y)[2],
                    crs = 3035)

# this gives me a plot of the posterior random field - the response version below may be better 
# as it shows the size of the effect - seems like london and the north west have a lower growth rate 
plot(xmean_ras)
points(xye, pch = 16, cex = 0.2)


# ------------- Fitting boundary models -------------

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

countries <- ne_countries(scale = "medium", returnclass = "sf",country = c("United Kingdom","Spain","Finland"))

# fix geometry
c2<-st_transform(countries,crs=CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
plot(c2$geometry)
points(xye)

# convert into spatial polygons
peurp<-as_Spatial(c2$geometry)
cmesh <- inla.mesh.2d(boundary = peurp,
                      loc=xye,
                      max.edge = c(15000,200000),cutoff = 5000)

Ac <- inla.spde.make.A(cmesh,loc=xye)

# make the barrier 
# set the geometry
d<-st_as_sf(c2$geometry)
d2<-as_Spatial(d) # have to turn into sp
# points of the mesh over the polygon
land.tri <- inla.over_sp_mesh(d2, y = cmesh, 
                              type = "centroid",ignore.CRS = TRUE)
numtri<- length(cmesh$graph$tv[,1]) # number of triangles
barrier.tri <- setdiff(1:numtri,land.tri) # those not in the land are barriers 
poly.barrier <- inla.barrier.polygon(cmesh, 
                                     barrier.triangles = barrier.tri)
plot(poly.barrier,col="red") # take a look at barriers 
points(xye) # add points


# Make the barrier spde - also used in the model 
barrier.model <- inla.barrier.pcmatern(cmesh, 
                                       barrier.triangles = barrier.tri, 
                                       prior.range = c(50000, 0.05), #  less than prob
                                       prior.sigma = c(10, 0.05))

# I rewrite the stack here in a slightly different way as the other way was causing trouble
stackb <- inla.stack(
  data = list(y = butterfly.dt$n_change), # response
  A = list(Ac,1,1),
  effects=list(s= 1:cmesh$n,Intercept=rep(1,1504),
               list(span=butterfly.dt$span, 
                    span_2 =butterfly.dt$span_2,
                    span_3 =butterfly.dt$span_3,
                    span_in=butterfly.dt$span_in, 
                    span_in_2 =butterfly.dt$span_in_2,
                    span_in_3 =butterfly.dt$span_in_3,
                    rng = butterfly.dt$rng,
                    species = butterfly.dt$species,
                    species.0 = butterfly.dt$species,
                    species.1 = butterfly.dt$species,
                    species.2 = butterfly.dt$species,
                    species.4 = butterfly.dt$species,
                    species.5 = butterfly.dt$species,
                    lognminus1 =butterfly.dt$lognminus1,
                    transect = butterfly.dt$transect,
                    year=butterfly.dt$year)) # spatial random effect
  
  ,
  tag= "data"
)


inla.space1b<-inla(y~0+Intercept + lognminus1+ span +span_2 + span_3 + rng + span_in + span_in_2 + span_in_3 +
                     f(s,model=barrier.model)+
                     f(transect,model="iid",hyper=rand_prior) + f(year,model="iid",hyper=rand_prior) +  # transect + year iid effects
                     f(species.0,model="iid",hyper=rand_prior) +
                     f(species,lognminus1,model="iid",hyper=rand_prior) +  # density per species
                     f(species.1,span,model="iid",hyper=rand_prior) + f(species.2,span_2,model="iid",hyper=rand_prior)+ # anomaly response per species
                     f(species.4,span_in,model="iid") + f(species.5,span_in_2,model="iid"),
                     data = inla.stack.data(stackb),family="Gaussian",
                     control.inla =list(int.strategy ="eb",tolerance=1e-6,reordering ="metis",strategy="gaussian"),
                     control.compute = list(dic = TRUE,waic = TRUE),
                     control.fixed = list(
                     mean=0, prec= 2, # fixed effect priors 
                     mean.intercept= 0, prec.intercept= 1 ),# priors on the intercept
                     control.predictor = list(A = inla.stack.A(stackb),compute=TRUE)
)


summary(inla.space1b)
bri.hyperpar.plot(inla.space1b,together = F)

# the spatial field looks roughly consistent with previous
book.plot.field(
  inla.space1b$summary.random$s$mean + inla.space1b$summary.fixed$mean[1] , 
  mesh = cmesh) # the north east south west lower intercept is mainta

# just spatial field 
book.plot.field(
  inla.space1b$summary.random$s$mean  , 
  mesh = cmesh) 

# Pred vs obs - similar to the first model 
plot(inla.space1b$summary.fitted.values$mean[1:1504],butterfly.dt$n_change,col = alpha("grey", 0.4),xlab="Pred",ylab="obs")
abline(0,1)

# plots
book.plot.field(
  inla.space1b$summary.random$s$mean + inla.space1b$summary.fixed$mean[1] , 
  mesh = cmesh,poly = poly.barrier )

# test some nice colors other than viridis
library(RColorBrewer)
library(viridis)

p1<-plasma(201)

#tiff("sptf1.tif", res=600, compression = "lzw", height=6, width=10, units="in")
book.plot.field(
  inla.space1b$summary.random$s$mean + inla.space1b$summary.fixed$mean[1] , 
  mesh = cmesh,poly = poly.barrier,col=p1)
#dev.off()

#### other components 
# hyper parameter variation
bri.hyperpar.plot(inla.space1b,together = F)

# fixed effects
# density dependence
plot(inla.space1b$marginals.fixed$lognminus1,type="l",ylab="Density",
      xlab="lognminus1")
abline(v=inla.zmarginal(inla.space1b$marginals.fixed$lognminus1)[3],col="red",lty="dashed")
abline(v=inla.zmarginal(inla.space1b$marginals.fixed$lognminus1)[7],col="red",lty="dashed")

# anomaly 
plot(inla.space1b$marginals.fixed$span,type="l",ylab="Density",
     xlab="anomaly")
abline(v=inla.zmarginal(inla.space1b$marginals.fixed$span)[3],col="red",lty="dashed")
abline(v=inla.zmarginal(inla.space1b$marginals.fixed$span)[7],col="red",lty="dashed")
abline(v=0,col="black",lty="dashed")

# anomaly ^ 2
plot(inla.space1b$marginals.fixed$span_2,type="l",ylab="Density",
      xlab="anomaly2")
abline(v=inla.zmarginal(inla.space1b$marginals.fixed$span_2)[3],col="red",lty="dashed")
abline(v=inla.zmarginal(inla.space1b$marginals.fixed$span_2)[7],col="red",lty="dashed")
abline(v=0,col="black",lty="dashed")

# anomaly ^ 3
plot(inla.space1b$marginals.fixed$span_3,type="l",ylab="Density",
     xlab="anomaly3")
abline(v=inla.zmarginal(inla.space1b$marginals.fixed$span_3)[3],col="red",lty="dashed")
abline(v=inla.zmarginal(inla.space1b$marginals.fixed$span_3)[7],col="red",lty="dashed")
abline(v=0,col="black",lty="dashed")

# range position
plot(inla.space1b$marginals.fixed$rng,type="l",ylab="Density",
     xlab="Range position")
abline(v=inla.zmarginal(inla.space1b$marginals.fixed$rng)[3],col="red",lty="dashed")
abline(v=inla.zmarginal(inla.space1b$marginals.fixed$rng)[7],col="red",lty="dashed")
abline(v=0,col="black",lty="dashed")

# interaction first order
plot(inla.space1b$marginals.fixed$span_in,type="l",ylab="Density",
     xlab="anomaly*range positon")
abline(v=inla.zmarginal(inla.space1b$marginals.fixed$span_in)[3],col="red",lty="dashed")
abline(v=inla.zmarginal(inla.space1b$marginals.fixed$span_in)[7],col="red",lty="dashed")
abline(v=0,col="black",lty="dashed")

# interaction second order
plot(inla.space1b$marginals.fixed$span_in_2,type="l",ylab="Density",
     xlab="anomaly2*range position")
abline(v=inla.zmarginal(inla.space1b$marginals.fixed$span_in_2)[3],col="red",lty="dashed")
abline(v=inla.zmarginal(inla.space1b$marginals.fixed$span_in_2)[7],col="red",lty="dashed")
abline(v=0,col="black",lty="dashed")

# interaction third order
plot(inla.space1b$marginals.fixed$span_in_3,type="l",ylab="Density",
     xlab="anomaly3*range position")
abline(v=inla.zmarginal(inla.space1b$marginals.fixed$span_in_3)[3],col="red",lty="dashed")
abline(v=inla.zmarginal(inla.space1b$marginals.fixed$span_in_3)[7],col="red",lty="dashed")
abline(v=0,col="black",lty="dashed")


# classic check of residuals 
hist(inla.space1b$summary.fitted.values[id_space2, "mean"]-butterfly.dt$n_change)

# plot the main fit of the model
inla.space1b$summary.fixed # use means and sds

# make data 
simspan<-seq(from=range(butterfly.dt$span)[1],to=range(butterfly.dt$span)[2],0.1)
simrng<-seq(-1,1,0.2)

simdt<-expand.grid(simrng,simspan)
names(simdt)<-c("simrng","simspan")
simdt$simspan2 <-simdt$simspan^2
simdt$simspan3 <-simdt$simspan^3

simdt$simin1 <- simdt$simspan*simdt$simrng
simdt$simin2 <- simdt$simspan2*simdt$simrng
simdt$simin3 <- simdt$simspan3*simdt$simrng

# order as in dataframe
simdt<-simdt[,c(2,3,4,1,5,6,7)]

# Function to produce fit from model
findfit<-function(row,model){
  
  outm<-row[1]*model$summary.fixed$mean[3] +
    row[2]*model$summary.fixed$mean[4] +
    row[3]*model$summary.fixed$mean[5] +
    row[4]*model$summary.fixed$mean[6] +
    row[5]*model$summary.fixed$mean[7] +
    row[6]*model$summary.fixed$mean[8] +
    row[7]*model$summary.fixed$mean[9] 
  
  outlw<-row[1]* model$summary.fixed$`0.025quant`[3] +
    row[2]*model$summary.fixed$`0.025quant`[4]  +
    row[3]*model$summary.fixed$`0.025quant`[5]  +
    row[4]*model$summary.fixed$`0.025quant`[6]  +
    row[5]*model$summary.fixed$`0.025quant`[7]  +
    row[6]*model$summary.fixed$`0.025quant`[8]  +
    row[7]*model$summary.fixed$`0.025quant`[9] 
  
  outup<-row[1]*model$summary.fixed$`0.975quant`[3] +
    row[2]*model$summary.fixed$`0.975quant`[4]  +
    row[3]*model$summary.fixed$`0.975quant`[5]  +
    row[4]*model$summary.fixed$`0.975quant`[6]  +
    row[5]*model$summary.fixed$`0.975quant`[7]  +
    row[6]*model$summary.fixed$`0.975quant`[8]  +
    row[7]*model$summary.fixed$`0.975quant`[9] 
  
  
  ans<-c(outm,outlw,outup)
  return(ans)

    
}

# add to original data 
fits<-apply(simdt,MARGIN = 1,FUN = findfit,model=inla.space1b)
fits<-data.frame(t(fits))
names(fits) <- c("outm","outlw","outup")
simdt<-cbind(simdt,fits)

library(cowplot)

cols <- brewer.pal(11, "Spectral")
cols<-rev(cols)

# not the same fit due to low sample 
#tiff("fit1.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
ggplot(simdt,aes(x=simspan,y=outm,group=as.factor(simrng),color=as.factor(simrng)))+
  geom_line(size=1.5,alpha=0.9)+
  theme_cowplot()+
  scale_color_manual(values =cols)+
  theme(legend.position = "none",axis.text=element_text(size=20),
        axis.title=element_text(size=20))+
  ylab("Impact on growth rate")+
  xlab("Temperature anomaly")
#dev.off()


#tiff("fit1Un.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
ggplot(simdt,aes(x=simspan,y=outm,group=as.factor(simrng),color=as.factor(simrng)))+
  geom_line(size=1.5,alpha=0.9)+
  theme_cowplot()+
  geom_line(aes(x=simspan,y=outlw,group=as.factor(simrng),color=as.factor(simrng)),lty="dashed")+
  geom_line(aes(x=simspan,y=outup,group=as.factor(simrng),color=as.factor(simrng)),lty="dashed")+
  scale_color_manual(values =cols)+
  theme(legend.position = "none",axis.text=element_text(size=20),
        axis.title=element_text(size=20))+
  ylab("Impact on growth rate")+
  xlab("Temperature anomaly") +
  facet_wrap(~as.factor(simrng))
#dev.off()


# color bar
ranges<-seq(-1,1,0.2)
lev<-rep(1,11)
colorbar<-cbind.data.frame(ranges,lev)

# How joined up these are depends on the size
#tiff("colorbar.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
ggplot(colorbar,aes(ranges,lev,color=as.factor(ranges),group=as.factor(ranges)))+
  geom_point(size=20,shape=15)+
  scale_color_manual(values =cols)+
  theme_cowplot()+
  ylim(0.9,1)+
  theme(axis.title.y=element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),legend.position = "none",axis.text=element_text(size=14),
        axis.title=element_text(size=18)) +
 annotate("text",x=ranges,y=rep(0.99,11),label=c("-1","-0.8","-0.6","-0.4","-0.2","0","0.2","0.4","0.6","0.8","1"))
#dev.off()


# ------------- Reviewers request for Hierarchical GAM approach -------------
# Appendix 3
library(mgcv)

butterflym1$SPECIES_NAME= as.factor(butterflym1$SPECIES_NAME)
butterflym1$TRANSECT_ID= as.factor(butterflym1$TRANSECT_ID)

# this fails due to using sample data rather than the full dataset - but code is as used for Appendix 3
SPgam<-bam(n_change ~ lognminus1 + s(TRANSECT_ID, bs="re") + te(siteanomalyspring,rangepositionspring,bs=c("tp", "tp"),k=c(5, 5), m=2) +
                                       t2(siteanomalyspring,rangepositionspring, SPECIES_NAME, bs=c("tp", "tp", "re"), m=2, full=TRUE),
           data=butterflym1, method="REML")

summary(SPgam)
plot(SPgam)

library(ggeffects)
mydf3 <- ggpredict(SPgam, terms = c( "siteanomalyspring[all]","rangepositionspring[-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1]"))

cols2 <- brewer.pal(11, "Spectral")
cols2<-rev(cols2)

#tiff("gamalt.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
ggplot(mydf3, aes(x, predicted, colour = group)) +
  geom_line(size=1.5,alpha=0.9)+
  ylab("Impact on growth rate")+
  xlab("Temperature anomaly")+
  theme_cowplot()+
  scale_color_manual(values =cols2)+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20))+
   labs(color = "Thermal range \n position")
#dev.off()


# simplified model with spatial effects 
stackR <- inla.stack(
  data = list(y = butterfly.dt$n_change), # response
  A = list(Ac,1,1),
  effects=list(s= 1:cmesh$n,Intercept=rep(1,1504),
               list(span=butterfly.dt$span, 
                    span_in=butterfly.dt$span_in, 
                    rng = butterfly.dt$rng,
                    species = butterfly.dt$species,
                    species.0 = butterfly.dt$species,
                    species.1 = butterfly.dt$species,
                    lognminus1 =butterfly.dt$lognminus1,
                    transect = butterfly.dt$transect,
                    year=butterfly.dt$year)) # spatial random effect
  ,
  tag= "data"
)

inla.space1R<-inla(y~0+Intercept + lognminus1+ span + rng + span_in +
                     f(s,model=barrier.model)+
                     f(transect,model="iid",hyper=rand_prior) + f(year,model="iid",hyper=rand_prior) +  # transect + year iid effects
                     f(species.0,model="iid",hyper=rand_prior) +
                     f(species,lognminus1,model="iid",hyper=rand_prior) +  # density per species
                     f(species.1,span,model="iid") ,
                     data = inla.stack.data(stackR),family="Gaussian",
                     control.inla =list(int.strategy ="eb",tolerance=1e-6,reordering ="metis",strategy="gaussian"),
                     control.compute = list(dic = TRUE,waic = TRUE),
                     control.fixed = list(
                     mean=0, prec= 2, # fixed effect priors 
                     mean.intercept= 0, prec.intercept= 1 ),# priors on the intercept
                     control.predictor = list(A = inla.stack.A(stackR),compute=TRUE)
) 
summary(inla.space1R)

# function to produce model fit
findfit2<-function(row,model){
  
  outm<-row[1]*model$summary.fixed$mean[3] +
    row[2]*model$summary.fixed$mean[4] +
    row[3]*model$summary.fixed$mean[5] 

  outlw<-row[1]* (model$summary.fixed$mean[3] - (2*model$summary.fixed$sd[3]))
  row[2]*(model$summary.fixed$mean[4] - (2*model$summary.fixed$sd[4]))
  row[3]*(model$summary.fixed$mean[5] - (2*model$summary.fixed$sd[5]))

  outup<-row[1]* (model$summary.fixed$mean[3] + (2*model$summary.fixed$sd[3]))
  row[2]*(model$summary.fixed$mean[4] + (2*model$summary.fixed$sd[4]))
  row[3]*(model$summary.fixed$mean[5] + (2*model$summary.fixed$sd[5]))


  ans<-c(outm,outlw,outup)
  return(ans)
  
}

simdt2<-simdt[,c(1,4,5)]

fitsismp<-apply(simdt2,MARGIN = 1,FUN = findfit2,inla.space1R)
fitsismp<-data.frame(t(fitsismp))
names(fitsismp) <- c("outm","outlw","outup")
simdt2<-cbind(simdt,fitsismp)

cols <- brewer.pal(11, "Spectral")
cols<-rev(cols)


#tiff("simplemodel.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
ggplot(simdt2,aes(x=simspan,y=outm,group=as.factor(simrng),color=as.factor(simrng)))+
  geom_line(size=1.5,alpha=0.9)+
  theme_cowplot()+
  ylim(-0.9,0.3)+
  scale_color_manual(values =cols)+
  theme(legend.position = "none",axis.text=element_text(size=20),
        axis.title=element_text(size=20))+
  ylab("Impact on growth rate")+
  xlab("Temperature anomaly")
#dev.off()


# ------------- Reviewers request for Forest style plot -------------

rtab=cbind.data.frame(inla.space1b$names.fixed, inla.space1b$summary.fixed$mean,inla.space1b$summary.fixed$`0.025quant`,inla.space1b$summary.fixed$`0.975quant`,1:9)

names(rtab) = c("param","mean","lw","up","id")

ggplot(rtab, aes(mean,id))+
  geom_point(shape = 15,size=2)+
  theme_cowplot()+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_errorbarh(aes(xmax = up, xmin = lw,height=0))+
  xlab("Effect size")
  
  
