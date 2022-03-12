# --------------- Code --------------
#  Bioclimatic context of species’ populations determines community stability'
#  H.2 Community stability from thermal range position

library(tidyverse)
library(INLA)
library(brinla)


# ----------- load and prep data -----------

# load data 
stability<-read.delim("/BioclimaticContextCommunityStability/DataFiles/H2CommunityData.txt", sep= " ")


# -------- Fit non spatial model ----------

# default model prior to spde
# range of stability is 0.5 to 8 and range score smallest variable with size 2, so if it explained the full
# effect gradient would be 4, including this value within 2 sds gives prec 0.25
inla.m4 <- inla(ICOV ~ 1 + rangescore + rangescore_2 + meanpopsize + sppnum,
                data=stability,
                control.compute = list(config=TRUE,waic = TRUE,dic=TRUE),
                control.fixed = list(
                  mean=0, prec= 0.25, # fixed effect priors
                  mean.intercept= 0, prec.intercept= 1 )# p# priors on the intercept
)

summary(inla.m4)
bri.fixed.plot(inla.m4,together = F)

# -------- Fit spatial barrier model ----------

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

countries <- ne_countries(scale = "medium", returnclass = "sf",country = c("United Kingdom","Spain","Finland"))

# fix geometry
c2<-st_transform(countries,crs=CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
plot(c2$geometry)
points(stability$transect_lon,stability$transect_lat)

xye<-cbind(stability$transect_lon,stability$transect_lat)

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


# create data stack
stackc <- inla.stack(
  data = list(y = stability$ICOV), # response
  A = list(Ac,1,1),
  effects=list(s= 1:cmesh$n,Intercept=rep(1,616),
               list(meanpopsize = stability$meanpopsize,
                    rangescore = stability$rangescore ,
                    rangescore_2 = stability$rangescore_2,
                    sppnum = stability$sppnum)) # spatial random effect
  #iidx=1:nrow(df), # iid random effect
  ,
  tag= "data"
)

# run model
inla.space1c<-inla(y~0+Intercept + rangescore + rangescore_2 + meanpopsize + sppnum +
                   f(s,model=barrier.model),
                   data = inla.stack.data(stackc),family="Gaussian",
                   control.compute = list(config=TRUE,waic = TRUE,dic=TRUE),
                   control.predictor = list(A = inla.stack.A(stackc),compute=TRUE),    
                   control.fixed = list(
                     mean=0, prec= 0.25, # fixed effect priors
                     mean.intercept= 0, prec.intercept= 1 )# p# priors on the intercept
)


summary(inla.space1c)



# -------- Plot spatial barrier model ----------


# code download https://www.r-inla.org/learnmore/books
source('INLAbookfunc/R/spde-book-functions.R')

# plot hyper params
bri.hyperpar.plot(inla.space1c,together = F)

p1<-plasma(201)
#tiff("sptf2.tif", res=600, compression = "lzw", height=6, width=10, units="in")
book.plot.field(
  inla.space1c$summary.random$s$mean + inla.space1c$summary.fixed$mean[1] , 
  mesh = cmesh,poly = poly.barrier,col=p1)
#dev.off()

# classic check of residuals - slight skew but probably OK at this sample size
hist(stability$ICOV-inla.space1c$summary.fitted.values[id_space2, "mean"])

# plot out effects 
points<-cmesh$loc
step <- 2000 # m grid

east.range <- diff(range(points[,1]))  # calculate the length of the Easting range
north.range <- diff(range(points[,2])) # calculate the length of the Northing range

nxy <- round(c(east.range, north.range)/step)  # Calculate the number of cells in the x and y ranges

# Project the spatial field on the mesh vertices using the inla.mesh.projector() function
library(raster)
projgrid <- inla.mesh.projector(cmesh,
                                xlim = range(points[,1]),
                                ylim = range(points[,2]),
                                dims = nxy)

xmean <- inla.mesh.project(projgrid,
                           inla.space1c$summary.random$s$mean)


xmean2 <- t(xmean)
xmean3 <- xmean2[rev(1:length(xmean2[,1])),]
xmean_ras <- raster(xmean3,
                    xmn = range(projgrid$x)[1], xmx = range(projgrid$x)[2],
                    ymn = range(projgrid$y)[1], ymx = range(projgrid$y)[2],
                    crs = 3035)

# this gives me a plot of the posterior random field - the response version below may be better 
# as it shows the size of the effect - seems like london and the north west have a lower growth rate 
plot(xmean_ras)

# no, the values from this method make more sense 
epoints<-cbind(raster::extract(xmean_ras,SpatialPoints(xye[,1:2])),xye[,1:2])

max(epoints[,1]) +  2.298 # 3.5
min(epoints[,1]) +  2.298 # 1.75



# compare this to the fitted intercept using this method - same range 
book.plot.field(
      inla.space1c$summary.random$s$mean + inla.space1c$summary.fixed$mean[1] , 
      mesh = cmesh,poly = poly.barrier)
  
epoints[,1] <- epoints[,1] + inla.space1c$summary.fixed$mean[1]



# -------- Plot spatial model marginal fit ----------

# need an row id for epoints and stability 
stability$id <- 1:616
epoints<-data.frame(epoints)
epoints$id <- 1:616

# plot the marginal effects of rangescore - Now I want marginal after controlling for other effects 
simrange<-seq(from=range(stability$rangescore)[1],to=range(stability$rangescore)[2],0.01)

# function to get marginal fits
findmargresid <- function(row) {
  
  id <- as.numeric(row[10])

  marg <- epoints$X1[id] + as.numeric(row[4]) * inla.space1c$summary.fixed$mean[4] +  as.numeric(row[8]) * inla.space1c$summary.fixed$mean[5] 
  
  marg <- as.numeric(row[3]) - marg
  
  return(marg)
}

mfits<-apply(stability,MARGIN = 1,FUN = findmargresid)

# quick check
plot(stability$rangescore,mfits)


# now draw nicely 
simout<-simrange * inla.space1c$summary.fixed$mean[2] +  simrange ^ 2 * inla.space1c$summary.fixed$mean[3] 
simoutu<-simrange * (inla.space1c$summary.fixed$mean[2] + (inla.space1c$summary.fixed$sd[2]*1.96) )   +  simrange ^ 2 * (inla.space1c$summary.fixed$mean[3] +  (inla.space1c$summary.fixed$sd[3]*1.96)) 
simoutl<-simrange * (inla.space1c$summary.fixed$mean[2] - (inla.space1c$summary.fixed$sd[2]*1.96) )   +  simrange ^ 2 * (inla.space1c$summary.fixed$mean[3] -  (inla.space1c$summary.fixed$sd[3]*1.96)) 

# line
simout<-cbind.data.frame(simrange, simout,simoutu,simoutl)

# points
ptdt<-cbind.data.frame(stability$rangescore,mfits, substr(stability$TRANSECT_ID,1,1))
names(ptdt) <-c("rangescore","resids","country")

library(ggsci)
library(cowplot)

#tiff("fit2.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
ggplot(ptdt,aes(rangescore,resids,color=country,group=country))+
  geom_point(size=4,alpha=0.7)+
  scale_color_npg()+
  theme_cowplot() +
  geom_line(data=simout,aes(simrange,simout), size=2,alpha=0.9,inherit.aes = F) +
  geom_line(data=simout,aes(simrange,simoutu), size=1,alpha=0.9,inherit.aes = F,lty="dashed") +
  geom_line(data=simout,aes(simrange,simoutl), size=1,alpha=0.9,inherit.aes = F,lty="dashed") +
  theme(legend.position = "none", axis.text=element_text(size=20),
        axis.title=element_text(size=20))+
  xlab("Mean range score")+
  ylab("Partial residual community stability")
#dev.off()

             

# ------------- Reviewers request for Forest style plot -------------

rtab=cbind.data.frame(inla.space1c$names.fixed, inla.space1c$summary.fixed$mean,inla.space1c$summary.fixed$`0.025quant`,inla.space1c$summary.fixed$`0.975quant`,1:5)
summary(inla.space1c)
names(rtab) = c("param","mean","lw","up","id")

# non-overlap
no= cbind.data.frame(c(3,3,3,3),c(5,3,2,1))
names(no) = c("x","y")             

#tiff("forestm2.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
ggplot(rtab, aes(mean,as.factor(id)))+
  geom_point(shape = 15,size=2)+
  theme_cowplot()+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_errorbarh(aes(xmax = up, xmin = lw,height=0))+
  xlab("Effect size")+
  ylab("")+
  scale_y_discrete(labels=c("1" = "Intercept", "2" = "Thermal\n range position",
                          "3" = "Thermal \n range positon","4"="Mean\n abundance",
                           "5"="Species \n richness"))+
  geom_point(inherit.aes = F,data=no, aes(x,y),shape=8)

#dev.off()
