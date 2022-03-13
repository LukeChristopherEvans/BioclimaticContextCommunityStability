# --------------- Code --------------
#  Bioclimatic context of species’ populations determines community stability'
#  H.3 Community stability from asynchrony
#  H.4 Community stability from average population stability 

# The code runs on a slice of the data as population data for a project must be requested from the 
# Finnish, Catalonian, BMS schemes. Consequently, results will not be identical to those presented in the manuscript

library(tidyverse)
library(INLA)
library(brinla)

# ----------- H4 -----------
# Note: this will not run due to the random slice, but code is kept in case it is useful

# ----------- load and prep data -----------

# load data 
stability<-read.delim("/BioclimaticContextCommunityStability/DataFiles/H2CommunityData.txt", sep= " ")
butterflym1<-read.delim("/BioclimaticContextCommunityStability/DataFiles/H1Slice.txt", sep= " ")


# get gamdens the density of the population 
butterflym1<-mutate(butterflym1, gamdens =  gam_index / transect_length)

# the inverse coefficient of variation
IVCOV<- function(datarange){
  cov<- ( sd(datarange,na.rm = T)  / mean(datarange,na.rm=T))
  cov <-  (1+ (1 / ( 4 * samplesize))) * cov
  incov1 <- 1 / cov
  return(incov1)
}


# work out the stability of individual species 
speciesstability<-butterflym1 %>% group_by(TRANSECT_ID,SPECIES_NAME) %>%
  summarise(n=n(),
            ICOV = IVCOV(gamdens,n),
            meanpop = mean(log(gamdens))
  )

speciesstability<-speciesstability[complete.cases(speciesstability),]

# get range 
IDCODE<-butterflym1 %>% distinct(TRANSECT_ID,SPECIES_NAME,.keep_all = T)

# add on range
speciesstability<-left_join(speciesstability,IDCODE,by=c("TRANSECT_ID","SPECIES_NAME"))

# get rid of too short timeseries 
speciesstability<-speciesstability %>% filter(n>=10)

speciesstability$rangepositionspring_2 <- speciesstability$rangepositionspring ^2


# -------- Fit non spatial model ----------

lmod<-lm(ICOV ~  rangepositionspring+ rangepositionspring_2 + meanpop ,data=speciesstability)
summary(lmod)
sdres <- sqrt(var(residuals(lmod)))

rand_prior <- list(prec = list(prior="pc.prec", param = c(3*sdres,0.01))) 

speciesstability$speciesname1 <- speciesstability$SPECIES_NAME
speciesstability$speciesname2 <- speciesstability$SPECIES_NAME
speciesstability$speciesname3 <- speciesstability$SPECIES_NAME
speciesstability$speciesname4 <- speciesstability$SPECIES_NAME

inla.m5<- inla(ICOV ~ 1 + rangepositionspring+ rangepositionspring_2 + meanpop +
                 f(TRANSECT_ID,model="iid",hyper=rand_prior) +  
                 f(speciesname1,model="iid",hyper=rand_prior) +
                 f(speciesname2,meanpop,model="iid",hyper=rand_prior) +
                 f(speciesname3,rangepositionspring,model="iid",hyper=rand_prior) + f(speciesname4,rangepositionspring_2,model="iid",hyper=rand_prior) 
                ,  
                 data=speciesstability,
                 control.compute = list(config=TRUE,waic = TRUE,dic=TRUE),
                 control.fixed = list(
                 mean=0, prec= 0.25, # fixed effect priors
                 mean.intercept= 0, prec.intercept= 1 )# p# priors on the intercept
)

summary(inla.m5)
bri.fixed.plot(inla.m5,together = F)
bri.hyperpar.plot(inla.m5,together = F)

# pretty normal
hist(speciesstability$ICOV-inla.m5$summary.fitted.values$mean)


# -------- Fit spatial barrier model ----------

## spatial field 
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

countries <- ne_countries(scale = "medium", returnclass = "sf",country = c("United Kingdom","Spain","Finland"))

# fix geometry
c2<-st_transform(countries,crs=CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
plot(c2$geometry)
points(speciesstability$transect_lon,speciesstability$transect_lat)

xye<-cbind(speciesstability$transect_lon,speciesstability$transect_lat)

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


stackd <- inla.stack(
  data = list(y = speciesstability$ICOV), # response
  A = list(Ac,1,1),
  effects=list(s= 1:cmesh$n,Intercept=rep(1,9053),
               list(rangepositionspring = speciesstability$rangepositionspring,
                    rangepositionspring_2 = speciesstability$rangepositionspring_2,
                    meanpop = speciesstability$meanpop,
                    TRANSECT_ID = speciesstability$TRANSECT_ID ,
                    YEAR= speciesstability$YEAR ,
                    speciesname4 = speciesstability$SPECIES_NAME,
                    speciesname1= speciesstability$speciesname1 ,
                    speciesname3= speciesstability$speciesname1 ,
                    speciesname2= speciesstability$speciesname2 )) # spatial random effect
  ,
  tag= "data"
)


inla.space1d<- inla(y~0+Intercept + rangepositionspring+ rangepositionspring_2+ meanpop+
                 f(s,model=barrier.model) +#+  f(YEAR,model="iid",hyper=rand_prior)+
                 f(TRANSECT_ID,model="iid",hyper=rand_prior) +
                 f(speciesname1,model="iid",hyper=rand_prior) +
                 f(speciesname2,meanpop,model="iid",hyper=rand_prior) +
                 f(speciesname3,rangepositionspring,model="iid",hyper=rand_prior) +f(speciesname4,rangepositionspring_2,model="iid",hyper=rand_prior)
                 , 
                 data = inla.stack.data(stackd),family="Gaussian",
               control.compute = list(config=TRUE,waic = TRUE,dic=TRUE),
               control.predictor = list(A = inla.stack.A(stackd),compute=TRUE),
               control.fixed = list(
                 mean=0, prec= 0.25, # fixed effect priors
                 mean.intercept= 0, prec.intercept= 1 )# p# priors on the intercept
)



summary(inla.space1d)


# -------- Plot spatial barrier model ----------


# code download https://www.r-inla.org/learnmore/books
source('INLAbookfunc/R/spde-book-functions.R')

bri.hyperpar.plot(inla.space1d,together = F)
bri.fixed.plot(inla.space1d,together = F)

book.plot.field(
  inla.space1d$summary.random$s$mean + inla.space1d$summary.fixed$mean[1] , 
  mesh = cmesh,poly = poly.barrier)

book.plot.field(
  inla.space1d$summary.random$s$mean, 
  mesh = cmesh,poly = poly.barrier)


# simulate the fit
simrange<-seq(from=range(speciesstability$rangepositionspring)[1],to=range(speciesstability$rangepositionspring)[2],0.01)
simrange_2 <- simrange^2

# estimate

estim<- inla.space1d$summary.fixed$mean[2] * simrange + inla.space1d$summary.fixed$mean[3] * simrange_2
estiml<- (inla.space1d$summary.fixed$mean[2] - (inla.space1d$summary.fixed$sd[2]*1.96)  ) * simrange + (inla.space1d$summary.fixed$mean[3]- (inla.space1d$summary.fixed$sd[2]*1.96)  ) * simrange_2
estimu<-(inla.space1d$summary.fixed$mean[2] + (inla.space1d$summary.fixed$sd[2]*1.96)  ) * simrange + (inla.space1d$summary.fixed$mean[3]+ (inla.space1d$summary.fixed$sd[2]*1.96)  ) * simrange_2

# make line dt
simrange<-cbind.data.frame(simrange,estim,estiml,estimu)
# make pt dt
ptdt<-cbind.data.frame(speciesstability$rangepositionspring,speciesstability$ICOV,speciesstability$meanpop, substr(speciesstability$TRANSECT_ID,1,1))
names(ptdt) <-c("rangescore","ICOV","meanpop","country")

ptdt<-mutate(ptdt, resid = ICOV -  (inla.space1d$summary.fixed$mean[1] + inla.space1d$summary.fixed$mean[4] * meanpop )  )

# ----------- Plot the fits  -----------

library(ggsci)
library(cowplot)

#tiff("fit3_1.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
ggplot(data=ptdt,aes(rangescore,resid,color=country,group=country))+
  geom_point(size=4,alpha=0.7)+
  scale_color_npg()+
  theme_cowplot() +
  geom_line(data=simrange,aes(simrange,estim), size=2,alpha=0.9,inherit.aes = F) +
  geom_line(data=simrange,aes(simrange,estiml), size=1,alpha=0.9,inherit.aes = F,lty="dashed") +
  geom_line(data=simrange,aes(simrange,estimu), size=1,alpha=0.9,inherit.aes = F,lty="dashed") +
  theme(legend.position = "none", axis.text=element_text(size=20),
        axis.title=element_text(size=20))+
  xlab("Range position")+
  ylab("Partial residual stability")+
  ylim(-1,1)
#dev.off()

# Fit of mean population

# residual dt
ptdt2<-cbind.data.frame(speciesstability$rangepositionspring,speciesstability$rangepositionspring_2,speciesstability$ICOV,speciesstability$meanpop, substr(speciesstability$TRANSECT_ID,1,1))
names(ptdt2) <- c("rng","rng2","ICOV","meanpop","country")

ptdt2<-mutate(ptdt2, resid = ICOV -  (inla.space1d$summary.fixed$mean[1] + inla.space1d$summary.fixed$mean[2] * rng+ inla.space1d$summary.fixed$mean[3] * rng2))             


# simulate the fit
simpop<-seq(from=range(speciesstability$meanpop)[1],to=range(speciesstability$meanpop)[2],0.02)


# estimate
estimP<-inla.space1d$summary.fixed$mean[4] * simpop
estimlP<-(inla.space1d$summary.fixed$mean[4] - inla.space1d$summary.fixed$sd[4]*1.96 ) * simpop
estimuP<-(inla.space1d$summary.fixed$mean[4] + inla.space1d$summary.fixed$sd[4]*1.96 ) * simpop 


simpop<-cbind.data.frame(simpop,estimP,estimlP,estimuP)


#tiff("fit3_2.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
ggplot(data=ptdt2,aes(meanpop,resid,color=country,group=country))+
  geom_point(size=4,alpha=0.7)+
  scale_color_npg()+
  theme_cowplot() +
  geom_line(data=simpop,aes(simpop,estimP), size=2,alpha=0.9,inherit.aes = F) +
  geom_line(data=simpop,aes(simpop,estimlP), size=1,alpha=0.9,inherit.aes = F,lty="dashed") +
  geom_line(data=simpop,aes(simpop,estimuP), size=1,alpha=0.9,inherit.aes = F,lty="dashed") +
  theme(legend.position = "none", axis.text=element_text(size=20),
        axis.title=element_text(size=20))+
  xlab("Log abundance")+
  ylab("Partial residual stability")
#dev.off()


# ----------- H3 -----------

### Function to calculate correlations 
findcor<-function(nams,dataf){
  
  name1 <- nams[1]
  name2<- nams[2]
  
  N1<-filter(dataf, SPECIES_NAME   == paste(name1))
  N2<- filter(dataf, SPECIES_NAME   == paste(name2))

  K1<-dplyr::select(N1, YEAR, n_change)
  K2<-dplyr::select(N2, YEAR, n_change)
  joined<-inner_join(K1,K2, by="YEAR") 
  joined<-joined[complete.cases(joined),]
  
  ANS<-ifelse(nrow(joined) >= 10,cor(joined$n_change.x,joined$n_change.y),NA) # use same limit
  
  return(ANS)
  
}

synchdat <- NULL

# run functions on data 
for(i in unique(butterflym1$TRANSECT_ID)){

  thisdata<-filter(butterflym1, TRANSECT_ID == i)
  
  specieshere<-unique(thisdata$SPECIES_NAME) 
  
  if(length(specieshere) > 1) {
  
  # get pair wise combinations
  pairwis<-combn(specieshere,2)
  pairwis<-t(pairwis)  
  combination<-cbind.data.frame(pairwis[,1],pairwis[,2])
  
  cors<-apply(combination,1,findcor,dataf=thisdata)
  
  namess<-combination[which(!is.na(cors)),]
  
  namess2<-unique(c(as.character(namess[,1]),as.character(namess[,2])))
  
  nspp <- length(namess2)
  avgcor<-mean(cors,na.rm=T)
  medcor<-median(cors,na.rm=T)
  sdcor<-sd(cors,na.rm=T)
  secor<- sdcor / sqrt(length(na.omit(cors)))
  varcor<- var(cors,na.rm=T)
  
  pos<-thisdata %>% 
    filter(gam_index > 0) %>% filter(!is.na(gam_index))%>% 
    distinct(SPECIES_NAME,.keep_all = T) %>%
    filter(SPECIES_NAME %in% namess2) %>%
    summarise(rangescore=mean(rangepositionspring,na.rm=T))
  
  if(nrow(pos) == 0){
    out<-c(NA,NA,NA,NA,NA,NA,NA)
  } else {
    
    out<-cbind.data.frame(i,pos,avgcor,medcor,nspp,sdcor,secor,varcor)
  }
  
  
  } else {out<-cbind.data.frame(i,pos,avgcor,medcor,nspp,sdcor,secor,varcor)}
  synchdat<-rbind(synchdat,out)
  
}  

# add on square
synchdat$rangescore_2 <- synchdat$rangescore^2
synchdat<-synchdat[complete.cases(synchdat),]


# -------- Fit model ----------
inla.m6<- inla(avgcor ~ 1 + rangescore+ rangescore_2 + nspp +varcor, 
               data=synchdat,
               control.compute = list(config=TRUE,waic = TRUE,dic=TRUE),
               control.fixed = list(
                 mean=0, prec= 0.25, # fixed effect priors
                 mean.intercept= 0, prec.intercept= 1 )# p# priors on the intercept
)

summary(inla.m6)
bri.fixed.plot(inla.m6)

inla.m6$summary.fixed


# use brms for weighted regression
library(brms)
synchdat$weights <- 1/synchdat$varcor

b6 <- brm(data = synchdat, family = gaussian,
      avgcor|weights(weights) ~ 1 +rangescore+ rangescore_2  + nspp ,
      prior = c(prior(normal(0, 1), class = Intercept),
                prior(normal(0, 2), class = b),
                prior(uniform(0, 10), class = sigma)),
      iter = 2000, warmup = 500, chains = 4, cores = 4,
      seed = 5)


plot(b6)
summary(b6)

bgsum<-summary(b6)

# get posterior 
sb6<-posterior_samples(b6)

# ----------- Plot the fits  -----------

# add country marker
synchdat$country  <- substr(synchdat$i,1,1)
library(ggsci)
library(cowplot)

ggplot(synchdat,aes(nspp,avgcor))+
  geom_point(aes(color=country,size=(weights)),alpha=0.7)+
  geom_abline(intercept = fixef(b6)[1,1],
              slope = fixef(b6)[4,1],size=2)+
  geom_abline(intercept = (fixef(b6)[1,3]),
              slope = fixef(b6)[4,3],lty="dashed")+
  geom_abline(intercept = (fixef(b6)[1,4]),
              slope = fixef(b6)[4,4],lty="dashed")+
  scale_color_npg()+
  theme_cowplot()+
  theme(legend.position = "none",axis.text=element_text(size=20),
        axis.title=element_text(size=20))+
  xlab("Species richness")+
  ylab("Mean synchrony")



# new version
simrangea<-seq(-1,1,0.1)
simrange_2a <- simrangea^2

# estimate
estima<- simrangea * mean(sb6$b_rangescore) +  simrange_2a * mean(sb6$b_rangescore_2)
estimal<- simrangea * (mean(sb6$b_rangescore)- 1.96*sd(sb6$b_rangescore)) + simrange_2a * (mean(sb6$b_rangescore_2)- 1.96*sd(sb6$b_rangescore_2))
estimau<- simrangea * (mean(sb6$b_rangescore)+ 1.96*sd(sb6$b_rangescore)) + simrange_2a * (mean(sb6$b_rangescore_2)+ 1.96*sd(sb6$b_rangescore_2))

# make line dt
simrangea<-cbind.data.frame(simrangea,estima,estimal,estimau)


# make pt dt
ptdta<-cbind.data.frame(synchdat$rangescore,synchdat$rangescore_2,synchdat$avgcor,synchdat$nspp ,substr(synchdat$i,1,1),synchdat$weights)
names(ptdta) <-c("rangescore","rangescore_2","avgcor","nspp","country","weights")

ptdta<-mutate(ptdta, resid = avgcor -  (mean(sb6$b_Intercept) + mean(sb6$b_nspp) *nspp )  )

#tiff("fit3_3.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
ggplot(data=ptdta,aes(rangescore,resid,color=country,group=country))+
  geom_point(aes(size=weights),alpha=0.7)+
  scale_color_npg()+
  theme_cowplot() +
  geom_line(data=simrangea,aes(simrangea,estima), size=2,alpha=0.9,inherit.aes = F) +
  geom_line(data=simrangea,aes(simrangea,estimal), size=1,alpha=0.9,inherit.aes = F,lty="dashed") +
  geom_line(data=simrangea,aes(simrangea,estimau), size=1,alpha=0.9,inherit.aes = F,lty="dashed") +
  theme(legend.position = "none", axis.text=element_text(size=20),
        axis.title=element_text(size=20))+
  xlab("Range position")+
  ylab("Partial residual synchrony")+
  ylim(-0.2,0.2)
#dev.off()



## now the species richness marginal
ptdta2<- ptdta 
ptdta2<-mutate(ptdta2, resid = avgcor -  (mean(sb6$b_Intercept) + mean(sb6$b_rangescore) *rangescore + mean(sb6$b_rangescore_2) *rangescore_2 )  )

# lines
simnspp<-seq(from=range(synchdat$nspp)[1],to=range(synchdat$nspp)[2],0.02)

# estimate
estimS<- mean(sb6$b_nspp) * simnspp
estimlS<-(mean(sb6$b_nspp) - 1.96 * sd(sb6$b_nspp))  * simnspp
estimuS<-(mean(sb6$b_nspp) + 1.96 * sd(sb6$b_nspp))  * simnspp

simnspp<-cbind.data.frame(simnspp,estimS,estimlS,estimuS)

#tiff("fit3_4.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
ggplot(data=ptdta2,aes(nspp,resid,color=country,group=country))+
  geom_point(aes(size=weights),alpha=0.7)+
  scale_color_npg()+
  theme_cowplot() +
  geom_line(data=simnspp,aes(simnspp,estimS), size=2,alpha=0.9,inherit.aes = F) +
  geom_line(data=simnspp,aes(simnspp,estimlS), size=1,alpha=0.9,inherit.aes = F,lty="dashed") +
  geom_line(data=simnspp,aes(simnspp, estimuS), size=1,alpha=0.9,inherit.aes = F,lty="dashed") +
  theme(legend.position = "none", axis.text=element_text(size=20),
        axis.title=element_text(size=20))+
  xlab("Species richness")+
  ylab("Partial residual synchrony")
#dev.off()



# ------------- Reviewers request for Forest style plot -------------

rtab=cbind.data.frame(names(bgsum$fixed[,1]),unname(bgsum$fixed[,1]),unname(bgsum$fixed[,3]),unname(bgsum$fixed[,4]),1:4)

names(rtab) = c("param","mean","lw","up","id")

# non-overlap
no= cbind.data.frame(c(0.2,0.2,0.2),c(1,3,4))
names(no) = c("x","y")             

#tiff("forestm3.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
ggplot(rtab, aes(mean,as.factor(id)))+
  geom_point(shape = 15,size=2)+
  theme_cowplot()+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_errorbarh(aes(xmax = up, xmin = lw,height=0))+
  xlab("Effect size")+
  ylab("")+
  scale_y_discrete(labels=c("1" = "Intercept", "2" = "Thermal\n range position",
                            "3" = "Thermal \n range position","4"="Species \n richness"))+
  geom_point(inherit.aes = F,data=no, aes(x,y),shape=8)
#dev.off()

