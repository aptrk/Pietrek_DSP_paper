#Author: Alejandro Pietrek
#Last updated 12/1/2024
#Using spOccupancy package to analyze diademed plover occupancy data


library(spOccupancy)
library(coda)
library(MCMCvis)
library(bayesplot)
library(stars)
library(tidyverse)
library(readxl)
library(patchwork)

#Load removal occupancy data for Phegornis and covariates

datar <- read_excel("datar.xlsx")
y1= select(datar,P1:P3) #Defining the Ys
zst <- apply(y1, 1, max,na.rm = TRUE) # Zst for initial values

#Corr matrix

covs <- datar[,4:18]

library(corrr)

tiff("correlation.jpeg", width = 7, height = 5, units = 'in', res = 300)
covs[,1:14] %>% correlate() %>% shave() %>% rplot(print_cor=TRUE)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()


#Add latent zst to the dataset

zst <- as.factor(zst)
datar <- cbind(datar,zst)
datar <- as_tibble(datar) #format for data input


####Visualize data####


#Pendiente (slope)

tiff("slope.jpeg", width = 7, height = 5, units = 'in', res = 300)
ggplot(datar, aes(OBJECTID,pendiente, colour =zst)) +    
  geom_point() +
  geom_text(aes(label = Vega), hjust = - 0.5)
dev.off()

#Area
tiff("area.jpeg", width = 7, height = 5, units = 'in', res = 300)
ggplot(datar, aes(OBJECTID, Shape_Area, label = Vega, colour =zst)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Vega), hjust = - 0.5)
dev.off()

#NDWI amp, NDWI min


p1 <- ggplot(datar, aes(OBJECTID, ndwi_amp, label = Vega,colour =zst)) +    
  geom_point() +
  geom_text(aes(label = Vega), hjust = - 0.5) #INTERESTING

p2 <- ggplot(datar, aes(OBJECTID, ndwiMean, label = Vega, colour =zst)) +    
  geom_point() +
  geom_text(aes(label = Vega), hjust = - 0.5) #INTERESTING

p3 <- ggplot(datar, aes(OBJECTID, ndwiMin, label = Vega, colour =zst)) +    
  geom_point() +
  geom_text(aes(label = Vega), hjust = - 0.5) #INTERESTING

p4 <- ggplot(datar, aes(OBJECTID, ndwiMax, label = Vega, colour =zst)) +    
  geom_point() +
  geom_text(aes(label = Vega), hjust = - 0.5) 

tiff("ndwi.jpeg", width = 12, height = 6, units = 'in', res = 300)
(p1+p2)/(p3+p4)
dev.off()

#TPI

ggplot(singleocc, aes(OBJECTID, tpi, label = Vega, colour =zst)) +    
  geom_point() +
  geom_text(aes(label = Vega), hjust = - 0.5)

#Grazing

tiff("grazing.jpeg", width = 7, height = 5, units = 'in', res = 300)
ggplot(datar, aes(OBJECTID, Grazing, label = Vega, colour =zst)) +    
  geom_point() +
  geom_text(aes(label = Vega), hjust = - 0.5)
dev.off()


#Plot grazing vs area

ggplot(datar, aes(Shape_Area, Grazing, label = Vega, colour =zst)) +    
  geom_point() +
  geom_text(aes(label = Vega), hjust = - 0.5)


####JUMP STRAIGHT TO THE MODELS, START WITH CONSTANT DETECTION####

#First define inits and priors, for details check the spOccupancy Vignette
#https://www.jeffdoser.com/files/spoccupancy-web/articles/modelfitting

#Data


datar <- as_tibble(datar)
str( sp.data <- list(y = y1) )


dsp.inits <- list(alpha = 0, 
                  beta = 0, 
                  z = apply(sp.data$y, 1, max, na.rm = TRUE))

dsp.priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
                   beta.normal = list(mean = 0, var = 2.72))

#MCMC parameters

n.burn <- 50000
n.thin <- 20
n.chains <- 3
n.samples <- 100000





#Constant detection

dsp.1 <- PGOcc(occ.formula = ~1 , 
                det.formula = ~1, 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors, 
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.1)
plot(dsp.1$beta.samples, density = FALSE)
plot(dsp.1$beta.samples, density = TRUE)


#Constant detection and occupancy by Area

str( sp.data <- list(y = y1, occ.covs = datar[,14]) )

dsp.2 <- PGOcc(occ.formula = ~ scale(Shape_Area) , 
                det.formula = ~1, 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors, 
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.2)
plot(dsp.2$beta.samples, density = FALSE)
plot(dsp.2$beta.samples, density = TRUE)



#Constant detection and occupancy by slope

str( sp.data <- list(y = y1, occ.covs = datar[,16]) )

dsp.3 <- PGOcc(occ.formula = ~scale(pendiente) , 
                det.formula = ~1, 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors, 
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.3)
plot(dsp.3$beta.samples, density = FALSE)
plot(dsp.3$beta.samples, density = TRUE)


#Constant detection and occupancy by grazing

str( sp.data <- list(y = y1, occ.covs = datar[,18]) )

dsp.4 <- PGOcc(occ.formula = ~scale(Grazing) , 
                det.formula = ~1, 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors, 
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.4)
plot(dsp.4$beta.samples, density = FALSE)
plot(dsp.4$beta.samples, density = TRUE)



#Constant detection and occupancy by NDWI

str( sp.data <- list(y = y1, occ.covs = datar[,5]) )

dsp.5 <- PGOcc(occ.formula = ~scale(ndwiMean) , 
                det.formula = ~1, 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors, 
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.5)
plot(dsp.5$beta.samples, density = FALSE)
plot(dsp.5$beta.samples, density = TRUE)

#K-fold validation

out.k.fold1 <- PGOcc(occ.formula = ~scale(ndwiMean) , 
                     det.formula = ~1, 
                     data = sp.data, 
                     inits = dsp.inits, 
                     n.samples = n.samples, 
                     priors = dsp.priors, 
                     n.omp.threads = 1, 
                     verbose = TRUE, 
                     n.report = 1000, 
                     n.burn = n.burn, 
                     n.thin = n.thin, 
                     n.chains = n.chains,
                     k.fold = 40) #This is a Leave one out validation

out.k.fold1$k.fold.deviance


#Constant detection and occupancy by area and NDWI

str( sp.data <- list(y = y1, occ.covs = datar[,c(5,14)]))

dsp.6 <- PGOcc(occ.formula = ~ scale(ndwiMean) + scale(Shape_Area) , 
                det.formula = ~1, 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors,
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.6)
plot(dsp.6$beta.samples, density = FALSE)
plot(dsp.6$beta.samples, density = TRUE)


#Constant detection and occupancy by area and slope

str( sp.data <- list(y = y1, occ.covs = datar[,c(14,16)]) )

dsp.7 <- PGOcc(occ.formula = ~ scale(pendiente) + scale(Shape_Area) , 
                det.formula = ~1, 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors,
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.7)
plot(dsp.7$beta.samples, density = FALSE)
plot(dsp.71$beta.samples, density = TRUE)




#Constant detection and occupancy by area and grazing

str( sp.data <- list(y = y1, occ.covs = datar[,c(14,18)]) )

dsp.8 <- PGOcc(occ.formula = ~ scale(Grazing) + scale(Shape_Area) , 
                det.formula = ~1, 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors,
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection


waicOcc(dsp.8)
plot(dsp.8$beta.samples, density = FALSE)
plot(dsp.8$beta.samples, density = TRUE)


#Constant detection and occupancy by ndwi and slope


str( sp.data <- list(y = y1, occ.covs = datar[,c(5,16)]))


dsp.9 <- PGOcc(occ.formula = ~ scale(ndwiMean) + scale(pendiente) , 
                det.formula = ~1, 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors,
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.9)
plot(dsp.9$beta.samples, density = FALSE)
plot(dsp.9$beta.samples, density = TRUE)



#Constant detection and occupancy by ndwi and grazing

str( sp.data <- list(y = y1, occ.covs = datar[,c(5,18)]) )


dsp.10 <- PGOcc(occ.formula = ~ scale(ndwiMean) + scale(Grazing) , 
                det.formula = ~1, 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors,
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.10)
plot(dsp.10$beta.samples, density = FALSE)
plot(dsp.10$beta.samples, density = TRUE)




#Constant detection and occupancy by slope and grazing

str( sp.data <- list(y = y1, occ.covs = datar[,c(16,18)] ) )


dsp.11 <- PGOcc(occ.formula = ~ scale(pendiente) + scale(Grazing), 
                det.formula = ~1, 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors,
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.11)
plot(dsp.11$beta.samples, density = FALSE)
plot(dsp.11$beta.samples, density = TRUE)


#Constant detection and occupancy with a an area-grazing interaction


str( sp.data <- list(y = y1, occ.covs = datar[,c(14,18)]) )


dsp.12 <- PGOcc(occ.formula = ~ scale(Shape_Area) * scale(Grazing), 
                det.formula = ~1, 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors,
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.12)
plot(dsp.12$beta.samples, density = FALSE)
plot(dsp.12$beta.samples, density = TRUE)



####INCORPORATING AREA TO MODEL DETECTION WITH REMOVAL DESIGN####



#occupancy by Area

str( sp.data <- list(y = y1, occ.covs = datar[,14] ,det.covs = datar[,14]) )

dsp.13 <- PGOcc(occ.formula = ~ scale(Shape_Area) , 
                det.formula = ~ scale(Shape_Area), 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors, 
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.13)
plot(dsp.13$beta.samples, density = FALSE)
plot(dsp.13$beta.samples, density = TRUE)



#occupancy by slope

str( sp.data <- list(y = y1, occ.covs = datar[,16] ,det.covs = datar[,14]) )

dsp.14 <- PGOcc(occ.formula = ~scale(pendiente) , 
                det.formula = ~scale(Shape_Area), 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors, 
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.14)
plot(dsp.14$beta.samples, density = FALSE)
plot(dsp.14$beta.samples, density = TRUE)



#occupancy by grazing

str( sp.data <- list(y = y1, occ.covs = datar[,18] , det.covs = datar[,14]) )

dsp.15 <- PGOcc(occ.formula = ~scale(Grazing) , 
                det.formula = ~scale(Shape_Area), 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors, 
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.15)
plot(dsp.15$beta.samples, density = FALSE)
plot(dsp.15$beta.samples, density = TRUE)




#occupancy by ndwi

str( sp.data <- list(y = y1, occ.covs = datar[,5] , det.covs = datar[,14]) )

dsp.16<- PGOcc(occ.formula = ~scale(ndwiMean) , 
               det.formula = ~scale(Shape_Area), 
               data = sp.data, 
               inits = dsp.inits, 
               n.samples = n.samples, 
               priors = dsp.priors, 
               n.omp.threads = 1, 
               verbose = TRUE, 
               n.report = 1000, 
               n.burn = n.burn, 
               n.thin = n.thin, 
               n.chains = n.chains)

#Model selection

waicOcc(dsp.16)
plot(dsp.16$beta.samples, density = FALSE)
plot(dsp.16$beta.samples, density = TRUE)
plot(dsp.16$alpha.samples, density = TRUE)


#occupancy by area and NDWI

str( sp.data <- list(y = y1, occ.covs = datar[,c(5,14)] ,det.covs = datar[,14]) )

dsp.17 <- PGOcc(occ.formula = ~ scale(ndwiMean) + scale(Shape_Area) , 
                det.formula = ~ scale(Shape_Area), 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors,
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.17)
plot(dsp.17$beta.samples, density = FALSE)
plot(dsp.17$beta.samples, density = TRUE)


#Constant detection and occupancy by area and slope

str( sp.data <- list(y = y1, occ.covs = datar[,c(14,16)] , det.covs = datar[,14]) )

dsp.18 <- PGOcc(occ.formula = ~ scale(pendiente) + scale(Shape_Area) , 
                det.formula = ~ scale(Shape_Area), 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors,
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.18)
plot(dsp.18$beta.samples, density = FALSE)
plot(dsp.18$beta.samples, density = TRUE)




#Detection  by area and occupancy by area and grazing

str( sp.data <- list(y = y1, occ.covs = datar[,c(14,18)] ,det.covs = datar[,14]) )

dsp.19<- PGOcc(occ.formula = ~ scale(Grazing) + scale(Shape_Area) , 
               det.formula = ~scale(Shape_Area), 
               data = sp.data, 
               inits = dsp.inits, 
               n.samples = n.samples, 
               priors = dsp.priors,
               n.omp.threads = 1, 
               verbose = TRUE, 
               n.report = 1000, 
               n.burn = n.burn, 
               n.thin = n.thin, 
               n.chains = n.chains)

#Model selection

waicOcc(dsp.19)
plot(dsp.19$beta.samples, density = FALSE)
plot(dsp.19$beta.samples, density = TRUE)



#occupancy by ndwi and slope

str( sp.data <- list(y = y1, occ.covs = datar[,c(5,16)] ,  det.covs = datar[,14]) )


dsp.20 <- PGOcc(occ.formula = ~ scale(ndwiMean) + scale(pendiente) , 
                det.formula = ~scale(Shape_Area), 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors,
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.20)
plot(dsp.20$beta.samples, density = FALSE)
plot(dsp.20$beta.samples, density = TRUE)
plot(dsp.20$alpha.samples, density = TRUE)

#Change parameter names for pretty plots

as.matrix(dsp.20$alpha.samples)
colnames(dsp.20$alpha.samples) <- c("Intercept", "Area")
as.matrix(dsp.20$beta.samples)
colnames(dsp.20$beta.samples) <- c("Intercept", "Mean NDWI", "Slope")

MCMCplot(dsp.20$beta.samples, ref_ovl = TRUE, ci = c(50, 95),sz_labels = 1.5)
MCMCplot(dsp.20$alpha.samples, ref_ovl = TRUE, ci = c(50, 95),sz_labels = 1.5)

#K-fold validation

out.k.fold2 <- PGOcc(occ.formula = ~ scale(ndwiMean) + scale(pendiente) , 
                    det.formula = ~scale(Shape_Area), 
                    data = sp.data, 
                    inits = dsp.inits, 
                    n.samples = n.samples, 
                    priors = dsp.priors,
                    n.omp.threads = 1, 
                    verbose = TRUE, 
                    n.report = 1000, 
                    n.burn = n.burn, 
                    n.thin = n.thin, 
                    n.chains = n.chains,
                    k.fold = 40)


out.k.fold2$k.fold.deviance # 81.34

#Goodness of fit, posterior predictive plots

ppc.sp.out1 <- ppcOcc(dsp.44, fit.stat = 'freeman-tukey', group = 1) #across sites
summary(ppc.sp.out1)

ppc.sp.out2 <- ppcOcc(dsp.44, fit.stat = 'freeman-tukey', group = 2) #across surveys
summary(ppc.sp.out2)


#occupancy by ndwi and grazing


str( sp.data <- list(y = y1, occ.covs = datar[,c(5,18)] , det.covs = datar[,14]))


dsp.21 <- PGOcc(occ.formula = ~ scale(ndwiMean) + scale(Grazing) , 
                det.formula = ~ scale(Shape_Area), 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors,
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.21)
plot(dsp.21$beta.samples,density = FALSE) 
plot(dsp.21$beta.samples, density = TRUE)


#occupancy by slope and grazing


str( sp.data <- list(y = y1, occ.covs = datar[,c(16,18)] ,  det.covs = datar[,14],
                     coords = coord) )


dsp.22 <- PGOcc(occ.formula = ~ scale(pendiente) + scale(Grazing), 
                det.formula = ~scale(Shape_Area), 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors,
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.22)
plot(dsp.22$beta.samples, density = FALSE)
plot(dsp.22$beta.samples, density = TRUE)


#occupancy with a an area-grazing interaction



str( sp.data <- list(y = y1, occ.covs = datar[,c(14,18)] ,det.covs = datar[,14]) )


dsp.23 <- PGOcc(occ.formula = ~ scale(Shape_Area) * scale(Grazing), 
                det.formula = ~scale(Shape_Area), 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors,
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.23)
plot(dsp.23$beta.samples, density = FALSE)
plot(dsp.23$beta.samples, density = TRUE)


#occupancy constant and detection by area


str( sp.data <- list(y = y1 ,  det.covs = datar[,14]) )


dsp.24 <- PGOcc(occ.formula = ~ 1, 
                det.formula = ~scale(Shape_Area), 
                data = sp.data, 
                inits = dsp.inits, 
                n.samples = n.samples, 
                priors = dsp.priors,
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = n.burn, 
                n.thin = n.thin, 
                n.chains = n.chains)

#Model selection

waicOcc(dsp.24)
plot(dsp.24$beta.samples, density = FALSE)
plot(dsp.24$beta.samples, density = TRUE)


####Add spatial autocorrelation to the the model with constant detection and NDWI####

#Adding space to data


coord <- datar[,19:20]

str( sp.data <- list(y = y1, occ.covs = datar[,5] ,
                     coords = coord) )



#Define priors and parameters

# Pair-wise distances between all sites
dist.hbef <- dist(sp.data$coords)
# Exponential covariance model
cov.model <- "exponential"
# Specify list of inits
sp.inits <- list(alpha = 0, 
                 beta = 0, 
                 z = apply(sp.data$y, 1, max, na.rm = TRUE), 
                 sigma.sq = 1.5, 
                 phi = 3 / mean(dist.hbef), 
                 w = rep(0, nrow(sp.data$y)))


#MCMC parameters

batch.length <- 25
n.batch <- 4000
n.burn <- 50000
n.thin <- 20
n.chains <- 3

sp.tuning <- list(phi = 1)
# accept.rate = 0.43 by default, so we do not specify it.


#Specifying priors

min.dist <- min(dist.hbef)
max.dist <- max(dist.hbef)
sp.priors <- list(beta.normal = list(mean = 0, var = 2.72), 
                  alpha.normal = list(mean = 0, var = 2.72), 
                  sigma.sq.ig = c(2, 1), 
                  phi.unif = c(3/max.dist, 3/min.dist))

#Computation parameters

n.omp.threads <- 1
verbose <- TRUE
n.report <- 100 # Report progress at every 100th batch.

#Fitting the model

dsp.sp <- spPGOcc(occ.formula = ~scale(ndwiMean), 
                  det.formula = ~1, 
                  data = sp.data, 
                  inits = sp.inits, 
                  n.batch = n.batch, 
                  batch.length = batch.length, 
                  priors = sp.priors, 
                  cov.model = cov.model, 
                  NNGP = TRUE, 
                  n.neighbors = 10,
                  tuning = sp.tuning, 
                  n.report = n.report, 
                  n.burn = n.burn, 
                  n.thin = n.thin, 
                  n.chains = n.chains)

#Model selection

waicOcc(dsp.sp)
plot(dsp.sp$beta.samples, density = TRUE)
plot(dsp.sp$theta.samples,density = TRUE)


#Change parameter names for pretty plots

as.matrix(dsp.sp$alpha.samples)
colnames(dsp.sp$alpha.samples) <- c("Intercept")
as.matrix(dsp.sp$beta.samples)
colnames(dsp.sp$beta.samples) <- c("Intercept", "Mean NDWI")

MCMCplot(dsp.sp$beta.samples, ref_ovl = TRUE, ci = c(50, 95),sz_labels = 1.5)
MCMCplot(dsp.sp$alpha.samples, ref_ovl = TRUE, ci = c(50, 95),sz_labels = 1.5)

#K-fold validation

out.k.fold3 <-  spPGOcc(occ.formula = ~scale(ndwiMean), 
                        det.formula = ~1, 
                        data = sp.data, 
                        inits = sp.inits, 
                        n.batch = n.batch, 
                        batch.length = batch.length, 
                        priors = sp.priors, 
                        cov.model = cov.model, 
                        NNGP = TRUE, 
                        n.neighbors = 10,
                        tuning = sp.tuning, 
                        n.report = n.report, 
                        n.burn = n.burn, 
                        n.thin = n.thin, 
                        n.chains = n.chains,
                        k.fold = 40)



out.k.fold3$k.fold.deviance # 81.34

#Goodness of fit, posterior predictive plots

ppc.sp.out1 <- ppcOcc(dsp.sp, fit.stat = 'freeman-tukey', group = 1) #across sites
summary(ppc.sp.out1)

ppc.sp.out2 <- ppcOcc(dsp.sp, fit.stat = 'freeman-tukey', group = 2) #across surveys
summary(ppc.sp.out2)


#Plot means of random effects


means <-colMeans(dsp.sp$w.samples)

tiff("plot5.tiff", width = 7, height = 5, units = 'in', res = 300)
plot(means, xlab="Bog", ylab= "Spatial random effect", col= 2,pch = 19)+
  abline(h=0)
dev.off()




#Posterior predictive checks

ppc.sp.out1 <- ppcOcc(dsp.sp, fit.stat = 'freeman-tukey', group = 1) #across sites
summary(ppc.sp.out1)

ppc.sp.out2 <- ppcOcc(dsp.sp, fit.stat = 'freeman-tukey', group = 2) #across surveys
summary(ppc.sp.out2)

#Posterior predictive plot

ppc.df <- data.frame(fit = ppc.sp.out1$fit.y, 
                     fit.rep = ppc.sp.out1$fit.y.rep, 
                     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True')
lines(ppc.df$fit, ppc.df$fit, col = 'black')


tiff("plot6.tiff", width = 7, height = 5, units = 'in', res = 300)
diff.fit <- ppc.sp.out1$fit.y.rep.group.quants[3, ] - ppc.sp.out1$fit.y.group.quants[3, ]
plot(diff.fit, pch = 19, xlab = 'Site ID', ylab = 'Replicate - True Discrepancy') #Vegas 16 y 27 desviadas
axis(1, at = seq(1, 40, by = 1), las=2)
dev.off()







# One of the paper figures!

p5 <-   ggplot(datar, aes(OBJECTID, ndwiMean, colour =zst)) +    
  geom_point(size = 2) +
  theme_bw()+
  scale_color_discrete(name= "Latent occurence") +
  theme(
    # Hide panel borders and remove grid lines
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change axis line
    axis.line = element_line(colour = "black")
    ) +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank() 
    ) +
    labs(y = "Mean NDWI") 
   

# Same thing with min NDWI

p6 <-   ggplot(datar, aes(OBJECTID, ndwiMin, colour =zst)) +    
  geom_point(size = 2) +
  theme_bw()+
  scale_color_discrete(name= "Latent occurence") +
  theme(
    # Hide panel borders and remove grid lines
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change axis line
    axis.line = element_line(colour = "black")
  ) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() 
  ) +
  labs(y = "Min NDWI") +
  theme(legend.position="none")


# Try boxplots


datar <- datar %>%
    mutate(Detection = recode_factor(zst, `1` = "Detected", `0` = "Not detected" ))
  
  

p7 <-   ggplot(datar, aes(Detection, ndwiMean, colour = Detection)) + 
  geom_boxplot(outlier.alpha = 0, width = 0.5) +
  geom_jitter()+
  theme_bw()+
  theme(
    # Hide panel borders and remove grid lines
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change axis line
    axis.line = element_line(colour = "black")
  ) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position="none"
  ) +
  labs(y = "NDWI mean") +
  labs(x = c("Detected", "Not detected"))




p8 <-   ggplot(datar, aes(Detection, ndwiMin, colour = Detection)) + 
  geom_boxplot(outlier.alpha = 0, width = 0.5) +
  geom_jitter()+
  theme_bw()+ 
  theme(
    # Hide panel borders and remove grid lines
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change axis line
    axis.line = element_line(colour = "black")
  ) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position="none"
  ) +
  labs(y = "NDWI min") +
  labs(x = c("Detected", "Not detected")) 


   

#Matrix format of p. samples and rename parameters

posterior <- as.matrix(dsp.sp$beta.samples)
colnames(posterior) <- c("Intercept", "NDWI mean")



#Caterpillar plots parameters  

p9 <- mcmc_intervals(posterior, pars = c("Intercept", "NDWI mean")) + 
  theme(text = element_text(family = "Arial", size = 12))




tiff("fig3.tiff", width = 7.5, height = 4.5, units = 'in', res = 300)
(p8/p7)|p9
dev.off()


# We'll try to map random effects



# Load necessary packages (if not already loaded)
library(leaflet)
library(sf)
library(mapview)



# Sample data with UTM coordinates and associated values

utm_x <- datar$coords.x1
utm_y <- datar$coords.x2
utm_zone <- 19

# Create an sf data frame with UTM coordinates and CRS (Zone 19S)
utm_data <- data.frame(x = utm_x, y = utm_y)
utm_crs <- st_crs(paste0("+proj=utm +zone=", utm_zone, " +south +datum=WGS84"))
utm_sf <- st_as_sf(utm_data, coords = c("x", "y"), crs = utm_crs)

# Convert UTM coordinates to geographic coordinates (decimal degrees)
geo_sf <- st_transform(utm_sf, 4326)  # EPSG 4326 is WGS84

# Extract latitude and longitude values
latitude <- st_coordinates(geo_sf)[, 2]
longitude <- st_coordinates(geo_sf)[, 1]


# Convert UTM coordinates to latitude and longitude
data_points <- cbind(latitude,longitude, means)
data_points <- as.data.frame(data_points)

# Create a data frame with your data points
# (Same data frame creation code as before)

# Create a leaflet map with satellite background and inset map


library(leaflet)



# Set the desired extent (bounding box)
lng1 <- -66.5  # Minimum longitude
lat1 <- -22.485    # Minimum latitude
lng2 <- -66.9   # Maximum longitude
lat2 <- -22.70   # Maximum latitude



p10 <-  leaflet(data = data_points) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(
    lng = ~longitude,
    lat = ~latitude,
    radius = ~sqrt(abs(means)) * 10,
    color = ~ifelse(means >= 0, "red", "blue"),
    stroke = TRUE,
    fillOpacity = 0.6,
    popup = ~paste("Latitude:", latitude, "<br>Longitude:", longitude, "<br>Mean:", means)
  ) %>%
  setView(lng = (lng1 + lng2) / 2, lat = (lat1 + lat2) / 2, zoom = 13) %>% 
   addMiniMap(
     tiles = providers$Jawg.Sunny,
     position = 'bottomright', 
     width = 315, height = 240,
     toggleDisplay = FALSE
   ) %>% 
  addScaleBar(position = "topright",options = scaleBarOptions( maxWidth = 300,
                                                                 metric = TRUE,
                                                                 imperial = FALSE)) %>% 
  addGraticule(interval = .2) %>% 
 fitBounds(lng1, lat1, lng2, lat2) 



# Capture and save the map as an image
 
 mapshot( p10,  # Capture the last plotted map
  file = "leaflet_map.png",  # Specify the output file name (use .png for image format)
  remove_controls = c("zoomControl", "layersControl", "homeButton",
                      "drawToolbar", "easyButton"))


 ######THE END##################################################################
 
 
 
 