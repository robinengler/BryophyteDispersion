######################################
#                                    #
#     Zanatta, F. et al 2019 NC      #
#                                    #
######################################

######################################
#           Part 1: SDMs             #
#     Part 2: kernel computation     #
######################################

setwd() # add your pathway where you want to work

# Libraries
if(!require("raster")){
  install.packages("raster")
  library(raster)}
if(!require("foreign")){
  install.packages("foreign")
  library(foreign)}
if(!require("rgdal")){
  install.packages("rgdal")
  library(rgdal)}
if(!require("maptools")){
  install.packages("maptools")
  library(maptools)}
if(!require("dismo")){
  install.packages("dismo")
  library(dismo)}
if(!require("biomod2")){
  install.packages("biomod2")
  library(biomod2)}
if(!require("ecospat")){
  install.packages("ecospat")
  library(ecospat)}

library(parallel)

########################################
# Part 0: Function used in this script #
########################################

### Fonction for computing the friction coefficient for the canopy height at the wind measurement ###
## Present in Part 2

calcFricS <- function(wind){
  return((wind * 0.4) / log(10/0.03))
}

### Fonction for computing the friction coefficient with the canopy height of each pixels from friction coefficient at the canopy height from the measurement ###
## Present in Part 2

calcFricF <- function(fricS,hc){
  return(fricS*(log(200)-log(0.03))/(log(200)-log(0.1*hc)))
}

### Fonction for kernel allowing to make an integration between the limits of a pixel ###
## Present in Part 2

## Parameters
# d1 and d2 = the lower and the greater distance from the pixel source to the target pixel pixel
# sigmaw = sigmaw raster
# wind = wind speed raster
# hc = canopy height raster
# zo = spore release height
# vset = the settling velocity

kernel <- function(sigmaw,wind,hc){
  
  d1 <- get("d1")
  d2 <- get("d2")
  vset <- get("vset")
  zo <- get("zo")
  
  if(is.na(sigmaw) || is.na(wind) ||is.na(hc)){return(NA)}
  else{
    return(integrate(Vectorize(function(x){return((zo/(sqrt(0.4*hc*2*(sigmaw/wind))*sqrt(2*pi*x^3)))*exp(-((x-(zo*wind/vset))^2)/(2*(0.4*hc*2*(sigmaw/wind))*((wind^2)/vset^2)*x)))}),lower=d1,upper=d2)$value)
  } #integrate based on the wald algorithm (Bullock et al. 2012 & Katul et al. 2005)
}


### Change the wind at the different release heights based on two formula ###
## Present in Part 2

## Parameters
# z = release height
# hc = canopy height raster
# fricF = friction cofficient calculated with the canopy height of the pixel 

changeWind <- function(hc,fricF,z){
  if(z < hc ){ # if the release if below hc, we apply the formula from .... with the alpha computation present in ....
    alpha <- 0.24 + 0.096 * log(0.1*hc) + 0.016 * log(0.1*hc)^2
    windHc <- (fricF/0.4) * log((hc-0.7*hc)/(0.1*hc))
    return(windHc*exp(alpha*((z/hc)-1)))
  }
  else{ # if not, the wind is calculated with the log-profile wind speed
    return((fricF/0.4)*log((z-0.7*hc)/(0.1*hc)))
  }
}

### End of the function part ###

#####################
#   Part 1: SDMs    #
#####################

# Stack bioclimatic data from WorldClim v1.4
IndVar <- stack("Variables/_Present/bio_18.tif","Variables/_Present/bio_2.tif", "Variables/_Present/bio_4.tif", "Variables/_Present/bio_10.tif", "Variables/_Present/bio_13.tif")
proj4string(IndVar) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

IndVar2 <- stack("Variables/_mp45bi50/bio_18.tif","Variables/_mp45bi50/bio_2.tif", "Variables/_mp45bi50/bio_4.tif", "Variables/_mp45bi50/bio_10.tif", "Variables/_mp45bi50/bio_13.tif")
proj4string(IndVar2) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

IndVar3 <- stack("Variables/_mp85bi50/bio_18.tif","Variables/_mp85bi50/bio_2.tif", "Variables/_mp85bi50/bio_4.tif", "Variables/_mp85bi50/bio_10.tif", "Variables/_mp85bi50/bio_13.tif")
proj4string(IndVar3) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

IndVar4 <- stack("Variables/_he45bi50/bio_18.tif","Variables/_he45bi50/bio_2.tif", "Variables/_he45bi50/bio_4.tif", "Variables/_he45bi50/bio_10.tif", "Variables/_he45bi50/bio_13.tif")
proj4string(IndVar4) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

IndVar5 <- stack("Variables/_he85bi50/bio_18.tif","Variables/_he85bi50/bio_2.tif", "Variables/_he85bi50/bio_4.tif", "Variables/_he85bi50/bio_10.tif", "Variables/_he85bi50/bio_13.tif")
proj4string(IndVar5) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

v1 <- subset (IndVar, 1)

biome<-shapefile("Variables/biomesEU2.shp")

# All the studied species
species <- c("Amphm","Plag","Amphl","Anas","Anom","Anth","Arct","Atri","Bart",
             "Bazz","Cors","Cyno","Cyrt","Dicr","Dipla","Diplt","Fabr","Foss","Frul",
             "Glyph","Grim","Gymn","Habr","Herb","Homa","Lept","Leuc","Mast","Metz","Myur",
             "Ortha","Orthl","Oxym","Palu","Ptych","Sacc","Scle","Scor","Spha","Ulot")

for(i in 1:length(species)){
  
  # Store species name
  spName <- species[i]
  
  # Read species presence points, remove duplicates and points in the water (due to the resolution) 
  species.pres.temp <- read.dbf(paste0("Distrib/", spName,"/", spName, ".dbf"))
  species.pres.filt <- species.pres.temp[,3:4]
  species.pres.xy <- unique(species.pres.filt [!is.na (extract (v1, species.pres.filt)),])
  colnames(species.pres.xy)<-c("x","y") 
  
  # Remove species occurrences that are closer than 0.1°
  species.pres.xy<-ecospat.occ.desaggregation(species.pres.xy,min.dist=0.1)
  colnames(species.pres.xy)<-c("lon","lat")
  
  # Sample pseudo-absences (PA)
  background.tmp<-data.frame(coordinates(spsample(biome,50000,type="random")))
  background.xy<-background.tmp[!is.na(extract(v1, background.tmp)),]
  names(background.xy)<-c("lon","lat")
  
  # Species data with PA and true occurrences
  myResp.xy <- rbind(species.pres.xy,background.xy);names(myResp.xy)<-c("x","y");row.names(myResp.xy)<-c(1:nrow(myResp.xy))
  myResp <- data.frame(c(rep(1,nrow(species.pres.xy)),rep(NA,nrow(background.xy))));names(myResp)<-"pa"; row.names(myResp)<-c(1:nrow(myResp.xy))
  
  # Environmental data used to build the models
  myExpl <- data.frame (extract (IndVar, myResp.xy))
  
  # We used GLM and GBM to build the models. Each model techniques were run separately because we adapted the number of PA in function of the technics
  
  ## GLM part ##
  
  # Formating the data
  myBiomodData <- BIOMOD_FormatingData (resp.var = myResp, resp.xy = myResp.xy, expl.var = myExpl, resp.name = paste0(spName,"GLM"), PA.nb.rep = 1, PA.nb.absences = 10000, PA.strategy = 'random')
  
  # Default options
  myBiomodOption <- BIOMOD_ModelingOptions()
  
  # Modeling part
  myBiomodModelOut <- BIOMOD_Modeling (
    myBiomodData,
    models = 'GLM',
    models.options = myBiomodOption,
    NbRunEval=10,
    DataSplit=80,
    Yweights=NULL,
    VarImport=3,
    models.eval.meth = c('TSS','ROC'),
    SaveObj = TRUE,
    rescal.all.models = TRUE,
    do.full.models=T)
  
  #Ensemble Modeling
  myBiomodEM <- BIOMOD_EnsembleModeling(
    modeling.output = myBiomodModelOut,
    chosen.models = 'all',
    eval.metric = c('TSS'),
    eval.metric.quality.threshold = NULL,
    prob.mean = F,
    prob.cv = F,
    prob.ci = F,
    prob.ci.alpha = 0.05,
    prob.median = F,
    committee.averaging = F,
    prob.mean.weight = T,
    prob.mean.weight.decay = 'proportional' )
  
  # Evaluation outputs
  TSS<-get_evaluations(myBiomodModelOut)[1,,,,]
  ROC<-get_evaluations(myBiomodModelOut)[2,,,,]
  myModelsVarImport<-get_variables_importance(myBiomodModelOut)
  VarImport<-myModelsVarImport[,1,1,]
  evalAss <- get_evaluations(myBiomodEM)
  
  # Store the evaluation outputs
  write.table(TSS,paste0("evaluation/",spName,"_TSS_GLM.txt"),sep="\t",row.names=T)
  write.table(ROC,paste0("evaluation/",spName,"_ROC_GLM.txt"),sep="\t",row.names=T)
  write.table(evalAss,paste0("evaluation/",spName,"_eval_assemble_GLM.txt"),sep="\t",row.names=T)
  write.table(VarImport,paste0("evaluation/",spName,"_VarImport_GLM.txt"),sep="\t",row.names=T)
  
  # Projection of the GLM models
  myBiomodProj <- BIOMOD_Projection (modeling.output = myBiomodModelOut, new.env = IndVar, 
                                     proj.name = 'current',selected.models = 'all', 
                                     binary.meth = 'TSS')
  
  # Projection of the consensus GLM model
  myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                           EM.output = myBiomodEM)
  
  PredEF<-raster(paste0(spName, "GLM/proj_current/proj_current_", spName, "GLM_ensemble.grd"))
  
  # Store the raster at present time
  writeRaster (PredEF, paste0(spName,"Pred_pres_GLM.tif"), datatype='INT2S', options="COMPRESS=LZW", overwrite=TRUE)
  
  # MPI 4.5 scenario
  myBiomomodProj2050MP45 <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = IndVar2,
    proj.name = '2050MP45',
    selected.models = 'all',
    binary.meth = 'TSS',
    clamping.mask = T)
  
  myBiomodEF <- BIOMOD_EnsembleForecasting(
    projection.output = myBiomomodProj2050MP45,
    EM.output = myBiomodEM )
  
  MP452050<-raster(paste0(spName, "GLM/proj_2050MP45/proj_2050MP45_", spName, "GLM_ensemble.grd"))
  
  writeRaster (MP452050, paste0(spName,"Pred50MP45_GLM.tif"), datatype='INT2S', options="COMPRESS=LZW", overwrite=TRUE)
  
  # MPI 8.5 scenario
  myBiomomodProj2050MP85 <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = stack(IndVar3),
    proj.name = '2050MP85',
    selected.models = 'all',
    binary.meth = 'TSS',
    clamping.mask = T)
  
  myBiomodEF <- BIOMOD_EnsembleForecasting(
    projection.output = myBiomomodProj2050MP85,
    EM.output = myBiomodEM )
  
  MP852050<-raster(paste0(spName, "GLM/proj_2050MP85/proj_2050MP85_", spName, "GLM_ensemble.grd"))
  
  writeRaster (MP852050, paste0(spName,"Pred50MP85_GLM.tif"), datatype='INT2S', options="COMPRESS=LZW", overwrite=TRUE)
  
  ## GBM part ##
  
  myBiomodData <- BIOMOD_FormatingData (resp.var = myResp, resp.xy = myResp.xy, expl.var = myExpl, resp.name = paste0(spName,"GBM"), PA.nb.rep = 1, PA.nb.absences = nrow(species.pres.xy), PA.strategy = 'random')
  
  myBiomodOption <- BIOMOD_ModelingOptions(
    GBM = list(n.trees = 5000,
               cv.folds = 3,
               interaction.depth = 5,
               shrinkage = 0.001,
               bag.fraction = 0.8,
               train.fraction = 1.0)
  )
  
  myBiomodModelOut <- BIOMOD_Modeling (
    myBiomodData,
    models = 'GBM',
    models.options = myBiomodOption,
    NbRunEval=10,
    DataSplit=80,
    Yweights=NULL,
    VarImport=3,
    models.eval.meth = c('TSS','ROC'),
    SaveObj = TRUE,
    rescal.all.models = TRUE,
    do.full.models=T)
  
  #Ensemble Modeling
  myBiomodEM <- BIOMOD_EnsembleModeling(
    modeling.output = myBiomodModelOut,
    chosen.models = 'all',
    eval.metric = c('TSS'),
    eval.metric.quality.threshold = NULL,
    prob.mean = F,
    prob.cv = F,
    prob.ci = F,
    prob.ci.alpha = 0.05,
    prob.median = F,
    committee.averaging = F,
    prob.mean.weight = T,
    prob.mean.weight.decay = 'proportional' )
  
  TSS<-get_evaluations(myBiomodModelOut)[1,,,,]
  ROC<-get_evaluations(myBiomodModelOut)[2,,,,]
  myModelsVarImport<-get_variables_importance(myBiomodModelOut)
  VarImport<-myModelsVarImport[,1,1,]
  evalAss <- get_evaluations(myBiomodEM)
  write.table(TSS,paste0("evaluation/",spName,"_TSS_GBM.txt"),sep="\t",row.names=T)
  write.table(ROC,paste0("evaluation/",spName,"_ROC_GBM.txt"),sep="\t",row.names=T)
  write.table(evalAss,paste0("evaluation/",spName,"_eval_assemble_GBM.txt"),sep="\t",row.names=T)
  write.table(VarImport,paste0("evaluation/",spName,"_VarImport_GBM.txt"),sep="\t",row.names=T)
  
  
  myBiomodProj <- BIOMOD_Projection (modeling.output = myBiomodModelOut, new.env = IndVar, 
                                     proj.name = 'current',selected.models = 'all', 
                                     binary.meth = 'TSS')
  
  myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                           EM.output = myBiomodEM)
  
  PredEF<-raster(paste0(spName, "GBM/proj_current/proj_current_", spName, "GBM_ensemble.grd"))
  
  writeRaster (PredEF, paste0(spName,"Pred_pres_GBM.tif"), datatype='INT2S', options="COMPRESS=LZW", overwrite=TRUE)
  
  myBiomomodProj2050MP45 <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = IndVar2,
    proj.name = '2050MP45',
    selected.models = 'all',
    binary.meth = 'TSS',
    clamping.mask = T)
  
  myBiomodEF <- BIOMOD_EnsembleForecasting(
    projection.output = myBiomomodProj2050MP45,
    EM.output = myBiomodEM )
  
  MP452050<-raster(paste0(spName, "GBM/proj_2050MP45/proj_2050MP45_", spName, "GBM_ensemble.grd"))
  
  writeRaster (MP452050, paste0(spName,"Pred50MP45_GBM.tif"), datatype='INT2S', options="COMPRESS=LZW", overwrite=TRUE)
  
  myBiomomodProj2050MP85 <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = stack(IndVar3),
    proj.name = '2050MP85',
    selected.models = 'all',
    binary.meth = 'TSS',
    clamping.mask = T)
  
  myBiomodEF <- BIOMOD_EnsembleForecasting(
    projection.output = myBiomomodProj2050MP85,
    EM.output = myBiomodEM )
  
  MP852050<-raster(paste0(spName, "GBM/proj_2050MP85/proj_2050MP85_", spName, "GBM_ensemble.grd"))
  
  writeRaster (MP852050, paste0(spName,"Pred50MP85_GBM.tif"), datatype='INT2S', options="COMPRESS=LZW", overwrite=TRUE)
  
  myBiomomodProj2050HE45 <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = IndVar4,
    proj.name = '2050HE45',
    selected.models = 'all',
    binary.meth = 'TSS',
    clamping.mask = T)
  
  myBiomodEF <- BIOMOD_EnsembleForecasting(
    projection.output = myBiomomodProj2050HE45,
    EM.output = myBiomodEM )
  
  He452050<-raster(paste0(species[i], "GBM/proj_2050HE45/proj_2050HE45_", species[i], "GBM_ensemble.grd"))
  
  writeRaster (He452050, paste0(species[i],"Pred50HE45_GBM.tif"), datatype='INT2S', options="COMPRESS=LZW", overwrite=TRUE)
  
  
  myBiomomodProj2050HE85 <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = IndVar5,
    proj.name = '2050HE85',
    selected.models = 'all',
    binary.meth = 'TSS',
    clamping.mask = T)
  
  myBiomodEF <- BIOMOD_EnsembleForecasting(
    projection.output = myBiomomodProj2050HE85,
    EM.output = myBiomodEM )
  
  He852050<-raster(paste0(species[i], "GBM/proj_2050HE85/proj_2050HE85_", species[i], "GBM_ensemble.grd"))
  
  writeRaster (He852050, paste0(species[i],"Pred50HE85_GBM.tif"), datatype='INT2S', options="COMPRESS=LZW", overwrite=TRUE)
  removeTmpFiles(h=0.1)
}

## Merging the two technics by averaging the consensus maps
## Computing rasters between 2010 (present time) and 2050 
wind <-raster("Variables_kernels/WindMax_2011-2020_HE85.tif")
e <- extent(wind) # Extent used for the study
valWind <- values(wind) # to remove areas where we don't have wind speed data

for(i in 1: length(species)){
  
  # Species Data #
  spName <- species[i] # Species Name
  
  # Crop the different SDMs #
  presentGLM <- crop(raster(paste0(spName,"Pred_Pres_GLM.tif")),e)
  presentGBM <- crop(raster(paste0(spName,"Pred_Pres_GBM.tif")),e)
  
  HE45GLM <- crop(raster(paste0(spName,"GLM/",spName,"Pred50HE45_GLM.tif")),e)
  HE45GBM <- crop(raster(paste0(spName,"GBM/",spName,"Pred50HE45_GBM.tif")),e)
  HE85GLM <- crop(raster(paste0(spName,"GLM/",spName,"Pred50HE85_GLM.tif")),e)
  HE85GBM <- crop(raster(paste0(spName,"GBM/",spName,"Pred50HE85_GBM.tif")),e)

  MP45GLM <- crop(raster(paste0(spName,"Pred50MP45_GLM.tif")),e)
  MP45GBM <- crop(raster(paste0(spName,"Pred50MP45_GBM.tif")),e)
  MP85GLM <- crop(raster(paste0(spName,"Pred50MP85_GLM.tif")),e)
  MP85GBM <- crop(raster(paste0(spName,"Pred50MP85_GBM.tif")),e)
  
  # New Final map Formation
  pres <- overlay(presentGLM,presentGBM,fun=function(x,y){return((x+y)/2)})
  writeRaster(pres,paste0(spName,"_Present.tif"),overwrite=T)
  values(pres)[is.na(valWind)]=NA
  writeRaster(pres,paste0("pour_MigClim/",spName,"_Present.tif"),overwrite=T)
  HE45 <- overlay(HE45GLM,HE45GBM,fun=function(x,y){return((x+y)/2)})
  writeRaster(HE45,paste0(spName,"_50HE45.tif"),overwrite=T)
  values(HE45)[is.na(valWind)]=NA
  writeRaster(HE45,paste0("pour_MigClim/",spName,"_50HE45.tif"),overwrite=T)
  HE85 <- overlay(HE85GLM,HE85GBM,fun=function(x,y){return((x+y)/2)})
  writeRaster(HE85,paste0(spName,"_50HE85.tif"),overwrite=T)
  values(HE85)[is.na(valWind)]=NA
  writeRaster(HE85,paste0("pour_MigClim/",spName,"_50HE85.tif"),overwrite=T)
  
  MP45 <- overlay(MP45GLM,MP45GBM,fun=function(x,y){return((x+y)/2)})
  writeRaster(MP45,paste0(spName,"_50MP45.tif"),overwrite=T)
  values(MP45)[is.na(valWind)]=NA
  writeRaster(MP45,paste0("pour_MigClim/",spName,"_50MP45.tif"),overwrite=T)
  MP85 <- overlay(MP85GLM,MP85GBM,fun=function(x,y){return((x+y)/2)})
  writeRaster(MP85,paste0(spName,"_50MP85.tif"),overwrite=T)
  values(MP85)[is.na(valWind)]=NA
  writeRaster(MP85,paste0("pour_MigClim/",spName,"_50MP85.tif"),overwrite=T)
  
  # Computing maps between 2010 and 2050 (with the two scenarios)
  c2020 <- overlay(pres,HE45, fun=function(x,y){return(0.75*x+0.25*y)},filename=paste0("pour_MigClim/",spName,"_20HE45.tif"),overwrite=T)
  c2030 <- overlay(pres,HE45,fun=function(x,y){return(0.5*x+0.5*y)},filename=paste0("pour_MigClim/",spName,"_30HE45.tif"),overwrite=T)
  c2040 <- overlay(pres,HE45,fun=function(x,y){return(0.25*x+0.75*y)},filename=paste0("pour_MigClim/",spName,"_40HE45.tif"),overwrite=T)
  
  c2020 <- overlay(pres,HE85, fun=function(x,y){return(0.75*x+0.25*y)},filename=paste0("pour_MigClim/",spName,"_20HE85.tif"),overwrite=T)
  c2030 <- overlay(pres,HE85,fun=function(x,y){return(0.5*x+0.5*y)},filename=paste0("pour_MigClim/",spName,"_30HE85.tif"),overwrite=T)
  c2040 <- overlay(pres,HE85,fun=function(x,y){return(0.25*x+0.75*y)},filename=paste0("pour_MigClim/",spName,"_40HE85.tif"),overwrite=T)
  
  
  # Computing maps between 2010 and 2050 (with the two scenarios)
  c2020 <- overlay(pres,MP45, fun=function(x,y){return(0.75*x+0.25*y)},filename=paste0("pour_MigClim/",spName,"_20MP45.tif"),overwrite=T)
  c2030 <- overlay(pres,MP45,fun=function(x,y){return(0.5*x+0.5*y)},filename=paste0("pour_MigClim/",spName,"_30MP45.tif"),overwrite=T)
  c2040 <- overlay(pres,MP45,fun=function(x,y){return(0.25*x+0.75*y)},filename=paste0("pour_MigClim/",spName,"_40MP45.tif"),overwrite=T)

  c2020 <- overlay(pres,MP85, fun=function(x,y){return(0.75*x+0.25*y)},filename=paste0("pour_MigClim/",spName,"_20MP85.tif"),overwrite=T)
  c2030 <- overlay(pres,MP85,fun=function(x,y){return(0.5*x+0.5*y)},filename=paste0("pour_MigClim/",spName,"_30MP85.tif"),overwrite=T)
  c2040 <- overlay(pres,MP85,fun=function(x,y){return(0.25*x+0.75*y)},filename=paste0("pour_MigClim/",spName,"_40MP85.tif"),overwrite=T)

}


rm(list = ls())

########################################################

#######################################################
# Part 2 : kernel computation and data preparation    #
# Exemple with the scenario HadGem2-ES rcp 8.5 (HE85) #
#######################################################

## Data loading and modification before the kernel computation ##
# Species Name & spore diameter data #

donnees <- read.table("tab_spore.csv",sep=";",dec=",",h=T) # A data.frame containing the species names, the diameter spores (µm), and the kernel category allowing to regroup several species that have the same diameter spore (name of the column: Kernel)
catKern <- levels(donnees$Kernel)

# Release height (m) #
# Here, we tested several realease heights 
ZO <- c(0.03,1,10) 

# Mean and max wind speed values (m/s) during a time period of 10 years #

# Modified Data from Euro-Cordex to allow 1km pixel resolution

# Mean wind speed
windMean10 <- raster("Variables_Kernels/WindMean201101-202012_HE85.tif") #2011-2020
windMean20 <- raster("Variables_Kernels/WindMean202101-203012_HE85.tif") #2021-2030
windMean30 <- raster("Variables_Kernels/WindMean203101-204012_HE85.tif") #2031-2040
windMean40 <- raster("Variables_Kernels/WindMean204101-205012_HE85.tif") #2041-2050

# Max wind speed
windMax10 <- raster("Variables_Kernels/WindMax_2011-2020_HE85.tif") #2011-2020
windMax20 <- raster("Variables_Kernels/WindMax_2021-2030_HE85.tif") #2021-2030
windMax30 <- raster("Variables_Kernels/WindMax_2031-2040_HE85.tif") #2031-2040
windMax40 <- raster("Variables_Kernels/WindMax_2041-2050_HE85.tif") #2041-2050

# Canopy height raster (m) #
# Data from the NASA 1km resolution
canop <- raster("Variables_Kernels/europe_canop/w001001.adf")
values(canop)[values(canop)==0]=0.3 # Replace 0 values by 0.3 
canop <- crop(canop,extent(windMean10))

## Computing sigmaw ##

# Step 1: determining friction velocity of our wind data with Z0=0.03 (canopy height = 30cm) #

fricSMean10 <- calcFricS(windMean10)
fricSMean20 <- calcFricS(windMean20)
fricSMean30 <- calcFricS(windMean30)
fricSMean40 <- calcFricS(windMean40)
fricSMax10 <- calcFricS(windMax10)
fricSMax20 <- calcFricS(windMax20)
fricSMax30 <- calcFricS(windMax30)
fricSMax40 <- calcFricS(windMax40)

# Step 2: Computing the friction velocity at the canopy level

# at 200m, u*_f = u*_s [log(200)-log(0.1*h_s)]/[log(200)-log(0.1*h_f)], where
# u*_f = friction coefficient in a forest area (canopy height different of 0.3), u*_s = friction coefficient for our wind data (canopy height of 0.3)
fricFMean10 <- calcFricF(fricSMean10,canop)
fricFMean20 <- calcFricF(fricSMean20,canop)
fricFMean30 <- calcFricF(fricSMean30,canop)
fricFMean40 <- calcFricF(fricSMean40,canop)
fricFMax10 <- calcFricF(fricSMax10,canop)
fricFMax20 <- calcFricF(fricSMax20,canop)
fricFMax30 <- calcFricF(fricSMax30,canop)
fricFMax40 <- calcFricF(fricSMax40,canop)


# Step 3: Depth-averaged vertical velocity standard deviation (m/s) #

sigmawMean10 <- 1.25*fricFMean10
sigmawMean20 <- 1.25*fricFMean20
sigmawMean30 <- 1.25*fricFMean30
sigmawMean40 <- 1.25*fricFMean40
sigmawMax10 <- 1.25*fricFMax10
sigmawMax20 <- 1.25*fricFMax20
sigmawMax30 <- 1.25*fricFMax30
sigmawMax40 <- 1.25*fricFMax40

# Save the maps #

writeRaster(sigmawMean10, paste0("Variables_Kernels/sigmawMean10_HE85.tif"), datatype='FLT4S', overwrite=TRUE)
writeRaster(sigmawMean20, paste0("Variables_Kernels/sigmawMean20_HE85.tif"), datatype='FLT4S', overwrite=TRUE)
writeRaster(sigmawMean30, paste0("Variables_Kernels/sigmawMean30_HE85.tif"), datatype='FLT4S', overwrite=TRUE)
writeRaster(sigmawMean40, paste0("Variables_Kernels/sigmawMean40_HE85.tif"), datatype='FLT4S', overwrite=TRUE)
writeRaster(sigmawMax10, paste0("Variables_Kernels/sigmawMax10_HE85.tif"), datatype='FLT4S', overwrite=TRUE)
writeRaster(sigmawMax20, paste0("Variables_Kernels/sigmawMax20_HE85.tif"), datatype='FLT4S', overwrite=TRUE)
writeRaster(sigmawMax30, paste0("Variables_Kernels/sigmawMax30_HE85.tif"), datatype='FLT4S', overwrite=TRUE)
writeRaster(sigmawMax40, paste0("Variables_Kernels/sigmawMax40_HE85.tif"), datatype='FLT4S', overwrite=TRUE)


## wind speed computation at the different release height
z.003 <- windMean10
values(z.003) = 0.03# release height
windMean10.003 <- overlay(canop,fricFMean10,z.003,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMean201101-202012_HE85-h03.tif", overwrite=TRUE)
windMean20.003 <- overlay(canop,fricFMean20,z.003,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMean202101-203012_HE85-h03.tif", overwrite=TRUE)
windMean30.003 <- overlay(canop,fricFMean30,z.003,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMean203101-204012_HE85-h03.tif", overwrite=TRUE)
windMean40.003 <- overlay(canop,fricFMean40,z.003,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMean204101-205012_HE85-h03.tif", overwrite=TRUE)
windMax10.003 <- overlay(canop,fricFMax10,z.003,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMax_2011-2020_HE85-h03.tif", overwrite=TRUE)
windMax20.003 <- overlay(canop,fricFMax20,z.003,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMax_2021-2030_HE85-h03.tif", overwrite=TRUE)
windMax30.003 <- overlay(canop,fricFMax30,z.003,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMax_2031-2040_HE85-h03.tif", overwrite=TRUE)
windMax40.003 <- overlay(canop,fricFMax40,z.003,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMax_2041-2050_HE85-h03.tif", overwrite=TRUE)

z.1 <- windMean10
values(z.1) = 1 # release height
windMean10.1 <- overlay(canop,fricFMean10,z.1,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMean201101-202012_HE85-h1.tif", overwrite=TRUE)
windMean20.1 <- overlay(canop,fricFMean20,z.1,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMean202101-203012_HE85-h1.tif", overwrite=TRUE)
windMean30.1 <- overlay(canop,fricFMean30,z.1,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMean203101-204012_HE85-h1.tif", overwrite=TRUE)
windMean40.1 <- overlay(canop,fricFMean40,z.1,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMean204101-205012_HE85-h1.tif", overwrite=TRUE)
windMax10.1 <- overlay(canop,fricFMax10,z.1,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMax_2011-2020_HE85-h1.tif", overwrite=TRUE)
windMax20.1 <- overlay(canop,fricFMax20,z.1,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMax_2021-2030_HE85-h1.tif", overwrite=TRUE)
windMax30.1 <- overlay(canop,fricFMax30,z.1,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMax_2031-2040_HE85-h1.tif", overwrite=TRUE)
windMax40.1 <- overlay(canop,fricFMax40,z.1,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMax_2041-2050_HE85-h1.tif", overwrite=TRUE)

## 10 
z.10 <- windMean10
values(z.10) = 10# release height
windMean10.10 <- overlay(canop,fricFMean10,z.10,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMean201101-202012_HE85-h10.tif",overwrite=T)
windMean20.10 <- overlay(canop,fricFMean20,z.10,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMean202101-203012_HE85-h10.tif",overwrite=T)
windMean30.10 <- overlay(canop,fricFMean30,z.10,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMean203101-204012_HE85-h10.tif",overwrite=T)
windMean40.10 <- overlay(canop,fricFMean40,z.10,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMean204101-205012_HE85-h10.tif",overwrite=T)
windMax10.10 <- overlay(canop,fricFMax10,z.10,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMax_2011-2020_HE85-h10.tif", overwrite=TRUE)
windMax20.10 <- overlay(canop,fricFMax20,z.10,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMax_2021-2030_HE85-h10.tif", overwrite=TRUE)
windMax30.10 <- overlay(canop,fricFMax30,z.10,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMax_2031-2040_HE85-h10.tif", overwrite=TRUE)
windMax40.10 <- overlay(canop,fricFMax40,z.10,fun=Vectorize(changeWind),filename="Variables_Kernels/WindMax_2041-2050_HE85-h10.tif", overwrite=TRUE)

# Remove unused objects
rm(windMean10.003,windMean10.1,windMean10.10,windMean20.003,windMean20.1,windMean20.10,windMean30.003,windMean30.1,windMean30.10,
   windMean40.003,windMean40.1,windMean40.10,windMax10.003,windMax10.1,windMax10.10,windMax20.003,windMax20.1,windMax20.10,windMax30.003,
   windMax30.1,windMax30.10,windMax40.003,windMax40.1,windMax40.10,fricSMean10,fricSMean20,fricSMean30,fricSMean40,fricSMax10,fricSMax20,
   fricSMax30,fricSMax40,fricFMean10,fricFMean20,fricFMean30,fricFMean40,fricFMax10,fricFMax20,fricFMax30,fricFMax40,windMean10,windMean20,
   windMean30,windMean40,windMax10,windMax20,windMax30,windMax40)


### Kernel computing part ###

for(i in 1:length(catKern)){ # for each category of spore diameter
  
  # Species having the same spore diameter
  
  spore <- unique(donnees[donnees$Kernel==catKern[i],2])
  
  # Mean spore diameter (m) #
  
  d3 <- spore*0.000001
  
  # Stokes? Law relation for determining the settling velocity of small spherical particles #
  # vset1 <- ((g*(d?)*(Pp-Pf))/(18*?)) #
  
  vset <- ((9.81*(d3*d3)*(1100-1.225))/(18*(0.0000178)))
  
  for(k in 1 : length(ZO)){ # for each release height
    # zo
    zo <-  ZO[k]
    d1 <- 500
    d2 <- 1500
    
    # Mean wind speed raster (m/s) at the release height #
    
    windMean10 <- raster(paste0("Variables_Kernels/WindMean201101-202012_HE85-",nomHauteur[k],".tif"))
    
    windMean20 <- raster(paste0("Variables_Kernels/WindMean202101-203012_HE85-",nomHauteur[k],".tif"))
    
    windMean30 <- raster(paste0("Variables_Kernels/WindMean203101-204012_HE85-",nomHauteur[k],".tif"))
    
    windMean40 <- raster(paste0("Variables_Kernels/WindMean204101-205012_HE85-",nomHauteur[k],".tif"))
    
    windMax10 <- raster(paste0("Variables_Kernels/WindMax_2011-2020_HE85-",nomHauteur[k],".tif"))
    
    windMax20 <- raster(paste0("Variables_Kernels/WindMax_2021-2030_HE85-",nomHauteur[k],".tif"))
    
    windMax30 <- raster(paste0("Variables_Kernels/WindMax_2031-2040_HE85-",nomHauteur[k],".tif"))
    
    windMax40 <- raster(paste0("Variables_Kernels/WindMax_2041-2050_HE85-",nomHauteur[k],".tif"))
    
    for(l in 1: 10){ # between 1km to 10km 
      
      surf1 <- ((pi*(d2^2))-(pi*(d1^2)))/1000000 #donuts surface
      
      ## Wind Mean
      # Computation for the different time period between 2020 and 2050 
      kMean10 <- overlay(sigmawMean10,windMean10,canop,fun=Vectorize(kernel))
      # Division of the probability with the diffusion parameter
      PDifMean10 <- kMean10/surf1
      rm(kMean10)
      # Successful Seeds calculation with coefficient derived from L?nnel et al. 2012 with WALD reversely applied
      PDMean10 <- 1-(1-PDifMean10)^764.53
      rm(PDifMean10)
      is.na(PDMean10) <- -9999
      NAvalue(PDMean10) <- -9999
      writeRaster(PDMean10, paste0("Kernels/HE85/MeanW/",zo,"/",catKern[i],"_", d1 ,"-",d2, "_1.tif"), datatype='FLT4S', overwrite=TRUE)
      rm(PDMean10)
      
      kMean20 <- overlay(sigmawMean20,windMean20,canop,fun=Vectorize(kernel))
      PDifMean20 <- kMean20/surf1
      rm(kMean20)
      PDMean20 <- 1-(1-PDifMean20)^764.53
      is.na(PDMean20) <- -9999
      NAvalue(PDMean20) <- -9999
      rm(PDifMean20)
      writeRaster(PDMean20, paste0("Kernels/HE85/MeanW/",zo,"/",catKern[i],"_",d1 ,"-",d2, "_2.tif"), datatype='FLT4S', overwrite=TRUE)
      rm(PDMean20)
      
      kMean30 <- overlay(sigmawMean30,windMean30,canop,fun=Vectorize(kernel))
      PDifMean30 <- kMean30/surf1
      rm(kMean30)
      PDMean30 <- 1-(1-PDifMean30)^764.53
      rm(PDifMean30)
      is.na(PDMean30) <- -9999
      NAvalue(PDMean30) <- -9999
      writeRaster(PDMean30, paste0("Kernels/HE85/MeanW/",zo,"/",catKern[i],"_", d1 ,"-",d2,"_3.tif"), datatype='FLT4S', overwrite=TRUE)
      rm(PDMean30)
      
      
      kMean40 <- overlay(sigmawMean40,windMean40,canop,fun=Vectorize(kernel))
      PDifMean40 <- kMean40/surf1
      rm(kMean40)
      PDMean40 <- 1-(1-PDifMean40)^764.53
      rm(PDifMean40)
      is.na(PDMean40) <- -9999
      NAvalue(PDMean40) <- -9999
      writeRaster(PDMean40, paste0("Kernels/HE85/MeanW/",zo,"/",catKern[i],"_",d1 ,"-",d2, "_4.tif"), datatype='FLT4S', overwrite=TRUE)
      rm(PDMean40)
      
      
      ## Wind MAx
      kMax10 <- overlay(sigmawMax10,windMax10,canop,fun=Vectorize(kernel))
      PDifMax10 <- kMax10/surf1
      rm(kMax10)
      PDMax10 <- 1-(1-PDifMax10)^764.53
      rm(PDifMax10)
      is.na(PDMax10) <- -9999
      NAvalue(PDMax10) <- -9999
      writeRaster(PDMax10, paste0("Kernels/HE85/MaxW/",zo,"/", catKern[i],"_", d1 ,"-",d2,"_1.tif"), datatype='FLT4S', overwrite=TRUE)
      
      
      kMax20 <- overlay(sigmawMax20,windMax20,canop,fun=Vectorize(kernel))
      PDifMax20 <- kMax20/surf1
      rm(kMax20)
      PDMax20 <- 1-(1-PDifMax20)^764.53
      rm(PDifMax20)
      is.na(PDMax20) <- -9999
      NAvalue(PDMax20) <- -9999
      writeRaster(PDMax20, paste0("Kernels/HE85/MaxW/",zo,"/", catKern[i],"_", d1 ,"-",d2,"_2.tif"), datatype='FLT4S', overwrite=TRUE)
      rm(PDMax20)
      
      kMax30 <- overlay(sigmawMax30,windMax30,canop,fun=Vectorize(kernel))
      PDifMax30 <- kMax30/surf1
      rm(kMax30)
      PDMax30 <- 1-(1-PDifMax30)^764.53
      rm(PDifMax30)
      is.na(PDMax30) <- -9999
      NAvalue(PDMax30) <- -9999
      writeRaster(PDMax30, paste0("Kernels/HE85/MaxW/",zo,"/", catKern[i],"_", d1 ,"-",d2,"_3.tif"), datatype='FLT4S', overwrite=TRUE)
      rm(PDMax30)
      
      kMax40 <- overlay(sigmawMax40,windMax40,canop,fun=Vectorize(kernel))
      PDifMax40 <- kMax40/surf1
      rm(kMax40)
      PDMax40 <- 1-(1-PDifMax40)^764.53
      rm(PDifMax40)
      is.na(PDMax40) <- -9999
      NAvalue(PDMax40) <- -9999
      writeRaster(PDMax40, paste0("Kernels/HE85/MaxW/",zo,"/", catKern[i],"_", d1 ,"-",d2,"_4.tif"), datatype='FLT4S', overwrite=TRUE)
      rm(PDMax40)
      
      removeTmpFiles(h=0.1)
      d1 =d1 +1000 # increase the distance with the pixel size
      d2 =d2+1000
    }
    
  }
  
} 
