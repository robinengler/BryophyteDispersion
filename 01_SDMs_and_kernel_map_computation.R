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


### Fonction for kernel allowing to make an integration between the limits of a pixel ###

kernel <- function(m,a,b,d,f,g,h,q){
  # m the indices of a pixel value
  # a corresponds to sigmaw values
  # b corresponds to wind data
  # d = canopy height
  # f = spore release height
  # g = settling velocity
  # h and q corresponds to min and max distances between the source
  
  if(is.na(a[m])||is.na(b[m])||is.na(d[m])){return(NA)}
  else{
    return(integrate(f=function(x){return((f/(sqrt(0.3*d[m]*2*(a[m]/b[m]))*sqrt(2*pi*x^3)))*exp(-((x-(f*b[m]/g))^2)/(2*(0.3*d[m]*2*(a[m]/b[m]))*((b[m]^2)/g^2)*x)))},lower=h,upper=q)$value)
  } #integrate based on the wald algorithm (Bullock et al. 2012)
}

### End of the fonction ###

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
v1 <- subset (IndVar, 1)

biome<-shapefile("Variables/biomesEU2.shp")

# All the studied species
species <- c("Leje","Amphm","Plag","Amphl","Anas","Anom","Anth","Arct","Atri","Bart",
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
  removeTmpFiles(h=0.1)
}

## Merging the two technics by averaging the consensus maps
## Computing rasters between 2010 (present time) and 2050 

e <- extent(raster("Variables/wind/WindMax_2011-2020_MPI45.tif")) # Extent used for the study

for(i in 1: length(species)){
  
  # Species Data #
  spName <- species[i] # Species Name
  
  # Crop the different SDMs #
  presentGLM <- crop(raster(paste0(spName,"Pred_Pres_GLM.tif")),e)
  presentGBM <- crop(raster(paste0(spName,"Pred_Pres_GBM.tif")),e)
  MP45GLM <- crop(raster(paste0(spName,"Pred50MP45_GLM.tif")),e)
  MP45GBM <- crop(raster(paste0(spName,"Pred50MP45_GBM.tif")),e)
  MP85GLM <- crop(raster(paste0(spName,"Pred50MP85_GLM.tif")),e)
  MP85GBM <- crop(raster(paste0(spName,"Pred50MP85_GBM.tif")),e)

  # New Final map Formation
  pres <- overlay(presentGLM,presentGBM,fun=function(x,y){return((x+y)/2)},filename= paste0(spName,"_Present.tif"),overwrite=T)
  Mp45 <- overlay(MP45GLM,MP85GBM,fun=function(x,y){return((x+y)/2)},filename= paste0(spName,"_50MP45.tif"),overwrite=T)
  MP85 <- overlay(MP85GLM,MP85GBM,fun=function(x,y){return((x+y)/2)},filename= paste0(spName,"_50MP85.tif"),overwrite=T)
  
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

###############################
# Part 2 : kernel computation #
###############################

## Data loading and modification before the kernel computation ##
# Species Name & spore diameter data #

donnees <- read.table("tab_spore.csv",sep=";",dec=",",h=T) # A table containing the species names, the diameter spores (µm), and the kernel category allowing to regroup several species that have the same diameter spore (name of the column: Kernel)
catKern <- levels(donnees$Kernel)

# Release height (m) #
# Here, we tested several realease heights 
ZO <- c(0.03,1,10) 

# Mean and max wind speed values (m/s) during a time period of 10 years #
# scenario MPI rcp 8.5 #

# Modified Data from Euro-Cordex to allow 1km pixel resolution

# Mean wind speed
windMean10 <- raster("Variables/wind/WindMean201101-202012_MPI85.tif") #2011-2020
wMean10Val <- values(windMean10)
windMean20 <- raster("Variables/wind/WindMean202101-203012_MPI85.tif") #2021-2030
wMean20Val <- values(windMean20)
windMean30 <- raster("Variables/wind/WindMean203101-204012_MPI85.tif") #2031-2040
wMean30Val <- values(windMean30)
windMean40 <- raster("Variables/wind/WindMean204101-205012_MPI85.tif") #2041-2050
wMean40Val <- values(windMean40)

# Max wind speed
windMax10 <- raster("Variables/wind/WindMax_2011-2020_MPI85.tif") #2011-2020
wMax10Val <- values(windMax10)
windMax20 <- raster("Variables/wind/WindMax_2021-2030_MPI85.tif") #2021-2030
wMax20Val <- values(windMax20)
windMax30 <- raster("Variables/wind/WindMax_2031-2040_MPI85.tif") #2031-2040
wMax30Val <- values(windMax30)
windMax40 <- raster("Variables/wind/WindMax_2041-2050_MPI85.tif") #2041-2050
wMax40Val <- values(windMax40)

# Canopy height raster (m) #
# Data from the NASA 1km resolution
canop <- raster("Variables/canopy/europe_canop/w001001.adf")
values(canop)[values(canop)==0]=0.3 # Replace 0 values by 0.3 
canopVal <- values(canop)

## Computing sigmaw ##

# Step 1: determining friction velocity of our wind data with Z0=0.03 (canopy height = 30cm) #

fricSMean10 <- overlay(windMean10,fun=function(x){return((x*0.4)/log(10/0.03))})
fricSMean20 <- overlay(windMean20,fun=function(x){return((x*0.4)/log(10/0.03))})
fricSMean30 <- overlay(windMean30,fun=function(x){return((x*0.4)/log(10/0.03))})
fricSMean40 <- overlay(windMean40,fun=function(x){return((x*0.4)/log(10/0.03))})
fricSMax10 <- overlay(windMax10,fun=function(x){return((x*0.4)/log(10/0.03))})
fricSMax20 <- overlay(windMax20,fun=function(x){return((x*0.4)/log(10/0.03))})
fricSMax30 <- overlay(windMax30,fun=function(x){return((x*0.4)/log(10/0.03))})
fricSMax40 <- overlay(windMax40,fun=function(x){return((x*0.4)/log(10/0.03))})

# Step 2: Computing the friction velocity at the canopy level

# at 200m, u*_f = u*_s [log(200)-log(0.1*h_s)]/[log(200)-log(0.1*h_f)], where
# u*_f = friction coefficient in a forest area (canopy height different of 0.3), u*_s = friction coefficient for our wind data (canopy height of 0.3)
fricFMean10 <- overlay(fricSMean10,canop,fun=function(x,y){return(x*(log(200)-log(0.03))/(log(200)-log(0.1*y)))})
fricFMean20 <- overlay(fricSMean20,canop,fun=function(x,y){return(x*(log(200)-log(0.03))/(log(200)-log(0.1*y)))})
fricFMean30 <- overlay(fricSMean30,canop,fun=function(x,y){return(x*(log(200)-log(0.03))/(log(200)-log(0.1*y)))})
fricFMean40 <- overlay(fricSMean40,canop,fun=function(x,y){return(x*(log(200)-log(0.03))/(log(200)-log(0.1*y)))})
fricFMax10 <- overlay(fricSMax10,canop,fun=function(x,y){return(x*(log(200)-log(0.03))/(log(200)-log(0.1*y)))})
fricFMax20 <- overlay(fricSMax20,canop,fun=function(x,y){return(x*(log(200)-log(0.03))/(log(200)-log(0.1*y)))})
fricFMax30 <- overlay(fricSMax30,canop,fun=function(x,y){return(x*(log(200)-log(0.03))/(log(200)-log(0.1*y)))})
fricFMax40 <- overlay(fricSMax40,canop,fun=function(x,y){return(x*(log(200)-log(0.03))/(log(200)-log(0.1*y)))})


# Step 3: Depth-averaged vertical velocity standard deviation (m/s) #

sigmawMean10 <- 1.25*fricFMean10
siMean10Val <- values(sigmawMean10)
sigmawMean20 <- 1.25*fricFMean20
siMean20Val <- values(sigmawMean20)
sigmawMean30 <- 1.25*fricFMean30
siMean30Val <- values(sigmawMean30)
sigmawMean40 <- 1.25*fricFMean40
siMean40Val <- values(sigmawMean40)
sigmawMax10 <- 1.25*fricFMax10
siMax10Val <- values(sigmawMax10)
sigmawMax20 <- 1.25*fricFMax20
siMax20Val <- values(sigmawMax20)
sigmawMax30 <- 1.25*fricFMax30
siMax30Val <- values(sigmawMax30)
sigmawMax40 <- 1.25*fricFMax40
siMax40Val <- values(sigmawMax40)

# Save the maps #

writeRaster(sigmawMean10, paste0("Variables/sigmaw/sigmawMean10_MP85.tif"), datatype='FLT4S', overwrite=TRUE)
writeRaster(sigmawMean20, paste0("Variables/sigmaw/sigmawMean20_MP85.tif"), datatype='FLT4S', overwrite=TRUE)
writeRaster(sigmawMean30, paste0("Variables/sigmaw/sigmawMean30_MP85.tif"), datatype='FLT4S', overwrite=TRUE)
writeRaster(sigmawMean40, paste0("Variables/sigmaw/sigmawMean40_MP85.tif"), datatype='FLT4S', overwrite=TRUE)
writeRaster(sigmawMax10, paste0("Variables/sigmaw/sigmawMax10_MP85.tif"), datatype='FLT4S', overwrite=TRUE)
writeRaster(sigmawMax20, paste0("Variables/sigmaw/sigmawMax20_MP85.tif"), datatype='FLT4S', overwrite=TRUE)
writeRaster(sigmawMax30, paste0("Variables/sigmaw/sigmawMax30_MP85.tif"), datatype='FLT4S', overwrite=TRUE)
writeRaster(sigmawMax40, paste0("Variables/sigmaw/sigmawMax40_MP85.tif"), datatype='FLT4S', overwrite=TRUE)


siMean10Val <- values(sigmawMean10)
siMean20Val <- values(sigmawMean20)
siMean30Val <- values(sigmawMean30)
siMean40Val <- values(sigmawMean40)
siMax10Val <- values(sigmawMax10)
siMax20Val <- values(sigmawMax20)
siMax30Val <- values(sigmawMax30)
siMax40Val <- values(sigmawMax40)

# BaseMap #

# Here we used a wind speed raster and replace all the values by -9999
baseMap <- raster("Variables/wind/WindMean201101-202012_MPI85.tif")
values(baseMap)[!is.na(wMean10Val)]=-9999
e <- extent(baseMap) # keep the extent of the map

# Remove all the variables that is not useful for the computation #

rm(sigmawMax10,sigmawMax20,sigmawMax30,sigmawMax40,sigmawMean10,sigmawMean20,sigmawMean30,sigmawMean40,
   windMax10,windMax20,windMax30,windMax40,windMean10,windMean20,windMean30,windMean40,canop)

### Part 2 : Kernel computing ###

for(i in 1:length(catKern)){
  
  # Species names of those having the same spore diameter
  spNames <- as.character(donnees[donnees$Kernel==catKern[i],1]) 
  
  # Convert mean spore diameter in µm to m #
  
  spore <- unique(donnees[donnees$Kernel==catKern[i],2]) 
  d3 <- spore*0.000001 
  
  # Stokes? Law relation for determining the settling velocity of small spherical particles #
  
  # vset1 <- ((g*(d?)*(Pp-Pf))/(18*?)) 
  vset1 <- ((9.81*(d3*d3)*(1100-1.225))/(18*(0.0000178)))
  
  for(k in 1 : length(ZO)){
    
    # release height
    zo <- ZO[k]
    d1 = 500 # the min distance between the source pixel and the nearest pixels (1km resolution)
    d2 = 1500 # the max distance between the source pixel and the nearest pixels (1km resolution)
    
    stopVentMean <- 0 # Allows to stop the computation when the probability values = 0 everywhere 
    stopVentMax <- 0 # Allows to stop the computation when the probability values = 0 everywhere 
    
    for(l in 1: 10){ #We consider that after 10km, we do not have SDD but a LDD
      
      # Mean Wind #
      
      if(stopVentMean==0){
        
        # Computation of wald algorithm integration #
        
        kMean10 <- baseMap
        # Here we used a parallelisation but if you do not have enough power compution, 
        # you can use sapply instead of parSapply
        cl <-makeCluster(detectCores()-1) 
        # Use the fonction kernel on each pixel values of the raster
        valMe10 <- parSapply(cl,1:length(siMean10Val),kernel,a=siMean10Val,b=wMean10Val,d=canopVal,f=zo,g=vset1,h=d1,q=d2)
        stopCluster(cl)
        values(kMean10)=valMe10 # Add the computed values in the basemap
        maxValMe10 <- max(valMe10,na.rm=T) #Allows to see if the probability if > 0
        rm(valMe10)
        
        kMean20 <- baseMap
        cl <-makeCluster(detectCores()-1) 
        valMe20 <- parSapply(cl,1:length(siMean20Val),kernel,a=siMean20Val,b=wMean20Val,d=canopVal,f=zo,g=vset1,h=d1,q=d2)
        stopCluster(cl)
        values(kMean20)=valMe20
        maxValMe20 <- max(valMe20,na.rm=T)
        rm(valMe20)
        
        kMean30 <- baseMap
        cl <-makeCluster(detectCores()-1) 
        valMe30 <- parSapply(cl,1:length(siMean30Val),kernel,a=siMean30Val,b=wMean30Val,d=canopVal,f=zo,g=vset1,h=d1,q=d2)
        stopCluster(cl)
        values(kMean30)=valMe30
        maxValMe30 <- max(valMe30,na.rm=T)
        rm(valMe30)
        
        kMean40 <- baseMap
        cl <-makeCluster(detectCores()-1) 
        valMe40 <- parSapply(cl,1:length(siMean40Val),kernel,a=siMean40Val,b=wMean40Val,d=canopVal,f=zo,g=vset1,h=d1,q=d2)
        stopCluster(cl)
        values(kMean40)=valMe40
        maxValMe40 <- max(valMe40,na.rm=T)
        rm(valMe40)
        
        # Diffusion parameter into account #
        
        surf1 <- ((pi*(d2^2))-(pi*(d1^2)))/1000000 #give the number of pixels present at the same distance
        
        PDifMean10 <- overlay(kMean10, fun=function(a){return(a/surf1)})
        PDifMean20 <- overlay(kMean20, fun=function(a){return(a/surf1)})
        PDifMean30 <- overlay(kMean30, fun=function(a){return(a/surf1)})
        PDifMean40 <- overlay(kMean40, fun=function(a){return(a/surf1)})
        rm(kMean10,kMean20,kMean30,kMean40)
        
        # Successful Seeds calculation with coefficient (79.31) derived from Lonnel et al. 2012 with WALD reversely applied
        
        PDMean10 <- overlay(PDifMean10, fun=function(a){return(1-(1-a)^79.31)})
        PDMean20 <- overlay(PDifMean20, fun=function(a){return(1-(1-a)^79.31)})
        PDMean30 <- overlay(PDifMean30, fun=function(a){return(1-(1-a)^79.31)})
        PDMean40 <- overlay(PDifMean40, fun=function(a){return(1-(1-a)^79.31)})
        rm(PDifMean10,PDifMean20,PDifMean30,PDifMean40)
        
        # Make sure that the maps have the same extent with Na values corresponding to -9999
        finMean10 <- crop(PDMean10, e)
        crs(finMean10) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
        is.na(finMean10) <- -9999
        NAvalue(finMean10) <- -9999
        
        finMean20 <- crop(PDMean20, e)
        crs(finMean20) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
        is.na(finMean20) <- -9999
        NAvalue(finMean20) <- -9999
        
        finMean30 <- crop(PDMean30, e)
        crs(finMean30) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
        is.na(finMean30) <- -9999
        NAvalue(finMean30) <- -9999
        
        finMean40 <- crop(PDMean40, e)
        crs(finMean40) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
        is.na(finMean40) <- -9999
        NAvalue(finMean40) <- -9999
        
        rm(PDMean10,PDMean40,PDMean30,PDMean20)
      }
      
      # Max wind speed #
      
      if(stopVentMax==0){
        
        kMax10 <- baseMap
        cl <-makeCluster(detectCores()-1) 
        valMa10 <- parSapply(cl,1:length(siMax10Val),kernel,a=siMax10Val,b=wMax10Val,d=canopVal,f=zo,g=vset1,h=d1,q=d2)
        stopCluster(cl)
        values(kMax10)=valMa10
        maxValMa10 <- max(valMa10,na.rm=T)
        rm(valMa10)
        
        kMax20 <- baseMap
        cl <-makeCluster(detectCores()-1) 
        valMa20 <- parSapply(cl,1:length(siMax20Val),kernel,a=siMax20Val,b=wMax20Val,d=canopVal,f=zo,g=vset1,h=d1,q=d2)
        stopCluster(cl)
        values(kMax20)=valMa20
        maxValMa20 <- max(valMa20,na.rm=T)
        rm(valMa20)
        
        kMax30 <- baseMap
        cl <-makeCluster(detectCores()-1) 
        valMa30 <- parSapply(cl,1:length(siMax30Val),kernel,a=siMax30Val,b=wMax30Val,d=canopVal,f=zo,g=vset1,h=d1,q=d2)
        stopCluster(cl)
        values(kMax30)=valMa30
        maxValMa30 <- max(valMa30,na.rm=T)
        rm(valMa30)
        
        kMax40 <- baseMap
        cl <-makeCluster(detectCores()-1) 
        valMa40 <- parSapply(cl,1:length(siMax40Val),kernel,a=siMax40Val,b=wMax40Val,d=canopVal,f=zo,g=vset1,h=d1,q=d2)
        stopCluster(cl)
        values(kMax40)=valMa40
        maxValMa40 <- max(valMa40,na.rm=T)
        rm(valMa40)
        
        # Diffusion parameter into account #
        
        surf1 <- ((pi*(d2^2))-(pi*(d1^2)))/1000000 
        
        PDifMax10 <- overlay(kMax10, fun=function(a){return(a/surf1)})
        PDifMax20 <- overlay(kMax20, fun=function(a){return(a/surf1)})
        PDifMax30 <- overlay(kMax30, fun=function(a){return(a/surf1)})
        PDifMax40 <- overlay(kMax40, fun=function(a){return(a/surf1)})
        rm(kMax10,kMax20,kMax30,kMax40)
        
        # Successful Seeds calculation with coefficient (79.31) derived from Lonnel et al. 2012 with WALD reversely applied
        
        PDMax10 <- overlay(PDifMax10, fun=function(a){return(1-(1-a)^79.31)})
        PDMax20 <- overlay(PDifMax20, fun=function(a){return(1-(1-a)^79.31)})
        PDMax30 <- overlay(PDifMax30, fun=function(a){return(1-(1-a)^79.31)})
        PDMax40 <- overlay(PDifMax40, fun=function(a){return(1-(1-a)^79.31)})
        rm(PDifMax10,PDifMax20,PDifMax30,PDifMax40)
        
        finMax10 <- crop(PDMax10, e)
        crs(finMax10) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
        is.na(finMax10) <- -9999
        NAvalue(finMax10) <- -9999
        
        finMax20 <- crop(PDMax20, e)
        crs(finMax20) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
        is.na(finMax20) <- -9999
        NAvalue(finMax20) <- -9999
        
        finMax30 <- crop(PDMax30, e)
        crs(finMax30) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
        is.na(finMax30) <- -9999
        NAvalue(finMax10) <- -9999
        
        finMax40 <- crop(PDMax40, e)
        crs(finMax40) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
        is.na(finMax40) <- -9999
        NAvalue(finMax10) <- -9999
        rm(PDMax10,PDMax20,PDMax30,PDMax40)
      }
      
      # Check if the max probability value is > 0 #
      
      if(maxValMe10==0 & maxValMe20==0 & maxValMe30==0 & maxValMe40==0){stopVentMean=1}
      if(maxValMa10==0 & maxValMa20==0 & maxValMa30==0 & maxValMa40==0){stopVentMax=1}
      
      # Save the created maps for all species of a same category #
      
      for(j in 1 : length(spNames)){
        
        writeRaster(finMean10, paste0("Kernels/MP85/MeanW/",zo,"/", spNames[j],"_", d1 ,"-",d2, "_1.tif"), datatype='FLT4S', overwrite=TRUE)
        writeRaster(finMean20, paste0("Kernels/MP85/MeanW/",zo,"/", spNames[j],"_",d1 ,"-",d2, "_2.tif"), datatype='FLT4S', overwrite=TRUE)
        writeRaster(finMean30, paste0("Kernels/MP85/MeanW/",zo,"/", spNames[j],"_", d1 ,"-",d2,"_3.tif"), datatype='FLT4S', overwrite=TRUE)
        writeRaster(finMean40, paste0("Kernels/MP85/MeanW/",zo,"/", spNames[j],"_",d1 ,"-",d2, "_4.tif"), datatype='FLT4S', overwrite=TRUE)
        
        writeRaster(finMax10, paste0("Kernels/MP85/MaxW/",zo,"/", spNames[j],"_", d1 ,"-",d2,"_1.tif"), datatype='FLT4S', overwrite=TRUE)
        writeRaster(finMax20, paste0("Kernels/MP85/MaxW/",zo,"/", spNames[j],"_", d1 ,"-",d2,"_2.tif"), datatype='FLT4S', overwrite=TRUE)
        writeRaster(finMax30, paste0("Kernels/MP85/MaxW/",zo,"/", spNames[j],"_", d1 ,"-",d2,"_3.tif"), datatype='FLT4S', overwrite=TRUE)
        writeRaster(finMax40, paste0("Kernels/MP85/MaxW/",zo,"/", spNames[j],"_", d1 ,"-",d2,"_4.tif"), datatype='FLT4S', overwrite=TRUE)
        
      }
      
      removeTmpFiles(h=0.1) # Delete temporary files
      d1 =d1 +1000 # Go to the next distance pixels
      d2 =d2+1000
      rm(finMax10,finMean10,finMean20,finMean30,finMean40,finMax30,finMax20,finMax40)
    }
    
  }
}  
