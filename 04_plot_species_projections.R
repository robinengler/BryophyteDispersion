#!/software/R/3.5.1/bin/R
########################################################################################################################
### Plot projections of suitable habitat for year 2010 and 2050.
### ***********************************************************
### This script creates plots (for each of the 40 species) of a species distribution under current (2010) and 
### future (2050) climatic conditions. It does so for two different reclassification thresholds: "tss" and "min"
###  -> "tss" = maximum TSS threshold, the threshold that maximises the TSS (True Skill Statistic) value.
###  -> "min" = threshold corresponding to the minimum projected value observed at the locations of species occurence.
###
### The idea is that we want to compare the projections when reclassified by "tss" and by "min" so that we can 
### possibly evaluate which one is the more "realistic" to use.
### 
########################################################################################################################

########################################################################################################################
### Data input and environment setup.
### ********************************
#
# Input/Output locations.
rootDir    = "/home/rengler/Documents/Projects/Ecospat/2018_MigClim_FZanatta/inputData"
hsmapDir   = file.path(rootDir, "hsMaps")
spFile     = file.path(rootDir, "speciesSummaryTable.txt")
ccScenario = 'mp45'

# Load species file.
spTable = read.table(spFile, h=TRUE, sep='\t', as.is=TRUE)
spList = spTable$speciesAbrev[which( !is.na(spTable$maxTSSthreshold) & !is.na(spTable$minProjValue)) ]

# Setup environment.
setwd(rootDir)
require(rgdal)
require(raster)
###
########################################################################################################################


########################################################################################################################
### Plot projected distributions
### ******************************************
### 
#
pdf("tss_vs_min.pdf", bg="white")

for(spName in spList){
    
    # Progress report.
    cat("### Plotting data for species", spName, ":\n")
    
    proj2010 = raster(paste0(hsmapDir,'/',spName,'/',spName,'_GLMGBM_2010.tif'))
    proj2050 = raster(paste0(hsmapDir,'/',spName,'/',spName,'_GLMGBM_2050_',ccScenario,'.tif'))
    NAvalue(proj2010) = -9999
    NAvalue(proj2050) = -9999
    
    # Reclass raster values with "tss" or "min" thresholds.
    tssThreshold = spTable[which(spTable$speciesAbrev==spName),'maxTSSthreshold'] 
    minThreshold = spTable[which(spTable$speciesAbrev==spName),'minProjValue']
    if(minThreshold > tssThreshold) stop("ERROR: tss threshold should not be smaller than min threshold.")
    diff2010 = calc(proj2010, fun=function(x){ifelse(x >= tssThreshold, 2, 0)}) +
               calc(proj2010, fun=function(x){ifelse(x >= minThreshold, 1, 0)}) 
    diff2050 = calc(proj2050, fun=function(x){ifelse(x >= tssThreshold, 2, 0)}) + 
               calc(proj2050, fun=function(x){ifelse(x >= minThreshold, 1, 0)})
    
    # Load shapefile containing species observed occurences and raster with projections for 
    # year 2010 (current climate).
    shp <- readOGR(dsn=file.path(rootDir,"distribution",spName), layer=spName)
    
    
    # Plot maps to PDF file.
    par(mfrow=c(1,2), oma=c(0, 0, 4, 0))
    
    # maps for 2010.
    minValue = cellStats(diff2010, 'min')
    if(minValue == 0) colorSubset = c(1:3)
    if(minValue == 1) colorSubset = c(2:3)
    if(minValue > 1) stop("error: minimum value should be <= 1.")
    colorScale = c("grey","lightblue","blue")
    plot(diff2010, col=colorScale[colorSubset], legend=F, main="present [2010]")
    points(shp, pch=19, col='red', cex=0.1)
    legend(x=-12, y=85, legend=c("unsuitable","min","tss"), fill=colorScale, horiz=F, inset=0.03, cex=0.8)
    
    # maps for 2050.
    minValue = cellStats(diff2050, 'min')
    if(minValue == 0) colorSubset = c(1:3)
    if(minValue == 1) colorSubset = c(2:3)
    if(minValue > 1) stop("error: minimum value should be <= 1.")
    colorScale = c("grey","orange","darkred")
    plot(diff2050, col=colorScale[colorSubset], legend=F, main=paste("future",ccScenario,"[2050]"))
    points(shp, pch=19, col='red', cex=0.1)
    legend(x=-12, y=85, legend=c("unsuitable","min","tss"), fill=colorScale, horiz=F, inset=0.03, cex=0.8)
    
    mainTitle = paste(spTable[which(spTable$speciesAbrev==spName),'speciesName'], " [",spName,"] \n", 
                      "Difference between reclassified projections with 'TSS' or 'min' thresholds.\n",
                      "min threshold=", minThreshold, "  tss threshold=", tssThreshold, "  points=sp occurences", sep='')
    mtext(mainTitle, side=3, outer=TRUE, cex=1.2)
    

    # Progress report and clear memory.
    cat("###  -> OK.\n")
    rm(proj2010, proj2050, tssThreshold, minThreshold, diff2010, diff2050, mainTitle, shp, minValue, colorSubset)
}

dev.off()

###
########################################################################################################################
    


