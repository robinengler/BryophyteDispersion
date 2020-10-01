####################################################################################################
### Setup input data.
### ****************
### Create input data structure for migclim runs of the 'bryophyte' project. In this structure, 
### each species has its own directory with all the input data needed for its migclim runs.
### Input raster data are taken from their origin location, converted to ascii format, renamed, 
### and then saved in the correct directory.
###
####################################################################################################

####################################################################################################
### Input and environment setup.
### ***************************
species_list = c('Amphl','Amphm','Anas','Anom','Anth','Arct','Atri','Bart','Bazz','Cors','Cyno',
                 'Cyrt','Dicr','Dipla','Diplt','Fabr','Foss','Frul','Glyph','Grim','Gymn','Habr',
                 'Herb','Homa','Lept','Leuc','Mast','Metz','Myur','Ortha','Orthl','Oxym','Palu',
                 'Plag','Ptych','Sacc','Scle','Scor','Spha','Ulot')

cc_scenario_list     = c('mp45', 'mp85', 'he45', 'he85')
wind_scenario_list   = c('maxW', 'meanW')
release_height_list  = c('h03', 'h1', 'h10')
env_change_year_list = c(2020, 2030, 2040, 2050)
kernel_distance_list = c(1:10)
ldd_freq_list        = c('0', '0.0001', '0.001', '0.01', '0.1')
test_mode            = TRUE

# Input and output directory locations.
root_dir              = '/scratch/temporary/rengler/migclim_run'
original_data_dir     = file.path(root_dir, 'data_original')
output_dir            = file.path(root_dir, 'data_input')
species_table_file    = file.path(output_dir, 'species_table.txt')
simulation_table_file = file.path(output_dir, 'simulation_table.txt')


# Command line agruments (optional).
# *********************************
# The number of the species to process and the test_mode value (TRUE/FALSE) can be passed to the
# script as command line arguments. If nothing is passed, then the script defaults to using all
# species and test_mode = TRUE.
args = commandArgs(trailingOnly=TRUE)    # commandArgs() retrieves all input arguments that were
                                         # passed to the script with 'Rscript'.
if(length(args) > 0){
    species_to_process = args[[1]]       # e.g. '1,2,5'
    test_mode = as.logical(args[[2]])    # either 'TRUE' or 'FALSE'.
    
    # Convert species to process input to numeric values.
    if(species_to_process != 'all'){
        species_to_process = as.numeric(unlist(strsplit(species_to_process, split=',')))
        if(any(is.na(species_to_process))) stop('invalid input species number passed.\n')
        if(any(species_to_process > length(species_list))) stop('input species nb out of range.\n')
        species_list = species_list[species_to_process]
    }
    
    # Verify test mode value is TRUE/FALSE.
    if(is.na(test_mode)) stop('invalid test mode value passed.\n')
    
}

    
# Progress report.
# ***************
cat('################################################################################ \n')
cat('### Setup input data script. \n')
cat('### ***********************  \n')
cat('### Number of species:', length(species_list), '\n')
cat('### original data dir:', original_data_dir, '\n')
cat('### output dir       :', output_dir, '\n')
cat('### Test mode        :', test_mode, '\n')
cat('### species list     :', paste(species_list, collapse=', '), '\n')
cat('### \n')
Sys.sleep(5)


# Verify all input/output directories exist.
if(!dir.exists(original_data_dir)) stop('original data dir does not exist')
if(!dir.exists(output_dir)) stop('output dir does not exist')

# Load required libraries.
require(raster, quietly=TRUE)
####################################################################################################



####################################################################################################
### Generate "species_table.txt" file.
### *********************************
### This text file contains the list of all species and their associated max TSS reclassification
### thresholds.
cat('### Generate species data table: \n')
if(file.exists(species_table_file)){
    cat('###  -> load existing file. \n')
    species_table = read.table(species_table_file, h=T, sep='\t', as.is=T)
    
} else{
    # Create new table.
    cat('###  -> create new file. \n')
    species_table = data.frame(species = species_list, 
                               max_TSS_threshold_GLM = NA, 
                               max_TSS_threshold_GBM = NA, 
                               max_TSS_threshold = NA)
    
    # Loop through all species and get their TSS reclassification thresholds.
    for(species in species_list){
        # Get TSS threshold for GLM model
        input_file = paste0(original_data_dir,'/sdm_evaluation/',species,'_TSS_GLM.txt')
        tss_glm = read.table(input_file, h=T, sep='\t', as.is=T)['Cutoff', 'Full']
        # Note: the same values can also be retrieved from different files:
        #input_file = paste0(original_data_dir,'/sdm_evaluation/',species,'_eval_assemble_GLM.txt')
        #cutoff_colname = paste0(species, 'GLM_EMwmeanByTSS_mergedAlgo_Full_PA1.Cutoff')
        #tss_glm = read.table(input_file, h=T, sep='\t', as.is=T)['TSS', cutoff_colname]
        
        # Get TSS threshold for GBM model
        input_file = paste0(original_data_dir,'/sdm_evaluation/',species,'_TSS_GBM.txt')
        tss_gbm = read.table(input_file, h=T, sep='\t', as.is=T)['Cutoff', 'Full']
        # Note: the same values can also be retrieved from different files:
        #input_file = paste0(original_data_dir,'/sdm_evaluation/',species,'_eval_assemble_GBM.txt')
        #cutoff_colname = paste0(species, 'GBM_EMwmeanByTSS_mergedAlgo_Full_PA1.Cutoff')
        #tss_gbm = read.table(input_file, h=T, sep='\t', as.is=T)['TSS', cutoff_colname]
        #rm(cutoff_colname)
        
        colnb = which(species_table$species == species)
        species_table[colnb, 'max_TSS_threshold_GLM'] = round(tss_glm)
        species_table[colnb, 'max_TSS_threshold_GBM'] = round(tss_gbm)
        species_table[colnb, 'max_TSS_threshold'] = round(mean(c(tss_glm, tss_gbm), na.rm=TRUE))
    }
    
    # Save species_table to disk.
    write.table(species_table, file=species_table_file, sep='\t',quote=F,row.names=F)
    rm(species_table_file, species, input_file, tss_glm, tss_gbm, colnb)
}
cat('###  -> completed. \n')
cat('### \n')
####################################################################################################



####################################################################################################
### Generate "simulation_table.txt" file.
### *************************************
### Table with the list of migclim simulations to run, i.e. all parameter combinations to test.
cat('### Generate species data table: \n')
simulation_table = data.frame('species'        = rep(species_table[,'species'], 
                                                     each=(length(cc_scenario_list) * 
                                                           length(wind_scenario_list) * 
                                                           length(release_height_list) * 
                                                           length(ldd_freq_list)) ),
                              'cc_scenario'    = rep(cc_scenario_list, 
                                                     each=(length(wind_scenario_list) * 
                                                           length(release_height_list) * 
                                                           length(ldd_freq_list)) ), 
                              'wind_speed'     = rep(wind_scenario_list, 
                                                     each=(length(release_height_list) * 
                                                           length(ldd_freq_list)) ),
                              'release_height' = rep(release_height_list, 
                                                     each=length(ldd_freq_list)), 
                              'ldd_freq'       = ldd_freq_list, 
                              stringsAsFactors = FALSE)
simulation_table[,'simulation_name'] = paste0(simulation_table$species, '_', 
                                              simulation_table$cc_scenario, '_', 
                                              simulation_table$wind_speed, '_', 
                                              simulation_table$release_height,
                                              '_ldd',
                                              sapply(simulation_table$ldd_freq, 
                                                     FUN=function(x) which(x == ldd_freq_list)))
stopifnot(nrow(simulation_table) == nrow(unique(simulation_table)))
stopifnot(nrow(simulation_table) == 4800)
write.table(simulation_table, file=simulation_table_file, sep='\t',quote=F,row.names=F)
rm(simulation_table_file, simulation_table)
cat('###  -> completed. \n')
cat('### \n')
####################################################################################################



####################################################################################################
### Prepare dispersal kernel rasters.
### ********************************
#
# Loop through all species,
species_kernel <- read.table(file.path(original_data_dir,"tab_spore.csv"),sep="\t",dec=",",h=T)
distance <- c("500-1500","1500-2500","2500-3500","3500-4500","4500-5500","5500-6500","6500-7500",
              "7500-8500","8500-9500","9500-10500")

cat('### Convert kernel rasters to ascii format:\n')
for(species in species_list){
  
  cat_kernel <- as.character(species_kernel[species_kernel$spName==species,"Kernel"]) #Category for the species kernel
  # Create a subdirectory
  # If needed, create a sub-directory for current species.
  species_dir = file.path(output_dir, species)
  if(!dir.exists(species_dir) & !test_mode) dir.create(species_dir)
  
  # Loop through climate change scenarios, wind scenarios and release heights.
  for(cc_scenario in cc_scenario_list){
    for(wind_scenario in wind_scenario_list){
      for(release_height in release_height_list){
        for(year in env_change_year_list){
          year_index = year/10 - 201
          for(i in 1:10){
            # Progress report.
            cat("\r###  ->", species, cc_scenario, wind_scenario, release_height, year, ':', i)
            
            # Define input/output file names. Check the input file exists.
            input_file = file.path(original_data_dir, 'kernels', cc_scenario, wind_scenario,
                                   release_height, paste0(cat_kernel, '_', distance[i], '_', year_index, '.tif'))
            output_file = paste0(species_dir, '/kernel_', species, '_', cc_scenario, '_',
                                 wind_scenario, '_', release_height, '_', year, '_', i,'.asc')
            if(!file.exists(input_file)) stop('missing input file:', input_file)
            
            # Convert original raster to ascii format and write it to the correct location.
            if(!file.exists(output_file) & !test_mode){
              # Write raster in ascii format to the output directory for the current migclim run.
              input_raster = raster(input_file)
              writeRaster(input_raster, filename=output_file,
                          format='ascii', datatype='FLT4S', overwrite=TRUE)
              rm(input_raster)
            }
          }
          # Print carriage return to start a new line the next time we print something.
          cat(' \n')
          rm(i, input_file, output_file)
        }
      }
    }
  }
}
rm(species, species_dir, cc_scenario, wind_scenario, release_height, year,year_index,species_kernel,cat_kernel)
cat('###  -> completed. \n')
cat('### \n')
####################################################################################################



####################################################################################################
### Prepare habitat suitability rasters.
### ***********************************
### Reclass present (2010) projection to generate 'iniDist', the initial distribution of the
### species (1 = species is present, 0 = species is absent).
### Convert future projections (2020, 2030, 2040, 2050) to ascii grid format without reclassifying
### them (reclassification will be done in the migclim run directly).

cat('### Convert habitat suitability rasters to ascii format:\n')
species_table <- read.table(file.path(output_dir,"species_table.txt"),sep="\t",h=T)
for(species in species_list){
  cat('### -> processing species:', species, '\n')
  species_dir = file.path(output_dir, species)

  # Reclass projections under current climate.
  # *****************************************
  cat('###  -> present \n')

  # Define input/output file names. Check the input file exists.
  input_file = paste0(original_data_dir, '/sdm/', species, '_Present.tif')
  output_file = paste0(species_dir, '/hsmap_', species, '_present.asc')
  if(!file.exists(input_file)) stop('missing input file:', input_file)

  # Get max TSS reclassification threshold for the current species.
  rc_threshold = species_table[which(species_table$species==species), 'max_TSS_threshold']
  stopifnot(length(rc_threshold) == 1)

  # Reclass raster.
  if(!file.exists(output_file) & !test_mode){
    input_raster = raster(input_file)
    tmp = calc(input_raster, fun=function(x){ifelse(x >= rc_threshold, 1, 0)},
               filename=output_file, format='ascii', datatype='INT2S', overwrite=TRUE)
    rm(input_raster, tmp)
  }


  # Reclass projections under climate change scenarios.
  # **************************************************
  for(cc_scenario in cc_scenario_list){
    for(year in env_change_year_list){
      cat('###  ->', cc_scenario, year, ' \n')

      # Define input/output file names. Check the input file exists.
      year_index = year/10 - 201
      input_file = paste0(original_data_dir, '/sdm/', species, '_',
                          year-2000, toupper(cc_scenario), '.tif')
      output_file = paste0(species_dir, '/hsmap_',
                           species, '_', cc_scenario, '_', year_index, '.asc')
      if(!file.exists(input_file)) stop('missing input file:', input_file)

      # If output file is missing, generate it.
      if(!file.exists(output_file) & !test_mode){
        input_raster = raster(input_file)
        NAvalue(input_raster) = -9999
        writeRaster(input_raster, filename=output_file,
                    format='ascii', datatype='INT2S', overwrite=TRUE)
        rm(input_raster)
      }
    }
  }
}
rm(species, species_dir, input_file, output_file, rc_threshold, cc_scenario, year, year_index)
cat('###  -> completed. \n')
cat('### \n')
####################################################################################################


####################################################################################################
### Generate symlinks for years 2060 - 2100.
### ***************************************
### To extend the simulation to from 2050 to 2100, but with no more change in habitat suitability,
### we create symlinks to the data for 2050.
cat('### Generate symlinks for years 2060-2100:\n')
for(species in species_list){
  cat('### -> create symlinks for:', species, '\n')   
  species_dir = file.path(output_dir, species)
  
  for(cc_scenario in cc_scenario_list){
    for(year in c(2060,2070,2080,2090,2100)){
      
      # Symlink for dispersal kernels.
      # *****************************
      for(wind_scenario in wind_scenario_list){
        for(release_height in release_height_list){
          for(i in 1:10){
            target_file = paste0(species_dir, '/kernel_', species, '_', cc_scenario, '_', 
                                 wind_scenario, '_', release_height, '_2050_', i, '.asc')
            symlink_file = paste0(species_dir, '/kernel_', species, '_', cc_scenario, '_', 
                                  wind_scenario, '_', release_height, '_', year, '_', i,'.asc')
            if(!file.exists(symlink_file) & !test_mode){
              if(!file.exists(target_file)) stop('missing input file:', target_file)
              file.symlink(basename(target_file), symlink_file)
            }
          }
        }
      }
      
      # Symlinks for habitat suitability rasters.
      # ****************************************
      year_index = year/10 - 201
      target_file = paste0(species_dir, '/hsmap_', species, '_', cc_scenario, '_4.asc')
      symlink_file = paste0(species_dir, '/hsmap_', 
                            species, '_', cc_scenario, '_', year_index, '.asc')
      if(!file.exists(symlink_file) & !test_mode){
        if(!file.exists(target_file)) stop('missing input file:', target_file)
        file.symlink(basename(target_file), symlink_file)
      }
    }
  }
}
rm(species, species_dir, cc_scenario, year, year_index, 
   wind_scenario, release_height, i, target_file, symlink_file)
cat('###  -> completed. \n')
cat('### \n')
cat('################################################################################ \n')
####################################################################################################

