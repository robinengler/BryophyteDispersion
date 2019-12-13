#!/software/R/3.5.1/bin/R
####################################################################################################
### Run migclim simulations on for the "Bryophyte" data (Florian, Flavien and Alain's data).
### ***************************************************************************************
###
###
####################################################################################################


####################################################################################################
### Input and environment setup.
### ***************************
# Input data location and range of parameters to test.
root_dir            = '/scratch/temporary/rengler/migclim_run'
input_dir           = file.path(root_dir, 'data_input')
if(startsWith(system('hostname', intern=T), 'rserv')) input_dir = '/scratch/local/monthly/rengler/data_input'
if(startsWith(system('hostname', intern=T), 'dee-serv')) input_dir = '/scratch/local/monthly/rengler/data_input'
output_dir          = file.path(root_dir, 'data_output')
output_dir_rasters  = file.path(output_dir, '01_raster_files')
species_file        = file.path(input_dir, 'species_table.txt')
simulation_file     = file.path(input_dir, 'simulation_table.txt')
test_mode           = TRUE

# migclim parameters.
env_change_steps   = 50     # number of environmental change steps.
dispersal_steps    = 10     # number of dispersal steps.
dispersal_distance = 10     # maximum dispersal distance, in cell units.
ldd_min_dist       = 11     # minimum long distance dispersal distance, in cell units.
ldd_max_dist       = 3000   # maximum long distance dispersal distance, in cell units.
replicate_nb       = 30     # number of replicates to perform for each migclim simulation.


# Load species list from species table.
# ************************************
if(!dir.exists(input_dir))  stop('root dir does not exist:', input_dir, '\n')
if(!dir.exists(output_dir)) stop('directory does not exist:', output_dir, '\n')
if(!dir.exists(output_dir_rasters)) stop('directory does not exist:', output_dir_rasters, '\n')
if(!file.exists(species_file)) stop('input file does not exist:', species_file, '\n')
if(!file.exists(simulation_file)) stop('input file does not exist:', species_file, '\n')
species_table = read.table(species_file, h=T, sep='\t', as.is=T, stringsAsFactors=F)
stopifnot(nrow(species_table) == 40)
simulation_table = read.table(simulation_file, h=T, sep='\t', 
                              as.is=T, stringsAsFactors=F, colClasses='character')

# Remove wind speed and release height columns from table, because we don't vary these parameters
# in the simple kernel simulations.
simulation_table = simulation_table[, c('species', 'cc_scenario', 'ldd_freq', 'simulation_name')]
simulation_table[,'simulation_name'] = gsub('_maxW|_meanW|_h03|_h1|_h10', 
                                            '', simulation_table$simulation_name)
simulation_table[,'simulation_name'] = paste0(simulation_table[,'simulation_name'],'_simplekernel')
simulation_table = unique(simulation_table)
simulation_list = 1:nrow(simulation_table)
stopifnot(length(simulation_list) == 400)


# Command line agruments (optional).
# *********************************
# The index of the simulation(s) to run (one or more numbers in the range 1 - 2400) and the test 
# mode value (TRUE/FALSE) can be passed to the script as command line arguments. If nothing is 
# passed, then the script defaults to running all simulations in test mode.
args = commandArgs(trailingOnly=TRUE)    # commandArgs() retrieves all input arguments that were
                                         # passed to the script with 'Rscript'.
if(length(args) > 0){
    simulation_index = args[[1]]        # e.g. '1,2,5'
    test_mode = as.logical(args[[2]])   # either 'TRUE' or 'FALSE'.
    
    # Convert species to process input to numeric values.
    if(simulation_index != 'all'){
        simulation_index = as.numeric(unlist(strsplit(simulation_index, split=',')))
        if(any(is.na(simulation_index))) stop('invalid input species number passed.\n')
        if(any(simulation_index < 1)) stop('input species nb out of range.\n')
        if(any(simulation_index > nrow(simulation_table))) stop('input species nb out of range.\n')
        simulation_list = simulation_list[simulation_index]
    }
    # Verify test mode value is TRUE/FALSE.
    if(is.na(test_mode)) stop('invalid test mode value passed.\n')
}


# Progress report.
# ***************
cat('################################################################################ \n')
cat('### Setup input data script. \n')
cat('### ***********************  \n')
cat('### Number of simulations:', length(simulation_list), '\n')
cat('### input data dir       :', input_dir, '\n')
cat('### Test mode            :', test_mode, '\n')
cat('### \n')
Sys.sleep(3)

# Load required libraries.
require(raster, quietly=TRUE)
require(SDMTools, quietly=TRUE, warn.conflicts=FALSE)
require(MigClim, lib.loc = '/Home/rengler/R/x86_64-koji-linux-gnu-library/3.5/', quietly=TRUE)
####################################################################################################



####################################################################################################
### Verify all input data for the current run are present.
### *****************************************************
cat('### Input data check: \n')
for(simulation_index in simulation_list){
    # Retrieve simulation parameter values.
    species        = simulation_table[simulation_index, 'species']
    cc_scenario    = simulation_table[simulation_index, 'cc_scenario']
    #wind_speed     = simulation_table[simulation_index, 'wind_speed']
    #release_height = simulation_table[simulation_index, 'release_height']
    cat('###  ->', species, cc_scenario, '\n')
    
    # Set working directory for current simulation.
    species_dir = file.path(input_dir, species)
    if(!dir.exists(species_dir)) stop('missing input directory:', species_dir, '\n')
    setwd(species_dir)
    
    # Check that species initial distribution is present.
    input_file = paste0('hsmap_', species, '_present.asc')
    if(!file.exists(input_file)) stop('missing input file:', input_file, '\n')
    
    # Check that habitat suitability maps are present.
    # Note that there are only 4 actual habitat suitability files, because for the other 
    # environmental change steps we just keep re-using the data for 2050.
    for(i in 1:4){
        input_file = paste0('hsmap_', species, '_', cc_scenario, '_', i, '.asc')
        if(!file.exists(input_file)) stop('missing input file:', input_file, '\n')
    }
    
    # Check that disersal kernel files are present.
    for(i in 1:4){
        input_file = paste0('simple_kernel_hsmap_', species, '_', cc_scenario, '_', i, '.txt')
        if(!file.exists(input_file)) stop('missing input file:', input_file, '\n')
    }
    rm(species, cc_scenario, species_dir, input_file, i)
}
cat('###  -> completed. \n')
cat('### \n')
####################################################################################################



####################################################################################################
### Run migclim.
### ***********
for(simulation_index in simulation_list){
    # Retrieve simulation parameter values.
    species         = simulation_table[simulation_index, 'species']
    cc_scenario     = simulation_table[simulation_index, 'cc_scenario']
    #wind_speed      = simulation_table[simulation_index, 'wind_speed']
    #release_height  = simulation_table[simulation_index, 'release_height']
    ldd_frequency   = simulation_table[simulation_index, 'ldd_freq']
    simulation_name = simulation_table[simulation_index, 'simulation_name']
    tss_threshold   = species_table[which(species_table$species==species),'max_TSS_threshold']
    
    # Set working directory for current simulation, and if needed, delete existing output 
    # directory for the current run.
    setwd(file.path(input_dir, species))
    if(dir.exists(simulation_name)) unlink(simulation_name, recursive=TRUE)
    
    
    # Run migclim simulation.
    # **********************
    cat('### Starting migclim simulation:', simulation_name, '\n')
    cat('### Test mode:', test_mode, '\n')
    N = MigClim.migrate(iniDist         = paste0('hsmap_', species, '_present.asc'), 
                        hsMap           = paste0('hsmap_', species, '_', cc_scenario, '_'), 
                        rcThreshold     = tss_threshold, 
                        envChgSteps     = env_change_steps, 
                        dispSteps       = dispersal_steps, 
                        dispKernel      = rep(1,10), 
                        barrier         = NULL, 
                        barrierType     = 'strong', 
                        iniMatAge       = 1, 
                        propaguleProd   = 1, 
                        lddFreq         = as.numeric(ldd_frequency), 
                        lddMinDist      = ldd_min_dist, 
                        lddMaxDist      = ldd_max_dist, 
                        simulName       = simulation_name, 
                        replicateNb     = replicate_nb, 
                        overWrite       = TRUE, 
                        testMode        = test_mode, 
                        fullOutput      = FALSE, 
                        keepTempFiles   = FALSE, 
                        randomGeneratorSeed = 10 + simulation_index)
    
    
    # Re-organize output files.
    # ************************
    if(!test_mode){
        # Delete some of the raster outputs to save space.
        cat('###  -> delete some of replicated rasters. \n')
        if(replicate_nb > 5){
            for(i in 6:replicate_nb) unlink(file.path(simulation_name, 
                                            paste0(simulation_name, '_', i, '_raster.asc')))
        }
        # Move outputs to output directory. Text and raster files go into separate directories.
        cat('###  -> move data to output directory. \n')
        output_dir_text   = file.path(output_dir, species, simulation_name)
        output_dir_raster = file.path(output_dir_rasters, species, simulation_name)
        if(!dir.exists(output_dir_text)) dir.create(output_dir_text, recursive=T)
        if(!dir.exists(output_dir_raster)) dir.create(output_dir_raster, recursive=T)
        for(f in list.files(simulation_name, full.names=F)){
            if(endsWith(f,'.asc')){
                file.copy(from=file.path(simulation_name, f), to=file.path(output_dir_raster, f))
            } else{
                file.copy(from=file.path(simulation_name, f), to=file.path(output_dir_text, f))
            }
            file.remove(file.path(simulation_name, f))
        } 
        # Delete now empty directory.
        file.remove(simulation_name)
        rm(output_dir_text, output_dir_raster, f)
    }
    
    # Clean environment.    
    rm(species, cc_scenario, ldd_frequency, tss_threshold, simulation_name, N)
}


# Progress report.
# ***************
cat('###  -> completed. \n')
cat('### \n')
cat('### migclim run completed successfully. \n')
cat('################################################################################ \n')
####################################################################################################


