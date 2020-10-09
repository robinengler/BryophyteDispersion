#!/software/R/3.5.1/bin/R
####################################################################################################
### Merge migclim simulations output.
### *********************************
### This script merges the migclim "_stats.txt" and "_summary.txt" files for all migclim 
### simulations into a single file.
###
###
####################################################################################################


####################################################################################################
### Input and environment setup.
### ***************************
# Input parameter values.
# **********************
root_dir         = '/scratch/temporary/rengler/migclim_run'
output_dir       = file.path(root_dir, 'data_output')
simulation_file  = file.path(root_dir, 'data_input', 'simulation_table.txt')
replicate_nb     = 30
env_change_steps = 50
dispersal_steps  = 10
year_list = 2009 + (1:((env_change_steps * dispersal_steps) + 1))


# Command line agruments (optional).
# *********************************
test_mode = FALSE
args = commandArgs(trailingOnly=TRUE)    # commandArgs() retrieves all input arguments that were
if(length(args) > 0){
    test_mode = as.logical(args[[1]])   # either 'TRUE' or 'FALSE'.
    if(is.na(test_mode)) stop('invalid test mode value passed.\n')
}


# Load simulation list.
# ********************
if(!dir.exists(root_dir))  stop('root dir does not exist:', input_dir, '\n')
if(!dir.exists(output_dir)) stop('directory does not exist:', output_dir, '\n')
if(!file.exists(simulation_file)) stop('input file does not exist:', species_file, '\n')
simulation_table = read.table(simulation_file, h=T, sep='\t', 
                              as.is=T, stringsAsFactors=F, colClasses='character')
stopifnot(nrow(simulation_table) == 4800)

####################################################################################################


#tmp_subset = c('Amphl','Amphm','Anas','Anom', 'Anth', 'Gymn', 'Habr', 'Herb', 'Palu', 'Plag', 'Ptych', 'Sacc')
#simulation_table = simulation_table[which(simulation_table$species %in% tmp_subset), ]



####################################################################################################
### Generate summary accross all simulations for the "_stats.txt" files.
### *******************************************************************
out_df = simulation_table[rep(1:nrow(simulation_table), each=length(year_list)), ]
out_df[,'year'] = year_list

out_df[,'univDispersal']        = NA
out_df[,'noDispersal']          = NA
out_df[,'occupied_mean']        = NA
out_df[,'occupied_sd']          = NA
out_df[,'absent_mean']          = NA
out_df[,'absent_sd']            = NA

out_df[,'sumColonized_mean']    = NA  
out_df[,'sumColonized_sd']      = NA
out_df[,'sumDecolonized_mean']  = NA  
out_df[,'sumDecolonized_sd']    = NA
out_df[,'sumLDDsuccess_mean']   = NA
out_df[,'sumLDDsuccess_sd']     = NA

out_df[,'stepColonized_mean']   = NA  
out_df[,'stepColonized_sd']     = NA
out_df[,'stepDecolonized_mean'] = NA  
out_df[,'stepDecolonized_sd']   = NA
out_df[,'stepLDDsuccess_mean']  = NA
out_df[,'stepLDDsuccess_sd']    = NA

# Loop through all simulations.
for(i in 1:nrow(simulation_table)){
    
    # Species and simulation name for the current loop.
    cat('### processing simulation', i, '/', nrow(simulation_table), '\n')
    species         = simulation_table[i, 'species']
    simulation_name = simulation_table[i, 'simulation_name']

        
    # Load "_stats" data for all replicates.
    file_list = file.path(output_dir, species, simulation_name, 
                          paste0(simulation_name, '_', 1:replicate_nb,'_stats.txt'))
    missing_count = sum(!sapply(file_list, FUN=file.exists))
    if(missing_count > 0) cat('ERROR: missing', missing_count, 'files for',simulation_name,'\n')
    if(test_mode) next
    stopifnot(missing_count == 0)
    
    tmp = lapply(file_list, FUN=function(x) read.table(x,h=T,sep='\t',as.is=T,stringsAsFactors=F))
    array_stats = array(unlist(tmp), 
                        dim=c(nrow(tmp[[1]]), ncol(tmp[[1]]), length(tmp)), 
                        dimnames=list(year_list, 
                                      colnames(tmp[[1]]), 
                                      paste0('replicate_', 1:replicate_nb)) )
        
    # Compute mean and standard deviation values for each year.
    start_index = (i - 1) * length(year_list) + 1
    end_index   = (i - 1) * length(year_list) + length(year_list)
    sub_df = out_df[start_index:end_index, ]
    for(j in 1:length(year_list)){
        #row_nb = (i - 1) * length(year_list) + j
        year = as.character(out_df$year[j])
        
        tmp = unique(array_stats[year, 'univDispersal', ])
        stopifnot(length(tmp) == 1)
        sub_df[j, 'univDispersal'] = tmp
        
        tmp = unique(array_stats[year, 'NoDispersal', ])
        stopifnot(length(tmp) == 1)
        sub_df[j, 'noDispersal'] = tmp
        
        # Mean and standard deviations of counts.
        for(stat_name in c('occupied', 'absent', 
                             'stepColonized', 'stepDecolonized', 'stepLDDsuccess')){
            tmp = array_stats[year, stat_name, ]
            sub_df[j, paste0(stat_name, '_mean')] = round(mean(tmp), 1)
            sub_df[j, paste0(stat_name, '_sd')]   = round(sd(tmp), 1)
            
        }
        
        # Mean and standard deviation of sums.
        for(stat_name in c('sumColonized', 'sumDecolonized', 'sumLDDsuccess')){
            tmp = array_stats[1:j, sub('sum','step',stat_name), ]
            if(j > 1) tmp = apply(tmp, MARGIN=2, FUN=sum)
            sub_df[j, paste0(stat_name, '_mean')] = round(mean(tmp), 1)
            sub_df[j, paste0(stat_name, '_sd')]   = round(sd(tmp), 1)
        }
        rm(year, tmp, stat_name)
    }
    
    out_df[start_index:end_index, ] = sub_df
    rm(species, simulation_name, file_list, array_stats, sub_df, start_index, end_index)
}

# Save output data frame to disk.
if(! test_mode){
    stopifnot(all(!is.na(out_df[,c(1,2,5:25)])))
    write.table(out_df, file=file.path(output_dir, 'output_summary_stats.txt'), 
                row.names=FALSE, sep='\t', quote=FALSE)
}
####################################################################################################

