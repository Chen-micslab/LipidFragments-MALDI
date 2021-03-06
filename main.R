# Demo R script for illustrating the workflow of pre-processing,
# processing and statistics of data generated by mass spectrometers


# Check for dependencies
if (!all(c('alsace', 'ggplot2', 'gridExtra', 'MALDIquant', "mzR") 
         %in% installed.packages()))
    stop(paste0('Please install the following packages before running:', '\n\t',
                 '"alsace", "ggplot2", "gridExtra", "MALDIquant","mzR"'))


# Process the 1-D mass spectra generated by MALDI-TOF
source('process/processMZ_20210910.R')  # DHB + Choline standards
source('process/processMZ_20211024.R')  # DHB + Carnitine standards,
                                        # DHB + d9-Choline + Serum
source('process/processMZ_20211025.R')  # DHB + d3-Carnitine + Serum
source('process/processMZ_20211107.R')  # TiO2 + PE
source('process/processMZ_20211108.R')  # TiO2 + PE(MS2), TiO2 + PG(MS2)
source('process/processMZ_20211109.R')  # TiO2 + PC, TiO2 + PG
source('process/processMZ_20211114.R')  # TiO2 + PE(MS2), TiO2 + PG(MS2)
source('process/processMZ_20211115.R')  # TiO2 + Serum(EtOH)(MS2)
source('process/processMZ_20211123.R')  # DHB + Serum(EtOH)
source('process/processMZ_20211124.R')  # TiO2 + PC(MS2), TiO2 + PE(MS2), 
                                        # TiO2 + PG(MS2)
source('process/processMZ_20211125.R')  # DHB + Serum(MeOH-CHCl3)(MS2),
                                        # TiO2 + PC(MS2),
                                        # TiO2 + Serum(EtOH)(MS2)
source('process/processMZ_20211126.R')  # TiO2 + Serum(MeOH-CHCl3)
source('process/processMZ_20211127.R')  # TiO2 + PC(MS2),
source('process/processMZ_20211128.R')  # DHB + Standards(MS2), DHB + Serum(MS2)
source('process/processMZ_20211207.R')  # DHB + Standards(MS2)
source('process/processMZ_20211213.R')  # TiO2 + Serum(EtOH)
source('process/processMZ_20211215.R')  # TiO2 + Serum(MeOH-CHCl3)
source('process/processMZ_20211221.R')  # TiO2 + Standards(MS2), 
                                        # TiO2 + Serum(EtOH)
source('process/processMZ_20220210.R')  # TiO2 + Serum(ACN, EtOH)
source('process/processMZ_20220211.R')  # TiO2 + PC
source('process/processMZ_20220212.R')  # TiO2 + PC(MS2)


# Pre-process and process the 2-D mass spectra generated by HPLC-MS
# Skip them if it takes a long time...
source('process/processLCMZ_20210917.R')  # SPCS + Choline standards, 
                                            # choline@Serum(EtOH)
source('process/processLCMZ_20210918.R')  # SPCS + Choline standards
source('process/processLCMZ_20211017.R')  # Carnitine@Serum
source('process/processLCMZ_20211019.R')  # SPCS + Carnitine standards
source('process/processLCMZ_20211025.R')  # Choline@Serum(MeOH-CHCl3)


# Summarize processed spectra and peak lists
source('summary/Summary-CommonMZ_20210910.R')  # Choline, by MALDI
source('summary/Summary-CommonMZ_20210917.R')  # Choline, by HPLC-MS
source('summary/Summary-CommonMZ_20211017.R')  # Carnitine, by HPLC-MS
source('summary/Summary-CommonMZ_20211024.R')  # Carnitine, by MALDI
source('summary/Summary-CommonMZ_20211107.R')  # TiO2 + PE
source('summary/Summary-CommonMZ_20211109.R')  # TiO2 + PC, TiO2 + PG
source('summary/Summary-CommonMZ_20211123.R')  # DHB + Serum(EtOH)
source('summary/Summary-CommonMZ_20211126.R')  # DHB + Serum(MeOH-CHCl3)
source('summary/Summary-CommonMZ_20211213.R')  # TiO2 + Serum(EtOH)
source('summary/Summary-CommonMZ_20211215.R')  # TiO2 + Serum(MeOH-CHCl3)
source('summary/Summary-CommonMZ_20220210.R')  # TiO2 + Serum(ACN, MeOH)
source('summary/Summary-CommonMZ_20220211.R')  # TiO2 + PC


# Establish calibration equations for target metabolites
source('summary/Summary-Linearity_20210910.R')  # Choline, by MALDI
source('summary/Summary-Linearity_20210917.R')  # Choline, by HPLC-MS
source('summary/Summary-Linearity_20211019_2.R')  # Carnitine, by HPLC-MS
source('summary/Summary-Linearity_20211024.R')  # Carnitine, by MALDI


# Calculate metabolite concentration using the calibration equations 
source('summary/Summary-Concentration_20210917.R')  # Choline, by HPLC-MS
source('summary/Summary-Concentration_20211017.R')  # Carnitine, by HPLC-MS
source('summary/Summary-Concentration_20211019.R')  # Carnitine, by MALDI
source('summary/Summary-Concentration_20211024.R')  # Choline, by MALDI


# Plot spectra
# MS/MS spectra make take longer time to plot
source('plot/Plot-Intensity-20211109.R')  # TiO2 + PC, TiO2 + PE, TiO2 + PG
source('plot/Plot-Intensity-20211126.R')  # DHB + Serum
source('plot/Plot-Intensity-20211213.R')  # TiO2 + PC, TiO2 + Serum(EtOH)
source('plot/Plot-Intensity-20211215.R')  # TiO2 + Serum
source('plot/Plot-Intensity-20220210.R')  # TiO2 + PC, TiO2 + PE, TiO2 + PG
source('plot/Plot-Intensity-20220211.R')  # TiO2 + PC, TiO2 + Serum(EtOH)
source('plot/Plot-Intensity-20220211-2.R')  # TiO2 + PC
                                            # TiO2 + Serum(ACN)
                                            # TiO2 + Serum(MeOH)
source('plot/Plot-Intensity_MS2-20211125.R')  # TiO2 + PE(MS2), TiO2 + PG(MS2)
source('plot/Plot-Intensity_MS2-20211207.R')  # DHB + Standards(MS2), 
                                              # DHB + Serum(MS2)
source('plot/Plot-Intensity_MS2-20211221.R')  # TiO2 + PC(MS2)
source('plot/Plot-Intensity_MS2-20211221_2.R')  # TiO2 + PC(MS2), 
                                                # TiO2 + Serum(EtOH)(MS2)
source('plot/Plot-Intensity_MS2-20220210.R')  # TiO2 + PC(MS2)
source('plot/Plot-Intensity_MS2-20220212.R')  # TiO2 + PC(MS2), 
                                              # TiO2 + Serum(EtOH)(MS2)
