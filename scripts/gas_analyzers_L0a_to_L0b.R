
# ---
# Authors: Camille Minaudo, Benjamin Misteli
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script loads a raw measurement file (data level L0a) from one of the gas 
# analyzers used in the project, and transform it into a unified harmonized csv 
# file (data level L0b) allowing for visualization of the data and further data processing.


# Directories
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"
datapath <- paste0(dropbox_root,"/GHG/RAW data")


setwd(datapath)




