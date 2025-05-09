#Selection criteria for CH4 water (aquaGHG)


#Chosing best flux will follow this logic: 

#1. If there is water, flux separator method will be calculated

#2. IF No bubles are found (i.e. length difusion = total length-1 )--> GOflux criteria (to determine, similar logic to CO2 criteria, after inspection for hard thresholds for g.fact)

#3. IF Bubles are found (difusion length < total length - 1): 

  #IF ebullition is significant (based on flux +- SD of total flux, difusive flux, ebullitive flux): Use total flux as best flux 

  #IF ebullition is non-significant, or method fails (difusion > total flux):

      #Check Goflux LM quality, if above threshold (to determine), best flux is taken from LM
          #IF LM quality is below threshold, best flux is total flux 


#STEPS: 
  #1. First of all, total flux must be calculated for the whole period of the incubation (now it seems it is only calculated from first point of detected difusive chunk)

  #2. Inpsect scripts aquaGHG to determine significance of ebullition

  #3. Determine g.fact threshold for goflux when no bubles are detected

  #4. Determine LM quality thresholds for LM vs total.flux for non-significant ebullition. 



