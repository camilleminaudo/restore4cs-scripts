##Funtions to calculate gas concentration in the water using headspace technique
##and headspace samples measured through gas cromatography

## Main function
#' Function to calculate the concentration in the water in uM
#' it is the same that GHG_water_uM but you can choose the gas (CO2, CH4 or N2O) for which you want to calculate
#' the concentration
#nGHG_water_uM
nGHG_water_uM <- function(Specie = "CH4", GHG_ppmv, Vol_H2O_HS, Vol_aire_HS, T_Celsius, Patm_eq = 1, R = 0.08206, GHG_atm_ppmv = 1.7){
  if(Specie == "CH4"){
    K <- Kh_CH4(T_Celsius)
  }
  if(Specie == "CO2"){
    K <- Kh_CO2(T_Celsius)
  }
  if(Specie == "N2O"){
    K <- Kh_N2O(T_Celsius)
  }
  nGHG_agua_eq_mol <-nGHG_agua_eq(Patm_eq, GHG_ppmv, Kh_GHG = K, Vol_H2O_HS)
  nGHG_gas_mol <- nGHG_gas(GHG_ppmv, Vol_aire_HS, T_Celsius, Patm_eq, R)
  nGHG_total_mol <- nGHG_total(nGHG_agua_eq = nGHG_agua_eq_mol, nGHG_gas = nGHG_gas_mol)
  nGHG_aire_inicial_mol <- nGHG_aire_inicial(GHG_atm_ppmv, Vol_aire_HS, T_Celsius, Patm_eq, R)
  GHG_uM_agua(nGHG_total_mol, nGHG_aire_inicial_mol, Vol_H2O_HS)
}

#' Arguments
#' Specie = character vector indicating for which gas you want to calculate the concentration. Options are "CO2", "CH4" and "N2O"
#' GHG_ppmv = numeric vector with the concentration in ppmv of the gas of interest in the headspace phase after equilibration. The one you measured with the GC
#' Vol_H2O_HS = A numeric vector with the volume of water (mL) used in the headspace equilibration
#' Vol_aire_HS = numeric vector with volume of air in the headspace
#' T_celsius = numeric vector with temperature in celsius during headspace equilibration
#' GHG_atm_ppmv = numeric vector with the concentration in ppmv of the gas of interest in the headspace phase before equilibration. Set by default as 1.7 ppm NOAA
#' Patm_eq = numeric vector with pressure during equilibration of headspace. Set by default to 1, atmospheric pressure  
#' R = numeric vector the constant of ideal gas

#' 
#' 

## Funtion 1
#' Calculate gas concentration in the water, units uM
#' Description
#' The GHG_uM_agua function is used to calculate gas concentration in the water

#' 
#' The function
GHG_uM_agua <- function(nGHG_total, nGHG_aire_inicial, Vol_H2O_HS){
  (nGHG_total - nGHG_aire_inicial)/ (Vol_H2O_HS/1000)*1000000
}

#' Arguments
#' nGHG_total = A numeric vector of the total number of moles of the gas of interest 
#' present in the sample (water + headspace)
#' nGHG_aire_inicial = A numeric vector of the total number of moles of gas of interest already present 
#' present in the sample before the equilibration is reached. If N2 is used for headspace equilibration then nGHG_aire_inicial will be 0
#' Vol_H2O_HS = A numeric vector with the volume of water (mL) used in the headspace equilibration

## Function 2

#' Calculate the total number of moles of the gas of intereset present in the sample (water + headspace)

#' The function
nGHG_total <- function(nGHG_agua_eq, nGHG_gas){
  nGHG_agua_eq + nGHG_gas
}

#' Arguments
#' nGHG_aqua_eq = numeric vector with the number of moles of the gas of interest present in the equilibrated water
#' nGHG_agua = numeric vector with the number of moles of the gas of interest present in the equilibrated headspace

##Funtion 3
#'  Calculate the number of moles of the gas of interest that there is in the water sample
nGHG_agua_eq <- function(Patm_eq, GHG_ppmv, Kh_GHG, Vol_H2O_HS){
  (Patm_eq * GHG_ppmv/1000000)*Kh_GHG*(Vol_H2O_HS/1000)
}

#' Arguments
#' Patm_eq = numeric vector with the pressure at which the headspace was equilibrated
#' GHG_ppmv = numeric vector with the concentration in ppmv of the gas of interest in the headspace phase
#' kh = numeric vector with the constant of solubility of the gas of interest in mol/l*atm
#' Vol_H2O_HS = numeric vector with the volume of water (ml) in the headspace equilibration
#' 
#' ##Funtions 4
#' #' This function calculate the constant of solubility for each gas (CO2, CH4 and N2O) based on temperature of equilibrium
  #' Kh_CH4----
  Kh_CH4 <- function(T_Celsius){
    (0.000014*exp(1600*(1/(T_Celsius+273.15)-1/298.15)))*101325/1000
  }
  
  #' Kh_CO2----
  Kh_CO2 <- function(T_Celsius){
    (0.00033*exp(2400*(1/(T_Celsius+273.15)-1/298.15)))*101325/1000
  }
  
  #' Kh_N2O----
  Kh_N2O <- function(T_Celsius){
    (0.00024*exp(1600*(1/(T_Celsius+273.15)-1/298.15)))*101325/1000}
  
  #' Arguments
  #' T_celsius = numeric vector with temperature in celsius during headspace equilibration

## Funtion 5
#' Funtion to calculate number of moles in the headspace phase after equilibration
#' nGHG_gas----
  nGHG_gas <- function(GHG_ppmv, Vol_aire_HS,T_Celsius, Patm_eq = 1, R = 0.08206){
    (Patm_eq*GHG_ppmv/1000000)*(Vol_aire_HS/1000)/R/(T_Celsius+273.15)
  }
#' Patm_eq = numeric vector with the pressure at which the headspace was equilibrated
#' GHG_ppmv = numeric vector with the concentration in ppmv of the gas of interest in the headspace phase after equilibration
#' Vol_aire_HS = numeric vector with volume of air in the headspace
#' T_Celsius = numeric vector with the temperature in Celisus during equilibration of headspace
#' Patm_eq = numeric vector with pressure during equilibration of headspace. Set by default to 1, atmospheric pressure  
#' R = numeric vector the constant of ideal gas
  
#' Function to calculate the number of moles in the headspace phase before starting equilibration (in case you air, not N2)  
  ##nGHG_aire_inicial----
  nGHG_aire_inicial <- function(GHG_atm_ppmv = 1.7, Vol_aire_HS, T_Celsius, Patm_eq = 1, R = 0.08206){
    (Patm_eq*GHG_atm_ppmv/1000000)*(Vol_aire_HS/1000)/R/(T_Celsius+273.15)
  }
#' Patm_eq = numeric vector with the pressure at which the headspace was equilibrated
#' GHG_atm_ppmv = numeric vector with the concentration in ppmv of the gas of interest in the headspace phase before equilibration. Set by default as 1.7 ppm NOAA
#' Vol_aire_HS = numeric vector with volume of air in the headspace
#' T_Celsius = numeric vector with the temperature in Celisus during equilibration of headspace
#' Patm_eq = numeric vector with pressure during equilibration of headspace. Set by default to 1, atmospheric pressure  
#' R = numeric vector the constant of ideal gas
  
