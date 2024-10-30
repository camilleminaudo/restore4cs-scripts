## Functions to calculate gas concentration in water using the headspace technique
## and headspace samples measured through gas chromatography

## Main Function
#' Calculates the concentration of a specified gas (CO2, CH4, or N2O) in water in µM
#' This function is similar to GHG_water_uM but allows the user to choose the gas.
nGHG_water_uM <- function(Specie = "CH4", GHG_ppmv, Vol_H2O_HS, Vol_air_HS, T_Celsius, Patm_eq = 1, R = 0.08206, GHG_atm_ppmv = 1.7) {
  if (Specie == "CH4") {
    K <- Kh_CH4(T_Celsius)
  }
  if (Specie == "CO2") {
    K <- Kh_CO2(T_Celsius)
  }
  if (Specie == "N2O") {
    K <- Kh_N2O(T_Celsius)
  }
  nGHG_water_eq_mol <- nGHG_water_eq(Patm_eq, GHG_ppmv, Kh_GHG = K, Vol_H2O_HS)
  nGHG_gas_mol <- nGHG_gas(GHG_ppmv, Vol_air_HS, T_Celsius, Patm_eq, R)
  nGHG_total_mol <- nGHG_total(nGHG_water_eq = nGHG_water_eq_mol, nGHG_gas = nGHG_gas_mol)
  nGHG_air_initial_mol <- nGHG_air_initial(GHG_atm_ppmv, Vol_air_HS, T_Celsius, Patm_eq, R)
  GHG_uM_water(nGHG_total_mol, nGHG_air_initial_mol, Vol_H2O_HS)
}

#' Arguments
#' Specie = character vector indicating which gas to calculate concentration for. Options: "CO2", "CH4", "N2O"
#' GHG_ppmv = numeric vector with the concentration in ppmv of the gas of interest in the headspace phase after equilibration (measured via GC)
#' Vol_H2O_HS = numeric vector with the volume of water (mL) in the headspace equilibration
#' Vol_air_HS = numeric vector with the volume of air in the headspace
#' T_Celsius = numeric vector with the temperature (Celsius) during equilibration
#' GHG_atm_ppmv = numeric vector with initial gas concentration in ppmv before equilibration (default 1.7 ppm, NOAA)
#' Patm_eq = numeric vector for atmospheric pressure during equilibration (default 1)
#' R = numeric vector for the ideal gas constant

## Function 1
#' Calculates gas concentration in water in µM
#' GHG_uM_water function is used to calculate gas concentration in water
GHG_uM_water <- function(nGHG_total, nGHG_air_initial, Vol_H2O_HS) {
  (nGHG_total - nGHG_air_initial) / (Vol_H2O_HS / 1000) * 1e6
}

#' Arguments
#' nGHG_total = numeric vector for the total moles of the gas of interest in the sample (water + headspace)
#' nGHG_air_initial = numeric vector for moles of gas present in the sample before equilibration (0 if N2 is used)
#' Vol_H2O_HS = numeric vector for the volume of water in headspace equilibration (mL)

## Function 2
#' Calculates the total moles of gas of interest in the sample (water + headspace)
nGHG_total <- function(nGHG_water_eq, nGHG_gas) {
  nGHG_water_eq + nGHG_gas
}

#' Arguments
#' nGHG_water_eq = numeric vector for moles of gas in the equilibrated water
#' nGHG_gas = numeric vector for moles of gas in the equilibrated headspace

## Function 3
#' Calculates moles of gas in the water sample
nGHG_water_eq <- function(Patm_eq, GHG_ppmv, Kh_GHG, Vol_H2O_HS) {
  (Patm_eq * GHG_ppmv / 1e6) * Kh_GHG * (Vol_H2O_HS / 1000)
}

#' Arguments
#' Patm_eq = numeric vector for pressure at equilibration
#' GHG_ppmv = numeric vector for concentration in ppmv of the gas of interest in headspace
#' Kh_GHG = solubility constant of the gas of interest in mol/L*atm
#' Vol_H2O_HS = numeric vector for water volume in headspace equilibration (mL)

## Functions 4
#' Calculates solubility constant for each gas based on equilibrium temperature
Kh_CH4 <- function(T_Celsius) {
  (0.000014 * exp(1600 * (1 / (T_Celsius + 273.15) - 1 / 298.15))) * 101325 / 1000
}

Kh_CO2 <- function(T_Celsius) {
  (0.00033 * exp(2400 * (1 / (T_Celsius + 273.15) - 1 / 298.15))) * 101325 / 1000
}

Kh_N2O <- function(T_Celsius) {
  (0.00024 * exp(1600 * (1 / (T_Celsius + 273.15) - 1 / 298.15))) * 101325 / 1000
}

#' Arguments
#' T_Celsius = numeric vector with the temperature in Celsius during equilibration

## Function 5
#' Calculates moles in headspace phase after equilibration
nGHG_gas <- function(GHG_ppmv, Vol_air_HS, T_Celsius, Patm_eq = 1, R = 0.08206) {
  (Patm_eq * GHG_ppmv / 1e6) * (Vol_air_HS / 1000) / R / (T_Celsius + 273.15)
}

#' Arguments
#' Patm_eq = pressure at equilibration
#' GHG_ppmv = gas concentration in ppmv after equilibration
#' Vol_air_HS = air volume in the headspace
#' T_Celsius = temperature in Celsius during equilibration
#' R = ideal gas constant

## Function 6
#' Calculates moles in headspace phase before equilibration
nGHG_air_initial <- function(GHG_atm_ppmv = 1.7, Vol_air_HS, T_Celsius, Patm_eq = 1, R = 0.08206) {
  (Patm_eq * GHG_atm_ppmv / 1e6) * (Vol_air_HS / 1000) / R / (T_Celsius + 273.15)
}

#' Arguments
#' GHG_atm_ppmv = initial gas concentration in ppmv before equilibration (default 1.7 ppm NOAA)
#' Vol_air_HS = volume of air in the headspace
#' T_Celsius = temperature in Celsius during equilibration
#' Patm_eq = atmospheric pressure during equilibration (default 1)
#' R = ideal gas constant
