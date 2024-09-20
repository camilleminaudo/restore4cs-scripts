#Methane functions====
  ##Concentraci??n de metano en la muestra de agua expresada en uM----
   GHG_uM_agua <- function(nGHG_total, nGHG_aire_inicial, Vol_H2O_HS){
     (nGHG_total - nGHG_aire_inicial)/ (Vol_H2O_HS/1000)*1000000
   }
   
  #nGHG_total n??mero de moles en la muestra de gas analizada.
  #nGHG_inicial n??mero de moles en la muestra de gas de aire.
  #Vol_H2O_HS es el volumen de agua en el HeadSpaces, mililitros.
  
  ##N??mero de moles de metano totales----
   nGHG_total <- function(nGHG_agua_eq, nGHG_gas){
     nGHG_agua_eq + nGHG_gas
   }
  #nGHG_agua_eq n??mero de moles de metano en el agua en equilibrio
  #nGHG_gas n??mero de moles de metano en el gas
  
  ##nGHG_agua_eq----
    nGHG_agua_eq <- function(Patm_eq, GHG_ppmv, Kh_GHG, Vol_H2O_HS){
      (Patm_eq * GHG_ppmv/1000000)*Kh_GHG*(Vol_H2O_HS/1000)
    }
    
  #Patm_eq presi??n atmosf??rica en equilibrio, atm. Generalmente 1 atm
  #CH4_ppm concentraci??n de metano en la muestra en ppmv
  #Kh_CH4 constante de solubilidad del metano, mol/l*atm
  #Vol_H2O_HS volumen de agua en el HeadSpace en ml.
  
  ##Kh_CH4----
    Kh_CH4 <- function(T_Celsius){
      (0.000014*exp(1600*(1/(T_Celsius+273.15)-1/298.15)))*101325/1000
    }
  
  ##Kh_CO2----
    Kh_CO2 <- function(T_Celsius){
      (0.00033*exp(2400*(1/(T_Celsius+273.15)-1/298.15)))*101325/1000
    }
    
  ##Kh_N2O----
    Kh_N2O <- function(T_Celsius){
      (0.00024*exp(1600*(1/(T_Celsius+273.15)-1/298.15)))*101325/1000}
    
  ##nCH4_gas----
      nGHG_gas <- function(GHG_ppmv, Vol_aire_HS,T_Celsius_eq, Patm_eq = 1, R = 0.08206){
      (Patm_eq*GHG_ppmv/1000000)*(Vol_aire_HS/1000)/R/(T_Celsius_eq+273.15)
    }
   
  ##nGHG_aire_inicial----
    nGHG_aire_inicial <- function(GHG_atm_ppmv = 1.7, Vol_aire_HS, T_Celsius_eq, Patm_eq = 1, R = 0.08206){
      (Patm_eq*GHG_atm_ppmv/1000000)*(Vol_aire_HS/1000)/R/(T_Celsius_eq+273.15)
    }

    
##Una funci??n que las recoja a todas====
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
 
 #nCO2_water_uM
 nCO2_water_uM <- function(GHG_ppmv, Vol_H2O_HS, Vol_aire_HS, T_Celsius, Patm_eq = 1, R = 0.08206, GHG_atm_ppmv = 1.7){
   K <- Kh_CO2(T_Celsius)
   nGHG_agua_eq_mol <-nGHG_agua_eq(Patm_eq, CO2_ppmv, Kh = K, Vol_H2O_HS)
   nGHG_gas_mol <- nGHG_gas(CO2_ppmv, Vol_aire_HS, T_Celsius, Patm_eq, R)
   nGHG_total_mol <- nGHG_total(nGHG_agua_eq = nGHG_agua_eq_mol, nGHG_gas = nGHG_gas_mol)
   nGHG_aire_inicial_mol <- nGHG_aire_inicial(GHG_atm_ppmv, Vol_aire_HS, T_Celsius, Patm_eq, R)
   GHG_uM_agua(nGHG_total_mol, nGHG_aire_inicial_mol, Vol_H2O_HS)
 }
 