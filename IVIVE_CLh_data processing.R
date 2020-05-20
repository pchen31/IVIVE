## Load packages 
pacman::p_load(pacman, dplyr, readxl, writexl, lubridate) 

## Import raw data (alternatively use GUI: Right hand panel > Environment > Import Dataset > Choose "from excel" > Browse
Rawdata <- read_excel("CLint_Raw data.xlsx", sheet = "Rawdata")
#View(Rawdata)

Df_SF_Qh <- read_excel("CLint_Raw data.xlsx", sheet = "SF_Qh")
#View(Df_SF_Qh)

Df_Exp_Unit <- read_excel("CLint_Raw data.xlsx", sheet = "Exp_Unit")
#View(Df_Exp_Unit)

### Set variables for experimental settings and columns in raw data
Val_protein_conc <- as.numeric(Df_Exp_Unit$Value[1])                 # Protein concentration used in microsomal stability assay
Val_cell_density <- as.numeric(Df_Exp_Unit$Value[2])                 # Cell density used in hepatocyte stability assay
Val_cpd_conc <- as.numeric(Df_Exp_Unit$Value[3])                     # Substrate concentration used in the assay
Val_fu_inc_digit <- as.numeric(Df_Exp_Unit$Value[4])                 # Number of decimal places for fu,inc values
Val_CLint_s_digit <- as.numeric(Df_Exp_Unit$Value[5])                # Number of decimal places for fu,inc values
Val_CLh_digit <- as.numeric(Df_Exp_Unit$Value[6])                    # Number of decimal places for fu,inc values
Val_pH_plasma <- as.numeric(Df_Exp_Unit$Value[11])                   # pH of plasma, for F1 calculation
Val_pH_heps <- as.numeric(Df_Exp_Unit$Value[12])                     # pH of hepatoctytes (inside the cells), for F1 calculation
Val_CLh_unit_conv <- ifelse(Df_Exp_Unit$Unit[10] == "L/hr/kg",       # Value to convert CLh to CLint based on desired CLh unit (L/hr/kg or mL/min/kg)
                            60/1000000, 1/1000)

Num_cpd = nrow(Rawdata)                                              # Calculate number of compounds in the Rawdata, for data processing

Col_ID <- Rawdata$ID
Col_species <- Rawdata$Species
Col_matrix <- Rawdata$Matrix
Col_clint <- Rawdata$CLint
Col_CLobs <- Rawdata$CLobs
Col_exp_fu_p <- Rawdata$exp_fup
Col_c_fup <- Rawdata$c_fup
Col_exp_fu_mic <- Rawdata$exp_fu_mic
Col_exp_fu_hep <- Rawdata$exp_fu_hep
Col_exp_logP <- Rawdata$exp_LogP
Col_c_logP <- Rawdata$c_LogP
Col_exp_logD <- Rawdata$exp_LogD
Col_c_logD <- Rawdata$c_LogD
Col_exp_pKa <- Rawdata$ex_pKa
Col_c_pKa <- Rawdata$c_pKa
Col_pKa_type <- Rawdata$pKa_Type

# Get week, month and year of the metabolic stability experiment
Rawdata <- mutate(Rawdata, Exp_Week = week(Rawdata$Exp_Date),
                           Exp_Month_digit = month(Rawdata$Exp_Date),
                           Exp_Month = months(Rawdata$Exp_Date, abbreviate = TRUE),
                           Exp_Year = year(Rawdata$Exp_Date))

Rawdata$Exp_Date <- as_date(Rawdata$Exp_Date)               # Convert exp date back to only yyyy-mm-dd

## Create function for different fu,inc (i.e. fu,mic & fu,hep) prediction methods, input = logP/D
cfu_hep_Austin <- function(LogP_D, Digits) { zval <- (1/(10^(0.4*LogP_D-1.38)+1)) 
                                             zval <- round(zval, digits = Digits)
                                             return(zval)} 

cfu_hep_Austin.fumic <- function(LogP_D, Protein_conc, Digits) { zval <- (1/(Protein_conc*10^(0.56*LogP_D-1.41)+1))
                                                   zval <- (1/(10^((log10((1-zval)/zval)-0.06)/1.52)+1))
                                                   zval <- round(zval, digits = Digits)
                                                   return(zval)} 

cfu_hep_Kilford <- function(LogP_D, Digits) { zval <- (1/(125*0.005*10^(0.072*LogP_D^2+0.067*LogP_D-1.126)+1))      # VR value is 0.005 for 1M cells/mL incubation
                                              zval <- round(zval, digits = Digits) 
                                              return(zval)} 

cfu_mic_Austin <- function(LogP_D, Protein_conc, Digits) { zval <- (1/(Protein_conc*10^(0.56*LogP_D-1.41)+1)) 
                                             zval <- round(zval, digits = Digits) 
                                             return(zval)} 

cfu_mic_Houston <- function(LogP_D, Protein_conc, Digits) { zval <- (1/(Protein_conc*10^(0.072*LogP_D^2+0.067*LogP_D-1.126)+1)) 
                                              zval <- round(zval, digits = Digits) 
                                              return(zval)} 

cfu_mic_Turner.basic <- function(LogP, Protein_conc, Digits) { zval <- (1/(Protein_conc*10^(0.58*LogP-2.02)+1)) 
                                                 zval <- round(zval, digits = Digits) 
                                                 return(zval)} 

cfu_mic_Turner.acidic <- function(LogP, Protein_conc, Digits) { zval <- (1/(Protein_conc*10^(0.20*LogP-1.54)+1)) 
                                                  zval <- round(zval, digits = Digits) 
                                                  return(zval)} 

cfu_mic_Turner.neutral <- function(LogP, Protein_conc, Digits) { zval <- (1/(Protein_conc*10^(0.46*LogP-1.51)+1)) 
                                                   zval <- round(zval, digits = Digits) 
                                                   return(zval)} 

#Test to see if the fuctions are working
#cfu_hep_Austin(3.59, 4)
#cfu_hep_Austin.fumic(3.59, 0.25, 4)
#cfu_hep_Kilford(3.59, 4)
#cfu_mic_Austin(3.59, 0.25, 4)
#cfu_mic_Houston(3.59, 0.25, 4)
#cfu_mic_Turner.basic(3.59, 0.25, 4)
#cfu_mic_Turner.acidic(3.59, 0.25, 4)
#cfu_mic_Turner.neutral(3.59, 0.25, 4)

## Prepare lists for fu,inc calculation

List_c_log_P_D <- 1                                             # Initialize the variable for loop
for (x in 1:Num_cpd) {                                        # Make a list to select logD for acidic compounds, and logP for others (i.e. basic, neutral or NA)
  List_c_log_P_D[x] <- ifelse(Col_pKa_type[x] %in% c("Basic", "Neutral", NA), Col_c_logP[x], 
                            Col_c_logD[x]) }

Rawdata <- mutate(Rawdata, LogP_D = List_c_log_P_D)           # Combine the logP/D list into data frame for fu,inc calculation

## Perform fu,inc calculation
List_cfu_hep_Austin <- sapply(List_c_log_P_D, 
                              cfu_hep_Austin, Digits = Val_fu_inc_digit)

List_cfu_hep_Austin_fumic <- sapply(List_c_log_P_D, 
                                    cfu_hep_Austin.fumic, Protein_conc = Val_protein_conc, Digits = Val_fu_inc_digit)

List_cfu_hep_Kilford <- sapply(List_c_log_P_D, 
                               cfu_hep_Kilford, Digits = Val_fu_inc_digit)

List_cfu_mic_Austin <- sapply(List_c_log_P_D, 
                              cfu_mic_Austin, Protein_conc = Val_protein_conc, Digits = Val_fu_inc_digit)

List_cfu_mic_Houston <- sapply(List_c_log_P_D, 
                               cfu_mic_Houston, Protein_conc = Val_protein_conc, Digits = Val_fu_inc_digit)

List_cfu_mic_Turner <- 1                                      # Initialize the variable for loop
for (x in 1:Num_cpd) {                                        # Calculate fu,mic using Turner method; logP only for all pKa classes
  zval <- Col_c_logP[x]
  List_cfu_mic_Turner[x] <- ifelse(Col_pKa_type[x] %in% "Basic", cfu_mic_Turner.basic(zval, Val_protein_conc, Val_fu_inc_digit),
                                  ifelse(Col_pKa_type[x] %in% c("Neutral", NA), cfu_mic_Turner.neutral(zval, Val_protein_conc, Val_fu_inc_digit), 
                                         cfu_mic_Turner.acidic(zval, Val_protein_conc, Val_fu_inc_digit)))}

## Combine fu,mic and fu,hep into data frame
Rawdata <- mutate(Rawdata, cfu_hep_Austin_logP_D = List_cfu_hep_Austin,
                           cfu_hep_Austin_fu_mic = List_cfu_hep_Austin_fumic, 
                           cfu_hep_Kilford_logP_D = List_cfu_hep_Kilford, 
                           cfu_mic_Austin_logP_D = List_cfu_mic_Austin,
                           cfu_mic_Houston_logP_D = List_cfu_mic_Houston,
                           cfu_mic_Turner_logP = List_cfu_mic_Turner)

## Create function to get SF1, SF2 and Qh based on species; need "" for argument (e.g. type "Rat" instead of Rat)
SF1 <- function(Species) { zchr <- paste(Species)
                           zval = grep(zchr, Df_SF_Qh$Species) 
                           zval <- Df_SF_Qh[zval, 2]
                           zval <- as.numeric(zval)
                           return(zval)}

SF2.hep <- function(Species) { zchr <- paste(Species)
                               zval = grep(zchr, Df_SF_Qh$Species) 
                               zval <- Df_SF_Qh[zval, 3]
                               zval <- as.numeric(zval)
                               return(zval)}

SF2.LM <- function(Species) { zchr <- paste(Species)
                              zval = grep(zchr, Df_SF_Qh$Species) 
                              zval <- Df_SF_Qh[zval, 4]
                              zval <- as.numeric(zval)
                              return(zval)}

Qh <- function(Species) { zchr <- paste(Species)
                          zval = grep(zchr, Df_SF_Qh$Species) 
                          zval <- Df_SF_Qh[zval, 5]
                          zval <- as.numeric(zval)
                          return(zval)}

#Test to see if the fuctions are working
#SF1("Rat")
#SF2.hep("Rat")
#SF2.LM("Rat")
#Qh("Rat")
#SF1("Human")
#SF2.hep("Human")
#SF2.LM("Human")
#Qh("Human")

# Calculate scaled CLint
CLint_scaled.hep <- function(CLint, Species, Cell_density, Digits) { zchr <- paste(Species)             
                                                                     zval = CLint*SF1(zchr)*SF2.hep(zchr)/Cell_density*Val_CLh_unit_conv     
                                                                     zval <- round(zval, digits = Digits)
                                                                     return(zval)}

CLint_scaled.LM <- function(CLint, Species, Protein_conc, Digits) { zchr <- paste(Species)             
                                                                    zval = CLint*SF1(zchr)*SF2.LM(zchr)/Protein_conc*Val_CLh_unit_conv     
                                                                    zval <- round(zval, digits = Digits)
                                                                    return(zval)}

#Test to see if the fuctions are working
#CLint_scaled.hep(4.74, "Rat", 1, 2)
#CLint_scaled.hep(4.74, "Rat", 0.5, 2)
#CLint_scaled.hep(4.74, "Human", 0.5, 2)
#CLint_scaled.LM(54, "Rat", 0.25, 2)
#CLint_scaled.LM(54, "Rat", 1, 2)
#CLint_scaled.LM(54, "Human", 1, 2)

List_CLint_s <- 1                                           # Initialize the variable for loop, please note that the unit of CLint,s is L/hr/kg 
for (x in 1:Num_cpd) { zspecies <- Col_species[x]                             
                       List_CLint_s[x] <- ifelse(Col_matrix[x] %in% "Hepatocytes", 
                                                 CLint_scaled.hep(Col_clint[x],zspecies, Val_cell_density, Val_CLint_s_digit), 
                                                 CLint_scaled.LM(Col_clint[x], Val_protein_conc, Val_CLint_s_digit))}

## Create function for different IVIVE methods
IVIVE <- function(CLint_s, Species, Digits) { zchr <- paste(Species)
                                              zval = Qh(zchr)*CLint_s/(Qh(zchr)+CLint_s)
                                              zval <- round(zval, digits = Digits)
                                              return(zval)}

IVIVE.fup <- function(CLint_s, Species, fu_p, Digits) { zchr <- paste(Species)
                                                        zval = Qh(zchr)*CLint_s*fu_p/(Qh(zchr)+CLint_s*fu_p)
                                                        zval <- round(zval, digits = Digits)
                                                        return(zval)}

IVIVE.fup.fuinc <- function(CLint_s, Species, fu_p, fu_inc, Digits) { zchr <- paste(Species)
                                                                      zval = Qh(zchr)*CLint_s*(fu_p/fu_inc)/(Qh(zchr)+CLint_s*(fu_p/fu_inc))
                                                                      zval <- round(zval, digits = Digits)
                                                                      return(zval)}

F1.basic <- function(pKa, pH_p, pH_hep) { Val1 <- 1/(1+10^(pKa - pH_p))
                                          Val2 <- 1/(1+10^(pKa - pH_hep))
                                          zval = Val1/Val2
                                          zval <- round(zval, digits = 1)
                                          return(zval)}

F1.acidic <- function(pKa, pH_p, pH_hep) { Val3 <- 1/(1+10^(pH_p - pKa))
                                           Val4 <- 1/(1+10^(pH_hep - pKa))
                                           zval = Val3/Val4
                                           zval <- round(zval, digits = 1)
                                           return(zval)}

IVIVE.F1 <- function(CLint_s, Species, F1, Digits) { zchr <- paste(Species)
                                                     zval = Qh(zchr)*CLint_s*F1/(Qh(zchr)+CLint_s*F1)
                                                     zval <- round(zval, digits = Digits)
                                                     return(zval)}

IVIVE.fup.F1 <- function(CLint_s, Species, fu_p, F1, Digits) { zchr <- paste(Species)
                                                               zval = Qh(zchr)*CLint_s*fu_p*F1/(Qh(zchr)+CLint_s*fu_p*F1)
                                                               zval <- round(zval, digits = Digits)
                                                               return(zval)}

IVIVE.fup.fuinc.F1 <- function(CLint_s, Species, fu_p, fu_inc, F1, Digits) { zchr <- paste(Species)
                                                                             zval = Qh(zchr)*CLint_s*(fu_p/fu_inc)*F1/(Qh(zchr)+CLint_s*(fu_p/fu_inc)*F1)
                                                                             zval <- round(zval, digits = Digits)
                                                                             return(zval)}
#Test to see if the fuctions are working
# IVIVE(1.54, "Rat", 3)
# IVIVE.fup(1.54, "Rat", 0.327, 2)
# IVIVE.fup.fuinc(1.54, "Rat", 0.327, 0.966, 2)
# F1.basic(9.7, 7.4, 7.0)
# F1.acidic(4, 7.4, 7.0)
# IVIVE.F1(8.514, "Human", 2.5, 2)
# IVIVE.fup.F1(8.514, "Human", 0.059, 2.5, 2)
# IVIVE.fup.fuinc.F1(8.514, "Human", 0.141, 2.5, 2)

# Calculate F1
List_F1_cpKa <- 1                                                            # Initialize the variable for loop 
for (x in 1:Num_cpd) { zpKatype <- Col_pKa_type[x]
                       zpKa <- Col_c_pKa[x]
                       List_F1_cpKa[x] <- ifelse(Col_pKa_type[x] %in% "Basic", F1.basic(zpKa, Val_pH_plasma, Val_pH_heps),
                                            ifelse(Col_pKa_type[x] %in% "Acidic", F1.acidic(zpKa, Val_pH_plasma, Val_pH_heps), 1))}

Rawdata <- mutate(Rawdata, F1 = List_F1_cpKa)

## IVIVE with no correction, fu,p and fu,p + fu,inc
List_CLh <- 1                                                                # Initialize the variable for loop
for (x in 1:Num_cpd) { zspecies <- Col_species[x]                             
                       List_CLh[x] <- IVIVE(List_CLint_s[x], zspecies, Val_CLh_digit)}

List_CLh_cfup <- 1                                                           # Initialize the variable for loop
for (x in 1:Num_cpd) { zspecies <- Col_species[x]                             
                       List_CLh_cfup[x] <- IVIVE.fup(List_CLint_s[x], zspecies, 
                                                     Col_c_fup[x], Val_CLh_digit)}

List_CLh_cfup_cfuhep_Austin <- 1                                             # Initialize the variable for loop
for (x in 1:Num_cpd) { zspecies <- Col_species[x]                             
                       List_CLh_cfup_cfuhep_Austin[x] <- IVIVE.fup.fuinc(List_CLint_s[x], zspecies, 
                                                                         Col_c_fup[x], List_cfu_hep_Austin[x], Val_CLh_digit)}

List_CLh_cfup_cfuhep_Kilford <- 1                                            # Initialize the variable for loop
for (x in 1:Num_cpd) { zspecies <- Col_species[x]                             
                       List_CLh_cfup_cfuhep_Kilford[x] <- IVIVE.fup.fuinc(List_CLint_s[x], zspecies, 
                                                                        Col_c_fup[x], List_cfu_hep_Kilford[x], Val_CLh_digit)}

List_CLh_cfup_cfuhep_Austin_fumic <- 1                                       # Initialize the variable for loop
for (x in 1:Num_cpd) { zspecies <- Col_species[x]                             
                       List_CLh_cfup_cfuhep_Austin_fumic[x] <- IVIVE.fup.fuinc(List_CLint_s[x], zspecies, 
                                                                               Col_c_fup[x], List_cfu_hep_Austin_fumic[x], Val_CLh_digit)}

List_CLh_cfup_cfumic_Austin <- 1                                             # Initialize the variable for loop
for (x in 1:Num_cpd) { zspecies <- Col_species[x]                             
                       List_CLh_cfup_cfumic_Austin[x] <- IVIVE.fup.fuinc(List_CLint_s[x], zspecies, 
                                                                         Col_c_fup[x], List_cfu_mic_Austin[x], Val_CLh_digit)}

List_CLh_cfup_cfumic_Houston <- 1                                            # Initialize the variable for loop
for (x in 1:Num_cpd) { zspecies <- Col_species[x]                             
                       List_CLh_cfup_cfumic_Houston[x] <- IVIVE.fup.fuinc(List_CLint_s[x], zspecies, 
                                                                          Col_c_fup[x], List_cfu_mic_Houston[x], Val_CLh_digit)}

List_CLh_cfup_cfumic_Turner <- 1                                             # Initialize the variable for loop
for (x in 1:Num_cpd) { zspecies <- Col_species[x]                             
                       List_CLh_cfup_cfumic_Turner[x] <- IVIVE.fup.fuinc(List_CLint_s[x], zspecies, 
                                                                         Col_c_fup[x], List_cfu_mic_Turner[x], Val_CLh_digit)}

# IVIVE using F1 or fup + F1 
List_CLh_F1 <- 1                                                             # Initialize the variable for loop
for (x in 1:Num_cpd) { zspecies <- Col_species[x]                             
                       List_CLh_F1[x] <- IVIVE.F1(List_CLint_s[x], zspecies, 
                                                  List_F1_cpKa[x], Val_CLh_digit)}

List_CLh_cfup_F1 <- 1                                                        # Initialize the variable for loop
for (x in 1:Num_cpd) { zspecies <- Col_species[x]                             
                       List_CLh_cfup_F1[x] <- IVIVE.fup.F1(List_CLint_s[x], zspecies, 
                                                           Col_c_fup[x], List_F1_cpKa[x], Val_CLh_digit)}

# Compile all CLh values using different IVIVE methods
Rawdata <- mutate(Rawdata, CLh = List_CLh,
                           CLh_cfup = List_CLh_cfup, 
                           CLh_cfup_cfuhep_Austin = List_CLh_cfup_cfuhep_Austin, 
                           CLh_cfup_cfuhep_Austin_fumic = List_CLh_cfup_cfuhep_Austin_fumic,
                           CLh_cfup_cfuhep_Kilford = List_CLh_cfup_cfuhep_Kilford,
                           CLh_cfup_cfumic_Austin = List_CLh_cfup_cfumic_Austin,
                           CLh_cfup_cfumic_Houston = List_CLh_cfup_cfumic_Houston,
                           CLh_cfup_cfumic_Turner = List_CLh_cfup_cfumic_Turner,
                           CLh_F1 = List_CLh_F1,
                           CLh_cfup_F1 = List_CLh_cfup_F1)

# Calculate % Extraction Ratio (ER)
List_Qh <- sapply(Col_species, Qh)
Rawdata <- mutate(Rawdata, ER = List_CLh/List_Qh*100,
                  Pct_ER_cfup = List_CLh_cfup/List_Qh*100, 
                  Pct_ER_cfup_cfuhep_Austin = List_CLh_cfup_cfuhep_Austin/List_Qh*100, 
                  Pct_ER_cfup_cfuhep_Austin_fumic = List_CLh_cfup_cfuhep_Austin_fumic/List_Qh*100,
                  Pct_ER_cfup_cfuhep_Kilford = List_CLh_cfup_cfuhep_Kilford/List_Qh*100,
                  Pct_ER_cfup_cfumic_Austin = List_CLh_cfup_cfumic_Austin/List_Qh*100,
                  Pct_ER_cfup_cfumic_Houston = List_CLh_cfup_cfumic_Houston/List_Qh*100,
                  Pct_ER_cfup_cfumic_Turner = List_CLh_cfup_cfumic_Turner/List_Qh*100,
                  Pct_ER_F1 = List_CLh_F1/List_Qh*100,
                  Pct_ER_cfup_F1 = List_CLh_cfup_F1/List_Qh*100)

# Calculate AFE: CLh(predicted)/CLh(observed)
Rawdata <- mutate(Rawdata, AFE_no_fu_correction = List_CLh/Col_CLobs,
                  AFE_fup = List_CLh_cfup/Col_CLobs, 
                  AFE_fup_fuhep_Austin = List_CLh_cfup_cfuhep_Austin/Col_CLobs, 
                  AFE_fup_fuhep_Austin_fumic = List_CLh_cfup_cfuhep_Austin_fumic/Col_CLobs,
                  AFE_fup_fuhep_Kilford = List_CLh_cfup_cfuhep_Kilford/Col_CLobs,
                  AFE_fup_fumic_Austin = List_CLh_cfup_cfumic_Austin/Col_CLobs,
                  AFE_fup_fumic_Houston = List_CLh_cfup_cfumic_Houston/Col_CLobs,
                  AFE_fup_fumic_Turner = List_CLh_cfup_cfumic_Turner/Col_CLobs,
                  AFE_F1 = List_CLh_F1/Col_CLobs,
                  AFE_fup_F1 = List_CLh_cfup_F1/Col_CLobs)

# Format values to desired decimal places
View(colnames(Rawdata))
Rawdata[ ,55:74] <- round(Rawdata[ ,55:74], digits = 1)     
Rawdata[ ,c(5:6)] <- round(Rawdata[ ,c(5:6)], digits = 2)     
Rawdata[ ,c(9, 38:43)] <- round(Rawdata[ ,c(9, 38:43)], digits = 3) 
       
IVIVE_CLh <- Rawdata[ ,c(1:6, 8:12, 37:74)]
View(IVIVE_CLh)

Val_current_date <- Sys.Date()                                                # Get the current date to attach to the file name
write_xlsx(IVIVE_CLh,                                                         # Export a short version + current data
           paste(Val_current_date, " IVIVE_CLh", ".xlsx", sep = ""))  
write_xlsx(Rawdata,                                                           # Export full version + current date
           paste(Val_current_date, " IVIVE_CLh_all", ".xlsx", sep = ""))
write_xlsx(Rawdata,                                                           # Export full version (will overwrite the old file)
           "IVIVE_CLh_all.xlsx")


##### CLEAN UP ######

rm(list = ls())   # Clear environment
p_unload(all)     # Remove all add-ons
cat("\014")       # ctrl+L # Clear console