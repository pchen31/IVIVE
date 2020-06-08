pacman::p_load(pacman, dplyr, readxl, writexl, lubridate, stringr) 

# Import Data
Rawdata <- read_excel("CLint_Raw data.xlsx", sheet = "Rawdata")
Df_SF_Qh <- read_excel("CLint_Raw data.xlsx", sheet = "SF_Qh")
Df_Exp_Unit <- read_excel("CLint_Raw data.xlsx", sheet = "Exp_Unit")

# Create Variables
Val_protein_conc <- as.numeric(Df_Exp_Unit$Value[1])                 # Protein concentration used in microsomal stability assay
Val_cell_density <- as.numeric(Df_Exp_Unit$Value[2])                 # Cell density used in hepatocyte stability assay
Val_cpd_conc <- as.numeric(Df_Exp_Unit$Value[3])                     # Substrate concentration used in the assay
Val_fu_digit <- as.numeric(Df_Exp_Unit$Value[4])                     # Number of decimal places for fu values
Val_CLint_s_digit <- as.numeric(Df_Exp_Unit$Value[5])                # Number of decimal places for CLint,s values
Val_CLh_digit <- as.numeric(Df_Exp_Unit$Value[6])                    # Number of decimal places for CLh values
Val_pH_plasma <- as.numeric(Df_Exp_Unit$Value[11])                   # pH of plasma, for F1 calculation
Val_pH_heps <- as.numeric(Df_Exp_Unit$Value[12])                     # pH of hepatoctytes (inside the cells), for F1 calculation
Val_CLint_u_LM_low_CL <- as.numeric(Df_Exp_Unit$Value[13])           # Cutoff of CLint,u as low CL for liver microsomes
Val_CLint_u_LM_high_CL <- as.numeric(Df_Exp_Unit$Value[14])          # Cutoff of CLint,u as high CL for liver microsomes
Val_CLint_u_hep_low_CL <- as.numeric(Df_Exp_Unit$Value[15])          # Cutoff of CLint,u as low CL for hepacoytes
Val_CLint_u_hep_high_CL <- as.numeric(Df_Exp_Unit$Value[16])         # Cutoff of CLint,u as high CL for hepatocytes
Val_ER_low_CL <- as.numeric(Df_Exp_Unit$Value[17])                   # Cutoff of ER as low CL
Val_ER_high_CL <- as.numeric(Df_Exp_Unit$Value[18])                  # Cutoff of ER as low CL
Val_CLh_unit_conv <- ifelse(Df_Exp_Unit$Unit[10] == "L/hr/kg",       # Value to convert CLh to CLint based on desired CLh unit (L/hr/kg or mL/min/kg)
                            60/1000000, 1/1000)

Num_cpd = nrow(Rawdata)                                              # Calculate number of compounds in the Rawdata, for data processing

Col_ID <- Rawdata$ID
Col_species <- Rawdata$Species
Col_clint_lm <- Rawdata$CLint_LM
Col_clint_hep <- Rawdata$CLint_hep
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
Col_exp_date <- Rawdata$Exp_Date

# Get Experiment Dates
Rawdata <- mutate(Rawdata, Exp_Week = week(Rawdata$Exp_Date),
                           Exp_Month_digit = month(Rawdata$Exp_Date),
                           Exp_Month = months(Rawdata$Exp_Date, abbreviate = TRUE),
                           Exp_Year = year(Rawdata$Exp_Date))

Rawdata$Exp_Date <- as_date(Rawdata$Exp_Date)               # Convert exp date back to only yyyy-mm-dd
head(Rawdata)

# Create Functions for fu,inc
cfu_mic_Austin <- function(LogP_D, Protein_conc, Digits) { zval <- (1/(Protein_conc*10^(0.56*LogP_D-1.41)+1)) 
                                                           zval <- round(zval, digits = Digits) 
                                                           return(zval)} 

cfu_mic_Houston <- function(LogP_D, Protein_conc, Digits) { zval <- (1/(Protein_conc*10^(0.072*LogP_D^2+0.067*LogP_D-1.126)+1)) 
                                                            zval <- round(zval, digits = Digits) 
                                                            return(zval)}

cfu_mic_Turner.acidic <- function(LogP, Protein_conc, Digits) { zval <- (1/(Protein_conc*10^(0.20*LogP-1.54)+1)) 
                                                                zval <- round(zval, digits = Digits) 
                                                                return(zval)} 

cfu_mic_Turner.neutral <- function(LogP, Protein_conc, Digits) { zval <- (1/(Protein_conc*10^(0.46*LogP-1.51)+1)) 
                                                                 zval <- round(zval, digits = Digits) 
                                                                 return(zval)} 

cfu_mic_Turner.basic <- function(LogP, Protein_conc, Digits) { zval <- (1/(Protein_conc*10^(0.58*LogP-2.02)+1)) 
                                                               zval <- round(zval, digits = Digits) 
                                                               return(zval)} 

cfu_hep_Austin <- function(LogP_D, Digits) { zval <- (1/(10^(0.4*LogP_D-1.38)+1)) 
                                             zval <- round(zval, digits = Digits)
                                             return(zval)} 

cfu_hep_Austin.fumic <- function(LogP_D, Protein_conc, Digits) { zval <- (1/(Protein_conc*10^(0.56*LogP_D-1.41)+1))
                                                                 zval <- (1/(10^((log10((1-zval)/zval)-0.06)/1.52)+1))
                                                                 zval <- round(zval, digits = Digits)
                                                                 return(zval)} 

cfu_hep_Kilford <- function(LogP_D, Digits) { zval <- (1/(125*0.005*10^(0.072*LogP_D^2+0.067*LogP_D-1.126)+1))      
                                              zval <- round(zval, digits = Digits) 
                                              return(zval)}                 # VR value is 0.005 for 1M cells/mL incubation

# Create a new column: cLogD for acidic, cLogP for basic and neutral compounds 
Rawdata <- Rawdata %>% mutate(c_LogP_D = ifelse(pKa_Type %in% "Acidic", c_LogD, c_LogP)) 
head(Rawdata)

# Calculate fu,inc
Rawdata <- Rawdata %>% mutate(cfu_mic_Austin_logP_D = cfu_mic_Austin(c_LogP_D, Val_protein_conc, Val_fu_digit),
                              cfu_mic_Houston_logP_D = cfu_mic_Houston(c_LogP_D, Val_protein_conc, Val_fu_digit),
                              cfu_mic_Turner_logP = ifelse(pKa_Type == "Acidic", cfu_mic_Turner.acidic(c_LogP, Val_protein_conc, Val_fu_digit),
                                                           ifelse(pKa_Type == "Basic", cfu_mic_Turner.basic(c_LogP, Val_protein_conc, Val_fu_digit),
                                                                  cfu_mic_Turner.neutral(c_LogP, Val_protein_conc, Val_fu_digit))),
                              cfu_hep_Austin_logP_D = cfu_hep_Austin(c_LogP_D, Val_fu_digit),
                              cfu_hep_Austin_fu_mic = cfu_hep_Austin.fumic(c_LogP_D, Val_protein_conc, Val_fu_digit), 
                              cfu_hep_Kilford_logP_D = cfu_hep_Kilford(c_LogP_D, Val_fu_digit) 
                                )

# Calculate CLint,u
Rawdata <- Rawdata %>% mutate(CLint_u_LM_Austin = CLint_LM/cfu_mic_Austin_logP_D,
                              CLint_u_LM_Houston = CLint_LM/cfu_mic_Houston_logP_D,
                              CLint_u_LM_Turner = CLint_LM/cfu_mic_Turner_logP,
                              CLint_u_hep_Austin = CLint_hep/cfu_hep_Austin_logP_D,
                              CLint_u_hep_Austin_fumic = CLint_hep/cfu_hep_Austin_fu_mic,
                              CLint_u_hep_Kilford = CLint_hep/cfu_hep_Kilford_logP_D)

# Create Functions to Get SF & Qh Based on Species, and Extraploate CLint to CLint,s 
SF1 <- function(Species) { zval = which(Df_SF_Qh$Species %in% Species) 
                           zval <- Df_SF_Qh[zval, 2]
                           zval <- as.numeric(unlist(zval))
                           return(zval)}

SF2.hep <- function(Species) { zval = which(Df_SF_Qh$Species %in% Species)
                               zval <- Df_SF_Qh[zval, 3]
                               zval <- as.numeric(unlist(zval))
                               return(zval)}

SF2.LM <- function(Species) { zval = which(Df_SF_Qh$Species %in% Species) 
                              zval <- Df_SF_Qh[zval, 4]
                              zval <- as.numeric(unlist(zval))
                              return(zval)}

Qh <- function(Species) { zval = which(Df_SF_Qh$Species %in% Species) 
                          zval <- Df_SF_Qh[zval, 5]
                          zval <- as.numeric(unlist(zval))
                          return(zval)}

CLint_scaled.LM <- function(CLint, Species, Digits) { zval = CLint*SF1(Species)*SF2.LM(Species)*Val_CLh_unit_conv     
                                                      zval <- round(zval, digits = Digits)
                                                      return(zval)}

CLint_scaled.hep <- function(CLint, Species, Digits) { zval = CLint*SF1(Species)*SF2.hep(Species)*Val_CLh_unit_conv     
                                                       zval <- round(zval, digits = Digits)
                                                       return(zval)}

# Calculate CLint,s
Rawdata <- Rawdata %>% mutate(CLint_s_LM = CLint_scaled.LM(Col_clint_lm, Col_species, Val_CLint_s_digit),
                              CLint_s_hep = CLint_scaled.hep(Col_clint_hep, Col_species, Val_CLint_s_digit)
                              )

# Create Functions for IVIVE
IVIVE <- function(CLint_s, Species, Digits) { zval = Qh(Species)*CLint_s/(Qh(Species)+CLint_s)
                                              zval <- round(zval, digits = Digits)
                                              return(zval)}

IVIVE.fup <- function(CLint_s, Species, fu_p, Digits) { zval = Qh(Species)*CLint_s*fu_p/(Qh(Species)+CLint_s*fu_p)
                                                        zval <- round(zval, digits = Digits)
                                                        return(zval)}

IVIVE.fup.fuinc <- function(CLint_s, Species, fu_p, fu_inc, Digits) { zval = Qh(Species)*CLint_s*(fu_p/fu_inc)/(Qh(Species)+CLint_s*(fu_p/fu_inc))
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

IVIVE.F1 <- function(CLint_s, Species, F1, Digits) { zval = Qh(Species)*CLint_s*F1/(Qh(Species)+CLint_s*F1)
                                                     zval <- round(zval, digits = Digits)
                                                     return(zval)}

IVIVE.fup.F1 <- function(CLint_s, Species, fu_p, F1, Digits) { zval = Qh(Species)*CLint_s*fu_p*F1/(Qh(Species)+CLint_s*fu_p*F1)
                                                               zval <- round(zval, digits = Digits)
                                                               return(zval)}

IVIVE.fup.fuinc.F1 <- function(CLint_s, Species, fu_p, fu_inc, F1, Digits) { zval = Qh(Species)*CLint_s*(fu_p/fu_inc)*F1/(Qh(Species)+CLint_s*(fu_p/fu_inc)*F1)
                                                                             zval <- round(zval, digits = Digits)
                                                                             return(zval)}

CL_class <- function(CL, Cutoff_low_CL, Cutoff_high_CL) { zval = ifelse(CL < Cutoff_low_CL, "Low",
                                                                     ifelse(CL > Cutoff_high_CL, "High", "Medium")
                                                                 )
                                                          return(zval)
                                                        }

# Calculate F1
Rawdata <- Rawdata %>% mutate(F1 = ifelse(pKa_Type == "Basic", F1.basic(c_pKa, Val_pH_plasma, Val_pH_heps),
                                          ifelse(pKa_Type == "Acidic", F1.acidic(c_pKa, Val_pH_plasma, Val_pH_heps), 1))
                              )

# Calculate CLh
Rawdata <- mutate(Rawdata, 
                  CLh_LM = IVIVE(CLint_s_LM, Species, Val_CLh_digit),
                  CLh_LM_cfup = IVIVE.fup(CLint_s_LM, Species, c_fup, Val_CLh_digit),
                  CLh_LM_cfup_cfumic_Austin = IVIVE.fup.fuinc(CLint_s_LM, Species,c_fup, cfu_mic_Austin_logP_D, Val_CLh_digit),
                  CLh_LM_cfup_cfumic_Houston = IVIVE.fup.fuinc(CLint_s_LM, Species,c_fup, cfu_mic_Houston_logP_D, Val_CLh_digit),
                  CLh_LM_cfup_cfumic_Turner = IVIVE.fup.fuinc(CLint_s_LM, Species,c_fup, cfu_mic_Turner_logP, Val_CLh_digit),
                  CLh_hep = IVIVE(CLint_s_hep, Species, Val_CLh_digit),
                  CLh_hep_cfup = IVIVE.fup(CLint_s_hep, Species,c_fup, Val_CLh_digit), 
                  CLh_hep_cfup_cfuhep_Austin = IVIVE.fup.fuinc(CLint_s_hep, Species,c_fup, cfu_hep_Austin_logP_D, Val_CLh_digit), 
                  CLh_hep_cfup_cfuhep_Austin_fumic = IVIVE.fup.fuinc(CLint_s_hep, Species,c_fup, cfu_hep_Austin_fu_mic, Val_CLh_digit),
                  CLh_hep_cfup_cfuhep_Kilford = IVIVE.fup.fuinc(CLint_s_hep, Species, c_fup, cfu_hep_Kilford_logP_D, Val_CLh_digit),
                  CLh_hep_F1 = IVIVE.F1(CLint_s_hep, Species, F1, Val_CLh_digit),
                  CLh_hep_cfup_F1 = IVIVE.fup.F1(CLint_s_hep, Species, c_fup, F1, Val_CLh_digit)
                  )

# Calculate %ER
List_Qh <- sapply(Col_species, Qh)
Rawdata <- mutate(Rawdata, 
                  Pct_ER_LM = CLh_LM/List_Qh*100,
                  Pct_ER_LM_cfup = CLh_LM_cfup/List_Qh*100,
                  Pct_ER_LM_cfup_cfumic_Austin = CLh_LM_cfup_cfumic_Austin/List_Qh*100,
                  Pct_ER_LM_cfup_cfumic_Houston = CLh_LM_cfup_cfumic_Houston/List_Qh*100,
                  Pct_ER_LM_cfup_cfumic_TurnER_LM = CLh_LM_cfup_cfumic_Turner/List_Qh*100,
                  Pct_ER_hep = CLh_hep/List_Qh*100,
                  Pct_ER_hep_cfup = CLh_hep_cfup/List_Qh*100, 
                  Pct_ER_hep_cfup_cfuhep_Austin = CLh_hep_cfup_cfuhep_Austin/List_Qh*100, 
                  Pct_ER_hep_cfup_cfuhep_Austin_fumic = CLh_hep_cfup_cfuhep_Austin_fumic/List_Qh*100,
                  Pct_ER_hep_cfup_cfuhep_Kilford = CLh_hep_cfup_cfuhep_Kilford/List_Qh*100,
                  Pct_ER_hep_F1 = CLh_hep_F1/List_Qh*100,
                  Pct_ER_hep_cfup_F1 = CLh_hep_cfup_F1/List_Qh*100
                  )

# Calculate AFE ( = predicted CL/Observed CL)
Rawdata <- mutate(Rawdata, 
                  AFE_LM = CLh_LM/Col_CLobs,
                  AFE_LM_cfup = CLh_LM_cfup/Col_CLobs,
                  AFE_LM_cfup_cfumic_Austin = CLh_LM_cfup_cfumic_Austin/Col_CLobs,
                  AFE_LM_cfup_cfumic_Houston = CLh_LM_cfup_cfumic_Houston/Col_CLobs,
                  AFE_LM_cfup_cfumic_Turner = CLh_LM_cfup_cfumic_Turner/Col_CLobs,
                  AFE_hep = CLh_hep/Col_CLobs,
                  AFE_hep_cfup = CLh_hep_cfup/Col_CLobs, 
                  AFE_hep_cfup_cfuhep_Austin = CLh_hep_cfup_cfuhep_Austin/Col_CLobs, 
                  AFE_hep_cfup_cfuhep_Austin_fumic = CLh_hep_cfup_cfuhep_Austin_fumic/Col_CLobs,
                  AFE_hep_cfup_cfuhep_Kilford = CLh_hep_cfup_cfuhep_Kilford/Col_CLobs,
                  AFE_hep_F1 = CLh_hep_F1/Col_CLobs,
                  AFE_hep_cfup_F1 = CLh_hep_cfup_F1/Col_CLobs
                  )

# CL Classification
Rawdata <- Rawdata %>% mutate(CL_Class_CLobs = CL_class(CLobs, List_Qh*Val_ER_low_CL, List_Qh*Val_ER_high_CL),
                              CL_Class_CLint_u_cfumic_Austin = CL_class(CLint_u_LM_Austin, Val_CLint_u_LM_low_CL, Val_CLint_u_LM_high_CL),
                              CL_Class_CLint_u_cfumic_Houston = CL_class(CLint_u_LM_Houston, Val_CLint_u_LM_low_CL, Val_CLint_u_LM_high_CL),
                              CL_Class_CLint_u_cfumic_Turner = CL_class(CLint_u_LM_Turner, Val_CLint_u_LM_low_CL, Val_CLint_u_LM_high_CL),
                              CL_Class_Clint_u_cfuhep_Austin = CL_class(CLint_u_hep_Austin, Val_CLint_u_hep_low_CL, Val_CLint_u_hep_high_CL),
                              CL_Class_Clint_u_cfuhep_Austin_fumic = CL_class(CLint_u_hep_Austin_fumic, Val_CLint_u_hep_low_CL, Val_CLint_u_hep_high_CL),
                              CL_Class_Clint_u_cfuhep_Kilford = CL_class(CLint_u_hep_Kilford, Val_CLint_u_hep_low_CL, Val_CLint_u_hep_high_CL),
                              )

# Round values to desired digits
# View(colnames(Rawdata))
Rawdata[ ,65:88] <- round(Rawdata[ ,65:88], digits = 1)     
Rawdata[ ,c(4:6,44:51,53:64)] <- round(Rawdata[ ,c(4:6,44:51,53:64)], digits = 2)     
Rawdata[ ,c(8:11, 38:43)] <- round(Rawdata[ ,c(8:11, 38:43)], digits = 3) 

# Export Processed Data to Excel File Format
Val_current_date <- Sys.Date()                                                # Get the current date to attach to the file name

write_xlsx(Rawdata,                                                           # Export full version + current date
           paste(Val_current_date, " IVIVE_CLh", ".xlsx", sep = ""))

write_xlsx(Rawdata,                                                           # Export full version (will overwrite the old file)
           "IVIVE_CLh.xlsx")

##### CLEAN UP ######

rm(list = ls())   # Clear environment
p_unload(all)     # Remove all add-ons
cat("\014")       # ctrl+L # Clear console