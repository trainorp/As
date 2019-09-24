############ Prereqs ############
## Begin always run
options(stringsAsFactors = FALSE, scipen = 600)
oldPar <- par()

library(tidyverse)
library(Hmisc)
library(survey)

os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive", "~/gdrive")

setwd(paste0(baseDir, "/NHANES/Arsenic/"))

########### Import datasets ###########
years <- c("2011", "2013", "2015")
df1 <- data.frame()
for(year in years){
  demo <- sasxport.get(paste0("Data/DEMO_", year, ".xpt"))
  hepE <- sasxport.get(paste0("Data/HEPE_", year, ".xpt"))
  if(year != "2011"){
    totalAs <- sasxport.get(paste0("Data/UTAS_", year, ".xpt"))
    specAs <- sasxport.get(paste0("Data/UAS_", year, ".xpt"))
    As <- totalAs %>% full_join(specAs %>% select(-wtsa2yr), by = "seqn")
  }else{
    As <- sasxport.get(paste0("Data/UAS_", year, ".xpt"))
  }
  
  hb <- sasxport.get(paste0("Data/GHB_", year, ".xpt"))
  hpv <- sasxport.get(paste0("Data/HPVSWR_", year, ".xpt"))
  
  if(year == "2015"){
    hpv2 <- sasxport.get(paste0("Data/HPVSWC_", year, ".xpt"))
    hpv <- hpv %>% full_join(hpv2, by = "seqn")
  }
  
  # Join all of the datasets together:
  df1year <- demo %>% full_join(hepE, by = "seqn") %>% 
    full_join(As, by = "seqn") %>% full_join(hb, by = "seqn") %>% full_join(hpv, by = "seqn")
  df1year$year <- year
  # labs <- label(df1year)
  
  df1 <-bind_rows(df1, df1year)
  # label(df1, self = FALSE) <- labs
}
names(df1)[names(df1) == "wtsa2yr.x"] <- "wtsa2yr"

# Adjust the weights:
df1$myWt <- 1/3 * df1$wtsa2yr

# IgG HepE
df1$lbdheg <- factor(df1$lbdheg)
levels(df1$lbdheg) <- c("Pos HepE IgG", "Neg HepE IgG")
df1$lbdheg <- factor(df1$lbdheg, levels = c("Neg HepE IgG", "Pos HepE IgG"))

# IgM HepE
df1$lbdhem <- factor(df1$lbdhem)
levels(df1$lbdhem) <- c("Pos HepE IgM", "Neg HepE IgM")
df1$lbdhem <- factor(df1$lbdhem, levels = c("Neg HepE IgM", "Pos HepE IgM"))

# High risk HPV:
df1$lbxhp2c <- factor(df1$lbxhp2c)
levels(df1$lbxhp2c) <- c("Positive", "Negative", "Inadequate")
df1$lbxhp2c <- factor(df1$lbxhp2c, levels = c("Negative", "Positive", "Inadequate"))
df1$lbxhp2c2 <- df1$lbxhp2c
df1$lbxhp2c2[!is.na(df1$lbxhp2c) & df1$lbxhp2c == "Inadequate"] <- NA
df1$lbxhp2c2 <- factor(df1$lbxhp2c2)

xtabs(~is.na(urxuas) + lbdheg, data = df1)
xtabs(~is.na(urxuas) + lbdhem, data = df1)
xtabs(~lbdhem + lbdheg, data = df1)

ggplot(df1, aes(x = lbdheg, y = log(urxuas))) + geom_boxplot() +
  theme_bw()

glmHepEIgG_Total <- glm(lbdheg ~ log(urxuas), data = df1, family = "binomial")
summary(glmHepEIgG_Total)

glmHepEIgM_Total <- glm(lbdhem ~ log(urxuas), data = df1, family = "binomial")
summary(glmHepEIgM_Total)

nhanesDesign <- svydesign(id = ~ sdmvpsu, strata = ~ sdmvstra, weights = ~ myWt, nest = TRUE, 
                          data = df1 %>% filter(!is.na(wtsa2yr)))

# HepE IgG:
glmEIgG_Total <- svyglm(lbdheg ~ log(urxuas), data = df1 %>% filter(!is.na(wtsa2yr)), 
                 family = "binomial", design = nhanesDesign)
glmEIgG_TotalCoef <- as.data.frame(summary(glmEIgG_Total)$coefficients)
glmEIgG_TotalCoef$var <- "Total"

glmEIgG_ArsenousAcid <- svyglm(lbdheg ~ log(urxuas3), data = df1 %>% filter(!is.na(wtsa2yr)), 
                               family = "binomial", design = nhanesDesign)
glmEIgG_ArsenousAcidCoef <- as.data.frame(summary(glmEIgG_ArsenousAcid)$coefficients)
glmEIgG_ArsenousAcidCoef$var <- "Arsenous Acid"

glmEIgG_ArsenicAcid <- svyglm(lbdheg ~ log(urxuas5), data = df1 %>% filter(!is.na(wtsa2yr)), 
                            family = "binomial", design = nhanesDesign)
glmEIgG_ArsenicAcidCoef <- as.data.frame(summary(glmEIgG_ArsenicAcid)$coefficients)
glmEIgG_ArsenicAcidCoef$var <- "Arsenic Acid"

glmEIgG_Arsenobetaine <- svyglm(lbdheg ~ log(urxuab), data = df1 %>% filter(!is.na(wtsa2yr)), 
                                family = "binomial", design = nhanesDesign)
glmEIgG_ArsenobetaineCoef <- as.data.frame(summary(glmEIgG_Arsenobetaine)$coefficients)
glmEIgG_ArsenobetaineCoef$var <- "Arsenobetaine"

glmEIgG_Arsenocholine <- svyglm(lbdheg ~ log(urxuac), data = df1 %>% filter(!is.na(wtsa2yr)), 
                                family = "binomial", design = nhanesDesign)
glmEIgG_ArsenocholineCoef <- as.data.frame(summary(glmEIgG_Arsenocholine)$coefficients)
glmEIgG_ArsenocholineCoef$var <- "Arsenocholine"

glmEIgG_Dimethylarsinic <- svyglm(lbdheg ~ log(urxudma), data = df1 %>% filter(!is.na(wtsa2yr)), 
                                  family = "binomial", design = nhanesDesign)
glmEIgG_DimethylarsinicCoef <- as.data.frame(summary(glmEIgG_Dimethylarsinic)$coefficients)
glmEIgG_DimethylarsinicCoef$var <- "Dimethylarsinic"

glmEIgG_Monomethylacrsonic <- svyglm(lbdheg ~ log(urxumma), data = df1 %>% filter(!is.na(wtsa2yr)), 
                                      family = "binomial", design = nhanesDesign)
glmEIgG_MonomethylacrsonicCoef <- as.data.frame(summary(glmEIgG_Monomethylacrsonic)$coefficients)
glmEIgG_MonomethylacrsonicCoef$var <- "Monomethylacrsonic"

write.csv(rbind(glmEIgG_TotalCoef, glmEIgG_ArsenousAcidCoef, glmEIgG_ArsenicAcidCoef, glmEIgG_ArsenobetaineCoef, 
      glmEIgG_ArsenocholineCoef, glmEIgG_DimethylarsinicCoef, glmEIgG_MonomethylacrsonicCoef), 
      file = "individualModels.csv",
      row.names = FALSE)

# Hemoglobin A1c:
lmGhb <- lm(log(lbxgh) ~ log(urxuas), data = df1 %>% filter(!is.na(wtsa2yr)))
summary(lmGhb)
svlm1 <- svyglm(log(lbxgh) ~ log(urxuas), data = df1 %>% filter(!is.na(wtsa2yr)), 
                family = "gaussian", design = nhanesDesign)
summary(svlm1)

# HR HPV:
table(df1$lbxhp2c2[!is.na(df1$wtsa2yr)])
glmHRHPV_Total <- glm(lbxhp2c2 ~ log(urxuas), data = df1, family = "binomial")
summary(glmHRHPV_Total)

glmHRHPV_ArsenousAcid <- glm(lbxhp2c2 ~ log(urxuas3), data = df1, family = "binomial")
summary(glmHRHPV_ArsenousAcid)

glmHRHPV_ArsenicAcid <- glm(lbxhp2c2 ~ log(urxuas5), data = df1, family = "binomial")
summary(glmHRHPV_ArsenicAcid)

glmHRHPV_Arsenobetaine <- glm(lbxhp2c2 ~ log(urxuab), data = df1, family = "binomial")
summary(glmHRHPV_Arsenobetaine)

glmHRHPV_Arsenocholine <- glm(lbxhp2c2 ~ log(urxuac), data = df1, family = "binomial")
summary(glmHRHPV_Arsenocholine)

glmHRHPV_Dimethylarsinic <- glm(lbxhp2c2 ~ log(urxudma), data = df1, family = "binomial")
summary(glmHRHPV_Dimethylarsinic)

glmHRHPV_Monomethylacrsonic <- glm(lbxhp2c2 ~ log(urxumma), data = df1, family = "binomial")
summary(glmHRHPV_Monomethylacrsonic)

svglmHRHPV <- svyglm(lbxhp2c2 ~ log(urxuas), data = df1 %>% filter(!is.na(wtsa2yr)),
                     family = "binomial", design = nhanesDesign)
summary(svglmHRHPV)
