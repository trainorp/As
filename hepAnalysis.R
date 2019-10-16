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

save.image("ArsenicData.RData")
load("ArsenicData.RData")

########### Some data processing ###########
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

xtabs(~is.na(urxuas), data = df1)

xtabs(~is.na(urxuas) + lbdheg, data = df1)
xtabs(~is.na(urxuas) + lbdhem, data = df1)
xtabs(~lbdhem + lbdheg, data = df1)

########### Demographics ###########
# P-value format function
pRound <- function(x){
  if(x > .01){
    x <- as.character(round(x, 2))
  }else{
    if(x < .0001){
      x <- as.character("<0.0001")
    }
    else{
      if(x < .001){
        x <- as.character(round(x, 4))
      }else{
        x <- as.character(round(x, 3))
      }
    }
  }
  return(x)
}

# Survey data:
dfDemo <- df1 %>% select(seqn, sdmvpsu, sdmvstra, ridageyr, riagendr, ridreth3, myWt, urxuas, lbdheg) %>% 
  filter(!is.na(myWt) & !is.na(urxuas) & !is.na(lbdheg))
nhanesDesign0 <- svydesign(id = ~ sdmvpsu, strata = ~ sdmvstra, weights = ~ myWt, nest = TRUE, 
                          data = dfDemo)

# Quantiles for total arsenic:
AsTotalQ <- quantile(dfDemo$urxuas, probs = c(.25, .5, .75))
dfDemo$AsQ <- ifelse(dfDemo$urxuas < AsTotalQ[1], "Q1",
                     ifelse(dfDemo$urxuas < AsTotalQ[2], "Q2", 
                            ifelse(dfDemo$urxuas < AsTotalQ[3], "Q3", "Q4")))

sumStatsFun <- function(var, isNum){
  if(isNum){
    tab0 <- dfDemo %>% group_by(AsQ) %>% dplyr::summarize(mean = mean((!!rlang::sym(var))),
                                                          sd = sd((!!rlang::sym(var))), median = median((!!rlang::sym(var))), 
                                                          Q1 = quantile((!!rlang::sym(var)), probs = .25), Q3 = quantile((!!rlang::sym(var)), probs = .75))
    dec <- 4 - ceiling(log10(max(tab0$mean)))
    # Format nice:
    tab <- rbind(paste0(format(round(tab0$mean, dec), nsmall = dec, trim = TRUE), " pm ", 
                        format(round(tab0$sd, dec), nsmall = dec, trim = TRUE)),
                 paste0(format(round(tab0$median, dec), nsmall = dec, trim = TRUE), " (", 
                        format(round(tab0$Q1, dec), nsmall = dec, trim = TRUE), ", ",
                        format(round(tab0$Q3, dec), nsmall = dec, trim = TRUE), ")"))
    tab <- as.data.frame(tab) 
    names(tab) <- paste0("Q", 1:4)
    tab$levels <- c("Mean pm SD", "Median (Q1, Q3)")
  }else{
    # Summary stats:
    form <- as.formula(paste0("~", var, "+", "AsQ"))
    tab <- as.data.frame(xtabs(form, data = dfDemo)) %>% spread(key = quantVar, value = Freq)
    pTab <- as.data.frame(prop.table(xtabs(form, data = dfDemo), 1)) %>% spread(key = quantVar, value = Freq)
    tab <- tab %>% left_join(pTab, by = var, suffix = c("_Freq", "_Prop"))
    # Format nice:
    tab$Q1 <- paste0(tab$Q1_Freq, " (", format(round(tab$Q1_Prop * 100, 1), nsmall = 1), "%)")
    tab$Q2 <- paste0(tab$Q2_Freq, " (", format(round(tab$Q2_Prop * 100, 1), nsmall = 1), "%)")
    tab$Q3 <- paste0(tab$Q3_Freq, " (", format(round(tab$Q3_Prop * 100, 1), nsmall = 1), "%)")
    tab$Q4 <- paste0(tab$Q4_Freq, " (", format(round(tab$Q4_Prop * 100, 1), nsmall = 1), "%)")
    tab <- tab[, !grepl("_", names(tab))]
    names(tab)[names(tab) == var] <- "levels"
  }
  
  form2 <- as.formula(paste0("log(urxuas)", "~", var))
  form2b <- as.formula(paste0("log(urxuas)", "~1"))
  lmRes <- svyglm(form2, nhanesDesign0)
  lmResb <- svyglm(form2b, nhanesDesign0)
  anova1 <- anova(lmRes, lmResb)
  tab$anova <- anova1$p
  
  tab$var <- var
  tab$totN <- nrow(data)
  return(tab)
}
sumStatsFun("riagendr", FALSE)
sumStatsFun("ridageyr", TRUE)
sumStatsFun("ridreth3", FALSE)


########### Proportion of measurements at LOD ###########
df2 <- as.matrix(df1[,names(df1) %in% c("urxuas", "urxuas3", "urxuas5", "urxuab", "urxuac", "urxudma", "urxumma")])
nonMissingAs <- apply(df2, 2, function(x) sum(!is.na(x)))
isLLOQ <- apply(apply(df2, 2, function(x) x == min(x, na.rm = TRUE)), 2, function(x) sum(x, na.rm = TRUE))
# write.csv(isLLOQ / nonMissingAs, file = "isLLOQ.csv", row.names = FALSE)

df1$urxuas3Above <- factor(ifelse(df1$urdua3lc == 0, 1, 0))
df1$urxuas5Above <- factor(ifelse(df1$urdua5lc == 0, 1, 0))
df1$urxuabAbove <- factor(ifelse(df1$urduablc == 0, 1, 0))
df1$urxuacAbove <- factor(ifelse(df1$urduaclc == 0, 1, 0))
df1$urxudmaAbove <- factor(ifelse(df1$urdudalc == 0, 1, 0))
df1$urxummaAbove <- factor(ifelse(df1$urdummal == 0, 1, 0))

########### HepE IgG logistic models ###########
glmHepEIgG_Total <- glm(lbdheg ~ log(urxuas), data = df1, family = "binomial")
summary(glmHepEIgG_Total)

glmHepEIgM_Total <- glm(lbdhem ~ log(urxuas), data = df1, family = "binomial")
summary(glmHepEIgM_Total)

# Survey design:
nhanesDesign <- svydesign(id = ~ sdmvpsu, strata = ~ sdmvstra, weights = ~ myWt, nest = TRUE, 
                          data = df1 %>% filter(!is.na(wtsa2yr)))

# HepE IgG:
glmEIgG_Total <- svyglm(lbdheg ~ log(urxuas), data = df1 %>% filter(!is.na(wtsa2yr)), 
                 family = "binomial", design = nhanesDesign)
glmEIgG_TotalCoef <- as.data.frame(summary(glmEIgG_Total)$coefficients)
glmEIgG_TotalCoef$var <- "Total"

glmEIgG_ArsenousAcid <- svyglm(lbdheg ~ log(urxuas3), data = df1 %>% filter(!is.na(wtsa2yr)), 
                               family = "binomial", design = nhanesDesign)
glmEIgG_ArsenousAcid2 <- update(glmEIgG_ArsenousAcid, . ~ . + log(urxuas3):urxuas3Above)
glmEIgG_ArsenousAcidCoef <- as.data.frame(summary(glmEIgG_ArsenousAcid)$coefficients)
glmEIgG_ArsenousAcidCoef$var <- "Arsenous Acid"
glmEIgG_ArsenousAcidCoef2 <- as.data.frame(summary(glmEIgG_ArsenousAcid2)$coefficients)
glmEIgG_ArsenousAcidCoef2$var <- "Arsenous Acid2"

glmEIgG_ArsenicAcid <- svyglm(lbdheg ~ log(urxuas5), data = df1 %>% filter(!is.na(wtsa2yr)), 
                            family = "binomial", design = nhanesDesign)
glmEIgG_ArsenicAcid2 <- update(glmEIgG_ArsenicAcid, . ~ . + log(urxuas5):urxuas5Above)
glmEIgG_ArsenicAcidCoef <- as.data.frame(summary(glmEIgG_ArsenicAcid)$coefficients)
glmEIgG_ArsenicAcidCoef$var <- "Arsenic Acid"
glmEIgG_ArsenicAcidCoef2 <- as.data.frame(summary(glmEIgG_ArsenicAcid2)$coefficients)
glmEIgG_ArsenicAcidCoef2$var <- "Arsenic Acid2"

glmEIgG_Arsenobetaine <- svyglm(lbdheg ~ log(urxuab), data = df1 %>% filter(!is.na(wtsa2yr)), 
                                family = "binomial", design = nhanesDesign)
glmEIgG_Arsenobetaine2 <- update(glmEIgG_Arsenobetaine, . ~ . + log(urxuab):urxuabAbove)
glmEIgG_ArsenobetaineCoef <- as.data.frame(summary(glmEIgG_Arsenobetaine)$coefficients)
glmEIgG_ArsenobetaineCoef$var <- "Arsenobetaine"
glmEIgG_ArsenobetaineCoef2 <- as.data.frame(summary(glmEIgG_Arsenobetaine2)$coefficients)
glmEIgG_ArsenobetaineCoef2$var <- "Arsenobetaine2"

glmEIgG_Arsenocholine <- svyglm(lbdheg ~ log(urxuac), data = df1 %>% filter(!is.na(wtsa2yr)), 
                                family = "binomial", design = nhanesDesign)
glmEIgG_Arsenocholine2 <- update(glmEIgG_Arsenocholine, . ~ . + log(urxuac):urxuacAbove)
glmEIgG_ArsenocholineCoef <- as.data.frame(summary(glmEIgG_Arsenocholine)$coefficients)
glmEIgG_ArsenocholineCoef$var <- "Arsenocholine"
glmEIgG_ArsenocholineCoef2 <- as.data.frame(summary(glmEIgG_Arsenocholine2)$coefficients)
glmEIgG_ArsenocholineCoef2$var <- "Arsenocholine2"

glmEIgG_Dimethylarsinic <- svyglm(lbdheg ~ log(urxudma), data = df1 %>% filter(!is.na(wtsa2yr)), 
                                  family = "binomial", design = nhanesDesign)
glmEIgG_Dimethylarsinic2 <- update(glmEIgG_Dimethylarsinic, . ~ . + log(urxudma):urxudmaAbove)
glmEIgG_DimethylarsinicCoef <- as.data.frame(summary(glmEIgG_Dimethylarsinic)$coefficients)
glmEIgG_DimethylarsinicCoef$var <- "Dimethylarsinic"
glmEIgG_DimethylarsinicCoef2 <- as.data.frame(summary(glmEIgG_Dimethylarsinic2)$coefficients)
glmEIgG_DimethylarsinicCoef2$var <- "Dimethylarsinic2"

glmEIgG_Monomethylacrsonic <- svyglm(lbdheg ~ log(urxumma), data = df1 %>% filter(!is.na(wtsa2yr)), 
                                      family = "binomial", design = nhanesDesign)
glmEIgG_Monomethylacrsonic2 <- update(glmEIgG_Monomethylacrsonic, . ~ . + log(urxumma):urxummaAbove)
glmEIgG_MonomethylacrsonicCoef <- as.data.frame(summary(glmEIgG_Monomethylacrsonic)$coefficients)
glmEIgG_MonomethylacrsonicCoef$var <- "Monomethylacrsonic"
glmEIgG_MonomethylacrsonicCoef2 <- as.data.frame(summary(glmEIgG_Monomethylacrsonic2)$coefficients)
glmEIgG_MonomethylacrsonicCoef2$var <- "Monomethylacrsonic2"

# write.csv(rbind(glmEIgG_TotalCoef, glmEIgG_ArsenousAcidCoef, glmEIgG_ArsenousAcidCoef2, 
#                 glmEIgG_ArsenicAcidCoef, glmEIgG_ArsenicAcidCoef2,
#                 glmEIgG_ArsenobetaineCoef, glmEIgG_ArsenobetaineCoef2,
#                 glmEIgG_ArsenocholineCoef, glmEIgG_ArsenocholineCoef2,
#                 glmEIgG_DimethylarsinicCoef, glmEIgG_DimethylarsinicCoef2,
#                 glmEIgG_MonomethylacrsonicCoef, glmEIgG_MonomethylacrsonicCoef2), 
#       file = "individualModels.csv",
#       row.names = FALSE)

########### Medians by group ###########
df3 <- df1 %>% filter(!is.na(urxuas) & !is.na(lbdheg))
xtabs(~lbdheg, data = df3)
round(prop.table(xtabs(~lbdheg, data = df3)),3)

medianAs <- df3 %>% group_by(lbdheg) %>% dplyr::summarize(`Total` = median(urxuas), `Arsenous acid` = median(urxuas3, na.rm = TRUE),
                  `Arsenic acid` = median(urxuas5, na.rm = TRUE), `Arsenobetaine` = median(urxuab, na.rm = TRUE),
                  `Arsenocholine` = median(urxuac, na.rm = TRUE), `Dimethylarsinic acid` = median(urxudma, na.rm = TRUE),
                  `Monomethylacrsonic acid` = median(urxumma, na.rm = TRUE))
medianAs$quant <- "Q2"
# write.csv(medianAs, file = "medianAs.csv")

q3As <- df3 %>% group_by(lbdheg) %>% dplyr::summarize(`Total` = quantile(urxuas, probs = .75), 
                                                      `Arsenous acid` = quantile(urxuas3, probs = .75, na.rm = TRUE),
                            `Arsenic acid` = quantile(urxuas5, probs = .75, na.rm = TRUE), 
                            `Arsenobetaine` = quantile(urxuab, probs = .75, na.rm = TRUE),
                            `Arsenocholine` = quantile(urxuac, probs = .75, na.rm = TRUE), 
                            `Dimethylarsinic acid` = quantile(urxudma, probs = .75, na.rm = TRUE),
                            `Monomethylacrsonic acid` = quantile(urxumma, probs = .75, na.rm = TRUE))
q3As$quant <- "Q3"

p90As <- df3 %>% group_by(lbdheg) %>% dplyr::summarize(`Total` = quantile(urxuas, probs = .9), 
                                                      `Arsenous acid` = quantile(urxuas3, probs = .9, na.rm = TRUE),
                                                      `Arsenic acid` = quantile(urxuas5, probs = .9, na.rm = TRUE), 
                                                      `Arsenobetaine` = quantile(urxuab, probs = .9, na.rm = TRUE),
                                                      `Arsenocholine` = quantile(urxuac, probs = .9, na.rm = TRUE), 
                                                      `Dimethylarsinic acid` = quantile(urxudma, probs = .9, na.rm = TRUE),
                                                      `Monomethylacrsonic acid` = quantile(urxumma, probs = .9, na.rm = TRUE))
p90As$quant <- "90"

write.csv(rbind(medianAs, q3As, p90As), file = "quants.csv", row.names = FALSE)

########### Quantile values ###########
quantile(df1$urxuab[df1$urxuabAbove == 1], probs = c(1/4, 1/2, 3/4), na.rm = TRUE)

quantile(df1$urxuab, probs = seq(0, 1, .10), na.rm = TRUE)
predict(glmEIgG_Arsenobetaine, type = "response",
        newdata = data.frame(urxuabAbove = factor(rep(1, 4), levels = c("0", "1")),
                urxuab = c(0.84, quantile(df1$urxuab[df1$urxuabAbove == 1], probs = c(1/4, 1/2, 3/4), na.rm = TRUE))))

predict(glmEIgG_Arsenobetaine, type = "link",
        newdata = data.frame(urxuabAbove = factor(rep(1, 4), levels = c("0", "1")),
                urxuab = c(0.84, quantile(df1$urxuab[df1$urxuabAbove == 1], probs = c(1/4, 1/2, 3/4), na.rm = TRUE))))

########### Some plots ###########
# png(filename = "TotalBP.png", height = 4, width = 5, units = "in", res = 300)
ggplot(df1 %>% filter(!is.na(lbdheg)), aes(x = lbdheg, y = log(urxuas))) + geom_boxplot() +
  theme_bw() + labs(x = "HepE IgG", y = "Log(Total)")
# dev.off()

# png(filename = "ArsenobetaineBP.png", height = 4, width = 5, units = "in", res = 300)
ggplot(df1 %>% filter(!is.na(lbdheg)), aes(x = lbdheg, y = log(urxuab))) + geom_boxplot() +
  theme_bw() + labs(x = "HepE IgG", y = "Log(Arsenobetaine)")
#dev.off()

########### Variance and correlation analysis ###########
labs <- data.frame(labs = c("Total Urinary As", "Arsenous acid", "Arsenic acid", "Arsenobetaine", 
          "Arsenocholine", "Dimethylarsinic acid", "Monomethylacrsonic acid"), 
          vars = c("urxuas", "urxuas3", "urxuas5", "urxuab", "urxuac", "urxudma", "urxumma"))
apply(df2, 2, function(x) var(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
df2Log <- log(df2)
colnames(df2Log) <- labs$labs
df2cor <- cor(df2Log, use = "pairwise.complete.obs")

# png(filename = "corplot.png", height = 8, width = 8, units = "in", res = 300)
corrplot::corrplot(df2cor, order = "hclust", addCoef.col = "black", diag = FALSE)
# dev.off()

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
