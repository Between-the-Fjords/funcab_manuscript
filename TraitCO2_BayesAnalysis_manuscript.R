# Required analysis packages
library(tidyverse)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(corrplot)
library(Hmisc)
library(ggfortify)
library(devtools)
library(factoextra)
library(nimble)
library(DHARMa)
library(coda)
library(cowplot)
library(ggridges)
library(eulerr)
library(parallel)  # Adds the library that allows to detect the number of cores you have
library(future)    # Adds the library that allows asynchronous calling of R code
plan(multisession) # Allows the future library to do parallel processing

# Define the location of the workspace where the script is to import data files and
# produce output files
curWorkspace <- tempdir()  # Change this to the relevant path on your system

CO2_mass_traits<-readRDS(file= file.path(curWorkspace, "CO2_traits06122021.rds"))%>%
  mutate(GPP = -1*GPP)%>% # turn GPP positive for modelling 
  filter(PAR >200)

# Plot predictors vs GPP

#colors for temperature gradient
#display.brewer.pal(n = 3, name = 'Reds')
#brewer.pal(n = 3, name = "BuGn")
#   "#FEE8C8", "#FDBB84", "#E34A33"

#GPP response to PAR
ggplot(CO2_mass_traits, aes(x=PAR, y= GPP, fill=as.factor(T_level), shape=as.factor(P_level)))+
         geom_point()+
         scale_fill_manual(values= c("#FEE0D2", "#FC9272", "#DE2D26"), name = "Temperature", labels = c("Alpine", "Sub-alpine", "Boreal"))+ 
         guides(fill=guide_legend(override.aes=list(shape=21)))+
         scale_shape_manual(values= c(24, 22, 21, 25), name = "Precipitation", labels = c("600mm", "1200mm", "2000mm", "2700mm"))+
         theme_classic()
         
GPP_predictors<- CO2_mass_traits%>%
  select(GPP, T_level, P_level, T_summer, P_annual, Richness, Evenness, Diversity, VegetationHeight, 
         CWM_N, CWM_C, CWM_CN, CWM_LDMC, CWM_LT, CWM_LA, CWM_SLA, CWM_VH, F.Dispersion)

library(plotly)
plot_ly(GPP_predictors, x= ~VegetationHeight, y= ~Richness, z= ~GPP, 
        marker = list(color = ~GPP, showscale = TRUE))

#display.brewer.pal(n = 3, name = 'Reds')
#brewer.pal(n = 3, name = "BuGn")
#   "#FEE8C8", "#FDBB84", "#E34A33"

GPP_predictors = melt(GPP_predictors, id.vars= c("GPP", "T_level", "P_level"))
ggplot(GPP_predictors) +
  geom_point(aes(value, GPP, fill=as.factor(T_level), shape=as.factor(P_level)))+ 
  #geom_smooth(aes(value, GPP), method=lm, se=FALSE, col="black") +
  scale_fill_manual(values= c("#FEE0D2", "#FC9272", "#DE2D26"), name = "Temperature", labels = c("Alpine", "Sub-alpine", "Boreal"))+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  scale_shape_manual(values= c(24, 22, 21, 25), name = "Precipitation", labels = c("600mm", "1200mm", "2000mm", "2700mm"))+
  facet_wrap(~variable, scales="free_x") +
  labs(x = "Predictor", y = "GPP")+
  theme_classic()


Reco_predictors<- CO2_mass_traits%>%
  select(Reco, T_level, P_level, T_summer, P_annual, Richness, Evenness, Diversity, VegetationHeight, 
         CWM_N, CWM_C, CWM_CN, CWM_LDMC, CWM_LT, CWM_LA, CWM_SLA, CWM_VH, F.Dispersion)

Reco_predictors = melt(Reco_predictors, id.vars= c("Reco", "T_level", "P_level"))
ggplot(Reco_predictors) +
  geom_point(aes(value, Reco, fill=as.factor(T_level), shape=as.factor(P_level)))+ 
  #geom_smooth(aes(value, Reco), method=lm, se=FALSE, col="black") +
  scale_fill_manual(values= c("#FEE0D2", "#FC9272", "#DE2D26"), name = "Temperature", labels = c("Alpine", "Sub-alpine", "Boreal"))+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  scale_shape_manual(values= c(24, 22, 21, 25), name = "Precipitation", labels = c("600mm", "1200mm", "2000mm", "2700mm"))+
  facet_wrap(~variable, scales="free_x") +
  labs(x = "Predictor", y = "Reco")+
  theme_classic()

NPP_predictors<- CO2_mass_traits%>%
  select(total.biomass, P_level, T_level, VegetationHeight) #Richness, Evenness, Diversity, 

lmNPP<- lm(total.biomass~VegetationHeight, data = NPP_predictors)
summary(lmNPP)

ggplot(NPP_predictors, aes(VegetationHeight, total.biomass, shape=factor(P_level)))+
  geom_point(aes(fill=factor(T_level)), size=4)+
  #geom_smooth(aes(value, total.biomass), method=lm, se=FALSE, col="black") +
  scale_fill_manual(values= c("#FEE0D2", "#FC9272", "#DE2D26"), name = "Temperature", labels = c("Alpine", "Sub-alpine", "Boreal"))+ 
  guides(fill=guide_legend(override.aes=list(shape=21, size= 4)))+
  scale_shape_manual(values= c(24, 22, 21, 25), name = "Precipitation", labels = c("600mm", "1200mm", "2000mm", "2700mm"))+
  labs(y = "NPP (g/m2)")+
  theme_classic()


NPP_predictors = melt(NPP_predictors, id.vars= c("VegetationHeight", "T_level", "P_level"))
ggplot(NPP_predictors) +
  geom_point(aes(value,VegetationHeight, fill=as.factor(T_level), shape=as.factor(P_level), shape =1)) + 
  #geom_smooth(aes(value, VegetationHeight), method=lm, se=FALSE, col="black") +
  scale_fill_manual(values= c("#FEE0D2", "#FC9272", "#DE2D26"), name = "Temperature", labels = c("Alpine", "Sub-alpine", "Boreal"))+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  scale_shape_manual(values= c(24, 22, 21, 25), name = "Precipitation", labels = c("600mm", "1200mm", "2000mm", "2700mm"))+
  facet_wrap(~variable, scales="free_x") +
  labs(x = "Predictor", y = "NPP")+
  theme_classic()

###### correlation matrix of predictors
###### Calculate median per plot per year to not have duplicate predictors for correlation matrix
CO2_mass_traits_median <- CO2_mass_traits%>%
  dplyr::rename(Site =site)%>%
  dplyr::group_by(turfID, Site, year, P_level, T_level, treatment) %>%
  #summarise(no_rows = length(turfID))%>%
  dplyr::summarise_all(median, na.rm=TRUE)%>%
  dplyr::rename(Year = year)%>%
  ungroup()

apply(CO2_mass_traits_median, 2, function(x) length(which(!is.na(x))))

C_VT_temp<- CO2_mass_traits_median%>%
  select(T_summer, T_level, Richness, Evenness, Diversity, VegetationHeight, 
         CWM_N, CWM_C, CWM_CN, CWM_LDMC, CWM_LT, CWM_LA, CWM_SLA, CWM_VH, F.Dispersion)

C_VT_temp = melt(C_VT_temp, id.vars= c("T_summer", "T_level"))
ggplot(C_VT_temp) +
  geom_jitter(aes(T_summer, value, col= T_level), shape =1) + 
  geom_smooth(aes(T_summer, value), method= "lm", se=FALSE, col= "black") +
  facet_wrap(~variable, scales="free") +
  labs(x = "Summer Temperature", y = "response")+
  theme_classic()

C_VT_prec<- CO2_mass_traits_median%>%
  select(P_annual, P_level, Richness, Evenness, Diversity, VegetationHeight, 
         CWM_N, CWM_C, CWM_CN, CWM_LDMC, CWM_LT, CWM_LA, CWM_SLA, CWM_VH, F.Dispersion)

C_VT_prec = melt(C_VT_prec, id.vars= c("P_annual", "P_level"))
ggplot(C_VT_prec) +
  geom_jitter(aes(P_annual, value, col= P_level), shape =1) + 
  geom_smooth(aes(P_annual, value, col = P_level), method= "lm", se=FALSE, col= "black") +
  facet_wrap(~variable, scales="free") +
  labs(x = "annual Precipitation", y = "response")+
  theme_classic()

# create correlation matrix for predictors Climate, Veg structure and traits 
cor_CVT <- CO2_mass_traits_median %>%
  ungroup()%>%
  select(T_summer, P_annual, Richness, Evenness, Diversity, VegetationHeight,
         CWM_N, CWM_C, CWM_CN, CWM_LDMC, CWM_LT, CWM_LA, CWM_SLA, CWM_VH, F.Dispersion)

#T_summer, P_annual, Richness, Evenness, Diversity, VegetationHeight, CWM_N, CWM_C, CWM_CN, CWM_LDMC, CWM_Lth, CWM_LA, CWM_SLA, CWM_VH, FDiv_Traits, FEve_Traits, FRic_Traits, FDis_Traits, RaoQ_Traits, FDis_C, FDis_N, FDis_CN, FDis_SLA, FDis_LA,  FDis_LDMC, FDis_Lth, FDis_VH

cor_CVT <- rcorr(as.matrix(cor_CVT))

cor_CVT$r #the correlation matrix
cor_CVT$n #the matrix of the number of observations used in analyzing each pair of variables
cor_CVT$P #the p-values corresponding to the significance levels of correlations.

par(mfrow=c(1,1))
#brewer.pal(11, "PiYG")
#brewer.pal(11, "PuOr")
#col1 <- colorRampPalette(c("#40004B", "#762A83", "#9970AB", "#C2A5CF", "#E7D4E8", "#F7F7F7", "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#00441B"))
#col2 <- colorRampPalette(c("#7F3B08", "#B35806", "#E08214", "#FDB863", "#FEE0B6", "#F7F7F7", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788", "#2D004B"))
col3 <- colorRampPalette(c("#8E0152", "#C51B7D", "#DE77AE", "#F1B6DA", "#FDE0EF", "#F7F7F7", "#E6F5D0", "#B8E186", "#7FBC41", "#4D9221", "#276419"))

corrplot.mixed(cor_CVT$r,  lower.col = "black", upper.col= col3(10), number.cex = .7, upper = "square", diag = "u", 
               tl.pos = "lt", tl.col = "black", p.mat = cor_CVT$P, sig.level = 0.05,  insig = "blank")


###########################################################################################################################
#PCA

Traitdata<- CO2_mass_traits_median%>%
  select(Richness, Evenness, Diversity, "Veg Height"= VegetationHeight, N = CWM_N , C =CWM_C, CN = CWM_CN, LDMC = CWM_LDMC, LT = CWM_LT, LA = CWM_LA, SLA = CWM_SLA, VH = CWM_VH, FDispersion = F.Dispersion)
TraitPCA <- princomp(Traitdata, cor= TRUE, scores=TRUE) #, Temperature = T_summer_longterm, Precipitation = P_annual_longterm

PCAplot<- autoplot(TraitPCA, data = CO2_mass_traits_median, fill= "T_level", shape = "P_level", size = 3,
         loadings = TRUE, loadings.colour = 'black', loadings.label.colour = "black", 
         loadings.label = TRUE, loadings.label.size = 5, loadings.label.vjust = -.6, loadings.label.hjust = 0.9)+
  scale_fill_manual(values= c("#FEE0D2", "#FC9272", "#DE2D26"), name = "Temperature", labels = c("Alpine", "Sub-alpine", "Boreal"))+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  scale_shape_manual(values= c(24, 22, 21, 25), name = "Precipitation", labels = c("600mm", "1200mm", "2000mm", "2700mm"))+
  theme_classic()+
  theme(legend.position = "right", title=element_text(size=14), axis.text = element_text(size=13), legend.text=element_text(size=14), legend.title = element_text(size=14))
PCAplot

#purple color pallete for climate #efedf5, #bcbddc, #756bb1
# red color pallete "#FEE0D2", "#FC9272", "#DE2D26"
  
Climdata<- CO2_mass_traits_median%>%
  select(T_summer, P_annual)
ClimPCA<- princomp(Climdata, cor= TRUE, scores=TRUE)
autoplot(ClimPCA, data = CO2_mass_traits_median, fill= "T_level", col = "T_level", shape = "P_level",
         loadings = TRUE, loadings.colour = 'green', 
         loadings.label = TRUE, loadings.label.size = 4, loadings.label.vjust = -1, loadings.label.hjust = 1)

unclass(TraitPCA$loadings)
screeplot(TraitPCA)

TraitPCA.scores<- TraitPCA$scores
TraitPCA.scores<- cbind(TraitPCA.scores, CO2_mass_traits_median$turfID)
TraitPCA.scores<- cbind(TraitPCA.scores, as.character(CO2_mass_traits_median$Year))
TraitPCA.scores<- as.data.frame(TraitPCA.scores)%>%
  rename(turfID = V10, year = V11)%>%
  mutate(Comp.1 = as.numeric(as.character(Comp.1)),
         Comp.2 = as.numeric(as.character(Comp.2)))


#### Climate effects on Vegetation properties ######
########################################## Separate models for indirect effects to not duplicate plot level data !! ########################################


# Import glmNIMBLE and gppLightCurveCorrection functions and their dependencies
# from the PaGAn archive folder
source("https://raw.githubusercontent.com/joechip90/PaGAn/master/inst/archive_R/mcmcInternals.R")
source("https://raw.githubusercontent.com/joechip90/PaGAn/master/inst/archive_R/glmNIMBLE.R")
source("https://raw.githubusercontent.com/joechip90/PaGAn/master/inst/archive_R/gppLightCurveCorrection.R")


vegHeightModel<- glmNIMBLE(modelFormula = VegetationHeight ~ T_summer * P_annual, 
                           errorFamily = Gamma(link = "log"), regCoeffs = "lasso", 
                           mcmcParams = list(numRuns = 100000, thinDensity = 40, predictThinDensity = 40),
                           inputData =CO2_mass_traits_median)
#plot(vegHeightModel$DHARMaResiduals)
#vegHeightModel$parameterFigure

RichnessModel<- glmNIMBLE(modelFormula = Richness ~ T_summer * P_annual, 
                          errorFamily = list(family = "negbinomial", link = "log") , regCoeffs = "lasso", 
                          mcmcParams = list(numRuns = 100000, thinDensity = 40, predictThinDensity = 40),
                          inputData = CO2_mass_traits_median)
plot(RichnessModel$DHARMaResiduals)
RichnessModel$parameterFigure
#!negative binomial gives inverse parameters

EvennessModel<- glmNIMBLE(modelFormula = Evenness ~ T_summer * P_annual, 
                          errorFamily = list(family = "beta", link = "logit"), regCoeffs = "lasso", 
                          inputData =CO2_mass_traits_median)
#plot(EvennessModel$DHARMaResiduals)
#EvennessModel$parameterFigure

DiversityModel<- glmNIMBLE(modelFormula = Diversity ~ T_summer * P_annual, 
                           errorFamily = Gamma(link = "log"), regCoeffs = "lasso", 
                           mcmcParams = list(numRuns = 100000, thinDensity = 40, predictThinDensity = 40),
                           inputData =CO2_mass_traits_median)
#plot(DiversityModel$DHARMaResiduals)
#DiversityModel$parameterFigure

VHModel<- glmNIMBLE(modelFormula = CWM_VH ~ T_summer * P_annual, 
                    errorFamily = Gamma(link = "log"), regCoeffs = "lasso", 
                    inputData =CO2_mass_traits_median)
#plot(VHModel$DHARMaResiduals)
#VHModel$parameterFigure

LAModel<- glmNIMBLE(modelFormula = CWM_LA ~ T_summer * P_annual, 
                    errorFamily = gaussian(), regCoeffs = "lasso", 
                    inputData =CO2_mass_traits_median)
#plot(LAModel$DHARMaResiduals)
#LAModel$parameterFigure

SLAModel<- glmNIMBLE(modelFormula = CWM_SLA ~ T_summer * P_annual, 
                     errorFamily = Gamma(link = "log"), regCoeffs = "lasso", 
                     inputData =CO2_mass_traits_median)
#plot(SLAModel$DHARMaResiduals)
#SLAModel$parameterFigure

LTModel<- glmNIMBLE(modelFormula = CWM_LT ~ T_summer * P_annual, 
                    errorFamily = Gamma(link = "log"), regCoeffs = "lasso", 
                    inputData =CO2_mass_traits_median)
#plot(LTModel$DHARMaResiduals)
#LTModel$parameterFigure

LDMCModel<- glmNIMBLE(modelFormula = CWM_LDMC ~ T_summer * P_annual, 
                      errorFamily = Gamma(link = "log"), regCoeffs = "lasso", 
                      inputData =CO2_mass_traits_median)
#plot(LDMCModel$DHARMaResiduals)
#LDMCModel$parameterFigure

NModel<- glmNIMBLE(modelFormula = CWM_N/100 ~ T_summer * P_annual, 
                   errorFamily = list(family = "beta", link = "logit"), regCoeffs = "lasso", 
                   inputData =CO2_mass_traits_median)
#plot(NModel$DHARMaResiduals)
#NModel$parameterFigure

CModel<- glmNIMBLE(modelFormula = CWM_C/100 ~ T_summer * P_annual, 
                   errorFamily = list(family = "beta", link = "logit"), regCoeffs = "lasso", 
                   inputData =CO2_mass_traits_median)
#errorFamily = Gamma(link = "log")
#plot(CModel$DHARMaResiduals)
#CModel$parameterFigure

CNModel<- glmNIMBLE(modelFormula = CWM_CN ~ T_summer * P_annual, 
                    errorFamily = Gamma(link = "log"), regCoeffs = "lasso", 
                    inputData =CO2_mass_traits_median)
#plot(CNModel$DHARMaResiduals)
#CNModel$parameterFigure

FDisModel<- glmNIMBLE(modelFormula = F.Dispersion ~ T_summer * P_annual, 
                      errorFamily = Gamma(link = "log"), regCoeffs = "lasso", 
                      inputData =CO2_mass_traits_median)
#plot(FDisModel$DHARMaResiduals)
#FDisModel$parameterFigure

# create datatable with all output of individual models
C.indirect<-rbind(vegHeightModel$parameterSummary, RichnessModel$parameterSummary, EvennessModel$parameterSummary, DiversityModel$parameterSummary, VHModel$parameterSummary, LAModel$parameterSummary, SLAModel$parameterSummary, LTModel$parameterSummary, LDMCModel$parameterSummary, CModel$parameterSummary, NModel$parameterSummary, CNModel$parameterSummary, FDisModel$parameterSummary)
parameter<- rownames(C.indirect)
C.indirect<- as.data.frame(cbind(C.indirect, parameter))%>%
  mutate(Mean = as.numeric(Mean),
         Median = as.numeric(Median),
         St.Dev. = as.numeric(St.Dev.),
         CIlow_95 = as.numeric(`95%CI_low`),
         CIupp_95 = as.numeric(`95%CI_upp`),
         significance= ifelse(CIlow_95 < 0 & CIupp_95 < 0, "significant", ifelse(CIlow_95 > 0 & CIupp_95 > 0, "significant", "non-significant")))%>%
  mutate_if(is.numeric, round, 3)%>%
  select(-`95%CI_low`, -`95%CI_upp`)

write.csv2(C.indirect, file.path(curWorkspace, "Ceffect_VT.csv"))

T_effect<- C.indirect%>%
  filter(parameter == "T_summerCoeff")%>%
  mutate(Climate = "Temperature",
         variable = c("VegetationHeight", "Richness", "Eveness", "Diversity", "VH", "LA", "SLA", "LT", "LDMC", "C", "N", "CN", "FDis"))%>%
  select(-parameter)

P_effect<- C.indirect%>%
  filter(parameter == "P_annualCoeff")%>%
  mutate(Climate = "Precipitation",
         variable = c("VegetationHeight", "Richness", "Eveness", "Diversity", "VH", "LA", "SLA", "LT", "LDMC", "C", "N", "CN", "FDis"))%>%
  select(-parameter)

TP_interaction<- C.indirect%>%
  filter(parameter == "T_summer_P_annualCoeff")%>%
  mutate(Climate = "Interaction",
         variable = c("VegetationHeight", "Richness", "Eveness", "Diversity", "VH", "LA", "SLA", "LT", "LDMC", "C", "N", "CN", "FDis"))%>%
  select(-parameter)

# correct inverse parameter values from negative binomial distribution for Richness 
Climate_effect<- rbind(T_effect, P_effect, TP_interaction)%>%
  mutate(Mean = ifelse(variable == "Richness", Mean*-1, Mean),
         Median = ifelse(variable == "Richness", Median*-1, Median),
         CIlow_95 = ifelse(variable == "Richness", CIlow_95*-1, CIlow_95),
         CIupp_95 = ifelse(variable == "Richness", CIupp_95*-1, CIupp_95)
  )
         
rects <- data.frame(xstart = seq(0.5,11.5,1), xend = seq(1.5,12.5,1), 
                    col = rep(c("a", "b")))


Ceffect_VT<-ggplot(Climate_effect)+
  geom_point(aes(x= variable, y = Mean, col= Climate, shape= significance), size=2.5, position=position_dodge(width = .5))+
  geom_errorbar(aes(x= variable, ymin = CIlow_95, ymax= CIupp_95, col= Climate), size= 1, width= 0, position=position_dodge(width = .5))+ 
  scale_color_manual(name = "Climate", labels= c("Temp:Prec", "Precipitation", "Temperature"), values = c("#2e4057", "#00798c", "#d1495b"))+
  scale_shape_manual(values=c(1,19), guide= "none")+
  geom_rect(data=rects, aes(ymin=-0.2, ymax=Inf, xmin=xstart, xmax=xend, fill=col), alpha =0.1)+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = c("white", "grey30"), guide ="none") +
  ylim(-0.2, 0.4)+
  coord_flip()+
  scale_y_continuous(limits = c(-0.2, 0.35), expand = c(0, 0))+
  scale_x_discrete(limits=c("FDis", "CN", "N", "C", "LDMC", "LT", "SLA", "LA", "VH", 
                            "VegetationHeight", "Diversity", "Eveness", "Richness"),
                   labels=c("F.Dispersion", "CN", "N", "C", "LDMC", "LT", "SLA", "LA", "VH", 
                            "Veg Height", "Diversity", "Evenness", "Richness"))+
  theme_classic()+
  theme(legend.position = c(0.8,0.1), strip.background = element_blank(), 
        strip.text.x= element_blank(), axis.title.x=element_blank(), 
        axis.text.x=element_text(size = 14), axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14), axis.title.y=element_blank(),
        title=element_text(size=14), legend.text=element_text(size=14))
Ceffect_VT

CVT_PCA <- plot_grid(Ceffect_VT, PCAplot, labels = c("A", "B"),   hjust = -1,   nrow = 1 )
CVT_PCA



####################################################################################################################
###### Direct vs Indirect effects of Climate on GPP and Reco

CVToutValues <- gppLightCurveCorrection(
  # Set the input data
  inputData = CO2_mass_traits,
  # Tell the model which column represents the light values (for the curve correction)
  lightValues = "PAR",
  # Set the three sub-models for the model components.  You probably want to keep the x-assymtote model as an intercept-only model
  # and turn off the multiplier model (by setting the multiplier to 1.0).  This means only the y-asymptote will change with the
  # environmental covariates + VegetationHeight + Richness + Evenness + Diversity + CWM_N + CWM_C + CWM_CN + CWM_LDMC + CWM_LT + CWM_LA + CWM_SLA + CWM_VH + FDis_Traits
  yAsymModel = GPP ~ T_summer * P_annual + Richness + Evenness + Diversity + CWM_N + CWM_C + CWM_CN + CWM_LDMC + CWM_LT + CWM_LA + CWM_SLA + CWM_VH + 
    F.Dispersion, 
  xAsymModel = ~ 1,
  multiplierModel = 1,
  # Tell the model that you want to do LASSO regularisation
  regCoeffs = "lasso" , 
  # Define the set of indirect models that you also want to fit (check how the sub-models are defined)
  indirectComponents = list(
    list(modelFormula = Reco ~ T_summer * P_annual + Richness + Evenness + Diversity + CWM_N + CWM_C + CWM_CN + CWM_LDMC + CWM_LT + CWM_LA + CWM_SLA + CWM_VH + F.Dispersion, errorFamily = Gamma(link = "log"),   
         regCoeffs = "lasso", 
         suffix = "_Reco_Model")
    ),
  # A vector of different PAR values you want to get predictions for (so that you can do PAR standardisation of the values)
  lightStandards = c(800),
  mcmcParams = list(numRuns = 1000000, thinDensity = 400, predictThinDensity = 400),
  numCores = 0)

CoutValues <- gppLightCurveCorrection(
  # Set the input data
  inputData = CO2_mass_traits,
  # Tell the model which column represents the light values (for the curve correction)
  lightValues = "PAR",
  # Set the three sub-models for the model components.  You probably want to keep the x-assymtote model as an intercept-only model
  # and turn off the multiplier model (by setting the multiplier to 1.0).  This means only the y-asymptote will change with the
  # environmental covariates + VegetationHeight + Richness + Evenness + Diversity + CWM_N + CWM_C + CWM_CN + CWM_LDMC + CWM_LT + CWM_LA + CWM_SLA + CWM_VH + FDis_Traits
  yAsymModel = GPP ~ T_summer * P_annual, 
  xAsymModel = ~ 1,
  multiplierModel = 1,
  # Tell the model that you want to do LASSO regularisation
  regCoeffs = "lasso" , 
  # Define the set of indirect models that you also want to fit (check how the sub-models are defined)
  indirectComponents = list(
    list(modelFormula = Reco ~ T_summer * P_annual, errorFamily = Gamma(link = "log"), 
         regCoeffs = "lasso", 
         suffix = "_Reco_Model")
  ),
  # A vector of different PAR values you want to get predictions for (so that you can do PAR standardisation of the values)
  lightStandards = c(800),
  mcmcParams = list(numRuns = 1000000, thinDensity = 400, predictThinDensity = 400),
  numCores = 0)

VoutValues <- gppLightCurveCorrection(
  # Set the input data
  inputData = CO2_mass_traits,
  # Tell the model which column represents the light values (for the curve correction)
  lightValues = "PAR",
  # Set the three sub-models for the model components.  You probably want to keep the x-assymtote model as an intercept-only model
  # and turn off the multiplier model (by setting the multiplier to 1.0).  This means only the y-asymptote will change with the
  # environmental covariates + VegetationHeight + Richness + Evenness + Diversity + CWM_N + CWM_C + CWM_CN + CWM_LDMC + CWM_LT + CWM_LA + CWM_SLA + CWM_VH + FDis_Traits
  yAsymModel = GPP ~ Richness + Evenness + Diversity, 
  xAsymModel = ~ 1,
  multiplierModel = 1,
  # Tell the model that you want to do LASSO regularisation
  regCoeffs = "lasso" , 
  # Define the set of indirect models that you also want to fit (check how the sub-models are defined)
  indirectComponents = list(
    list(modelFormula = Reco ~ Richness + Evenness + Diversity, errorFamily = Gamma(link = "log"),
         regCoeffs = "lasso", 
         suffix = "_Reco_Model")
  ),
  # A vector of different PAR values you want to get predictions for (so that you can do PAR standardisation of the values)
  lightStandards = c(800),
  mcmcParams = list(numRuns = 1000000, thinDensity = 400, predictThinDensity = 400),
  numCores = 0)

ToutValues <- gppLightCurveCorrection(
  # Set the input data
  inputData = CO2_mass_traits,
  # Tell the model which column represents the light values (for the curve correction)
  lightValues = "PAR",
  # Set the three sub-models for the model components.  You probably want to keep the x-assymtote model as an intercept-only model
  # and turn off the multiplier model (by setting the multiplier to 1.0).  This means only the y-asymptote will change with the
  # environmental covariates + VegetationHeight + Richness + Evenness + Diversity + CWM_N + CWM_C + CWM_CN + CWM_LDMC + CWM_LT + CWM_LA + CWM_SLA + CWM_VH + FDis_Traits
  yAsymModel = GPP ~ CWM_N + CWM_C + CWM_CN + CWM_LDMC + CWM_LT + CWM_LA + CWM_SLA + CWM_VH + F.Dispersion, 
  xAsymModel = ~ 1,
  multiplierModel = 1,
  # Tell the model that you want to do LASSO regularisation
  regCoeffs = "lasso" , 
  # Define the set of indirect models that you also want to fit (check how the sub-models are defined)
  indirectComponents = list(
    list(modelFormula = Reco ~ CWM_N + CWM_C + CWM_CN + CWM_LDMC + CWM_LT + CWM_LA + CWM_SLA + CWM_VH + F.Dispersion,
         errorFamily = Gamma(link = "log"),
         regCoeffs = "lasso", 
         suffix = "_Reco_Model")
  ),
  # A vector of different PAR values you want to get predictions for (so that you can do PAR standardisation of the values)
  lightStandards = c(800),
  mcmcParams = list(numRuns = 1000000, thinDensity = 400, predictThinDensity = 400),
  numCores = 0)

CToutValues <- gppLightCurveCorrection(
  # Set the input data
  inputData = CO2_mass_traits,
  # Tell the model which column represents the light values (for the curve correction)
  lightValues = "PAR",
  # Set the three sub-models for the model components.  You probably want to keep the x-assymtote model as an intercept-only model
  # and turn off the multiplier model (by setting the multiplier to 1.0).  This means only the y-asymptote will change with the
  # environmental covariates + VegetationHeight + Richness + Evenness + Diversity + CWM_N + CWM_C + CWM_CN + CWM_LDMC + CWM_LT + CWM_LA + CWM_SLA + CWM_VH + FDis_Traits
  yAsymModel = GPP ~ T_summer * P_annual + CWM_N + CWM_C + CWM_CN + CWM_LDMC + CWM_LT + CWM_LA + CWM_SLA + CWM_VH + F.Dispersion, 
  xAsymModel = ~ 1,
  multiplierModel = 1,
  # Tell the model that you want to do LASSO regularisation
  regCoeffs = "lasso" , 
  # Define the set of indirect models that you also want to fit (check how the sub-models are defined)
  indirectComponents = list(
    list(modelFormula = Reco ~ T_summer * P_annual + CWM_N + CWM_C + CWM_CN + CWM_LDMC + CWM_LT + CWM_LA + CWM_SLA + CWM_VH + F.Dispersion, errorFamily = Gamma(link = "log"),
         regCoeffs = "lasso", 
         suffix = "_Reco_Model")
  ),
  # A vector of different PAR values you want to get predictions for (so that you can do PAR standardisation of the values)
  lightStandards = c(800),
  mcmcParams = list(numRuns = 1000000, thinDensity = 400, predictThinDensity = 400),
  numCores = 0)

VToutValues <- gppLightCurveCorrection(
  # Set the input data
  inputData = CO2_mass_traits,
  # Tell the model which column represents the light values (for the curve correction)
  lightValues = "PAR",
  # Set the three sub-models for the model components.  You probably want to keep the x-assymtote model as an intercept-only model
  # and turn off the multiplier model (by setting the multiplier to 1.0).  This means only the y-asymptote will change with the
  # environmental covariates + VegetationHeight + Richness + Evenness + Diversity + CWM_N + CWM_C + CWM_CN + CWM_LDMC + CWM_LT + CWM_LA + CWM_SLA + CWM_VH + FDis_Traits
  yAsymModel = GPP ~ Richness + Evenness + Diversity + CWM_N + CWM_C + CWM_CN + CWM_LDMC + CWM_LT 
  + CWM_LA + CWM_SLA + CWM_VH + F.Dispersion, 
  xAsymModel = ~ 1,
  multiplierModel = 1,
  # Tell the model that you want to do LASSO regularisation
  regCoeffs = "lasso" , 
  # Define the set of indirect models that you also want to fit (check how the sub-models are defined)
  indirectComponents = list(
    list(modelFormula = Reco ~ Richness + Evenness + Diversity + CWM_N + CWM_C + CWM_CN + CWM_LDMC
         + CWM_LT + CWM_LA + CWM_SLA + CWM_VH + F.Dispersion, errorFamily = Gamma(link = "log"),
         regCoeffs = "lasso", 
         suffix = "_Reco_Model")
  ),
  # A vector of different PAR values you want to get predictions for (so that you can do PAR standardisation of the values)
  lightStandards = c(800),
  mcmcParams = list(numRuns = 1000000, thinDensity = 400, predictThinDensity = 400),
  numCores = 0)



CVoutValues <- gppLightCurveCorrection(
  # Set the input data
  inputData = CO2_mass_traits,
  # Tell the model which column represents the light values (for the curve correction)
  lightValues = "PAR",
  # Set the three sub-models for the model components.  You probably want to keep the x-assymtote model as an intercept-only model
  # and turn off the multiplier model (by setting the multiplier to 1.0).  This means only the y-asymptote will change with the
  # environmental covariates + VegetationHeight + Richness + Evenness + Diversity + CWM_N + CWM_C + CWM_CN + CWM_LDMC + CWM_LT + CWM_LA + CWM_SLA + CWM_VH + FDis_Traits
  yAsymModel = GPP ~ T_summer * P_annual + Richness + Evenness + Diversity, 
  xAsymModel = ~ 1,
  multiplierModel = 1,
  # Tell the model that you want to do LASSO regularisation
  regCoeffs = "lasso" , 
  # Define the set of indirect models that you also want to fit (check how the sub-models are defined)
  indirectComponents = list(
    list(modelFormula = Reco ~ T_summer * P_annual + Richness + Evenness + Diversity, errorFamily = Gamma(link = "log"),
         regCoeffs = "lasso", 
         suffix = "_Reco_Model")
  ),
  # A vector of different PAR values you want to get predictions for (so that you can do PAR standardisation of the values)
  lightStandards = c(800),
  mcmcParams = list(numRuns = 1000000, thinDensity = 400, predictThinDensity = 400),
  numCores = 0)

# Save the model files locally
saveRDS(CoutValues, file = file.path(curWorkspace, "GPP_RECO_C_11102022.rds"))
saveRDS(VoutValues, file = file.path(curWorkspace, "GPP_RECO_V_11102022.rds"))
saveRDS(ToutValues, file = file.path(curWorkspace, "GPP_RECO_T_19122021.rds"))
saveRDS(CVoutValues, file = file.path(curWorkspace, "GPP_RECO_CV_11102022.rds"))
saveRDS(VToutValues, file = file.path(curWorkspace, "GPP_RECO_VT_11102022.rds"))
saveRDS(CVToutValues, file = file.path(curWorkspace, "GPP_RECO_CVT_11102022.rds"))
saveRDS(CToutValues, file = file.path(curWorkspace, "GPP_RECO_CT_11102022.rds"))

CVToutValues$standardSummary

GPPpredicted<- as.data.frame(CVToutValues$standardSummary)
GPPpredicted<- cbind(GPPpredicted, CO2_mass_traits)

Reco_predictors<- GPPpredicted%>%
  select(lightStandard800.Mean, T_level, P_level, T_summer, P_annual, Richness, Evenness, Diversity, VegetationHeight, 
         CWM_N, CWM_C, CWM_CN, CWM_LDMC, CWM_LT, CWM_LA, CWM_SLA, CWM_VH, F.Dispersion)

Reco_predictors = melt(Reco_predictors, id.vars= c("lightStandard800.Mean", "T_level", "P_level"))
ggplot(Reco_predictors) +
  geom_point(aes(value, lightStandard800.Mean, fill=as.factor(T_level), shape=as.factor(P_level)))+ 
  #geom_smooth(aes(value, Reco), method=lm, se=FALSE, col="black") +
  scale_fill_manual(values= c("#FEE0D2", "#FC9272", "#DE2D26"), name = "Temperature", labels = c("Alpine", "Sub-alpine", "Boreal"))+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  scale_shape_manual(values= c(24, 22, 21, 25), name = "Precipitation", labels = c("600mm", "1200mm", "2000mm", "2700mm"))+
  facet_wrap(~variable, scales="free_x") +
  labs(x = "Predictor", y = "GPP")+
  theme_classic()



# Get convergence plots models
plot(CoutValues$parameterSamples) 
plot(CVToutValues$parameterSamples) 

# Retrieve the DHARMa plot of the model
plot(ToutValues$DHARMaResiduals) # GPP
plot(CVToutValues$DHARMaResiduals) # GPP
# Retrieve the DHARMa plot of the first indirect sub-model (Reco)
plot(ToutValues$indirectModelOutputs[[1]]$DHARMaResiduals)
plot(CVToutValues$indirectModelOutputs[[1]]$DHARMaResiduals)

# Retrieve the predictions of the Cmodel
CoutValues$predictionSummary
# Get the predictions standardised to the first light standard
CoutValues$standardSummary[[1]]

# Retrieve the predictions of the CVTmodel
CVToutValues$predictionSummary
# Get the predictions standardised to the first light standard
CVToutValues$standardSummary[[1]]

plotResiduals(CoutValues$DHARMaResiduals, CO2_mass_traits$Diversity)

obs_pred<-as.data.frame(cbind(CVToutValues$DHARMaResiduals$fittedPredictedResponse, CVToutValues$predictionSummary[,1], CoutValues$DHARMaResiduals$fittedPredictedResponse, CoutValues$predictionSummary[,1]))
obs_pred <- cbind(CO2_mass_traits, obs_pred)

p5<-ggplot(obs_pred, aes(GPP, V3, col=PAR))+
  geom_point()+
  scale_color_gradient(low="grey", high="red")+
  ylab("GPP_fitted")+
  geom_abline(slope=1, intercept = 0)+
  theme_classic()

p6<-ggplot(obs_pred, aes(GPP, V1, col=PAR))+
  geom_point()+
  scale_color_gradient(low="grey", high="red")+
  ylab("GPP fitted")+
  geom_abline(slope=1, intercept = 0)+
  theme_classic()


#Paramater values GPP models
GPPsummary_C<-CoutValues$parameterSummary #summary of GPP model
GPPsummary_C<- as.data.frame(GPPsummary_C)%>% 
  mutate_if(is.numeric, round, 3)
write.csv2(GPPsummary_C, file= file.path(curWorkspace, "GPP_Coutput_summary_20102022.csv"))
CoutValues$WAIC
CoutValues$rSquared
  
GPPsummary_V<-VoutValues$parameterSummary #summary of GPP model
GPPsummary_V<- as.data.frame(GPPsummary_V)%>% 
  mutate_if(is.numeric, round, 3)
write.csv2(GPPsummary_V, file= file.path(curWorkspace, "GPP_Voutput_summary_20102022.csv"))
VoutValues$WAIC
VoutValues$rSquared

GPPsummary_T<-ToutValues$parameterSummary #summary of GPP model
GPPsummary_T<- as.data.frame(GPPsummary_T)%>% 
  mutate_if(is.numeric, round, 3)
write.csv2(GPPsummary_T, file= file.path(curWorkspace, "GPP_Toutput_summary_20102022.csv"))
ToutValues$WAIC
ToutValues$rSquared

GPPsummary_CV<-CVoutValues$parameterSummary #summary of GPP model
GPPsummary_CV<- as.data.frame(GPPsummary_CV)%>% 
  mutate_if(is.numeric, round, 3)
write.csv2(GPPsummary_CV, file= file.path(curWorkspace, "GPP_CVoutput_summary_20102022.csv"))
CVoutValues$WAIC
CVoutValues$rSquared

GPPsummary_CT<-CToutValues$parameterSummary #summary of GPP model
GPPsummary_CT<- as.data.frame(GPPsummary_CT)%>% 
  mutate_if(is.numeric, round, 3)
write.csv2(GPPsummary_CT, file= file.path(curWorkspace, "GPP_CToutput_summary_20102022.csv"))
CToutValues$WAIC
CToutValues$rSquared


GPPsummary_VT<-VToutValues$parameterSummary #summary of GPP model
GPPsummary_VT<- as.data.frame(GPPsummary_VT)%>% 
  mutate_if(is.numeric, round, 3)
write.csv2(GPPsummary_VT, file= file.path(curWorkspace, "GPP_VToutput_summary_20102022.csv"))
VToutValues$WAIC
VToutValues$rSquared

GPPsummary_CVT<-CVToutValues$parameterSummary #summary of GPP model
GPPsummary_CVT<- as.data.frame(GPPsummary_CVT)%>% 
  mutate_if(is.numeric, round, 3)
write.csv2(GPPsummary_CVT, file= file.path(curWorkspace, "GPP_CVToutput_summary_20102022.csv"))
CVToutValues$WAIC
CVToutValues$rSquared

# Parameter values Reco models
Recosummary_C<-as.data.frame(CoutValues$indirectModelOutputs[[1]]$parameterSummary) #summary of Reco model
Recosummary_C<- Recosummary_C%>% 
  mutate_if(is.numeric, round, 3)
write.csv2(Recosummary_C, file= file.path(curWorkspace, "Reco_Coutput_summary_20102022.csv"))
CoutValues$indirectModelOutputs[[1]]$WAIC
CoutValues$indirectModelOutputs[[1]]$rSquared

Recosummary_V<-as.data.frame(VoutValues$indirectModelOutputs[[1]]$parameterSummary) #summary of Reco model
Recosummary_V<- Recosummary_V%>% 
  mutate_if(is.numeric, round, 3)
write.csv2(Recosummary_V, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_Voutput_summary_20102022.csv")
VoutValues$indirectModelOutputs[[1]]$WAIC
VoutValues$indirectModelOutputs[[1]]$rSquared

Recosummary_T<-as.data.frame(ToutValues$indirectModelOutputs[[1]]$parameterSummary) #summary of Reco model
Recosummary_T<- Recosummary_T%>% 
  mutate_if(is.numeric, round, 3)
write.csv2(Recosummary_T, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_Toutput_summary_20102022.csv")
ToutValues$indirectModelOutputs[[1]]$WAIC
ToutValues$indirectModelOutputs[[1]]$rSquared

Recosummary_CV<-as.data.frame(CVoutValues$indirectModelOutputs[[1]]$parameterSummary) #summary of Reco model
Recosummary_CV<- Recosummary_CV%>% 
  mutate_if(is.numeric, round, 3)
write.csv2(Recosummary_CV, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_CVoutput_summary_20102022.csv")
CVoutValues$indirectModelOutputs[[1]]$WAIC
CVoutValues$indirectModelOutputs[[1]]$rSquared

Recosummary_CT<-as.data.frame(CToutValues$indirectModelOutputs[[1]]$parameterSummary) #summary of Reco model
Recosummary_CT<- Recosummary_CT%>% 
  mutate_if(is.numeric, round, 3)
write.csv2(Recosummary_CT, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_CToutput_summary_20102022.csv")
CToutValues$indirectModelOutputs[[1]]$WAIC
CToutValues$indirectModelOutputs[[1]]$rSquared

Recosummary_VT<-as.data.frame(VToutValues$indirectModelOutputs[[1]]$parameterSummary) #summary of Reco model
Recosummary_VT<- Recosummary_VT%>% 
  mutate_if(is.numeric, round, 3)
write.csv2(Recosummary_VT, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_VToutput_summary_20102022.csv")
VToutValues$indirectModelOutputs[[1]]$WAIC
VToutValues$indirectModelOutputs[[1]]$rSquared

Recosummary_CVT<-as.data.frame(CVToutValues$indirectModelOutputs[[1]]$parameterSummary) #summary of Reco model
Recosummary_CVT<- Recosummary_CVT%>% 
  mutate_if(is.numeric, round, 3)
write.csv2(Recosummary_CVT, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_CVToutput_summary_20102022.csv")
CVToutValues$indirectModelOutputs[[1]]$WAIC
CVToutValues$indirectModelOutputs[[1]]$rSquared