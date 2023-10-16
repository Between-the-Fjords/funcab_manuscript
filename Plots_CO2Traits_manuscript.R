library(tidyverse)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

CO2_mass_traits<-readRDS(file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\CO2_traits06122021.rds")%>%
  mutate(GPP = -1*GPP)%>% # turn GPP positive for modelling 
  filter(PAR >200)

#GPP response to PAR, Supplementary figure 1
ggplot(CO2_mass_traits, aes(x=PAR, y= GPP, fill=as.factor(T_level), shape=as.factor(P_level)), size=2)+
  geom_point()+
  scale_fill_manual(values= c( "#FEE0D2", "#FC9272", "#DE2D26"), name = "Temperature", labels = c("Alpine", "Sub-alpine", "Boreal"))+ 
  guides(fill=guide_legend(override.aes=list(shape=21, size =3)))+
  scale_shape_manual(values= c(24, 22, 21, 25), name = "Precipitation", labels = c("600mm", "1200mm", "2000mm", "2700mm"))+
  xlab(expression(PAR~(micromol/m^{2}/s))) +
  ylab(expression(GPP~(micromol/m^{2}/s))) +
  theme_classic()+
  theme( axis.title = element_text(size = 14), axis.text = element_text(size =12), legend.text = element_text(size =11) )

# Figure 2a
C.indirect<-read.csv2("C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Ceffect_VT.csv")

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
  scale_color_manual(name = "Climate", labels= c("Temp:Prec", "Precipitation", "Temperature"), values = c("#2e4057", "#0099cc", "#DE2D26"))+
  scale_shape_manual(values=c(1,19), guide= "none")+
  geom_rect(data=rects, aes(ymin=-0.2, ymax=Inf, xmin=xstart, xmax=xend, fill=col), alpha =0.1)+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = c("white", "grey30"), guide ="none") +
  ylim(-0.2, 0.4)+
  coord_flip()+
  scale_y_continuous(limits = c(-0.2, 0.35), expand = c(0, 0))+
  scale_x_discrete(limits=c("FDis", "CN", "N", "C", "LDMC", "LT", "SLA", "LA", "VH", 
                             "Diversity", "Eveness", "Richness", "VegetationHeight"),
                   labels=c("F.Dispersion", "CN", "N", "C", "LDMC", "LT", "SLA", "LA", "VH", 
                             "Diversity", "Evenness", "Richness", "Veg Height"))+
  theme_classic()+
  theme(legend.position = c(0.8,0.3), strip.background = element_blank(), 
        strip.text.x= element_blank(), axis.title.x=element_blank(), 
        axis.text.x=element_text(size = 14), axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14), axis.title.y=element_blank(),
        title=element_text(size=14), legend.text=element_text(size=14))
Ceffect_VT



#PCA plot figure 2B
# calculate means to not duplicate values for plots due to multiple flux measurements on plot per year
CO2_mass_traits_median <- CO2_mass_traits%>%
  dplyr::rename(Site =site)%>%
  dplyr::group_by(turfID, Site, year, P_level, T_level, treatment) %>%
  #summarise(no_rows = length(turfID))%>%
  dplyr::summarise_all(median, na.rm=TRUE)%>%
  dplyr::rename(Year = year)%>%
  ungroup()

apply(CO2_mass_traits_median, 2, function(x) length(which(!is.na(x)))) # check number of plots

library(ggfortify)
library(devtools)
library(factoextra)
Traitdata<- CO2_mass_traits_median%>%
  select(Richness, Evenness, Diversity, "Veg Height"= VegetationHeight, N = CWM_N , C =CWM_C, CN = CWM_CN, LDMC = CWM_LDMC, LT = CWM_LT, LA = CWM_LA, SLA = CWM_SLA, VH = CWM_VH, FDispersion = F.Dispersion)
TraitPCA <- princomp(Traitdata, cor= TRUE, scores=TRUE) #, Temperature = T_summer_longterm, Precipitation = P_annual_longterm

PCAplot<- autoplot(TraitPCA, data = CO2_mass_traits_median, fill= "T_level", shape = "P_level", size = 3,
                   loadings = TRUE, loadings.colour = 'black', loadings.label.colour = "black", 
                   loadings.label = TRUE, loadings.label.size = 5, loadings.label.vjust = 1.2, loadings.label.hjust = 0.6)+
  scale_fill_manual(values= c( "#FEE0D2", "#FC9272", "#DE2D26"), name = "Temperature", labels = c("Alpine", "Sub-alpine", "Boreal"))+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  scale_shape_manual(values= c(24, 22, 21, 25), name = "Precipitation", labels = c("600mm", "1200mm", "2000mm", "2700mm"))+
  theme_classic()+
  theme(legend.position = "left", title=element_text(size=14), axis.text = element_text(size=13), legend.text=element_text(size=14), legend.title = element_text(size=14))
PCAplot

# combine figure 2A and 2B
library(cowplot)
CVT_PCA <- plot_grid(Ceffect_VT, PCAplot, labels = c("A", "B"),   hjust = -1,   ncol = 1 )
CVT_PCA

# load dataframes to create figure 3A and 3B
# density plots made from full CVT model output
# points to indicate effect size and signifance from Climate (C) only, Vegetation structure (V) only, and Traits only (T) model
CoutValues<-readRDS("C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_RECO_C_20102022.rds")
VoutValues<-readRDS("C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_RECO_V_20102022.rds")
ToutValues<-readRDS("C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_RECO_T_20102022.rds")
CVoutValues<-readRDS("C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_RECO_CV_20102022.rds")
CToutValues<-readRDS("C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_RECO_CT_20102022.rds")
VToutValues<-readRDS("C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_RECO_VT_20102022.rds")
CVToutValues<-readRDS("C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_RECO_CVT_20102022.rds")

#Paramater values GPP models for creating supplementary material tables 
GPPsummary_C<-CoutValues$parameterSummary #summary of GPP model
GPPsummary_C<- as.data.frame(GPPsummary_C)%>% 
  mutate(model = "C")
#write.csv2(GPPsummary_C, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_Coutput_summary_20102022.csv")
CoutValues$WAIC
CoutValues$rSquared

GPPsummary_V<-VoutValues$parameterSummary #summary of GPP model
GPPsummary_V<- as.data.frame(GPPsummary_V)%>% 
  mutate(model = "V")
#write.csv2(GPPsummary_V, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_Voutput_summary_20102022.csv")
VoutValues$WAIC
VoutValues$rSquared

GPPsummary_T<-ToutValues$parameterSummary #summary of GPP model
GPPsummary_T<- as.data.frame(GPPsummary_T)%>% 
  mutate(model = "T")
#write.csv2(GPPsummary_T, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_Toutput_summary_20102022.csv")
ToutValues$WAIC
ToutValues$rSquared

GPPsummary_CV<-CVoutValues$parameterSummary #summary of GPP model
GPPsummary_CV<- as.data.frame(GPPsummary_CV)%>% 
  mutate(model = "CV")
#write.csv2(GPPsummary_CV, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_CVoutput_summary_20102022.csv")
CVoutValues$WAIC
CVoutValues$rSquared

GPPsummary_CT<-CToutValues$parameterSummary #summary of GPP model
GPPsummary_CT<- as.data.frame(GPPsummary_CT)%>% 
  mutate(model = "CT")
#write.csv2(GPPsummary_CT, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_CToutput_summary_20102022.csv")
CToutValues$WAIC
CToutValues$rSquared

GPPsummary_VT<-VToutValues$parameterSummary #summary of GPP model
GPPsummary_VT<- as.data.frame(GPPsummary_VT)%>% 
  mutate(model = "VT")
#write.csv2(GPPsummary_VT, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_VToutput_summary_20102022.csv")
VToutValues$WAIC
VToutValues$rSquared

GPPsummary_CVT<-CVToutValues$parameterSummary #summary of GPP model
GPPsummary_CVT<- as.data.frame(GPPsummary_CVT)%>% 
  mutate(model = "CVT")
#write.csv2(GPPsummary_CVT, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_CVToutput_summary_20102022.csv")
CVToutValues$WAIC
CVToutValues$rSquared

# Parameter values Reco models
Recosummary_C<-as.data.frame(CoutValues$indirectModelOutputs[[1]]$parameterSummary) #summary of Reco model
Recosummary_C<- Recosummary_C%>% 
  mutate(model = "C")
#write.csv2(Recosummary_C, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_Coutput_summary_20102022.csv")
CoutValues$indirectModelOutputs[[1]]$WAIC
CoutValues$indirectModelOutputs[[1]]$rSquared

Recosummary_V<-as.data.frame(VoutValues$indirectModelOutputs[[1]]$parameterSummary) #summary of Reco model
Recosummary_V<- Recosummary_V%>% 
  mutate(model = "V")
#write.csv2(Recosummary_V, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_Voutput_summary_20102022.csv")
VoutValues$indirectModelOutputs[[1]]$WAIC
VoutValues$indirectModelOutputs[[1]]$rSquared

Recosummary_T<-as.data.frame(ToutValues$indirectModelOutputs[[1]]$parameterSummary) #summary of Reco model
Recosummary_T<- Recosummary_T%>% 
  mutate(model = "T")
#write.csv2(Recosummary_T, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_Toutput_summary_20102022.csv")
ToutValues$indirectModelOutputs[[1]]$WAIC
ToutValues$indirectModelOutputs[[1]]$rSquared

Recosummary_CV<-as.data.frame(CVoutValues$indirectModelOutputs[[1]]$parameterSummary) #summary of Reco model
Recosummary_CV<- Recosummary_CV%>% 
  mutate(model = "CV")
#write.csv2(Recosummary_CV, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_CVoutput_summary_20102022.csv")
CVoutValues$indirectModelOutputs[[1]]$WAIC
CVoutValues$indirectModelOutputs[[1]]$rSquared

Recosummary_CT<-as.data.frame(CToutValues$indirectModelOutputs[[1]]$parameterSummary) #summary of Reco model
Recosummary_CT<- Recosummary_CT%>% 
  mutate(model = "CT")
#write.csv2(Recosummary_CT, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_CToutput_summary_20102022.csv")
CToutValues$indirectModelOutputs[[1]]$WAIC
CToutValues$indirectModelOutputs[[1]]$rSquared

Recosummary_VT<-as.data.frame(VToutValues$indirectModelOutputs[[1]]$parameterSummary) #summary of Reco model
Recosummary_VT<- Recosummary_VT%>% 
  mutate(model = "VT")
write.csv2(Recosummary_VT, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_VToutput_summary_20102022.csv")
VToutValues$indirectModelOutputs[[1]]$WAIC
VToutValues$indirectModelOutputs[[1]]$rSquared

Recosummary_CVT<-as.data.frame(CVToutValues$indirectModelOutputs[[1]]$parameterSummary) #summary of Reco model
Recosummary_CVT<- Recosummary_CVT%>% 
  mutate(model = "CVT")
#write.csv2(Recosummary_CVT, file= "C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_CVToutput_summary_20102022.csv")
CVToutValues$indirectModelOutputs[[1]]$WAIC
CVToutValues$indirectModelOutputs[[1]]$rSquared


#### Bind together output of Climate (C), Vegetation structure (V), Traits (T) model for GPP and Reco and determine significance for different predictors for these submodels
GPPsummary_C<-read.csv2("C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_Coutput_summary_20102022.csv")%>% 
  mutate(model = "C")
GPPsummary_V<-read.csv2("C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_Voutput_summary_20102022.csv") %>% 
  mutate(model = "V")
GPPsummary_T<-read.csv2("C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_Toutput_summary_20102022.csv")%>% 
  mutate(model = "T")
GPPsummary_CVT<-read.csv2("C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\GPP_CVToutput_summary_20102022.csv")%>% 
  mutate(model = "CVT")

GPP_C_V_T_separate<-rbind(GPPsummary_C, GPPsummary_V, GPPsummary_T, GPPsummary_CVT)%>%
  mutate(predictor = X)%>%
  filter(!grepl("intercept", predictor))%>%
  filter(!grepl("lasso", predictor))%>%
  filter(!grepl("SD", predictor))%>%
  mutate(covName=  sub(pattern = "logyAssymCoeff", "", predictor))%>%
  #mutate(covName=  sub(pattern = "1", "", covName))%>%
  mutate(CIlow_95 = as.numeric(X95.CI_low),
         CIupp_95 = as.numeric(X95.CI_upp),
         significance= ifelse(CIlow_95 < 0 & CIupp_95 < 0, "significant", 
                              ifelse(CIlow_95 > 0 & CIupp_95 > 0, "significant", "non-significant")))%>%
  select(covName, Mean, significance, CIlow_95, CIupp_95, model)

Recosummary_C<-read.csv2("C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_Coutput_summary_20102022.csv")%>% 
  mutate(model = "C")
Recosummary_V<-read.csv2("C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_Voutput_summary_20102022.csv") %>% 
  mutate(model = "V")
Recosummary_T<-read.csv2("C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_Toutput_summary_20102022.csv")%>% 
  mutate(model = "T")
Recosummary_CVT<-read.csv2("C:\\Users\\ialt\\OneDrive - NORCE\\FunCab\\Data\\FunCaB2\\TraitCO2\\Analysis\\Reco_CVToutput_summary_20102022.csv")%>%   mutate(model = "CVT")


Reco_C_V_T_separate<-rbind(Recosummary_C, Recosummary_V, Recosummary_T, Recosummary_CVT)%>%
  mutate(predictor = X)%>%
  filter(!grepl("intercept", predictor))%>%
  filter(!grepl("lasso", predictor))%>%
  filter(!grepl("SD", predictor))%>%
  mutate(covName=  sub(pattern = "Coeff", "", predictor))%>%
  mutate(covName=  sub(pattern = "1", "", covName))%>%
  mutate(CIlow_95 = as.numeric(X95.CI_low),
         CIupp_95 = as.numeric(X95.CI_upp),
         significance= ifelse(CIlow_95 < 0 & CIupp_95 < 0, "significant", 
                              ifelse(CIlow_95 > 0 & CIupp_95 > 0, "significant", "non-significant")))%>%
  select(covName, Mean, significance, CIlow_95, CIupp_95,  model)


## NEW FIGURE
# create dataframe for bands in figure 2A
rects <- data.frame(ystart = seq(0.5,15.5,1), yend = seq(1.5,16.5,1), 
                    col = rep(c("a", "b")))

#combine CVT model output with C, V and T only models
GPP_CVT_figure<- ggplot() +
  geom_point(data = GPP_C_V_T_separate, aes(x = Mean, y = covName, shape = significance, color= model), 
             position = position_dodge(width = 0.8), size = 3)+
  geom_errorbar(data = GPP_C_V_T_separate, aes(y= covName, xmin = CIlow_95, xmax= CIupp_95, shape = significance, color= model), 
                position = position_dodge(width = 0.8), alpha= 0.6, width= 0, size = 1)+
  scale_shape_manual(values=c(1,19), guide = "none")+
  scale_color_manual(values= c("#fc8d62","#66c2a5",  "#8da0cb", "#636363"), 
                       breaks= c("C", "V", "T", "CVT"),
                       labels= c("C only", "V only", "T only", "CVT combined"))+
  scale_fill_manual(values = c("white", "grey30"), guide ="none") +
  geom_vline(xintercept = 0.0, colour = "grey", size = 1) + 
  geom_hline(yintercept = c(9.5, 13.5), linetype = "dotted", size =1) +
  geom_rect(data=rects, aes(ymin=ystart, ymax=yend, xmin=-Inf, xmax=Inf, fill=col), alpha =0.1)+
  scale_y_discrete(limits=c("F_Dispersion", "CWM_CN", "CWM_N", "CWM_C", "CWM_LDMC", "CWM_LT", "CWM_SLA", "CWM_LA",
                            "CWM_VH",  "Diversity", "Evenness",
                            "Richness","VegetationHeight", "T_summer_P_annual", "P_annual", "T_summer"),
                   labels = c("F_Dispersion" = "F.Dispersion", "CWM_CN" = "CN", "CWM_N" = "N", "CWM_C" = "C", "CWM_LDMC" = "LDMC",
                              "CWM_LT" = "LT", "CWM_SLA" ="SLA", "CWM_LA" = "LA", "CWM_VH" = "VH", "Diversity", "Evenness",
                              "Richness","Vegetation Height", "T_summer_P_annual" = "Temp:Prec", 
                              "P_annual" = "Precipitation", "T_summer" = "Temperature")) +
  ggtitle("GPP")+
  theme_classic() + xlab("Scaled coefficient value") +
  theme(strip.background = element_blank(), strip.text.x= element_blank(), axis.title.x=element_text(size = 14), 
        axis.text.x=element_text(size = 14), axis.title = element_text(size = 12), axis.text.y = element_text(size = 14), 
        axis.title.y=element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5))

Reco_CVT_figure<- ggplot() +
  geom_point(data = Reco_C_V_T_separate, aes(x = Mean, y = covName, shape = significance, color= model), 
             position = position_dodge(width = 0.8), size = 3)+
  geom_errorbar(data = Reco_C_V_T_separate, aes(y= covName, xmin = CIlow_95, xmax= CIupp_95, shape = significance, color= model), 
                position = position_dodge(width = 0.8), alpha= 0.6, width= 0, size = 1)+
  scale_shape_manual(values=c(1,19))+
  scale_color_manual(values= c("#fc8d62","#66c2a5",  "#8da0cb", "#636363"), 
                     breaks= c("C", "V", "T", "CVT"),
                     labels= c("C only", "V only", "T only", "CVT combined"))+
  scale_fill_manual(values = c("white", "grey30"), guide ="none") +
  geom_vline(xintercept = 0.0, colour = "grey", size = 1) + 
  geom_hline(yintercept = c(9.5, 13.5), linetype = "dotted", size =1) +
  geom_rect(data=rects, aes(ymin=ystart, ymax=yend, xmin=-Inf, xmax=Inf, fill=col), alpha =0.1)+
  scale_y_discrete(limits=c("F_Dispersion", "CWM_CN", "CWM_N", "CWM_C", "CWM_LDMC", "CWM_LT", "CWM_SLA", "CWM_LA",
                            "CWM_VH",  "Diversity", "Evenness",
                            "Richness", "VegetationHeight", "T_summer_P_annual", "P_annual", "T_summer"),
                   labels = c("F_Dispersion" = "F.Dispersion", "CWM_CN" = "CN", "CWM_N" = "N", "CWM_C" = "C", "CWM_LDMC" = "LDMC",
                              "CWM_LT" = "LT", "CWM_SLA" ="SLA", "CWM_LA" = "LA", "CWM_VH" = "VH", "Diversity", "Evenness",
                              "Richness", "Vegetation Height","T_summer_P_annual" = "Temp:Prec", 
                              "P_annual" = "Precipitation", "T_summer" = "Temperature")) +
  ggtitle("Reco")+
  theme_classic() + xlab("Scaled coefficient value") +
  theme(strip.background = element_blank(), strip.text.x= element_blank(), axis.title.x=element_text(size = 14), 
        axis.text.x=element_text(size = 14), axis.title = element_text(size = 14), axis.text.y=element_blank(), 
        axis.title.y=element_blank(), legend.position = "right", plot.title = element_text(hjust = 0.5))

# Combine GPP and Reco density plots
library(cowplot)
prow <- plot_grid(GPP_CVT_figure+ theme(legend.position="none"),
                  Reco_CVT_figure,
                  align = 'vh',  labels = c("A", "B"),   hjust = -1,   nrow = 1 )


