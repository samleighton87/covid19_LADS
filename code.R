library(readr)
library(mice)
library(doParallel)
library(VIM)
library(corrplot)
library(plyr)
library(geofacet)
library(ggplot2)
library(viridis) # colour blind friendly palette, works in B&W also
library(dotwhisker)
library(broom)
library(dplyr)
library(car)
library(systemfit) #for seemingly unrelated regression
library(multcomp)  #for hypothesis tests of models coefficients
library(RCurl)

#Robust standard errors in summary function (only works single imputation)
url_robust <- "https://raw.githubusercontent.com/IsidoreBeautrelet/economictheoryblog/master/robust_summary.R"
eval(parse(text = getURL(url_robust, ssl.verifypeer = FALSE)),
     envir=.GlobalEnv)

#enable multicore (windows) which roughly halfs time for analysis runs
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

options(scipen=999)

#################################################################

covid_all = read_csv("CovidDataNew.csv")

covid_all$Country = as.factor(covid_all$Country)
covid_all$Local_Authority = as.factor(covid_all$Local_Authority)
covid_all$Area_Code = as.factor(covid_all$Area_Code)

#multiple imputation https://www.r-bloggers.com/imputing-missing-data-with-r-mice-package/

covid_death_all = covid_all[ ,!(colnames(covid_all) %in% c("Death_Registrations_Covid_Week_16","Death_Registrations_Care_Home_Covid_Week_16",
                                                               "Death_Registrations_Care_Home_Covid_Week_16_Per_100000",
                                                               "Death_Registrations_Hospital_Covid_Week_16",
                                                               "Death_Registrations_Hospital_Covid_Week_16_Per_100000",
                                                               "Death_Registrations_Home_Covid_Week_16",
                                                               "Death_Registrations_Home_Covid_Week_16_Per_100000",
                                                               "Death_Registrations_OoHosp_Week_16","Death_Registrations_OoHosp_Week_16_Per_100000"
                                                           ))]

#summary of missing data - max missing 5.5%
aggr_plot <- aggr(covid_death_all, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(covid_death_all), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

covid_death_all_data = covid_death_all[,-c(1,5:8)]


#standardise everything except outcome
covid_death_all_data[,-29] = scale(covid_death_all_data[,-29])
summary(covid_death_all_data)

tempData_all <- mice(covid_death_all_data,seed=987)
summary(tempData_all)

#Visualise predictor correlation
Data_imp = mice::complete(tempData_all,1)
Data_imp = plyr::rename(Data_imp, c("Density_People_sq_km"="Population Density","Local_Proportion_Decile_1_Z_Ind_IMD"="Deprivation",
                     "Females_2018"="Women","Lives_In_Communal"="Communal Living","BAME_Including_Jews"="BAME",
                     "Median_Age_2018"="Age","Over_1_5_Per_Room"="Crowded Living","Obese_Z_Ind"="Obesity",
                    "Current_Smoking_Z_Ind"="Smoking (current)","Ex_Smoker_Z_Ind"="Ex-Smokers","Tonnes_Port_Over_750"="Port Activity",
                   "Passengers_To_From_Pop_Stand_Normalised"="Flight Passengers","Age_Standardised_CVD_Deaths_U75_100000"="CVD deaths (<75)",
                     "Diabetes_Prevalence_Z_Ind"="Diabetes","RA_Z_Ind"="Rheumatoid Arthritis","HTN_Z_Ind"="Hypertension","COPD_Z_Ind"="COPD","Cancer_Z_Ind"="Cancer",
                   "CKD_Z_Ind"="CKD", "Average_Humid_23_2_1_3_8_3_15_3"="Humidity","Average_tempHighCel_23_2_1_3_8_3_15_3"="Temperature", 
                   "Coronavirus_Google_010220_190320"="Google Searching","Dementia_Z_Ind"="Dementia","Self_Funders_2017"="Self-Funding Care Home", 
                   "PM_2_5_Total_2018"= "Air Pollution"))
M = cor(Data_imp[,-29], method = "pearson")
diag(M) = NA
res1 = cor.mtest(Data_imp[,-29], conf.level = .95)
pdf("corr_predictors.pdf",height=10,width=10)
corrplot(M, p.mat = res1$p, order = "hclust", sig.level = c(.001, .01, .05), pch.cex = .9,
         insig = "label_sig", pch.col = "white", addrect = 4, method = "color",na.label = " ", col = viridis(100), tl.col = "black")
dev.off()

#not care home funding excecpt for care home model
tempData_all <- mice(covid_death_all_data,seed=987)
summary(tempData_all)


modelFit_covid_death_all = with(tempData_all, lm(Death_Registrations_Covid_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
                                                  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                                                   Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                                                  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                                                  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                                                   Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                                                  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018
                                                ))
View(summary(pool(modelFit_covid_death_all)))
pool.r.squared(modelFit_covid_death_all, adjusted = T)

#post regression diagnostics and sensitivity analyses with robust standard errors
Data_imp = mice::complete(tempData_all,1)

modTestRes = lm(Death_Registrations_Covid_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
                  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018, Data_imp)
plot(modTestRes)
summary(modTestRes, robust = F)
summary(modTestRes, robust = T)

data = summary(pool(modelFit_covid_death_all))
data <- tibble::rownames_to_column(data, "term")

dwplot(data)
#look at 95CIs
View(data.frame(data$term,data$estimate,(data$estimate-(1.96*data$std.error)),(data$estimate+(1.96*data$std.error)),data$p.value))

# Create list of brackets (label, topmost included predictor, bottommost included predictor)
five_brackets <- list(c("Country", "Scotland cf. England", "Wales cf. England"), 
                      c("Weather/Air", "Air Pollution", "Temperature"),
                      c("Demographics/Social", "Age", "Self Funding Care Home"),
                      c("Population Health Rates","Cancer","Ex-Smokers"),
                      c("Other Factors","Google Searching","Port Activity"))
pdf("Model_All_final3.pdf",11,7)
{dwplot(data, 
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2),
        dot_args = list(color = "#F8766D"), # color for the dot
        whisker_args = list(color = "Grey47")   # color for the whisker
        ) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(Scotland = "Scotland cf. England",                       # relabel predictors
                         Wales = "Wales cf. England",
                         PM_2_5_Total_2018 = "Air Pollution",
                         Average_Humid_23_2_1_3_8_3_15_3 = "Humidity", 
                         Average_tempHighCel_23_2_1_3_8_3_15_3 = "Temperature", 
                         Median_Age_2018 = "Age", 
                         BAME_Including_Jews = "BAME",
                         Lives_In_Communal = "Communal Living",
                         Over_1_5_Per_Room = "Crowded Living",
                         Local_Proportion_Decile_1_Z_Ind_IMD = "Deprivation",
                         Density_People_sq_km = "Population Density",
                         Females_2018 = "Women",
                         Self_Funders_2017 = "Self Funding Care Home",
                         Cancer_Z_Ind = "Cancer",
                         CKD_Z_Ind = "CKD",
                         COPD_Z_Ind = "COPD",
                         Age_Standardised_CVD_Deaths_U75_100000 = "CVD deaths (<75)",
                         Dementia_Z_Ind = "Dementia",
                         Diabetes_Prevalence_Z_Ind = "Diabetes",
                         HTN_Z_Ind = "Hypertension",
                         Obese_Z_Ind = "Obesity",
                         RA_Z_Ind = "Rheumatoid Arthritis",
                         Current_Smoking_Z_Ind = "Smoking (current)",
                         Ex_Smoker_Z_Ind = "Ex-Smokers",
                         Coronavirus_Google_010220_190320 = "Google Searching",
                         Passengers_To_From_Pop_Stand_Normalised = "Flight Passengers",
                         Tonnes_Port_Over_750 = "Port Activity"
                         )) +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    ggtitle("Local Authority Covid-19 Deaths per 100,000 in All Settings") + geom_text(aes(label = "Adjusted R² = 0.46\n(0.38, 0.53)"), x=9,y=2, size = 3)+
    theme(plot.title = element_text(face="bold", hjust = 0.5), legend.position = "none",plot.margin = unit(c(0,5,1,1), "cm"))+
    geom_text(x = 13.1,
              y=27, label = " 0.98 (-1.28, 3.25)", #Scotland
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=26, label = " 0.89  (-0.66, 2.45)", #Wales
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=25, label = " 8.18  (4.51, 11.84)   ***", #pm2.5
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=24, label = "-0.13 (-2.42, 2.16)", #Humidity
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=23, label = "-5.60 (-9.12, -2.08) **", #temp
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=22, label = " 1.43 (-2.78, 5.63)", #age
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=21, label = "  3.82 (-0.57, 8.21)", #BAME
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=20, label = "-1.64 (-3.28, 0.01)", #communal
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=19, label = "-1.62 (-5.49, 2.24)", #crowded
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=18, label = "  3.94 (1.70, 6.18)     **", #depriv
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=17, label = "  1.62 (-1.05, 4.29)", #density
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=16, label = "  3.42 (1.64, 5.20)     ***", #women
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=15, label = "-2.89 (-5.01, -0.78) **", #self funding
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=14, label = "  0.92 (-2.76, 4.60)", #cancer
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=13, label = "-0.07 (-3.64, 3.49)", #ckd
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=12, label = "-1.88 (-5.06, 1.29)", #copd
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=11, label = "-2.72 (-5.57, 0.13)", #cvd
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=10, label = "-2.75 (-5.02, -0.47) *", #dementia
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=9, label = "  1.21 (-1.05, 3.47)", #diabetes
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=8, label = "  1.65 (-1.83, 5.13)", #BP
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=7, label = "  0.65 (-1.96, 3.27)", #obesity
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=6, label = "-0.36 (-4.50, 3.78)", #RA
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=5, label = "-1.67 (-4.31, 0.98)", #smoking current
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=4, label = "-3.40 (-5.80, -1.00) **", #ex-smoking
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=3, label = "-0.26 (-1.89, 1.37)", #google
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=2, label = "  1.93 (-0.12, 3.98)", #flight passengers
              hjust = 0,
              size = 3) +
    geom_text(x = 13.1,
              y=1, label = "-0.45 (-1.69, 0.79)", #port activity
              hjust = 0,
              size = 3) +
    coord_cartesian(clip = 'off')  # This keeps the labels from disappearing
  
  } %>% 
  add_brackets(five_brackets)
dev.off()

################################################

covid_death_care = covid_all[ ,!(colnames(covid_all) %in% c("Death_Registrations_Covid_Week_16","Death_Registrations_Care_Home_Covid_Week_16",
                                                           "Death_Registrations_Covid_Week_16_Per_100000",
                                                           "Death_Registrations_Hospital_Covid_Week_16",
                                                           "Death_Registrations_Hospital_Covid_Week_16_Per_100000",
                                                           "Death_Registrations_Home_Covid_Week_16",
                                                           "Death_Registrations_Home_Covid_Week_16_Per_100000",
                                                           "Death_Registrations_OoHosp_Week_16","Death_Registrations_OoHosp_Week_16_Per_100000"))]

#summary of missing data - max missing 5.5%

covid_death_care_data = covid_death_care[,-c(1,5:8)]


#standardise everything except outcome
covid_death_care_data[,-29] = scale(covid_death_care_data[,-29])
summary(covid_death_care_data)

tempData_care <- mice(covid_death_care_data,seed=987)
modelFit_covid_death_care = with(tempData_care, lm(Death_Registrations_Care_Home_Covid_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
                                                     Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+BAME_Including_Jews+Lives_In_Communal+
                                                     Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                                                     Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                                                     Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                                                     Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                                                     Dementia_Z_Ind +Self_Funders_2017+PM_2_5_Total_2018
                                                 ))
summary(pool(modelFit_covid_death_care))
pool.r.squared(modelFit_covid_death_care, adjusted = T)

#post regression diagnostics and sensitivity analyses with robust standard errors
Data_imp = mice::complete(tempData_care,1)

modTestRes = lm(Death_Registrations_Care_Home_Covid_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
                  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018, Data_imp)
plot(modTestRes)
summary(modTestRes, robust = F)
summary(modTestRes, robust = T)


data = summary(pool(modelFit_covid_death_care))
data <- tibble::rownames_to_column(data, "term")

dwplot(data)
#look at 95CIs
View(data.frame(data$term,data$estimate,(data$estimate-(1.96*data$std.error)),(data$estimate+(1.96*data$std.error)),data$p.value))

# Create list of brackets (label, topmost included predictor, bottommost included predictor)
five_brackets <- list(c("Country", "Scotland cf. England", "Wales cf. England"), 
                      c("Weather/Air", "Air Pollution", "Temperature"),
                      c("Demographics/Social", "Age", "Self Funding Care Home"),
                      c("Population Health Rates","Cancer","Ex-Smokers"),
                      c("Other Factors","Google Searching","Port Activity"))
pdf("Model_Care_final3.pdf",11,7)
{dwplot(data, 
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2),
        dot_args = list(color = "#F8766D"), # color for the dot
        whisker_args = list(color = "Grey47")   # color for the whisker
) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(Scotland = "Scotland cf. England",                       # relabel predictors
                         Wales = "Wales cf. England",
                         PM_2_5_Total_2018 = "Air Pollution",
                         Average_Humid_23_2_1_3_8_3_15_3 = "Humidity", 
                         Average_tempHighCel_23_2_1_3_8_3_15_3 = "Temperature", 
                         Median_Age_2018 = "Age", 
                         BAME_Including_Jews = "BAME",
                         Lives_In_Communal = "Communal Living",
                         Over_1_5_Per_Room = "Crowded Living",
                         Local_Proportion_Decile_1_Z_Ind_IMD = "Deprivation",
                         Density_People_sq_km = "Population Density",
                         Females_2018 = "Women",
                         Self_Funders_2017 = "Self Funding Care Home",
                         Cancer_Z_Ind = "Cancer",
                         CKD_Z_Ind = "CKD",
                         COPD_Z_Ind = "COPD",
                         Age_Standardised_CVD_Deaths_U75_100000 = "CVD deaths (<75)",
                         Dementia_Z_Ind = "Dementia",
                         Diabetes_Prevalence_Z_Ind = "Diabetes",
                         HTN_Z_Ind = "Hypertension",
                         Obese_Z_Ind = "Obesity",
                         RA_Z_Ind = "Rheumatoid Arthritis",
                         Current_Smoking_Z_Ind = "Smoking (current)",
                         Ex_Smoker_Z_Ind = "Ex-Smokers",
                         Coronavirus_Google_010220_190320 = "Google Searching",
                         Passengers_To_From_Pop_Stand_Normalised = "Flight Passengers",
                         Tonnes_Port_Over_750 = "Port Activity"
    )) +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    ggtitle("Local Authority Covid-19 Deaths per 100,000 in Care Homes") + geom_text(aes(label = "Adjusted R² = 0.25\n(0.18, 0.34)"), x=1.75,y=2, size = 3)+
    theme(plot.title = element_text(face="bold", hjust = 0.5), legend.position = "none",plot.margin = unit(c(0,5,1,1), "cm"))+
    geom_text(x = 2.45,
              y=27, label = "  1.20 (0.38, 2.02)     **", #Scotland
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=26, label = "-0.11 (-0.68, 0.46)", #Wales
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=25, label = "  0.71 (-0.66, 2.08)", #air
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=24, label = "  0.71 (-0.09, 1.52)", #Humidity
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=23, label = "-1.08 (-2.34, 0.18)", #temp
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=22, label = "-0.14 (-1.69, 1.41)", #age
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=21, label = "-1.19 (-2,82, 0.44)", #BAME
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=20, label = "-0.22 (-0.79, 0.35)", #communal
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=19, label = "  0.73 (-0.70, 2.16)", #crowded
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=18, label = "  1.12 (0.31, 1.93)     **", #depriv
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=17, label = "  0.61 (-0.38, 1.60)", #density
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=16, label = "  1.03 (0.36, 1.69)     **", #women
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=15, label = "  0.22 (-0.55, 1.00)", #self funding
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=14, label = "-1.10 (-2.52, 0.32)", #cancer
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=13, label = "  0.37 (-0.94, 1.68", #ckd
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=12, label = "-0.25 (-1.37, 0.87)", #copd
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=11, label = "-0.63 (-1.66, 0.40)", #cvd
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=10, label = "-0.32 (-1.16, 0.51)", #dementia
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=9, label = "-0.87 (-1.70, -0.05) *", #diabetes
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=8, label = "  0.08 (-1.18, 1.34)", #BP
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=7, label = "  0.88 (-0.08, 1.83)", #obesity
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=6, label = "  0.49 (-1.11, 2.08)", #RA
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=5, label = "-0.84 (-1.76, 0.08)", #smoking current
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=4, label = "-1.47 (-2.38, -0.56) **", #ex-smoking
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=3, label = "  0.03 (-0.57, 0.62)", #google
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=2, label = "  0.15 (-0.60, 0.91)", #flight passengers
              hjust = 0,
              size = 3) +
    geom_text(x = 2.45,
              y=1, label = "-0.35 (-0.80, 0.11)", #port activity
              hjust = 0,
              size = 3) +
    coord_cartesian(clip = 'off')  # This keeps the labels from disappearing
  } %>% 
  add_brackets(five_brackets)
dev.off()


################################################

covid_death_hosp = covid_all[ ,!(colnames(covid_all) %in% c("Death_Registrations_Covid_Week_16","Death_Registrations_Care_Home_Covid_Week_16",
                                                            "Death_Registrations_Covid_Week_16_Per_100000","Death_Registrations_Care_Home_Covid_Week_16_Per_100000",
                                                            "Death_Registrations_Hospital_Covid_Week_16",
                                                            "Death_Registrations_Home_Covid_Week_16",
                                                            "Death_Registrations_Home_Covid_Week_16_Per_100000",
                                                            "Death_Registrations_OoHosp_Week_16","Death_Registrations_OoHosp_Week_16_Per_100000",
                                                            "Death_Registrations_OoHosp_Week_16","Death_Registrations_OoHosp_Week_16_Per_100000"))]

covid_death_hosp_data = covid_death_hosp[,-c(1,5:8)]


#standardise everything except outcome
covid_death_hosp_data[,-29] = scale(covid_death_hosp_data[,-29])
summary(covid_death_hosp_data)

tempData_hosp <- mice(covid_death_hosp_data,seed=987)
modelFit_covid_death_hosp = with(tempData_hosp, lm(Death_Registrations_Hospital_Covid_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
                                                     Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+BAME_Including_Jews+Lives_In_Communal+
                                                     Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                                                     Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                                                     Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                                                     Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                                                     Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018
                                                 ))
summary(pool(modelFit_covid_death_hosp))
pool.r.squared(modelFit_covid_death_hosp, adjusted = T)

#post regression diagnostics and sensitivity analyses with robust standard errors
Data_imp = mice::complete(tempData_hosp,1)

modTestRes = lm(Death_Registrations_Hospital_Covid_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
                  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018, Data_imp)
plot(modTestRes)
summary(modTestRes, robust = F)
summary(modTestRes, robust = T)

data = summary(pool(modelFit_covid_death_hosp))
data <- tibble::rownames_to_column(data, "term")

View(data.frame(data$term,data$estimate,(data$estimate-(1.96*data$std.error)),(data$estimate+(1.96*data$std.error)),data$p.value))
dwplot(data)

# Create list of brackets (label, topmost included predictor, bottommost included predictor)
five_brackets <- list(c("Country", "Scotland cf. England", "Wales cf. England"), 
                      c("Weather/Air", "Air Pollution", "Temperature"),
                      c("Demographics/Social", "Age", "Self Funding Care Home"),
                      c("Population Health Rates","Cancer","Ex-Smokers"),
                      c("Other Factors","Google Searching","Port Activity"))
pdf("Model_Hosp_final3.pdf",11,7)
{dwplot(data, 
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2),
        dot_args = list(color = "#F8766D"), # color for the dot
        whisker_args = list(color = "Grey47")   # color for the whisker
) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(Scotland = "Scotland cf. England",                       # relabel predictors
                         Wales = "Wales cf. England",
                         PM_2_5_Total_2018 = "Air Pollution",
                         Average_Humid_23_2_1_3_8_3_15_3 = "Humidity", 
                         Average_tempHighCel_23_2_1_3_8_3_15_3 = "Temperature", 
                         Median_Age_2018 = "Age", 
                         BAME_Including_Jews = "BAME",
                         Lives_In_Communal = "Communal Living",
                         Over_1_5_Per_Room = "Crowded Living",
                         Local_Proportion_Decile_1_Z_Ind_IMD = "Deprivation",
                         Density_People_sq_km = "Population Density",
                         Females_2018 = "Women",
                         Self_Funders_2017 = "Self Funding Care Home",
                         Cancer_Z_Ind = "Cancer",
                         CKD_Z_Ind = "CKD",
                         COPD_Z_Ind = "COPD",
                         Age_Standardised_CVD_Deaths_U75_100000 = "CVD deaths (<75)",
                         Dementia_Z_Ind = "Dementia",
                         Diabetes_Prevalence_Z_Ind = "Diabetes",
                         HTN_Z_Ind = "Hypertension",
                         Obese_Z_Ind = "Obesity",
                         RA_Z_Ind = "Rheumatoid Arthritis",
                         Current_Smoking_Z_Ind = "Smoking (current)",
                         Ex_Smoker_Z_Ind = "Ex-Smokers",
                         Coronavirus_Google_010220_190320 = "Google Searching",
                         Passengers_To_From_Pop_Stand_Normalised = "Flight Passengers",
                         Tonnes_Port_Over_750 = "Port Activity"
    )) +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    ggtitle("Local Authority Covid-19 Deaths per 100,000 in Hospitals") + geom_text(aes(label = "Adjusted R² = 0.50\n(0.43, 0.57)"), x=7.5,y=2, size = 3)+
    theme(plot.title = element_text(face="bold", hjust = 0.5), legend.position = "none",plot.margin = unit(c(0,5,1,1), "cm"))+
    geom_text(x = 11.15,
              y=27, label = "-0.58 (-2.35, 1.19)", #Scotland
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=26, label = "  0.93 (-0.28, 2.14)", #Wales
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=25, label = "  7.22 (4.33, 10.11)   ***", #air
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=24, label = "-0.83 (-2.56, 0.90)", #Humidity
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=23, label = "-3.75 (-6.39, -1.10) **", #temp
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=22, label = "  1.41 (-1.89, 4.70)", #age
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=21, label = "  4.88 (1.41, 8.34)     **", #BAME
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=20, label = "-1.59 (-2.88, -0.30) *", #communal
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=19, label = "-2.99 (-6.03, 0.05)", #crowded
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=18, label = "  2.37 (0.63, 4.11)     **", #depriv
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=17, label = "  0.58 (-1.52, 2.69)", #density
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=16, label = "  1.84 (0.42, 3.26)     *", #women
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=15, label = "-3.10 (-4.76, -1.45) ***", #self funding
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=14, label = "  1.90 (-0.98, 4.78)", #cancer
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=13, label = "-0.26 (-3.07, 2.54)", #ckd
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=12, label = "-1.27 (-3.67, 1.13)", #copd
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=11, label = "-1.54 (-3.80, 0.73)", #cvd
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=10, label = "-2.23 (-4.01, -0.46) *", #dementia
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=9, label = "  1.84 (0.11, 3.58)     *", #diabetes
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=8, label = "  1.39 (-1.28, 4.05)", #BP
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=7, label = "-0.12 (-1.86, 1.62)", #obesity
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=6, label = "-0.79 (-4.05, 2.47)", #RA
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=5, label = "-0.78 (-2.87, 1.30)", #smoking current
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=4, label = "-1.72 (-3.56, 0.12)", #ex-smoking
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=3, label = "-0.12 (-1.40, 1.16)", #google
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=2, label = "  1.94 (0.35, 3.53)     *", #flight passengers
              hjust = 0,
              size = 3) +
    geom_text(x = 11.15,
              y=1, label = "-0.11 (-1.08, 0.87)", #port activity
              hjust = 0,
              size = 3) +
    coord_cartesian(clip = 'off')  # This keeps the labels from disappearing
  } %>% 
  add_brackets(five_brackets)
dev.off()

################################################

covid_death_home = covid_all[ ,!(colnames(covid_all) %in% c("Death_Registrations_Covid_Week_16","Death_Registrations_Care_Home_Covid_Week_16",
                                                            "Death_Registrations_Covid_Week_16_Per_100000","Death_Registrations_Care_Home_Covid_Week_16_Per_100000",
                                                            "Death_Registrations_Hospital_Covid_Week_16", "Death_Registrations_Hospital_Covid_Week_16_Per_100000",
                                                            "Death_Registrations_Home_Covid_Week_16","Death_Registrations_OoHosp_Week_16","Death_Registrations_OoHosp_Week_16_Per_100000"))]

covid_death_home_data = covid_death_home[,-c(1,5:8)]

#standardise everything except outcome
covid_death_home_data[,-29] = scale(covid_death_home_data[,-29])
summary(covid_death_home_data)

tempData_home <- mice(covid_death_home_data,seed=987)
summary(tempData_home)
modelFit_covid_death_home = with(tempData_home, lm(Death_Registrations_Home_Covid_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
                                                     Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+BAME_Including_Jews+Lives_In_Communal+
                                                     Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                                                     Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                                                     Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                                                     Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                                                     Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018
                                                 ))


summary(pool(modelFit_covid_death_home))
pool.r.squared(modelFit_covid_death_home, adjusted = T)

#post regression diagnostics and sensitivity analyses with robust standard errors
Data_imp = mice::complete(tempData_home,1)

modTestRes = lm(Death_Registrations_Home_Covid_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
                  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018, Data_imp)
plot(modTestRes)
summary(modTestRes, robust = F)
summary(modTestRes, robust = T)

data = summary(pool(modelFit_covid_death_home))
data <- tibble::rownames_to_column(data, "term")

View(data.frame(data$term,data$estimate,(data$estimate-(1.96*data$std.error)),(data$estimate+(1.96*data$std.error)),data$p.value))
dwplot(data)

# Create list of brackets (label, topmost included predictor, bottommost included predictor)
five_brackets <- list(c("Country", "Scotland cf. England", "Wales cf. England"), 
                      c("Weather/Air", "Air Pollution", "Temperature"),
                      c("Demographics/Social", "Age", "Self Funding Care Home"),
                      c("Population Health Rates","Cancer","Ex-Smokers"),
                      c("Other Factors","Google Searching","Port Activity"))
pdf("Model_Home_final3.pdf",11,7)
{dwplot(data, 
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2),
        dot_args = list(color = "#F8766D"), # color for the dot
        whisker_args = list(color = "Grey47")   # color for the whisker
) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(Scotland = "Scotland cf. England",                       # relabel predictors
                         Wales = "Wales cf. England",
                         PM_2_5_Total_2018 = "Air Pollution",
                         Average_Humid_23_2_1_3_8_3_15_3 = "Humidity", 
                         Average_tempHighCel_23_2_1_3_8_3_15_3 = "Temperature", 
                         Median_Age_2018 = "Age", 
                         BAME_Including_Jews = "BAME",
                         Lives_In_Communal = "Communal Living",
                         Over_1_5_Per_Room = "Crowded Living",
                         Local_Proportion_Decile_1_Z_Ind_IMD = "Deprivation",
                         Density_People_sq_km = "Population Density",
                         Females_2018 = "Women",
                         Self_Funders_2017 = "Self Funding Care Home",
                         Cancer_Z_Ind = "Cancer",
                         CKD_Z_Ind = "CKD",
                         COPD_Z_Ind = "COPD",
                         Age_Standardised_CVD_Deaths_U75_100000 = "CVD deaths (<75)",
                         Dementia_Z_Ind = "Dementia",
                         Diabetes_Prevalence_Z_Ind = "Diabetes",
                         HTN_Z_Ind = "Hypertension",
                         Obese_Z_Ind = "Obesity",
                         RA_Z_Ind = "Rheumatoid Arthritis",
                         Current_Smoking_Z_Ind = "Smoking (current)",
                         Ex_Smoker_Z_Ind = "Ex-Smokers",
                         Coronavirus_Google_010220_190320 = "Google Searching",
                         Passengers_To_From_Pop_Stand_Normalised = "Flight Passengers",
                         Tonnes_Port_Over_750 = "Port Activity"
    )) +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    ggtitle("Local Authority Covid-19 Deaths per 100,000 at Home") + geom_text(aes(label = "Adjusted R² = 0.36\n(0.28, 0.44)"), x=0.6,y=2, size = 3)+
    theme(plot.title = element_text(face="bold", hjust = 0.5), legend.position = "none",plot.margin = unit(c(0,5,1,1), "cm"))+
    geom_text(x = 0.95,
              y=27, label = "  0.41 (0.17, 0.65)     **", #Scotland
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=26, label = "  0.12 (-0.05, 0.28)", #Wales
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=25, label = "  0.05 (-0.35, 0.44)", #air
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=24, label = "-0.08 (-0.31, 0.16)", #Humidity
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=23, label = "-0.48 (-0.85, -0.12) *", #temp
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=22, label = "  0.27 (-0.18, 0.71)", #age
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=21, label = "  0.21 (-0.27, 0.68)", #BAME
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=20, label = "  0.05 (-0.12, 0.22)", #communal
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=19, label = "  0.44 (0.02, 0.85)     *", #crowded
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=18, label = "  0.33 (0.09, 0.57)     *", #depriv
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=17, label = "  0.24 (-0.05, 0.53)", #density
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=16, label = "  0.31 (0.12, 0.51)     **", #women
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=15, label = "-0.09 (-0.32, 0.14)", #self funding
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=14, label = "  0.02 (-0.38, 0.41)", #cancer
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=13, label = "-0.21 (-0.61, 0.18)", #ckd
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=12, label = "-0.13 (-0.46, 0.20)", #copd
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=11, label = "-0.41 (-0.73, -0.10) *", #cvd
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=10, label = "-0.23 (-0.48, 0.01)", #dementia
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=9, label = "  0.10 (-0.15, 0.34)", #diabetes
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=8, label = "  0.06 (-0.32, 0.44)", #BP
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=7, label = "-0.12 (-0.38, 0.14)", #obesity
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=6, label = "  0.12 (-0.33, 0.57)", #RA
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=5, label = "-0.05 (-0.32, 0.23)", #smoking current
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=4, label = "-0.21 (-0.32, 0.23)", #ex-smoking
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=3, label = "-0.09 (-0.27, 0.08)", #google
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=2, label = "-0.01 (-0.23, 0.21)", #flight passengers
              hjust = 0,
              size = 3) +
    geom_text(x = 0.95,
              y=1, label = "  0.03 (-0.10, 0.17)", #port activity
              hjust = 0,
              size = 3) +
    coord_cartesian(clip = 'off')  # This keeps the labels from disappearing
  } %>% 
  add_brackets(five_brackets)
dev.off()


######################################################################################################
# 
#
# R code for use with the geofacet library
# https://hafen.github.io/geofacet/ 
#


#
# Great Britain
# https://github.com/OrdnanceSurvey/equal-area-cartogram

great_britain = read_csv("lad-cartogram-square-2019-great-britain.csv")

View(great_britain)
covid_deaths = read_csv("covid_deaths_all.csv")

pdf("deathsPer100000LADHeat.pdf",width = 50,height = 75)
  ggplot(covid_deaths, aes(`Location`, `Covid Deaths Per 100,000`, fill = `Covid Deaths Per 100,000`)) +
  geom_col() + scale_fill_viridis() +
  facet_geo(~ name, grid = great_britain, labeller = label_wrap_gen(width=15)) +
  labs(title = "COVID-19 Deaths Per 100,000 People",
       caption = "Data Source: https://www.ons.gov.uk/ & https://www.nrscotland.gov.uk/",
       x = NULL,
       y = "Deaths Per 100,000") +
  theme(strip.text.x = element_text(size = 15), plot.title = element_text(size = 75),
        plot.caption = element_text(size = 30), legend.title = element_blank())
dev.off()

pdf("deathsPer100000ExceptLADHeat.pdf",width = 50,height = 75)
ggplot(covid_deaths, aes(`Location`, `All Deaths Per 100000`, fill = `All Deaths Per 100000`)) +
  geom_col() + scale_fill_viridis() +
  facet_geo(~ name, grid = great_britain, labeller = label_wrap_gen(width=15)) +
  labs(title = "All Deaths (Except COVID-19) Per 100,000 People",
       caption = "Data Source: https://www.ons.gov.uk/ & https://www.nrscotland.gov.uk/",
       x = NULL,
       y = "Deaths Per 100,000") +
  theme(strip.text.x = element_text(size = 15), plot.title = element_text(size = 75),
        plot.caption = element_text(size = 30), legend.title = element_blank())
dev.off()

###################################################################################################
#Week 16 all minus covid

minus_all = read_csv("AllMinusDataNew.csv")

minus_all$Country = as.factor(minus_all$Country)
minus_all$Local_Authority = as.factor(minus_all$Local_Authority)
minus_all$Area_Code = as.factor(minus_all$Area_Code)

#multiple imputation https://www.r-bloggers.com/imputing-missing-data-with-r-mice-package/

minus_death_all = minus_all[ ,!(colnames(minus_all) %in% c("Death_Registrations_Week_16","Death_Registrations_Care_Home_Week_16",
                                                           "Death_Registrations_Care_Home_Week_16_Per_100000",
                                                           "Death_Registrations_Hospital_Week_16",
                                                           "Death_Registrations_Hospital_Week_16_Per_100000",
                                                           "Death_Registrations_Home_Week_16",
                                                           "Death_Registrations_Home_Week_16_Per_100000",
                                                           "Death_Registrations_OoHosp_Week_16","Death_Registrations_OoHosp_Week_16_Per_100000"
))]

minus_death_all_data = minus_death_all[,-c(1,5:8)]


#standardise everything except outcome
minus_death_all_data[,-29] = scale(minus_death_all_data[,-29])
summary(minus_death_all_data)

tempData_all <- mice(minus_death_all_data,seed=987)
summary(tempData_all)

modelFit_minus_death_all = with(tempData_all, lm(Death_Registrations_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
                                                   Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                                                   Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                                                   Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                                                   Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                                                   Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                                                   Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018
))
View(summary(pool(modelFit_minus_death_all)))
pool.r.squared(modelFit_minus_death_all, adjusted = T)

#post regression diagnostics and sensitivity analyses with robust standard errors
Data_imp = mice::complete(tempData_all,1)

modTestRes = lm(Death_Registrations_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
                  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018, Data_imp)
plot(modTestRes)
summary(modTestRes, robust = F)
summary(modTestRes, robust = T)

data = summary(pool(modelFit_minus_death_all))
data <- tibble::rownames_to_column(data, "term")

dwplot(data)
#look at 95CIs
View(data.frame(data$term,data$estimate,(data$estimate-(1.96*data$std.error)),(data$estimate+(1.96*data$std.error)),data$p.value))

# Create list of brackets (label, topmost included predictor, bottommost included predictor)
five_brackets <- list(c("Country", "Scotland cf. England", "Wales cf. England"), 
                      c("Weather/Air", "Air Pollution", "Temperature"),
                      c("Demographics/Social", "Age", "Self Funding Care Home"),
                      c("Population Health Rates","Cancer","Ex-Smokers"),
                      c("Other Factors","Google Searching","Port Activity"))
pdf("Model_All_Minus_final3.pdf",11,7)
{dwplot(data, 
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2),
        dot_args = list(color = "#F8766D"), # color for the dot
        whisker_args = list(color = "Grey47")   # color for the whisker
) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(Scotland = "Scotland cf. England",                       # relabel predictors
                         Wales = "Wales cf. England",
                         PM_2_5_Total_2018 = "Air Pollution",
                         Average_Humid_23_2_1_3_8_3_15_3 = "Humidity", 
                         Average_tempHighCel_23_2_1_3_8_3_15_3 = "Temperature", 
                         Median_Age_2018 = "Age", 
                         BAME_Including_Jews = "BAME",
                         Lives_In_Communal = "Communal Living",
                         Over_1_5_Per_Room = "Crowded Living",
                         Local_Proportion_Decile_1_Z_Ind_IMD = "Deprivation",
                         Density_People_sq_km = "Population Density",
                         Females_2018 = "Women",
                         Self_Funders_2017 = "Self Funding Care Home",
                         Cancer_Z_Ind = "Cancer",
                         CKD_Z_Ind = "CKD",
                         COPD_Z_Ind = "COPD",
                         Age_Standardised_CVD_Deaths_U75_100000 = "CVD deaths (<75)",
                         Dementia_Z_Ind = "Dementia",
                         Diabetes_Prevalence_Z_Ind = "Diabetes",
                         HTN_Z_Ind = "Hypertension",
                         Obese_Z_Ind = "Obesity",
                         RA_Z_Ind = "Rheumatoid Arthritis",
                         Current_Smoking_Z_Ind = "Smoking (current)",
                         Ex_Smoker_Z_Ind = "Ex-Smokers",
                         Coronavirus_Google_010220_190320 = "Google Searching",
                         Passengers_To_From_Pop_Stand_Normalised = "Flight Passengers",
                         Tonnes_Port_Over_750 = "Port Activity"
    )) +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    ggtitle("Local Authority All Deaths (except COVID-19) per 100,000 in All Settings") + geom_text(aes(label = "Adjusted R² = 0.87\n(0.84, 0.89)"), x=34,y=2, size = 3)+
    theme(plot.title = element_text(face="bold", hjust = 0.5), legend.position = "none",plot.margin = unit(c(0,5,1,1), "cm"))+
    geom_text(x = 45,
              y=27, label = "  8.54 (3.43, 13.66)    **", #Scotland
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=26, label = "  7.75 (4.26, 11.23)    ***", #Wales
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=25, label = "  13.16 (4.83, 21.48)  **", #air
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=24, label = "-2.66 (-7.58, 2.26)", #Humidity
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=23, label = "-5.48 (-13.20, 2.25)", #temp
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=22, label = "  32.10 (22.59, 41.62)***", #age
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=21, label = "  5.32 (-4.52, 15.16)", #BAME
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=20, label = "  9.35 (5.83, 12.86)    ***", #communal
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=19, label = "-6.41 (-15.17, 2.35)", #crowded
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=18, label = "  9.41 (4.40, 14.42)    ***", #depriv
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=17, label = "-2.04 (-8.08, 3.99)", #density
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=16, label = "  7.42 (3.40, 11.43)    ***", #women
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=15, label = "-4.30 (-9.07, 0.46)", #self funding
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=14, label = "-7.44 (-15.95, 1.07)", #cancer
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=13, label = "  19.96 (11.87, 28.05)***", #ckd
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=12, label = "  8.73 (1.69, 15.77)    *", #copd
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=11, label = "  12.72 (6.27, 19.17)  ***", #cvd
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=10, label = "  9.10 (3.98, 14.21)    **", #dementia
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=9, label = "-1.23 (-6.22, 3.75)", #diabetes
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=8, label = "  1.65 (-5.98, 9.29)", #BP
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=7, label = "-2.45 (-7.32, 2.42)", #obesity
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=6, label = "  8.94 (-0.35, 18.23)", #RA
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=5, label = "  5.29 (-1.33, 11.90)", #smoking current
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=4, label = "  2.74 (-2.68, 8.16)", #ex-smoking
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=3, label = "-0.35 (-4.06, 3.36)", #google
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=2, label = "-6.92 (-11.48, -2.36)**", #flight passengers
              hjust = 0,
              size = 3) +
    geom_text(x = 45,
              y=1, label = "-0.49 (-3.28, 2.30)", #port activity
              hjust = 0,
              size = 3) +
    coord_cartesian(clip = 'off')  # This keeps the labels from disappearing
  
} %>% 
  add_brackets(five_brackets)
dev.off()


################################################################################

minus_care_all = minus_all[ ,!(colnames(minus_all) %in% c("Death_Registrations_Week_16","Death_Registrations_Week_16_Per_100000",
                                                          "Death_Registrations_Care_Home_Week_16",
                                                          "Death_Registrations_Hospital_Week_16",
                                                           "Death_Registrations_Hospital_Week_16_Per_100000",
                                                           "Death_Registrations_Home_Week_16",
                                                           "Death_Registrations_Home_Week_16_Per_100000",
                                                           "Death_Registrations_OoHosp_Week_16","Death_Registrations_OoHosp_Week_16_Per_100000"
))]

minus_care_all_data = minus_care_all[,-c(1,5:8)]

#standardise everything except outcome
minus_care_all_data[,-29] = scale(minus_care_all_data[,-29])
summary(minus_care_all_data)

tempData_all <- mice(minus_care_all_data,seed=987)
summary(tempData_all)


modelFit_minus_care_all = with(tempData_all, lm(Death_Registrations_Care_Home_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
                                                   Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                                                   Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                                                   Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                                                   Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                                                   Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                                                   Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018
))
summary(pool(modelFit_minus_care_all))
pool.r.squared(modelFit_minus_care_all, adjusted = T)

#post regression diagnostics and sensitivity analyses with robust standard errors
Data_imp = mice::complete(tempData_all,1)

modTestRes = lm(Death_Registrations_Care_Home_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
                  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018, Data_imp)
plot(modTestRes)
summary(modTestRes, robust = F)
summary(modTestRes, robust = T)

data = summary(pool(modelFit_minus_care_all))
data <- tibble::rownames_to_column(data, "term")

dwplot(data)
#look at 95CIs
View(data.frame(data$term,data$estimate,(data$estimate-(1.96*data$std.error)),(data$estimate+(1.96*data$std.error)),data$p.value))

# Create list of brackets (label, topmost included predictor, bottommost included predictor)
five_brackets <- list(c("Country", "Scotland cf. England", "Wales cf. England"), 
                      c("Weather/Air", "Air Pollution", "Temperature"),
                      c("Demographics/Social", "Age", "Self Funding Care Home"),
                      c("Population Health Rates","Cancer","Ex-Smokers"),
                      c("Other Factors","Google Searching","Port Activity"))
pdf("Model_Care_Minus_final3.pdf",11,7)
{dwplot(data, 
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2),
        dot_args = list(color = "#F8766D"), # color for the dot
        whisker_args = list(color = "Grey47")   # color for the whisker
) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(Scotland = "Scotland cf. England",                       # relabel predictors
                         Wales = "Wales cf. England",
                         PM_2_5_Total_2018 = "Air Pollution",
                         Average_Humid_23_2_1_3_8_3_15_3 = "Humidity", 
                         Average_tempHighCel_23_2_1_3_8_3_15_3 = "Temperature", 
                         Median_Age_2018 = "Age", 
                         BAME_Including_Jews = "BAME",
                         Lives_In_Communal = "Communal Living",
                         Over_1_5_Per_Room = "Crowded Living",
                         Local_Proportion_Decile_1_Z_Ind_IMD = "Deprivation",
                         Density_People_sq_km = "Population Density",
                         Females_2018 = "Women",
                         Self_Funders_2017 = "Self Funding Care Home",
                         Cancer_Z_Ind = "Cancer",
                         CKD_Z_Ind = "CKD",
                         COPD_Z_Ind = "COPD",
                         Age_Standardised_CVD_Deaths_U75_100000 = "CVD deaths (<75)",
                         Dementia_Z_Ind = "Dementia",
                         Diabetes_Prevalence_Z_Ind = "Diabetes",
                         HTN_Z_Ind = "Hypertension",
                         Obese_Z_Ind = "Obesity",
                         RA_Z_Ind = "Rheumatoid Arthritis",
                         Current_Smoking_Z_Ind = "Smoking (current)",
                         Ex_Smoker_Z_Ind = "Ex-Smokers",
                         Coronavirus_Google_010220_190320 = "Google Searching",
                         Passengers_To_From_Pop_Stand_Normalised = "Flight Passengers",
                         Tonnes_Port_Over_750 = "Port Activity"
    )) +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    ggtitle("Local Authority All Deaths (except COVID-19) per 100,000 in Care Homes") + geom_text(aes(label = "Adjusted R² = 0.66\n(0.60, 0.71)"), x=15,y=2, size = 3)+
    theme(plot.title = element_text(face="bold", hjust = 0.5), legend.position = "none",plot.margin = unit(c(0,5,1,1), "cm"))+
    geom_text(x = 20.5,
              y=27, label = "  0.17 (-3.12, 3.46)", #Scotland
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=26, label = "-4.78 (-7.05, -2.52) ***", #Wales
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=25, label = "  0.08 (-5.37, 5.52)", #air
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=24, label = "  1.02 (-2.23, 4.27)", #Humidity
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=23, label = "  2.12 (-2.90, 7.14)", #temp
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=22, label = "  12.49 (6.28, 18.71) ***", #age
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=21, label = "  6.15 (-0.28, 12.58)", #BAME
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=20, label = "  6.07 (3.75, 8.40)     ***", #communal
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=19, label = "  1.21 (-4.52, 6.94)", #crowded
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=18, label = "  2.50 (-0.77, 5.77)", #depriv
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=17, label = "-1.12 (-5.08, 2.84)", #density
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=16, label = "  5.70 (3.04, 8.35)     ***", #women
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=15, label = "-1.08 (-4.19, 2.03)", #self funding
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=14, label = "-5.02 (-10.44, 0.39)", #cancer
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=13, label = "  13.05 (7.66, 18.43) ***", #ckd
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=12, label = "-1.86 (-6.48, 2.76)", #copd
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=11, label = "  3.98 (-0.17, 8.14)", #cvd
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=10, label = "  7.45 (4.07, 10.82)   ***", #dementia
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=9, label = "-4.76 (-8.06, -1.46) **", #diabetes
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=8, label = "-4.05 (-9.20, 1.10)", #BP
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=7, label = "  2.41 (-0.88, 5.70)", #obesity
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=6, label = "  4.55 (-1.63, 10.72)", #RA
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=5, label = "  0.92 (-3.02, 4.86)", #smoking current
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=4, label = "  2.56 (-1.03, 6.14)", #ex-smoking
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=3, label = "  0.52 (-1.90, 2.93)", #google
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=2, label = "-1.34 (-4.34, 1.67)", #flight passengers
              hjust = 0,
              size = 3) +
    geom_text(x = 20.5,
              y=1, label = "  0.33 (-1.52, 2.18)", #port activity
              hjust = 0,
              size = 3) +

    coord_cartesian(clip = 'off')  # This keeps the labels from disappearing
  
} %>% 
  add_brackets(five_brackets)
dev.off()


################################################

minus_hosp_all = minus_all[ ,!(colnames(minus_all) %in% c("Death_Registrations_Week_16","Death_Registrations_Week_16_Per_100000",
                                                          "Death_Registrations_Care_Home_Week_16","Death_Registrations_Care_Home_Week_16_Per_100000",
                                                          "Death_Registrations_Hospital_Week_16",
                                                          "Death_Registrations_Home_Week_16",
                                                          "Death_Registrations_Home_Week_16_Per_100000",
                                                          "Death_Registrations_OoHosp_Week_16","Death_Registrations_OoHosp_Week_16_Per_100000"
))]

minus_hosp_all_data = minus_hosp_all[,-c(1,5:8)]

#standardise everything except outcome
minus_hosp_all_data[,-29] = scale(minus_hosp_all_data[,-29])
summary(minus_hosp_all_data)

tempData_all <- mice(minus_hosp_all_data,seed=987)
summary(tempData_all)

modelFit_minus_hosp_all = with(tempData_all, lm(Death_Registrations_Hospital_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
                                                  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                                                  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                                                  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                                                  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                                                  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                                                  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018
))
summary(pool(modelFit_minus_hosp_all))
pool.r.squared(modelFit_minus_hosp_all, adjusted = T)

#post regression diagnostics and sensitivity analyses with robust standard errors
Data_imp = mice::complete(tempData_all,1)

modTestRes = lm(Death_Registrations_Hospital_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
                  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018, Data_imp)
plot(modTestRes)
summary(modTestRes, robust = F)
summary(modTestRes, robust = T)

data = summary(pool(modelFit_minus_hosp_all))
data <- tibble::rownames_to_column(data, "term")

dwplot(data)
#look at 95CIs
View(data.frame(data$term,data$estimate,(data$estimate-(1.96*data$std.error)),(data$estimate+(1.96*data$std.error)),data$p.value))

# Create list of brackets (label, topmost included predictor, bottommost included predictor)
five_brackets <- list(c("Country", "Scotland cf. England", "Wales cf. England"), 
                      c("Weather/Air", "Air Pollution", "Temperature"),
                      c("Demographics/Social", "Age", "Self Funding Care Home"),
                      c("Population Health Rates","Cancer","Ex-Smokers"),
                      c("Other Factors","Google Searching","Port Activity"))
pdf("Model_Hosp_Minus_final3.pdf",11,7)
{dwplot(data, 
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2),
        dot_args = list(color = "#F8766D"), # color for the dot
        whisker_args = list(color = "Grey47")   # color for the whisker
) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(Scotland = "Scotland cf. England",                       # relabel predictors
                         Wales = "Wales cf. England",
                         PM_2_5_Total_2018 = "Air Pollution",
                         Average_Humid_23_2_1_3_8_3_15_3 = "Humidity", 
                         Average_tempHighCel_23_2_1_3_8_3_15_3 = "Temperature", 
                         Median_Age_2018 = "Age", 
                         BAME_Including_Jews = "BAME",
                         Lives_In_Communal = "Communal Living",
                         Over_1_5_Per_Room = "Crowded Living",
                         Local_Proportion_Decile_1_Z_Ind_IMD = "Deprivation",
                         Density_People_sq_km = "Population Density",
                         Females_2018 = "Women",
                         Self_Funders_2017 = "Self Funding Care Home",
                         Cancer_Z_Ind = "Cancer",
                         CKD_Z_Ind = "CKD",
                         COPD_Z_Ind = "COPD",
                         Age_Standardised_CVD_Deaths_U75_100000 = "CVD deaths (<75)",
                         Dementia_Z_Ind = "Dementia",
                         Diabetes_Prevalence_Z_Ind = "Diabetes",
                         HTN_Z_Ind = "Hypertension",
                         Obese_Z_Ind = "Obesity",
                         RA_Z_Ind = "Rheumatoid Arthritis",
                         Current_Smoking_Z_Ind = "Smoking (current)",
                         Ex_Smoker_Z_Ind = "Ex-Smokers",
                         Coronavirus_Google_010220_190320 = "Google Searching",
                         Passengers_To_From_Pop_Stand_Normalised = "Flight Passengers",
                         Tonnes_Port_Over_750 = "Port Activity"
    )) +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    ggtitle("Local Authority All Deaths (except COVID-19) per 100,000 in Hospitals") + geom_text(aes(label = "Adjusted R² = 0.73\n(0.68, 0.77)"), x=15,y=2, size = 3)+
    theme(plot.title = element_text(face="bold", hjust = 0.5), legend.position = "none",plot.margin = unit(c(0,5,1,1), "cm"))+
    geom_text(x = 19,
              y=27, label = "  5.13 (1.71, 8.56)      **", #Scotland
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=26, label = "  10.95 (8.64, 13.27)  ***", #Wales
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=25, label = "  8.56 (3.02, 14.10)    **", #ayr
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=24, label = "-4.10 (-7.38, -0.82)  *", #Humidity
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=23, label = "-5.88 (-10.91, -0.86)*", #temp
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=22, label = "  8.81 (2.51, 15.11)    **", #age
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=21, label = "-0.22 (-6.81, 6.38)", #BAME
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=20, label = "-0.37 (-2.71, 1.98)", #communal
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=19, label = "-8.11 (-13.92, -2.30)**", #crowded
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=18, label = "  2.07 (-1.25, 5.39)", #depriv
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=17, label = "-0.93 (-4.97, 3.12)", #density
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=16, label = "-0.62 (-3.29, 2.06)", #women
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=15, label = "-2.27 (-5.44, 0.89)", #self funding
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=14, label = "  3.34 (-2.26, 8.94)", #cancer
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=13, label = "  11.80 (6.39, 17.21)  ***", #ckd
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=12, label = "  6.15 (1.43, 10.87)    *", #copd
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=11, label = "  5.57 (1.37, 9.77)      *", #cvd
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=10, label = "  0.03 (-3.34, 3.41)", #dementia
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=9, label = "  0.90 (-2.55, 4.34)", #diabetes
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=8, label = "-0.75 (-5.92, 4.42)", #BP
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=7, label = "  0.40 (-2.84, 3.64)", #obesity
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=6, label = "-2.13 (-8.42, 4.16)", #RA
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=5, label = "  6.00 (1.96, 10.05)    **", #smoking current
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=4, label = "-3.38 (-6.90, 0.15)", #ex-smoking
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=3, label = "-0.68 (-3.13, 1.77)", #google
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=2, label = "-0.50 (-2.55, 3.55)", #flight passengers
              hjust = 0,
              size = 3) +
    geom_text(x = 19,
              y=1, label = "-0.97 (-2.84, 0.90)", #port activity
              hjust = 0,
              size = 3) +
    
    coord_cartesian(clip = 'off')  # This keeps the labels from disappearing
  
} %>% 
  add_brackets(five_brackets)
dev.off()

######################################################

minus_home_all = minus_all[ ,!(colnames(minus_all) %in% c("Death_Registrations_Week_16","Death_Registrations_Week_16_Per_100000",
                                                          "Death_Registrations_Care_Home_Week_16","Death_Registrations_Care_Home_Week_16_Per_100000",
                                                          "Death_Registrations_Hospital_Week_16", "Death_Registrations_Hospital_Week_16_Per_100000",
                                                          "Death_Registrations_Home_Week_16",
                                                          "Death_Registrations_OoHosp_Week_16","Death_Registrations_OoHosp_Week_16_Per_100000"
))]

minus_home_all_data = minus_home_all[,-c(1,5:8)]

#standardise everything except outcome
minus_home_all_data[,-29] = scale(minus_home_all_data[,-29])
summary(minus_home_all_data)

tempData_all <- mice(minus_home_all_data,seed=987)
summary(tempData_all)

modelFit_minus_home_all = with(tempData_all, lm(Death_Registrations_Home_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
                                                  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                                                  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                                                  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                                                  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                                                  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                                                  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018
))
summary(pool(modelFit_minus_home_all))
pool.r.squared(modelFit_minus_home_all, adjusted = T)

#post regression diagnostics and sensitivity analyses with robust standard errors
Data_imp = mice::complete(tempData_all,1)

modTestRes = lm(Death_Registrations_Home_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
                  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018, Data_imp)
plot(modTestRes)
summary(modTestRes, robust = F)
summary(modTestRes, robust = T)

data = summary(pool(modelFit_minus_home_all))
data <- tibble::rownames_to_column(data, "term")

dwplot(data)
#look at 95CIs
View(data.frame(data$term,data$estimate,(data$estimate-(1.96*data$std.error)),(data$estimate+(1.96*data$std.error)),data$p.value))

# Create list of brackets (label, topmost included predictor, bottommost included predictor)
five_brackets <- list(c("Country", "Scotland cf. England", "Wales cf. England"), 
                      c("Weather/Air", "Air Pollution", "Temperature"),
                      c("Demographics/Social", "Age", "Self Funding Care Home"),
                      c("Population Health Rates","Cancer","Ex-Smokers"),
                      c("Other Factors","Google Searching","Port Activity"))
pdf("Model_Home_Minus_final3.pdf",11,7)
{dwplot(data, 
        vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2),
        dot_args = list(color = "#F8766D"), # color for the dot
        whisker_args = list(color = "Grey47")   # color for the whisker
) %>% # plot line at zero _behind_ coefs
    relabel_predictors(c(Scotland = "Scotland cf. England",                       # relabel predictors
                         Wales = "Wales cf. England",
                         PM_2_5_Total_2018 = "Air Pollution",
                         Average_Humid_23_2_1_3_8_3_15_3 = "Humidity", 
                         Average_tempHighCel_23_2_1_3_8_3_15_3 = "Temperature", 
                         Median_Age_2018 = "Age", 
                         BAME_Including_Jews = "BAME",
                         Lives_In_Communal = "Communal Living",
                         Over_1_5_Per_Room = "Crowded Living",
                         Local_Proportion_Decile_1_Z_Ind_IMD = "Deprivation",
                         Density_People_sq_km = "Population Density",
                         Females_2018 = "Women",
                         Self_Funders_2017 = "Self Funding Care Home",
                         Cancer_Z_Ind = "Cancer",
                         CKD_Z_Ind = "CKD",
                         COPD_Z_Ind = "COPD",
                         Age_Standardised_CVD_Deaths_U75_100000 = "CVD deaths (<75)",
                         Dementia_Z_Ind = "Dementia",
                         Diabetes_Prevalence_Z_Ind = "Diabetes",
                         HTN_Z_Ind = "Hypertension",
                         Obese_Z_Ind = "Obesity",
                         RA_Z_Ind = "Rheumatoid Arthritis",
                         Current_Smoking_Z_Ind = "Smoking (current)",
                         Ex_Smoker_Z_Ind = "Ex-Smokers",
                         Coronavirus_Google_010220_190320 = "Google Searching",
                         Passengers_To_From_Pop_Stand_Normalised = "Flight Passengers",
                         Tonnes_Port_Over_750 = "Port Activity"
    )) +
    theme_bw() + xlab("Coefficient Estimate") + ylab("") +
    ggtitle("Local Authority All Deaths (except COVID-19) per 100,000 At Home") + geom_text(aes(label = "Adjusted R² = 0.65\n(0.59, 0.71)"), x=10.5,y=2, size = 3)+
    theme(plot.title = element_text(face="bold", hjust = 0.5), legend.position = "none",plot.margin = unit(c(0,5,1,1), "cm"))+
    geom_text(x = 16,
              y=27, label = "  8.27 (5.12, 11.49)     ***", #Scotland
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=26, label = "  2.44 (0.34, 4.54)       *", #Wales
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=25, label = "-0.52 (-5.55, 4.51)", #air
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=24, label = "  3.54 (0.58, 6.50)       *", #Humidity
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=23, label = "  1.92 (-2.67, 6.50)", #temp
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=22, label = "  8.76 (3.05, 14.48)     **", #age
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=21, label = "  0.15 (-5.66, 5.96)", #BAME
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=20, label = "  2.22 (0.09, 4.35)       *", #communal
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=19, label = "-1.40 (-6.52, 3.72)", #crowded
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=18, label = "  4.38 (1.43, 7.32)       **", #depriv
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=17, label = "  2.54 (-1.04, 6.11)", #density
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=16, label = "  1.49 (-0.91, 3.89)", #women
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=15, label = "-2.33 (-5.13, 0.47)", #self funding
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=14, label = "-8.17 (-13.17, -3.17) **", #cancer
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=13, label = "-4.51 (-9.32, 0.30)", #ckd
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=12, label = "  3.43 (-0.62, 7.49)", #copd
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=11, label = "  1.00 (-3.06, 5.05)", #cvd
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=10, label = "  0.13 (-2.93, 3.19)", #dementia
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=9, label = "  2.18 (-0.87, 5.22)", #diabetes
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=8, label = "  7.78 (2.84, 12.71)     **", #BP
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=7, label = "-2.21 (-5.34, 0.92)", #obesity
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=6, label = "  7.15 (1.47, 12.83)     *", #RA
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=5, label = "-2.60 (-6.25, 1.06)", #smoking current
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=4, label = "  0.90 (-2.18, 3.98)", #ex-smoking
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=3, label = "-0.27 (-2.45, 1.92)", #google
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=2, label = "-5.28 (-8.03, -2.53)   ***", #flight passengers
              hjust = 0,
              size = 3) +
    geom_text(x = 16,
              y=1, label = "  0.41 (-1.25, 2.06)", #port activity
              hjust = 0,
              size = 3) +
    coord_cartesian(clip = 'off')  # This keeps the labels from disappearing
  
} %>% 
  add_brackets(five_brackets)
dev.off()

#######################################################

#Multivariate multiple regression all settings covid, vs all settings except covid

both_death_all_data = covid_death_all[,-c(1,5:8)]

#standardise everything except outcome
both_death_all_data[,-29] = scale(both_death_all_data[,-29])

both_death_all_data$Death_Registrations_Week_16_Per_100000 = minus_death_all$Death_Registrations_Week_16_Per_100000
summary(both_death_all_data)

tempData_all <- mice(both_death_all_data,seed=987)

#multiple mulvariate regression doesn't work with multiple datasets
Data_imp = mice::complete(tempData_all,1)

modelFit_both_death_all = lm(cbind(Death_Registrations_Covid_Week_16_Per_100000, Death_Registrations_Week_16_Per_100000)~ Scotland+Wales+Density_People_sq_km+
                                                   Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                                                   Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                                                   Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                                                   Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                                                   Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                                                   Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018,data = Data_imp)
sum = summary(modelFit_both_death_all)
data_cov = sum$`Response Death_Registrations_Covid_Week_16_Per_100000`$coefficients
data_all = sum$`Response Death_Registrations_Week_16_Per_100000`$coefficients
#look at 95CIs
View(data.frame("est"=data_cov[,1],"low"=(data_cov[,1]-(1.96*data_cov[,2])),"high"=(data_cov[,1]+(1.96*data_cov[,2])),"p"=data_cov[,4]))

View(data.frame("est"=data_all[,1],"low"=(data_all[,1]-(1.96*data_all[,2])),"high"=(data_all[,1]+(1.96*data_all[,2])),"p"=data_all[,4]))

#Type II Manova - variables still different from null given both dependent outcomes
Anova(modelFit_both_death_all)

#scale outcomes for walds test of coefficient equality
Data_imp =  as.data.frame(scale(Data_imp))
summary(Data_imp)

#significant Coefficients from Manova, equality test across seemingly unrelated regressions
#https://andrewpwheeler.com/2017/06/12/testing-the-equality-of-coefficients-same-independent-different-dependent-variables/?unapproved=9632&moderation-hash=0c84f568615a42d280c784fc5a66379c#comment-9632
covidAllmod = Death_Registrations_Covid_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018

allAllmod = Death_Registrations_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018

fitsur <- systemfit(list(covidAllmod = covidAllmod, allAllmod = allAllmod), data=Data_imp, method="SUR")
summary(fitsur)

#testing whether income effect is equivalent for both models
Covid_more_scot <- glht(fitsur,linfct = c("covidAllmod_Scotland - allAllmod_Scotland= 0"))
summary(Covid_more_scot) # p = 0.296
Covid_more_wales <- glht(fitsur,linfct = c("covidAllmod_Wales - allAllmod_Wales= 0"))
summary(Covid_more_wales) # p = 0.274
Covid_more_air <- glht(fitsur,linfct = c("covidAllmod_PM_2_5_Total_2018 - allAllmod_PM_2_5_Total_2018= 0"))
summary(Covid_more_air) # p = 0.0067
Covid_more_temp <- glht(fitsur,linfct = c("covidAllmod_Average_tempHighCel_23_2_1_3_8_3_15_3 - allAllmod_Average_tempHighCel_23_2_1_3_8_3_15_3= 0"))
summary(Covid_more_temp) # p = 0.00358
Covid_more_age <- glht(fitsur,linfct = c("covidAllmod_Median_Age_2018 - allAllmod_Median_Age_2018= 0"))
summary(Covid_more_age) #p = 0.0067
Covid_more_comm <- glht(fitsur,linfct = c("covidAllmod_Lives_In_Communal - allAllmod_Lives_In_Communal= 0"))
summary(Covid_more_comm) # p = 0.00000396
Covid_more_dep <- glht(fitsur,linfct = c("covidAllmod_Local_Proportion_Decile_1_Z_Ind_IMD - allAllmod_Local_Proportion_Decile_1_Z_Ind_IMD = 0"))
summary(Covid_more_dep) # p = 0.0682
Covid_more_fem <- glht(fitsur,linfct = c("covidAllmod_Females_2018 - allAllmod_Females_2018 = 0"))
summary(Covid_more_fem) # p = 0.0355
Covid_more_ckd <- glht(fitsur,linfct = c("covidAllmod_CKD_Z_Ind - allAllmod_CKD_Z_Ind = 0"))
summary(Covid_more_ckd) # p = 0.017
Covid_more_copd <- glht(fitsur,linfct = c("covidAllmod_COPD_Z_Ind - allAllmod_COPD_Z_Ind = 0"))
summary(Covid_more_copd) # p =  0.00321
Covid_more_cvd <- glht(fitsur,linfct = c("covidAllmod_Age_Standardised_CVD_Deaths_U75_100000 - allAllmod_Age_Standardised_CVD_Deaths_U75_100000 = 0"))
summary(Covid_more_cvd) # 0.000061
Covid_more_dem <- glht(fitsur,linfct = c("covidAllmod_Dementia_Z_Ind - allAllmod_Dementia_Z_Ind = 0"))
summary(Covid_more_dem) # 0.00000329
Covid_more_exsm <- glht(fitsur,linfct = c("covidAllmod_Ex_Smoker_Z_Ind - allAllmod_Ex_Smoker_Z_Ind = 0"))
summary(Covid_more_exsm) # 0.00159
Covid_more_flight <- glht(fitsur,linfct = c("covidAllmod_Passengers_To_From_Pop_Stand_Normalised - allAllmod_Passengers_To_From_Pop_Stand_Normalised = 0"))
summary(Covid_more_flight) # 0.000836
Covid_more_self <- glht(fitsur,linfct = c("covidAllmod_Self_Funders_2017 - allAllmod_Self_Funders_2017 = 0"))
summary(Covid_more_self) # 0.0972

View(p.adjust(p = c(0.296,0.274,0.0067,0.00358,0.0067,0.00000396,0.0682,0.0355,0.017,0.00321,0.000061,0.00000329,0.00159,0.000836,0.0972), method = "fdr"))

############################################################
#Multivariate multiple regression care home setting covid, vs care home setting except covid

both_death_care_data = covid_death_care[,-c(1,5:8)]

#standardise everything except outcome
both_death_care_data[,-29] = scale(both_death_care_data[,-29])

both_death_care_data$Death_Registrations_Care_Home_Week_16_Per_100000 = minus_care_all$Death_Registrations_Care_Home_Week_16_Per_100000
summary(both_death_care_data)

tempData_care <- mice(both_death_care_data,seed=987)

#multiple mulvariate regression doesn't work with multiple datasets
Data_imp = mice::complete(tempData_care,1)

modelFit_both_death_care = lm(cbind(Death_Registrations_Care_Home_Covid_Week_16_Per_100000, Death_Registrations_Care_Home_Week_16_Per_100000)~ Scotland+Wales+Density_People_sq_km+
                               Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                               Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                               Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                               Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                               Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                               Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018,data = Data_imp)
summary(modelFit_both_death_care)

#Type II Manova - variables still different from null given both dependent outcomes
Anova(modelFit_both_death_care)

#scale outcomes for walds test of coefficient equality
Data_imp =  as.data.frame(scale(Data_imp))
summary(Data_imp)

#significant Coefficients from Manova, equality test across seemingly unrelated regressions
#https://andrewpwheeler.com/2017/06/12/testing-the-equality-of-coefficients-same-independent-different-dependent-variables/?unapproved=9632&moderation-hash=0c84f568615a42d280c784fc5a66379c#comment-9632
covidCaremod = Death_Registrations_Care_Home_Covid_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018

allCaremod = Death_Registrations_Care_Home_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018

fitsur2 <- systemfit(list(covidCaremod = covidCaremod, allCaremod = allCaremod), data=Data_imp, method="SUR")
summary(fitsur2)

#testing whether effect is equivalent for both models
Covid_more_2_scot <- glht(fitsur2,linfct = c("covidCaremod_Scotland - allCaremod_Scotland= 0"))
summary(Covid_more_2_scot) #p = 0.00765
Covid_more_2_wales <- glht(fitsur2,linfct = c("covidCaremod_Wales - allCaremod_Wales= 0"))
summary(Covid_more_2_wales) #p = 0.0277
Covid_more_2_temp <- glht(fitsur2,linfct = c("covidCaremod_Average_tempHighCel_23_2_1_3_8_3_15_3 - allCaremod_Average_tempHighCel_23_2_1_3_8_3_15_3= 0"))
summary(Covid_more_2_temp) #p = 0.0067
Covid_more_2_age <- glht(fitsur2,linfct = c("covidCaremod_Median_Age_2018 - allCaremod_Median_Age_2018= 0"))
summary(Covid_more_2_age) #p = 0.00526
Covid_more_2_bame <- glht(fitsur2,linfct = c("covidCaremod_BAME_Including_Jews - allCaremod_BAME_Including_Jews = 0"))
summary(Covid_more_2_bame) #p = 0.00487
Covid_more_2_comm <- glht(fitsur2,linfct = c("covidCaremod_Lives_In_Communal - allCaremod_Lives_In_Communal= 0"))
summary(Covid_more_2_comm) #p = 0.0000211
Covid_more_2_dep <- glht(fitsur2,linfct = c("covidCaremod_Local_Proportion_Decile_1_Z_Ind_IMD - allCaremod_Local_Proportion_Decile_1_Z_Ind_IMD = 0"))
summary(Covid_more_2_dep) #p =   0.146
Covid_more_2_fem <- glht(fitsur2,linfct = c("covidCaremod_Females_2018 - allCaremod_Females_2018 = 0"))
summary(Covid_more_2_fem) #p = 0.756
Covid_more_2_ckd <- glht(fitsur2,linfct = c("covidCaremod_CKD_Z_Ind - allCaremod_CKD_Z_Ind = 0"))
summary(Covid_more_2_ckd) #p = 0.0264
Covid_more_2_dem <- glht(fitsur2,linfct = c("covidCaremod_Dementia_Z_Ind - allCaremod_Dementia_Z_Ind = 0"))
summary(Covid_more_2_dem) #p = 0.000266
Covid_more_2_diab <- glht(fitsur2,linfct = c("covidCaremod_Diabetes_Prevalence_Z_Ind - allCaremod_Diabetes_Prevalence_Z_Ind = 0"))
summary(Covid_more_2_diab) #p = 0.81
Covid_more_2_exsm <- glht(fitsur2,linfct = c("covidCaremod_Ex_Smoker_Z_Ind - allCaremod_Ex_Smoker_Z_Ind = 0"))
summary(Covid_more_2_exsm) #p = 0.0000142

View(p.adjust(p = c(0.00765,0.0277,0.0067,0.00526,0.00487,0.0000211,0.146,0.756,0.0264,0.000266,0.81,0.0000142), method = "fdr"))

############################################################
#Multivariate multiple regression hospital setting covid, vs hospital setting except covid

both_death_hosp_data = covid_death_hosp[,-c(1,5:8)]

#standardise everything except outcome
both_death_hosp_data[,-29] = scale(both_death_hosp_data[,-29])

both_death_hosp_data$Death_Registrations_Hospital_Week_16_Per_100000 = minus_hosp_all$Death_Registrations_Hospital_Week_16_Per_100000
summary(both_death_hosp_data)

tempData_hosp <- mice(both_death_hosp_data,seed=987)

#multiple mulvariate regression doesn't work with multiple datasets
Data_imp = mice::complete(tempData_hosp,1)

modelFit_both_death_hosp = lm(cbind(Death_Registrations_Hospital_Covid_Week_16_Per_100000, Death_Registrations_Hospital_Week_16_Per_100000)~ Scotland+Wales+Density_People_sq_km+
                                Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                                Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                                Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                                Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                                Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                                Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018,data = Data_imp)
summary(modelFit_both_death_hosp)

#Type II Manova - variables still different from null given both dependent outcomes
Anova(modelFit_both_death_hosp)

#scale outcomes for walds test of coefficient equality
Data_imp =  as.data.frame(scale(Data_imp))
summary(Data_imp)

#significant Coefficients from Manova, equality test across seemingly unrelated regressions
#https://andrewpwheeler.com/2017/06/12/testing-the-equality-of-coefficients-same-independent-different-dependent-variables/?unapproved=9632&moderation-hash=0c84f568615a42d280c784fc5a66379c#comment-9632
covidHospmod = Death_Registrations_Hospital_Covid_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018

allHospmod = Death_Registrations_Hospital_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018

fitsur3 <- systemfit(list(covidHospmod = covidHospmod, allHospmod = allHospmod), data=Data_imp, method="SUR")
summary(fitsur3)

#testing whether effect is equivalent for both models
Covid_more_3_scot <- glht(fitsur3,linfct = c("covidHospmod_Scotland - allHospmod_Scotland= 0"))
summary(Covid_more_3_scot) #p = 0.0117
Covid_more_3_wales <- glht(fitsur3,linfct = c("covidHospmod_Wales - allHospmod_Wales= 0"))
summary(Covid_more_3_wales) #p = 0.00000041
Covid_more_3_air <- glht(fitsur3,linfct = c("covidHospmod_PM_2_5_Total_2018 - allHospmod_PM_2_5_Total_2018 = 0"))
summary(Covid_more_3_air) #p = 0.015
Covid_more_3_temp <- glht(fitsur3,linfct = c("covidHospmod_Average_tempHighCel_23_2_1_3_8_3_15_3 - allHospmod_Average_tempHighCel_23_2_1_3_8_3_15_3 = 0"))
summary(Covid_more_3_temp) #p = 0.289
Covid_more_3_age <- glht(fitsur3,linfct = c("covidHospmod_Median_Age_2018 - allHospmod_Median_Age_2018 = 0"))
summary(Covid_more_3_age) #p = 0.359
Covid_more_3_bame <- glht(fitsur3,linfct = c("covidHospmod_BAME_Including_Jews - allHospmod_BAME_Including_Jews = 0"))
summary(Covid_more_3_bame) #p = 0.00615
Covid_more_3_comm <- glht(fitsur3,linfct = c("covidHospmod_Lives_In_Communal - allHospmod_Lives_In_Communal= 0"))
summary(Covid_more_3_comm) #p = 0.0211
Covid_more_3_crowd <- glht(fitsur3,linfct = c("covidHospmod_Over_1_5_Per_Room - allHospmod_Over_1_5_Per_Room= 0"))
summary(Covid_more_3_crowd) #p = 0.948
Covid_more_3_dep <- glht(fitsur3,linfct = c("covidHospmod_Local_Proportion_Decile_1_Z_Ind_IMD - allHospmod_Local_Proportion_Decile_1_Z_Ind_IMD = 0"))
summary(Covid_more_3_dep) #p =  0.0579
Covid_more_3_fem <- glht(fitsur3,linfct = c("covidHospmod_Females_2018 - allHospmod_Females_2018 = 0"))
summary(Covid_more_3_fem) #p = 0.00606
Covid_more_3_self <- glht(fitsur3,linfct = c("covidHospmod_Self_Funders_2017 - allHospmod_Self_Funders_2017 = 0"))
summary(Covid_more_3_self) #p = 0.0191
Covid_more_3_ckd <- glht(fitsur3,linfct = c("covidHospmod_CKD_Z_Ind - allHospmod_CKD_Z_Ind = 0"))
summary(Covid_more_3_ckd) #p = 0.00187
Covid_more_3_copd <- glht(fitsur3,linfct = c("covidHospmod_COPD_Z_Ind - allHospmod_COPD_Z_Ind = 0"))
summary(Covid_more_3_copd) #p = 0.00418
Covid_more_3_cvd <- glht(fitsur3,linfct = c("covidHospmod_Age_Standardised_CVD_Deaths_U75_100000 - allHospmod_Age_Standardised_CVD_Deaths_U75_100000 = 0"))
summary(Covid_more_3_cvd) #p = 0.000123
Covid_more_3_dem <- glht(fitsur3,linfct = c("covidHospmod_Dementia_Z_Ind - allHospmod_Dementia_Z_Ind = 0"))
summary(Covid_more_3_dem) #p =  0.0162
Covid_more_3_smok <- glht(fitsur3,linfct = c("covidHospmod_Current_Smoking_Z_Ind - allHospmod_Current_Smoking_Z_Ind = 0"))
summary(Covid_more_3_smok) #p = 0.0195
Covid_more_3_exs <- glht(fitsur3,linfct = c("covidHospmod_Ex_Smoker_Z_Ind - allHospmod_Ex_Smoker_Z_Ind = 0"))
summary(Covid_more_3_exs) #p = 0.241

View(p.adjust(p = c(0.0117,0.00000041,0.015,0.289,0.359,0.00615,0.0211,0.948,0.0579,0.00606,0.0191,0.00187,0.00418,0.000123,0.0162,0.0195,0.241), method = "fdr"))

############################################################
#Multivariate multiple regression at home setting covid, vs at home setting except covid

both_death_home_data = covid_death_home[,-c(1,5:8)]

#standardise everything except outcome
both_death_home_data[,-29] = scale(both_death_home_data[,-29])

both_death_home_data$Death_Registrations_Home_Week_16_Per_100000 = minus_home_all$Death_Registrations_Home_Week_16_Per_100000
summary(both_death_home_data)

tempData_home <- mice(both_death_home_data,seed=987)

#multiple mulvariate regression doesn't work with multiple datasets
Data_imp = mice::complete(tempData_home,1)

modelFit_both_death_home = lm(cbind(Death_Registrations_Home_Covid_Week_16_Per_100000, Death_Registrations_Home_Week_16_Per_100000)~ Scotland+Wales+Density_People_sq_km+
                                Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
                                Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
                                Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
                                Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
                                Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
                                Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018,data = Data_imp)
summary(modelFit_both_death_home)

#Type II Manova - variables still different from null given both dependent outcomes
Anova(modelFit_both_death_home)

#scale outcomes for walds test of coefficient equality
Data_imp =  as.data.frame(scale(Data_imp))
summary(Data_imp)

#significant Coefficients from Manova, equality test across seemingly unrelated regressions
#https://andrewpwheeler.com/2017/06/12/testing-the-equality-of-coefficients-same-independent-different-dependent-variables/?unapproved=9632&moderation-hash=0c84f568615a42d280c784fc5a66379c#comment-9632
covidHomemod = Death_Registrations_Home_Covid_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018

allHomemod = Death_Registrations_Home_Week_16_Per_100000 ~ Scotland+Wales+Density_People_sq_km+
  Local_Proportion_Decile_1_Z_Ind_IMD+Females_2018+Lives_In_Communal+BAME_Including_Jews+
  Obese_Z_Ind+Over_1_5_Per_Room+Median_Age_2018+Current_Smoking_Z_Ind+Ex_Smoker_Z_Ind+
  Tonnes_Port_Over_750+Passengers_To_From_Pop_Stand_Normalised+Age_Standardised_CVD_Deaths_U75_100000+
  Diabetes_Prevalence_Z_Ind+RA_Z_Ind+HTN_Z_Ind+COPD_Z_Ind+Cancer_Z_Ind+CKD_Z_Ind+
  Average_Humid_23_2_1_3_8_3_15_3+Average_tempHighCel_23_2_1_3_8_3_15_3+Coronavirus_Google_010220_190320+
  Dementia_Z_Ind+Self_Funders_2017+PM_2_5_Total_2018

fitsur4 <- systemfit(list(covidHomemod = covidHomemod, allHomemod = allHomemod), data=Data_imp, method="SUR")
summary(fitsur4)

#testing whether effect is equivalent for both models
Covid_more_4_scot <- glht(fitsur4,linfct = c("covidHomemod_Scotland - allHomemod_Scotland= 0"))
summary(Covid_more_4_scot) #p = 0.342
Covid_more_4_wales <- glht(fitsur4,linfct = c("covidHomemod_Wales - allHomemod_Wales= 0"))
summary(Covid_more_4_wales) #p = 0.668
Covid_more_4_hum <- glht(fitsur4,linfct = c("covidHomemod_Average_Humid_23_2_1_3_8_3_15_3 - allHomemod_Average_Humid_23_2_1_3_8_3_15_3 = 0"))
summary(Covid_more_4_hum) #p = 0.0196
Covid_more_4_temp <- glht(fitsur4,linfct = c("covidHomemod_Average_tempHighCel_23_2_1_3_8_3_15_3 - allHomemod_Average_tempHighCel_23_2_1_3_8_3_15_3 = 0"))
summary(Covid_more_4_temp) #p = 0.00108
Covid_more_4_age <- glht(fitsur4,linfct = c("covidHomemod_Median_Age_2018 - allHomemod_Median_Age_2018 = 0"))
summary(Covid_more_4_age) #p = 0.294
Covid_more_4_dep <- glht(fitsur4,linfct = c("covidHomemod_Local_Proportion_Decile_1_Z_Ind_IMD - allHomemod_Local_Proportion_Decile_1_Z_Ind_IMD = 0"))
summary(Covid_more_4_dep) #p =  0.595
Covid_more_4_fem <- glht(fitsur4,linfct = c("covidHomemod_Females_2018 - allHomemod_Females_2018 = 0"))
summary(Covid_more_4_fem) #p = 0.0229
Covid_more_4_cancer <- glht(fitsur4,linfct = c("covidHomemod_Cancer_Z_Ind - allHomemod_Cancer_Z_Ind = 0"))
summary(Covid_more_4_cancer) #p = 0.0163
Covid_more_4_cvd <- glht(fitsur4,linfct = c("covidHomemod_Age_Standardised_CVD_Deaths_U75_100000 - allHomemod_Age_Standardised_CVD_Deaths_U75_100000 = 0"))
summary(Covid_more_4_cvd) #p = 0.00376
Covid_more_4_bp <- glht(fitsur4,linfct = c("covidHomemod_HTN_Z_Ind - allHomemod_HTN_Z_Ind = 0"))
summary(Covid_more_4_bp) #p =  0.0245
Covid_more_4_ra <- glht(fitsur4,linfct = c("covidHomemod_RA_Z_Ind - allHomemod_RA_Z_Ind = 0"))
summary(Covid_more_4_ra) #p =  0.186
Covid_more_4_flight <- glht(fitsur4,linfct = c("covidHomemod_Passengers_To_From_Pop_Stand_Normalised - allHomemod_Passengers_To_From_Pop_Stand_Normalised = 0"))
summary(Covid_more_4_flight) #p = 0.0128

View(p.adjust(p = c(0.342,0.668,0.0196,0.00108,0.294,0.595,0.0229,0.0163,0.00376,0.0245,0.186,0.0128), method = "fdr"))

#####################################################################
#END
