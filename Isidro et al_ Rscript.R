setwd("C:/Users/isidrom/Desktop/PhD/Halimeda_Enrichment/Stats/For paper submition")

#---------Libraries used----------------
library(ggplot2)
library(extrafont)
library(dplyr)
library(stringr)
library(lme4)
library(lmerTest)
library(lmtest)
library(car)
library(tidyr)
library(stats)
library(ggsignif)
library(ggpubr)
library(lsmeans)
library(devtools)
library(permuco)
#install_github("onofriandreapg/aomisc")
library(aomisc)
library(statforbiology)
library(broom)
library(purrr)
library(nls2)
library(tibble)


#--------------------------Extra fonts for graphs-----------------------------
extrafont::loadfonts()

#-------------------------------- Theme for plots-----------------------------------
newtheme <- theme_classic() + theme(text = element_text(size=14,family="Arial"))+
  theme(axis.text.x = element_text(size=12,colour="black"), axis.text.y = element_text(size=12,colour="black"))+
  theme(plot.margin = unit(c(5.5,5.5,5.5,20), "pt"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(panel.background = element_blank())+
  theme(axis.line = element_blank())+theme()

#------------------------------Set sampling days for plotting----------------------------
day_breaks=c(0,7,10,14,20,26,36,42)
day_breaks2=c(7,10,14,20,26,36,42)

#-----------------------------------------------------------------------------------------
#   Temporal changes in d15N values of Halimeda opuntia in response to nutrient enrichment
#-----------------------------------------------------------------------------------------

#Open dataframes
hal_data<-read.csv("hal_enrichment.csv")
colnames(hal_data) <- c("Treatment", "Time_point","Days", "Plot", "Replicate","ID", "d15N", "d13C", "perctN", "perctC", "ratioCN")

#Delete sketchy values from Control T42
hal_data <- hal_data[-(239:249),]


#-----------------------Calculate means, sd, se, and sample size----------------------------
#Per day
hal_means<-hal_data %>% group_by(Treatment,Time_point,Days) %>% 
  summarize(N=sum(!is.na(d15N)),
            mean_d15N=mean(d15N,na.rm=T),
            mean_CN=mean(ratioCN,na.rm=T),
            sd_d15N=sd(d15N),
            sd_CN= sd(ratioCN),
            se_d15N=sd(d15N)/sqrt(N),
            se_CN=sd(ratioCN)/sqrt(N)) %>% ungroup()

#Per day and per replicate plot
hal_plot_mean <- hal_data %>%  group_by(Treatment, Time_point, Days, Plot) %>%
  summarize(N = sum(!is.na(d15N)),
            sd_d15N = sd(d15N, na.rm = TRUE),
            sd_ratioCN = sd(ratioCN, na.rm = TRUE),
            mean_d15N = mean(d15N, na.rm = TRUE),
            mean_CN = mean(ratioCN, na.rm = TRUE),
            se_d15N = sd(d15N, na.rm = TRUE) / sqrt(sum(!is.na(d15N))),
            se_CN = sd(ratioCN, na.rm = TRUE) / sqrt(sum(!is.na(ratioCN)))) %>%  ungroup()


#------------Mixed effect model for Halimeda opuntia d15N values across time--------------

#Write up model for Halimeda opuntia d15N values across time. Plot and Replicate is a random effect
d15N_lmer<-lmer(d15N~Treatment*as.factor(Days)+(1|Plot)+(1|Replicate),data=hal_data)
summary(d15N_lmer) 
anova(d15N_lmer,type=2)

# Look at Plot as random effect
ranova(d15N_lmer)   ## Plot as random effect is significant--> δ15N varies between different plots, even after accounting for Treatment, Days, and their interaction
sum(resid(d15N_lmer)^2)

#Check residuals 
res <- resid(d15N_lmer)
hist(res, breaks = 20)
qqnorm(res)
qqline(res, col = "red")
shapiro.test(res)   #0.2006--> residuals are normal

#Levene's Test 
fixed_model <- lm(d15N ~ Treatment * as.factor(Days), data = subset(hal_data, !Days %in% c(1, 2, 3, 5)))
leveneTest(residuals(fixed_model) ~ Treatment * as.factor(Days), 
           data = subset(hal_data, !Days %in% c(1, 2, 3, 5)))  #0.945 assumption of homogeneity is met


##Plot line plot for Halimeda opuntia d15N values across time 
# All days in your data (for tick marks)
all_days <- seq(min(hal_means$Days), max(hal_means$Days), by = 1)

# Define X axis lable
main_labels <- c(0, 7, 10, 14, 20, 26, 36, 42)
asterisk_days <- c(14, 20, 26, 36, 42)

d15N_enrichment <- ggplot() +
  geom_line(data = hal_means, aes(x = Days, y = mean_d15N, group = Treatment, color = Treatment)) +
  geom_errorbar(data = hal_means,
                aes(x = Days, y = mean_d15N,
                    ymin = mean_d15N - (se_d15N * 1.96), 
                    ymax = mean_d15N + (se_d15N * 1.96),
                    color = Treatment),
                linewidth = 2, width = 0) +
  geom_point(data = hal_means,
             aes(x = Days, y = mean_d15N, fill = Treatment, shape = Treatment),
             col = 'black', size = 5) +
  geom_point(data = hal_data,
             aes(x = Days, y = d15N, color = Treatment, fill = Treatment, shape = Treatment),
             col = 'black',
             position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2),
             alpha = 0.5) +
  scale_shape_manual(values = c(21, 24)) +
  geom_hline(yintercept = 0.3, lty = 2) +
  geom_hline(yintercept = 0.3 + 0.153, linetype = "dotted", color = "grey85", linewidth = 0.5) +
  geom_hline(yintercept = 0.3 - 0.153, linetype = "dotted", color = "grey85", linewidth = 0.5) +
  scale_x_continuous(breaks = all_days,labels = ifelse(all_days %in% main_labels, all_days, "")) +
  annotate("text", x = asterisk_days,y = 3.9,  label = "*",  size = 8, vjust = 1) +
  scale_color_manual(values = c("#b9c9cc", "#5fbfaf")) +
  scale_fill_manual(values = c("#b9c9cc", "#5fbfaf")) +
  xlab("Days") +
  ylab(expression({delta}^15*N~'(\u2030)')) +
  theme(axis.ticks.length = unit(0.25, "cm"),axis.ticks = element_line(color = "black") ) +newtheme
d15N_enrichment


#----------------Mixed effect model for Halimeda opuntia C:N values across time-----------------

#write model for Halimeda opuntia C:N ratios across time. Plot is a ramdom effect
CN_lmer<-lmer(ratioCN~Treatment*as.factor(Days)+(1|Plot)+(1|Replicate),data=hal_data)
summary(CN_lmer)
anova(CN_lmer,type=2)

# Look at Plot as random effect
ranova(CN_lmer)       #C:N varies between plots but not replicates 
sum(resid(CN_lmer)^2)

#Check residuals 
res <- resid(CN_lmer)
hist(res, breaks = 20)
qqnorm(res)
qqline(res, col = "red")
shapiro.test(res)   #0.01--> residuals are not normal

#Levene's Test
fixed_model2 <- lm(ratioCN ~ Treatment * as.factor(Days), data = hal_data)
leveneTest(residuals(fixed_model2) ~ Treatment * as.factor(Days), 
           data = hal_data)  #p<0.05 assumption of homogeneity is violeted

#Breusch-Pagan test
bptest(fixed_model2) #p> 0.03 homogeneous

##Plot line plot for Halimeda opuntia C:N values across time 
# All days for tick marks
all_days <- seq(min(day_breaks), max(day_breaks), by = 1)

main_labels <- day_breaks
asterisk_days <- c(7, 10, 14, 20, 26, 36, 42)  

ratioCN_enrichment <- ggplot() +
  geom_line(data = hal_means, aes(x = Days, y = mean_CN, group = Treatment, color = Treatment)) +
  geom_errorbar(data = hal_means,
                aes(x = Days, y = mean_CN,
                    ymin = mean_CN - (se_CN * 1.96),
                    ymax = mean_CN + (se_CN * 1.96),
                    color = Treatment),
                linewidth = 2, width = 0) +
  geom_point(data = hal_means,
             aes(x = Days, y = mean_CN, fill = Treatment, shape = Treatment),
             col = 'black', size = 5) +
  geom_point(data = hal_data,
             aes(x = Days, y = ratioCN, color = Treatment, fill = Treatment, shape = Treatment),
             col = 'black',
             position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2),
             alpha = 0.5) +
  scale_shape_manual(values = c(21, 24)) +
  scale_x_continuous(breaks = all_days,
                     labels = ifelse(all_days %in% main_labels, all_days, "")) +
    annotate("text", x = asterisk_days, y = 17.8, label = "*", size = 8, vjust = 1) +
  expand_limits(y = asterisk_days) +
  scale_color_manual(values = c("#b9c9cc", "#5fbfaf")) +
  scale_fill_manual(values = c("#b9c9cc", "#5fbfaf")) +
  labs(x = "Days", y = "C:N") +
  theme(axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(color = "black")) + newtheme
ratioCN_enrichment

#------------------------------------------------------------------------------------
#          Regression curve for calculating turnover times & Bootstrap
#---------------------------------------------------------------------------------------

#Temporal changes for each ID for the ones with 4 or more time points
hal_data_filtered <- hal_data %>%
  group_by(ID) %>%
  filter(n_distinct(Days) >= 4) %>%   
  ungroup()

hal_data_filtered <- hal_data %>%
  group_by(ID, Treatment) %>%
  filter(n_distinct(Days) >= 4) %>%   
  ungroup()

# All days for tick marks
all_days <- seq(min(day_breaks), max(day_breaks), by = 1)
main_labels <- day_breaks

reps_time <- ggplot(hal_data_filtered, aes(x = Days, y = d15N, group = ID, color = ID)) +
  geom_line(show.legend = FALSE) +
  geom_point(size = 1.5, show.legend = FALSE) +
  facet_wrap(~Treatment) +
  scale_x_continuous( breaks = all_days,labels = ifelse(all_days %in% main_labels, all_days,"")) +
  labs(x = "Days", y = expression(delta^{15}*N)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +  newtheme
reps_time


#Calculate difference between d15N in enriched Halimeda and the mean of fertilizer(=0.3‰)
df1 <- hal_data %>%
  dplyr::select(Treatment, Days, Plot, ID, d15N) %>%       
  tidyr::pivot_wider(names_from = Treatment, values_from = d15N) %>% 
  dplyr::mutate(d15N_diff = Enriched - 0.3) %>%      
  tidyr::drop_na(Enriched) %>%
  dplyr::select(Days, Plot, ID, d15N_diff)


model2<- nls(d15N_diff ~ NLS.expoDecay(Days, a,c),
             data =df1)

# Extract estimated model parameters
params <- coef(model2)
a_est <- as.numeric(params["a"])
c_est <- as.numeric(params["c"])

c_se  <- coef(summary(model2))["c", "Std. Error"]

# Half-life (t0.5)
half_life <- log(2) /c_est    #11.23 days 
half_life

#Half-life with SE
half_life_se <- (log(2) / (c_est^2)) * c_se
half_life_se           #1.18


#Full Turnover (t95)
time_to_turnover <- log(1 - 0.95) / (-c_est)   #48.52 days 
time_to_turnover

#Full turnover with SE
time_turover_se <- (log(1 - 0.95) / (c_est^2)) * c_se
time_turover_se         #5.11


##Plot regression model 

# Generate Days for plotting
curve_days <- data.frame(Days = seq(min(df1$Days), max(df1$Days), length.out = 200))

# Predicted values
curve_days$predicted <- predict(model2, newdata = curve_days)

# Partial derivatives for SE
X <- cbind(
  a = 1 - exp(-coef(model2)["c"] * curve_days$Days),  
  c = coef(model2)["a"] * curve_days$Days * exp(-coef(model2)["c"] * curve_days$Days))

# Standard errors and 95% CI
V <- vcov(model2)
se <- sqrt(rowSums((X %*% V) * X))
curve_days$lower <- curve_days$predicted - 1.96 * se
curve_days$upper <- curve_days$predicted + 1.96 * se

# Add predictions to raw data 
df1$predicted <- predict(model2, newdata = df1)

#----------------------------Plot Turnover rate with 95% CI----------------------
# All days for tick marks
all_days <- seq(min(day_breaks), max(day_breaks), by = 1)
main_labels <- day_breaks

model2_turnover <- ggplot() +
  geom_point(data = df1, aes(x = Days, y = d15N_diff),
             size = 3, fill = "deepskyblue4", pch = 21, col = "black") +
  geom_line(data = curve_days, aes(x = Days, y = predicted),
            color = "black", linewidth = 1.2) +
  geom_ribbon(data = curve_days, aes(x = Days, ymin = lower, ymax = upper),
              linetype = "dotted", fill = "grey", alpha = 0.2) +
  scale_x_continuous(breaks = all_days,
                     labels = ifelse(all_days %in% main_labels, all_days, "")) +
  labs(x = "Days") +ylab(expression(Delta^{15}*N ~ (delta^{15}*N[Halimeda] - delta^{15}*N[Fertilizer]))) +
  newtheme
model2_turnover

#---------------------------------Bootstrap-----------------------------------
# Get number of observations per Days group
nested <- df1 %>%
  group_by(Days) %>%
  nest() %>%
  ungroup() %>%
  mutate(n = map_int(data, nrow))

#Bootstrap settings
nboot <- 100 #final with 10000
set.seed(123)

all <- NULL 

# To store results
all_preds <- list()
all_params <- list()

# Prediction grid for smooth curves
pred_days <- seq(min(df1$Days), max(df1$Days), length.out = 100)

# Bootstrap loop
for (i in 1:nboot) {
  
  # Sample with replacement within each Days group
  samps <- nested %>%
    mutate(samp = map2(data, n, ~ sample_n(.x, size = .y, replace = TRUE)))
  
  # Combine sampled data
  samps2 <- samps %>%
    dplyr::select(-data) %>%
    tidyr::unnest(samp) %>%
    dplyr::select(n, Days, d15N_diff)
  
  
  # Fit model on sampled data
  model_boot <- nls(d15N_diff ~ NLS.expoDecay(Days, a, c),
                    data = samps2)
  
  
  # Extract coefficients
  a_est <- coef(model_boot)[["a"]]
  c_est <- coef(model_boot)[["c"]]
  
  # Derived metrics
  half_life <- log(2) /c_est   
  time_to_turnover <- log(1 - 0.95) / (-c_est)   
  
  # Store predictions
  preds <- tibble(
    Days = pred_days,
    pred = predict(model_boot, newdata = tibble(Days = pred_days)),
    Iteration = i)
  
  # Store parameters
  params_tb <- tibble(
    Iteration = i,
    a = a_est,
    c = c_est,
    half_life = half_life,
    time_to_turnover = time_to_turnover)
  
  # Append to list
  all_preds[[i]] <- preds
  all_params[[i]] <- params_tb
 
}


#Check if n matches with df_diff_reps
check<-samps2%>% group_by(Days) %>% summarise(N=sum(!is.na(d15N_diff)))


# Combine results
preds <- bind_rows(all_preds)
params <- bind_rows(all_params)


# Histograms to see distribution of data
ggplot(params, aes(x = half_life)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(x = "Half-life", y = "Count") 

ggplot(preds, aes(x = Days, y = pred, group = Iteration)) +
  geom_line(alpha = 0.1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(x = "Days", y = "Predicted δ15N difference")


# Summaries bootstrap predictions by Days
ci_boot <- preds %>%
  group_by(Days) %>%
  summarise(
    median_pred = median(pred, na.rm = TRUE),
    mean_pred =  mean(pred, na.rm = TRUE),
    lower_95  = quantile(pred, 0.025, na.rm = TRUE),
    upper_95  = quantile(pred, 0.975, na.rm = TRUE),
    lower_75  = quantile(pred, 0.125, na.rm = TRUE),
    upper_75  = quantile(pred, 0.875, na.rm = TRUE),
    .groups = "drop"  )

#Estimate mean half-life and CI
half_life_summary <- params %>%
  summarise(
    mean_hl = mean(half_life, na.rm = TRUE), 
    mediab_hl= median(half_life, na.rm = TRUE),
    lower_95 = quantile(half_life, 0.025, na.rm = TRUE),
    upper_95 = quantile(half_life, 0.975, na.rm = TRUE))
half_life_summary

#mean_hl mediab_hl lower_95 upper_95
# 11.3      11.2     9.27     13.7

##Estimate mean time to turover and CI
time_to_turnover_summary <- params %>%
  summarise(
    mean_ft = mean(time_to_turnover, na.rm = TRUE), 
    mediannft= median(time_to_turnover, na.rm = TRUE),
    lower_95 = quantile(time_to_turnover, 0.025, na.rm = TRUE),
    upper_95 = quantile(time_to_turnover, 0.975, na.rm = TRUE))
time_to_turnover_summary

#mean_ft median_ft lower_95 upper_95
# 48.9      48.6     40.0     59.3

#------------------------Plot Bootstrap with 95% CI-------------------------
boot.d15N<- ggplot(ci_boot, aes(x = Days, y = median_pred)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95),    # 95% CI ribbon
              fill = "grey", alpha = 0.25) +
  geom_line(color = "black", linewidth = 1) +
  geom_jitter(data = df1,
              aes(x = Days, y = d15N_diff),
              width = 0.2, height = 0,
              fill = "deepskyblue4", size = 2, pch=21, col='black') +
  scale_x_continuous(breaks = day_breaks) +
  labs(x = "Days",y = expression(Delta^{15}*N ~ (delta^{15}*N[Halimeda] - 
        delta^{15}*N[Fertilizer]))) +  newtheme
boot.d15N

#All individual bootstrap------- no 95% intervals 
# All days for tick marks
all_days <- seq(min(day_breaks), max(day_breaks), by = 1)
main_labels <- day_breaks

bootall.d15N <- ggplot() +
  geom_line(data = preds,            # individual bootstraps
            aes(x = Days, y = pred, group = Iteration),
            color = "grey", alpha = 0.25, linewidth = 0.3) +
  geom_line(data = ci_boot,         #solid line- median 
            aes(x = Days, y = median_pred),
            color = "black", linewidth = 1.2) +
  geom_jitter(data = df1,
              aes(x = Days, y = d15N_diff),
              width = 0.2, height = 0,
              fill = "deepskyblue4", size = 2, pch=21, col='black') +
  scale_x_continuous( breaks = all_days,labels = ifelse(all_days %in% main_labels, all_days,"")) +
  labs(x = "Days",y = expression(Delta^{15}*N ~ 
  (delta^{15}*N[Halimeda] - delta^{15}*N[Fertilizer]))) + newtheme
bootall.d15N


#-------------------------------------------------------------------------------------
#         Spatial variation in d15N values across Halimeda opuntia thallus
#-------------------------------------------------------------------------------------

#Open dataset
hal_thallus<-read.csv("hal_thallus.csv")
colnames(hal_thallus) <- c("Location", "Time_point","Plant", "Replicate", "d15N", "d13C", "perctN", "perctC", "ratioCN")

#Location and time point as factors
hal_thallus$Time_point <- as.factor(hal_thallus$Time_point)
hal_thallus$Location <- factor(hal_thallus$Location , levels=c("Top", "Mid", "Bottom", "Holdfast"))

#Delete sketchy values from plant E
hal_thallus <- hal_thallus[-c(43,44,121,122,123),]

#Calculate mean, sd, se by location and time point 
thallus_means_location <-hal_thallus %>% group_by(Location,Time_point) %>% 
  summarize(N=sum(!is.na(d15N)),
            mean_d15N=mean(d15N,na.rm=T),
            mean_CN=mean(ratioCN,na.rm=T),
            se_d15N=sd(d15N)/sqrt(N),
            se_CN=sd(ratioCN)/sqrt(N),
            sd_d15N=sd(d15N),
            sd_CN=sd(ratioCN)) %>% ungroup()

#Remove holdfast due to to low replication and lack of considerable difference in means after 15 days
hal_thallus<- hal_thallus %>%
  filter(Location != "Holdfast")

thallus_delta <- thallus_means_location %>%
  dplyr::select(Location, Time_point, mean_d15N) %>%
  tidyr::pivot_wider(names_from = Time_point, values_from = mean_d15N) %>%
  dplyr::mutate(delta_d15N = `15` - `0`)

thallus_delta <- thallus_delta %>%
  mutate(delta_d15N = `15` - `0`)


#----------Linear Mixed model-differences between different locations in thallus-----------------------------------------

#Write up model for d15N values in different locations across time point. Replicate is a random effect
d15N_thallus <- lmer(d15N ~ Location * Time_point + (1 | Plant), data = hal_thallus)
summary(d15N_thallus)
anova(d15N_thallus,type=2)

# Look at Plot as random effect
ranova(d15N_thallus)       #Plot varies
sum(resid(d15N_thallus)^2)

#Check residuals 
res <- resid(d15N_thallus)
hist(res, breaks = 20)
qqnorm(res)
qqline(res, col = "red")
shapiro.test(res)   #0.06--> residuals are not normal

#Levene's Test
fixed_model2 <- lm(d15N ~ Location * Time_point, data = hal_thallus)
leveneTest(d15N ~ Location * Time_point, data = hal_thallus)   #0.29-> homogeneous


#----------Boxplot of Halimeda opuntia d15N across diferent locations across time ------------------
thallus_d15N<- ggplot(hal_thallus, aes(x=Time_point, y=d15N, fill=Location))+
  geom_boxplot(outlier.colour = NA)+ geom_point(aes(fill=Location,group=Location),shape=21,position=position_jitterdodge(.2), alpha=0.5)+
  scale_color_manual(values = c("#b9c9cc", "#68A4A5", "#4C8055","#31473A"))+
  scale_fill_manual(values = c("#b9c9cc", "#68A4A5", "#4C8055","#31473A"))+
  labs (x= "Day", y=expression({delta}^15*N~'(\u2030)'))+newtheme
thallus_d15N


#--------------------------------------------------------------------------------
#           Effect of decalcification in Halimeda opuntia stable isotopic values
#-----------------------------------------------------------------------------

#Open dataset
hal_t0<-read.csv("halimeda_t0_caltest.csv")
colnames(hal_t0) <- c("Time_point", "Treatment","Plots", "Replicate", "d15N", "d13C", "perctN", "perctC", "ratioCN")


#-----------------------------Calculate means, sd, se, ----------------------------

#Create new column Calc which indicates if samples are calcified or decalcified 
#and another column Enriched which indicates treatment level(enriched or control)  
hal_t0<-hal_t0 %>% mutate(Calc=ifelse(str_detect(Treatment, "_"), "Calcified", "Non-Calcified"),
                          Enriched=ifelse(str_detect(Treatment,"Enriched"),"Enriched","Control"))

#Calculate mean, SD, SE d15N values grouped by calcification status (Calc), replication plot, and enrichment treatment
hal_t0_means<-hal_t0 %>% group_by(Calc,Plots,Enriched ) %>% 
  summarize(N=sum(!is.na(d15N)),
            mean_d15N=mean(d15N,na.rm=T),
            mean_CN=mean(ratioCN,na.rm=T),
            sd_d15N=sd(d15N),
            sd_CN= sd(ratioCN),
            se_d15N=sd(d15N)/sqrt(N),
            se_CN=sd(ratioCN)/sqrt(N)) %>% ungroup()


#---------------------------------- ANOVA Test --------------------------------------
#Check normality of d15N values
hist(hal_t0$d15N)
ggqqplot(hal_t0$d15N)
shapiro.test(hal_t0$d15N)  #0.34, normal distribution
leveneTest(d15N ~ Calc, data = hal_t0)   #0.03 -- unequal variances 

#2 way ANOVA
d15N_calc <- aov(d15N~Calc*Plots, data = hal_t0)
summary(d15N_calc)

#Residuals 
res <- residuals(d15N_calc)
qqnorm(res)
qqline(res, col = "red", lwd = 2)
hist(res, col = "lightblue", border = "black")
shapiro.test(res)        #0.10 - normal distributed 
leveneTest(res ~ hal_t0$Calc)  # 0.17--- equal variances


#------------Regression line with d15N values for calcified and decalcified-------------

#Create column with d15N values per calcification level
hal_wide<-hal_t0 %>% dplyr::select(Plots,Replicate,Treatment,d15N) %>%
  pivot_wider(names_from=Treatment,values_from=d15N)  

#Graph with regression line 
lr.rep <- ggplot() +
  geom_point(data = hal_wide, mapping = aes(x = Enriched_Calcified, y = Enriched, color = "Enriched"), size = 3) +
  geom_abline() +
  geom_point(data = hal_wide, mapping = aes(x = Control_Calcified, y = Control, color = "Control"), size = 3) +
  labs(x = "Calcified", y = "De-calcified") +
  scale_color_manual(values = c("Enriched" = "#b9c9cc", "Control" = "#5fbfaf")) +
  guides(fill = guide_legend(title = NULL))+  newtheme
lr.rep


# ----------Boxplot with δ15N across plots, grouped by calcification status and faceted by Enrichment treatment------------------

#d15N values
d15N_calc <- ggplot(hal_t0, aes(x = Plots, y = d15N, fill = Calc)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_point(aes(group = Calc), shape = 21, 
             position = position_jitterdodge(0.2), alpha = 0.5) +
  scale_fill_manual(values = c("darkseagreen", "darkseagreen2")) +
  newtheme +
  xlab("Plots") +
  ylab(expression({delta}^15*N~'(\u2030)')) +
  guides(fill = guide_legend(title = NULL))
d15N_calc


#------------------------------------------------------------------------------------
#                        Nutrient Concentration in water samples
#--------------------------------------------------------------------------------------

#Open dataframes
nutrient<-read.csv("nutrient_water.csv")
colnames(nutrient) <- c("Treatment", "Plot", "Time_point", "Days","Replicate", "Silica", "Nitrite", "Nitrate", "Phosphate", "Silica_mol", "Nitrite_mol", "Nitrate_mol", "Phosphate_mol")

#Remove Silica data 
nutrient <- nutrient[, !(names(nutrient) %in% c("Silica", "Silica_mol"))]

#Remove negative values in Phosphate 
nutrient <- nutrient %>%
  filter(!is.na(Phosphate),
         !is.na(Phosphate_mol), Phosphate >= 0,Phosphate_mol >= 0)

#Calculate DIN (Nitrite+ Nitrate) in a new column
nutrient$DIN_mol <- nutrient$Nitrate_mol + nutrient$Nitrite_mol

#Average replicates
nutrient_average<-nutrient %>% 
  group_by(Treatment, Plot, Time_point,Days) %>% 
  summarize(DIN=mean(DIN_mol,na.rm=T),
            Phosphate=mean(Phosphate_mol,na.rm=T))

#Reorganize data for plot, different nutrients as rows  
nutrient_rows <- nutrient_average %>%
  pivot_longer(cols = c(DIN, Phosphate),
    names_to = "Nutrient",values_to = "umol_l" )

#Calculate mean, sd, se for the different nutrients at each sampling day
nutrient_summary <- nutrient_rows %>%
  filter(!is.na(umol_l), 
         !is.na(Treatment), 
         !is.na(Nutrient), 
         !is.na(Days)) %>%
  group_by(Treatment, Nutrient, Days) %>%   
  summarise(mean = mean(umol_l, na.rm = TRUE),
            sd = sd(umol_l, na.rm = TRUE),
           n = n(),
           se = sd / sqrt(n), .groups = "drop") %>%
  arrange(Nutrient, Treatment, Days)


#----------------------Linear mixed model for nutrient values -----------------------------------
##DIN
leveneTest(nutrient_average$ DIN~nutrient_average$Treatment) 
qqp(nutrient_average$DIN)

DIN.mod<-lmer(DIN~Treatment*as.factor(Days)+(1|Plot),data=nutrient_average)
summary(DIN.mod)
anova(DIN.mod,type=2)
#unlike Phosphate - fit of this model not singular because are some quantifable differences among plots 

emmN = emmeans(DIN.mod, ~ Treatment|Days)
pairs(emmN,adjust = "bonferroni") #differences on days 20, 26, 36

#Calculate % difference between control and enriched plots during significantly different time points
means<-nutrient_average %>% group_by(Treatment, Days) %>% summarize(n=mean(DIN))
diffs<-means %>% pivot_wider(names_from=Treatment,values_from=n) %>% 
  mutate(Pdiff=((Enriched-Control)/Control)*100,diff=Enriched-Control)

##Phosphate
leveneTest(nutrient_average$Phosphate~nutrient_average$Treatment) 
qqp(nutrient_average$Phosphate)

P.mod<-lmer(Phosphate~Treatment*as.factor(Days)+(1|Plot),data=nutrient_average,REML=F)
summary(P.mod)
anova(P.mod,type=2)

emmN = emmeans(P.mod, ~ Days)
pairs(emmN,adjust = "bonferroni") # only significant timepoint is 36


# ------------------Plot Nutriment values across time -------------------------
# All days for tick marks
all_days <- seq(min(day_breaks), max(day_breaks), by = 1)
main_labels <- day_breaks 

nutrient_water <- ggplot(data = nutrient_summary, mapping = aes(Days, mean, col = Treatment)) +  
  geom_line(data = nutrient_summary, mapping = aes(x = Days, y = mean, group = Treatment, color = Treatment)) +
  geom_point(data = nutrient_summary, mapping = aes(x = Days, y = mean, fill = Treatment, shape = Treatment),
             col = 'black', size = 3.5) +
  geom_errorbar(data = nutrient_summary, mapping = aes(x = Days, y = mean,
               ymin = mean - (se * 1.96), ymax = mean + (se * 1.96), color = Treatment),
                linewidth = 1, width = 0, alpha = 0.5) +
  geom_blank( data = data.frame(Nutrient = "DIN", Days = 0, mean = 1.95,
      Treatment = NA) ) +
  facet_grid(rows = vars(Nutrient), scales = "free") +  
  labs(x = "Days", y = "μmol/l") +
  scale_color_manual(values = c("#b9c9cc", "#5fbfaf")) +
  scale_fill_manual(values = c("#b9c9cc", "#5fbfaf")) +
  scale_shape_manual(values = c(21, 24)) +
  scale_x_continuous(breaks = all_days,
  labels = ifelse(all_days %in% main_labels, all_days, "")) +
  scale_y_continuous(labels = function(x) format(round(x, 2), nsmall = 2)) +
  geom_text(data = data.frame(
      Nutrient = c("DIN", "DIN", "DIN", "Phosphate"),
      Days = c(20, 26, 36, 36), y = c(1.85, 1.85, 1.85, 0.10),  
      label = "*"), aes(x = Days, y = y, label = label),inherit.aes = FALSE,size = 6) +newtheme
nutrient_water
