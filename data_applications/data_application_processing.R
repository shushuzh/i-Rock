library(foreach)
library(conquer)
library(fGarch)
library(RColorBrewer)
library(dplyr)
library(compare)
library(doParallel)
library(MonoInc)
library(Rcpp)
library(RcppArmadillo)
library(fastkmedoids)
library(FNN)
library(SparseM)
library(quantreg)
library(splines)
library(zoo)
library(extraDistr)

us <- readLines("/home/shushuz/M_rock/birthweight_2022_us.txt")
ps <- readLines("/home/shushuz/M_rock/birthweight_2022_ps.txt")
all_lines = c(us,ps)
### Extract variables ###
data = data.frame(
  birthweight = as.numeric(unlist(lapply(all_lines,substr,start=504,stop=507))),
  race_mom = unlist(lapply(all_lines,substr,start=117,stop=117)),
  prenatal_visit = as.numeric(unlist(lapply(all_lines,substr,start=238,stop=239))),
  sex = unlist(lapply(all_lines,substr,start=475,stop=475)),
  age = as.numeric(unlist(lapply(all_lines,substr,start=75,stop=76))),
  diabete = unlist(lapply(all_lines,substr,start=314,stop=314)),
  hypertension = unlist(lapply(all_lines,substr,start=316,stop=316)),
  cigarette = as.numeric(unlist(lapply(all_lines,substr,start=264,stop=264))),
  plurality = as.numeric(unlist(lapply(all_lines,substr,start=454,stop=454))),
  WIC = unlist(lapply(all_lines,substr,start=251,stop=251)),
  education = as.numeric(unlist(lapply(all_lines,substr,start=124,stop=124))),
  marital = as.numeric(unlist(lapply(all_lines,substr,start=120,stop=120)))
  )
### only focus on male and single babies ###
data = data[data$sex=="M",]
data = data[data$plurality==1,]
### remove missing data ###
data = data[data$birthweight!=9999,]
data = data[data$race_mom!="8",]
data = data[data$WIC!="U",]
data = data[(data$diabete!="U"),]
data = data[(data$hypertension!="U"),]
data = data[(data$cigarette!=6),]
data = data[(data$education!=9),]
data = data[(data$marital!=9)&(!is.na(data$marital)),]
data = data[(data$prenatal_visit!=99),]

data <- data %>%
        mutate(marital_cat = case_when(
                                    marital == 1 ~ "Married",
                                    marital >1 ~ "Unmarried"
                                    )) %>%
        mutate(age_cat = case_when(
    age <= 19  ~ "<=19",
    age >= 20 & age <=34 ~ "20-34",
    age >= 35 ~ ">=35"
  )) %>%
        mutate(education_cat = case_when(
    education >= 1 & education <= 3  ~ "<=HS",
    education == 4 | education == 5 ~ "some college",
    education >= 6 ~ "college"
  )) %>%
        mutate(cigarette = ifelse(cigarette > 0, 1, 0))%>%
       mutate(prenatal_visit_category = case_when(
    prenatal_visit >= 0 & prenatal_visit <= 5 ~ "0-5",
    prenatal_visit >= 6 & prenatal_visit <= 10 ~ "6-10",
    prenatal_visit > 10 ~ ">10"
  ))
data$marital_cat <- relevel(factor(data$marital_cat),ref="Married")
data$age_cat <- relevel(factor(data$age_cat),ref="20-34")
data$education_cat <- relevel(factor(data$education_cat),ref = "<=HS")
data$cigarette <- relevel(factor(data$cigarette),ref = "0")
data = subset(data,race_mom %in% c(1,2,4,7))
data$race_mom <- factor(data$race_mom)
data$race_mom <- relevel(data$race_mom, ref = 1)
data$prenatal_visit_category <- relevel(factor(data$prenatal_visit_category),ref="0-5")
# write.csv(data,file="./application_data.csv")
design_matrix = model.matrix(~race_mom+prenatal_visit_category + diabete + hypertension + cigarette + age_cat +  WIC + marital_cat + education_cat,data=data)
data_design <- cbind(data$birthweight,design_matrix[,-1])
write.csv(data_design,file="./design_matrix_main.csv")
