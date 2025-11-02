.libPaths("/home/shushuz/R_libs")
library(foreach, lib.loc = "/home/shushuz/R_libs")
library(conquer, lib.loc = "/home/shushuz/R_libs")
library(fGarch, lib.loc = "/home/shushuz/R_libs")
library(RColorBrewer, lib.loc = "/home/shushuz/R_libs")
library(dplyr, lib.loc = "/home/shushuz/R_libs")
library(compare, lib.loc = "/home/shushuz/R_libs")
library(doParallel, lib.loc = "/home/shushuz/R_libs")
library(MonoInc, lib.loc = "/home/shushuz/R_libs")
library(Rcpp, lib.loc = "/home/shushuz/R_libs")
library(RcppArmadillo, lib.loc = "/home/shushuz/R_libs")
library(fastkmedoids, lib.loc = "/home/shushuz/R_libs")
library(FNN, lib.loc = "/home/shushuz/R_libs")
library(SparseM, lib.loc = "/home/shushuz/R_libs")
library(quantreg, lib.loc = "/home/shushuz/R_libs")
library(splines, lib.loc = "/home/shushuz/R_libs")
library(zoo, lib.loc = "/home/shushuz/R_libs")
#library(drf, lib.loc = "/home/shushuz/R_libs")
library(extraDistr, lib.loc = "/home/shushuz/R_libs")

#all_lines = readLines("/home/shushuz/M_rock/birthweight_2022_ps.txt")
us <- readLines("/home/shushuz/M_rock/birthweight_2022_us.txt")
ps <- readLines("/home/shushuz/M_rock/birthweight_2022_ps.txt")
all_lines = c(us,ps)
### Extract variables ###
data = data.frame(
  birthweight = as.numeric(unlist(lapply(all_lines,substr,start=504,stop=507))),
  race_mom = unlist(lapply(all_lines,substr,start=117,stop=117)),
  prenatal_visit = as.numeric(unlist(lapply(all_lines,substr,start=238,stop=239))),
  #weight_mom = as.numeric(unlist(lapply(all_lines,substr,start=292,stop=294))),
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
#data = data[data$weight_mom!=999,]
#data <- data %>%
#  mutate(diabete = case_when(
#    diabete == "N" ~ 0,
#    diabete == "Y" ~ 1
#  )) %>%
#  mutate(hypertension = case_when(
#    hypertension == "N" ~ 0,
#    hypertension == "Y" ~ 1
#  )) %>%
#  mutate(WIC = case_when(
#    WIC == "N" ~ 0,
#    WIC == "Y" ~ 1
#  ))
#write.csv(data,file="/home/shushuz/M_rock/data.csv")

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
write.csv(data,file="/home/shushuz/M_rock/application_data.csv")
design_matrix = model.matrix(~race_mom+prenatal_visit_category + diabete + hypertension + cigarette + age_cat +  WIC + marital_cat + education_cat,data=data)
data_design <- cbind(data$birthweight,design_matrix[,-1])
#print(head(data_design))
write.csv(data_design,file="/home/shushuz/M_rock/design_matrix_main.csv")

if(FALSE){
print("overall")
#data = read.csv("/home/shushuz/M_rock/data.csv")
print(summary(lm(prenatal_visit~age+diabete+hypertension+cigarette+WIC+education+marital,data=data)))
print("white")
data_white = data[data$race==1,]
print(summary(lm(prenatal_visit~age+diabete+hypertension+cigarette+WIC+education+marital,data=data_white)))
print("black")
data_black = data[data$race==2,]
print(summary(lm(prenatal_visit~age+diabete+hypertension+cigarette+WIC+education+marital,data=data_black)))
print("asian")
data_asian = data[data$race==4,]
print(summary(lm(prenatal_visit~age+diabete+hypertension+cigarette+WIC+education+marital,data=data_asian)))
print("hispanic")
data_hispanic = data[data$race==7,]
print(summary(lm(prenatal_visit~age+diabete+hypertension+cigarette+WIC+education+marital,data=data_hispanic)))
#print(apply(data,2,function(x) sum(is.na(x))))
#print(summary(data))
#corr_white = cor(data[data$race_mom==1,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital")])
#corr_black = cor(data[data$race_mom==2,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital")])
#corr_asian = cor(data[data$race_mom==4,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital")])
#corr_hispanic = cor(data[data$race_mom==7,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital")])
print("0")
  data_10 = data[data$prenatal_visit==0,]
print(nrow(data_10[data_10$race_mom==1,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital","education")]))
print(nrow(data_10[data_10$race_mom==2,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital","education")]))
print(nrow(data_10[data_10$race_mom==4,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital","education")]))
print(nrow(data_10[data_10$race_mom==7,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital","education")]))

print("<=5")
data_10 = data[data$prenatal_visit<=5,]
print(nrow(data_10[data_10$race_mom==1,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital","education")]))
print(nrow(data_10[data_10$race_mom==2,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital","education")]))
print(nrow(data_10[data_10$race_mom==4,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital","education")]))
print(nrow(data_10[data_10$race_mom==7,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital","education")]))

print("6-10")
  data_10 = data[(data$prenatal_visit>5)&(data$prenatal_visit<=10),]
print(nrow(data_10[data_10$race_mom==1,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital","education")]))
print(nrow(data_10[data_10$race_mom==2,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital","education")]))
print(nrow(data_10[data_10$race_mom==4,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital","education")]))
print(nrow(data_10[data_10$race_mom==7,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital","education")]))

print(">10")
data_10 = data[data$prenatal_visit>10,]
print(nrow(data_10[data_10$race_mom==1,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital","education")]))
print(nrow(data_10[data_10$race_mom==2,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital","education")]))
print(nrow(data_10[data_10$race_mom==4,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital","education")]))
print(nrow(data_10[data_10$race_mom==7,c("prenatal_visit","diabete","hypertension","cigarette","age","WIC","marital","education")]))

}
