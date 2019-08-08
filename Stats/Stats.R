library(R.matlab)

path_code = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Stats'

path_data_1 = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/Model_6_10_CoV_50_Ur_Rec_3/'
path_data_2 = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/Model_6_10_CoV_50_Ur_Rec_3_shortTendon/'

setwd(path_code)

#============================================================
data.temp <- readMat(paste(path_data_1,'cov_Force.mat',sep = ""))
data_1 <- data.temp$cov.Force

data.temp <- readMat(paste(path_data_2,'cov_Force.mat',sep = ""))
data_2 <- data.temp$cov.Force

p <-numeric(dim(data_1)[2])

for(i in 1:dim(data_1)[2]){ #for each simulated experiment
  z <- t.test(data_1[,i], data_2[,i],alternative = "two.sided", var.equal = FALSE)
  p[i]<-z$p.value #get the p-value and store it
}

round(p.adjust(p, method = "BH"),3)
