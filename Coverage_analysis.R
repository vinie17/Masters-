library(readxl)
library(ggvenn)
library(stringr)


setwd("C:\\Users\\victo\\Desktop\\speciale\\Data\\gpmaw\\5-10-21")
getwd()
RTC4<-list()
RTC2<-list()
#IMPORT XLS data from gpmaw
files<-c("20211005-1756_EXP3_1000147_20211005_PH_VGC_PH_IgG_01-04.xls",
            "20211005-2340_EXP3_1000152_20211005_PH_VGC_PH_IgG_05-8.xls",
            "20211006-1401_EXP3_1000160_20211005_PH_VGC_PH_IgG_09-12.xls",
            "20211006-1946_EXP3_1000165_20211005_PH_VGC_PH_IgG_13-16.xls",
            "20211007-0130_EXP3_1000170_20211005_PH_VGC_PH_IgG_17-20.xls",
            "20211007-0716_EXP3_1000175_20211005_PH_VGC_PH_IgG_21-24.xls",
            "20211008-2033_EXP3_1000209_20211005_PH_VGC_PH_IgG_26-29.xls",
            "20211009-0219_EXP3_1000214_20211005_PH_VGC_PH_IgG_30-33.xls",
            "20211009-0805_EXP3_1000219_20211005_PH_VGC_PH_IgG_34-37.xls",
            "20211009-1534_EXP3_1000227_20211005_PH_VGC_PH_IgG_39-42.xls",
            "20211009-2120_EXP3_1000232_20211005_PH_VGC_PH_IgG_43-46.xls",
            "20211011-1527_EXP3_1000243_20211005_PH_VGC_PH_IgG_47-50.xls"
            )
for(i in 1:6){
RTC4[[i]]<-read_xls(files[i+6],col_names = F)
RTC2[[i]]<-read_xls(files[i],col_names = F)
}
set.seed(321)


#Count how many times each peptide is seen
#remove NA

for(i in 1:length(RTC4)){
  counts<-aggregate(...2~...11, data = RTC4[[i]], summary)
  for(j in 1:dim(counts)[1]){
     if(counts$...2[j] == "1"){
       RTC4[[i]]<-RTC4[[i]][!is.na(RTC4[[i]]$...11), ]
       for(k in 1:dim(RTC4[[i]][1])){
         if(RTC4[[i]]$...11[k] == counts$...11[j]){
           RTC4[[i]]$...11[k]<- NA
         }
       }
     }
  }
}
RTC4PEP<-RTC4
# Count the number of missed cleavages
for(i in 1:length(RTC4)){
  RTC4[[i]]$number.of.r <- str_count(RTC4[[i]]$...11, "R")
  RTC4[[i]]$number.of.k <- str_count(RTC4[[i]]$...11, "K")
  RTC4[[i]]$missedcleavage <- (RTC4[[i]]$number.of.r+RTC4[[i]]$number.of.k)
}
 
#Remove missed cleavages

for(i in 1:length(RTC4)){
  for(j in 1:dim(RTC4[[i]])[1]){
    
    try(if(RTC4[[1]]$missedcleavage[j] == 3){
      RTC4[[i]]$missedcleavage[j] <- NA
    },silent = T)
    
  }
  RTC4[[i]]<-RTC4[[i]][!is.na(RTC4[[i]]$missedcleavage), ]
}



#C4
O<-list()
O<- list(
  A = RTC4[[1]]$...11, 
  B = RTC4[[2]]$...11, 
  C = RTC4[[3]]$...11
)        

TY<-list()
TY <- list(
  A = RTC4[[4]]$...11, 
  B = RTC4[[5]]$...11, 
  C = RTC4[[6]]$...11
)



#C2
O<-list()
O<- list(
  A = RTC2[[1]]$...11, 
  B = RTC2[[2]]$...11, 
  C = RTC2[[3]]$...11
)

TY<-list()
TY <- list(
  A = RTC2[[4]]$...11, 
  B = RTC2[[5]]$...11, 
  C = RTC2[[6]]$...11
)
#Ven diagram
ggvenn(
  O, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 1, set_name_size = 3
)


