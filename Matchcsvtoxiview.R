#MeroxAdjust #
#Date: 06-01-2022
#Author: Victor Chrone
#Import fucntions

setwd("C:\\Users\\victo\\Desktop\\speciale\\R\\MEROXSORT\\MeroxAdjust")
library(plyr)
library(stringr)
source("sort_merox_functions_v2.R")


RTC4<-MeroxAdjust(RTC4)
RTC4D<-MeroxAdjust(RTC4D)
RTC_TEMP<-MeroxAdjust(RTC_TEMP)
RTC_OVERNIGHT<-MeroxAdjust(RTC_OVERNIGHT)
RTC_FILTER<-MeroxAdjust(RTC_FILTER)
RTC_HOMEOTYPIC<-MeroxAdjust(RTC_HOMEOTYPIC)
RTCINF<-MeroxAdjust(RTCINF)
RTCINF30<-MeroxAdjust(RTCINF30)
RTCINF250<-MeroxAdjust(RTCINF250)
RTC2<-MeroxAdjust(RTC2)
RTCab<-MeroxAdjust(RTCab)
RTCabcd<-MeroxAdjust(RTCabcd)
RTC_SI<-MeroxAdjust(RTC_SI)
RTC_SO<-MeroxAdjust(RTC_SO)
RTC_ST<-MeroxAdjust(RTC_ST)
RTC410_40<-MeroxAdjust(RTC410_40)
RTC_PG<-MeroxAdjust(RTC_PG)
RTC_IGA<-MeroxAdjust(RTC_IGA)

#import csv
setwd("C:\\Users\\victo\\Desktop\\speciale\\Data\\MeroX\\5-10-21\\xiview")
setwd("C:\\Users\\victo\\Desktop\\speciale\\Data\\MeroX\\16-12-21\\xiview")
setwd("C:\\Users\\victo\\Desktop\\speciale\\Data\\MeroX\\25-03-22\\xiview")
ocre1<-read.csv("20211009-1534_EXP3_1000227_20211005_PH_VGC_PH_IgG_39-42_XL.csv_MeroxAdjust.csv",header = T)
ocre2<-read.csv("20211009-2120_EXP3_1000232_20211005_PH_VGC_PH_IgG_43-46_XL.csv_MeroxAdjust.csv",header = T)
ocre3<-read.csv("20211011-1527_EXP3_1000243_20211005_PH_VGC_PH_IgG_47-50_XL.csv_MeroxAdjust.csv",header = T)
ocre<-rbind(ocre1,ocre2,ocre3)



ocrelizumab<-tysabri

ocre$conpep<-str_c(str_c("[",ocre$PepSeq1,"]",""),str_c("[",ocre$PepSeq2,"]",""),"")
ocre$conresid<-str_c(str_c((as.numeric(ocre$PepPos1)+as.numeric(ocre$LinkPos1)-1),(as.numeric(ocre$PepPos2)+as.numeric(ocre$LinkPos2)-1),sep = "|"),str_c(str_c(">",ocre$Protein1,sep=""),str_c(">",ocre$Protein2,sep=""),sep="|"),sep = "|")

ocre$corrected<-"False"
for(i in 1:dim(ocre)[1]){
  for(j in 1:dim(ocrelizumab)[1]){
    if(ocre$conpep[i]==ocrelizumab$concpep[j]){
      ocre$corrected[i] = "True"
      if(length(strsplit(ocrelizumab$best.linkage.position.peptide.1[j],"")[[1]])==2){
        ocre$LinkPos1[i] =strsplit(ocrelizumab$best.linkage.position.peptide.1[j],"")[[1]][2]
      }else{ocre$LinkPos1[i] =str_c(strsplit(ocrelizumab$best.linkage.position.peptide.1[j],"")[[1]][2:3],"")
      }
      if(length(strsplit(ocrelizumab$best.linkage.position.peptide.2[j],"")[[1]])==2){
        ocre$LinkPos2[i] =strsplit(ocrelizumab$best.linkage.position.peptide.2[j],"")[[1]][2]
      }else{ocre$LinkPos2[i] =str_c(strsplit(ocrelizumab$best.linkage.position.peptide.2[j],"")[[1]][2:3],"")
      }
    }
}
}
`%!in%` <- Negate(`%in%`)

#first filtration
forxiview<-filter(ocre, conresid %in% ocrelizumab$conresid)
unique(forxiview$conpep)

ocrelizumab<-filter(ocrelizumab,conresid %!in% forxiview$conresid)

#part 2 to get extra
ocrelizumab$Peptide.11<-"hej"
ocrelizumab$Peptide22<-ocrelizumab$Peptide2

for(i in 1:dim(ocrelizumab)[1]){
  pep<-str_split(ocrelizumab$Peptide.1[i],"")
  
  for(j in 1:(length(pep[[1]])-1)){
    if(pep[[1]][j]==pep[[1]][j+1]){
      pep[[1]][j]<-NA
    }
  }
  pep<-na.omit(pep[[1]])
  pep<-as.character(pep)
  res <- rle(pep)
  compressedString <- paste(res$values, collapse = "", sep = "")
  ocrelizumab$Peptide.11[i]<-compressedString
}
for(i in 1:dim(ocrelizumab)[1]){
  pep2<-str_split(ocrelizumab$Peptide2[i],"")
  if(ocrelizumab$Peptide2[i] != 0){
    if(ocrelizumab$Peptide2[i] != 1){
  for(j in 1:(length(pep2[[1]])-1)){
    if(pep2[[1]][j]==pep2[[1]][j+1]){
      pep2[[1]][j]<-NA
    }
  }
  pep2<-na.omit(pep2[[1]])
  pep2<-as.character(pep2)
  res2 <- rle(pep2)
  compressedString2 <- paste(res2$values, collapse = "", sep = "")
  ocrelizumab$Peptide22[i]<-compressedString2
    }}
  
}
ocrelizumab$conpep2<-str_c(ocrelizumab$Peptide.1,ocrelizumab$Peptide2,"")

#second filtration
forxiview2<-filter(ocre, conpep %in% ocrelizumab$conpep2)

forxiview<-rbind(forxiview,forxiview2)

#thrid filtration
ocrelizumab<-filter(ocrelizumab,conpep2 %!in% forxiview$conpep)
ocrelizumab<-filter(ocrelizumab,Peptide22 != 0)
for(i in 1:dim(ocrelizumab)[1]){
  if(ocrelizumab$Peptide2[i]==1){
    ocrelizumab$Peptide22[i] <-ocrelizumab$Peptide.11[i]
  }
}

ocrelizumab$conpep2<-str_c(ocrelizumab$Peptide.11,ocrelizumab$Peptide22,"")

forxiview3<-filter(ocre, conpep %in% ocrelizumab$conpep2)

forxiview<-rbind(forxiview,forxiview3)

#fourth filtration
ocrelizumab<-filter(ocrelizumab,conpep2 %!in% forxiview$conpep)
ocrelizumab$conpep2<-str_c(ocrelizumab$Peptide.11,ocrelizumab$Peptide22,"")
forxiview4<-filter(ocre, conpep %in% ocrelizumab$conpep2)
forxiview<-rbind(forxiview,forxiview4)

#fith filtration
ocrelizumab<-filter(ocrelizumab,conpep2 %!in% forxiview$conpep)
forxiview5<-filter(ocre,ExpMz %in% ocrelizumab$m.z )
forxiview<-rbind(forxiview,forxiview5)

#sixth filtration
ocrelizumab<-filter(ocrelizumab,m.z %!in% forxiview5$ExpMz)
ocre$conpep3<-str_c(str_c("{",ocre$PepSeq1,"]",""),str_c("[",ocre$PepSeq2,"]",""),"")
ocre$conpep4<-str_c(str_c("[",ocre$PepSeq1,"]",""),str_c("{",ocre$PepSeq2,"]",""),"")
forxiview6<-filter(ocre,conpep3 %in% ocrelizumab$conpep2)
ocrelizumab<-filter(ocrelizumab,conpep2 %!in% forxiview6$conpep3)

forxiview7<-filter(ocre,conpep4 %in% ocrelizumab$conpep2 )
ocrelizumab<-filter(ocrelizumab,conpep2 %!in% forxiview7$conpep4)

forxiview<-rbind(forxiview,forxiview6[,1:27])
forxiview<-rbind(forxiview,forxiview7[,1:27])

forxiview<-forxiview[,1:25]
setwd("C:\\Users\\victo\\Desktop")
write.csv(x = forxiview,file = "filtered_IgG_tys.csv")
write.csv(x = ocre,file = "all_IgG_nat.csv")
#write.csv(x = ocre,file = "ocre.csv")



tys1<-read.csv("20211009-1534_EXP3_1000227_20211005_PH_VGC_PH_IgG_39-42_XL.csv",header = T)
tys2<-read.csv("20211009-2120_EXP3_1000232_20211005_PH_VGC_PH_IgG_43-46_XL.csv",header = T)
tys3<-read.csv("20211011-1527_EXP3_1000243_20211005_PH_VGC_PH_IgG_47-50_XL.csv",header = T)
tys<-rbind(tys1,tys2,tys3)

tys$conpep<-str_c(str_c("[",tys$PepSeq1,"]",""),str_c("[",tys$PepSeq2,"]",""),"")
tys$corrected<-""
tys$corrected<-"False"
for(i in 1:dim(tys)[1]){
  for(j in 1:dim(tysabri)[1]){
    if(tys$conpep[i]==tysabri$concpep[j]){
      tys$corrected[i] = "True"
      if(length(strsplit(tysabri$best.linkage.position.peptide.1[j],"")[[1]])==2){
        tys$LinkPos1[i] =strsplit(tysabri$best.linkage.position.peptide.1[j],"")[[1]][2]
      }else{tys$LinkPos1[i] =str_c(strsplit(tysabri$best.linkage.position.peptide.1[j],"")[[1]][2:3],"")
      }
      if(length(strsplit(tysabri$best.linkage.position.peptide.2[j],"")[[1]])==2){
        tys$LinkPos2[i] =strsplit(tysabri$best.linkage.position.peptide.2[j],"")[[1]][2]
      }else{tys$LinkPos2[i] =str_c(strsplit(tysabri$best.linkage.position.peptide.2[j],"")[[1]][2:3],"")
      }
    }
  }
}

forxiview<-filter(tys, conpep %in% tysabri$concpep)
#forxiview<-filter(ocre, corrected == "True")
unique(forxiview$conpep)
forxiview<-forxiview[,1:25]
setwd("C:\\Users\\victo\\Desktop")
write.csv(x = forxiview,file = "filtered_IgGtys.csv")
write.csv(x = tys,file = "tys.csv")
