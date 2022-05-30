
#open structure #####
open<-filter(keepsafe,model=="open")

for(i in 1:dim(open)[1]){
  if(open$type1[i]=="A"){
    open$type1[i]=">light"}
  if(open$type2[i]=="A"){
    open$type2[i]=">light"}
  if(open$type1[i]=="C"){
    open$type1[i]=">light"}
  if(open$type2[i]=="C"){
    open$type2[i]=">light"}
  if(open$type1[i]=="B"){
    open$type1[i]=">heavy"}
  if(open$type2[i]=="B"){
    open$type2[i]=">heavy"}
  if(open$type1[i]=="D"){
    open$type1[i]=">heavy"}
  if(open$type2[i]=="D"){
    open$type2[i]=">heavy"}
}

open$conpep<-str_c(open$Lys1,open$Lys2,open$type1,open$type2,sep = "|")

setwd("C:\\Users\\victo\\Desktop\\speciale\\Data\\MeroX\\5-10-21\\xiview")
meroxcsv<-read.csv("20211008-2033_EXP3_1000209_20211005_PH_VGC_PH_IgG_26-29_XL.csv")
meroxcsv$newpos1<-""
meroxcsv$newpos2<-""
for(i in 1:dim(meroxcsv)[1]){
  if(meroxcsv$Protein1[i]=="light"){
    meroxcsv$newpos1[i]=">light"
  }
  if(meroxcsv$Protein2[i]=="light"){
    meroxcsv$newpos2[i]=">light"
  }
  if(meroxcsv$Protein1[i]=="heavy"){
    meroxcsv$newpos1[i]=">heavy"
  }
  if(meroxcsv$Protein2[i]=="heavy"){
    meroxcsv$newpos2[i]=">heavy"
  }
}
meroxcsv$lys1<-(meroxcsv$PepPos1+meroxcsv$LinkPos1-1)
meroxcsv$lys2<-(meroxcsv$PepPos2+meroxcsv$LinkPos2-1)

meroxcsv$conpep<-str_c(meroxcsv$lys1,meroxcsv$lys2,meroxcsv$newpos1,meroxcsv$newpos2,sep = "|")

openfinal<-meroxcsv[meroxcsv$conpep %in% open$conpep,]
openfinal<-openfinal[,1:25]
#openfinal<-checkrev(test)
#test2<-table(test$conresid)  
#ocrelizmab_33_resid_count<-ldply(test2,data.frame)
#ocrelizmab_33_merox<-test[!duplicated(openfinal$conresid),]
setwd("C:\\Users\\victo\\Desktop")
write.csv(x = openfinal,file = "open.csv")


#closed structure #####
closed<-filter(keepsafe,model=="closed")

for(i in 1:dim(closed)[1]){
  if(closed$type1[i]=="A"){
    closed$type1[i]=">light"}
  if(closed$type2[i]=="A"){
    closed$type2[i]=">light"}
  if(closed$type1[i]=="C"){
    closed$type1[i]=">light"}
  if(closed$type2[i]=="C"){
    closed$type2[i]=">light"}
  if(closed$type1[i]=="B"){
    closed$type1[i]=">heavy"}
  if(closed$type2[i]=="B"){
    closed$type2[i]=">heavy"}
  if(closed$type1[i]=="D"){
    closed$type1[i]=">heavy"}
  if(closed$type2[i]=="D"){
    closed$type2[i]=">heavy"}
}

closed$conpep<-str_c(closed$Lys1,closed$Lys2,closed$type1,closed$type2,sep = "|")

setwd("C:\\Users\\victo\\Desktop\\speciale\\Data\\MeroX\\5-10-21\\xiview")
meroxcsv<-read.csv("20211008-2033_EXP3_1000209_20211005_PH_VGC_PH_IgG_26-29_XL.csv")
meroxcsv$newpos1<-""
meroxcsv$newpos2<-""
for(i in 1:dim(meroxcsv)[1]){
  if(meroxcsv$Protein1[i]=="light"){
    meroxcsv$newpos1[i]=">light"
  }
  if(meroxcsv$Protein2[i]=="light"){
    meroxcsv$newpos2[i]=">light"
  }
  if(meroxcsv$Protein1[i]=="heavy"){
    meroxcsv$newpos1[i]=">heavy"
  }
  if(meroxcsv$Protein2[i]=="heavy"){
    meroxcsv$newpos2[i]=">heavy"
  }
}
meroxcsv$lys1<-(meroxcsv$PepPos1+meroxcsv$LinkPos1-1)
meroxcsv$lys2<-(meroxcsv$PepPos2+meroxcsv$LinkPos2-1)

meroxcsv$conpep<-str_c(meroxcsv$lys1,meroxcsv$lys2,meroxcsv$newpos1,meroxcsv$newpos2,sep = "|")

closedfinal<-meroxcsv[meroxcsv$conpep %in% closed$conpep,]
closedfinal<-closedfinal[,1:25]
#openfinal<-checkrev(test)
#test2<-table(test$conresid)  
#ocrelizmab_33_resid_count<-ldply(test2,data.frame)
#ocrelizmab_33_merox<-test[!duplicated(openfinal$conresid),]
setwd("C:\\Users\\victo\\Desktop")
write.csv(x = closedfinal,file = "closed.csv")





#once the open and closed structures have been made #####
#For ocrelizumab open and closed
setwd("C:\\Users\\victo\\Desktop\\speciale\\Data\\OpenClosed structure")
setwd("C:\\Users\\victo\\Desktop\\speciale\\Data\\Anti- murine CD163")
open<-read.csv("open.csv")
closed<-read.csv("closed.csv")
open$newpos1<-""
open$newpos2<-""
closed$newpos1<-""
closed$newpos2<-""

for(i in 1:dim(open)[1]){
  if(open$Protein1[i]=="light"){
    open$newpos1[i]=">light"
  }
  if(open$Protein2[i]=="light"){
    open$newpos2[i]=">light"
  }
  if(open$Protein1[i]=="heavy"){
    open$newpos1[i]=">heavy"
  }
  if(open$Protein2[i]=="heavy"){
    open$newpos2[i]=">heavy"
  }
}
for(i in 1:dim(closed)[1]){
  if(closed$Protein1[i]=="light"){
    closed$newpos1[i]=">light"
  }
  if(closed$Protein2[i]=="light"){
    closed$newpos2[i]=">light"
  }
  if(closed$Protein1[i]=="heavy"){
    closed$newpos1[i]=">heavy"
  }
  if(closed$Protein2[i]=="heavy"){
    closed$newpos2[i]=">heavy"
  }
}
open$lys1<-(open$PepPos1+open$LinkPos1-1)
open$lys2<-(open$PepPos2+open$LinkPos2-1)
closed$lys1<-(closed$PepPos1+closed$LinkPos1-1)
closed$lys2<-(closed$PepPos2+closed$LinkPos2-1)

open$conpep<-str_c(open$lys1,open$lys2,open$newpos1,open$newpos2,sep = "|")
closed$conpep<-str_c(closed$lys1,closed$lys2,closed$newpos1,closed$newpos2,sep = "|")
binded<-ldply(RTC_SO[5],data.frame)

closedfinal5<-binded[binded$conresid %in% closed$conpep,]

openfinal5<-binded[binded$conresid %in% open$conpep,]


open1<-c(dim(openfinal11)[1],dim(openfinal12)[1],dim(openfinal13)[1],dim(openfinal14)[1],dim(openfinal15)[1])
closed1<-c(dim(closedfinal11)[1],dim(closedfinal12)[1],dim(closedfinal13)[1],dim(closedfinal14)[1],dim(closedfinal15)[1])


ratio<-data.frame(open1,closed1)
ratio$ratio<-open1/closed1
<
SLICESratio<-ratio
#openfinal<-openfinal[,1:25]
#openfinal<-checkrev(test)
#test2<-table(test$conresid)  
#ocrelizmab_33_resid_count<-ldply(test2,data.frame)
#ocrelizmab_33_merox<-test[!duplicated(openfinal$conresid),]
setwd("C:\\Users\\victo\\Desktop")
write.csv(x = closedfinal,file = "closed.csv")

test2<-table(openfinal5$conresid)
test2<-ldply(test2,data.frame)
test2<-filter(test2,X..i..>15)

openfinal15<-openfinal5[openfinal5$conresid %in% test2$.id,]
min30final<-min30final[,1:25]
setwd("C:\\Users\\victo\\Desktop")
write.csv(x = min30final,file = "min30.csv")

