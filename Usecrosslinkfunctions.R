#set working directory to active file location
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path)) 
getwd()
#Import fucntions
source("crosslinkFunctionsV2.R")

setwd("C:\\Users\\victo\\Desktop\\speciale\\Data\\Modelsforfigures")
setwd("C:\\Users\\victo\\Desktop\\speciale\\Data\\OpenClosed structure")
setwd("C:\\Users\\victo\\Desktop\\speciale\\Data\\ocrelizumab\\model")
setwd("C:\\Users\\victo\\Desktop\\speciale\\Data\\ocrelizumab\\links")
setwd("C:\\Users\\victo\\Desktop\\speciale\\Data\\Tysabri\\Model\\Sheila")
#open2<-read.csv("open_strict.csv",header=T)

#import pdbfile
linkImporter(
  file = "openmodel2.pdb",  #pdb file name located in dir or txt file of links
  modeltypes = "Model 2.1",              #model used for distance generated
  filenames = "closed", #filename without pdb
  links = "shortest",                        #all possible links or shortest
  antibodys = "tys",                   #which antibody such as ocre or tys or iga
  directlyfromPDB = TRUE,
  correctcordinates = T)               #file being uploaded is pdb?

rm(distanalysis)
rm(distanalysis1)
rm(keepsafe)
#generate plots
crosslinkplotspecific(                 
  data = keepsafe,                 #data being used. recomend keepsafe or a filtered version
  filteringafter = "tys",              #filtering after tys ocre or igA
  IgGIgA = "IgG",                      #is it IgG or IgA
  lightchain =  lightChain,            #Single letter sequence of light chain
  heavychain =  heavyChain,            #Single letter sequence of heavy chain
  jitx = 2,                           #How much jitter for scatter plot 
  jity = 2,
  multiple = F)                           #same as above 

#classfication(crosslinkdist,"ocre")
#colorclass(distanalysis1)
#chainsplit(distanalysis1,"IgG",lightChain,heavyChain)
#crosslinkplot(scatterplotdist,lightChain,heavyChain,2,2)

#Import csv file with data for linkImporter
import_from_csv("linkdataimport.xlsx")

#distribution plot
linkdist
#scatter plots

linkscatter

#open structure #####
open<-filter(keepsafe,model=="open")

for(i in 1:dim(open)[1]){
  if(open$type1[i]=="A"){
    open$type1[i]=">E10B10-L"}
  if(open$type2[i]=="A"){
    open$type2[i]=">E10B10-L"}
  if(open$type1[i]=="C"){
    open$type1[i]=">E10B10-L"}
  if(open$type2[i]=="C"){
    open$type2[i]=">E10B10-L"}
  if(open$type1[i]=="B"){
    open$type1[i]=">E10B10-H"}
  if(open$type2[i]=="B"){
    open$type2[i]=">E10B10-H"}
  if(open$type1[i]=="D"){
    open$type1[i]=">E10B10-H"}
  if(open$type2[i]=="D"){
    open$type2[i]=">E10B10-H"}
  }

open$conpep<-str_c(open$Lys1,open$Lys2,open$type1,open$type2,sep = "|")

setwd("C:\\Users\\victo\\Desktop\\speciale\\Data\\MeroX\\5-10-21\\xiview")
setwd("C:\\Users\\victo\\Desktop\\speciale\\Data\\MeroX\\25-03-22\\xiview")
meroxcsv<-read.csv("20220325-2100_EXP3_1001571_PH_VGC_20220324_Ab1cd1_XL.csv")
meroxcsv$newpos1<-""
meroxcsv$newpos2<-""
for(i in 1:dim(meroxcsv)[1]){
  if(meroxcsv$Protein1[i]=="E10B10-L"){
    meroxcsv$newpos1[i]=">E10B10-L"
  }
  if(meroxcsv$Protein2[i]=="E10B10-L"){
    meroxcsv$newpos2[i]=">E10B10-L"
  }
  if(meroxcsv$Protein1[i]=="E10B10-H"){
    meroxcsv$newpos1[i]=">E10B10-H"
  }
  if(meroxcsv$Protein2[i]=="E10B10-H"){
    meroxcsv$newpos2[i]=">E10B10-H"
  }
}
meroxcsv$lys1<-(meroxcsv$PepPos1+meroxcsv$LinkPos1-1)
meroxcsv$lys2<-(meroxcsv$PepPos2+meroxcsv$LinkPos2-1)

meroxcsv$conpep<-str_c(meroxcsv$lys1,meroxcsv$lys2,meroxcsv$newpos1,meroxcsv$newpos2,sep = "|")

openfinal<-meroxcsv[meroxcsv$conpep %in% open$conpep,]
#openfinal<-openfinal[,1:25]
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
    closed$type1[i]=">E10B10-L"}
  if(closed$type2[i]=="A"){
    closed$type2[i]=">E10B10-L"}
  if(closed$type1[i]=="C"){
    closed$type1[i]=">E10B10-L"}
  if(closed$type2[i]=="C"){
    closed$type2[i]=">E10B10-L"}
  if(closed$type1[i]=="B"){
    closed$type1[i]=">E10B10-H"}
  if(closed$type2[i]=="B"){
    closed$type2[i]=">E10B10-H"}
  if(closed$type1[i]=="D"){
    closed$type1[i]=">E10B10-H"}
  if(closed$type2[i]=="D"){
    closed$type2[i]=">E10B10-H"}
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

test2<-table(meroxcsv$conpep)
test2<-ldply(test2,data.frame)
test2<-filter(test2,X..i..>15)

min30final<-meroxcsv[meroxcsv$conpep %in% test2$.id,]
min30final<-min30final[,1:25]
setwd("C:\\Users\\victo\\Desktop")
write.csv(x = min30final,file = "min30.csv")
