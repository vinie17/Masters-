#set seed
set.seed(123)
####### For Anti-Murine CD163  ####
setwd("C:\\Users\\victo\\Desktop\\speciale\\R\\MEROXSORT\\MeroxAdjust")
E10B10<-read_xlsx("final_E10B10_dataset_type0(+)3.xlsx")
E10B10CD163<-read_xlsx("final_E10B10+CD163_dataset_type0(+)3.xlsx")
dataset<-list(E10B10,E10B10CD163)
linklist<-c("E10B10","E10B10 + CD163")

RTCab<-dataset



conpepres2(RTCab)
dataset<-listdata

conpepres2(RTCabcd)
RTCabcd<-listdata
unique(E10B10CD163$Protein.1)
for(i in 1:length(RTCab)){
  filt1<-filter(RTCab[[i]],Protein.1 %in% c(">E10B10-H",">E10B10-L"))
  filt2<-filter(filt1,Protein.2 %in% c(">E10B10-H",">E10B10-L","intrapeptidal","NH3","H2O"))
  RTCab[[i]]<-filt2
}

for(i in 1:length(RTCabcd)){
  filt1<-filter(RTCabcd[[i]],Protein.1 %in% c(">E10B10-H",">E10B10-L"))
  filt2<-filter(filt1,Protein.2 %in% c(">E10B10-H",">E10B10-L","intrapeptidal","NH3","H2O"))
  RTCabcd[[i]]<-filt2
}


for(i in 1:length(RTCab)){
  for(j in 1:dim(RTCab[[i]])[1]){
    if(RTCab[[i]]$Protein.1[j]==">sp|Q2VLH6|C163A_MOUSE Scavenger receptor cysteine-rich type 1 protein M130 OS=Mus musculus OX=10090 GN=Cd163 PE=1 SV=2"){
      RTCab[[i]]$Protein.1[j] = "CD163"
      print("replaced")
    }
    if(RTCab[[i]]$Protein.2[j]==">sp|Q2VLH6|C163A_MOUSE Scavenger receptor cysteine-rich type 1 protein M130 OS=Mus musculus OX=10090 GN=Cd163 PE=1 SV=2"){
      RTCab[[i]]$Protein.2[j] = "CD163"
      print("replaced 2")
    }
  }
  
}


for(i in 1:length(RTCabcd)){
  for(j in 1:dim(RTCabcd[[i]])[1]){
    if(RTCabcd[[i]]$Protein.1[j]==">sp|Q2VLH6|C163A_MOUSE Scavenger receptor cysteine-rich type 1 protein M130 OS=Mus musculus OX=10090 GN=Cd163 PE=1 SV=2"){
      RTCabcd[[i]]$Protein.1[j] = "CD163"
      print("replaced")
    }
    if(RTCabcd[[i]]$Protein.2[j]==">sp|Q2VLH6|C163A_MOUSE Scavenger receptor cysteine-rich type 1 protein M130 OS=Mus musculus OX=10090 GN=Cd163 PE=1 SV=2"){
      RTCabcd[[i]]$Protein.2[j] = "CD163"
      print("replaced 2")
    }
  }
  
}





dataset<-list(RTCab[[1]],RTCab[[2]],RTCab[[3]],RTCabcd[[1]],RTCabcd[[2]],RTCabcd[[3]],RTCabcd[[4]])
linklist<-c("E10B10 1","E10B10 2","E10B10 3","E10B10 + 2CD163","E10B10 + 1CD163_1","E10B10 + 1CD163_2","E10B10 + 1CD163_3")

linktypelistzero<-list(NA,NA,NA,NA,NA,NA)
linktypelistone<-list(NA,NA,NA,NA,NA,NA)
linktypelistrest<-list(NA,NA,NA,NA,NA,NA)

linktypelistzerores<-list(NA,NA,NA,NA,NA,NA)
linktypelistoneres<-list(NA,NA,NA,NA,NA,NA)
linktypelistrestres<-list(NA,NA,NA,NA,NA,NA)


crosslinkclassifyresid2(dataset = dataset,
                        linklist = linklist,
                        linktypelistzero,
                        linktypelistone,
                        linktypelistrest)

#filter linktypelistrestres filtering 
linktypelisthomeotypic<-list(NA,NA,NA,NA,NA,NA)
datares<-list(NA,NA,NA,NA,NA)
homeotypiclist<-list(NA,NA,NA,NA,NA,NA)

checkrest(linktypelistrestres)

linktypelistrestres<-datares  

NAremoveres(linktypelistrestres,linktypelisthomeotypic)

#filter linktypelistoneres filtering 
linktypelisthomeotypicone<-list(NA,NA,NA,NA,NA,NA,NA)
dataone<-list(NA,NA,NA,NA,NA,NA,NA)
homeotypiclistone<-list(NA,NA,NA,NA,NA,NA,NA)

checkones(linktypelistoneres)

linktypelistoneres<-dataone

NAremoveresone(linktypelistoneres,linktypelisthomeotypicone)


#filter the homeotypic peptides from original data
#can be used for either rest(homeotypic) and one(to find S and Y links not filtered by meroxsort)
findhomeo(dataset,linktypelisthomeotypic)
for(i in 1:length(homeotypiclist)){
  homeotypiclist[[i]]<-filter(homeotypiclist[[i]], Peptide2 != 1)
  linktypelisthomeotypic[[i]]<-homeotypiclist[[i]]$conresid
}

#Run each time to remove false homeotypic (S and Y)
findhomeoone(dataset,linktypelisthomeotypicone)
for(i in 1:length(homeotypiclistone)){
  homeotypiclistone[[i]]<-filter(homeotypiclistone[[i]], Peptide2 == 1)
  linktypelisthomeotypicone[[i]]<-homeotypiclistone[[i]]$conresid
}

#combine remaining type one and type rest
oneandrest<-list()
oneandrest<-mapply(c,linktypelistoneres[1:2], linktypelistrestres[1:2], SIMPLIFY=FALSE)
homeo<-mapply(c,linktypelisthomeotypic,linktypelisthomeotypicone)
#barplot of crosslink types
p<-list(NA)
crosslinktypeplot(p,linktypesres)


#look at different types
#C4

#Used to filter data using T5
O2<- list(
  A = setmincount(na.omit(oneandrest[[1]]),2), 
  B = setmincount(na.omit(oneandrest[[2]]),2), 
  C = setmincount(na.omit(oneandrest[[3]]),2)
)

TY2 <- list(
  A = setmincount(na.omit(oneandrest[[5]]),2), 
  B = setmincount(na.omit(oneandrest[[6]]),2), 
  C = setmincount(na.omit(oneandrest[[7]]),2)
)
#Ven diagram
o<-list(
  "E10B10" = oneandrest[[1]],
  "E10B10 + CD163" = oneandrest[[2]]
)
forven<-list(TY2,O2)
forven<-list(o)

crosslinkven(forven)

#uniquelists<-list()
barplotsres<-list()
linktypeslistres<-list()
crosslinklistres<-list()
venslistres<-list()
savedatares(p,linktypesres,Crosslinkclasses,VENS,1)


joinedCD163<-c(secs[[1]]$`2`,secs[[2]]$`3`,secs[[3]]$`5`,secs[[4]]$`1`)
binded<-ldply(RTCabcd[2:4],data.frame)
joinedCD163<-binded[binded$conresid %in% joinedCD163,]

joinedE10B10<-c(secs[[5]]$`2`,secs[[6]]$`3`,secs[[7]]$`5`,secs[[8]]$`1`)
binded<-ldply(RTCab[1:3],data.frame)
joinedE10B10<-binded[binded$conresid %in% joinedE10B10,]
setwd("C:\\Users\\victo\\Desktop\\speciale\\R\\MEROXSORT\\MeroxAdjust")
writexl::write_xlsx(joinedE10B10,path = "final_E10B10_dataset_type0(+)3.xlsx")
writexl::write_xlsx(joinedCD163,path = "final_E10B10+CD163_dataset_type0(+)3.xlsx")



####### For PROTEIN G AND PROTEIN A ####
conpepres2(RTC_PG)
RTC_PG<-listdata


unique(RTC_PG[[6]]$Protein.2)
unique(RTC_PG[[6]]$Protein.1)

#Replace protein A id with new short
for(i in 1:length(RTC_PG)){
  for(j in 1:dim(RTC_PG[[i]])){
  if(RTC_PG[[i]]$Protein.1[j]==">tr|Q70AB8|Q70AB8_STAAU Protein A  Fragment  OS=Staphylococcus aureus OX=1280 GN=spa PE=1 SV=1" || 
     RTC_PG[[i]]$Protein.1[j]==">tr|Q70AB8|Q70AB8_STAAU Protein A  Fragment  OS=Staphylococcus aureus OX=1280 GN=spa PE=1 SV=1(>tr|Q70AB8|Q70AB8_STAAU Protein A  Fragment  OS=Staphylococcus aureus OX=1280 GN=spa PE=1 SV=1)"){
    RTC_PG[[i]]$Protein.1[j] = ">Protein A"
    
  }
    if(RTC_PG[[i]]$Protein.2[j]==">tr|Q70AB8|Q70AB8_STAAU Protein A  Fragment  OS=Staphylococcus aureus OX=1280 GN=spa PE=1 SV=1" || 
       RTC_PG[[i]]$Protein.2[j]==">tr|Q70AB8|Q70AB8_STAAU Protein A  Fragment  OS=Staphylococcus aureus OX=1280 GN=spa PE=1 SV=1(>tr|Q70AB8|Q70AB8_STAAU Protein A  Fragment  OS=Staphylococcus aureus OX=1280 GN=spa PE=1 SV=1)"){
      RTC_PG[[i]]$Protein.2[j] = ">Protein A"
      
    }
  }
}
#Replace protein G with new short
for(i in 1:length(RTC_PG)){
  for(j in 1:dim(RTC_PG[[i]])){
    if(str_detect(RTC_PG[[i]]$Protein.1[j],"protein G")){
      RTC_PG[[i]]$Protein.1[j] = ">Protein G"
    }
    if(str_detect(RTC_PG[[i]]$Protein.2[j],"protein G")){
      RTC_PG[[i]]$Protein.2[j] = ">Protein G"
    }
  }
}
    
    
    
    
  



#CONTAINS pA : 4,11
#Contains pG: 5(should be pA),6,8(pl),9,
#conresid 
dataset<-list(RTC4[[1]],RTC4[[2]],RTC4[[3]],RTC4[[4]],RTC4[[5]],RTC4[[6]])
linklist<-c("OA","OB","OC","TYA","TYB","TYC")

linktypelistzero<-list(NA,NA,NA,NA,NA,NA)
linktypelistone<-list(NA,NA,NA,NA,NA,NA)
linktypelistrest<-list(NA,NA,NA,NA,NA,NA)

linktypelistzerores<-list(NA,NA,NA,NA,NA,NA)
linktypelistoneres<-list(NA,NA,NA,NA,NA,NA)
linktypelistrestres<-list(NA,NA,NA,NA,NA,NA)


crosslinkclassifyresid2(dataset = dataset,
                        linklist = linklist,
                        linktypelistzero,
                        linktypelistone,
                        linktypelistrest)

#filter linktypelistrestres filtering 
linktypelisthomeotypic<-list(NA,NA,NA,NA,NA,NA)
datares<-list(NA,NA,NA,NA,NA)
homeotypiclist<-list(NA,NA,NA,NA,NA,NA)

checkrest(linktypelistrestres)

linktypelistrestres<-datares  

NAremoveres(linktypelistrestres,linktypelisthomeotypic)

#filter linktypelistoneres filtering 
linktypelisthomeotypicone<-list(NA,NA,NA,NA,NA,NA)
dataone<-list(NA,NA,NA,NA,NA)
homeotypiclistone<-list(NA,NA,NA,NA,NA,NA)

checkones(linktypelistoneres)

linktypelistoneres<-dataone

NAremoveresone(linktypelistoneres,linktypelisthomeotypicone)


#filter the homeotypic peptides from original data
#can be used for either rest(homeotypic) and one(to find S and Y links not filtered by meroxsort)
findhomeo(dataset,linktypelisthomeotypic)
for(i in 1:length(homeotypiclist)){
  homeotypiclist[[i]]<-filter(homeotypiclist[[i]], Peptide2 != 1)
  linktypelisthomeotypic[[i]]<-homeotypiclist[[i]]$conresid
}

#Run each time to remove false homeotypic (S and Y)
findhomeoone(dataset,linktypelisthomeotypicone)
for(i in 1:length(homeotypiclistone)){
  homeotypiclistone[[i]]<-filter(homeotypiclistone[[i]], Peptide2 == 1)
  linktypelisthomeotypicone[[i]]<-homeotypiclistone[[i]]$conresid
}

#combine remaining type one and type rest
oneandrest<-list()
oneandrest<-mapply(c,linktypelistoneres, linktypelistrestres, SIMPLIFY=FALSE)
homeo<-mapply(c,linktypelisthomeotypic,linktypelisthomeotypicone)
#barplot of crosslink types
p<-list(NA)
crosslinktypeplot(p,linktypesres)


#look at different types
#C4

#combined with one(homeotypic removed)
O<- list(
  A = na.omit(oneandrest[[1]]), 
  B = na.omit(oneandrest[[2]]), 
  C = na.omit(oneandrest[[3]])
)

TY <- list(
  A = na.omit(oneandrest[[4]]), 
  B = na.omit(oneandrest[[5]]), 
  C = na.omit(oneandrest[[6]])
)
#Ven diagram

forven<-list(TY3,O3)


crosslinkven(forven)

#uniquelists<-list()
barplotsres<-list()
linktypeslistres<-list()
crosslinklistres<-list()
venslistres<-list()
savedatares(p,linktypesres,Crosslinkclasses,VENS,1)


####### FOR Multimeric OCRELIZUMAB AND TYSABRI ########################
#concatenate peptide 1 and peptide 2 so we known where links are between
conpepres2(RTC4T)
RTC4T<-listdata

#conresid 
dataset<-list(RTC4T[[1]],RTC4T[[2]])
linklist<-c("Ocrelizumab","Tysabri")

linktypelistzero<-list(NA,NA,NA,NA,NA,NA)
linktypelistone<-list(NA,NA,NA,NA,NA,NA)
linktypelistrest<-list(NA,NA,NA,NA,NA,NA)

linktypelistzerores<-list(NA,NA,NA,NA,NA,NA)
linktypelistoneres<-list(NA,NA,NA,NA,NA,NA)
linktypelistrestres<-list(NA,NA,NA,NA,NA,NA)


crosslinkclassifyresid2(dataset = dataset,
                        linklist = linklist,
                        linktypelistzero,
                        linktypelistone,
                        linktypelistrest)

#filter linktypelistrestres filtering 
linktypelisthomeotypic<-list(NA,NA,NA,NA,NA,NA)
datares<-list(NA,NA,NA,NA,NA)
homeotypiclist<-list(NA,NA,NA,NA,NA,NA)

checkrest(linktypelistrestres)

linktypelistrestres<-datares  

NAremoveres(linktypelistrestres,linktypelisthomeotypic)

#filter linktypelistoneres filtering 
linktypelisthomeotypicone<-list(NA,NA,NA,NA,NA,NA)
dataone<-list(NA,NA,NA,NA,NA)
homeotypiclistone<-list(NA,NA,NA,NA,NA,NA)

checkones(linktypelistoneres)

linktypelistoneres<-dataone

NAremoveresone(linktypelistoneres,linktypelisthomeotypicone)


#filter the homeotypic peptides from original data
#can be used for either rest(homeotypic) and one(to find S and Y links not filtered by meroxsort)
findhomeo(dataset,linktypelisthomeotypic)
for(i in 1:length(homeotypiclist)){
  homeotypiclist[[i]]<-filter(homeotypiclist[[i]], Peptide2 != 1)
  linktypelisthomeotypic[[i]]<-homeotypiclist[[i]]$conresid
}

#Run each time to remove false homeotypic (S and Y)
findhomeoone(dataset,linktypelisthomeotypicone)
for(i in 1:length(homeotypiclistone)){
  homeotypiclistone[[i]]<-filter(homeotypiclistone[[i]], Peptide2 == 1)
  linktypelisthomeotypicone[[i]]<-homeotypiclistone[[i]]$conresid
}

#combine remaining type one and type rest
oneandrest<-list()
oneandrest<-mapply(c,linktypelistoneres, linktypelistrestres, SIMPLIFY=FALSE)
homeo<-mapply(c,linktypelisthomeotypic,linktypelisthomeotypicone)
#barplot of crosslink types
p<-list(NA)
crosslinktypeplot(p,linktypesres)


#look at different types
#C4

#combined with one(homeotypic removed)
O<- list(
  A = na.omit(oneandrest[[1]]), 
  B = na.omit(oneandrest[[2]]), 
  C = na.omit(oneandrest[[3]])
)

TY <- list(
  A = na.omit(oneandrest[[4]]), 
  B = na.omit(oneandrest[[5]]), 
  C = na.omit(oneandrest[[6]])
)

#Ven diagram

forven<-list(TY3,O3)


crosslinkven(forven)

#uniquelists<-list()
barplotsres<-list()
linktypeslistres<-list()
crosslinklistres<-list()
venslistres<-list()
savedatares(p,linktypesres,Crosslinkclasses,VENS,1)

####### FOR HOMEOTYPIC ANALYSIS ######################### ####

#concatenate peptide 1 and peptide 2 so we known where links are between
conpepres2(RTC_HOMEOTYPIC)
dataset<-listdata


linklist<-c("Natalizumab_light","Natalizumab_heavy","Natalizumab_heavy+light","Ocrelizumab_light","Ocrelizumab_heavy","Ocrelizumab_light+heavy")

linktypelistzero<-list(NA,NA,NA,NA,NA,NA)
linktypelistone<-list(NA,NA,NA,NA,NA,NA)
linktypelistrest<-list(NA,NA,NA,NA,NA,NA)
linktypelistzerores<-list(NA,NA,NA,NA,NA,NA)
linktypelistoneres<-list(NA,NA,NA,NA,NA,NA)
linktypelistrestres<-list(NA,NA,NA,NA,NA,NA)


crosslinkclassifyresid2(dataset = dataset,
                       linklist = linklist,
                       linktypelistzero,
                       linktypelistone,
                       linktypelistrest)


#filter linktypelistrestres filtering 
linktypelisthomeotypic<-list(NA,NA,NA,NA,NA,NA)
datares<-list(NA,NA,NA,NA,NA)
homeotypiclist<-list(NA,NA,NA,NA,NA,NA)

checkrest(linktypelistrestres)

linktypelistrestres<-datares  

NAremoveres(linktypelistrestres,linktypelisthomeotypic)

#filter linktypelistoneres filtering 
linktypelisthomeotypicone<-list(NA,NA,NA,NA,NA,NA)
dataone<-list(NA,NA,NA,NA,NA)
homeotypiclistone<-list(NA,NA,NA,NA,NA,NA)

checkones(linktypelistoneres)

linktypelistoneres<-dataone

NAremoveresone(linktypelistoneres,linktypelisthomeotypicone)


#filter the homeotypic peptides from original data
#can be used for either rest(homeotypic) and one(to find S and Y links not filtered by meroxsort)
findhomeo(dataset,linktypelisthomeotypic)
for(i in 1:length(homeotypiclist)){
  homeotypiclist[[i]]<-filter(homeotypiclist[[i]], Peptide2 != 1)
  linktypelisthomeotypic[[i]]<-homeotypiclist[[i]]$conresid
}

#Run each time to remove false homeotypic (S and Y)
findhomeoone(dataset,linktypelisthomeotypicone)
for(i in 1:length(homeotypiclistone[[1:2]])){
  homeotypiclistone[[i]]<-filter(homeotypiclistone[[i]], Peptide2 == 1)
  linktypelisthomeotypicone[[i]]<-homeotypiclistone[[i]]$conresid
}

#combine remaining type one and type rest
oneandrest<-list()
oneandrest<-mapply(c,linktypelistoneres, linktypelistrestres, SIMPLIFY=FALSE)
homeo<-mapply(c,linktypelisthomeotypic,linktypelisthomeotypicone)

#barplot of crosslink types
p<-list(NA,NA,NA)
crosslinktypeplot(p,linktypesres2)


#linktypesres2<-rbind(linktypesres[3,],linktypesres[6,])
#look at different types
#C4

#combined with one(homeotypic removed)
TY<-list(
  "Tysabri_light" = oneandrest[[1]],
  "Tysabri_heavy" = oneandrest[[2]],
  "tysabri_heavy+light" = oneandrest[[3]]
)
O<-list("Ocrelizumab_light" = oneandrest[[4]],
         "Ocrelizumab_heavy" = oneandrest[[5]],
         "Ocrelizumab_light+heavy" = oneandrest[[6]])

#Ven diagram

forven<-list(TY,O)


crosslinkven4(forven)
crosslinkven(forven)

savedatares(p,linktypesres,Crosslinkclasses,VENS,1)

#EXTRACT INTERSECTS AND OUTERSECTS

TOPP<-Crosslinkclasses[[1]]$..values..[2]
BOT<-Crosslinkclasses[[1]]$..values..[3]
MID<-Crosslinkclasses[[1]]$..values..[1]

####### FOR OVERNIGHT TYSABRI ANALYSIS ######################### ####
RTC110<-list(RTC_OVERNIGHT[[1]],RTC_OVERNIGHT[[3]],RTC_OVERNIGHT[[4]])
RTC1250<-list(RTC_OVERNIGHT[[5]],RTC_OVERNIGHT[[7]],RTC_OVERNIGHT[[9]],RTC_OVERNIGHT[[11]])
RTC1250<-rbind(RTC_OVERNIGHT[[5]],RTC_OVERNIGHT[[7]],RTC_OVERNIGHT[[9]],RTC_OVERNIGHT[[11]])
RTC1250T<-list(RTC_OVERNIGHT[[6]],RTC_OVERNIGHT[[8]],RTC_OVERNIGHT[[10]],RTC_OVERNIGHT[[12]])
RTC1250T<-rbind(RTC_OVERNIGHT[[6]],RTC_OVERNIGHT[[8]],RTC_OVERNIGHT[[10]],RTC_OVERNIGHT[[12]])
RTC<-list(RTC1250,RTC1250T)
#concatenate peptide 1 and peptide 2 so we known where links are between
RTC<-list(RTC_OVERNIGHT[[6]],RTC_OVERNIGHT[[8]],RTC_OVERNIGHT[[10]],RTC_OVERNIGHT[[12]],RTC_OVERNIGHT[[5]],RTC_OVERNIGHT[[7]],RTC_OVERNIGHT[[9]],RTC_OVERNIGHT[[11]])
conpepres2(RTC)
dataset<-listdata

#split depending on type of cross-link and count total number
linklist<-c("Bottom","Top")
linklist<-c("pH 4 T","pH 5 T","pH 6 T","pH 7 T","pH 4 B","pH 5 B","pH 6 B","pH 7 B")
linktypelistzero<-list(NA,NA)
linktypelistone<-list(NA,NA)
linktypelistrest<-list(NA,NA)
linktypelistzerores<-list(NA,NA)
linktypelistoneres<-list(NA,NA)
linktypelistrestres<-list(NA,NA)


crosslinkclassifyresid2(dataset = dataset,
                       linklist = linklist,
                       linktypelistzero,
                       linktypelistone,
                       linktypelistrest)


#filter linktypelistrestres filtering 
linktypelisthomeotypic<-list(NA,NA,NA,NA,NA,NA)
datares<-list(NA,NA,NA,NA,NA)
homeotypiclist<-list(NA,NA,NA,NA,NA,NA)

checkrest(linktypelistrestres[1:3])

linktypelistrestres<-datares  

NAremoveres(linktypelistrestres,linktypelisthomeotypic)

#filter linktypelistoneres filtering 
linktypelisthomeotypicone<-list(NA,NA,NA,NA,NA,NA)
dataone<-list(NA,NA,NA,NA,NA)
homeotypiclistone<-list(NA,NA,NA,NA,NA,NA)

checkones(linktypelistoneres[1:8])

linktypelistoneres<-dataone

NAremoveresone(linktypelistoneres,linktypelisthomeotypicone)


#filter the homeotypic peptides from original data
#can be used for either rest(homeotypic) and one(to find S and Y links not filtered by meroxsort)
findhomeo(dataset,linktypelisthomeotypic)
for(i in 1:length(homeotypiclist)){
  homeotypiclist[[i]]<-filter(homeotypiclist[[i]], Peptide2 != 1)
  linktypelisthomeotypic[[i]]<-homeotypiclist[[i]]$conresid
}

#Run each time to remove false homeotypic (S and Y)
findhomeoone(dataset,linktypelisthomeotypicone)
for(i in 1:length(homeotypiclistone)){
  homeotypiclistone[[i]]<-filter(homeotypiclistone[[i]], Peptide2 == 1)
  linktypelisthomeotypicone[[i]]<-homeotypiclistone[[i]]$conresid
}

#combine remaining type one and type rest
oneandrest<-list()
oneandrest<-mapply(c,linktypelistoneres[1:2], linktypelistrestres[1:2], SIMPLIFY=FALSE)
for(i in 1:length(oneandrest)){
  oneandrest[[i]] = setmincount(na.omit(oneandrest[[i]]),2)
  }
homeo<-mapply(c,linktypelisthomeotypic,linktypelisthomeotypicone)


#barplot of crosslink types
p<-list(NA,NA)
crosslinktypeplot(p,linktypesres)


#look at different types
#C4

#combined with one(homeotypic removed)
O<- list(
  "Bottom" = oneandrest[[1]],
  "Top" = oneandrest[[2]]
  #"pH 6 T"= oneandrest[[3]],
  #"pH 7 T" =oneandrest[[4]]
  #"pH 7 T"=oneandrest[[4]]
)

#Ven diagram

forven<-list(O)


crosslinkven4(forven)
crosslinkven(forven)

savedatares(p,linktypesres,Crosslinkclasses,VENS,1)

#EXTRACT INTERSECTS AND OUTERSECTS

TOPP<-Crosslinkclasses[[1]]$..values..[2]
BOT<-Crosslinkclasses[[1]]$..values..[3]
MID<-Crosslinkclasses[[1]]$..values..[1]


####### FOR FILTERED PH TYSABRI ################################# #################
RTC110<-list(RTC_FILTER[[1]],RTC_FILTER[[2]],RTC_FILTER[[3]])
RTC1250<-list(RTC_FILTER[[5]],RTC_FILTER[[7]],RTC_FILTER[[9]],RTC_FILTER[[10]])
RTC1250T<-list(RTC_FILTER[[4]],RTC_FILTER[[6]],RTC_FILTER[[8]],RTC_FILTER[[11]])
RTC1250<-rbind(RTC1250[[1]],RTC1250[[2]],RTC1250[[3]],RTC1250[[4]])
RTC1250T<-rbind(RTC1250T[[1]],RTC1250T[[2]],RTC1250T[[3]],RTC1250T[[4]], SIMPLIFY=FALSE)
BT<-list(RTC1250,RTC1250T)
#concatenate peptide 1 and peptide 2 so we known where links are between
RTC_TEMP<-list(RTC_FILTER[[5]],RTC_FILTER[[7]],RTC_FILTER[[9]],RTC_FILTER[[10]],RTC_FILTER[[4]],RTC_FILTER[[6]],RTC_FILTER[[8]],RTC_FILTER[[11]])
conpepres2(BT)
dataset<-listdata

#split depending on type of cross-link and count total number
linklist<-c("B","T")#"pH 4 T","pH 5 T","pH 6 T","pH 7 T")
linktypelistzero<-list(NA,NA,NA,NA)
linktypelistone<-list(NA,NA,NA,NA)
linktypelistrest<-list(NA,NA,NA,NA)
linktypelistzerores<-list(NA,NA,NA,NA)
linktypelistoneres<-list(NA,NA,NA,NA)
linktypelistrestres<-list(NA,NA,NA,NA)


crosslinkclassifyresid2(dataset = dataset,
                        linklist = linklist,
                        linktypelistzero,
                        linktypelistone,
                        linktypelistrest)


#filter linktypelistrestres filtering 
linktypelisthomeotypic<-list(NA,NA,NA,NA,NA,NA)
datares<-list(NA,NA,NA,NA,NA)
homeotypiclist<-list(NA,NA,NA,NA,NA,NA)

checkrest(linktypelistrestres)

linktypelistrestres<-datares  

NAremoveres(linktypelistrestres,linktypelisthomeotypic)

#filter linktypelistoneres filtering 
linktypelisthomeotypicone<-list(NA,NA,NA,NA,NA,NA)
dataone<-list(NA,NA,NA,NA,NA)
homeotypiclistone<-list(NA,NA,NA,NA,NA,NA)

checkones(linktypelistoneres)

linktypelistoneres<-dataone

NAremoveresone(linktypelistoneres,linktypelisthomeotypicone)


#filter the homeotypic peptides from original data
#can be used for either rest(homeotypic) and one(to find S and Y links not filtered by meroxsort)
findhomeo(dataset,linktypelisthomeotypic)
for(i in 1:length(homeotypiclist)){
  homeotypiclist[[i]]<-filter(homeotypiclist[[i]], Peptide2 != 1)
  linktypelisthomeotypic[[i]]<-homeotypiclist[[i]]$conresid
}

#Run each time to remove false homeotypic (S and Y)
findhomeoone(dataset,linktypelisthomeotypicone)
for(i in 1:length(homeotypiclistone)){
  homeotypiclistone[[i]]<-filter(homeotypiclistone[[i]], Peptide2 == 1)
  linktypelisthomeotypicone[[i]]<-homeotypiclistone[[i]]$conresid
}

#combine remaining type one and type rest
oneandrest<-list()
oneandrest<-mapply(c,linktypelistoneres, linktypelistrestres, SIMPLIFY=FALSE)
homeo<-mapply(c,linktypelisthomeotypic,linktypelisthomeotypicone)


#barplot of crosslink types
p<-list(NA,NA,NA)
crosslinktypeplot(p,linktypesres)


#look at different types
#C4

#combined with one(homeotypic removed)
O<- list(
  "B" = oneandrest[[1]],
  "T" = oneandrest[[2]]
  
  
)

#Ven diagram

forven<-list(O)


crosslinkven4(forven)
crosslinkven(forven)

savedatares(p,linktypesres,Crosslinkclasses,VENS,1)


#EXTRACT INTERSECTS AND OUTERSECTS

TOPP<-Crosslinkclasses[[1]]$..values..[2]
BOT<-Crosslinkclasses[[1]]$..values..[3]
MID<-Crosslinkclasses[[1]]$..values..[1]




####### FOR INFLIXIMAB TRIPLICATE 1:250 C4 ##################### ##########################
setwd("C:\\Users\\victo\\Desktop\\speciale\\R\\MEROXSORT\\MeroxAdjust")
INF30<-read_xlsx("final_INFLIX130_dataset.xlsx")
INF250<-read_xlsx("final_INFLIX1250_dataset.xlsx")
dataset<-list(INF30,INF250)
linklist<-c("Infliximab 1:30","Infliximab 1:250")


for(i in 1:dim(INF250)[1]){
  if(INF250$Protein.1[j] == ">light "){
    INF250$Protein.1[j] = ">light"
  }
}
#concatenate peptide 1 and peptide 2 so we known where links are between
conpepres2(RTCINF30)
RTCINF<-listdata


#conresid 
dataset<-list(RTCINF[[1]],RTCINF[[2]],RTCINF[[3]],RTCINF250[[1]],RTCINF250[[2]],RTCINF250[[3]])
linklist<-c("INF_30_A","INF_30_B","INF_30_C","INF_250_A","INF_250_B","INF_250_C")

linktypelistzero<-list(NA,NA)
linktypelistone<-list(NA,NA,NA)
linktypelistrest<-list(NA,NA,NA)
linktypelistzerores<-list(NA,NA,NA)
linktypelistoneres<-list(NA,NA,NA)
linktypelistrestres<-list(NA,NA,NA)


crosslinkclassifyresid2(dataset = dataset,
                       linklist = linklist,
                       linktypelistzero,
                       linktypelistone,
                       linktypelistrest)


#filter linktypelistrestres filtering 
linktypelisthomeotypic<-list(NA,NA,NA,NA,NA,NA)
datares<-list(NA,NA,NA,NA,NA)
homeotypiclist<-list(NA,NA,NA,NA,NA,NA)

checkrest(linktypelistrestres)

linktypelistrestres<-datares  

NAremoveres(linktypelistrestres,linktypelisthomeotypic)

#filter linktypelistoneres filtering 
linktypelisthomeotypicone<-list(NA,NA,NA,NA,NA,NA)
dataone<-list(NA,NA,NA,NA,NA)
homeotypiclistone<-list(NA,NA,NA,NA,NA,NA)

checkones(linktypelistoneres)

linktypelistoneres<-dataone

NAremoveresone(linktypelistoneres,linktypelisthomeotypicone)


#filter the homeotypic peptides from original data
#can be used for either rest(homeotypic) and one(to find S and Y links not filtered by meroxsort)
findhomeo(dataset,linktypelisthomeotypic)
for(i in 1:length(homeotypiclist)){
  homeotypiclist[[i]]<-filter(homeotypiclist[[i]], Peptide2 != 1)
  linktypelisthomeotypic[[i]]<-homeotypiclist[[i]]$conresid
}

#Run each time to remove false homeotypic (S and Y)
findhomeoone(dataset,linktypelisthomeotypicone)
for(i in 1:length(homeotypiclistone)){
  homeotypiclistone[[i]]<-filter(homeotypiclistone[[i]], Peptide2 == 1)
  linktypelisthomeotypicone[[i]]<-homeotypiclistone[[i]]$conresid
}

#combine remaining type one and type rest
oneandrest<-list()
oneandrest<-mapply(c,linktypelistoneres[1:2], linktypelistrestres[1:2], SIMPLIFY=FALSE)
homeo<-mapply(c,linktypelisthomeotypic,linktypelisthomeotypicone)


#barplot of crosslink types
p<-list(NA,NA,NA)
crosslinktypeplot(p,linktypesres)


#look at different types
#C4

#combined with one(homeotypic removed)
#Used to filter data using T5
O2<- list(
  A = setmincount(na.omit(oneandrest[[1]]),2), 
  B = setmincount(na.omit(oneandrest[[2]]),2), 
  C = setmincount(na.omit(oneandrest[[3]]),2)
)

TY2 <- list(
  A = setmincount(na.omit(oneandrest[[4]]),2), 
  B = setmincount(na.omit(oneandrest[[5]]),2), 
  C = setmincount(na.omit(oneandrest[[6]]),2)
)

#Ven diagram

forven<-list(TY2,O2)


crosslinkven(forven)

#uniquelists<-list()
barplotsres<-list()
linktypeslistres<-list()
crosslinklistres<-list()
venslistres<-list()
savedatares(p,linktypesres,Crosslinkclasses,VENS,1)


#get unqiue list for comparison

INF<-ldply(dataset[1:6],data.frame)
INF<-INF[!duplicated(INF$conresid),]
writexl::write_xlsx(INFUNIQUE,path = "uniqueinf.xlsx")





####### FOR OCRELIZUMAB AND TYSABRI C4 1:250 #################### #########################
#concatenate peptide 1 and peptide 2 so we known where links are between
conpepres2(RTC4)
RTC4<-listdata

#conresid 
#after filtration import data
#ocrelizumab<-read_xlsx("final_ocrelizumab_dataset.xlsx")
#tysabri<-read_xlsx("final_tysabri_dataset.xlsx")
setwd("C:\\Users\\victo\\Desktop\\speciale\\R\\MEROXSORT\\MeroxAdjust")
ocrelizumab<-read_xlsx("final_ocrelizumab_dataset_type0(+)2.xlsx")
tysabri<-read_xlsx("final_tysabri_dataset_type0(+)2.xlsx")
dataset<-list(ocrelizumab,tysabri)
linklist<-c("Ocrelizumab","Natalizumab")

#pca
#library(factoextra)
#library("RColorBrewer")
#library(pca3d)
#library(plyr)
dataset<-list(RTC4[[1]],RTC4[[2]],RTC4[[3]],RTC4[[4]],RTC4[[5]],RTC4[[6]])

for(i in 1:length(dataset)){
dataset[[i]]<-cbind(dataset[[i]][,1:6],dataset[[i]][,9:10],dataset[[i]][,13:14],dataset[[i]][,17],dataset[[i]][,19],dataset[[i]][,28:30],dataset[[i]][,33],linklist[i])
}

datapca2<-ldply(dataset[1:6], data.frame)
for(i in 1:16){
  datapca2[,i]<-as.numeric(datapca2[,i])
}
datapca2$linklist.i.<-factor(datapca2$linklist.i.)

res.pca<-prcomp(datapca2[,-17],scale =T)

fviz_eig(res.pca)
fviz_pca_biplot(res.pca)

fviz_pca_ind(res.pca, 
             geom.ind = "point",
             col.ind = datapca2$linklist.i.,
             palette = brewer.pal(n=6,name="Dark2"),
             addEllipses = TRUE,
             legend.title = "Groups")

#pca3d(res.pca,group = datapca2$linklist.i., show.ellipses = T,
 #    ellipse.ci = 0.75, show.plane =F)
linklist<-c("OA","OB","OC","NA","NB","NC")

linktypelistzero<-list(NA,NA,NA,NA,NA)
linktypelistone<-list(NA,NA,NA,NA,NA,NA)
linktypelistrest<-list(NA,NA,NA,NA,NA,NA)

linktypelistzerores<-list(NA,NA)
linktypelistoneres<-list(NA,NA)
linktypelistrestres<-list(NA,NA)


crosslinkclassifyresid2(dataset = dataset,
                       linklist = linklist,
                       linktypelistzero,
                       linktypelistone,
                       linktypelistrest)

#filter linktypelistrestres filtering 
linktypelisthomeotypic<-list(NA,NA)
datares<-list(NA,NA)
homeotypiclist<-list(NA,NA)

checkrest(linktypelistrestres[1:2])

linktypelistrestres<-datares  

NAremoveres(linktypelistrestres,linktypelisthomeotypic)

#filter linktypelistoneres filtering 
linktypelisthomeotypicone<-list(NA,NA,NA,NA,NA,NA)
dataone<-list(NA,NA,NA,NA,NA)
homeotypiclistone<-list(NA,NA,NA,NA,NA,NA)

checkones(linktypelistoneres[1:2])

linktypelistoneres<-dataone

NAremoveresone(linktypelistoneres,linktypelisthomeotypicone)


#filter the homeotypic peptides from original data
#can be used for either rest(homeotypic) and one(to find S and Y links not filtered by meroxsort)
findhomeo(dataset,linktypelisthomeotypic)
for(i in 1:length(homeotypiclist)){
homeotypiclist[[i]]<-filter(homeotypiclist[[i]], Peptide2 != 1)
linktypelisthomeotypic[[i]]<-homeotypiclist[[i]]$conresid
}

#Run each time to remove false homeotypic (S and Y)
findhomeoone(dataset,linktypelisthomeotypicone)
for(i in 1:length(homeotypiclistone)){
  homeotypiclistone[[i]]<-filter(homeotypiclistone[[i]], Peptide2 == 1)
  linktypelisthomeotypicone[[i]]<-homeotypiclistone[[i]]$conresid
}

#combine remaining type one and type rest

oneandrest<-list()
oneandrest<-mapply(c,linktypelistoneres, linktypelistrestres,linktypelistzerores ,SIMPLIFY=FALSE)
oneandrest<-mapply(c,linktypelistoneres[1:2], linktypelistrestres[1:2],SIMPLIFY=FALSE)



#DOUBLE CHECK TO SEE IF HOMEOTYPIC STILL ARE PRESENT
for(i in 1:length(oneandrest[[1]])){
  if(oneandrest[[1]][i]=="106|106|>light|>light" ){
    print("YES")
  }
}

homeo<-mapply(c,linktypelisthomeotypic,linktypelisthomeotypicone)
#barplot of crosslink types
p<-list(NA)

crosslinktypeplot(p,linktypesres)


#look at different types
#C4

#combined with one(homeotypic removed)
O<- list(
  A = na.omit(oneandrest[[1]]), 
  B = na.omit(oneandrest[[2]]), 
  C = na.omit(oneandrest[[3]])
)

TY <- list(
  A = na.omit(oneandrest[[4]]), 
  B = na.omit(oneandrest[[5]]), 
  C = na.omit(oneandrest[[6]])
)

#Used to filter data using T5
O2<- list(
  A = setmincount(na.omit(oneandrest[[1]]),2), 
  B = setmincount(na.omit(oneandrest[[2]]),2), 
  C = setmincount(na.omit(oneandrest[[3]]),2)
)

TY2 <- list(
  A = setmincount(na.omit(oneandrest[[4]]),2), 
  B = setmincount(na.omit(oneandrest[[5]]),2), 
  C = setmincount(na.omit(oneandrest[[6]]),2)
)


O3<-ldply(O,data.frame)
TY3<-ldply(TY,data.frame)
O3include<-setmincount_merged(O3,3)
TY3include<-setmincount_merged(TY3,3)

#5 = 82.7 68.7
O3<- list(
  A = na.omit(oneandrest[[1]][oneandrest[[1]] %in% O3include$AA]), 
  B = na.omit(oneandrest[[2]][oneandrest[[2]] %in% O3include$AA]), 
  C = na.omit(oneandrest[[3]][oneandrest[[3]] %in% O3include$AA])
)

TY3 <- list(
  A = na.omit(oneandrest[[4]][oneandrest[[4]] %in% TY3include$AA]), 
  B = na.omit(oneandrest[[5]][oneandrest[[5]] %in% TY3include$AA]), 
  C = na.omit(oneandrest[[6]][oneandrest[[6]] %in% TY3include$AA])
)


#combined with one(homeotypic removed)
Oh<- list(
  A = na.omit(homeo[[1]]), 
  B = na.omit(homeo[[2]]), 
  C = na.omit(homeo[[3]])
)

TYh <- list(
  A = na.omit(homeo[[4]]), 
  B = na.omit(homeo[[5]]), 
  C = na.omit(homeo[[6]])
)



setmincount<-function(data,min_count){
dat<-as.data.frame(data)
colnames(dat)<-c("AA")
test<-as.data.frame(table(dat))

test2<-dplyr::filter(test, Freq > min_count)

test3<-dat[dat$AA %in% test2$dat,]
return(test3)
}

setmincount_merged<-function(data,min_count){
  dat<-as.data.frame(data)
  colnames(dat)<-c("ID","AA")
  test<-as.data.frame(table(dat$AA))
  
  test2<-dplyr::filter(test, Freq > min_count)
  
  test3<-unique(dat[dat$AA %in% test2$Var1,])
  return(test3)
}

#Ven diagram

forven<-list(TY2,O2)


crosslinkven(forven)


#uniquelists<-list()
barplotsres<-list()
linktypeslistres<-list()
crosslinklistres<-list()
venslistres<-list()
savedatares(p,linktypesres,Crosslinkclasses,VENS,1)

joinedtysabri<-c(secs[[1]]$`2`,secs[[2]]$`3`,secs[[3]]$`5`,secs[[4]]$`1`)
binded<-ldply(dataset[4:6],data.frame)
joinedtysabri<-binded[binded$conresid %in% joinedtysabri,]

joinedocrelizumab<-c(secs[[5]]$`2`,secs[[6]]$`3`,secs[[7]]$`5`,secs[[8]]$`1`)
binded<-ldply(RTC2[1:3],data.frame)
joinedocrelizumab<-binded[binded$conresid %in% joinedocrelizumab,]

writexl::write_xlsx(joinedocrelizumab,path = "final_INFLIX130_datasettype0.xlsx")
writexl::write_xlsx(joinedtysabri,path = "final_INFLIX1250_datasettype0.xlsx")


#Tysabri
ABC<-list(NA,NA)
trip<-list(NA,NA)
secs<-list(NA,NA)
#extract 1/3

  ABC[[1]]<-Crosslinkclasses[[1]]$..values..[7]
  ABC[[2]]<-Crosslinkclasses[[1]]$..values..[6]
  ABC[[3]]<-Crosslinkclasses[[1]]$..values..[4]

#extract 3/3
trip[[1]]<-Crosslinkclasses[[1]]$..values..[1]

#extract 2/3 + 3/3

  secs[[1]]<-Crosslinkclasses[[1]]$..values..[2]
  secs[[2]]<-Crosslinkclasses[[1]]$..values..[3]
  secs[[3]]<-Crosslinkclasses[[1]]$..values..[5]
  secs[[4]]<-Crosslinkclasses[[1]]$..values..[1]


#Ocrelizumab
#extract 1/3

  ABC[[4]]<-Crosslinkclasses[[2]]$..values..[7]
  ABC[[5]]<-Crosslinkclasses[[2]]$..values..[6]
  ABC[[6]]<-Crosslinkclasses[[2]]$..values..[4]

#extract 3/3
trip[[2]]<-Crosslinkclasses[[2]]$..values..[1]

#extract 2/3 + 3/3

  secs[[5]]<-Crosslinkclasses[[2]]$..values..[2]
  secs[[6]]<-Crosslinkclasses[[2]]$..values..[3]
  secs[[7]]<-Crosslinkclasses[[2]]$..values..[5]
  secs[[8]]<-Crosslinkclasses[[2]]$..values..[1]


test<-list()
test2<-list()
#for(i in 1:3){
binded<-ldply(dataset[1:3],data.frame)
test<-binded[binded$conresid %in% trip[[2]][["1"]],]
test<-checkrev(test)
test2<-table(test$conresid)  

ocrelizmab_33_resid_count<-ldply(test2,data.frame)
ocrelizmab_33_merox<-test[!duplicated(test$conresid),]

#get the noise filtered data for ocrelizumab 2/3+3/3
arrange(ocrelizmab_33_resid_count, by = X..i..)
  binded<-ldply(listdata[1:3],data.frame)
  test<-binded[binded$conresid %in% c(secs[[5]][[1]]),]
  test.1<-binded[binded$conresid %in% secs[[6]][[1]],]
  test.2<-binded[binded$conresid %in% secs[[7]][[1]],]
  test.3<-binded[binded$conresid %in% secs[[8]][[1]],]
  test<-checkrev(test)
  test.1<-checkrev(test.1)
  test.2<-checkrev(test.2)
  test.3<-checkrev(test.3)
  test2<-table(test$conresid)
  test2.1<-table(test.1$conresid)
  test2.2<-table(test.2$conresid)
  test2.3<-table(test.3$conresid)
  
  ocrelizmab_23_resid_count<-ldply(c(test2,test2.1,test2.2,test2.3),data.frame)
  binded<-rbind(test,test.1,test.2,test.3)
  ocrelizmab_23_merox<-binded[!duplicated(binded$conresid),]
  
  
  binded<-ldply(RTC4[1:3],data.frame)
  test<-binded[binded$conresid %in% c(ABC[[4]][[1]]),]
  test.1<-binded[binded$conresid %in% ABC[[5]][[1]],]
  test.2<-binded[binded$conresid %in% ABC[[6]][[1]],]
  test<-checkrev(test)
  test.1<-checkrev(test.1)
  test.2<-checkrev(test.2)
  test2<-table(test$conresid)
  test2.1<-table(test.1$conresid)
  test2.2<-table(test.2$conresid)
  
  ocrelizmab_13_resid_count<-ldply(c(test2,test2.1,test2.2),data.frame)
  binded<-rbind(test,test.1,test.2)
  ocrelizmab_13_merox<-binded[!duplicated(binded$conresid),]

  
  binded<-ldply(RTC4[4:6],data.frame)
  test<-binded[binded$conresid %in% trip[[1]][["1"]],]
  test<-checkrev(test)
  test2<-table(test$conresid)
  
  tysabri_33_resid_count<-ldply(test2,data.frame)
  tysabri_33_merox<-test[!duplicated(test$conresid),]
  

  binded<-ldply(listdata[4:6],data.frame)
  test<-binded[binded$conresid %in% c(secs[[1]][[1]]),]
  test.1<-binded[binded$conresid %in% secs[[2]][[1]],]
  test.2<-binded[binded$conresid %in% secs[[3]][[1]],]
  test.3<-binded[binded$conresid %in% secs[[4]][[1]],]
  test<-checkrev(test)
  test.1<-checkrev(test.1)
  test.2<-checkrev(test.2)
  test.3<-checkrev(test.3)
  test2<-table(test$conresid)
  test2.1<-table(test.1$conresid)
  test2.2<-table(test.2$conresid)
  test2.3<-table(test.3$conresid)
  
  tysabri_23_resid_count<-ldply(c(test2,test2.1,test2.2,test2.3),data.frame)
  binded<-rbind(test,test.1,test.2,test.3)
  tysabri_23_merox<-binded[!duplicated(binded$conresid),]
  
  binded<-ldply(RTC4[4:6],data.frame)
  test<-binded[binded$conresid %in% c(ABC[[1]][[1]]),]
  test.1<-binded[binded$conresid %in% ABC[[2]][[1]],]
  test.2<-binded[binded$conresid %in% ABC[[3]][[1]],]
  test<-checkrev(test)
  test.1<-checkrev(test.1)
  test.2<-checkrev(test.2)
  test2<-table(test$conresid)
  test2.1<-table(test.1$conresid)
  test2.2<-table(test.2$conresid)
  
  tysabri_13_resid_count<-ldply(c(test2,test2.1,test2.2),data.frame)
  binded<-rbind(test,test.1,test.2)
  tysabri_13_merox<-binded[!duplicated(binded$conresid),]

#write.csv2(binded,file = "ocrleizumab13")
setwd("C:\\Users\\victo\\Desktop\\speciale\\artikler\\IgG\\Data")
writexl::write_xlsx(list("ocrelizmab_33_resid_count" = ocrelizmab_33_resid_count ,"ocrelizmab_23_resid_count"=ocrelizmab_23_resid_count,"ocrelizmab_13_resid_count"=ocrelizmab_13_resid_count),path = "Ocrelizumab_resid_count_false_homeo4.xlsx")
writexl::write_xlsx(list("ocrelizmab_33_merox" = ocrelizmab_33_merox ,"ocrelizmab_23_merox"=ocrelizmab_23_merox,"ocrelizmab_13_merox"=ocrelizmab_13_merox),path = "Ocrelizumab_merox_false_homeo(-)4.xlsx")

writexl::write_xlsx(list("tysabri_33_resid_count" = tysabri_33_resid_count ,"tysabri_23_resid_count"=tysabri_23_resid_count,"tysabri_13_resid_count"= tysabri_13_resid_count),path = "tysabri_resid_count_false_homeo(-)4.xlsx")
writexl::write_xlsx(list("tysabri_33_merox" = tysabri_33_merox ,"tysabri_23_merox"=tysabri_23_merox,"tysabri_13_merox"= tysabri_13_merox),path = "tysabri_merox_false_homeo(-)4.xlsx")

writexl::write_xlsx(binded,path = "final_tysabri_dataset_type0(+).xlsx")


####### For OCRELIZUMAB AND TYSABRI C2 1:250 #################### ####
#concatenate peptide 1 and peptide 2 so we known where links are between
conpepres2(RTC2)
RTC2<-listdata

dataset<-list(joinedocrelizumab,joinedtysabri)
linklist<-c("Ocrelizumab C2","Tysabri C2")

#split depending on type of cross-link and count total number
dataset<-list(RTC2[[1]],RTC2[[2]],RTC2[[3]],RTC2[[4]],RTC2[[5]],RTC2[[6]])
linklist<-c("OAC2","OBC2","OCC2","TYAC2","TYBC2","TYCC2")
linktypelistzero<-list(NA,NA,NA,NA,NA,NA)
linktypelistone<-list(NA,NA,NA,NA,NA,NA)
linktypelistrest<-list(NA,NA,NA,NA,NA,NA)
linktypelistzerores<-list(NA,NA,NA,NA,NA,NA)
linktypelistoneres<-list(NA,NA,NA,NA,NA,NA)
linktypelistrestres<-list(NA,NA,NA,NA,NA,NA)

crosslinkclassifyresid2(dataset = dataset,
                       linklist = linklist,
                       linktypelistzero,
                       linktypelistone,
                       linktypelistrest)

#filter linktypelistoneres filtering 
linktypelisthomeotypicone<-list(NA,NA,NA,NA,NA,NA)
dataone<-list(NA,NA,NA,NA,NA)
homeotypiclistone<-list(NA,NA,NA,NA,NA,NA)

checkones(linktypelistoneres)

linktypelistoneres<-dataone

NAremoveresone(linktypelistoneres,linktypelisthomeotypicone)

#combine remaining type one and type rest
oneandrest<-list()
oneandrest<-mapply(c,linktypelistoneres, linktypelistrestres, SIMPLIFY=FALSE)

#filter the homeotypic peptides from original data
findhomeo(RTC2,linktypelisthomeotypic)

#barplot of crosslink types
p<-list(NA,NA,NA)
crosslinktypeplot(p,linktypesres)


#look at different types
#C4

#combined with one(homeotypic removed)
O2<- list(
  A = setmincount(na.omit(oneandrest[[1]]),2), 
  B = setmincount(na.omit(oneandrest[[2]]),2), 
  C = setmincount(na.omit(oneandrest[[3]]),2)
)

TY2 <- list(
  A = setmincount(na.omit(oneandrest[[4]]),2), 
  B = setmincount(na.omit(oneandrest[[5]]),2), 
  C = setmincount(na.omit(oneandrest[[6]]),2)
)

#Ven diagram

forven<-list(TY2,O2)


crosslinkven(forven)


savedatares(p,linktypesres,Crosslinkclasses,VENS,2)


####### For OCRELIZUMAB AND TYSABRI C4 1:250 DIMER ################ ####
#concatenate peptide 1 and peptide 2 so we known where links are between
conpepres2(RTC4D)
RTC4D<-listdata

#conresid 
dataset<-list(RTC4D[[1]],RTC4D[[2]],RTC4D[[3]],RTC4D[[4]],RTC4D[[5]],RTC4D[[6]])
linklist<-c("OAD","OBD","OCD","TYAD","TYBD","TYCD")

linktypelistzero<-list(NA,NA,NA,NA,NA,NA)
linktypelistone<-list(NA,NA,NA,NA,NA,NA)
linktypelistrest<-list(NA,NA,NA,NA,NA,NA)

linktypelistzerores<-list(NA,NA,NA,NA,NA,NA)
linktypelistoneres<-list(NA,NA,NA,NA,NA,NA)
linktypelistrestres<-list(NA,NA,NA,NA,NA,NA)


crosslinkclassifyresid2(dataset = dataset,
                        linklist = linklist,
                        linktypelistzero,
                        linktypelistone,
                        linktypelistrest)



#filter linktypelistrestres filtering 
linktypelisthomeotypic<-list(NA,NA,NA,NA,NA,NA)
datares<-list(NA,NA,NA,NA,NA)
homeotypiclist<-list(NA,NA,NA,NA,NA,NA)

checkrest(linktypelistrestres)

linktypelistrestres<-datares  

NAremoveres(linktypelistrestres,linktypelisthomeotypic)

#filter linktypelistoneres filtering 
linktypelisthomeotypicone<-list(NA,NA,NA,NA,NA,NA)
dataone<-list(NA,NA,NA,NA,NA)
homeotypiclistone<-list(NA,NA,NA,NA,NA,NA)

checkones(linktypelistoneres)

linktypelistoneres<-dataone

NAremoveresone(linktypelistoneres,linktypelisthomeotypicone)


#filter the homeotypic peptides from original data
#can be used for either rest(homeotypic) and one(to find S and Y links not filtered by meroxsort)
findhomeo(dataset,linktypelisthomeotypic)
for(i in 1:length(homeotypiclist)){
  homeotypiclist[[i]]<-filter(homeotypiclist[[i]], Peptide2 != 1)
  linktypelisthomeotypic[[i]]<-homeotypiclist[[i]]$conresid
}

#Run each time to remove false homeotypic (S and Y)
findhomeoone(dataset,linktypelisthomeotypicone)
for(i in 1:length(homeotypiclistone)){
  homeotypiclistone[[i]]<-filter(homeotypiclistone[[i]], Peptide2 == 1)
  linktypelisthomeotypicone[[i]]<-homeotypiclistone[[i]]$conresid
}

#combine remaining type one and type rest
oneandrest<-list()
oneandrest<-mapply(c,linktypelistoneres, linktypelistrestres, SIMPLIFY=FALSE)
homeo<-mapply(c,linktypelisthomeotypic,linktypelisthomeotypicone)
#barplot of crosslink types
p<-list(NA)
crosslinktypeplot(p,linktypesres)


#look at different types
#C4

#combined with one(homeotypic removed)
O<- list(
  A = na.omit(oneandrest[[1]]), 
  B = na.omit(oneandrest[[2]]), 
  C = na.omit(oneandrest[[3]])
)

TY <- list(
  A = na.omit(oneandrest[[4]]), 
  B = na.omit(oneandrest[[5]]), 
  C = na.omit(oneandrest[[6]])
)
#Used to filter data using T5
O2<- list(
  A = setmincount(na.omit(oneandrest[[1]]),2), 
  B = setmincount(na.omit(oneandrest[[2]]),2), 
  C = setmincount(na.omit(oneandrest[[3]]),2)
)

TY2 <- list(
  A = setmincount(na.omit(oneandrest[[4]]),2), 
  B = setmincount(na.omit(oneandrest[[5]]),2), 
  C = setmincount(na.omit(oneandrest[[6]]),2)
)

Oh<- list(
  A = na.omit(homeo[[1]]), 
  B = na.omit(homeo[[2]]), 
  C = na.omit(homeo[[3]])
)

TYh <- list(
  A = na.omit(homeo[[4]]), 
  B = na.omit(homeo[[5]]), 
  C = na.omit(homeo[[6]])
)
#Ven diagram

forven<-list(TY2,O2)


crosslinkven(forven)

savedatares(p,linktypesres,Crosslinkclasses,VENS,1)





####### COMPARE DIMER AND MONOMER OCRELIZUMAB AND TYSABRI ########### ####
C4O<-list()
C4O[[1]]<-bind_rows(RTC4[1],RTC4[2],RTC4[3])
C4O[[2]]<-bind_rows(RTC4D[1],RTC4D[2],RTC4D[3])
C4TY<-list()
C4TY[[1]]<-bind_rows(RTC4[4],RTC4[5],RTC4[6])
C4TY[[2]]<-bind_rows(RTC4D[4],RTC4D[5],RTC4D[6])
#concatenate peptide 1 and peptide 2 so we known where links are between
conpepres2(C4O)
C4O<-listdata
conpepres2(C4TY)
C4TY<-listdata

#conresid 
dataset<-list(C4O[[1]],C4O[[2]],C4TY[[1]],C4TY[[2]])
linklist<-c("OM","OD","TYM","TYD")

linktypelistzero<-list(NA,NA,NA,NA)
linktypelistone<-list(NA,NA,NA,NA)
linktypelistrest<-list(NA,NA,NA,NA)

linktypelistzerores<-list(NA,NA,NA,NA)
linktypelistoneres<-list(NA,NA,NA,NA)
linktypelistrestres<-list(NA,NA,NA,NA)


crosslinkclassifyresid2(dataset = dataset,
                        linklist = linklist,
                        linktypelistzero,
                        linktypelistone,
                        linktypelistrest)



#filter linktypelistrestres filtering 
linktypelisthomeotypic<-list(NA,NA,NA,NA)
datares<-list(NA,NA,NA,NA)
homeotypiclist<-list(NA,NA,NA,NA)

checkrest(linktypelistrestres)

linktypelistrestres<-datares  

NAremoveres(linktypelistrestres,linktypelisthomeotypic)

#filter linktypelistoneres filtering 
linktypelisthomeotypicone<-list(NA,NA,NA,NA)
dataone<-list(NA,NA,NA,NA)
homeotypiclistone<-list(NA,NA,NA,NA)

checkones(linktypelistoneres)

linktypelistoneres<-dataone

NAremoveresone(linktypelistoneres,linktypelisthomeotypicone)


#filter the homeotypic peptides from original data
#can be used for either rest(homeotypic) and one(to find S and Y links not filtered by meroxsort)
findhomeo(dataset,linktypelisthomeotypic)
for(i in 1:length(homeotypiclist)){
  homeotypiclist[[i]]<-filter(homeotypiclist[[i]], Peptide2 != 1)
  linktypelisthomeotypic[[i]]<-homeotypiclist[[i]]$conresid
}

#Run each time to remove false homeotypic (S and Y)
findhomeoone(dataset,linktypelisthomeotypicone)
for(i in 1:length(homeotypiclistone)){
  homeotypiclistone[[i]]<-filter(homeotypiclistone[[i]], Peptide2 == 1)
  linktypelisthomeotypicone[[i]]<-homeotypiclistone[[i]]$conresid
}

#combine remaining type one and type rest
oneandrest<-list()
oneandrest<-mapply(c,linktypelistoneres, linktypelistrestres, SIMPLIFY=FALSE)
homeo<-mapply(c,linktypelisthomeotypic,linktypelisthomeotypicone)
#barplot of crosslink types
p<-list(NA)
crosslinktypeplot(p,linktypesres)


#look at different types
#C4

#combined with one(homeotypic removed)
O<- list(
  A = na.omit(oneandrest[[1]]), 
  B = na.omit(oneandrest[[2]])
)

TY <- list(
  A = na.omit(oneandrest[[3]]), 
  B = na.omit(oneandrest[[4]])
)

Oh<- list(
  A = na.omit(homeo[[1]]), 
  B = na.omit(homeo[[2]])
)

TYh <- list(
  A = na.omit(homeo[[3]]), 
  B = na.omit(homeo[[4]])
)
#Ven diagram

forven<-list(TY,O)


crosslinkven(forven)

savedatares(p,linktypesres,Crosslinkclasses,VENS,1)



####### TYSABRI C4 1:250 TEMPERATURE ################################# ####
#concatenate peptide 1 and peptide 2 so we known where links are between
conpepres2(RTC_TEMP)
RTC_TEMP<-listdata

#conresid 
dataset<-list(RTC_TEMP[[1]],RTC_TEMP[[2]],RTC_TEMP[[3]],RTC_TEMP[[4]])
linklist<-c("N37","N47","N57","N67")

linktypelistzero<-list(NA,NA,NA,NA)
linktypelistone<-list(NA,NA,NA,NA)
linktypelistrest<-list(NA,NA,NA,NA)

linktypelistzerores<-list(NA,NA,NA,NA)
linktypelistoneres<-list(NA,NA,NA,NA)
linktypelistrestres<-list(NA,NA,NA,NA)


crosslinkclassifyresid2(dataset = dataset,
                       linklist = linklist,
                       linktypelistzero,
                       linktypelistone,
                       linktypelistrest)


#filter linktypelistrestres filtering 
linktypelisthomeotypic<-list(NA,NA,NA,NA)
datares<-list(NA,NA,NA,NA)
homeotypiclist<-list(NA,NA,NA,NA)

checkrest(linktypelistrestres)

linktypelistrestres<-datares  

NAremoveres(linktypelistrestres,linktypelisthomeotypic)

#filter linktypelistoneres filtering 
linktypelisthomeotypicone<-list(NA,NA,NA,NA)
dataone<-list(NA,NA,NA,NA)
homeotypiclistone<-list(NA,NA,NA,NA)

checkones(linktypelistoneres)

linktypelistoneres<-dataone

NAremoveresone(linktypelistoneres,linktypelisthomeotypicone)


#filter the homeotypic peptides from original data
#can be used for either rest(homeotypic) and one(to find S and Y links not filtered by meroxsort)
findhomeo(dataset,linktypelisthomeotypic)
for(i in 1:length(homeotypiclist)){
  homeotypiclist[[i]]<-filter(homeotypiclist[[i]], Peptide2 != 1)
  linktypelisthomeotypic[[i]]<-homeotypiclist[[i]]$conresid
}

#Run each time to remove false homeotypic (S and Y)
findhomeoone(dataset,linktypelisthomeotypicone)
for(i in 1:length(homeotypiclistone)){
  homeotypiclistone[[i]]<-filter(homeotypiclistone[[i]], Peptide2 == 1)
  linktypelisthomeotypicone[[i]]<-homeotypiclistone[[i]]$conresid
}

#combine remaining type one and type rest
oneandrest<-list()
oneandrest<-mapply(c,linktypelistoneres, linktypelistrestres, SIMPLIFY=FALSE)
homeo<-mapply(c,linktypelisthomeotypic,linktypelisthomeotypicone)
#barplot of crosslink types
p<-list(NA)
crosslinktypeplot(p,linktypesres)

#combined with one(homeotypic removed)
O<- list(
  "37" = oneandrest[[1]], 
  "47" = oneandrest[[2]], 
  "57" = oneandrest[[3]],
  "67" = oneandrest[[4]]
)



#Ven diagram

forven<-list(O)


crosslinkven4(forven)

#uniquelists<-list()
barplotsres<-list()
linktypeslistres<-list()
crosslinklistres<-list()
venslistres<-list()
savedatares(p,linktypesres,Crosslinkclasses,VENS,1)

test1<-RTC_TEMP[[4]] %>% drop_na(conres)
test2<-checkrev(test1)
test67<-table(test2$conresid) 
getwd()
setwd("C:\\Users\\victo\\Desktop\\speciale\\Data\\Data for artikel")
writexl::write_xlsx(list("37" =as.data.frame(test37)  ,"47"=as.data.frame(test47),"57"=as.data.frame(test57),"67"=as.data.frame(test67)),path = "tempexperiment.xlsx")

####### RTC4 1_40 ###############################################
#concatenate peptide 1 and peptide 2 so we known where links are between
conpepres2(RTC410_40)
RTC410_40<-listdata


dataset<-list(RTC410_40[[1]],RTC410_40[[2]],RTC410_40[[3]])


linklist<-c("1:10","1:20","1:40")

linktypelistzero<-list(NA,NA,NA,NA,NA)
linktypelistone<-list(NA,NA,NA,NA,NA,NA)
linktypelistrest<-list(NA,NA,NA,NA,NA,NA)

linktypelistzerores<-list(NA,NA)
linktypelistoneres<-list(NA,NA)
linktypelistrestres<-list(NA,NA)


crosslinkclassifyresid2(dataset = dataset,
                        linklist = linklist,
                        linktypelistzero,
                        linktypelistone,
                        linktypelistrest)

#filter linktypelistrestres filtering 
linktypelisthomeotypic<-list(NA,NA)
datares<-list(NA,NA)
homeotypiclist<-list(NA,NA)

checkrest(linktypelistrestres[1:2])

linktypelistrestres<-datares  

NAremoveres(linktypelistrestres,linktypelisthomeotypic)

#filter linktypelistoneres filtering 
linktypelisthomeotypicone<-list(NA,NA,NA,NA,NA,NA)
dataone<-list(NA,NA,NA,NA,NA)
homeotypiclistone<-list(NA,NA,NA,NA,NA,NA)

checkones(linktypelistoneres[1:3])

linktypelistoneres<-dataone

NAremoveresone(linktypelistoneres,linktypelisthomeotypicone)


#filter the homeotypic peptides from original data
#can be used for either rest(homeotypic) and one(to find S and Y links not filtered by meroxsort)
findhomeo(dataset,linktypelisthomeotypic)
for(i in 1:length(homeotypiclist)){
  homeotypiclist[[i]]<-filter(homeotypiclist[[i]], Peptide2 != 1)
  linktypelisthomeotypic[[i]]<-homeotypiclist[[i]]$conresid
}

#Run each time to remove false homeotypic (S and Y)
findhomeoone(dataset,linktypelisthomeotypicone)
for(i in 1:length(homeotypiclistone)){
  homeotypiclistone[[i]]<-filter(homeotypiclistone[[i]], Peptide2 == 1)
  linktypelisthomeotypicone[[i]]<-homeotypiclistone[[i]]$conresid
}

#combine remaining type one and type rest

oneandrest<-list()
oneandrest<-mapply(c,linktypelistoneres, linktypelistrestres,linktypelistzerores ,SIMPLIFY=FALSE)
oneandrest<-mapply(c,linktypelistoneres, linktypelistrestres,SIMPLIFY=FALSE)

#DOUBLE CHECK TO SEE IF HOMEOTYPIC STILL ARE PRESENT
for(i in 1:length(oneandrest[[1]])){
  if(oneandrest[[1]][i]=="106|106|>light|>light" ){
    print("YES")
  }
}

homeo<-mapply(c,linktypelisthomeotypic,linktypelisthomeotypicone)
#barplot of crosslink types
p<-list(NA)
crosslinktypeplot(p,linktypesres)


#look at different types
#C4

#combined with one(homeotypic removed)
O<- list(
  A = na.omit(oneandrest[[1]]), 
  B = na.omit(oneandrest[[2]]), 
  C = na.omit(oneandrest[[3]])
)

TY <- list(
  A = na.omit(oneandrest[[4]]), 
  B = na.omit(oneandrest[[5]]), 
  C = na.omit(oneandrest[[6]])
)

#Used to filter data using T5
O2<- list(
  A = setmincount(na.omit(oneandrest[[1]]),2), 
  B = setmincount(na.omit(oneandrest[[2]]),2), 
  C = setmincount(na.omit(oneandrest[[3]]),2)
)

TY2 <- list(
  A = setmincount(na.omit(oneandrest[[4]]),2), 
  B = setmincount(na.omit(oneandrest[[5]]),2), 
  C = setmincount(na.omit(oneandrest[[6]]),2)
)


O3<-ldply(O,data.frame)
TY3<-ldply(TY,data.frame)
O3include<-setmincount_merged(O3,3)
TY3include<-setmincount_merged(TY3,3)

#5 = 82.7 68.7
O3<- list(
  A = na.omit(oneandrest[[1]][oneandrest[[1]] %in% O3include$AA]), 
  B = na.omit(oneandrest[[2]][oneandrest[[2]] %in% O3include$AA]), 
  C = na.omit(oneandrest[[3]][oneandrest[[3]] %in% O3include$AA])
)

TY3 <- list(
  A = na.omit(oneandrest[[4]][oneandrest[[4]] %in% TY3include$AA]), 
  B = na.omit(oneandrest[[5]][oneandrest[[5]] %in% TY3include$AA]), 
  C = na.omit(oneandrest[[6]][oneandrest[[6]] %in% TY3include$AA])
)


#combined with one(homeotypic removed)
Oh<- list(
  A = na.omit(homeo[[1]]), 
  B = na.omit(homeo[[2]]), 
  C = na.omit(homeo[[3]])
)

TYh <- list(
  A = na.omit(homeo[[4]]), 
  B = na.omit(homeo[[5]]), 
  C = na.omit(homeo[[6]])
)


setmincount<-function(data,min_count){
  dat<-as.data.frame(data)
  colnames(dat)<-c("AA")
  test<-as.data.frame(table(dat))
  
  test2<-dplyr::filter(test, Freq > min_count)
  
  test3<-dat[dat$AA %in% test2$dat,]
  return(test3)
}

setmincount_merged<-function(data,min_count){
  dat<-as.data.frame(data)
  colnames(dat)<-c("ID","AA")
  test<-as.data.frame(table(dat$AA))
  
  test2<-dplyr::filter(test, Freq > min_count)
  
  test3<-unique(dat[dat$AA %in% test2$Var1,])
  return(test3)
}

#Ven diagram

forven<-list(TY2,O2)


crosslinkven(forven)


#uniquelists<-list()
barplotsres<-list()
linktypeslistres<-list()
crosslinklistres<-list()
venslistres<-list()
savedatares(p,linktypesres,Crosslinkclasses,VENS,1)




####### IGA ##########################################################
#concatenate peptide 1 and peptide 2 so we known where links are between
conpepres2(dataset)
dataset<-listdata


setwd("C:\\Users\\victo\\Desktop\\speciale\\R\\MEROXSORT\\MeroxAdjust")
IGA250<-read_xlsx("final_IGA250_dataset.xlsx")
IGA500<-read_xlsx("final_IGA500_dataset.xlsx")
IGA750<-read_xlsx("final_IGA750_dataset.xlsx")
dataset<-list(IGA250,IGA500,IGA750)
linklist<-c("IGA 1:250","IGA 1:500","IGA 1:750")

for(i in 1:length(dataset)){
  for(j in 1:dim(dataset[[i]])[1]){
    if(dataset[[i]]$Protein.1[j] == ">1IGA_1|Chains A, B|IGA1|Homo sapiens  9606"){
      dataset[[i]]$Protein.1[j] = ">heavy"
    }
    if(dataset[[i]]$Protein.2[j] == ">1IGA_1|Chains A, B|IGA1|Homo sapiens  9606"){
      dataset[[i]]$Protein.2[j] = ">heavy"
    }
    if(dataset[[i]]$Protein.1[j] == ">1IGA_2|Chains C, D|IGA1|Homo sapiens  9606"){
      dataset[[i]]$Protein.1[j] = ">light"
    }
    if(dataset[[i]]$Protein.2[j] == ">1IGA_2|Chains C, D|IGA1|Homo sapiens  9606"){
      dataset[[i]]$Protein.2[j] = ">light"
    }
  }
}

#linklist<-c("1:250A","1:250B","1:250C","1:500A","1:500B","1:500C","1:750A","1:750B","1:750C")

linktypelistzero<-list(NA,NA,NA)
linktypelistone<-list(NA,NA,NA)
linktypelistrest<-list(NA,NA,NA)

linktypelistzerores<-list(NA,NA,NA)
linktypelistoneres<-list(NA,NA,NA)
linktypelistrestres<-list(NA,NA,NA)




crosslinkclassifyresid2(dataset = dataset,
                        linklist = linklist,
                        linktypelistzero,
                        linktypelistone,
                        linktypelistrest)

#filter linktypelistrestres filtering 
linktypelisthomeotypic<-list(NA,NA)
datares<-list(NA,NA)
homeotypiclist<-list(NA,NA)

checkrest(linktypelistrestres[1:2])

linktypelistrestres<-datares  

NAremoveres(linktypelistrestres,linktypelisthomeotypic)

#filter linktypelistoneres filtering 
linktypelisthomeotypicone<-list(NA,NA,NA,NA,NA,NA,NA,NA,NA)
dataone<-list(NA,NA,NA,NA,NA,NA,NA,NA,NA)
homeotypiclistone<-list(NA,NA,NA,NA,NA,NA,NA,NA,NA)

checkones(linktypelistoneres)

linktypelistoneres<-dataone

NAremoveresone(linktypelistoneres,linktypelisthomeotypicone)


#filter the homeotypic peptides from original data
#can be used for either rest(homeotypic) and one(to find S and Y links not filtered by meroxsort)
findhomeo(dataset,linktypelisthomeotypic)
for(i in 1:length(homeotypiclist)){
  homeotypiclist[[i]]<-filter(homeotypiclist[[i]], Peptide2 != 1)
  linktypelisthomeotypic[[i]]<-homeotypiclist[[i]]$conresid
}

#Run each time to remove false homeotypic (S and Y)
findhomeoone(dataset,linktypelisthomeotypicone)
for(i in 1:length(homeotypiclistone)){
  homeotypiclistone[[i]]<-filter(homeotypiclistone[[i]], Peptide2 == 1)
  linktypelisthomeotypicone[[i]]<-homeotypiclistone[[i]]$conresid
}

#combine remaining type one and type rest

oneandrest<-list()
oneandrest<-mapply(c,linktypelistoneres, linktypelistrestres,linktypelistzerores ,SIMPLIFY=FALSE)
oneandrest<-mapply(c,linktypelistoneres, linktypelistrestres,SIMPLIFY=FALSE)



homeo<-mapply(c,linktypelisthomeotypic,linktypelisthomeotypicone)
#barplot of crosslink types
p<-list(NA)
crosslinktypeplot(p,linktypesres)


#look at different types
#C4

#combined with one(homeotypic removed)
A250<- list(
  A = na.omit(oneandrest[[1]]), 
  B = na.omit(oneandrest[[2]]), 
  C = na.omit(oneandrest[[3]])
)

B500 <- list(
  A = na.omit(oneandrest[[4]]), 
  B = na.omit(oneandrest[[5]]), 
  C = na.omit(oneandrest[[6]])
)

C750 <- list(
  A = na.omit(oneandrest[[7]]), 
  B = na.omit(oneandrest[[8]]), 
  C = na.omit(oneandrest[[9]])
)

#Used to filter data using T5
A250<- list(
  A = setmincount(na.omit(oneandrest[[1]]),2), 
  B = setmincount(na.omit(oneandrest[[2]]),2), 
  C = setmincount(na.omit(oneandrest[[3]]),2)
)

B500<- list(
  A = setmincount(na.omit(oneandrest[[4]]),2), 
  B = setmincount(na.omit(oneandrest[[5]]),2), 
  C = setmincount(na.omit(oneandrest[[6]]),2)
)

C750<- list(
  A = setmincount(na.omit(oneandrest[[7]]),2), 
  B = setmincount(na.omit(oneandrest[[8]]),2), 
  C = setmincount(na.omit(oneandrest[[9]]),2)
)


COMPARE<- list(
  A = na.omit(oneandrest[[1]]), 
  B = na.omit(oneandrest[[2]]), 
  C = na.omit(oneandrest[[3]])
)
#Ven diagram

forven<-list(A250,B500,C750)
forven<-list(COMPARE)


crosslinkven(forven)


#uniquelists<-list()
barplotsres<-list()
linktypeslistres<-list()
crosslinklistres<-list()
venslistres<-list()
savedatares(p,linktypesres,Crosslinkclasses,VENS,1)



secs<-list(NA,NA,NA)


secs[[1]]<-Crosslinkclasses[[1]]$..values..[2]
secs[[2]]<-Crosslinkclasses[[1]]$..values..[3]
secs[[3]]<-Crosslinkclasses[[1]]$..values..[5]
secs[[4]]<-Crosslinkclasses[[1]]$..values..[1]


secs[[5]]<-Crosslinkclasses[[2]]$..values..[2]
secs[[6]]<-Crosslinkclasses[[2]]$..values..[3]
secs[[7]]<-Crosslinkclasses[[2]]$..values..[5]
secs[[8]]<-Crosslinkclasses[[2]]$..values..[1]


secs[[9]]<-Crosslinkclasses[[3]]$..values..[2]
secs[[10]]<-Crosslinkclasses[[3]]$..values..[3]
secs[[11]]<-Crosslinkclasses[[3]]$..values..[5]
secs[[12]]<-Crosslinkclasses[[3]]$..values..[1]


joined250<-c(secs[[1]]$`2`,secs[[2]]$`3`,secs[[3]]$`5`,secs[[4]]$`1`)
binded<-ldply(dataset[1:3],data.frame)
joined250<-binded[binded$conresid %in% joined250,]

joined500<-c(secs[[5]]$`2`,secs[[6]]$`3`,secs[[7]]$`5`,secs[[8]]$`1`)
binded<-ldply(dataset[4:6],data.frame)
joined500<-binded[binded$conresid %in% joined500,]

joined750<-c(secs[[9]]$`2`,secs[[10]]$`3`,secs[[11]]$`5`,secs[[12]]$`1`)
binded<-ldply(dataset[7:9],data.frame)
joined750<-binded[binded$conresid %in% joined750,]

writexl::write_xlsx(joined250,path = "final_IGA250_dataset.xlsx")
writexl::write_xlsx(joined500,path = "final_IGA500_dataset.xlsx")
writexl::write_xlsx(joined750,path = "final_IGA750_dataset.xlsx")




