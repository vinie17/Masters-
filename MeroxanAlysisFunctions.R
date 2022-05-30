library(readxl)
library(ggvenn)
library(tidyverse)
library(ggplot2)
library(stringr)
library(VennDiagram)
library(plyr)
library(gridExtra)
#library(ggVennDiagram)
#library(Rcpp)
#library(eulerr)
#library(QCA)

#get the unique peptides in a,b,c
#uniquelist<-list(NA,NA)
#getsetdiff(forven[2])

#Get uniqe for each set
getsetdiff<-function(listdat){
 for(i in 1:length(listdat)){
  DIFF1<-outersect(listdat[[i]][[1]],listdat[[i]][[2]])
  #A<-outersect(DIFF1,listdat[[i]][[3]])
  DIFF2<-outersect(listdat[[i]][[1]],listdat[[i]][[3]])
  #C<-outersect(DIFF2,listdat[[i]][[1]])
  DIFF3<-outersect(listdat[[i]][[2]],listdat[[i]][[3]])
  #B<-setdiff(DIFF3,listdat[[i]][[2]])
  #abc<-ldply(A,B,C, data.frame)
  #uniquelist[i]<-abc
  DIFF12<-outersect(DIFF1,DIFF2)
  DIFF13<-outersect(DIFF1,DIFF3)
  DIFF23<-outersect(DIFF2,DIFF3)
  
  #DIFF1213<-outersect(DIFF12,DIFF13)
  
  
 }
  A<<-DIFF1
  B<<-DIFF2
  CC<<-DIFF3
  #return(uniquelist)
}

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}


#Extract rows containing peptides with C-terminal K best linkage position
checkk<-function(datlist){
  Klist<-list()
  for(i in 1:length(datlist)){
    CtermK<-data.frame()
    for(j in 1:dim(datlist[[i]])[1]){
      if((str_split(datlist[[i]]$Peptide.1[j],""))[[1]][str_count(datlist[[i]]$Peptide.1[j])-1] == "K"){
        if((str_count(datlist[[i]]$Peptide.1[j])-2==(as.numeric(str_split(datlist[[i]]$best.linkage.position.peptide.1[j],"")[[1]][2])))){
          CtermK[j,1:33]<-datlist[[i]][j,]
          datlist[[i]]$Peptide.1[j]<-NA
        }
      }else{
        if(datlist[[i]]$Peptide2[j] == "0"||datlist[[i]]$Peptide2[j] == "1"){next}
        
        if((str_split(datlist[[i]]$Peptide2[j],""))[[1]][str_count(datlist[[i]]$Peptide2[j])-1] == "K"){
          if((str_count(datlist[[i]]$Peptide2[j])-2==(as.numeric(str_split(datlist[[i]]$best.linkage.position.peptide.2[j],"")[[1]][2])))){
            CtermK[j,]<-datlist[[i]][j,]
            datlist[[i]]$Peptide2[j]<-NA
          }
        }
      }
    }
    CtermK<-CtermK[!is.na(CtermK$Score),]
    datlist[[i]]<-datlist[[i]][!is.na(datlist[[i]]$Peptide.1),]
    datlist[[i]]<-datlist[[i]][!is.na(datlist[[i]]$Peptide2),]
    Klist[[i]]<-CtermK
  }
  Klist<<-Klist
  return(datlist)
} 

#check pepties for K position
#checkk<-function(datlist){
 # Klist<-list()
  #for(i in 1:length(datlist)){
#NtermK<-data.frame()
 #   for(j in 1:dim(datlist[[i]])[1]){
  #    if((str_split(datlist[[i]]$Peptide.1[j],""))[[1]][str_count(datlist[[i]]$Peptide.1[j])-1] == "K"){
   #     NtermK[j,1:33]<-datlist[[i]][j,]
    #    datlist[[i]]$Peptide.1[j]<-NA
     # }else{
      #if(datlist[[i]]$Peptide2[j] == "0"||datlist[[i]]$Peptide2[j] == "1"){next}
      #
#      if((str_split(datlist[[i]]$Peptide2[j],""))[[1]][str_count(datlist[[i]]$Peptide2[j])-1] == "K"){
 #       NtermK[j,]<-datlist[[i]][j,]
  #      datlist[[i]]$Peptide2[j]<-NA
   #   }
    #  }
#    }
#NtermK<-NtermK[!is.na(NtermK$Score),]
#datlist[[i]]<-datlist[[i]][!is.na(datlist[[i]]$Peptide.1),]
#datlist[[i]]<-datlist[[i]][!is.na(datlist[[i]]$Peptide2),]
#Klist[[i]]<-NtermK
#  }
 # Klist<<-Klist
#  datlist<<-datlist
#} 
  

#concatenate peptide 1 and peptide 2 so we known where links are between ####

conpep<-function(listdata){
  for(i in 1:length(listdata)){
    for(j in 1:dim(listdata[[i]])[1]){
      listdata[[i]]$concpep[j]<-str_c(listdata[[i]]$Peptide.1[j],listdata[[i]]$Peptide2[j],collapse = ",")
      listdata[[i]]$concpeprev[j]<-str_c(listdata[[i]]$Peptide2[j],listdata[[i]]$Peptide.1[j],collapse = ",")
      id1<-str_split(listdata[[i]]$best.linkage.position.peptide.1[j],"")[[1]][2]
      id2<-str_split(listdata[[i]]$best.linkage.position.peptide.2[j],"")[[1]][2]
      if(listdata[[i]]$Peptide2[j] == 1){
        listdata[[i]]$conresid[j]<-str_c(str_c(as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,as.numeric(listdata[[i]]$From[j])+as.numeric(id2)-1,collapse = ",",sep = "|"),
                                         str_c(listdata[[i]]$Protein.1[j],
                                               listdata[[i]]$Protein.1[j],
                                               collapse = ",",sep = "|"),collapse = "," ,sep= "|")
      }else if(listdata[[i]]$Peptide2[j] == 0){
        listdata[[i]]$conresid[j]<-str_c(as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,"0",collapse = ",",sep = "|")
         }else{
               listdata[[i]]$conresid[j]<-str_c(
                 str_c(
                   as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,
                   as.numeric(listdata[[i]]$From.1[j])+as.numeric(id2)-1,
                   collapse = ",",sep = "|"),
                 str_c(listdata[[i]]$Protein.1[j],
                 listdata[[i]]$Protein.2[j],
                 collapse = ",",sep = "|"),collapse = "," ,sep= "|")
               }
    }
  }
  listdata<<-listdata
}
#str_c(RTC4[[1]]$Peptide2$best.linkage.position.peptide.1[1],RTC4[[1]]$best.linkage.position.peptide.2[1],collapse = ",")
#check for reversed resids 

#tes<-test[!duplicated(test$conresid),]

#checkrev(tes)
checkrev<-function(binded1){
#count=0
  for(j in 1:dim(binded1)[1]){
    for(i in j:dim(binded1)[1]){
      if(is.na(binded1$conresrev[i])){next}
      else if(j == i){next}
      else if(binded1$Peptide2[j] == 0){next}
      else if(binded1$conres[j] == binded1$conresrev[i]){
        #count=count+1
        #listrev[count]<-binded$conresid[i]
        binded1$conresid[i]<-binded1$conresid[j]
        #binded1$conresrev[i]<-NA
        
      }
    }}

return(binded1)
}


checkrev2<-function(binded1){
  #count=0
  for(j in 1:dim(binded1)[1]){
    for(i in j:dim(binded1)[1]){
      if(is.na(binded1$conresrev[i])){next}
      else if(j == i){next}
      else if(binded1$conres[j] == binded1$conresrev[i]){
        #count=count+1
        #listrev[count]<-binded$conresid[i]
        binded1$conresrev[i]<-binded1$conres[j]
        #binded1$conresrev[i]<-NA
        #print("yes")
        
      }
    }}
  
  return(binded1)
}



#listdata<-testlist
conpepres<-function(listdata){
  for(i in 1:length(listdata)){
    for(j in 1:dim(listdata[[i]])[1]){
      listdata[[i]]$concpep[j]<-str_c(listdata[[i]]$Peptide.1[j],listdata[[i]]$Peptide2[j],collapse = ",")
      listdata[[i]]$concpeprev[j]<-str_c(listdata[[i]]$Peptide2[j],listdata[[i]]$Peptide.1[j],collapse = ",")
      if(length(str_split(listdata[[i]]$best.linkage.position.peptide.1[j],"")[[1]])>2){
        id1<-str_c(str_split(listdata[[i]]$best.linkage.position.peptide.1[j],"")[[1]][2:3],collapse = "")
      }else{
        id1<-str_split(listdata[[i]]$best.linkage.position.peptide.1[j],"")[[1]][2]
      }
      if(length(str_split(listdata[[i]]$best.linkage.position.peptide.2[j],"")[[1]])>2){
        id2<-str_c(str_split(listdata[[i]]$best.linkage.position.peptide.2[j],"")[[1]][2:3],collapse = "")}
      else{
      id2<-str_split(listdata[[i]]$best.linkage.position.peptide.2[j],"")[[1]][2]
      }
      if(listdata[[i]]$Peptide2[j] == 1){
        listdata[[i]]$conres[j]<-str_c(as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,as.numeric(listdata[[i]]$From[j])+as.numeric(id2)-1,collapse = ",",sep = "|")
        listdata[[i]]$conresrev[j]<-str_c(as.numeric(listdata[[i]]$From[j])+as.numeric(id2)-1,as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,collapse = ",",sep = "|")
        listdata[[i]]$conresid[j]<-str_c(str_c(as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,as.numeric(listdata[[i]]$From[j])+as.numeric(id2)-1,collapse = ",",sep = "|"),
                                         str_c(listdata[[i]]$Protein.1[j],
                                               listdata[[i]]$Protein.1[j],
                                               collapse = ",",sep = "|"),collapse = "," ,sep= "|")
      }else if(listdata[[i]]$Peptide2[j] == 0){
        listdata[[i]]$conresid[j]<-str_c(as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,"0",listdata[[i]]$Protein.1[j],collapse = ",",sep = "|")
        listdata[[i]]$conres[j]<-NA
        listdata[[i]]$conresrev[j]<-NA
        
      }else{
        listdata[[i]]$conresid[j]<-str_c(
          str_c(
            as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,
            as.numeric(listdata[[i]]$From.1[j])+as.numeric(id2)-1,
            collapse = ",",sep = "|"),
          str_c(listdata[[i]]$Protein.1[j],
                listdata[[i]]$Protein.2[j],
                collapse = ",",sep = "|"),collapse = "," ,sep= "|")
        listdata[[i]]$conres[j]<-str_c(
          as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,
          as.numeric(listdata[[i]]$From.1[j])+as.numeric(id2)-1,
          collapse = ",",sep = "|")
        listdata[[i]]$conresrev[j]<-str_c(
          as.numeric(listdata[[i]]$From.1[j])+as.numeric(id2)-1,
          as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,
          collapse = ",",sep = "|")
      }
    }
  }
  listdata<<-listdata
}

conpepres2<-function(listdata){
  for(i in 1:length(listdata)){
    for(j in 1:dim(listdata[[i]])[1]){
      listdata[[i]]$concpep[j]<-str_c(listdata[[i]]$Peptide.1[j],listdata[[i]]$Peptide2[j],collapse = ",")
      listdata[[i]]$concpeprev[j]<-str_c(listdata[[i]]$Peptide2[j],listdata[[i]]$Peptide.1[j],collapse = ",")
        id1<-str_c(substring(listdata[[i]]$best.linkage.position.peptide.1[j],2),collapse = "")
        id2<-str_c(substring(listdata[[i]]$best.linkage.position.peptide.2[j],2),collapse = "")
    
      if(listdata[[i]]$Peptide2[j] == 1){
        listdata[[i]]$conres[j]<-str_c(as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,as.numeric(listdata[[i]]$From[j])+as.numeric(id2)-1,collapse = ",",sep = "|")
        listdata[[i]]$conresrev[j]<-str_c(as.numeric(listdata[[i]]$From[j])+as.numeric(id2)-1,as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,collapse = ",",sep = "|")
        listdata[[i]]$conresid[j]<-str_c(str_c(as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,as.numeric(listdata[[i]]$From[j])+as.numeric(id2)-1,collapse = ",",sep = "|"),
                                         str_c(listdata[[i]]$Protein.1[j],
                                               listdata[[i]]$Protein.1[j],
                                               collapse = ",",sep = "|"),collapse = "," ,sep= "|")
      }else if(listdata[[i]]$Peptide2[j] == 0){
        listdata[[i]]$conresid[j]<-str_c(as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,"0",listdata[[i]]$Protein.1[j],collapse = ",",sep = "|")
        listdata[[i]]$conres[j]<-NA
        listdata[[i]]$conresrev[j]<-NA
        
      }else{
        listdata[[i]]$conresid[j]<-str_c(
          str_c(
            as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,
            as.numeric(listdata[[i]]$From.1[j])+as.numeric(id2)-1,
            collapse = ",",sep = "|"),
          str_c(listdata[[i]]$Protein.1[j],
                listdata[[i]]$Protein.2[j],
                collapse = ",",sep = "|"),collapse = "," ,sep= "|")
        listdata[[i]]$conres[j]<-str_c(
          as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,
          as.numeric(listdata[[i]]$From.1[j])+as.numeric(id2)-1,
          collapse = ",",sep = "|")
        listdata[[i]]$conresrev[j]<-str_c(
          as.numeric(listdata[[i]]$From.1[j])+as.numeric(id2)-1,
          as.numeric(listdata[[i]]$From[j])+as.numeric(id1)-1,
          collapse = ",",sep = "|")
      }
    }
  }
  listdata<<-listdata
  #return(listdata)
  
}


#split depending on type of cross-link and count total number ####

#dataset1<-list(RTC4[[1]],RTC4[[2]],RTC4[[3]],RTC4[[4]],RTC4[[5]],RTC4[[6]])
#linklist1<-c("OA","OB","OC","TYA","TYB","TYC")
#linktypelistzero<-list(NA,NA,NA,NA,NA,NA)
#linktypelistone<-list(NA,NA,NA,NA,NA,NA)
#inktypelistrest<-list(NA,NA,NA,NA,NA,NA)
#crosslinkclassify(dataset = dataset1, linklist = linklist1,linktypelistzero, linktypelistone,linktypelistrest)
#input:
#dataset: list of lists to dataframes to be compared
#Linklist: list of what is in each dataframe same order
#Linktypelistzero/one/rest: list of NA length same as linklist

crosslinkclassify<-function(dataset1,linklist1,linktypelistzero,linktypelistone,linktypelistrest){
  linktypes<-data.frame()
for(i in 1:length(dataset)){
  #RTC4[[i]]<-RTC4[[i]][!is.na(RTC4[[i]]$Peptide2), ]
  zero<-0
  one<-0
  rest<-0
  for(k in 1:dim(dataset[[i]])[1]){
    try(if(dataset[[i]]$Peptide2[k] == 0){
      zero<-zero+1
      linktypelistzero[[i]][k]<-dataset[[i]]$concpep[k]
    } else if(dataset[[i]]$Peptide2[k] == 1){
      one<-one+1
      linktypelistone[[i]][k]<-dataset[[i]]$concpep[k]
    } else{
      rest<-rest+1
      linktypelistrest[[i]][k]<-dataset[[i]]$concpep[k]
    },silent = F)
    
  }
  linktypes[i,1]<-zero
  linktypes[i,2]<-one
  linktypes[i,3]<-rest
  linktypes[i,4]<-linklist[i]
  
}
#colnames for linktypes 
colnames(linktypes)<-c("zero","one","rest","File")

linktypelistone<<-linktypelistone
linktypelistrest<<-linktypelistrest
linktypelistzero<<-linktypelistzero
linktypes<<-linktypes
}


#FOR RESID
crosslinkclassifyresid<-function(dataset1,linklist1,linktypelistzero,linktypelistone,linktypelistrest){
  linktypes<-data.frame()
  for(i in 1:length(dataset)){
    #RTC4[[i]]<-RTC4[[i]][!is.na(RTC4[[i]]$Peptide2), ]
    zero<-0
    one<-0
    rest<-0
    for(k in 1:dim(dataset[[i]])[1]){
      if(dataset[[i]]$Peptide2[k] == 0){
        zero<-zero+1
        linktypelistzero[[i]][k]<-dataset[[i]]$conresid[k]
      } else if(dataset[[i]]$Peptide2[k] == 1){
        one<-one+1
        linktypelistone[[i]][k]<-dataset[[i]]$conresid[k]
      } else{
        rest<-rest+1
        linktypelistrest[[i]][k]<-dataset[[i]]$conresid[k]
      }
      
    }
    linktypes[i,1]<-zero
    linktypes[i,2]<-one
    linktypes[i,3]<-rest
    linktypes[i,4]<-linklist[i]
    
  }
  #colnames for linktypes 
  colnames(linktypes)<-c("zero","one","rest","File")
  
  linktypelistoneres<<-linktypelistone
  linktypelistrestres<<-linktypelistrest
  linktypelistzerores<<-linktypelistzero
  linktypesres<<-linktypes
}

crosslinkclassifyresid2<-function(dataset,linklist,linktypelistzero,linktypelistone,linktypelistrest){
  linktypes<-data.frame()
  for(i in 1:length(dataset)){
        linktypelistzero[[i]]<-filter(dataset[[i]], Peptide2 ==0)$conresid
        linktypelistone[[i]]<-filter(dataset[[i]], Peptide2 ==1)$conresid
        holdrest<-filter(dataset[[i]], Peptide2 !=1)
        linktypelistrest[[i]]<-filter(holdrest, Peptide2 !=0)$conresid
    
    linktypes[i,1]<-length(linktypelistzero[[i]])
    linktypes[i,2]<-length(linktypelistone[[i]])
    linktypes[i,3]<-length(linktypelistrest[[i]])
    linktypes[i,4]<-linklist[i]
    
  }
  #colnames for linktypes 
  colnames(linktypes)<-c("zero","one","rest","File")
  
  linktypelistoneres<<-linktypelistone
  linktypelistrestres<<-linktypelistrest
  linktypelistzerores<<-linktypelistzero
  linktypesres<<-linktypes
}

#remove NA from linktype lists ####
NAremove<-function(linktypelistone,linktypelistzero,linktypelistrest){
for(i in 1:length(linktypelistone)){
  linktypelistone[[i]]<-na.omit(linktypelistone[[i]])
}
for(i in 1:length(linktypelistzero)){
  linktypelistzero[[i]]<-na.omit(linktypelistzero[[i]])
}
for(i in 1:length(linktypelistrest)){
  linktypelistrest[[i]]<-na.omit(linktypelistrest[[i]])
}
linktypelistoneres<<-linktypelistone
linktypelistzerores<<-linktypelistzero
linktypelistrestres<<-linktypelistrest
}

#remove NA from linktype RES lists ####
NAremoveres<-function(linktypelistrestres,linktypelisthomeotypic){
  for(i in 1:length(linktypelistrestres)){
    linktypelistrestres[[i]]<-na.omit(linktypelistrestres[[i]])
  }
  for(i in 1:length(linktypelisthomeotypic)){
    linktypelisthomeotypic[[i]]<-na.omit(linktypelisthomeotypic[[i]])
  }
 
  linktypelistrestres<<-linktypelistrestres
  linktypelisthomeotypic<<-linktypelisthomeotypic
  
}

NAremoveresone<-function(linktypelistoneres,linktypelisthomeotypicone){
  for(i in 1:length(linktypelistoneres)){
    linktypelistoneres[[i]]<-na.omit(linktypelistoneres[[i]])
  }
  for(i in 1:length(linktypelisthomeotypicone)){
    linktypelisthomeotypicone[[i]]<-na.omit(linktypelisthomeotypicone[[i]])
  }
  
  linktypelistoneres<<-linktypelistoneres
  linktypelisthomeotypicone<<-linktypelisthomeotypicone
  
}

#check for homeotypic peptides in linkone ####
checkones<-function(dataone){
  for(i in 1:length(dataone)){
    for(j in 1:length(dataone[[i]])){
      if(str_split(dataone[[i]][j],"\\|")[[1]][1] == str_split(dataone[[i]][j],"\\|")[[1]][2]){
        linktypelisthomeotypicone[[i]][j]<-dataone[[i]][j]
        dataone[[i]][j]<-NA
      }else{
        linktypelisthomeotypicone[[i]][j]<-NA
      }
    }
  }
  linktypelisthomeotypicone<<-linktypelisthomeotypicone
  dataone<<-dataone
}
checkrest<-function(datares){
  for(i in 1:length(datares)){
    for(j in 1:length(datares[[i]])){
      #print(datares[[i]][j])
      if(str_split(datares[[i]][j],"\\|")[[1]][1] == str_split(datares[[i]][j],"\\|")[[1]][2]){
        linktypelisthomeotypic[[i]][j]<-datares[[i]][j]
        datares[[i]][j]<-NA
      }else{
        linktypelisthomeotypic[[i]][j]<-NA
      }
    }
  }
  linktypelisthomeotypic<<-linktypelisthomeotypic
  datares<<-datares
}


#filter homeotypic from original data ####
findhomeo<-function(listdat,reflist){
  testlist<-list()
  for(j in 1:length(listdat)){
    if(is_empty(reflist[[j]])){}else{
    
    test<-data.frame()
    test2<-data.frame()
    uniqueones<-unique(reflist[[j]])
    for(i in 1:length(uniqueones)){
      test<-filter(listdat[[j]],conresid == uniqueones[i])
      test2<-rbind(test2,test)
      
    }
    testlist[[j]]<-test2
    }
  }
  homeotypiclist<<-testlist
}



findhomeoone<-function(listdat,reflist){
  testlist<-list()
  for(j in 1:length(listdat)){
    for(j in 1:length(listdat)){
      if(is_empty(reflist[[j]])){}else{
    
    test<-data.frame()
    test2<-data.frame()
    uniqueones<-unique(reflist[[j]])
    for(i in 1:length(uniqueones)){
      test<-filter(listdat[[j]],conresid == uniqueones[i])
      test2<-rbind(test2,test)
      
    }
    testlist[[j]]<-test2
      }
    }
  }
  homeotypiclistone<<-testlist
}



findrest<-function(listdat,reflist){
  testlist<-list()
  for(j in 1:length(listdat)){
    for(j in 1:length(listdat)){
      if(is_empty(reflist[[j]])){}else{
        
        test<-data.frame()
        test2<-data.frame()
        uniqueones<-unique(reflist[[j]])
        for(i in 1:length(uniqueones)){
          test<-filter(listdat[[j]],conresid == uniqueones[i])
          test2<-rbind(test2,test)
          
        }
        testlist[[j]]<-test2
      }
    }
  }
  restlist<<-testlist
}
  
#barplot of crosslink types ####
#list of NA,NA,NA
#p<-list(NA,NA,NA)


crosslinktypeplot<-function(p,linktypes){
p[[1]]<-ggplot(data = linktypes, aes(x=File,y=zero, fill = File))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("skyblue","lightgreen","skyblue","lightgreen","lightgreen","lightgreen","grey70","grey80","grey90"))+
  ylab("Type 0")+theme(legend.position = "none")
#("skyblue","lightgreen","coral","#C77CFF","yellow","lightgreen","grey70","grey80","grey90")

p[[2]]<-ggplot(data = linktypes, aes(x=File,y=one, fill = File))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("skyblue","lightgreen","skyblue","lightgreen","lightgreen","lightgreen","grey70","grey80","grey90"))+
  ylab("Type 1")+theme(legend.position = "none")

p[[3]]<-ggplot(data = linktypes, aes(x=File,y=rest, fill = File))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("skyblue","lightgreen","skyblue","lightgreen","lightgreen","lightgreen","grey70","grey80","grey90"))+
  ylab("Type 2")+theme(legend.position = "none")
pp<-grid.arrange(p[[1]],p[[2]],p[[3]],ncol = 3)
p<<-pp
}

getwd()


setmincount<-function(data,min_count){
  dat<-as.data.frame(data)
  colnames(dat)<-c("AA")
  test<-as.data.frame(table(dat))
  
  test2<-dplyr::filter(test, Freq > min_count)
  
  test3<-dat[dat$AA %in% test2$dat,]
  return(test3)
}


#Ven diagrams ####
#input:
#list of lists to be compared
#O<- list(
#  A = linktypelistrest[[1]], 
#  B = linktypelistrest[[2]], 
#  C = linktypelistrest[[3]]
#)

#TY <- list(
#  A = linktypelistrest[[4]], 
#  B = linktypelistrest[[5]], 
#  C = linktypelistrest[[6]]
#)
#list of lists to be compared
#forven<-list(TY,O)

crosslinkven<-function(forven){
VENS<-list()
Crosslinkclasses<-list()
for(i in 1:length(forven)){
  VENS[[i]]<-ggvenn(
   forven[[i]], 
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
    stroke_size = 1, set_name_size = 3)
  Crosslinkclasses[[i]]<-get.venn.partitions(forven[[i]])
    #calculate.overlap(forven[[i]])
}
VENS<<-VENS
Crosslinkclasses<<-Crosslinkclasses
}

crosslinkven4<-function(forven){
  VENS<-list()
  Crosslinkclasses<-list()
  for(i in 1:length(forven)){
    VENS[[i]]<-ggvenn(
      forven[[i]], 
      fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF","#108686FF"),
      stroke_size = 1, set_name_size = 4)
    Crosslinkclasses[[i]]<-get.venn.partitions(forven[[i]])
    #calculate.overlap(forven[[i]])
  }
  VENS<<-VENS
  Crosslinkclasses<<-Crosslinkclasses
}

crosslinkven4<-function(forven){
  VENS<-list()
  Crosslinkclasses<-list()
  for(i in 1:length(forven)){
    VENS[[i]]<-ggvenn(
      forven[[i]], 
      fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF","#108686FF"),
      stroke_size = 1, set_name_size = 4)
    Crosslinkclasses[[i]]<-get.venn.partitions(forven[[i]])
    #calculate.overlap(forven[[i]])
  }
  VENS<<-VENS
  Crosslinkclasses<<-Crosslinkclasses
}
library(venn)
crosslinkven5<-function(forven){
  VENS<-list()
  Crosslinkclasses<-list()
  for(i in 1:length(forven)){
    VENS[[i]]<-venn(
      forven[[i]])
    Crosslinkclasses[[i]]<-get.venn.partitions(forven[[i]])
    #calculate.overlap(forven[[i]])
  }
  VENS<<-VENS
  Crosslinkclasses<<-Crosslinkclasses
}

        
#save data ####

savedata<-function(p,linktypes,Crosslinkclasses,VENS,rep){
  barplots[[rep]]<<-p
  linktypeslist[[rep]]<<-linktypes
  crosslinklist[[rep]]<<-Crosslinkclasses
  venslist[[rep]]<<-VENS
  #uniquelists[[rep]]<<-uniquelists
}


savedatares<-function(p,linktypes,Crosslinkclasses,VENS,rep){
  barplotsres[[rep]]<<-p
  linktypeslistres[[rep]]<<-linktypes
  crosslinklistres[[rep]]<<-Crosslinkclasses
  venslistres[[rep]]<<-VENS
  #uniquelists[[rep]]<<-uniquelists
}


#match conpep with crosslinks identified ####
conpepmatch<-function(meroxdata,listoflinks){
  linkspos<-data.frame()
  leftlist<-list()
  rightlist<-list()
  leftidlist<-list()
  rightidlist<-list()
  for(i in 1:length(listoflinks)){
    for(j in 1:length(meroxdata))
      for(k in 1:dim(meroxdata[[j]][1]))
        if(o[[i]] == meroxdata[[j]]$concpep[k]){
          if(meroxdata[[j]]$Peptide2[k] == "1"){
            leftlist[i]<-(meroxdata[[j]]$From[k]+as.numeric(gsub("K","",meroxdata[[j]]$best.linkage.position.peptide.1[k]))-1)
            rightlist[i]<-(meroxdata[[j]]$From[k]+as.numeric(gsub("K","",meroxdata[[j]]$best.linkage.position.peptide.2[k]))-1)
            leftidlist[i]<-meroxdata[[j]]$Protein.1[k]
            rightidlist[i]<-meroxdata[[j]]$Protein.2[k]
          }else{
          leftlist[i]<-(meroxdata[[j]]$From[k]+as.numeric(gsub("K","",meroxdata[[j]]$best.linkage.position.peptide.1[k]))-1)
          rightlist[i]<-(meroxdata[[j]]$From.1[k]+as.numeric(gsub("K","",meroxdata[[j]]$best.linkage.position.peptide.2[k]))-1)
          leftidlist[i]<-meroxdata[[j]]$Protein.1[k]
          rightidlist[i]<-meroxdata[[j]]$Protein.2[k]
        }}
  }
  
  linksleft<-ldply(leftlist, data.frame)
  linksright<-ldply(rightlist, data.frame)
  linksleftid<-ldply(leftidlist,data.frame)
  linksrightid<-ldply(rightidlist,data.frame)
  linkspos<-cbind(linksleftid,linksleft,linksrightid,linksright)
  colnames(linkspos)<-c("Protein.1","leftpos","Protein.2","rightpos")
  linkspos<<-linkspos
  
}




#Split crosslink list to two columns ####
splitlinks<-function(listoflinks){
  links<-list()
  linksright<-list()
  linksleft<-list()
  for(i in 1:length(listoflinks)){
    links[i]<-strsplit(listoflinks[i],"\\]")
    linksleft[i]<-gsub("\\[","",links[[i]][1])
    linksright[i]<-gsub("\\[","",links[[i]][2])
  
  }
  linksleft<-ldply(linksleft, data.frame)
  linksright<-ldply(linksright, data.frame)
  links<-cbind(linksleft,linksright)
  colnames(links)<-c("left","right")
  links<<-links
}

#check if proper split
propsplit<-function(linkframe){
  for(i in 1:dim(linkframe)[1]){
    if(is.na(linkframe$right[i])){
      link<-list()
      
      link[1]<-strsplit(linkframe$left[i],"}")
      linkframe$left[i]<-link[[1]][1]
      linkframe$right[i]<-link[[1]][2]
      
    }
  }
  links<<-linkframe
}



removerest<-function(linkframe){
for(i in 1:dim(linkframe)[1]){
  if(str_detect(links$left[i],"\\{")){
    linkframe$left[i]<-gsub("\\{","",linkframe$left[i])
  }
  if(str_detect(links$left[i],"\\}")){
    linkframe$left[i]<-gsub("\\}","",linkframe$left[i])
  }
  if(str_detect(links$right[i],"\\{")){
    linkframe$right[i]<-gsub("\\{","",linkframe$right[i])
  }
  if(str_detect(links$right[i],"\\}")){
    linkframe$right[i]<-gsub("\\}","",linkframe$right[i])
  }
}  
 links<<-linkframe 
}

#count occurences


counts<-function(datalist,countframe){
  countframe$counts<-0
  for(i in 1:dim(countframe)[1]){
    for(j in 1:length(datalist)){
      for(k in 1:dim(datalist[[j]])[1]){
        if(countframe$concpep[i] == datalist[[j]]$concpep[k]){
          countframe$counts[i]<-countframe$counts[i]+1
        }
        
      }
    }
    
  }
  forcount<<-countframe
}



