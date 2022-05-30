#Functions to for MeroxAdjust ####
#Written by: Victor Chrone


#Declare variables to be used 
interonek<-data.frame()
onek<-data.frame()
threek<-data.frame()
validated<-data.frame()

#Library
library(plyr)
library(readxl)
library(ggvenn)
library(tidyverse)
library(ggplot2)
library(stringr)
library(VennDiagram)
library(plyr)
library(gridExtra)

#Functions
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


#returns number of Ks in peptide
kcount1<-function(dat,i){
  counter = 0
  for(k in 1:length(str_split(dat$Peptide.1[i],"")[[1]])){
    if(str_split(dat$Peptide.1[i],"")[[1]][k]=="K"){
      counter = counter+1
    }
  }
  return(counter)
}
  
kcount2<-function(dat,i){
  counter = 0
  if(dat$Peptide2[i]=="1"){
    for(k in 1:length(str_split(dat$Peptide.1[i],"")[[1]])){
      if(str_split(dat$Peptide.1[i],"")[[1]][k]=="K"){
        counter = counter+1
      }
      }
    }else{
  for(k in 1:length(str_split(dat$Peptide2[i],"")[[1]])){
    if(str_split(dat$Peptide2[i],"")[[1]][k]=="K"){
      counter = counter+1
    }
  }
      }
  return(counter)
}

#Returns position of first internal K 
classone1<-function(dat,i){
    counterone<-0
    for(k in 1:length(str_split(dat$Peptide.1[i],"")[[1]])){
      counterone = counterone+1
      if(str_split(dat$Peptide.1[i],"")[[1]][k]=="K"){
        counterone = counterone-1
        break
      }
    }
  return(counterone)
}
    
classone2<-function(dat,i){
  counterone<-0
  for(k in 1:length(str_split(dat$Peptide2[i],"")[[1]])){
    counterone = counterone+1
    if(str_split(dat$Peptide2[i],"")[[1]][k]=="K"){
      counterone = counterone-1
      break
    }
  }
  return(counterone)
}


#Moves best linkage position to first internal K
kvalidation1<-function(dat,i){
  #if(str_count(dat$Peptide.1[i],"K")==1){
   # print(dat[i,])
    
  #}
  if(str_count(dat$Peptide.1[i])-2==(as.numeric(substring(dat$best.linkage.position.peptide.1[i],2)))){
    for(s in 2:(str_count(dat$Peptide.1[i])-1)){
      if(str_split(dat$Peptide.1[i],"")[[1]][s]=="K"){
        dat$best.linkage.position.peptide.1[i]<-paste("K",(s-1),sep="")
       break
      }
    }
  }
  return(dat[i,1:33])
}

#corrects peptide2
kvalidation2<-function(dat,i){
  if(dat$Peptide2[i] == "1"){
    if(str_count(dat$Peptide.1[i])-2==(as.numeric(substring(dat$best.linkage.position.peptide.2[i],2)))){
      for(s in 2:(str_count(dat$Peptide.1[i])-1)){
        if(str_split(dat$Peptide.1[i],"")[[1]][s]=="K"){
          dat$best.linkage.position.peptide.2[i]<-paste("K",(s-1),sep="")
          break
        }
      }
    }
  }else if(str_count(dat$Peptide2[i])-2==(as.numeric(substring(dat$best.linkage.position.peptide.2[i],2)))){
    for(s in 2:(str_count(dat$Peptide2[i])-1)){
      if(str_split(dat$Peptide2[i],"")[[1]][s]=="K"){
        dat$best.linkage.position.peptide.2[i]<-paste("K",(s-1),sep="")
        break
      }
    }
  }
  return(dat[i,1:33])
}
 

#Extract rows containing peptides with C-terminal K best linkage position
checkklist<-function(datlist){
  Klist<-list()
  for(i in 1:length(datlist)){
    CtermK<-data.frame()
    for(j in 1:dim(datlist[[i]])[1]){
      if((str_split(datlist[[i]]$Peptide.1[j],""))[[1]][str_count(datlist[[i]]$Peptide.1[j])-1] == "K"){
        if((str_count(datlist[[i]]$Peptide.1[j])-2==(as.numeric(substring(datlist[[i]]$best.linkage.position.peptide.1[i],2))))){
          CtermK[j,1:33]<-datlist[[i]][j,]
          datlist[[i]]$Peptide.1[j]<-NA
        }
      }else{if(datlist[[i]]$Peptide2[j] == "0"||datlist[[i]]$Peptide2[j] == "1"){next}
        else if((str_split(datlist[[i]]$Peptide2[j],""))[[1]][str_count(datlist[[i]]$Peptide2[j])-1] == "K"){
          if((str_count(datlist[[i]]$Peptide2[j])-2==(as.numeric(substring(datlist[[i]]$best.linkage.position.peptide.2[i],2))))){
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
  klist<<-Klist
    return(datlist)
  }
  


#Extract rows containing peptides with C-terminal K best linkage position

checkk<-function(datlist){
  Klist<-list()
  for(i in 1:length(datlist)){
    CtermK<-data.frame()
    for(j in 1:dim(datlist[[i]])[1]){
      if((str_split(datlist[[i]]$Peptide.1[j],""))[[1]][str_count(datlist[[i]]$Peptide.1[j])-1] == "K"){
        if((str_count(datlist[[i]]$Peptide.1[j])-2==(as.numeric(substring(datlist[[i]]$best.linkage.position.peptide.1[j],2))))){
          CtermK[j,1:33]<-datlist[[i]][j,]
          datlist[[i]]$Score[j]<-NA
          print("Peptide 1 wrong")
        }
      }
        if(datlist[[i]]$Peptide2[j] == "1"){
          if((str_split(datlist[[i]]$Peptide.1[j],""))[[1]][str_count(datlist[[i]]$Peptide.1[j])-1] == "K"){
            if((str_count(datlist[[i]]$Peptide.1[j])-2==((as.numeric(substring(datlist[[i]]$best.linkage.position.peptide.2[j],2)))))){
              CtermK[j,1:33]<-datlist[[i]][j,]
              datlist[[i]]$Score[j]<-NA
              print("inter peptide 2 worng")
            }
          }
        }else if(datlist[[i]]$Peptide2[j] != "0" ){
                if((str_split(datlist[[i]]$Peptide2[j],""))[[1]][str_count(datlist[[i]]$Peptide2[j])-1] == "K"){
                  if((str_count(datlist[[i]]$Peptide2[j])-2==((as.numeric(substring(datlist[[i]]$best.linkage.position.peptide.2[j],2)))))){
                    CtermK[j,1:33]<-datlist[[i]][j,]
                    datlist[[i]]$Score[j]<-NA
                    print("peptide 2 wrong")
                  }
                }
        }
    }
    CtermK<-CtermK[!is.na(CtermK$Score),]
    datlist[[i]]<-datlist[[i]][!is.na(datlist[[i]]$Score),]
    Klist[[i]]<-CtermK
  }
  Klist<<-Klist
  return(datlist)
}

#takes dataframe of cterm K best linkages and corrects them to internal K if the peptide only contains 2 K's
correct2k<-function(dat){
  for(i in 1:dim(dat)[1]){
    counter<-kcount1(dat,i)
    if(counter == 2){
      if(str_count(dat$Peptide.1[i])-2==((as.numeric(substring(dat$best.linkage.position.peptide.1[i],2))))){
        dat[i,1:33]<-kvalidation1(dat,i)
      }
    }
    counter<-kcount2(dat,i)
    if(counter == 2){
      if(dat$Peptide2[i]=="1"){
        if(str_count(dat$Peptide.1[i])-2==((as.numeric(substring(dat$best.linkage.position.peptide.2[i],2))))){
          dat[i,1:33]<-kvalidation2(dat,i)
        } 
      }else if(str_count(dat$Peptide2[i])-2==((as.numeric(substring(dat$best.linkage.position.peptide.2[i],2))))){
        dat[i,1:33]<-kvalidation2(dat,i)
      }
    }
  }
  return(dat)
}

#takes corrected dataframe of cterm Ks and sorts them based on the number of Ks and position
#correct2kdat
correct1k<-function(dat,interonek,onek,threek,validated){
  for(i in 1:dim(dat)[1]){
    counter<-kcount1(dat,i)
    if(counter == 1){
      #get position of k
      counter1<-classone1(dat,i)
      if(counter1 < (length(str_split(dat$Peptide.1[i],"")[[1]])-2)){
        if(kcount2(dat,i)>2){
          threek[i,1:33]<-dat[i,]
        }else if(kcount2(dat,i) == 1 && classone2(dat,i)== (length(str_split(dat$Peptide2[i],"")[[1]])-2)){
          onek[i,1:33]<-dat[i,]
        }
        else{
         interonek[i,1:33]<-dat[i,]
        }
      }else if(counter1 == (length(str_split(dat$Peptide.1[i],"")[[1]])-2)){
        onek[i,1:33]<-dat[i,]
      }
      dat[i,]<-NA
    }else if(counter == 2){
      validated[i,1:33]<-dat[i,]
      dat[i,]<-NA
    }else if(counter >=3){
            if(str_count(dat$Peptide.1[i])-2==((as.numeric(substring(dat$best.linkage.position.peptide.1[i],2))))){
              threek[i,1:33]<-dat[i,]
              dat[i,]<-NA
            }
    }
  }
  interonek<<-interonek
  onek<<-onek
  threek<<-threek
  validated<<-validated
  return(dat)
}

#resorts the new dataframes based on peptide 2

secondsort2k<-function(dat,interonek,onek,threek){
  #if(empty(dat)){}else{
  for(i in 1:dim(dat)[1]){
    if(is.na(dat$Peptide2[i])|dat$Peptide2[i] == 0){
    }else{
      if(dat$Peptide2[i]== 1){
        counter<-kcount2(dat,i)
        if(counter == 1){
          #get position of k
          counter1<-classone2(dat,i)
          if(counter1 < (length(str_split(dat$Peptide.1[i],"")[[1]])-2)){
            interonek[i,1:33]<-dat[i,]
            dat[i,]<-NA
          }else if(counter1 == (length(str_split(dat$Peptide.1[i],"")[[1]])-2)){
            onek[i,1:33]<-dat[i,]
            dat[i,]<-NA
          }
          
        }else if(counter >=3){
          if(str_count(dat$Peptide.1[i])-2==((as.numeric(substring(dat$best.linkage.position.peptide.2[i],2))))){
            threek[i,1:33]<-dat[i,]
            dat[i,]<-NA
          }
        }
      }else{
      counter<-kcount2(dat,i)
      if(counter == 1){
        #get position of k
        counter1<-classone2(dat,i)
        if(counter1 < (length(str_split(dat$Peptide2[i],"")[[1]])-2)){
          interonek[i,1:33]<-dat[i,]
          dat[i,]<-NA
        }else if(counter1 == (length(str_split(dat$Peptide2[i],"")[[1]])-2)){
          onek[i,1:33]<-dat[i,]
          dat[i,]<-NA
        }
        
      }else if(counter >=3){
        if(str_count(dat$Peptide2[i])-2==((as.numeric(substring(dat$best.linkage.position.peptide.2[i],2))))){
          threek[i,1:33]<-dat[i,]
          dat[i,]<-NA
        }
      }
      
      }
    }
  }
  #}
  interonek<<-interonek
  onek<<-onek
  threek<<-threek
  return(dat)
   
}


secondsort2krest<-function(dat,interonek,onek){
  for(i in 1:dim(dat)[1]){
    if(is.na(dat$Peptide2[i])|dat$Peptide2[i] == 0|dat$Peptide2[i] == 1){
      
    }else{
      counter<-kcount2(dat,i)
      if(counter == 1){
        #get position of k
        counter1<-classone2(dat,i)
        if(counter1 < (length(str_split(dat$Peptide2[i],"")[[1]])-2)){
          interonek[i,1:33]<-dat[i,]
          dat[i,]<-NA
        }else if(counter1 == (length(str_split(dat$Peptide2[i],"")[[1]])-2)){
          onek[i,1:33]<-dat[i,]
          dat[i,]<-NA
        }
      }
    }
  }
  interonek<<-interonek
  onek<<-onek
  return(dat)
}

mergedata<-function(datt,validatedlist,interoneklist){
  for(i in 1:length(datt)){
    datt[[i]]<-rbind(datt[[i]],validatedlist[[i]],interoneklist[[i]])
  }
  return(datt)
}





MeroxAdjust<-function(Merox_data){
  
  datlist<-Merox_data
  #Declare variables
  RTC<-list()
  Klist<-list()
  validatedlist<-list(NA,NA,NA)
  interoneklist<-list(NA,NA,NA)
  oneklist<-list(NA,NA,NA)
  threeklist<-list(NA,NA,NA)
  
  #Sort merox data 
  for(i in 1:length(datlist)){
    CtermK<-data.frame()
    for(j in 1:dim(datlist[[i]])[1]){
      if((str_split(datlist[[i]]$Peptide.1[j],""))[[1]][str_count(datlist[[i]]$Peptide.1[j])-1] == "K"){
        if((str_count(datlist[[i]]$Peptide.1[j])-2==(as.numeric(substring(datlist[[i]]$best.linkage.position.peptide.1[j],2))))){
          CtermK[j,1:33]<-datlist[[i]][j,]
          datlist[[i]]$Score[j]<-NA
          print("Peptide 1 wrong")
        }
      }
      if(datlist[[i]]$Peptide2[j] == "1"){
        if((str_split(datlist[[i]]$Peptide.1[j],""))[[1]][str_count(datlist[[i]]$Peptide.1[j])-1] == "K"){
          if((str_count(datlist[[i]]$Peptide.1[j])-2==((as.numeric(substring(datlist[[i]]$best.linkage.position.peptide.2[j],2)))))){
            CtermK[j,1:33]<-datlist[[i]][j,]
            datlist[[i]]$Score[j]<-NA
            print("inter peptide 2 worng")
          }
        }
      }else if(datlist[[i]]$Peptide2[j] != "0" ){
        if((str_split(datlist[[i]]$Peptide2[j],""))[[1]][str_count(datlist[[i]]$Peptide2[j])-1] == "K"){
          if((str_count(datlist[[i]]$Peptide2[j])-2==((as.numeric(substring(datlist[[i]]$best.linkage.position.peptide.2[j],2)))))){
            CtermK[j,1:33]<-datlist[[i]][j,]
            datlist[[i]]$Score[j]<-NA
            print("peptide 2 wrong")
          }
        }
      }
    }
    CtermK<-CtermK[!is.na(CtermK$Score),]
    RTC[[i]]<-datlist[[i]][!is.na(datlist[[i]]$Score),]
    Klist[[i]]<-CtermK
  }

#Sort and fix K  
  for(p in 1:length(Klist[1:length(Klist)])){
    interonek<-data.frame()
    onek<-data.frame()
    threek<-data.frame()
    validated<-data.frame()
    
   #Fix k in first colum if kcounter =2 
    for(i in 1:dim(Klist[[p]])[1]){
      counter<-kcount1(Klist[[p]],i)
      if(counter == 2){
        if(str_count(Klist[[p]]$Peptide.1[i])-2==((as.numeric(substring(Klist[[p]]$best.linkage.position.peptide.1[i],2))))){
          Klist[[p]][i,1:33]<-kvalidation1(Klist[[p]],i)
        }
      }
      counter<-kcount2(Klist[[p]],i)
      if(counter == 2){
        if(Klist[[p]]$Peptide2[i]=="1"){
          if(str_count(Klist[[p]]$Peptide.1[i])-2==((as.numeric(substring(Klist[[p]]$best.linkage.position.peptide.2[i],2))))){
            Klist[[p]][i,1:33]<-kvalidation2(Klist[[p]],i)
          } 
        }else if(str_count(Klist[[p]]$Peptide2[i])-2==((as.numeric(substring(Klist[[p]]$best.linkage.position.peptide.2[i],2))))){
          Klist[[p]][i,1:33]<-kvalidation2(Klist[[p]],i)
        }
      }
    }
    
    for(i in 1:dim(Klist[[p]])[1]){
      counter<-kcount1(Klist[[p]],i)
      if(counter == 1){
        #get position of k
        counter1<-classone1(Klist[[p]],i)
        if(counter1 < (length(str_split(Klist[[p]]$Peptide.1[i],"")[[1]])-2)){
          if(kcount2(Klist[[p]],i)>2){
            threek[i,1:33]<-Klist[[p]][i,]
          }else if(kcount2(Klist[[p]],i) == 1 && classone2(Klist[[p]],i)== (length(str_split(Klist[[p]]$Peptide2[i],"")[[1]])-2)){
            onek[i,1:33]<-Klist[[p]][i,]
          }
          else{
            interonek[i,1:33]<-Klist[[p]][i,]
          }
        }else if(counter1 == (length(str_split(Klist[[p]]$Peptide.1[i],"")[[1]])-2)){
          onek[i,1:33]<-Klist[[p]][i,]
        }
        Klist[[p]][i,]<-NA
      }else if(counter == 2){
        validated[i,1:33]<-Klist[[p]][i,]
        Klist[[p]][i,]<-NA
      }else if(counter >=3){
        if(str_count(Klist[[p]]$Peptide.1[i])-2==((as.numeric(substring(Klist[[p]]$best.linkage.position.peptide.1[i],2))))){
          threek[i,1:33]<-Klist[[p]][i,]
          Klist[[p]][i,]<-NA
        }
      }
    }
    
    
    #validated<-secondsort2k(dat = validated,interonek = interonek,onek = onek,threek = threek)
    
    for(i in 1:dim(validated)[1]){
      if(is.na(validated$Peptide2[i])|validated$Peptide2[i] == 0|validated$Peptide2[i] == 1){
        
      }else{
        counter<-kcount2(validated,i)
        if(counter == 1){
          #get position of k
          counter1<-classone2(validated,i)
          if(counter1 < (length(str_split(validated$Peptide2[i],"")[[1]])-2)){
            interonek[i,1:33]<-validated[i,]
            validated[i,]<-NA
          }else if(counter1 == (length(str_split(validated$Peptide2[i],"")[[1]])-2)){
            onek[i,1:33]<-validated[i,]
            validated[i,]<-NA
          }
        }
      }
    }
  
    
    #Klist[[i]]<-secondsort2krest(dat = Klist[[i]],interonek = interonek,onek = onek)
    for(i in 1:dim(Klist[[p]])[1]){
      if(is.na(Klist[[p]]$Peptide2[i])|Klist[[p]]$Peptide2[i] == 0|Klist[[p]]$Peptide2[i] == 1){
      }else{
        counter<-kcount2(Klist[[p]],i)
        if(counter == 1){
          #get position of k
          counter1<-classone2(Klist[[p]],i)
          if(counter1 < (length(str_split(Klist[[p]]$Peptide2[i],"")[[1]])-2)){
            interonek[i,1:33]<-Klist[[p]][i,]
            Klist[[p]][i,]<-NA
          }else if(counter1 == (length(str_split(Klist[[p]]$Peptide2[i],"")[[1]])-2)){
            onek[i,1:33]<-Klist[[p]][i,]
            Klist[[p]][i,]<-NA
          }
        }
      }
    }
    
    if(!empty(onek)){
      oneklist[[p]]<-onek %>% drop_na(Score)
    }else{oneklist[[p]]<-na.omit(onek)}
    if(!empty(interonek)){
      interoneklist[[p]]<-interonek %>% drop_na(Score)
    }else{interoneklist[[p]]<-na.omit(interonek)}
    if(!empty(threek)){
      threeklist[[p]]<-threek %>% drop_na(Score)
    }else{threeklist[[p]]<-na.omit(threek)}
    if(!empty(validated)){
      validatedlist[[p]]<-validated %>% drop_na(Score)
    }else{validated[[p]]<-na.omit(validated)}
    Klist[[p]]<-na.omit(Klist[[p]])
  }
  #Merge data
  for(i in 1:length(RTC)){
    RTC[[i]]<-rbind(RTC[[i]],validatedlist[[i]],interoneklist[[i]])
  }
  
  return(RTC)
}



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










