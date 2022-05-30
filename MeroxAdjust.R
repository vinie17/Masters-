library(shiny)
library(tidyr)
library(plyr)
library(writexl)
library(stringr)
#library(rstudioapi)
library(DT)
library(dplyr)
library(shinyWidgets)
library(bslib)
library(seqinr)
#setwd("C:\\Users\\victo\\Desktop\\speciale\\R\\MEROXSORT\\MeroxAdjust")
source("sort_merox_functions_v2.R")
#setwd("C:\\Users\\victo\\Desktop\\speciale\\R")
source("meroxfunctions2.R")

MeroxAdjust1<-function(Merox_data){
  
  datlist<-Merox_data
  
  #check which data set we have
  
  if(colnames(datlist)[1]== "PepSeq1"){
    
    colnames(datlist)<- c("Peptide.1","Peptide2","From","From.1","best.linkage.position.peptide.1","best.linkage.position.peptide.2","Protein.1","Protein.2",colnames(datlist)[9:25])
    datlist$Peptide.1<-str_c("[",datlist$Peptide.1,"]",sep = "")
    datlist$Peptide2<-str_c("[",datlist$Peptide2,"]",sep = "")
    datlist$best.linkage.position.peptide.1<-str_c("K",datlist$best.linkage.position.peptide.1,sep="")
    datlist$best.linkage.position.peptide.2<-str_c("K",datlist$best.linkage.position.peptide.2,sep="")
    datlist$Protein.1<-str_c(">",datlist$Protein.1,sep="")
    datlist$Protein.2<-str_c(">",datlist$Protein.2,sep="")
    
  }
  
  
  
  
  #Sort merox data 
  CtermK<-data.frame()
  for(j in 1:dim(datlist)[1]){
    if((str_split(datlist$Peptide.1[j],""))[[1]][str_count(datlist$Peptide.1[j])-1] == "K"){
      if((str_count(datlist$Peptide.1[j])-2==(as.numeric(substring(datlist$best.linkage.position.peptide.1[j],2))))){
        CtermK[j,1:33]<-datlist[j,]
        datlist$Score[j]<-NA
        #print("Peptide 1 wrong")
      }
    }
    if(datlist$Peptide2[j] == "1"){
      if((str_split(datlist$Peptide.1[j],""))[[1]][str_count(datlist$Peptide.1[j])-1] == "K"){
        if((str_count(datlist$Peptide.1[j])-2==((as.numeric(substring(datlist$best.linkage.position.peptide.2[j],2)))))){
          CtermK[j,1:33]<-datlist[j,]
          datlist$Score[j]<-NA
          #print("inter peptide 2 worng")
        }
      }
    }else if(datlist$Peptide2[j] != "0" ){
      if((str_split(datlist$Peptide2[j],""))[[1]][str_count(datlist$Peptide2[j])-1] == "K"){
        if((str_count(datlist$Peptide2[j])-2==((as.numeric(substring(datlist$best.linkage.position.peptide.2[j],2)))))){
          CtermK[j,1:33]<-datlist[j,]
          datlist$Score[j]<-NA
          #print("peptide 2 wrong")
        }
      }
    }
  }
  
  CtermK<-CtermK[!is.na(CtermK$Score),]
  RTC<-datlist[!is.na(datlist$Score),]
  RTC$corrected<-"FALSE"
  Klist<-CtermK
  
  #Sort and fix K  
  interonek<-data.frame()
  onek<-data.frame()
  threek<-data.frame()
  validated<-data.frame()
  
  #Fix k in first colum if kcounter =2 
  for(i in 1:dim(Klist)[1]){
    counter<-kcount1(Klist,i)
    if(counter == 2){
      if(str_count(Klist$Peptide.1[i])-2==((as.numeric(substring(Klist$best.linkage.position.peptide.1[i],2))))){
        Klist[i,1:33]<-kvalidation1(Klist,i)
      }
    }
    counter<-kcount2(Klist,i)
    if(counter == 2){
      if(Klist$Peptide2[i]=="1"){
        if(str_count(Klist$Peptide.1[i])-2==((as.numeric(substring(Klist$best.linkage.position.peptide.2[i],2))))){
          Klist[i,1:33]<-kvalidation2(Klist,i)
        } 
      }else if(str_count(Klist$Peptide2[i])-2==((as.numeric(substring(Klist$best.linkage.position.peptide.2[i],2))))){
        Klist[i,1:33]<-kvalidation2(Klist,i)
      }
    }
  }
  
  for(i in 1:dim(Klist)[1]){
    counter<-kcount1(Klist,i)
    if(counter == 1){
      #get position of k
      counter1<-classone1(Klist,i)
      if(counter1 < (length(str_split(Klist$Peptide.1[i],"")[[1]])-2)){
        if(kcount2(Klist,i)>2){
          threek[i,1:33]<-Klist[i,]
        }else if(kcount2(Klist,i) == 1 && classone2(Klist,i)== (length(str_split(Klist$Peptide2[i],"")[[1]])-2)){
          onek[i,1:33]<-Klist[i,]
        }
        else{
          interonek[i,1:33]<-Klist[i,]
        }
      }else if(counter1 == (length(str_split(Klist$Peptide.1[i],"")[[1]])-2)){
        onek[i,1:33]<-Klist[i,]
      }
      Klist[i,]<-NA
    }else if(counter == 2){
      validated[i,1:33]<-Klist[i,]
      Klist[i,]<-NA
    }else if(counter >=3){
      if(str_count(Klist$Peptide.1[i])-2==((as.numeric(substring(Klist$best.linkage.position.peptide.1[i],2))))){
        threek[i,1:33]<-Klist[i,]
        Klist[i,]<-NA
      }
    }
  }
  
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
  
  if(is.na(Klist$Peptide2[i])|Klist$Peptide2[i] == 0|Klist$Peptide2[i] == 1){
  }else{
    counter<-kcount2(Klist,i)
    if(counter == 1){
      #get position of k
      counter1<-classone2(Klist,i)
      if(counter1 < (length(str_split(Klist$Peptide2[i],"")[[1]])-2)){
        interonek[i,1:33]<-Klist[i,]
        Klist[i,]<-NA
      }else if(counter1 == (length(str_split(Klist$Peptide2[i],"")[[1]])-2)){
        onek[i,1:33]<-Klist[i,]
        Klist[i,]<-NA
      }
    }
  }
  
  if(!empty(onek)){
    onek<-onek %>% drop_na(Peptide.1)
  }else{onek<-na.omit(onek)}
  if(!empty(interonek)){
    interonek<-interonek %>% drop_na(Peptide.1)
  }else{interonek<-na.omit(interonek)}
  if(!empty(threek)){
    threek<-threek %>% drop_na(Peptide.1)
  }else{threek<-na.omit(threek)}
  if(!empty(validated)){
    validated<-validated %>% drop_na(Peptide.1)
  }else{validated<-na.omit(validated)}
  Klist<-na.omit(Klist)
  
  
  
  
  #Merge data
  
  
  
  if(dim(RTC)[2]==26){
    #print(validated)
    validated<-validated[,1:25]
    interonek<-interonek[,1:25]
    
    validated$corrected<-"TRUE"
    interonek$corrected<-"TRUE"
    
    RTC<-rbind(RTC,validated,interonek)
  }else{
    RTC<-rbind(RTC,validated,interonek)
    validated$corrected<-"TRUE"
    interonek$corrected<-"TRUE"
  }
  threek<<-threek
  onek<<-onek
  
  #for plot
  if(dim(RTC)[2]==26){
    correctedf1<-list(RTC)
    
    conpepres2(correctedf1)
    for(i in 1:dim(listdata[[1]])[1]){
      if(listdata[[1]]$Peptide2[i]== 1){
        if(str_split(listdata[[1]]$conres[i],"\\|")[[1]][1] == str_split(listdata[[1]]$conres[i],"\\|")[[1]][2]){
          #print(listdata[[1]]$conres[i])
          #print(RTC$Peptide2[i])
          RTC[i,]<-NA
          #print("removed")
          
        }
      }
    }
  }else{
    correctedf<-list(RTC)
    conpepres2(correctedf)
    
    #for(i in 1:dim(listdata)[1]){
    # if(listdata$Peptide2[i]== 1){
    #  if(str_split(listdata[i],"\\|")[[1]][1] == str_split(dataone[i],"\\|")[[1]][2]){
    #   RTC[i,]<-NA
    #  }
    #}
    #}
    
    linklist<-c("dat")
    linktypelistzero<-list(NA)
    linktypelistone<-list(NA)
    linktypelistrest<-list(NA)
    linktypelistzerores<-list(NA)
    linktypelistoneres<-list(NA)
    linktypelistrestres<-list(NA)
    linktypes<-data.frame()
    
    for(i in 1:length(listdata)){
      linktypelistzero[[i]]<-filter(listdata[[i]], Peptide2 ==0)$conresid
      linktypelistone[[i]]<-filter(listdata[[i]], Peptide2 ==1)$conresid
      holdrest<-filter(listdata[[i]], Peptide2 !=1)
      linktypelistrest[[i]]<-filter(holdrest, Peptide2 !=0)$conresid
      
      linktypes[i,1]<-length(linktypelistzero[[i]])
      linktypes[i,2]<-length(linktypelistone[[i]])
      linktypes[i,3]<-length(linktypelistrest[[i]])
      
    }
    
    #colnames for linktypes 
    colnames(linktypes)<-c("zero","one","rest")
    countunique<-unique(RTC$Protein.1)
    
    
    
    barlist<-list()
    for(i in 1:length(countunique)){
      barlist[[countunique[i]]]<-filterdata(linktypelistzero[[1]],
                                            linktypelistrest[[1]],
                                            linktypelistone[[1]],
                                            countunique[i])
    }
    oneandrest<-list()
    oneandrest<<-mapply(c,linktypelistone, linktypelistrest, SIMPLIFY=FALSE)
    linktypes<<-linktypes
    countunique<<-countunique
    barlist<<-barlist
    
  }
  if(!empty(RTC)){
    RTC<-RTC %>% drop_na(Peptide.1)
  }else{RTC<-na.omit(RTC)}
  
  
  #final check of wrongly corrected internal oneks
  for(i in 1:dim(RTC)[1]){
    
    if(RTC$best.linkage.position.peptide.1[i] == "K1"){
      if(str_count(RTC$Peptide.1[i]==1)){
        print("K1 detected")
        counterone<-1
        while(str_split(RTC$Peptide.1[i],"")[[1]][counterone]!="K"){
          counterone<-counterone+1
          print("K at position: ")
          print((counterone+1))
        }
        if(str_split(RTC$Peptide.1[i],"")[[1]][counterone+1]=="K"){
          RTC$best.linkage.position.peptide.1[i] = str_c("K",(counterone+1),"")
          print("corrected")
          print(RTC$best.linkage.position.peptide.1[i])
          print(RTC$Peptide.1[i])
        }
        
      }
    }
  }
  
  
  return(RTC)
}
filterdata<-function(zero,one,rest,chain){
  ocreAzero<-as.data.frame(na.omit(zero))
  ocreAone<-as.data.frame(na.omit(rest))
  ocreArest<-as.data.frame(na.omit(one))
  
  #for zero
  for(i in 1:dim(ocreAzero)[1]){
    ocreAzero$chainA[i]<-str_split(ocreAzero[i,1],"\\|")[[1]][3]
    ocreAzero$resid1[i]<-str_split(ocreAzero[i,1],"\\|")[[1]][1]
    ocreAzero$color[i]<-"zero"
  }
  #for one
  for(i in 1:dim(ocreAone)[1]){
    ocreAone$chainA[i]<-str_split(ocreAone[i,1],"\\|")[[1]][3]
    ocreAone$chainB[i]<-str_split(ocreAone[i,1],"\\|")[[1]][4]
    ocreAone$redid1[i]<-str_split(ocreAone[i,1],"\\|")[[1]][1]
    ocreAone$resid2[i]<-str_split(ocreAone[i,1],"\\|")[[1]][2]
    ocreAone$color[i]<-"one"
  }
  #for rest
  for(i in 1:dim(ocreArest)[1]){
    ocreArest$chainA[i]<-str_split(ocreArest[i,1],"\\|")[[1]][3]
    ocreArest$chainB[i]<-str_split(ocreArest[i,1],"\\|")[[1]][4]
    ocreArest$redid1[i]<-str_split(ocreArest[i,1],"\\|")[[1]][1]
    ocreArest$resid2[i]<-str_split(ocreArest[i,1],"\\|")[[1]][2]
    ocreArest$color[i]<-"rest"
  }
  #filter and create list
  ocreArest1<-cbind(ocreArest$chainA,ocreArest$redid1,ocreArest$color)
  ocreArest2<-cbind(ocreArest$chainB,ocreArest$resid2,ocreArest$color)
  ocreArest<-rbind(ocreArest1,ocreArest2)
  
  ocreAone1<-cbind(ocreAone$chainA,ocreAone$redid1,ocreAone$color)
  ocreAone2<-cbind(ocreAone$chainB,ocreAone$resid2,ocreAone$color)
  ocreAone<-rbind(ocreAone1,ocreAone2)
  
  ocreAzero<-cbind(ocreAzero$chainA,ocreAzero$resid1,ocreAzero$color)
  ocreAtot<-as.data.frame(rbind(ocreArest,ocreAone,ocreAzero))
  colnames(ocreAtot)<-c("chainA","resid","color")
  ocreAlight<-data.frame()
  ocreAlight<-filter(ocreAtot, chainA == chain)
  
  ocreAlight1<-as.data.frame(ocreAlight[order(ocreAlight$resid),])
  
  return(ocreAlight1)
}
lightchainplot<-function(dat,title,yaxisval){
  a<-ggplot(data = dat,aes(x=factor(resid), fill = color))+
    geom_bar(position = "stack")+labs(title=title,x= "Residue number")+
    theme(axis.text.x = element_text(angle = 90,size = 5, face = "bold"))+
    scale_y_continuous(limits=c(0,yaxisval))+
    scale_color_manual(values = c("one"="grey20",
                                  "rest"="grey40",
                                  "zero"="grey60",
                                  "No hits"="grey100"))+
    scale_fill_manual(values = c("one"="grey20",
                                 "rest"="grey40",
                                 "zero"="grey60",
                                 "No hits"="grey100"))
  return(a)
}
scatplotSLICE1<-function(dat, jitx,jity,lab,fasta,protein1,protein2){
  TOP<-dat
  TOPbefore <- as.data.frame(TOP[1])
  TOP<-list(strsplit(as.character(TOPbefore[,1]),'\\|'))
  TOP<-as.data.frame(do.call(rbind, TOP[[1]]))
  TOP$conres<-str_c(TOP$V1,TOP$V3,TOP$V2,TOP$V4,sep = "|")
  TOP$conresrev<-str_c(TOP$V2,TOP$V4,TOP$V1,TOP$V3,sep = "|")
  TOP<-checkrev2(TOP)
  TOP<-list(strsplit(as.character(TOP$conres),'\\|'))
  TOP<-as.data.frame(do.call(rbind, TOP[[1]]))
  TOP$V1<-as.numeric(TOP$V1)
  TOP$V3<-as.numeric(TOP$V3)
  
  #fasta<-TEST
  #if(as.double(length(fasta[[1]]))>as.double(length(fasta[[2]]))){
  # heavych<-as.double(length(fasta[[1]]))
  #lightch<-as.double(length(fasta[[2]]))
  #light<-attr(fasta[[2]],"name")
  #}else{
  # heavyhc<-as.double(length(fasta[[2]]))
  #lightch<-as.double(length(fasta[[1]]))
  #light<-attr(fasta[[1]],"name")
  #}
  
  heavych<-as.double(length(fasta[str_split(protein1,">")[[1]][2]][[1]]))
  lightch<-as.double(length(fasta[str_split(protein2,">")[[1]][2]][[1]]))
  
  
  
  #lightchC<-c(strsplit(lightChainC,""))
  #heavych<-c(strsplit(heavyChain,""))
  #heavychFc<-c(strsplit(heavyChainFc,""))
  #lightch<-as.double(length(lightch[[1]]))
  #lightchC<-as.double(length(lightchC[[1]]))
  #heavych<-as.double(length(heavych[[1]]))
  #heavychFc<-as.double(length(heavychFc[[1]]))
  TOP<-na.omit(TOP)
  
  if(protein1 != protein2){
    for(i in 1:dim(TOP)[1]){
      if(TOP$V2[i] == protein2){
        TOP$V1[i] = TOP$V1[i]+(heavych+1) 
      }
      if(TOP$V4[i] == protein2){
        TOP$V3[i] = TOP$V3[i]+(heavych+1) 
      }
    }
  }
  for(i in 1:dim(TOP)[1]){
    if(TOP$V2[i] == TOP$V4[i]){
      if(TOP$V1[i]==TOP$V3[i]){
        TOP$col[i]= "Homeotypic"
        TOP$V1[i]= TOP$V1[i]+0.1
      }
      else{
        TOP$col[i] = "Inter"
      }
    }
    else{TOP$col[i]="Intra"}
    
  }
  
  for(i in 1:dim(TOP)[1]){
    if(TOP$V3[i]<TOP$V1[i]){
      one<-TOP$V1[i]
      TOP$V1[i]<-TOP$V3[i]
      TOP$V3[i]<-one
      
    }
  }
  splitdat<<-TOP
  jitter<-position_jitter(width = jitx,height=jity)
  
  #scatter plot os distaces
  TOP$col<-as.factor(TOP$col)
  if(protein1 != protein2){
    ppp<-ggplot(TOP, aes(x=V1, y=V3, group=col))+
      geom_point(aes(colour=col),position =jitter, alpha=0.2,size=2)+
      scale_size(range = c(1,11 ), name="Distance (?)") +
      ylab("HC                                                                       HV                                                        LC                                          LV")+
      xlab("HC                                                                       HV                                                        LC                                          LV")+
      ggtitle(lab)+
      scale_color_manual(values = setNames(c('#619CFF',"#00BA38","#F8766D"),c("Intra","Inter","Homeotypic")))+
      
      scale_x_continuous(expand = c(0,0), breaks = seq(0,(heavych+lightch+1), by =30), limits = c(0,(1+lightch+heavych+10)))+
      scale_y_continuous(expand = c(0,0), breaks= seq(0, (heavych+lightch+1), by =30), limits = c(0,(1+lightch+heavych+10)))+
      theme(axis.text.x = element_text(angle = 90, face = "bold",size = 7),
            axis.text.y = element_text(face = "bold", size = 7),
            #axis.title.y = element_text(size = 7, hjust = 1), #for big plot
            #axis.title.x = element_text(hjust = 1),
            axis.title.y = element_text(size = 6, hjust = 1, face = "bold"), #for grouped
            axis.title.x = element_text(hjust = 1, size = 7, face = "bold"),
            axis.text = element_text(
              color =  c(
                #rep("grey1",length(seq(1,(heavychFc), by =30))),
                rep("grey40",length(seq(0,(heavych), by =30))),
                #rep("grey15",length(seq(heavych,(lightchC+heavych+1), by =30))),
                rep("grey55",length(seq((heavych+1),(lightch+heavych+1), by =30))))))
  }else{
    ppp<-ggplot(TOP, aes(x=V1, y=V3, group=col))+
      geom_point(aes(colour=col),position =jitter, alpha=0.2,size=2)+
      scale_size(range = c(1,11 ), name="Distance (?)") +
      ylab("beta peptide")+
      xlab("alpha peptide")+
      ggtitle(lab)+
      scale_color_manual(values = setNames(c('#619CFF',"#00BA38","#F8766D"),c("Intra","Inter","Homeotypic")))+
      annotate(geom="rect",xmin = 120, xmax=451,
               ymin = 541, ymax=661,
               color = NA,
               fill = "grey60",alpha = 0.3)+
      scale_x_continuous(expand = c(0,0), breaks = seq(0,(heavych+1), by =30), limits = c(0,(1+heavych+10)))+
      scale_y_continuous(expand = c(0,0), breaks= seq(0, (heavych+1), by =30), limits = c(0,(1+heavych+10)))+
      theme(axis.text.x = element_text(angle = 90, face = "bold",size = 7),
            axis.text.y = element_text(face = "bold", size = 7),
            #axis.title.y = element_text(size = 7, hjust = 1), #for big plot
            #axis.title.x = element_text(hjust = 1),
            axis.title.y = element_text(size = 6, hjust = 1, face = "bold"), #for grouped
            axis.title.x = element_text(hjust = 1, size = 7, face = "bold"),
            axis.text = element_text(
              color =  c(
                #rep("grey1",length(seq(1,(heavychFc), by =30))),
                rep("grey40",length(seq(0,(heavych), by =30))))))
    
  }
  return(ppp)
}

homeocheck<-function(correcteddata){
  
  conpepdata<-list(correcteddata)
  conpepres2(conpepdata)
  dataset<-listdata
  
  linklist<-c("Data")
  
  linktypelistzero<-list(NA)
  linktypelistone<-list(NA)
  linktypelistrest<-list(NA)
  
  linktypelistzerores<-list(NA)
  linktypelistoneres<-list(NA)
  linktypelistrestres<-list(NA)
  
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
  
  linktypelistoneres<-linktypelistone
  linktypelistrestres<-linktypelistrest
  linktypelistzerores<-linktypelistzero
  linktypesres<-linktypes
  
  
  #filter linktypelistrestres filtering 
  linktypelisthomeotypic<-list(NA)
  datares<-linktypelistrestres
  homeotypiclist<-list(NA)
  
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
  
  linktypelistrestres<-datares
  linktypelisthomeotypic<-linktypelisthomeotypic
  
  for(i in 1:length(linktypelistrestres)){
    linktypelistrestres[[i]]<-na.omit(linktypelistrestres[[i]])
  }
  for(i in 1:length(linktypelisthomeotypic)){
    linktypelisthomeotypic[[i]]<-na.omit(linktypelisthomeotypic[[i]])
  }
  
  
  #filter linktypelistoneres filtering 
  linktypelisthomeotypicone<-list(NA)
  dataone1<-linktypelistoneres
  homeotypiclistone<-list(NA)
  
  for(i in 1:length(dataone1)){
    for(j in 1:length(dataone1[[i]])){
      if(str_split(dataone1[[i]][j],"\\|")[[1]][1] == str_split(dataone1[[i]][j],"\\|")[[1]][2]){
        linktypelisthomeotypicone[[i]][j]<-dataone1[[i]][j]
        dataone1[[i]][j]<-NA
      }else{
        linktypelisthomeotypicone[[i]][j]<-NA
      }
    }
  }
  linktypelisthomeotypicone<-linktypelisthomeotypicone
  linktypelistoneres<-dataone1
  
  
  for(i in 1:length(linktypelistoneres)){
    linktypelistoneres[[i]]<-na.omit(linktypelistoneres[[i]])
  }
  for(i in 1:length(linktypelisthomeotypicone)){
    linktypelisthomeotypicone[[i]]<-na.omit(linktypelisthomeotypicone[[i]])
  }
  
  
  #filter the homeotypic peptides from original data
  #can be used for either rest(homeotypic) and one(to find S and Y links not filtered by meroxsort)
  if(is_empty(linktypelisthomeotypic[[1]])==FALSE){
    testlist<-list()
    if(is_empty(linktypelisthomeotypic[[1]])){}else{
      test<-data.frame()
      test2<-data.frame()
      uniqueones<-unique(linktypelisthomeotypic[[1]])
      #for(i in 1:length(uniqueones)){
      test<-filter(dataset[[1]],conresid %in% uniqueones)
      test2<-rbind(test2,test)
      
      testlist[[1]]<-test2
      homeotypiclist<-testlist
    }
    homeotypiclist[[1]]<-filter(homeotypiclist[[1]], Peptide2 != 1)
    
    linktypelisthomeotypic[[1]]<-homeotypiclist[[1]]$conresid
  }
  
  #un each time to remove false homeotypic (S and Y)
  if(is_empty(linktypelisthomeotypicone[[1]])==FALSE){
    findhomeoone(dataset,linktypelisthomeotypicone[[1]])
    testlist<-list()
    if(is_empty(linktypelisthomeotypicone[[1]])){}else{
      test<-data.frame()
      test2<-data.frame()
      uniqueones<-unique(linktypelisthomeotypicone[[1]])
      
      test<-filter(dataset[[1]],conresid %in% uniqueones)
      test2<-rbind(test2,test)
      
      testlist[[1]]<-test2
      homeotypiclistone<-testlist
    }
    homeotypiclistone[[1]]<-filter(homeotypiclistone[[1]], Peptide2 == 1)
    linktypelisthomeotypicone[[1]]<-homeotypiclistone[[1]]$conresid
  }
  `%!in%` <- Negate(`%in%`)
  if(is_empty(linktypelisthomeotypic[[1]])==FALSE){
    cleaned1<-filter(dataset[[1]],conresid %!in% linktypelisthomeotypic[[1]])
  }else{cleaned1<-dataset[[1]]}
  if(is_empty(linktypelisthomeotypicone[[1]])== FALSE){
    cleaned2<-filter(cleaned1,conresid %!in% linktypelisthomeotypicone[[1]])
  }else{cleaned2<-cleaned1}
  
  
  #cleaned3<-filter(cleaned2, conresid != homeotypiclist[[1]]$conresid)
  #cleaned4<-filter(cleaned3, conresid != homeotypiclistone[[1]]$conresid)
  
  
  
  correcteddfneg<<-cleaned2
  homeotypicdata<<-homeotypiclist[[1]] #[!duplicated(homeotypiclist[[1]][,c("conresid")]),]
  serinethreonine<<-homeotypiclistone[[1]] #[!duplicated(homeotypiclistone[[1]][,c("conresid")]),]
  
  
}

ui <- navbarPage(
  #theme = bs_theme(version = 4, bootswatch ="slate"),
  title = "Merox Adjust", 
  main_page <- tabPanel(
    title = "Adjust",
    sidebarLayout(
      sidebarPanel(
        title = "Inputs",
        fileInput("filedata","Select Merox CSV File to Import",accept=".csv"),
        checkboxInput("checkbox",label = "xiVIEW data",value = FALSE),
        fileInput("fasta","Select FastA file to import", accept = ".fastA"),
        actionButton("submit","Submit data"),
        downloadButton("downloadData","Download Adjusted data"),
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            title = "Adjusted data",
            textOutput('zero'),
            textOutput('one'),
            textOutput('rest'),
            textOutput('total'),
            textOutput('sandt'),
            textOutput('homeo2'),
            tabsetPanel(
              tabPanel("Validated", DT::dataTableOutput('corrected')),
              tabPanel("Serine/Threonine linked", DT::dataTableOutput('onek')),
              tabPanel("Needs further validation", DT::dataTableOutput('threek')),
              tabPanel("Validated Homeotypic(-)", DT::dataTableOutput('corrected2')),
              tabPanel("False homeotypic", DT::dataTableOutput('serine')),
              tabPanel("Homeotypic", DT::dataTableOutput('homeo'))
              
            )
          ),
          tabPanel(
            title = "Plots",
            tabsetPanel(
              tabPanel("Bar chart",
                       sidebarPanel(
                         title = "Inputs",
                         selectInput("Chain",
                                     label = h5("Select Protein/chain"),
                                     ""),
                         sliderInput("yaxis","Y-axis max:",
                                     min = 50, max= 500, value = 200, step = 5,
                                     animate = animationOptions(100)
                         ),
                       ),
                       mainPanel(
                         plotOutput("plot_1"),
                         #plotOutput("plot_2")
                       ),
              ),
              tabPanel("Scatterplot",
                       sidebarPanel(
                         title = "Choose chains/proteins",
                         selectInput("protein1",
                                     label = h5("Select Protein/chain 1"),
                                     ""),
                         selectInput("protein2",
                                     label = h5("Select Protein/chain 2"),
                                     "")),
                       plotOutput("plot_3")
              )
            )
          )
        )
      )
    )
  ),
  about_page <- tabPanel(title = "About",
                         titlePanel("About"),
                         "Merox Adjust is a tool developed to correct output from Merox. Merox systematically
                            wrongly assigns cross-links to the C-termianl Lysine in peptides eventhough this 
                            is not possible if working with trypic peptides as cross-linking with an NHS-ester,
                            reacts with the Lysine and removes the charge which is essential for the enzymatic
                            process of trypsin, another factor is that trypsin no longer can access the lysine 
                            due to the cross-linker sterically being in the way. This tool should however not 
                            be used if you have cross-linked after digestion with trypsin. Or if you havent 
                            quenched the cross-linker.",
                         br(),
                         "",
                         br(),
                         "If the peptide is C-terminally linked, the peptide is assigned one of three different
                            categories.",
                         br(),
                         "",
                         br(),
                         "The first category is double Lysine peptides (Includes C-termianl K), and 
                            the vast majority of tryptic cross-linked peptides in need of adjustment will be in
                            this category. Here the assigned link is simply reassigned to the correct and internal
                            Lysine.",
                         br(),
                         "",
                         br(),
                         "The second category is triple or more Lysine peptides. Here more than one possibility 
                            is available for the correction. To ensure correct assignment the original spectra in
                            Merox should be considered.",
                         br(),
                         "",
                         br(),
                         "The last category is Serine/Threonine cross-linked peptides. They do not contain any
                            internal Lysines. For this reason they can not be corrected.",
                         br(),
                         "",
                         br(),
                         "Author: Victor Chrone
                            ")
)
server <- function(input, output, session){
  options(shiny.maxRequestSize=2000*1024^2)
  
  
  observeEvent(input$submit,{
    
    if(input$checkbox == FALSE){
      data <-reactive({
        inFile <- input$filedata
        if (is.null(inFile))
          return(NULL)
        df <- read.csv2(inFile$datapath, header = TRUE)
        return(df)
      }) 
    }else{
      data <- reactive({
        inFile <- input$filedata
        if (is.null(inFile))
          return(NULL)
        df <- read.csv(inFile$datapath, header = TRUE)
        return(df)
      })
    }
    fasta <- reactive({
      inFile2 <- input$fasta
      if (is.null(inFile2))
        return(NULL)
      fastadat<-seqinr::read.fasta(inFile2$datapath)
      return(fastadat)
    })
    
    
    
    
    #MeroxAdjust
    
    
    correctedf<-MeroxAdjust1(data())
    
    countunique<-unique(correctedf$Protein.1)
    if(input$checkbox == FALSE){
      homeocheck(correctedf)
    }
    
    if(input$checkbox == FALSE){
      output$onek <-renderDT(onek[,c(1,7,8,11,12,21,22)])
      
      output$threek <-renderDT(threek[,c(1,7,8,11,12,21,22)])
      
      output$corrected<-renderDT(correctedf[,c(1,7,8,11,12,21,22,34)])
      
      output$corrected2<-renderDT(correcteddfneg[,c(1,7,8,11,12,21,22,34)])
      
      output$homeo<-renderDT(homeotypicdata[,c(1,7,8,11,12,21,22)])
      
      output$serine<-renderDT(serinethreonine[,c(1,7,8,11,12,21,22)])
      
      output$zero<-renderText(str_c("Dead-end linkers: ",linktypes$zero," "))
      output$one<-renderText(str_c("Inter peptide links: ",linktypes$one," "))
      output$rest<-renderText(str_c("Cross-links: ",linktypes$rest," "))
      output$total<-renderText(str_c("Total corrected: ", dim(filter(correctedf, corrected== "TRUE"))[1]," "))
      output$sandt<-renderText(str_c("Serine/Threonine linked: ",(dim(onek)[1]+dim(serinethreonine)[1])," "))
      output$homeo2<-renderText(str_c("Homeotypic links: ", dim(homeotypicdata)[1]," "))
      
    }else{
      output$onek <-renderDT(onek[,c(1,2,3,4,5,6,7,8,26)])
      
      output$threek <-renderDT(threek[,c(1,2,3,4,5,6,7,8,26)])
      
      output$corrected<-renderDT(correctedf[,c(1,2,3,4,5,6,7,8,26)])
      
      output$zero<-renderText(str_c("Dead-end linkers: ",linktypes$zero," "))
      output$one<-renderText(str_c("Inter peptide links: ",linktypes$one," "))
      output$rest<-renderText(str_c("Cross-links: ",linktypes$rest," "))
      output$total<-renderText(str_c("Total corrected: ", dim(filter(correctedf, Corrected== "TRUE"))[1]," "))
      output$sandt<-renderText(str_c("Serine/Threonine linked: ",(dim(onek)[1])," "))
      output$homeo2<-renderText(str_c("Homeotypic links: ", "NA"," "))
      
    }
    
    #for download
    
    
    if(input$checkbox == TRUE){
      
      if(is_empty(onek)==FALSE){
        onek<-onek[,1:26]
        colnames(onek)<-c("PepSeq1","PepSeq2","PepPos1","PepPos2","LinkPos1","LinkPos2","Protein1","Protein2","Charge","CrossLinkerModMass","ScanId","Scan","PeakListFileName","Rank","FragmentTolerance","IonTypes","ExpMz","CalcMz","Decoy1","Decoy2","Score","isDecoy","isTT","isTD","isDD","Corrected")
      }
      if(is_empty(threek)==FALSE){
        threek<-threek[,1:26]
        colnames(threek)<-c("PepSeq1","PepSeq2","PepPos1","PepPos2","LinkPos1","LinkPos2","Protein1","Protein2","Charge","CrossLinkerModMass","ScanId","Scan","PeakListFileName","Rank","FragmentTolerance","IonTypes","ExpMz","CalcMz","Decoy1","Decoy2","Score","isDecoy","isTT","isTD","isDD","Corrected")
      }
      
      correctedf<-correctedf[,1:26]
      
      colnames(correctedf)<-c("PepSeq1","PepSeq2","PepPos1","PepPos2","LinkPos1","LinkPos2","Protein1","Protein2","Charge","CrossLinkerModMass","ScanId","Scan","PeakListFileName","Rank","FragmentTolerance","IonTypes","ExpMz","CalcMz","Decoy1","Decoy2","Score","isDecoy","isTT","isTD","isDD","Corrected")
      
      for(i in 1:dim(correctedf)[1]){
        
        correctedf$PepSeq1[i]<-str_split(correctedf$PepSeq1[i],"\\[")[[1]][2]
        correctedf$PepSeq1[i]<-str_split(correctedf$PepSeq1[i],"\\]")[[1]][1]
        correctedf$PepSeq2[i]<-str_split(correctedf$PepSeq2[i],"\\[")[[1]][2]
        correctedf$PepSeq2[i]<-str_split(correctedf$PepSeq2[i],"\\]")[[1]][1]
        
        
        correctedf$LinkPos1[i]<-str_split(correctedf$LinkPos1[i],"")[[1]][2]
        correctedf$LinkPos2[i]<-str_split(correctedf$LinkPos2[i],"")[[1]][2]
        
        correctedf$Protein1[i]<-str_split(correctedf$Protein1[i],">")[[1]][2]
        correctedf$Protein2[i]<-str_split(correctedf$Protein2[i],">")[[1]][2]
        
      }
      
      
      output$downloadData<-downloadHandler(
        filename = function(){ paste(input$filedata,"_MeroxAdjust",".csv",sep = "")
        }, 
        content = function(file){write.table(correctedf[,1:25],file,sep=",",row.names = F,quote=F)
        }
      )
    }else{
      output$downloadData<-downloadHandler(
        filename = function(){ paste(input$filedata,"_MeroxAdjust",".xlsx",sep = "")
        }, 
        content = function(file){write_xlsx(list("Corrected" = correctedf ,"One_K"=onek,"Three_K"= threek),file)
        }
      )
      
      
      updateSelectInput(
        session,
        "Chain",
        choices = unique(correctedf$Protein.1)
      )
      
      updateSelectInput(
        session,
        "protein1",
        choices = unique(correctedf$Protein.1)
      )
      updateSelectInput(
        session,
        "protein2",
        choices = unique(correctedf$Protein.1)
      )
      
      
      
      output$plot_1<-renderPlot(lightchainplot(barlist[[input$Chain]],input$Chain,input$yaxis))
      
      output$plot_3<-renderPlot(scatplotSLICE1(oneandrest[1],0,0,"Scatterplot",fasta(),input$protein1,input$protein2))
    }
  })
  
  
  
}
shinyApp(ui = ui, server = server)




