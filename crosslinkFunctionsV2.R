####
#Data transformation functions for masters thesis
####
# Libraries
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library("tidyverse")
library(egg)
library("rstudioapi")
library("readxl")


lightChain<-c('DIQMTQSPSSLSASVGDRVTITCKTSQDINKYMAWYQQTPGKAPRLLIHYTSALQPGIPSRFSGSGSGRDYTFTISSLQPEDIATYYCLQYDNLWTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
heavyChain<-c('QVQLVQSGAEVKKPGASVKVSCKASGFNIKDTYIHWVRQAPGQRLEWMGRIDPANGYTKYDPKFQGRVTITADTSASTAYMELSSLRSEDTAVYYCAREGYYGNYGVYAMDYWGQGTLVTVSSASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTKTYTCNVDHKPSNTKVDKRVESKYGPPCPSCPAPEFLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPEVQFNWYVDGVEVHNAKTKPREEQFNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKGLPSSIEKTISKAKGQPREPQVYTLPPSQEEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSRLTVDKSRWQEGNVFSCSVMHEALHNHYTQKSLSLSLGK')

linkImporter<-function(file,modeltypes,filenames,links,antibodys,directlyfromPDB,correctcordinates){
  dataimport(file,modeltypes,filenames,links,antibodys,directlyfromPDB,correctcordinates)
  print("----------------------------")
  print("Data succesfully imported!  ")
  print("----------------------------")
  
  
  
  crosslinkclassification(processedlinks= crosslinkdist)
  print("----------------------------")
  print("Data succesfully classified!")
  print("----------------------------")
  

}


crosslinkplotspecific<-function(data,filteringafter,IgGIgA,lightchain,heavychain,jitx,jity,multiple){
  classfication(data,filteringafter)
  
  print("----------------------------")
  print("Data imported and filtered! ")
  print("----------------------------")
  
  Constraints(distanalysis1)
 
  print("----------------------------")
  print("Data color coded!           ")
  print("----------------------------")
  
  chainsplit(distanalysis1,IgGIgA, lightchain, heavychain)
  
  print("----------------------------")
  print("Light and heavy chain split!")
  print("----------------------------")
  
  if(isTRUE(multiple)){
    crosslinkplotMultiple(scatterplotdist, lightchain,heavychain,jitx,jity)
  }else{
    crosslinkplotSingle(scatterplotdist, lightchain,heavychain,jitx,jity)
    
  }
  print("----------------------------")
  print("Printing plot!              ")
  print("----------------------------")
}


conpeppdb<-function(datadis){
  datadis$Lys1<-as.character(datadis$Lys1)
  datadis$Lys2<-as.character(datadis$Lys2)
  for(i in 1:dim(datadis)[1]){
    datadis$concpepnr[i]<-str_c(datadis$Lys1[i],datadis$Lys2[i],collapse = ",")
  }
  datadis<<-datadis
}

#save for unique pdb
pdbmatch<-function(links,pdbs){
  pdb<-data.frame()
  forpdb<-data.frame()
  for(i in 1:length(links)){
    for(j in 1:dim(pdbs)){
      if(links[i] == pdbs$concpepnr[j]){
        pdb[i,1:11]<-pdbs[j,1:11]
        pdb$link[i]<-"LINK"
        pdb$CA[i]<-"CA"
        pdb$LYS[i]<-"LYS"
        pdb$blank[i]<-"        "
        
      }
    }
  }
  forpdb<-cbind(pdb[12],pdb[15],pdb[13],pdb[14],pdb[1:2],pdb[13],pdb[14],pdb[3:4],pdb[15],pdb[5])
  write.table(forpdb,file = "forpdb.txt", row.names = F,col.names = F,sep = "\t")
  pdb<<-pdb
  forpdb<<-forpdb
}

#import all info for linkimporter using xlsx 
import_from_csv<-function(csvfile){
  filedata<-read_excel(csvfile,col_names = T)
  
  for(i in 1:dim(filedata)[1]){
  linkImporter((filedata$dir[i]),filedata$model[i],filedata$filename[i],filedata$type[i],filedata$antibody[i],TRUE)
  }
  filedata<<-filedata
}

#links<-crosslinkdist
coordinatecorrect<-function(bool,oldlinks){
if(isTRUE(bool)){
for(i in 1:dim(oldlinks)[1]){
  if(oldlinks$type1[i] == "A"){
    oldlinks$type1[i] = "D"
  }else if(oldlinks$type1[i] == "D"){
    oldlinks$type1[i] = "A"
  }
  if(oldlinks$type2[i] == "A"){
    oldlinks$type2[i] = "D"
  }else if(oldlinks$type2[i] == "D"){
    oldlinks$type2[i] = "A"
  }
}
}
  return(oldlinks)
}

#read.csv2("linkdataimport.csv",header =T)
processFile = function(filepath){
  con = file(filepath,"r")
  linkrecord<-c()
  while(TRUE){
    line = readLines(con,n=1)
    if(length(line)==0){
      break
    }
    if(startsWith(line,"LINK")){
      linkrecord<-append(linkrecord,line)
    }
    if(startsWith(line,"ATOM")){
      break
    }
  }
  close(con)
  linkrecord<<-linkrecord
}



toDataFrame<-function(linkrecorddata){
  
  linkrecord<-as.data.frame(linkrecorddata)
  links<-data.frame()
  tt<-data.frame()
  
  for(i in 1:dim(linkrecord)[1]){
    links<-as.data.frame(strsplit(linkrecord[i,]," "))
    #remove all empty slots
    t<-as.data.frame(links[!apply(links =="",1,all),])
    #transpose
    t<-as.data.frame(transpose(t))
    tt<-rbind(tt,t)
  }
  tt<-cbind(tt[,4],tt[,5],tt[,8],tt[,9],tt[,10])
  colnames(tt)<-c("type1","Lys1","type2","Lys2","Distance")
  linkout<<-tt
}

dataimport<-function(file,modeltypes,filenames,link,antibodys,directlyfromPDB,correctcordinates){
  if(directlyfromPDB == F){
    #import file
    linkdata<-read.table(file, header = F)
  
    #Take the columns we need
    crosslinkdist<-cbind(linkdata[4:5], linkdata[8:9],linkdata[10])
    #change colnames to what we need
    colnames(crosslinkdist)<-c("type1","Lys1","type2","Lys2","Distance")
    #define model such as (Native or inflix. put inside "")
    crosslinkdist$model<-modeltypes
    #define filename put answer inside ""
    crosslinkdist$filename<-filenames
    #define link type such as short or all put answer inside ""
    crosslinkdist$links<-link
    #define the type of antibody such as ocrelizumab, tysabri or igA
    crosslinkdist$antibody<-antibodys
    #Classification of crosslink type
    crosslinkdist$Classification<-"none"
    #returs
    crosslinkdist<<-crosslinkdist
  }
  if(directlyfromPDB == T){
    #import pdb
    processFile(file)
    #transform to dataframe
    toDataFrame(linkrecord)
    
    #same as above
    crosslinkdist<-as.data.frame(linkout)
    #define model such as (Native or inflix. put inside "")
    crosslinkdist$model<-modeltypes
    #define filename put answer inside ""
    crosslinkdist$filename<-filenames
    #define link type such as short or all put answer inside ""
    crosslinkdist$links<-link
    #define the type of antibody such as ocrelizumab, tysabri or igA
    crosslinkdist$antibody<-antibodys
    #Classification of crosslink type
    crosslinkdist$Classification<-"none"
    
    #returs
    crosslinkdist<<-coordinatecorrect(correctcordinates,crosslinkdist)
    
  }
  crosslinkdist$Lys1<-as.double(crosslinkdist$Lys1)
  crosslinkdist$Lys2<-as.double(crosslinkdist$Lys2)
  crosslinkdist$Distance<-as.double(crosslinkdist$Distance)
  
}


crosslinkclassification<-function(processedlinks){
  #Classification of crosslink type
  processedlinks$Classification<-"none"
  
  for(i in 1:dim(processedlinks)[1]){
    if(processedlinks$type1[i] == processedlinks$type2[i] && processedlinks$Lys1[i] == processedlinks$Lys2[i]){
      processedlinks$Classification[i]<-"type 1 (same chain and same K (aritfact))"
    }
    else  if(processedlinks$type1[i] == processedlinks$type2[i]){
      processedlinks$Classification[i]<-"type 1 (same chain)"
    }
    else if(processedlinks$type1[i] == "A" && processedlinks$type2[i] == "C" && processedlinks$Lys1[i] == processedlinks$Lys2[i]){
      processedlinks$Classification[i]<-"type 2h"
    }
    else if(processedlinks$type1[i] == "A" && processedlinks$type2[i] == "C"){
      processedlinks$Classification[i]<-"type 1 (different chains)"
    }
    else if(processedlinks$type1[i] == "C" && processedlinks$type2[i] == "A" && processedlinks$Lys1[i] == processedlinks$Lys2[i]){
      processedlinks$Classification[i]<-"type 2h"
    }
    else if(processedlinks$type1[i] == "C" && processedlinks$type2[i] == "A"){
      processedlinks$Classification[i]<-"type 1 (different chains)"
    }
    else if(processedlinks$type1[i] == "B" && processedlinks$type2[i] == "D" && processedlinks$Lys1[i] == processedlinks$Lys2[i]){
      processedlinks$Classification[i]<-"type 2h"
    }
    else if(processedlinks$type1[i] == "B" && processedlinks$type2[i] == "D"){
      processedlinks$Classification[i]<-"type 1 (different chains)"
    } 
    else if(processedlinks$type1[i] == "D" && processedlinks$type2[i] == "B" && processedlinks$Lys1[i] == processedlinks$Lys2[i]){
      processedlinks$Classification[i]<-"type 2h"
    }
    else if(processedlinks$type1[i] == "D" && processedlinks$type2[i] == "B"){
      processedlinks$Classification[i]<-"type 1 (different chains)"
    }
    else{
      processedlinks$Classification[i]<-"type 2"
    }
  }
  #combine to the other sets
  if(exists("distanalysis")){
    distanalysis<-processedlinks
    keepsafe<-rbind(keepsafe,distanalysis)
  }else{
    distanalysis<-processedlinks
    keepsafe<-processedlinks
  }
  
  distanalysis$Lys1<-as.double(distanalysis$Lys1)
  distanalysis$Lys2<-as.double(distanalysis$Lys2)
  distanalysis$Distance<-as.double(distanalysis$Distance)

  #retuns
  distanalysis<<-distanalysis
  keepsafe<<-keepsafe
  
}

  classfication<-function(data,filteringafter){
    if(filteringafter == "tys"|| filteringafter =="ocre"||filteringafter =="iga"){
    forplot1<-filter(data, antibody == filteringafter)
    }
   #more can be put in if needed 
 
  #Histogram of the number of inter intra
  a<-ggplot(forplot1, aes(x=Classification, type=model, color =antibody)) + 
    geom_histogram(position = "dodge", stat = "count")
  #Histogram of distances
  forplot2<-forplot1
  forplot2$Distance<-as.numeric(forplot2$Distance)
  b<-ggplot(forplot2, aes(x=Distance, fill=model)) +
    geom_histogram(position = "dodge")
    
  
  
  distanalysis1<<-forplot1
  linkdist<<-b
  
  
  
  }
 
  #geat mean
  meanmodel<-function(model){
    xx<-data.frame("AVG")
    
    for(i in 1:dim(scatterplotdist)[1]){
      if(scatterplotdist$model[i] == model){
        xx[i,]<-scatterplotdist$Distance[i]
      }
      else{ xx[i,]<-NA}
    }
    xx<-na.omit(xx)
    mean(as.numeric(xx$X.AVG.))
  }
  
  #get mean of 80%
  meanmodel80<-function(model){
    xx<-data.frame("AVG")
    
    for(i in 1:dim(scatterplotdist)[1]){
      if(scatterplotdist$model[i] == model){
        xx[i,]<-scatterplotdist$Distance[i]
      }
      else{ xx[i,]<-NA}
    }
    xx<-na.omit(xx)
    xx<-xx[order(xx),]
    xx<-xx[xx <= quantile(as.numeric(xx),0.8)]
    mean(as.numeric(xx))
  }
  
 
  Constraints<-function(distanalysis1){
  #classification of distance by color
  distanalysis1$Constraints<-"none"
  for(i in 1:dim(distanalysis1)[1]){
    if(0 <= (distanalysis1$Distance[i]) && (distanalysis1$Distance[i])  <= 30){
      distanalysis1$Constraints[i]<-"Green"
    }
    else if(30.01 <= (distanalysis1$Distance[i]) && (distanalysis1$Distance[i])  <= 40){
      distanalysis1$Constraints[i]<-"Orange"
    }else{
      distanalysis1$Constraints[i]<-"Red"
    }
    
  }
  distanalysis1<<-distanalysis1
  }

  chainsplit<-function(distanalysis1, IgGIgA, lightchain, heavychain){
    #allowing for distiinciton between chains
    scatterplotdist<-distanalysis1
    scatterplotdist$Lys1<-as.double(scatterplotdist$Lys1)
    scatterplotdist$Lys2<-as.double(scatterplotdist$Lys2)
    scatterplotdist$Distance<-as.double(scatterplotdist$Distance)
    lightch<-c(strsplit(lightchain,""))
    heavych<-c(strsplit(heavychain,""))
    
    for(i in 1:dim(scatterplotdist)[1]){
      if(IgGIgA == "IgG" & scatterplotdist$type1[i] == "A" || scatterplotdist$type1[i]== "C"){
        if(scatterplotdist$Lys1[i] <= length(lightch[[1]])){
          scatterplotdist$Lys1[i] <- (scatterplotdist$Lys1[i] + length(heavych[[1]]))
        }
        #if(scatterplotdist$Lys2[i] <= length(lightch[[1]])){
         # scatterplotdist$Lys2[i] <- (scatterplotdist$Lys2[i] + length(heavych[[1]]))
        #}
      }
      if(IgGIgA == "IgG" & scatterplotdist$type2[i] == "C" || scatterplotdist$type2[i] == "A"){
        if(scatterplotdist$Lys2[i] <= length(lightch[[1]])){
          scatterplotdist$Lys2[i] <- (scatterplotdist$Lys2[i] + length(heavych[[1]]))
        }
        #if(scatterplotdist$Lys2[i] <= length(lightch[[1]])){
         # scatterplotdist$Lys2[i] <- (scatterplotdist$Lys2[i] + length(heavych[[1]]))
        #}
      }
      if(IgGIgA == "IgA" & scatterplotdist$type1[i] == "A" || scatterplotdist$type1[i] == "C"){
        if(scatterplotdist$Lys1[i] <= length(lightch[[1]])){
          scatterplotdist$Lys1[i] <- (scatterplotdist$Lys1[i] + length(heavych[[1]]))
        }
        #if(scatterplotdist$Lys2[i] <= length(lightch[[1]])){
         # scatterplotdist$Lys2[i] <- (scatterplotdist$Lys2[i] + length(heavych[[1]]))
        #}
      }
      if(IgGIgA == "IgA" & scatterplotdist$type2 == "C" & scatterplotdist$type2 == "A"){
        if(scatterplotdist$Lys2[i] <= length(lightch[[1]])){
          scatterplotdist$Lys2[i] <- (scatterplotdist$Lys2[i] + length(heavych[[1]]))
        }
        #if(scatterplotdist$Lys2[i] <= length(lightch[[1]])){
         # scatterplotdist$Lys2[i] <- (scatterplotdist$Lys2[i] + length(heavych[[1]]))
        #}
      }
      
    }
    scatterplotdist<<-scatterplotdist
  }
  
  crosslinkplotMultiple<-function(scatterplotdist,lightchain,heavychain,jitx,jity){
    lightChain<-c('DIQMTQSPSSLSASVGDRVTITCKTSQDINKYMAWYQQTPGKAPRLLIHYTSALQPGIPSRFSGSGSGRDYTFTISSLQPEDIATYYCLQYDNLWTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
    heavyChain<-c('QVQLVQSGAEVKKPGASVKVSCKASGFNIKDTYIHWVRQAPGQRLEWMGRIDPANGYTKYDPKFQGRVTITADTSASTAYMELSSLRSEDTAVYYCAREGYYGNYGVYAMDYWGQGTLVTVSSASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTKTYTCNVDHKPSNTKVDKRVESKYGPPCPSCPAPEFLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPEVQFNWYVDGVEVHNAKTKPREEQFNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKGLPSSIEKTISKAKGQPREPQVYTLPPSQEEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSRLTVDKSRWQEGNVFSCSVMHEALHNHYTQKSLSLSLGK')
    #OCRE
    #heavyChain<-c('EVQLVESGGGLVQPGGSLRLSCAASGYTFTSYNMHWVRQAPGKGLEWVGAIYPGNGDTSYNQKFKGRFTISVDKSKNTLYLQMNSLRAEDTAVYYCARVVYYSNSYWYFDVWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK')
    #lightChain<-c('DIQMTQSPSSLSASVGDRVTITCRASSSVSYMHWYQQKPGKAPKPLIYAPSNLASGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQWSFNPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
    #INF
    #lightChain<-c('DILLTQSPAILSVSPGERVSFSCRASQFVGSSIHWYQQRTNGSPRLLIKYASESMSGIPSRFSGSGSGTDFTLSINTVESEDIADYYCQQSHSWPFTFGSGTNLEVKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
    lightChainC<-c('RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
    #heavyChain<-c('EVKLEESGGGLVQPGGSMKLSCVASGFIFSNHWMNWVRQSPEKGLEWVAEIRSKSINSATHYAESVKGRFTISRDDSKSAVYLQMTDLRTEDTGVYYCSRNYYGSTYDYWGQGTTLTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLS')
    heavyChainFc<-c('ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK')
    
    lightch<-c(strsplit(lightChain,""))
    lightchC<-c(strsplit(lightChainC,""))
    heavych<-c(strsplit(heavyChain,""))
    heavychFc<-c(strsplit(heavyChainFc,""))
    lightch<-as.double(length(lightch[[1]]))
    lightchC<-as.double(length(lightchC[[1]]))
    heavych<-as.double(length(heavych[[1]]))
    heavychFc<-as.double(length(heavychFc[[1]]))
    
    
    
    jitter<-position_jitter(width = jitx,height=jity)
    #scatter plot os distaces
    scatterplotdist<-arrange(scatterplotdist, by = model)
    scatterplotdist$model<-as.factor(scatterplotdist$model)
    p<-ggplot(scatterplotdist, aes(x=Lys1, y=Lys2, group=Classification))+
      geom_point(aes(shape=Classification, colour=Constraints),position =jitter)+
      scale_color_manual(values = c("#00BA38",'#619CFF',"#F8766D"), label = c("Below 30Å", "Between 30-40Å","Above 40Å"))+
      scale_size(range = c(1,11 ), name="Distance (Å)") +
      ylab("Hv                                  Hc                         Lv                         Lc                                         HC                                  HV                         LC                         LV")+
      xlab("Hv                                             Hc                         Lv                         Lc                                  HC                                             HV                           LC                         LV")+
      #scale_x_binned(n.breaks = (lightch2+heavych2) )+
      #scale_y_binned(n.breaks = (lightch2+heavych2))+
      annotate(geom="rect",xmin = 300, xmax=451,
               ymin = 460, ymax=630,
               color = NA,
               fill = "grey60",alpha = 0.3)+
      scale_x_continuous(expand = c(0,0), breaks = seq(0,(heavych+lightch+1+49), by =30), limits = c(0,(1+lightch+heavych+50)))+
      scale_y_continuous(expand = c(0,0), breaks= seq(0, (heavych+lightch+1+49), by =30), limits = c(0,(1+lightch+heavych+50)))+
      theme(axis.text.x = element_text(angle = 90, face = "bold",size = 7),
            axis.text.y = element_text(face = "bold", size = 7),
            #axis.title.y = element_text(size = 7, hjust = 1), #for big plot
            #axis.title.x = element_text(hjust = 1),
            axis.title.y = element_text(size = 6, hjust = 1, face = "bold"), #for grouped
            axis.title.x = element_text(hjust = 1, size = 7, face = "bold"),
            axis.text = element_text(
              color =  c(
                rep("grey1",length(seq(1,(heavych-heavychFc), by =30))),
                rep("grey40",length(seq((heavych-heavychFc),(heavych), by =30))),
                rep("grey15",length(seq(heavych,((lightch-lightchC)+heavych), by =30))),
                rep("grey55",length(seq(((lightch-lightchC)+heavych),(lightch+heavych+1+49), by =30))))))
    my_tag_models<-as.data.frame(unique(as.character(scatterplotdist$model)))
    my_tag_means<-as.data.frame(apply(my_tag_models,1 , meanmodel))
    my_tag<-c(apply(my_tag_means, 1 , round))
    #my_tag<-c(my_tag[1:3],my_tag[4])
    
    
    pp<-p + facet_wrap(~model)+theme(panel.spacing = unit(4,"lines"))                 
    
    ppp<-tag_facet(pp, 
                   x = 500, y = -Inf, 
                   vjust = -1, hjust = -0.5,
                   fontface = 4,
                   size = 4,
                   family = "serif",
                   tag_pool = my_tag)  
    
    my_tag2<-c(my_tag_models[1:3,],my_tag_models[4,])
    pppp<-tag_facet(ppp,
                    x = 400, y = -Inf, 
                    vjust = -1, hjust = -0.25,
                    open = "", close = "",
                    fontface = 4,
                    size = 4,
                    family = "serif",
                    tag_pool = my_tag2)
    linkscatter<<-pppp
  }
  


  
  crosslinkplotSingle<-function(scatterplotdist,lightchain,heavychain,jitx,jity){
    lightChain<-c('DIQMTQSPSSLSASVGDRVTITCKTSQDINKYMAWYQQTPGKAPRLLIHYTSALQPGIPSRFSGSGSGRDYTFTISSLQPEDIATYYCLQYDNLWTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
    heavyChain<-c('QVQLVQSGAEVKKPGASVKVSCKASGFNIKDTYIHWVRQAPGQRLEWMGRIDPANGYTKYDPKFQGRVTITADTSASTAYMELSSLRSEDTAVYYCAREGYYGNYGVYAMDYWGQGTLVTVSSASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTKTYTCNVDHKPSNTKVDKRVESKYGPPCPSCPAPEFLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPEVQFNWYVDGVEVHNAKTKPREEQFNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKGLPSSIEKTISKAKGQPREPQVYTLPPSQEEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSRLTVDKSRWQEGNVFSCSVMHEALHNHYTQKSLSLSLGK')
    #OCRE
    #heavyChain<-c('EVQLVESGGGLVQPGGSLRLSCAASGYTFTSYNMHWVRQAPGKGLEWVGAIYPGNGDTSYNQKFKGRFTISVDKSKNTLYLQMNSLRAEDTAVYYCARVVYYSNSYWYFDVWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK')
    #lightChain<-c('DIQMTQSPSSLSASVGDRVTITCRASSSVSYMHWYQQKPGKAPKPLIYAPSNLASGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQWSFNPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
    #INF
    #lightChain<-c('DILLTQSPAILSVSPGERVSFSCRASQFVGSSIHWYQQRTNGSPRLLIKYASESMSGIPSRFSGSGSGTDFTLSINTVESEDIADYYCQQSHSWPFTFGSGTNLEVKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
    lightChainC<-c('RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
    #heavyChain<-c('EVKLEESGGGLVQPGGSMKLSCVASGFIFSNHWMNWVRQSPEKGLEWVAEIRSKSINSATHYAESVKGRFTISRDDSKSAVYLQMTDLRTEDTGVYYCSRNYYGSTYDYWGQGTTLTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLS')
    heavyChainFc<-c('ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK')
    
    lightch<-c(strsplit(lightChain,""))
    lightchC<-c(strsplit(lightChainC,""))
    heavych<-c(strsplit(heavyChain,""))
    heavychFc<-c(strsplit(heavyChainFc,""))
    lightch<-as.double(length(lightch[[1]]))
    lightchC<-as.double(length(lightchC[[1]]))
    heavych<-as.double(length(heavych[[1]]))
    heavychFc<-as.double(length(heavychFc[[1]]))
    
    
    
    jitter<-position_jitter(width = jitx,height=jity)
    #scatter plot os distaces
    scatterplotdist<-arrange(scatterplotdist, by = model)
    scatterplotdist$model<-as.factor(scatterplotdist$model)
    p<-ggplot(scatterplotdist, aes(x=Lys1, y=Lys2, group=Classification))+
      geom_point(aes(shape=Classification, colour=Constraints),position =jitter)+
      scale_color_manual(values = c("#00BA38",'#619CFF',"#F8766D"), label = c("Below 30Å", "Between 30-40Å","Above 40Å"))+
      scale_size(range = c(1,11 ), name="Distance (Å)") +
      ylab("Hv                                                                                                     Hc                                                                     Lv                                                                  Lc        ")+
      xlab("Hv                                                                                                     Hc                                                                     Lv                                                                  Lc        ")+
      #scale_x_binned(n.breaks = (lightch2+heavych2) )+
      #scale_y_binned(n.breaks = (lightch2+heavych2))+
      annotate(geom="rect",xmin = 270, xmax=451,
               ymin = 460, ymax=660,
               color = NA,
               fill = "grey60",alpha = 0.3)+
      scale_x_continuous(expand = c(0,0), breaks = seq(0,(heavych+lightch+1+49), by =30), limits = c(0,(1+lightch+heavych+50)))+
      scale_y_continuous(expand = c(0,0), breaks= seq(0, (heavych+lightch+1+49), by =30), limits = c(0,(1+lightch+heavych+50)))+
      theme(axis.text.x = element_text(angle = 90, face = "bold",size = 7),
            axis.text.y = element_text(face = "bold", size = 7),
            #axis.title.y = element_text(size = 7, hjust = 1), #for big plot
            #axis.title.x = element_text(hjust = 1),
            axis.title.y = element_text(size = 6, hjust = 1, face = "bold"), #for grouped
            axis.title.x = element_text(hjust = 1, size = 7, face = "bold"),
            axis.text = element_text(
              color =  c(
                rep("grey1",length(seq(1,(heavych-heavychFc), by =30))),
                rep("grey40",length(seq((heavych-heavychFc),(heavych), by =30))),
                rep("grey15",length(seq(heavych,((lightch-lightchC)+heavych), by =30))),
                rep("grey55",length(seq(((lightch-lightchC)+heavych),(lightch+heavych+1+49), by =30))))))
                                 my_tag_models<-as.data.frame(unique(as.character(scatterplotdist$model)))
                                 my_tag_means<-as.data.frame(apply(my_tag_models,1 , meanmodel))
                                 my_tag<-c(apply(my_tag_means, 1 , round))
                                           #my_tag<-c(my_tag[1:3],my_tag[4])
                                           
                                 
                                           pp<-p + facet_wrap(~model)+theme(panel.spacing = unit(4,"lines"))                 
                                           
                                           ppp<-tag_facet(pp, 
                                                          x = 500, y = -Inf, 
                                                          vjust = -1, hjust = -0.5,
                                                         fontface = 4,
                                                           size = 4,
                                                        family = "serif",
                                                         tag_pool = my_tag)  
                                           
                                           my_tag2<-c(my_tag_models[1:3,],my_tag_models[4,])
                                           pppp<-tag_facet(ppp,
                                                           x = 400, y = -Inf, 
                                                           vjust = -1, hjust = -0.25,
                                                           open = "", close = "",
                                                           fontface = 4,
                                                           size = 4,
                                                           family = "serif",
                                                           tag_pool = my_tag2)
                                           linkscatter<<-pppp
                                           
  }
  
  