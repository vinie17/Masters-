require(gridExtra)

#Functions ####
filterdata<-function(zero,one,rest,chain){
  ocreAzero<-as.data.frame(na.omit(zero))
  ocreAone<-as.data.frame(na.omit(rest))
  ocreArest<-as.data.frame(na.omit(one))
  
  #for zero
  for(i in 1:dim(ocreAzero)[1]){
    ocreAzero$chainA[i]<-str_split(ocreAzero[i,1],"\\|")[[1]][3]
    ocreAzero$resid1[i]<-str_split(ocreAzero[i,1],"\\|")[[1]][1]
    ocreAzero$color[i]<-"type 0"
  }
  #for one
  for(i in 1:dim(ocreAone)[1]){
    ocreAone$chainA[i]<-str_split(ocreAone[i,1],"\\|")[[1]][3]
    ocreAone$chainB[i]<-str_split(ocreAone[i,1],"\\|")[[1]][4]
    ocreAone$redid1[i]<-str_split(ocreAone[i,1],"\\|")[[1]][1]
    ocreAone$resid2[i]<-str_split(ocreAone[i,1],"\\|")[[1]][2]
    ocreAone$color[i]<-"type 1"
  }
  #for rest
  for(i in 1:dim(ocreArest)[1]){
    ocreArest$chainA[i]<-str_split(ocreArest[i,1],"\\|")[[1]][3]
    ocreArest$chainB[i]<-str_split(ocreArest[i,1],"\\|")[[1]][4]
    ocreArest$redid1[i]<-str_split(ocreArest[i,1],"\\|")[[1]][1]
    ocreArest$resid2[i]<-str_split(ocreArest[i,1],"\\|")[[1]][2]
    ocreArest$color[i]<-"type 2"
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
lightchainplot<-function(dat,title){
  a<-ggplot(data = dat,aes(x=factor(resid), fill = color))+
    geom_bar(position = "stack")+labs(title=title,x= "Residue number")+
    scale_y_continuous(limits=c(0,300))+
    theme(axis.text.x = element_text(angle = 90,size = 8, face = "bold"))+
    scale_color_manual(values = c("type 1"="skyblue",
                                  "type 2"="lightgreen",
                                  "type 0"="coral",
                                  "No hits"="grey100"))+
    scale_fill_manual(values = c("type 1"="skyblue",
                                 "type 2"="lightgreen",
                                 "type 0"="coral",
                                 "No hits"="grey100"))
  return(a)
}
heavychainplot<-function(dat,title){
  a<-ggplot(data = dat,aes(x=factor(resid), fill = color))+
    geom_bar(position = "stack")+labs(title=title,x= "Residue number")+
    scale_y_continuous(limits=c(0,300))+
    theme(axis.text.x = element_text(angle = 90,size = 8, face = "bold"))+
    
    scale_color_manual(values = c("type 1"="skyblue",
                                  "type 2"="lightgreen",
                                  "type 0"="coral",
                                  "No hits"="grey100"))+
    scale_fill_manual(values = c("type 1"="skyblue",
                                 "type 2"="lightgreen",
                                 "type 0"="coral",
                                 "No hits"="grey100"))
  return(a)
}
filloutmissingres<-function(unique,dat,chain){
  resid<-unique[!unique %in% dat$resid]
  chainA<-c(rep(chain,length(resid)))
  color<-c(rep("No hits",length(resid)))
  remaining<-data.frame(chainA,resid,color)
  dat<-droplevels(dat)
  dat<-rbind(dat,remaining)
  dat$resid<-as.numeric(as.character(dat$resid))
  dat<-as.data.frame(dat[order(dat$resid),])
  return(dat)
}

#barplot ####
#c4 ####

ocreAlight1<-filterdata(linktypelistzerores[[1]],
                        linktypelistrestres[[1]],
                        linktypelistoneres[[1]],
                        ">light")
ocreAlight11<-filterdata(linktypelistzerores[[1]],
                        linktypelistrestres[[1]],
                        linktypelistoneres[[1]],
                        ">light ")
ocreAlight1<-rbind(ocreAlight1,ocreAlight11)


ocreAlight2<-filterdata(linktypelistzerores[[2]],
                        linktypelistrestres[[2]],
                        linktypelistoneres[[2]],
                        ">light")
ocreAlight22<-filterdata(linktypelistzerores[[2]],
                        linktypelistrestres[[2]],
                        linktypelistoneres[[2]],
                        ">light ")
ocreAlight2<-rbind(ocreAlight2,ocreAlight22)


ocreAlight3<-filterdata(linktypelistzerores[[3]],
                        linktypelistrestres[[3]],
                        linktypelistoneres[[3]],
                        ">light")
ocreAlight4<-filterdata(linktypelistzerores[[4]],
                        linktypelistrestres[[4]],
                        linktypelistoneres[[4]],
                        ">light")
ocreAlight5<-filterdata(linktypelistzerores[[5]],
                        linktypelistrestres[[5]],
                        linktypelistoneres[[5]],
                        ">E10B10-L")
ocreAlight6<-filterdata(linktypelistzerores[[6]],
                        linktypelistrestres[[6]],
                        linktypelistoneres[[6]],
                        ">E10B10-L")
ocreAlight7<-filterdata(linktypelistzerores[[7]],
                        linktypelistrestres[[7]],
                        linktypelistoneres[[7]],
                        ">E10B10-L")

tysAlight1<-filterdata(linktypelistzerores[[4]],
                        linktypelistrestres[[4]],
                        linktypelistoneres[[4]],
                        ">light")
tysAlight2<-filterdata(linktypelistzerores[[5]],
                        linktypelistrestres[[5]],
                        linktypelistoneres[[5]],
                        ">light")
tysAlight3<-filterdata(linktypelistzerores[[6]],
                        linktypelistrestres[[6]],
                        linktypelistoneres[[6]],
                        ">light")
#heavy
ocreAheavy1<-filterdata(linktypelistzerores[[1]],
                        linktypelistrestres[[1]],
                        linktypelistoneres[[1]],
                        ">heavy")
ocreAheavy2<-filterdata(linktypelistzerores[[2]],
                        linktypelistrestres[[2]],
                        linktypelistoneres[[2]],
                        ">heavy")
ocreAheavy3<-filterdata(linktypelistzerores[[3]],
                        linktypelistrestres[[3]],
                        linktypelistoneres[[3]],
                        ">1IGA_2|Chains C, D|IGA1|Homo sapiens  9606")
ocreAheavy4<-filterdata(linktypelistzerores[[4]],
                        linktypelistrestres[[4]],
                        linktypelistoneres[[4]],
                        ">E10B10-H")
ocreAheavy5<-filterdata(linktypelistzerores[[5]],
                        linktypelistrestres[[5]],
                        linktypelistoneres[[5]],
                        ">E10B10-H")
ocreAheavy6<-filterdata(linktypelistzerores[[6]],
                        linktypelistrestres[[6]],
                        linktypelistoneres[[6]],
                        ">E10B10-H")
ocreAheavy7<-filterdata(linktypelistzerores[[7]],
                        linktypelistrestres[[7]],
                        linktypelistoneres[[7]],
                        ">E10B10-H")
tysAheavy1<-filterdata(linktypelistzerores[[4]],
                        linktypelistrestres[[4]],
                        linktypelistoneres[[4]],
                        ">heavy")
tysAheavy2<-filterdata(linktypelistzerores[[5]],
                        linktypelistrestres[[5]],
                        linktypelistoneres[[5]],
                        ">heavy")
tysAheavy3<-filterdata(linktypelistzerores[[6]],
                        linktypelistrestres[[6]],
                        linktypelistoneres[[6]],
                        ">heavy")

#for light chain ####
ocreAlight1<-filterdata(linktypelistzerores[[1]],
                        linktypelistrestres[[1]],
                        linktypelistoneres[[1]],
                        ">light")
ocreAlight2<-filterdata(linktypelistzerores[[2]],
                        linktypelistrestres[[2]],
                        linktypelistoneres[[2]],
                        ">light")
ocreAlight3<-filterdata(linktypelistzerores[[3]],
                        linktypelistrestres[[3]],
                        linktypelistoneres[[3]],
                        ">light")
ocreAlight4<-filterdata(linktypelistzerores[[4]],
                        linktypelistrestres[[4]],
                        linktypelistoneres[[4]],
                        ">light")
ocreAlight5<-filterdata(linktypelistzerores[[5]],
                        linktypelistrestres[[5]],
                        linktypelistoneres[[5]],
                        ">light")
ocreAlight6<-filterdata(linktypelistzerores[[6]],
                        linktypelistrestres[[6]],
                        linktypelistoneres[[6]],
                        ">light")
ocreAlight7<-filterdata(linktypelistzerores[[7]],
                        linktypelistrestres[[7]],
                        linktypelistoneres[[7]],
                        ">light")
ocreAlight8<-filterdata(linktypelistzerores[[8]],
                        linktypelistrestres[[8]],
                        linktypelistoneres[[8]],
                        ">light")
#for heavy chain
ocreAheavy1<-filterdata(linktypelistzerores[[1]],
                        linktypelistrestres[[1]],
                        linktypelistoneres[[1]],
                        ">heavy")
ocreAheavy2<-filterdata(linktypelistzerores[[2]],
                        linktypelistrestres[[2]],
                        linktypelistoneres[[2]],
                        ">heavy")
ocreAheavy3<-filterdata(linktypelistzerores[[3]],
                        linktypelistrestres[[3]],
                        linktypelistoneres[[3]],
                        ">heavy")
ocreAheavy4<-filterdata(linktypelistzerores[[4]],
                        linktypelistrestres[[4]],
                        linktypelistoneres[[4]],
                        ">heavy")
ocreAheavy5<-filterdata(linktypelistzerores[[5]],
                        linktypelistrestres[[5]],
                        linktypelistoneres[[5]],
                        ">heavy")
ocreAheavy6<-filterdata(linktypelistzerores[[6]],
                        linktypelistrestres[[6]],
                        linktypelistoneres[[6]],
                        ">heavy")
ocreAheavy7<-filterdata(linktypelistzerores[[7]],
                        linktypelistrestres[[7]],
                        linktypelistoneres[[7]],
                        ">heavy")
ocreAheavy8<-filterdata(linktypelistzerores[[8]],
                        linktypelistrestres[[8]],
                        linktypelistoneres[[8]],
                        ">heavy")


#for early plots
#find unique for light chain
finduni<-rbind(ocreAlight1,ocreAlight2,ocreAlight3,ocreAlight4)
finduni<-rbind(ocreAlight5,ocreAlight6,ocreAlight7,ocreAlight8)
#,ocreAlight5)#,ocreAlight6)#,ocreAlight7)
finduni<-unique(finduni$resid)

#finduni<-rbind(tysAlight1,tysAlight2,tysAlight3)
#finduni<-unique(finduni$resid)
#fill out missing for light
ocreAlight1<-filloutmissingres(finduni,ocreAlight1,">light")
ocreAlight2<-filloutmissingres(finduni,ocreAlight2,">light")
ocreAlight3<-filloutmissingres(finduni,ocreAlight3,">light")
ocreAlight4<-filloutmissingres(finduni,ocreAlight4,">light")
ocreAlight5<-filloutmissingres(finduni,ocreAlight5,">light")
ocreAlight6<-filloutmissingres(finduni,ocreAlight6,">light")
ocreAlight7<-filloutmissingres(finduni,ocreAlight7,">light")
ocreAlight8<-filloutmissingres(finduni,ocreAlight8,">light")

#tysAlight1<-filloutmissingres(finduni,tysAlight1,">light")
#tysAlight2<-filloutmissingres(finduni,tysAlight2,">light")
#tysAlight3<-filloutmissingres(finduni,tysAlight3,">light")

#find unique for heavy chain 
finduni<-rbind(ocreAheavy1,ocreAheavy2,ocreAheavy3,ocreAheavy4)
finduni<-rbind(ocreAheavy5,ocreAheavy6,ocreAheavy7,ocreAheavy8)#,ocreAheavy5)#,ocreAheavy6)#,ocreAheavy7)
finduni<-unique(finduni$resid)

#finduni<-rbind(tysAheavy1,tysAheavy2,tysAheavy3)#,ocreAheavy5)
#finduni<-unique(finduni$resid)

#fill out missing for heavy chain
ocreAheavy1<-filloutmissingres(finduni,ocreAheavy1,">heavy")
ocreAheavy2<-filloutmissingres(finduni,ocreAheavy2,">heavy")
ocreAheavy3<-filloutmissingres(finduni,ocreAheavy3,">heavy")
ocreAheavy4<-filloutmissingres(finduni,ocreAheavy4,">heavy")
ocreAheavy5<-filloutmissingres(finduni,ocreAheavy5,">heavy")
ocreAheavy6<-filloutmissingres(finduni,ocreAheavy6,">heavy")
ocreAheavy7<-filloutmissingres(finduni,ocreAheavy7,">heavy")
ocreAheavy8<-filloutmissingres(finduni,ocreAheavy8,">heavy")

#tysAheavy1<-filloutmissingres(finduni,tysAheavy1,">heavy")
#tysAheavy2<-filloutmissingres(finduni,tysAheavy2,">heavy")
#tysAheavy3<-filloutmissingres(finduni,tysAheavy3,">heavy")


#bar plot for light chain
a<-lightchainplot(ocreAlight1,linklist[1])
b<-lightchainplot(ocreAlight2,linklist[2])
c<-lightchainplot(ocreAlight3,linklist[3])
d<-lightchainplot(ocreAlight4,linklist[4])
e<-lightchainplot(ocreAlight5,linklist[5])
f<-lightchainplot(ocreAlight6,linklist[6])
g<-lightchainplot(ocreAlight7,linklist[7])
h<-lightchainplot(ocreAlight8,linklist[8])

#a<-lightchainplot(tysAlight1,linklist[4])
#b<-lightchainplot(tysAlight2,linklist[5])
#c<-lightchainplot(tysAlight3,linklist[6])

#arrange in grid light chain
gg<-grid.arrange(a,b,c,d,nrow=4)
ggg<-grid.arrange(e,f,g,h,nrow=4)#,e,f,g)#,c,d)#,d,e)
gg<-grid.arrange(e,d,c,b,a,ncol=1)
gg<-grid.arrange(a,b,c,d,e,f,g,h,ncol=1)
#bar plot for heavy chain
a<-heavychainplot(ocreAheavy1,linklist[1])
b<-heavychainplot(ocreAheavy2,linklist[2])
c<-heavychainplot(ocreAheavy3,linklist[3])
d<-heavychainplot(ocreAheavy4,linklist[4])
e<-heavychainplot(ocreAheavy5,linklist[5])
f<-heavychainplot(ocreAheavy6,linklist[6])
g<-heavychainplot(ocreAheavy7,linklist[7])
h<-heavychainplot(ocreAheavy8,linklist[8])

a<-heavychainplot(tysAheavy1,linklist[4])
b<-heavychainplot(tysAheavy2,linklist[5])
c<-heavychainplot(tysAheavy3,linklist[6])
#arrange in grid heavy chain
grid.arrange(gg,ggg,ncol=2)#,d,e)



#For HEATMAPS ####
###for light chain
ocreAlight1<-as.data.frame(table(ocreAlight1$resid))
ocreAlight2<-as.data.frame(table(ocreAlight2$resid))
ocreAlight3<-as.data.frame(table(ocreAlight3$resid))
ocreAlight4<-as.data.frame(table(ocreAlight4$resid))
ocreAlight5<-as.data.frame(table(ocreAlight5$resid))
ocreAlight6<-as.data.frame(table(ocreAlight6$resid))
ocreAlight7<-as.data.frame(table(ocreAlight7$resid))
### for heavy chain 
ocreAlight1<-as.data.frame(table(ocreAheavy1$resid))
ocreAlight2<-as.data.frame(table(ocreAheavy2$resid))
ocreAlight3<-as.data.frame(table(ocreAheavy3$resid))
ocreAlight4<-as.data.frame(table(ocreAheavy4$resid))
ocreAlight5<-as.data.frame(table(ocreAheavy5$resid))

ocreAlight1$Var1<-as.numeric(as.character(ocreAlight1$Var1))
ocreAlight2$Var1<-as.numeric(as.character(ocreAlight2$Var1))
ocreAlight3$Var1<-as.numeric(as.character(ocreAlight3$Var1))
ocreAlight4$Var1<-as.numeric(as.character(ocreAlight4$Var1))
ocreAlight5$Var1<-as.numeric(as.character(ocreAlight5$Var1))
ocreAlight6$Var1<-as.numeric(as.character(ocreAlight6$Var1))
ocreAlight7$Var1<-as.numeric(as.character(ocreAlight7$Var1))


ocreAlight1<-ocreAlight1[order(ocreAlight1$Var1),]
ocreAlight2<-ocreAlight2[order(ocreAlight2$Var1),]
ocreAlight3<-ocreAlight3[order(ocreAlight3$Var1),]
ocreAlight4<-ocreAlight4[order(ocreAlight4$Var1),]
ocreAlight5<-ocreAlight5[order(ocreAlight5$Var1),]
ocreAlight6<-ocreAlight6[order(ocreAlight6$Var1),]
ocreAlight7<-ocreAlight7[order(ocreAlight7$Var1),]

ocreAlight1$Var1<-factor(ocreAlight1$Var1)
ocreAlight2$Var1<-factor(ocreAlight2$Var1)
ocreAlight3$Var1<-factor(ocreAlight3$Var1)
ocreAlight4$Var1<-factor(ocreAlight4$Var1)
ocreAlight5$Var1<-factor(ocreAlight5$Var1)
ocreAlight6$Var1<-factor(ocreAlight6$Var1)
ocreAlight7$Var1<-factor(ocreAlight7$Var1)



ocreAlight<-cbind(ocreAlight1[,2],ocreAlight2[,2],ocreAlight3[,2],ocreAlight4[,2],ocreAlight5[,2])#,ocreAlight6[,2],ocreAlight7[,2])
rownames(ocreAlight)<-ocreAlight1$Var1
library(gplots)
colnames(ocreAlight)<-c(linklist[1],linklist[2],linklist[3],linklist[4],linklist[5])#,linklist[6],linklist[7])#,linklist[3])
dat<-as.matrix(ocreAlight)
heatmap.2(dat,scale = "column",trace = "none", density.info = "none",
          col = bluered(100),xlab = "Condition",ylab = "Resid",cexRow = 0.5,
          cexCol = 0.8,srtCol = 30)





