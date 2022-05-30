require(gridExtra)
####### FOR SLICES #################### #########################
#concatenate peptide 1 and peptide 2 so we known where links are between
#RTC_ST2<-RTC_ST
conpepres2(RTC_ST)
RTC_ST<-listdata

#conresid 
dataset<-list(RTC_SO[[1]],RTC_SO[[2]],RTC_SO[[3]],RTC_SO[[4]],RTC_SO[[5]])
dataset<-list(RTC_ST[[1]],RTC_ST[[2]],RTC_ST[[3]],RTC_ST[[4]],RTC_ST[[5]])
linklist<-c("S1(BOT)","S2","S3","S4","S5(TOP)")

linktypelistzero<-list(NA,NA,NA,NA,NA)
linktypelistone<-list(NA,NA,NA,NA,NA)
linktypelistrest<-list(NA,NA,NA,NA,NA)

linktypelistzerores<-list(NA,NA,NA,NA,NA)
linktypelistoneres<-list(NA,NA,NA,NA,NA)
linktypelistrestres<-list(NA,NA,NA,NA,NA)


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
linktypelisthomeotypicone<-list(NA,NA,NA,NA,NA)
dataone<-list(NA,NA,NA,NA,NA)
homeotypiclistone<-list(NA,NA,NA,NA,NA)
restlist<-list(NA,NA,NA,NA,NA)

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
for(i in 1:length(oneandrest)){
  oneandrest[[i]] = setmincount(na.omit(oneandrest[[i]]),2)
}
#find amount of other links so we can compare across
restlist<-list(NA,NA,NA,NA,NA)
oneandrest2<-oneandrest
findrest(dataset,oneandrest2)

#ratio of links to homeo
#inflix
174/1
916/24
559/34
478/13
174/1

#ocrelizumab
654/44
927/58
860/44
576/33
449/26

#tysabri
173/3
936/33
1072/23
501/15
206/1

#barplot of crosslink types
p<-list(NA,NA,NA)
crosslinktypeplot(p,linktypesres)


#combined with one(homeotypic removed)
O<- list(
  A = na.omit(oneandrest[[1]]), 
  B = na.omit(oneandrest[[2]]), 
  C = na.omit(oneandrest[[3]]),
  D = na.omit(oneandrest[[4]]),
  E = na.omit(oneandrest[[5]])
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

#tysabri #########################################

#normal
A<-list(oneandrest[[1]])
B<-list(oneandrest[[2]])
C<-list(oneandrest[[3]])
D<-list(oneandrest[[4]])
E<-list(oneandrest[[5]])

plot1<- scatplotSLICE(A,0,0,"1")+theme(legend.position = "none")
plot2<- scatplotSLICE(B,0,0,"2")+theme(legend.position ="none")
plot3<-scatplotSLICE(C,0,0,"3")+theme(legend.position = "none")
plot4<-scatplotSLICE(D,0,0,"4")+theme(legend.position = "none")
plot5<-scatplotSLICE(E,0,0,"5")+theme(legend.position = "none")

grid.arrange(plot1, plot2, plot3,plot4,plot5, ncol=3,nrow=2)


scatplotSLICE<-function(dat, jitx,jity,lab){
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

  lightChain<-c('DILLTQSPAILSVSPGERVSFSCRASQFVGSSIHWYQQRTNGSPRLLIKYASESMSGIPSRFSGSGSGTDFTLSINTVESEDIADYYCQQSHSWPFTFGSGTNLEVKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
  lightChainC<-c('RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
  heavyChain<-c('EVKLEESGGGLVQPGGSMKLSCVASGFIFSNHWMNWVRQSPEKGLEWVAEIRSKSINSATHYAESVKGRFTISRDDSKSAVYLQMTDLRTEDTGVYYCSRNYYGSTYDYWGQGTTLTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLS')
  heavyChainFc<-c('ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK')
  lightch<-c(strsplit(lightChain,""))
  lightchC<-c(strsplit(lightChainC,""))
  heavych<-c(strsplit(heavyChain,""))
  heavychFc<-c(strsplit(heavyChainFc,""))
  lightch<-as.double(length(lightch[[1]]))
  lightchC<-as.double(length(lightchC[[1]]))
  heavych<-as.double(length(heavych[[1]]))
  heavychFc<-as.double(length(heavychFc[[1]]))
  
  TOP<-na.omit(TOP)
  for(i in 1:dim(TOP)[1]){
    if(TOP$V2[i] == ">light "){
      TOP$V1[i] = TOP$V1[i]+(heavych+1) 
    }
    if(TOP$V4[i] == ">light "){
      TOP$V3[i] = TOP$V3[i]+(heavych+1) 
    }
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
      print("swapped")
    }
  }
  splitdat<<-TOP
  jitter<-position_jitter(width = jitx,height=jity)

  #scatter plot os distaces
  TOP$col<-as.factor(TOP$col)
  ppp<-ggplot(TOP, aes(x=V1, y=V3, group=col))+
    geom_point(aes(colour=col),position =jitter, alpha=0.2,size=2)+
    scale_size(range = c(1,11 ), name="Distance (?)") +
    ylab("HC                        HV               LC            LV")+
    xlab("HC                        HV               LC            LV")+
    ggtitle(lab)+
    annotate(geom="rect",xmin = 120, xmax=451,
                  ymin = 541, ymax=661,
              color = "grey60",
              fill = "grey60",alpha = 0.3)+
    scale_x_continuous(expand = c(0,0), breaks = seq(1,(heavych+lightch+1), by =30), limits = c(0,(1+lightch+heavych+10)))+
    scale_y_continuous(expand = c(0,0), breaks= seq(1, (heavych+lightch+1), by =30), limits = c(0,(1+lightch+heavych+10)))+
    theme(axis.text.x = element_text(angle = 90, face = "bold",size = 7),
          axis.text.y = element_text(face = "bold", size = 7),
          #axis.title.y = element_text(size = 7, hjust = 1), #for big plot
          #axis.title.x = element_text(hjust = 1),
          axis.title.y = element_text(size = 6, hjust = 1, face = "bold"), #for grouped
          axis.title.x = element_text(hjust = 1, size = 7, face = "bold"),
          axis.text = element_text(
            color =  c(
              rep("grey1",length(seq(1,(heavychFc), by =30))),
              rep("grey40",length(seq(heavychFc,(heavych), by =30))),
              rep("grey15",length(seq(heavych,(lightchC+heavych+1), by =30))),
              rep("grey55",length(seq((lightchC+heavych+1),(lightch+heavych+1), by =30))))))
  return(ppp)
}

