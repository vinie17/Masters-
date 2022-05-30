
# Libraries
#library(tidyverse)
#library(viridis)
#library(patchwork)
#library(hrbrthemes)
#library(igraph)
#library(ggraph)
#library(colormap)


require(gridExtra)
plot1 <- scatplot(oneandrest[[1]],1,1,"Lower band unique")+theme(legend.position = "none")
plot2 <- scatplot(oneandrest[[2]],1,1,"Intersect")+theme(legend.position ="none")
plot3<-scatplot(oneandrest[[3]],1,1,"Upper band unique")+theme(legend.position = "none")
grid.arrange(plot1, plot2, plot3, ncol=3)

#outersect
out<-Map(c,ABC[[4]],ABC[[5]],ABC[[6]])
out<-list(out)
#intersect
intersects<-list(trip[[2]])


data<-ldply(oneandrest[1:3],data.frame)
data<-list(data)

A<-list(oneandrest[[1]])
B<-list(oneandrest[[2]])
C<-list(oneandrest[[3]])
D<-list(oneandrest[[4]])
E<-list(oneandrest[[5]])
FF<-list(oneandrest[[6]])
G<-list(oneandrest[[7]])

plot1<-scatplotSLICE_CD163(A,0,0,linklist[1])+theme(legend.position = "none")
plot2<-scatplotSLICE_CD163(B,0,0,linklist[2])+theme(legend.position = "none")
plot3<-scatplotSLICE_CD163(G,0,0,linklist[7])+theme(legend.position = "none")








plot1<-scatplotSLICE(BOT,0,0,linklist[1])+theme(legend.position = "none")
plot2<- scatplotSLICE(MID,0,0,"MID")+theme(legend.position ="none")
plot3<-scatplotSLICE(TOPP,0,0,linklist[2])+theme(legend.position = "none")

plot1<-scatplotSLICE(A,0,0,linklist[1])+theme(legend.position ="none")
plot2<- scatplotSLICE(B,0,0,linklist[2])+theme(legend.position ="none")
plot3<-scatplotSLICE(C,0,0,linklist[3])+theme(legend.position = "none")
plot4<-scatplotSLICE(D,0,0,linklist[4])+theme(legend.position = "none")
plot5<-scatplotSLICE(E,0,0,linklist[5])+theme(legend.position = "none")
plot6<-scatplotSLICE(FF,0,0,linklist[6])+theme(legend.position = "none")
plot7<-scatplotSLICE(G,0,0,linklist[7])+theme(legend.position = "none")



grid.arrange(plot1,plot2 ,ncol = 2)
grid.arrange(plot3, plot6, ncol=2)

#scatplot(oneandrest[[6]],5,5,"Ocrelizumab light+heavy")

scatplotSLICE<-function(dat, jitx,jity,lab){
  #dat<-A
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
  #e10b10
  #lightChain<-c('DIVMTQTPSSQAVSAGERVTMRCKSSQSLLYSENKKNYLAWYQQKPGQSPKLLISWASTRESGVPDRFIGSGSGTDFTLTISSVQAEDLAVYYCDQYYDPPFTFGSGTKLEIKRADAAPTVSIFPPSTEQLATGGASVVCLMNNFYPRDISVKWKIDGTERRDGVLDSVTDQDSKDSTYSMSSTLSLTKADYESHNLYTCEVVHKTSSSPVVKSFNRNEC')
  #heavyChain<-c('QVKLQESGGGLVQPKESLKISCAASGFTFSTAAMYWVRQAPGKGLDWVARIRTKPDNYATYYPASVKGRFTISRDDSKGMVYLQMDNLKTEDTAIYYCTAAYYYDGRFDYWGQGVMVTVASAETTPKLVYPLAPGKHSTALKSNSMVTLGCLVKGYFPEPVTVTWNSGALSSGVHTFPAVLQSGLYTLTSSVTVPSSTWSSQAVTCNVAHPASSTKVDKKIVPRECNPCGCTGSEVSSVFIFPPKTKDVLTITLTPKVTCVVVDISQNDPEVRFSWFIDDVEVHTAQTHAPEKQSNSTLRSVSELPIVHRDWLNGKTFKCKVNSGAFPAPIEKSISKPEGTPRGPQVYTMAPPKEEMTQSQVSITCMVKGFYPPDIYTEWKMNGQPQENYKNTPPTMDTDGSYFLYSKLNVKKETWQQGNTFTCSVLHEGLHNHHTEKSLSHSPGK')
  #tys
  #lightChain<-c('DIQMTQSPSSLSASVGDRVTITCKTSQDINKYMAWYQQTPGKAPRLLIHYTSALQPGIPSRFSGSGSGRDYTFTISSLQPEDIATYYCLQYDNLWTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
  #heavyChain<-c('QVQLVQSGAEVKKPGASVKVSCKASGFNIKDTYIHWVRQAPGQRLEWMGRIDPANGYTKYDPKFQGRVTITADTSASTAYMELSSLRSEDTAVYYCAREGYYGNYGVYAMDYWGQGTLVTVSSASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTKTYTCNVDHKPSNTKVDKRVESKYGPPCPSCPAPEFLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPEVQFNWYVDGVEVHNAKTKPREEQFNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKGLPSSIEKTISKAKGQPREPQVYTLPPSQEEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSRLTVDKSRWQEGNVFSCSVMHEALHNHYTQKSLSLSLGK')
  #OCRE
  #heavyChain<-c('EVQLVESGGGLVQPGGSLRLSCAASGYTFTSYNMHWVRQAPGKGLEWVGAIYPGNGDTSYNQKFKGRFTISVDKSKNTLYLQMNSLRAEDTAVYYCARVVYYSNSYWYFDVWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK')
  #lightChain<-c('DIQMTQSPSSLSASVGDRVTITCRASSSVSYMHWYQQKPGKAPKPLIYAPSNLASGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQWSFNPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
  #INF
  lightChain<-c('DILLTQSPAILSVSPGERVSFSCRASQFVGSSIHWYQQRTNGSPRLLIKYASESMSGIPSRFSGSGSGTDFTLSINTVESEDIADYYCQQSHSWPFTFGSGTNLEVKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
  lightChainC<-c('RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
  heavyChain<-c('EVKLEESGGGLVQPGGSMKLSCVASGFIFSNHWMNWVRQSPEKGLEWVAEIRSKSINSATHYAESVKGRFTISRDDSKSAVYLQMTDLRTEDTGVYYCSRNYYGSTYDYWGQGTTLTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLS')
  heavyChainFc<-c('ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK')
  #iga
  #heavyChain<-c('QVKLLEQSGAEVKKPGASVKVSCKASGYSFTSYGLHWVRQAPGQRLEWMGWISAGTGNTKYSQKFRGRVTFTRDTSATTAYMGLSSLRPEDTAVYYCARDPYGGGKSEFDYWGQGTLVTVSSASPTSPKVFPLSLCSTQPDGNVVIACLVQGFFPQEPLSVTWSESGQGVTARNFPPSQDASGDLYTTSSQLTLPATQCLAGKSVTCHVKHYTNPSQDVTVPCPVPSTPPTPSPSTPPTPSPSCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGVTFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAEPWNHGKTFTCTAAYPESKTPLTATLSKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTFAVTSILRVAAEDWKKGDTFSCMVGHEALPLAFTQKTIDRLAGKPTHVNVSVVMAEVDGTCY')
  #lightChain<-c('ELVMTQSPSSLSASVGDRVNIACRASQGISSALAWYQQKPGKAPRLLIYDASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFAIYYCQQFNSYPLTFGGGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
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
    if(TOP$V2[i] == ">E10B10-L"){
      TOP$V1[i] = TOP$V1[i]+(heavych+1) 
    }
    if(TOP$V4[i] == ">E10B10-L"){
      TOP$V3[i] = TOP$V3[i]+(heavych+1) 
    }
    if(TOP$V2[i] == ">light "){
      TOP$V1[i] = TOP$V1[i]+(heavych+1) 
    }
    if(TOP$V4[i] == ">light "){
      TOP$V3[i] = TOP$V3[i]+(heavych+1) 
    }
    if(TOP$V2[i] == ">light"){
      TOP$V1[i] = TOP$V1[i]+(heavych+1) 
    }
    if(TOP$V4[i] == ">light"){
      TOP$V3[i] = TOP$V3[i]+(heavych+1) 
    }
    
    if(TOP$V2[i] == TOP$V4[i]){
      if(TOP$V1[i]==TOP$V3[i]){
        TOP$col[i]= "Homeotypic"
        TOP$V1[i]= TOP$V1[i]+0.1
      }
      else{
        TOP$col[i] = "Type 1"
      }
    }
    else{TOP$col[i]="Type 2"}
    
  }
  
  for(i in 1:dim(TOP)[1]){
    if(TOP$V3[i]<TOP$V1[i]){
      one<-TOP$V1[i]
      TOP$V1[i]<-TOP$V3[i]
      TOP$V3[i]<-one
      #print("swapped")
    }
  }
  splitdat<<-TOP
  jitter<-position_jitter(width = jitx,height=jity)
  
  #scatter plot os distaces
  TOP$col<-as.factor(TOP$col)
  ppp<-ggplot(TOP, aes(x=V1, y=V3, group=col))+
    geom_point(aes(colour=col),position =jitter, alpha=1,size=2)+
    scale_size(range = c(1,11 ), name="Distance (?)") +
    ylab("Hv                                                                                         Hc                                                         Lv                                      Lc        ")+
    xlab("Hv                                                             Hc                                                   Lv                            Lc  ")+
    ggtitle(lab)+
    scale_color_manual(values = setNames(c('#619CFF',"#00BA38","#F8766D"),c("Type 1","Type 2","Homeotypic")))+
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
  return(ppp)
}

scatplotSLICE_CD163<-function(dat, jitx,jity,lab){
  #dat<-A
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
  #e10b10
  #lightChain<-c('DIVMTQTPSSQAVSAGERVTMRCKSSQSLLYSENKKNYLAWYQQKPGQSPKLLISWASTRESGVPDRFIGSGSGTDFTLTISSVQAEDLAVYYCDQYYDPPFTFGSGTKLEIKRADAAPTVSIFPPSTEQLATGGASVVCLMNNFYPRDISVKWKIDGTERRDGVLDSVTDQDSKDSTYSMSSTLSLTKADYESHNLYTCEVVHKTSSSPVVKSFNRNEC')
  #heavyChain<-c('QVKLQESGGGLVQPKESLKISCAASGFTFSTAAMYWVRQAPGKGLDWVARIRTKPDNYATYYPASVKGRFTISRDDSKGMVYLQMDNLKTEDTAIYYCTAAYYYDGRFDYWGQGVMVTVASAETTPKLVYPLAPGKHSTALKSNSMVTLGCLVKGYFPEPVTVTWNSGALSSGVHTFPAVLQSGLYTLTSSVTVPSSTWSSQAVTCNVAHPASSTKVDKKIVPRECNPCGCTGSEVSSVFIFPPKTKDVLTITLTPKVTCVVVDISQNDPEVRFSWFIDDVEVHTAQTHAPEKQSNSTLRSVSELPIVHRDWLNGKTFKCKVNSGAFPAPIEKSISKPEGTPRGPQVYTMAPPKEEMTQSQVSITCMVKGFYPPDIYTEWKMNGQPQENYKNTPPTMDTDGSYFLYSKLNVKKETWQQGNTFTCSVLHEGLHNHHTEKSLSHSPGK')
  #tys
  #lightChain<-c('DIQMTQSPSSLSASVGDRVTITCKTSQDINKYMAWYQQTPGKAPRLLIHYTSALQPGIPSRFSGSGSGRDYTFTISSLQPEDIATYYCLQYDNLWTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
  #heavyChain<-c('QVQLVQSGAEVKKPGASVKVSCKASGFNIKDTYIHWVRQAPGQRLEWMGRIDPANGYTKYDPKFQGRVTITADTSASTAYMELSSLRSEDTAVYYCAREGYYGNYGVYAMDYWGQGTLVTVSSASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTKTYTCNVDHKPSNTKVDKRVESKYGPPCPSCPAPEFLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPEVQFNWYVDGVEVHNAKTKPREEQFNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKGLPSSIEKTISKAKGQPREPQVYTLPPSQEEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSRLTVDKSRWQEGNVFSCSVMHEALHNHYTQKSLSLSLGK')
  #OCRE
  #heavyChain<-c('EVQLVESGGGLVQPGGSLRLSCAASGYTFTSYNMHWVRQAPGKGLEWVGAIYPGNGDTSYNQKFKGRFTISVDKSKNTLYLQMNSLRAEDTAVYYCARVVYYSNSYWYFDVWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK')
  #lightChain<-c('DIQMTQSPSSLSASVGDRVTITCRASSSVSYMHWYQQKPGKAPKPLIYAPSNLASGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQWSFNPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC')
  #INF
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
    if(TOP$V2[i] == ">E10B10-L"){
      TOP$V1[i] = TOP$V1[i]+(heavych+1) 
    }
    if(TOP$V4[i] == ">E10B10-L"){
      TOP$V3[i] = TOP$V3[i]+(heavych+1) 
    }
    if(TOP$V2[i] == "CD163"){
      
      TOP$V1[i] = (TOP$V1[i]+(heavych+lightch+lightchC+1)-44) 
      print(TOP$V1[i])
    
      }
    if(TOP$V4[i] == "CD163"){
      
      TOP$V3[i] = (TOP$V3[i]+(heavych+lightch+lightchC+1)) 
      print(TOP$V3[i])
       }
    
    
    if(TOP$V2[i] == TOP$V4[i]){
      if(TOP$V1[i]==TOP$V3[i]){
        TOP$col[i]= "Homeotypic"
        TOP$V1[i]= TOP$V1[i]+0.1
      }
      else{
        TOP$col[i] = "Type 1"
      }
    }
    else{TOP$col[i]="Type 2"}
    if(TOP$V2[i] == "CD163"){
      TOP$col[i]= "CD163"}
    if(TOP$V4[i] == "CD163"){
      TOP$col[i]= "CD163"
    }
    
  }
  
  for(i in 1:dim(TOP)[1]){
    if(TOP$V3[i]<TOP$V1[i]){
      one<-TOP$V1[i]
      TOP$V1[i]<-TOP$V3[i]
      TOP$V3[i]<-one
      #print("swapped")
    }
  }
  splitdat<<-TOP
  jitter<-position_jitter(width = jitx,height=jity)
  
  #scatter plot os distaces
  TOP$col<-as.factor(TOP$col)
  ppp<-ggplot(TOP, aes(x=V1, y=V3, group=col))+
    geom_point(aes(colour=col),position =jitter, alpha=1,size=2)+
    scale_size(range = c(1,11 ), name="Distance (?)") +
    ylab("Hv                                                                                         Hc                                                         Lv                                      Lc        ")+
    xlab("Hv                                                             Hc                                                   Lv                            Lc  ")+
    ggtitle(lab)+
    scale_color_manual(values = setNames(c('#619CFF',"#00BA38","#F8766D","orange"),c("Type 1","Type 2","Homeotypic","CD163")))+
    annotate(geom="rect",xmin = 300, xmax=451,
             ymin = 460, ymax=630,
             color = NA,
             fill = "grey60",alpha = 0.3)+
    scale_x_continuous(expand = c(0,0), breaks = seq(0,(heavych+lightch+700), by =30), limits = c(0,(1+lightch+heavych+700)))+
    scale_y_continuous(expand = c(0,0), breaks= seq(0, (heavych+lightch+700), by =30), limits = c(0,(1+lightch+heavych+700)))+
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
              rep("grey55",length(seq(((lightch-lightchC)+heavych),(lightch+heavych+1+49), by =30))),
              rep("grey20",length(seq((lightch+heavych+1),(heavych+lightch+lightchC+700), by =30))))))
  return(ppp)
}
