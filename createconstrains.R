#create constraints file using oneandrest data
ocre<-as.data.frame(table(oneandrest[[1]]))
ocre$Var1<-droplevels(ocre$Var1)
for(i in 1:dim(ocre)[1]){
  if(ocre$Freq[i] > 100000){
    ocre$Var2[i]<-str_c("1",ocre$Var1[i],sep = ",")
  }else{
    ocre$Var2[i]<-str_c("0",ocre$Var1[i],sep = ",")
  }
  
}
ocre$Var2<-str_replace_all(ocre$Var2,">heavy","AB")
ocre$Var2<-str_replace_all(ocre$Var2,">light","CD")

write.csv(ocre$Var2,file= "constraints112.csv",quote = F,row.names = F)
getwd()



write.csv2(ocre,file= "constraints111.csv")
