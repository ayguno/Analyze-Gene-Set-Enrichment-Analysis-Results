
#This function is intended to summarize the results of a ssGSEA file
# After filtering by a desired FDR p-value, the program clusters the resulting normalized ssGSEA enrichment scores (NES)into a heatmap for visualization
# results (chr)= the Results.gct output from ssGSEA analysis containing NESs (converted to csv file)
# fdr (chr) = the fdr-pvalues.gct output from ssGSEA analysis containing NESs (converted to csv file)
# cutoff (num) = the desired fdr p-value cut off to be applied

# Last update 06/29/2016

results_ssGSEA<-function (results,fdr,cutoff) {
  
  require(gplots)
  require(RColorBrewer)
  
  NES<-read.csv(results,header=TRUE,skip=2) #Read the results.gct file   "CE1_CE2_CE3_Combined_.gct_Results.csv"
  FDR<-read.csv(fdr,header=TRUE,skip=2) #Read the fdr-pvalues.gct file    "CE1_CE2_CE3_Combined_.gct_Results-fdr-pvalues.csv"
  p<-cutoff
  
 
  
  tempcolnames<-substr(colnames(NES)[3:ncol(NES)],2,nchar(colnames(NES)[3:ncol(NES)]))
  colnames(NES)[3:ncol(NES)]<-tempcolnames
  colnames(FDR)[3:ncol(FDR)]<-tempcolnames # Droping the X in front of the relevant experiment column names
  
  NES<-NES[,-2] # No need to keep the description columns
  FDR<-FDR[,-2]
  
  
          w<-which((FDR[,2]<p & FDR[,3]<p) & NES[,2]*NES[,3]>0) #For the first experiment, get the gene set names that is reproducibly enriched or depleted AND reproducibly lower than the fdr cut off in both replicates
          
          FDRt<-FDR[w,]
          NESt<-NES[w,]
          
          write.csv(file=paste(p,"FDR_cutoff","FDR",colnames(FDR)[2],".csv",sep= "_"),FDRt) # Print the list for the experiment
          write.csv(file=paste(p,"FDR_cutoff","Results_NES",colnames(NES)[2],".csv",sep= "_"),NESt) # Print the list for the experiment
          
          Union_Name<-as.list(NULL)
          
          if(length(FDRt[,"Name"])> 0) 
                  
              {Union_Name<-as.list(as.character(FDRt[,"Name"]))}
          
          
            n<-ncol(FDR) 
            
            n1<-(n-1)/2 
            
            n2<-n1-1 
            
         #The next experiment column to be analyzed
            
            a<-4 
            
            for (i in 1:n2) { # Looping over the remaining individual experiments in the FDR and NES tables to print and collect union of gene set names
              
              w<-which((FDR[,a]<p & FDR[,a+1]<p) & NES[,a]*NES[,a+1]>0) #Get the gene set names that is reproducibly lower than the fdr cut off in both replicates
              
              FDRtemp<-FDR[w,]
              NEStemp<-NES[w,]
              
                    write.csv(file=paste(p,"FDR_cutoff","FDR",colnames(FDRtemp)[a],".csv",sep= "_"),FDRtemp) # Print the list for the next experiment
                    write.csv(file=paste(p,"FDR_cutoff","Results_NES",colnames(NEStemp)[a],".csv",sep= "_"),NEStemp) # Print the list for the next experiment
              
                    if(length(FDRtemp[,"Name"])> 0)    {Union_Name<-union(Union_Name,as.character(FDRtemp[,"Name"]))} #Get the concatenated union of the selected Gene list names for future access
              
              a<-a+2 # Move to the next experiment 
              
            } # continue until all the experiment columns 
  
           
  
                      ##Next prepare the union matrix and the heatmap
              
            if(length(Union_Name)>0) {
                    
                                      temp<-data.frame(unlist(Union_Name))     # convert the list into a data frame    
                                      colnames(temp)<-c("Name")   
                                      
                                      tempFDR<-temp  #The two matrices we will generate
                                      tempNES<-temp 
                      
                     
                      
                                       #  construct the union matrices
                                          
                                          tempFDR<-merge(tempFDR,FDR,all.x=T,by.x="Name",by.y="Name")
                                          tempNES<-merge(tempNES,NES,all.x=T,by.x="Name",by.y="Name")
                                        
                              
                                          write.csv(file=paste(p,"FDR_cutoff","Merged_Union_FDR",".csv",sep= "_"), tempFDR) # Print the two matrices for reference
                                          write.csv(file=paste(p,"FDR_cutoff","Merged_Union_Results_NES",".csv",sep= "_"),tempNES)
                              
                      
                      
                      
                      
                                              #Prepare the heatmaps
                                    
                                              xFDR<-as.matrix(tempFDR[2:ncol(tempFDR)],rownames.force = T)
                                              rownames(xFDR)<-tempFDR[,1]
                                              
                                              xNES<-as.matrix(tempNES[2:ncol(tempNES)],rownames.force = T)
                                              rownames(xNES)<-tempNES[,1]
                                              
                                                    pdf(file=paste(p,"FDR_pvalue_ssGSEA_clustered.pdf",sep="_"),width =11, height = 11)
                                                    
                                                    colfuncUPDN <- colorRampPalette(c("red","white")) 
                                                    heatmap.2(xFDR, Rowv = T, Colv = T, trace='n', 
                                                              dendrogram = 'both', notecol='black',
                                                              margins=c(13,15), density.info = 'none',col=colfuncUPDN(10), 
                                                              symbreaks=F,symkey = F,keysize = 1,
                                                              srtCol = 90, main=paste("FDR_pvalue<",p,sep=""),
                                                              key.title=F, key.xlab = "FDR_p_value",na.color = "grey",symm = F,
                                                              lwid=c(2,8),lhei = c(2,9) ) 
                                                    dev.off()
                                
                                                    pdf(file=paste(p,"NES_ssGSEA_clustered.pdf",sep="_"),width =11, height = 11)
                                                    colfuncUPDN <- colorRampPalette(c("blue","white","red")) 
                                                    heatmap.2(xNES, Rowv = T, Colv = F, trace='n', 
                                                              dendrogram = 'row', notecol='black',
                                                              margins=c(13,25), density.info = 'none',col=colfuncUPDN(1000), 
                                                              symbreaks=T,symkey = T,keysize = 1,
                                                              srtCol = 90,main=paste("FDR_pvalue<",p,sep=""),
                                                              key.title=F, key.xlab = "Normalized Enrichment Score",na.color = "grey",symm = F,
                                                              lwid=c(2,8),lhei = c(2,9) ) 
                                                    
                                                    dev.off()
                                                    
            }
            
           
              
            if(length(Union_Name)==0)
              
              {print("Error: analysis stopped! No gene sets were found based on the given FDR p value cut off")}
                      
                     
}