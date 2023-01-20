#################################
## Getting between TAD regions
#################################
getBetweenTADmap<-function(TADmap){
  
  ##Chr sizes, 
  ##So far only hg19 data since hg38 liftovered to hg19 so only load this
  ##Used to determine space between last TAD and chr end
  
  chrSizes<-read.delim(file = "data/hg19.chrom.sizes",
                       header = FALSE,
                       stringsAsFactors = FALSE,
                       sep="\t")
  
  colnames(chrSizes)<-c("chr","size")
  
  ##scaffold, used on each iteration to add new data before rbind with boundary map
  BtwnTADsScafold<-as.data.frame(matrix(data = NA, nrow = 1, ncol = 3))
  colnames(BtwnTADsScafold)<-c("chr","start","end")
  
  ##The final object to return
  BtwnTADs<-BtwnTADsScafold
  
  #hacer subset matriz por chromosoma
  for(studiedChr in unique(TADmap$chr)){
    subsetMapaTADchr<-subset(TADmap, chr==studiedChr)
    
    ContadorProceso<-0
    for(nfila in 1: (nrow(subsetMapaTADchr) + 1)) {
      ##El +1 de iteracion adicional es para, especificamente, pillar el hueco entre el ultimo TAD y el final del chromosoma
      
      ContadorProceso<-ContadorProceso + 1
      # print(ContadorProceso)
      # print(nfila)
      
      endTAD1<-(-1)##to initialize them
      startTAD2<-(-1)
      
      ##In case there is a space between chr start and first TAD
      if((nfila==1) && (subsetMapaTADchr[1,"start"]!=0)) {
        ##if nfila igual a 1 i la primera pos no es 0, ello indica que hay un hueco entre el inicio del chr y el primer TAD
        #cogemos ese hueco
        endTAD1<-0
        startTAD2<-subsetMapaTADchr[1,"start"]
        
      }else if((nfila > 1) && (nfila != (nrow(subsetMapaTADchr) + 1))){
        ##For the common scenario, space between TADs
        ##AND WE ARE NOT IN THE LAST ITERATION, to capture space between last TAD and chr end
        #empezamos a partir de la segunda fila
        endTAD1<-subsetMapaTADchr[nfila-1,"end"]
        startTAD2<-subsetMapaTADchr[nfila,"start"]
        
      }else if (nfila == (nrow(subsetMapaTADchr) + 1)){
        ##Last & extra iteration, capture space, if exists, between last TAD and chr end
        ##In case there is space between last TAD and chr end (most of it will be unknown)
        ##But get also this regions, because was rising error for some isolated and rare patient (U111 Cardiovascular)
        lastTAD<-subsetMapaTADchr[nrow(subsetMapaTADchr),]
        
        ##Get space between TAD end and chr end, if exists
        targetChrSize<-subset(chrSizes, chr == studiedChr)$size
        if(lastTAD$end < targetChrSize){
          ##There is a space between TAD end and chr end, capture it
          endTAD1<-lastTAD$end
          startTAD2<-targetChrSize
        }
      }
      
      ###############################################################################
      ## Storing between TADs - TAD external - regions
      ##If already found a region between TADs
      ##pq en la primera iteracion si el TAD empieza en 0 igual no encuentra space
      if(endTAD1!=-1 && startTAD2!=-1){
        ##generate limits
        ##We add the limit if the TADs are not 100%consecutive, 
        ##it means the next one start at the bp of the last one
        ##or if there is just one bp of diifference btwn the previous and the next one
        if((startTAD2-endTAD1)>1){
          BtwnTADsLimits<-BtwnTADsScafold
          BtwnTADsLimits[1,"chr"]<-studiedChr
          BtwnTADsLimits[1,"start"]<-endTAD1
          BtwnTADsLimits[1,"end"]<-startTAD2
          
          ##Anyadimos btwnTADlimits
          BtwnTADs<-rbind(BtwnTADs,BtwnTADsLimits)
          
        }
      }
    }
  }
  ##remove first row used to initialize the TAD map matrix
  BtwnTADs<-BtwnTADs[-1,]
  
  ##Returning Between TAD coordinates
  return(BtwnTADs)
}