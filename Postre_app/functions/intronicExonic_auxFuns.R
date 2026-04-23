#############################################################################################
## Functions introduced in POSTRE logic as of Sept2025 to improve intronic-exonic handling
#############################################################################################

#In this script code for functions:
# isGeneContainingBreakpoint()
# bedIntersect()
# isVariantAlteringExonicSequence()

######################################################################
######################################################################

isGeneContainingBreakpoint<-function(gtf, chrBreakP, coordBreakP){
  ##Extracting/Parsing breakp coord info as done in: AffectedRegions.R
  breakpointCoords<-unique(as.integer(unlist(strsplit(coordBreakP,split = ",", fixed = TRUE))))
  target_BP_start<-min(breakpointCoords)
  target_BP_end<-max(breakpointCoords)
  
  #gtf only modified inside the function environment, so no need to copy it
  gtf$IsBPInSameChr<-gtf$chr==chrBreakP
  #Start is simply the left most coord regardless of whether is TSS or TTS
  gtf$IsBPstartContainedInside <- ((gtf$start <=  target_BP_start) & (gtf$end >= target_BP_start))
  gtf$IsBPendContainedInside <- ((gtf$start <=  target_BP_end) & (gtf$end >= target_BP_end))
  
  ##Even if just one of the breakpoint limits (start-end if uncertainty region) is contained inside of the gene
  #Report as broken
  filt_gtf<-subset(gtf,IsBPInSameChr==TRUE)
  genesContainingBreakpoint<-c(filt_gtf$geneSymbol[filt_gtf$IsBPstartContainedInside],
                               filt_gtf$geneSymbol[filt_gtf$IsBPendContainedInside])
  
  genesContainingBreakpoint<-unique(genesContainingBreakpoint)
  
  #There might be multiple if large uncertainty regions, or genes overlapping in both strands
  return(genesContainingBreakpoint)
}




#################################################################################################
## bedIntersect: Function defined to select those A regions that have an overlap with other B regions
## A regions can be extended Xkb so that region is broader
## Adaptation of overlappedRegions originally defined in Tomas Project (not as faster as genomic ranges more recent approaches but avoiding extra package installation, reducing thus the number of dependencies)
## When running in POSTRE A_region corresponds to SV coord, and B_regions the exons coords

## NOTE, to load the function to POSTRE environment
## I've added to Postre/Postre_app/scripts_To_Load_Data/metaFunctionLoad.R the path to this script 
#################################################################################################
bedIntersect<-function(A_regions, B_regions, kb_extension_A_regions){
  #print("Check that chr in both files present the same namings. That the characters are not factor codified and that positions are in numeric class.")
  #print("It is required a bed file with the region coordinates. With column names: chr, start, end")
  
  ###
  kb_extension_A_regions<-kb_extension_A_regions*1000
  
  selectedRegions<-integer() ##To track the positions of the maintained regions
  
  for(nRegion in 1:nrow(A_regions)){
    
    ##nRegion will also be used to select and track the mantained rows
    
    #to know the process status
    #we print each thousand multiple
    # if(nRegion %% 1000 == 0){
    #   cat("nStudiedREgions: ",nRegion,"\n")
    # }
    
    #######
    ##Take A region coordinates and extend them the kb distance defined
    targetRegionStart<-A_regions[nRegion,"start"] - kb_extension_A_regions
    targetRegionEnd<-A_regions[nRegion,"end"] + kb_extension_A_regions
    targetRegionChr<-A_regions[nRegion,"chr"]
    
    
    ###Check if any B_region overlapps with the extended A_regions ((targetRegionStart))
    ##If there is an overlap, either the start or the end of the B region (or both) must be between the limits of the A extended region
    B_regions$isChr<-B_regions$chr==targetRegionChr
    
    ##Start B region between A extended borders
    aux1<-B_regions$start >= targetRegionStart
    aux2<-B_regions$start <= targetRegionEnd
    auxMat<-rbind(aux1,aux2)
    B_regions$isStartBetweenBorders<-colSums(auxMat)==2
    
    ##End B region between A extended borders
    aux1<-B_regions$end >= targetRegionStart
    aux2<-B_regions$end <= targetRegionEnd
    auxMat<-rbind(aux1,aux2)
    B_regions$isEndBetweenBorders<-colSums(auxMat)==2
    
    ##It can also happen that B region contains A region (so there is also an intersection)
    aux1<-B_regions$end >= targetRegionEnd
    aux2<-B_regions$start <= targetRegionStart
    auxMat<-rbind(aux1,aux2)
    B_regions$isB_HoldingA<-colSums(auxMat)==2
    
    logicalInfo<-B_regions[,c("isChr","isStartBetweenBorders","isEndBetweenBorders","isB_HoldingA")]
    
    ##let's subset per chromosome
    ##because we are not gonna work now with the == 3, but with the >=2. 
    ##And it could be that the start and end are between borders, without having a real overlap
    ##due to being in different chromosomes
    logicalInfo<-subset(logicalInfo, isChr==TRUE)
    
    ##It can be that both limits of the B regions are inside the A extended region, or just one to have the overlap, or A is entirely included inside of B
    ##that is why >=2
    B_region_mapped<-any(rowSums(logicalInfo)>=2) ##if TRUE, at least there is an overlap between the region A and a region B
    
    if(B_region_mapped == TRUE){
      ##Remember this row, will be selected afterwards
      selectedRegions<-c(selectedRegions, nRegion)
    }
  }
  
  ##Selecting A regions with at least an overlap with a B region
  resultMatrix<-A_regions[selectedRegions,]
  return(resultMatrix)
}

#####################################################
#####################################################
isVariantAlteringExonicSequence<-function(gtf, targetGene, chrBreakP_1, coordBreakP_1, chrBreakP_2, coordBreakP_2){
  
  ##The objective is to identify if a SV is stricly occuring inside of one intron.
  #se HA DE EVALUAR PRIMERO SI LA VARIANTE (region comprendida entre LEFT MOST COORD BREAKP1 AND RIGHT MOST COORD BREAKP 2)
  #SI MAPEAN CON UN EXON, SI LA RESPUESTA ES SI, tenemos que altera la secuencia
  #SINO, ES PQ ES INTRONICA, NO HAY OTRA
  
  #PODRIA SER QUE LA ALTERACIÓN SEA UNA INVERSIÓN OCURRIENDO ESTRICTAMENTE EN INTRONES
  #E.G ENTRE INTRONES 1 Y 3 DE MANERA QUE SE INTERCAMBIAN LAS POSICIONES DE LOS EXONES (no tengo claro el grado de patogenicidad de esto)
  #Pero en cualquier caso ahi tambien se altera el producto resultante, no va a ser la misma secuencia en orden 5'-3' aunque reubiquemos tripletes enteros
  #SI ESTO OCCURE DE MOMENTO LO MARCO COMO ALGO TRUNCATING/DISRUPTING YA QUE HAY UNA ALTERACIÓN REALMENTE DE LA SECUENCIA EXONICA, IGUAL AFECTA AL PLEGAMIENTO PROTEINA (SI HAY MUCHO FALSO POSITIVO, REFINAR)
  #Posibilidad de refinar podria ser via evaluacion de los exones son simetricos (divisibles por 3 o no) ya que eso en ppio es menos probable que altere pautas de lectura. Y lo contrario.
  #"Only symmetrical exons can be duplicated in tandem or deleted without affecting the reading frame. Duplication or deletion of asymmetrical exons would disrupt the reading frame downstream."
  #"Symmetric exons are the only ones that can be inserted into introns, undergo duplication, or be deleted without changing the reading frame"
  
  ## More than half of Primate Specific (PS) exons (56.27%) have a length multiple of three, also called symmetric... It has been previously reported that conserved alternative exons present a bias towards symmetry [6,37,38]. As most of the PS exons are alternative, these numbers could just reflect a relationship between reading frame preservation and inclusion levels, 
  ## Al final no te interesa que al hacer un splicing alternativo cambies toda la secuencia de aminoacidos
  #  No obstante se ve como un numero no despreciable de exones tiene un numero que no es multiple de 3
  
  ##Y lo dicho en ultima instancia, si se sigue observando numero alto de FALSO POSITIVO, REFINAR.
  
  ##TENER EN CUENTA QUE HAY GENES COMO H3C15 QUE SOLO TIENEN UN EXON!!
  #POR TANTO SI NO HAY EXONES LA VARIANTE NO PUEDE SER STRICTLY INTRONIC
  
  # #Chr has to be the same for both, OTHERWISE RISE ERROR, since this code should not be reachable
  # #Skipping for now as there are different checkpoints upstream, to avoid excessive verifications
  # #Thus taking first one as reference
  # ##CONTROL DE MISMO CHR, won't happen so minimizing code executions
  # if(chrBreakP_1 != chrBreakP_2){
  #   stop("Not possible both bp should be in the same chr")
  # }
  
  ##PARSEAR COORDENADAS COMO TOCA CONVIRTIENDO DE STRING A NUMERIC Y PILLANDO LEFTMOST AND RIGHT MOST
  ##Extracting/Parsing breakp coord info as done in: AffectedRegions.R
  ##ESTO PENDIENTE DE AJUSTAR
  breakpoint_1_Coords<-unique(as.integer(unlist(strsplit(coordBreakP_1,split = ",", fixed = TRUE))))
  breakpoint_2_Coords<-unique(as.integer(unlist(strsplit(coordBreakP_2,split = ",", fixed = TRUE))))
  
  #Covering the widest possible range
  target_region_start<-min(breakpoint_1_Coords)
  target_region_end<-max(breakpoint_2_Coords)
  
  #Filtrar gtf para el gen
  filt_gtf<-subset(gtf, geneSymbol==targetGene)
  
  ##CREAR META_Exons bed para el transcrit (same strategy as in RefSeq annotation processing script)
  exonStarts<-as.numeric(unlist(strsplit(filt_gtf$exonStarts, split = ",")))
  exonEnds<-as.numeric(unlist(strsplit(filt_gtf$exonEnds, split = ",")))
  exonMat<-data.frame("chr"=filt_gtf$chr, "start"=exonStarts, "end"=exonEnds)
  
  #Hacer bedtools intersect
  bed_variantRegion<-data.frame("chr"=filt_gtf$chr, "start"=target_region_start, "end"=target_region_end)
  overlaps_variant_exons<-bedIntersect(A_regions = bed_variantRegion, B_regions = exonMat, kb_extension_A_regions = 0)#0bp extension, looking for direct overlap
  
  if(nrow(overlaps_variant_exons)==0){
    ##NO EXON DISRUPTION
    return(FALSE)
  }else if(nrow(overlaps_variant_exons)>0){
    #The variant overlaps with at least one exon (this is principle is only 1,either the variant overlaps or not)
    return(TRUE)
  }else{
    stop("ERROR, THIS POINT SHOULD NOT BE REACHABLE ON EXONIC ALTERATION EVALUATION")
  }

}
