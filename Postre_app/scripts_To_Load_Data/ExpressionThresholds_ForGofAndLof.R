################################################
##Thresholds de expression para LOF y GOFs
################################################

threshold_MaxExpresion<-10##5##10##If it is above this threshold, considerably expressed gene
threshold_MinExpresion<-1##If a gene expresion below this threshold, not considerably expressed gene

### For ratios Enhancer Balance
# maxRatioEnhBalance<-1.75#2
##Calculado con regla de 3 inversa para el LOF makes sense, ya que 1 -> tener 100%, 2 -> 50%
maxRatioEnhBalance_LOF<-1.11 ##10% loss, pathogenic, as seen in PRS patient with SOX9 enh deletion. Dejar aqui a ver, y si no se disparan predicciones, dejar aqui
maxRatioEnhBalance_GOF<-2.5 ##For GOF I'm going to leave it to 1.5 which implies 50% gain, and is already much more permissive than positive controls, which go for 700% or from nothing to X (so infinite gain) this way we can expect intra tad Dup to have a chance to be pathogenic by upregulationo of intraTAD gene (as with SOX9 duplication of enhancers or SHH duplications), eventhough we do not have patients of that kind so far 

minRatioEnhBalance<-1