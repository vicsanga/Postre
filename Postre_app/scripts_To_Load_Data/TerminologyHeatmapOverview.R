###########################################################################
# Terminology software
# Can be seen in ranking script, or masterSummary, and heatmapOverview script
###########################################################################

##Terminologia interna
#Long-range pathogenic effects
"LongRange"
#Gene Truncation
"Direct_geneTruncation"
#Gene Deletion
"Direct_geneDeletion"
#For Duplications
#If impact long-range (Neo-TAD)
"LongRange_geneDuplication"
#If impact by direct duplication
"Direct_geneDuplication"
#If impact by both, either geneDuplication or long-range
"Direct_LongRange_geneDuplication"
#If gene in the uncertainty region (space between coordinates of 1 breakpoint)
"Direct_uncertaintyRegion"

##BUT! For duplications and heatmap we may have various definitions along stages:
##So more than one definition of gene impact
##it is expected to occur for gene duplications if we predict sth by LongRange eg "LongRange_geneDuplication" or sth direct "Direct_geneDuplication"
##For the overall geneImpact summary column we put the broader concept. But for the stage specific analysis, we work with the more specific condition, for the reports
# Matrix_integratedResults[gene,"GeneImpact"]<-"Direct_LongRange_geneDuplication"


