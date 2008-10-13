#
# Test functionality of isoClust under different scenarios
###############################################################################

library(ORCME)
data(geneData)
data(doseData)
data(esData)

##--  clustering for gene with monotone increasing profiles
incIso <- isoClust(expression=geneData,dose=doseData,alpha=1.2,lambda=0.05,phi=5,isoDir ='up',
    includeObserved=FALSE)


##--  clustering for genes with monotone decreasing profiles
decIso <- isoClust(expression=geneData,dose=doseData,alpha=1.2,lambda=0.05,phi=5,isoDir ='down',
    includeObserved = FALSE)

##--  clustering under both directions
Iso <- isoClust(expression=geneData,dose=doseData,alpha=1.2,lambda=0.05,phi=5,isoDir ='both',
    includeObserved = FALSE)

################################################################################
################# clustering based on isotonic means and observed data  ########
################################################################################

##--  clustering for gene with monotone increasing profiles
incIsoObs <- isoClust(expression=geneData,dose=doseData,alpha=1.2,lambda=0.2,phi=5,isoDir ='up',
    includeObserved = TRUE)
##--  clustering for genes with monotone decreasing profiles
decIsoObs <- isoClust(expression=geneData,dose=doseData,alpha=1.2,lambda=0.2,phi=5,isoDir ='down',
    includeObserved = TRUE)
##--  clustering under both directions
IsoObs <- isoClust(expression=geneData,dose=doseData,alpha=1.2,lambda=0.2,phi=5,isoDir ='both',
    includeObserved = TRUE)

############################################
### using ExpressionSet,character method ###
############################################

IsoObs <- isoClust(expression=esData,dose="dose",alpha=1.2,lambda=0.2,phi=5,isoDir ='both',
    includeObserved = TRUE)

