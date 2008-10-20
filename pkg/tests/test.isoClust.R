#
# Test functionality of isoClust under different scenarios
###############################################################################

library(ORCME)
data(geneData)
data(doseData)
data(timeGeneData)
data(timeData)
data(esData)

##--  clustering for gene with monotone increasing profiles
incIso <- isoClust(expression=geneData,dose=doseData,alpha=1.2,lambda=0.05,phi=5,isoDir ='up',
    includeObserved=FALSE,doseResponse=TRUE)


##--  clustering for genes with monotone decreasing profiles
decIso <- isoClust(expression=geneData,dose=doseData,alpha=1.2,lambda=0.05,phi=5,isoDir ='down',
    includeObserved = FALSE,doseResponse=TRUE)

##--  clustering under both directions
Iso <- isoClust(expression=geneData,dose=doseData,alpha=1.2,lambda=0.05,phi=5,isoDir ='both',
    includeObserved = FALSE,doseResponse=TRUE)

################################################################################
################# clustering based on isotonic means and observed data  ########
################################################################################

##--  clustering for gene with monotone increasing profiles
incIsoObs <- isoClust(expression=geneData,dose=doseData,alpha=1.2,lambda=0.2,phi=5,isoDir ='up',
    includeObserved = TRUE,doseResponse=TRUE)
##--  clustering for genes with monotone decreasing profiles
decIsoObs <- isoClust(expression=geneData,dose=doseData,alpha=1.2,lambda=0.2,phi=5,isoDir ='down',
    includeObserved = TRUE,doseResponse=TRUE)
##--  clustering under both directions
IsoObs <- isoClust(expression=geneData,dose=doseData,alpha=1.2,lambda=0.2,phi=5,isoDir ='both',
    includeObserved = TRUE,doseResponse=TRUE)

############################################
### using ExpressionSet,character method ###
############################################

IsoObs <- isoClust(expression=esData,dose="dose",alpha=1.2,lambda=0.2,phi=5,isoDir ='both',
    includeObserved = TRUE,doseResponse=TRUE)
    

################################################################################
################# clustering of time trend Data  ########
################################################################################

###--clustering of time trend Data with replicates
timeobs <- isoClust(expression=timeGeneData,dose=timeData,alpha=1.2,lambda=0.2,phi=5,isoDir ='up',
    includeObserved = TRUE,doseResponse=FALSE)
    
### clustering of time trend Data without replicates
timeobs <- isoClust(expression=timeGeneData,dose=timeData,alpha=1.2,lambda=0.2,phi=5,isoDir ='up',
    includeObserved = FALSE,doseResponse=FALSE)

