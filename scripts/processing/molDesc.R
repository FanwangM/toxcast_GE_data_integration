
## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved 
## 18/01/15
## Compute molecular descriptors and analyze chemical space of dataset

library(plyr)
library(Rcpi)
library(rcdk)
CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'

#chemicals
load(file.path(CALC_DATA_DIR,"/chemicals.Rdata"))
writeLines(chemicals$SMILES,file.path(CALC_DATA_DIR,'/SmilesOnly.smi'),sep='\n')

#convert smiles to mol objects
x.mol = readMolFromSmi(file.path(CALC_DATA_DIR,'/SmilesOnly.smi'),type='mol')
x.smi = readMolFromSmi(file.path(CALC_DATA_DIR,'/SmilesOnly.smi'),type='text')

# modify mol descriptor computation function
# so that errors are handled properly and
# dailing descriptors ignored 
extractDrugAIO = function (molecules, silent = TRUE, warn = TRUE) 
{
    if (warn == TRUE) {
        warning("Note that we need 3-D coordinates of the molecules to calculate some of the descriptors, if not provided, these descriptors will be NA")
    }
    descNames = c("org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.AminoAcidCountDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.AromaticBondsCountDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorCharge", 
        "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorMass", 
        "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorPolarizability", 
        "org.openscience.cdk.qsar.descriptors.molecular.BCUTDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.BPolDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.BondCountDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.CPSADescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.ChiChainDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.ChiClusterDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.ChiPathClusterDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.ChiPathDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.EccentricConnectivityIndexDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.FMFDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.FragmentComplexityDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.GravitationalIndexDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor", 
        "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.HybridizationRatioDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.IPMolecularLearningDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.KappaShapeIndicesDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.KierHallSmartsDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.LargestChainDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.LargestPiSystemDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.LengthOverBreadthDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.MDEDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.MannholdLogPDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.MomentOfInertiaDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.PetitjeanShapeIndexDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.VABCDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.VAdjMaDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.WHIMDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.WeightedPathDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor",
        "org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor")
    x = llply(descNames, function(molD) return(try(eval.desc(molecules, molD, verbose = !silent))))
    return(x)
} 


#extract all molecular descriptors
molDesc=extractDrugAIO(x.mol, silent = T, warn = F)
molDesc = molDesc[which(llply(molDesc,class) != 'try-error')]

# bind these in a unique data frame
md.df = do.call('cbind',molDesc)

# remove molecular descriptors with only NAs
nNAs = unlist(alply(md.df, 2, function(df) length(which(is.na(df)))))
md.df = md.df[,which(nNAs < nrow(md.df))]

# replace smiles by their corresponding lincs id
rownames(md.df) = chemicals[match(rownames(md.df),chemicals$SMILES),'LINCS_ID']

save(md.df,file=file.path(CALC_DATA_DIR,'datasets/molDescriptors.Rdata'))
