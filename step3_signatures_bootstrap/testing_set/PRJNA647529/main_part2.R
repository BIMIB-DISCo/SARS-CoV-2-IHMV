# load required libraries
library("harmonicmeanp")

# set the seed
set.seed(8765933)

# load results
load("results/signatures_significance.RData")

# NOTE: 
# S1 --> C>A
# S2 --> T>C
# S3 --> C>T

# goodness of fit only significant signatures
print(fivenum(signatures_significance$goodness_fit))
print(round(length(which(signatures_significance$goodness_fit>0.90))/length(signatures_significance$goodness_fit),2))

# goodness of fit bootstrap
print(fivenum(signatures_significance$bootstrap_estimates$goodness_fit))
print(round(length(which(signatures_significance$bootstrap_estimates$goodness_fit>0.90))/length(signatures_significance$bootstrap_estimates$goodness_fit),2))

# consider pvalues
pvalues = signatures_significance$bootstrap_estimates$pvalues
pvalues_raw = pvalues
pvalues[,"S1"] = p.adjust(pvalues[,"S1"],method="fdr")
pvalues[,"S2"] = p.adjust(pvalues[,"S2"],method="fdr")
pvalues[,"S3"] = p.adjust(pvalues[,"S3"],method="fdr")

# S1 (C>A)
all_samples = length(pvalues[,"S1"])
with_S_samples = length(which(!is.na(pvalues[,"S1"])))
significant_S_samples = length(which(pvalues[which(!is.na(pvalues[,"S1"])),"S1"]<0.05))
print(all_samples)
print(with_S_samples)
print(significant_S_samples)
print(round(with_S_samples/all_samples,2))
print(round(significant_S_samples/all_samples,2))
print(round(significant_S_samples/with_S_samples,2))
print(p.hmp(as.numeric(pvalues_raw[which(!is.na(pvalues_raw[,"S1"])),"S1"]),L=all_samples))

# S2 (T>C)
all_samples = length(pvalues[,"S2"])
with_S_samples = length(which(!is.na(pvalues[,"S2"])))
significant_S_samples = length(which(pvalues[which(!is.na(pvalues[,"S2"])),"S2"]<0.05))
print(all_samples)
print(with_S_samples)
print(significant_S_samples)
print(round(with_S_samples/all_samples,2))
print(round(significant_S_samples/all_samples,2))
print(round(significant_S_samples/with_S_samples,2))
print(p.hmp(as.numeric(pvalues_raw[which(!is.na(pvalues_raw[,"S2"])),"S2"]),L=all_samples))

# S3 (C>T)
all_samples = length(pvalues[,"S3"])
with_S_samples = length(which(!is.na(pvalues[,"S3"])))
significant_S_samples = length(which(pvalues[which(!is.na(pvalues[,"S3"])),"S3"]<0.05))
print(all_samples)
print(with_S_samples)
print(significant_S_samples)
print(round(with_S_samples/all_samples,2))
print(round(significant_S_samples/all_samples,2))
print(round(significant_S_samples/with_S_samples,2))
print(p.hmp(as.numeric(pvalues[which(!is.na(pvalues[,"S3"])),"S3"]),L=all_samples))
