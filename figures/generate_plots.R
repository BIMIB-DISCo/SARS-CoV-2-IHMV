##### USA DATASET ANALYSIS 
##### 1133 SAMPLES 
##################################

# To run the following code you need to install these CRAN packages.

# install.package('dplyr')
# install.package('ggplot2')
# install.package('ggforce')
# install.package('gridExtra')
# install.package('pheatmap')
# install.package('RColorBrewer')
# install.package('seqinr')
# install.package('viridis')

#RESET WORKSPACE
rm(list = ls())

library('seqinr')

unzip('data.zip',exdir='data')

#LOAD DATA
processed_variants <- readRDS("data/processed_variants.rds")
statistics <- readRDS("data/statistics.rds")
signatures_clustering <- readRDS("data/signatures_clustering.rds")
assignments <- readRDS("data/assignments.rds")

reference_Seq <- read.table(file = "data/REF_ANC_annotation.txt", header = T, sep = '\t', stringsAsFactors = F)
fileID <- file("data/REF-ANC.fasta", "r")
REFANC <- readLines(fileID)
close(fileID)
REFANC <- unlist(strsplit(REFANC[-1], ''))

if(!dir.exists('plots')){
  dir.create('plots')
}

#################
# PREPROCESSING #
#################

#SET-UP PARAMETERS
thrCLONAL <- 0.90
thrDETECTION <- 0.05
thrTRANS <- 0.20

processed_variants_na <- processed_variants

processed_variants[which(is.na(processed_variants))] <- 0


variant_position <- sapply(colnames(processed_variants), 
                           FUN = function(x) {strsplit(x, '_')[[1]][1]}, 
                           USE.NAMES = F)

# get mutation with more than one alternative allele
more_than_one_alt_alleles <- colnames(processed_variants)[!ave(variant_position, 
                                                               variant_position, FUN = length) == 1]
# Compute outliers and filtering
p_nMinor <- apply(processed_variants, MARGIN = 1, FUN = function(x){sum(x > thrDETECTION & x <= thrCLONAL, na.rm = T)})
p_nClonal <- apply(processed_variants, MARGIN = 1, FUN = function(x){sum(x > thrCLONAL, na.rm = T)})


var_freq_matrix_filt <- processed_variants

var_freq_matrix_filt <- var_freq_matrix_filt[,colSums(var_freq_matrix_filt, na.rm = T) > thrDETECTION]

processed_variants_na_filt <- processed_variants_na[match(rownames(var_freq_matrix_filt), rownames(processed_variants_na)),
                                                    match(colnames(var_freq_matrix_filt), colnames(processed_variants_na))]

p_nMinor_outFilt <- apply(var_freq_matrix_filt, MARGIN = 1, FUN = function(x){sum(x>=thrDETECTION & x<=thrCLONAL, na.rm = T)})

p_nClonal_outFilt <- apply(var_freq_matrix_filt, MARGIN = 1, FUN = function(x){sum(x>thrCLONAL, na.rm = T)})

p_VAFmean_Minor_outFilt <- apply(var_freq_matrix_filt, MARGIN = 1, FUN = function(x){mean(x[which(x>=thrDETECTION & x<=thrCLONAL, !is.na(x))], na.rm = T)})
p_VAFmean_Minor_outFilt[which(p_nMinor_outFilt==0)] <- 0 

p_nMinor_outFilt <- as.data.frame(p_nMinor_outFilt)
p_nMinor_outFilt$Run <- rownames(p_nMinor_outFilt)
colnames(p_nMinor_outFilt) <- c("nMinor", "Run")

p_nClonal_outFilt <- as.data.frame(p_nClonal_outFilt)
p_nClonal_outFilt$Run <- rownames(p_nClonal_outFilt)
colnames(p_nClonal_outFilt) <- c("nClonal", "Run")

p_VAFmean_Minor_outFilt <- as.data.frame(p_VAFmean_Minor_outFilt)
p_VAFmean_Minor_outFilt$Run <- rownames(p_VAFmean_Minor_outFilt)
colnames(p_VAFmean_Minor_outFilt) <- c("meanVF", "Run")

statistics$CollectionDate <- as.Date(as.character(statistics$CollectionDate), 
                                     format = "%Y-%m-%d")

statistics <- merge(x = statistics, y = p_nMinor_outFilt, by = "Run")
statistics <- merge(x = statistics, y = p_nClonal_outFilt, by = "Run")
statistics <- merge(x = statistics, y = p_VAFmean_Minor_outFilt, by = "Run")

# Get week number
statistics$Week <- as.numeric(strftime(statistics$CollectionDate, format = "%V"))
statistics$WeekMerged <- NA
statistics$WeekMerged[statistics$Week < 16] <- "<16"

for(i in seq(16, 38, by = 3)) {
  statistics$WeekMerged[statistics$Week >= i & statistics$Week <= (i+2)] <- paste0(i, '-', (i+2))
}

statistics$UMAPgroup <- 1

#Add signature information
clusters_signatures <- as.data.frame(signatures_clustering)
clusters_signatures$Run <- rownames(clusters_signatures)
colnames(clusters_signatures) <- c("SignatureCluster", "Run")
clusters_signatures$SignatureCluster <- gsub(pattern = 'Cluster ', replacement = '', x =  clusters_signatures$SignatureCluster)

statistics <- merge(x = statistics, y = clusters_signatures, by = "Run", all.x = T)
statistics$SignatureCluster[which(is.na(statistics$SignatureCluster))] <- 'SC#4'
statistics$SignatureCluster[statistics$nMinor == 0] <- NA

signature_label <- c("SC#1", "SC#2","SC#3", "SC#4")
statistics$SignatureCluster <- factor(statistics$SignatureCluster, levels = signature_label)

#Add VERSO cluster
clades <- as.data.frame(assignments)
clades$Run <-rownames(clades)
colnames(clades) <- c('VERSO', 'Run')

statistics <- merge(x = statistics, y = clades, by = 'Run')

SignatureColor = c("#d4a155","#738762","#a14d3f","#acaeb8")
names(SignatureColor) <- c("SC#1","SC#2","SC#3","SC#4")




## TESTING DISTRIBUTION OF NUMBER CLONAL OR MINOR SNV AMONG CLUSTERS

nClonalMinor_signature <- split(statistics[,c("nMinor", "nClonal", "SignatureCluster")], statistics$SignatureCluster)

minor_significance <- matrix(data = NA, nrow = 4, ncol = 4, dimnames = list(c('SC#1','SC#2','SC#3','SC#4'), c('SC#1','SC#2','SC#3','SC#4')))
clonal_significance <- matrix(data = NA, nrow = 4, ncol = 4, dimnames = list(c('SC#1','SC#2','SC#3','SC#4'), c('SC#1','SC#2','SC#3','SC#4')))

for(i in 1:length(nClonalMinor_signature)){
  for(j in 1:length(nClonalMinor_signature)){
    
    minor_significance[i,j] <- ks.test(x = unlist(nClonalMinor_signature[[i]][1]), y = unlist(nClonalMinor_signature[[j]][1]))$p.value
    clonal_significance[i,j] <- ks.test(x = unlist(nClonalMinor_signature[[i]][2]), y = unlist(nClonalMinor_signature[[j]][2]))$p.value
  }
}
print(clonal_significance)
print(minor_significance)

ggplot2::ggplot(data = statistics[!is.na(statistics$SignatureCluster) & statistics$SignatureCluster %in% c('SC#3','SC#4'), ], ggplot2::aes(x = nClonal, fill = SignatureCluster)) + 
  ggplot2::geom_histogram(ggplot2::aes(y = ..density..), alpha = 0.4, binwidth = 0.5, position = "identity")

##########################
# MUTATION TYPE METADATA #
##########################
nPatients <- 1

clonal_matrix <- var_freq_matrix_filt > thrCLONAL
minor_matrix <- ((var_freq_matrix_filt <= thrCLONAL) & (var_freq_matrix_filt > thrDETECTION))
low_matrix <- ((var_freq_matrix_filt <= thrTRANS) & (var_freq_matrix_filt > thrDETECTION))

MUTATIONdf <- data.frame(MuID = colnames(var_freq_matrix_filt))
MUTATIONdf$type <- NA
MUTATIONdf$TrueTrans <- 'no'
MUTATIONdf$HighMinor <- 'no'

CLONALvar <- colnames(var_freq_matrix_filt)[colSums(clonal_matrix, na.rm = T)>0]
MINORvar <- colnames(var_freq_matrix_filt)[colSums(minor_matrix, na.rm = T)>0]
LOWvar <- colnames(var_freq_matrix_filt)[colSums(low_matrix, na.rm = T)>=nPatients]

TRANSvar <- intersect(CLONALvar, MINORvar)
NOTLOWvar <- setdiff(TRANSvar, LOWvar)

MUTATIONdf$type[MUTATIONdf$MuID %in% MINORvar] <- "always minor"

MUTATIONdf$type[MUTATIONdf$MuID %in% CLONALvar] <- "always clonal"

MUTATIONdf$type[MUTATIONdf$MuID %in% TRANSvar] <- "mixed"

MUTATIONdf <- MUTATIONdf[!is.na(MUTATIONdf$type),]

MUTATIONdf$type <- factor(MUTATIONdf$type, levels = c('always clonal', 'mixed', 'always minor'))

MUTATIONdf$TrueTrans[MUTATIONdf$MuID %in% NOTLOWvar] <- "yes"

MUTATIONdf$HighMinor[!(MUTATIONdf$MuID %in% LOWvar)] <- "yes"

##############################
# PATIENT VARIATION METADATA #
##############################

matrixCLONAL <- var_freq_matrix_filt[, colnames(var_freq_matrix_filt) %in% MUTATIONdf$MuID[MUTATIONdf$type == "always clonal"]]
CLONALdf <- reshape2::melt(matrixCLONAL)
colnames(CLONALdf) <- c("Sample", "MuID", "VF")
CLONALdf <- CLONALdf[CLONALdf$VF > thrDETECTION,]
CLONALdf$type <- "always clonal"

matrixTRANSITION <- var_freq_matrix_filt[, colnames(var_freq_matrix_filt) %in% MUTATIONdf$MuID[MUTATIONdf$type == "mixed"]]
TRANSITIONdf <- reshape2::melt(matrixTRANSITION)
colnames(TRANSITIONdf) <- c("Sample", "MuID", "VF")
TRANSITIONdf <- TRANSITIONdf[TRANSITIONdf$VF > thrDETECTION,]
TRANSITIONdf$type <- "mixed"

matrixMINOR <- var_freq_matrix_filt[, colnames(var_freq_matrix_filt) %in% MUTATIONdf$MuID[MUTATIONdf$type == "always minor"]]
MINORdf <- reshape2::melt(matrixMINOR)
colnames(MINORdf) <- c("Sample", "MuID", "VF")
MINORdf <- MINORdf[MINORdf$VF > thrDETECTION,]
MINORdf$type <- "always minor"

VARIANTdf <- rbind.data.frame(CLONALdf,TRANSITIONdf,MINORdf)
VARIANTdf$type <- factor(VARIANTdf$type, levels = c("always clonal", "mixed", "always minor"))
VARIANTdf$pos <- sapply(as.character(VARIANTdf$MuID), FUN = function(x){
  strsplit(x, split = '_')[[1]][1]
})
VARIANTdf$ref <- sapply(as.character(VARIANTdf$MuID), FUN = function(x){
  strsplit(x, split = '_')[[1]][2]
})
VARIANTdf$alt <- sapply(as.character(VARIANTdf$MuID), FUN = function(x){
  strsplit(x, split = '_')[[1]][3]
})

VARIANTdf$pos <- factor(VARIANTdf$pos, levels = 1:max(as.numeric(VARIANTdf$pos)))
VARIANTdf <- VARIANTdf[order(as.numeric(VARIANTdf$pos)),]

################################
# MUTATION FUNCTIONAL ANALYSIS #
################################

MUTATIONdf$pos <- NA
MUTATIONdf$ref <- NA
MUTATIONdf$alt <- NA

MUTATIONdf$genomicRegion <- "intergenic" # 3UTR | 5UTR | intergenic | il nome l'orf
MUTATIONdf$codonPosition <- NA # 1 | 2 | 3
MUTATIONdf$effect <- NA # S | NS
MUTATIONdf$AAvar <- NA

for(i in 1:nrow(MUTATIONdf)) {
  
  pos_ref_alt_i <- unlist(strsplit(MUTATIONdf$MuID[i],split = '_'))
  
  MUTATIONdf$pos[i] <- as.numeric(pos_ref_alt_i[1])
  MUTATIONdf$ref[i] <- pos_ref_alt_i[2] 
  MUTATIONdf$alt[i] <- pos_ref_alt_i[3]
  
}

MUTATIONdf$VF_median <- sapply(X = MUTATIONdf$MuID, 
                               FUN = function(x){
                                 m_VF <- var_freq_matrix_filt[,x]
                                 return(median(m_VF[m_VF > thrDETECTION & !is.na(m_VF)]))
                               }
)


for(i in 1:nrow(reference_Seq)) {
  mut_idx <- which(MUTATIONdf$pos >= reference_Seq$START[i] & MUTATIONdf$pos <= reference_Seq$END[i])
  
  if(reference_Seq$genomicRegion[i] == "5UTR") {
    
    MUTATIONdf$genomicRegion[mut_idx] <- "5UTR"
    
  } else if(reference_Seq$genomicRegion[i] == "3UTR") {
    
    MUTATIONdf$genomicRegion[mut_idx] <- "3UTR"
    
  } else {
    
    MUTATIONdf$genomicRegion[mut_idx] <- reference_Seq$ORFname[i]
    
    ref_seq_i <- strsplit(reference_Seq$sequence[i], '')[[1]]
    
    if((length(ref_seq_i) %% 3) != 0) {
      stop("Error! CDS is not a triplette sequence.")
    }
    
    mod_seq_i_effect <- ref_seq_i
    
    
    mod_seq_i_kaks_SigClst <- list()
    
    for(sc in unique(statistics$SignatureCluster)) {
      mod_seq_i_kaks_SigClst[[sc]] <- ref_seq_i
    }
    
    ref_proteine_i <- seqinr::translate(seq = ref_seq_i, frame = 0)
    
    if((length(ref_proteine_i) * 3) != length(ref_seq_i)) {
      stop("Error! Translated proteine is truncated.")
    }
    
    cat(paste0('\n', reference_Seq$ORFname[i], '\n'))
    pb <- txtProgressBar(min = 0, max = length(mut_idx), style = 3)
    j = 0
    for(m_idx_i in mut_idx){
      
      j = j + 1
      setTxtProgressBar(pb, j)
      
      # handling frameshift in position 13468 (ref. NCBI)
      if(reference_Seq$ORFname[i] == "ORF1ab" & MUTATIONdf$pos[m_idx_i] > 13468) {
        pos_i <- MUTATIONdf$pos[m_idx_i] - reference_Seq$START[i] + 2
      } else {
        pos_i <- MUTATIONdf$pos[m_idx_i] - reference_Seq$START[i] + 1
      }
      
      if(MUTATIONdf$ref[m_idx_i]!=ref_seq_i[pos_i]) {
        stop("Error! reference allele and reference sequence don't match")
      }
      
      if(MUTATIONdf$pos[m_idx_i] == 13468) {
        MUTATIONdf$codonPosition[m_idx_i] <- NA
        
        mod_seq_i_effect[c(pos_i, pos_i +1)] <- MUTATIONdf$alt[m_idx_i]
        
      } else {
        MUTATIONdf$codonPosition[m_idx_i] <- ifelse(test = (pos_i %% 3) == 0, yes = 3, no = (pos_i %% 3))
        
        mod_seq_i_effect[pos_i] <- MUTATIONdf$alt[m_idx_i]
      }
      
      
      mod_proteine_i <- seqinr::translate(seq = mod_seq_i_effect, frame = 0)
      
      IDXsameAA <- mod_proteine_i == ref_proteine_i
      
      posAAmut <- NA
      if(sum(!IDXsameAA) == 0) {
        MUTATIONdf$effect[m_idx_i] <- "S"
        
      } else {
        
        MUTATIONdf$effect[m_idx_i] <- "NS"
        
        posAAmut <- which(IDXsameAA == FALSE)
        
        MUTATIONdf$AAvar[m_idx_i] <- paste0(reference_Seq$ORFname[i], ".", posAAmut,ref_proteine_i[posAAmut],'>',mod_proteine_i[posAAmut] )
        
      }
      
      mod_seq_i_effect <- ref_seq_i
      
      # kaks
      
      nucl_alt <- MUTATIONdf$alt[m_idx_i]
      mut_type <- MUTATIONdf$type[m_idx_i]
      
      idxMutDup <- which(MUTATIONdf$pos == MUTATIONdf$pos[m_idx_i])
      
      if(length(idxMutDup) > 1) {
        
        # Position with >2 alleles -> compute frequency
        VFsum <-colSums(var_freq_matrix_filt[, MUTATIONdf$MuID[idxMutDup]], na.rm = T)
        
        idVFmax <- which(VFsum == max(VFsum, na.rm = T))[1]
        
        nucl_alt <- MUTATIONdf$alt[idxMutDup[idVFmax]]
        mut_type <- MUTATIONdf$type[idxMutDup[idVFmax]]
        
      }
      
      for(sc in unique(VARIANTdf$SignatureCluster[VARIANTdf$MuID == MUTATIONdf$MuID[m_idx_i]])) {
        
        if(MUTATIONdf$pos[m_idx_i] == 13468) {  
          
          mod_seq_i_kaks_SigClst[[sc]][c(pos_i, pos_i+1)] <- nucl_alt
          
        } else {
          mod_seq_i_kaks_SigClst[[sc]][pos_i] <- nucl_alt
        }
      }
      
    }
  
  }
}

mutation_multAllele <- MUTATIONdf[MUTATIONdf$pos %in% as.numeric(names(which(table(MUTATIONdf$pos) > 1))),]

MUTATIONdf$effect[is.na(MUTATIONdf$effect)] <- "NC"

MUTATIONdf$effect <- factor(MUTATIONdf$effect, levels = c('S', 'NS', 'NC'))

MUTATIONdf$Variation <- paste0(MUTATIONdf$ref, ">", MUTATIONdf$alt)

VARIANTdf$effect <- MUTATIONdf$effect[match(VARIANTdf$MuID, MUTATIONdf$MuID)]
VARIANTdf$effect <- factor(as.character(VARIANTdf$effect), levels = c("S", "NS", "NC"))

VARIANTdf <- merge(x = VARIANTdf, y = statistics[,c('Run', 'SignatureCluster', 'nMinor', 'VERSO', 'Week')], by.x = 'Sample', by.y = 'Run')

# CODING - NON CODING

MUTATIONdf$Coding <- ifelse(MUTATIONdf$genomicRegion %in% c("5UTR", "3UTR", "intergenic"), yes = "NC", no = "C")
tbl_variation <- as.data.frame(table(MUTATIONdf$Variation, MUTATIONdf$type, MUTATIONdf$Coding))
colnames(tbl_variation) <- c("variation", "type", "coding", "count")

# SYNONIM - NON SYNONIM

freq_nucleotide <- table(REFANC)

#######################
# VARIATION LANDSCAPE #
#######################

VARIANTdf <- merge(x = VARIANTdf, y = MUTATIONdf[, c('MuID', 'Variation', 'TrueTrans')], by = 'MuID')

nPatient <- as.data.frame(table(statistics$nMinor))

VARIANTdf <- merge(x = VARIANTdf, y = nPatient, by.x = 'nMinor', by.y = 'Var1')
colnames(VARIANTdf)[ncol(VARIANTdf)] <- 'nPatientWnMinor'

VARIANTdf$Equivalence <- NA
VARIANTdf$Equivalence[VARIANTdf$Variation == 'G>T' | VARIANTdf$Variation == 'C>A'] <- 'C>A | G>T'
VARIANTdf$Equivalence[VARIANTdf$Variation == 'C>G' | VARIANTdf$Variation == 'G>C'] <- 'C>G | G>C'
VARIANTdf$Equivalence[VARIANTdf$Variation == 'C>T' | VARIANTdf$Variation == 'G>A'] <- 'C>T | G>A'
VARIANTdf$Equivalence[VARIANTdf$Variation == 'T>A' | VARIANTdf$Variation == 'A>T'] <- 'T>A | A>T'
VARIANTdf$Equivalence[VARIANTdf$Variation == 'T>C' | VARIANTdf$Variation == 'A>G'] <- 'T>C | A>G'
VARIANTdf$Equivalence[VARIANTdf$Variation == 'T>G' | VARIANTdf$Variation == 'A>C'] <- 'T>G | A>C'

nAlwaysMinor <- as.data.frame(table(VARIANTdf$Sample, VARIANTdf$type))
nAlwaysMinor <- nAlwaysMinor[nAlwaysMinor$Var2 == "always minor",c(1,3)]
colnames(nAlwaysMinor) <- c("Sample","nALWmin")

VARIANTdf <- merge(x = VARIANTdf, y = nAlwaysMinor, by = "Sample")

# backup workspace relevant objects for generating images from processed data
save(thrDETECTION, 
     thrCLONAL, 
     thrTRANS, 
     reference_Seq, 
     REFANC,
     freq_nucleotide, 
     var_freq_matrix_filt, 
     processed_variants_na_filt, 
     MUTATIONdf, 
     VARIANTdf,
     TRANSITIONdf,
     statistics,
     SignatureColor,
     file = "data/WorkSpace4plots.RData")

###########
# FIGURES #
###########

# cleaning
rm(list = ls())

#LOAD LIBRERIES
library(ggplot2)
library(gridExtra)

load("data/WorkSpace4plots.RData")
#################################
# FIGURE 1 - AGGREGATE ANALYSIS #
#################################

nObs_table <- table(statistics$nClonal, statistics$nMinor)
nObs <- c()
for(i in 1:nrow(statistics)){
  nObs <- c(nObs, nObs_table[as.numeric(rownames(nObs_table)) == statistics$nClonal[i],
                             as.numeric(colnames(nObs_table)) == statistics$nMinor[i]])
}

# Scatter clonal vs minor colorate per week
pdf(file = "plots/Fig1_A.pdf", width = 9, height = 6)
hist_top <- ggplot(statistics, aes(nClonal)) + 
  scale_x_continuous(breaks = seq(0,max(statistics$nClonal),by = 2)) +
  geom_histogram(binwidth = 0.5) + 
  xlab('') +
  ylab('count') +
  theme_bw() + 
  theme(text = element_text(size = 15))    
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank())

scatter <- ggplot(statistics, aes(x=nClonal, y=nMinor)) +
  geom_point(aes(size = nObs)) +
  xlab('# clonal SNVs') + 
  ylab('# minor SNVs') +
  scale_x_continuous(breaks = seq(0,max(statistics$nClonal),by = 2)) +
  scale_y_continuous(breaks = seq(0,max(statistics$nMinor),by = 10)) +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 20))

hist_right <- ggplot(statistics, aes(y = nMinor)) +
  geom_histogram(binwidth = 0.75)+
  scale_y_continuous(breaks = seq(0,max(statistics$nMinor),by = 10)) +
  ylab('') +
  xlab('count') +
  theme_bw() + 
  theme(text = element_text(size = 12), axis.text.x = element_text(angle = -45, vjust = -0.5)) 
grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 2), heights=c(1, 4))
dev.off()

stat_box_data <- function(y, upper_limit = max(statistics$nClonal) * 1.15) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('n =', length(y), '\n')
    )
  )
}

pdf(file = "plots/Fig1_B.pdf", width = 7, height = 6)
ggplot(statistics[!is.na(statistics$CollectionDate),], aes(x=as.factor(WeekMerged), y=nClonal)) + 
  geom_boxplot() + 
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.5
  ) + 
  xlab("Week") +
  ylab('# clonal SNVs') +
  theme_bw() + 
  theme(text = element_text(size = 20))
dev.off()

stat_box_data <- function(y, upper_limit = max(statistics$nMinor) * 1.15) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('n =', length(y), '\n')
    )
  )
}

pdf(file = "plots/Fig1_C.pdf", width = 7, height = 6)
ggplot(statistics[!is.na(statistics$CollectionDate),], aes(x=as.factor(WeekMerged), y=nMinor)) + 
  geom_boxplot() + 
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.5
  ) +
  xlab("Week") +
  ylab('# minor SNVs') +
  theme_bw() + 
  theme(text = element_text(size = 20))
dev.off()

stat_box_data <- function(y, upper_limit = max(VARIANTdf$VF, na.rm = T) * 1.1) {
  return(
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('m =', length(y), '\n')
    )
  )
}
pdf(file = "plots/Fig1_G.pdf", width = 6, height = 5)
ggplot(data = VARIANTdf, aes(x = type, fill = type, y = VF)) +
  geom_violin(scale = "width", draw_quantiles = 0.5)+
  stat_summary(
    fun.data = stat_box_data,
    geom = "text",
    hjust = 0.5,
    vjust = 0.5
  )+ 
  scale_y_continuous(breaks = seq(0.0,1,0.2), limits = c(0,1.05)) + 
  geom_hline(yintercept = thrTRANS, linetype = "dashed") +
  xlab("") +
  ylab('Variant Frequency') +
  theme_bw() + 
  scale_fill_manual(values = c("#CA644F", "#EBABA7", "#D3B3E3"))+ 
  theme(text = element_text(size = 20), legend.position="none")
dev.off()


library(dplyr)
library(ggforce) # for 'geom_arc_bar'

# INCLUDE POLYMORPHIC SITES

TableCountMut <- data.frame(hit = factor(c('Non mutated', 'Single hit', 'Multiple hits'), 
                                         levels = c('Non mutated', 'Single hit', 'Multiple hits')))
TableCountMut$count <- 0

TableCountMut$count[TableCountMut$hit == 'Single hit'] <- sum(!(duplicated(MUTATIONdf$pos) | duplicated(MUTATIONdf$pos, fromLast = T)))
TableCountMut$count[TableCountMut$hit == 'Multiple hits'] <- length(unique(MUTATIONdf$pos)) - TableCountMut$count[TableCountMut$hit == 'Single hit']
TableCountMut$count[TableCountMut$hit == 'Non mutated'] <- max(reference_Seq$END) - length(unique(MUTATIONdf$pos))

pdf(file = "plots/Fig1_D.pdf", width = 3, height = 6)
ggplot(data = TableCountMut) + 
  geom_bar(mapping = aes(x = 'Genome Position', fill = hit, weight = count), width = 0.33, position = 'fill') +
  facet_wrap(facets = ~ '') +
  ylab('ratio') + 
  scale_fill_manual(values = c("#aaa9c2", "#85dac8", "#c0a5a3"))+
  theme_bw() + 
  theme(text = element_text(size = 20), legend.position = 'none')
dev.off()

svg("plots/Fig1_E.svg", width = 3, height = 6)
ggplot(data = MUTATIONdf, mapping = aes(x = '', fill =type)) + 
  geom_bar(position = "fill", width = 0.33) +
  scale_fill_manual(values = c("#CA644F", "#EBABA7", "#D3B3E3"))+
  ylab('ratio') + 
  theme_bw() + 
  theme(text = element_text(size = 20), legend.position = 'none')
dev.off()

svg("plots/Fig1_F.svg", width = 6, height = 6)
ggplot(data = MUTATIONdf, mapping = aes(x = type, fill = effect)) + 
  geom_bar(position = "fill", width = 0.5)  +
  scale_fill_manual(values = c("#B1D8A0", "#8EC3ED", "#DBCD9D"))+
  ylab('ratio') + 
  theme_bw() + 
  theme(text = element_text(size = 20))
dev.off()


nMinorLM <- lm(nMinor~MedianCoverage, data = statistics)
pdf(file = "plots/FigS1_A.pdf", width = 6, height = 6)
ggplot(data = statistics, mapping = aes(x = MedianCoverage, y = nMinor)) + 
  geom_point() + 
  annotate("text", x=Inf, y = Inf, label = paste0("R2 = ", round(summary(nMinorLM)$adj.r.squared, digits = 4)), vjust=2, hjust=1.5) + 

  xlab("Median Coverage") +
  ylab('# minor SNVs') +
  theme_bw() + 
  theme(text = element_text(size = 20))
dev.off()


nMinorLM <- lm(nMinor~TotalCoverage, data = statistics)
pdf(file = "plots/FigS1_B.pdf", width = 6, height = 6)
ggplot(data = statistics, mapping = aes(x = TotalCoverage, y = nMinor)) + 
  geom_point() + 
  annotate("text", x=Inf, y = Inf, label = paste0("R2 = ", round(summary(nMinorLM)$adj.r.squared, digits = 4)), vjust=2, hjust=1.5) + 
  xlab("Total coverage") +
  ylab('# minor SNVs') +
  theme_bw() + 
  theme(text = element_text(size = 20))
dev.off()

###################
# FIGURE S2 - FUNCTIONAL ANALYSIS
##################
############
# boxlot - VF - position
######
pos_orf <- c()
orf_name <- c()

for(i in 1:nrow(reference_Seq)){
  start_orf <- min(as.numeric(VARIANTdf$pos[as.numeric(VARIANTdf$pos) >= reference_Seq$START[i]]))
  middle_orf <- min(as.numeric(VARIANTdf$pos[as.numeric(VARIANTdf$pos) >= ((reference_Seq$END[i] - reference_Seq$START[i]) / 2) + reference_Seq$START[i]]))
  end_orf <- max(as.numeric(VARIANTdf$pos[as.numeric(VARIANTdf$pos) <= reference_Seq$END[i]]))
  pos_orf <- c(pos_orf, middle_orf)
  orf_name <- c(orf_name, c(reference_Seq$ORFname[i]))
}

pos_union <- seq(0,max(reference_Seq$END),by = 2000)

pos_labels <- c()
for(i in 1:length(pos_union)){
  if(pos_union[i] %in% pos_orf){
    pos_labels <- c(pos_labels, orf_name[pos_orf == pos_union[i]])
  } else {
    pos_union[i] <- min(as.numeric(VARIANTdf$pos[as.numeric(VARIANTdf$pos) >= pos_union[i]]))
    
    pos_labels <- c(pos_labels, as.character(pos_union[i]))
  }
}

svg(filename = "plots/FigS2.svg", width = 19, height = 18)
ggplot(data = VARIANTdf, aes(x = pos, y = VF, color=effect)) + 
  geom_boxplot(outlier.size = 0.1, outlier.shape = 4) + 
  stat_summary(fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }, geom="crossbar", 
               size = 0.6, fatten=0, width=0.7, aes(group=1), color = "#00000050") +
  geom_hline(yintercept = 0.2, linetype = 'dashed') +
  facet_grid(type ~ .) + 
  scale_x_discrete("Genome position", breaks = pos_union, labels=pos_labels) + 
  ylab('Variant frequency') +
  scale_color_manual(values = c("#B1D8A0", "#8EC3ED", "#DBCD9D"))+
  theme_bw() + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


###################
# FIGURE 2 
##################


tbl_variation <- as.data.frame(table(MUTATIONdf$Variation, MUTATIONdf$type, MUTATIONdf$effect))
colnames(tbl_variation) <- c("variation", "type", "effect", "count")

tbl_variation$nucl_freq <- sapply(as.character(tbl_variation$variation), FUN = function(x){
  freq_nucleotide[strsplit(x, split = '>')[[1]][1]]
})
tbl_variation$type <- factor(as.character(tbl_variation$type), levels = c("always clonal", "mixed", "always minor"))
tbl_variation$effect <- factor(as.character(tbl_variation$effect), levels = c("S", "NS", "NC"))

pdf(file = "plots/Fig2_B.pdf", width = 7, height = 7)
ggplot(tbl_variation, aes(fill=type, y=count/nucl_freq, x=variation)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values =  c("#CA644F", "#EBABA7", "#D3B3E3"))+ 
  ylab('Frequency') + 
  xlab('Substitution')+
  theme_bw() + 
  theme(text = element_text(size = 20), legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

###################
# FIGURE 3 
##################

library(dplyr)
library(ggforce) # for 'geom_arc_bar'


TablePatientMut <- as.data.frame(table(statistics[,c('SignatureCluster')], useNA = 'ifany'), stringsAsFactors = F)
colnames(TablePatientMut) <- c("SignatureCluster", "count")
TablePatientMut$SignatureCluster[is.na(TablePatientMut$SignatureCluster)] <- '0 minors'
TablePatientMut$SignatureCluster <- factor(TablePatientMut$SignatureCluster, levels = c("SC#1","SC#2","SC#3","SC#4","0 minors"))
TablePatientMut <- TablePatientMut[TablePatientMut$count !=0,]

df <- TablePatientMut %>%
  mutate(end = 2 * pi * cumsum(count)/sum(count),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))

pdf("plots/Fig3_C.pdf", width = 10, height = 6)
ggplot(df) +
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                   start = start, end = end, fill = SignatureCluster)) +
  geom_text(aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = round((count/sum(count)*100), digits = 2),
                hjust = hjust, vjust = vjust, angle = 0), size = 4) +
  coord_fixed() +
  scale_x_continuous(limits = c(-1.5, 1.5),  # Adjust so labels are not cut off
                     name = "mixed", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1.5, 1.5),      # Adjust so labels are not cut off
                     name = "", breaks = NULL, labels = NULL) +
  scale_fill_manual(values = c("#ee746c","#7aaf2b", "#26b6bd","#c77aff", "grey90"))+
  theme_minimal() +
  theme(text = element_text(size = 15))
dev.off()


###################
# FIGURE 4 - SIGNATURE CLUSTER
##################

pdf(file = "plots/Fig4_A.pdf", width = 6, height = 6)
ggplot(data = statistics[!is.na(statistics$SignatureCluster),], aes(y = nClonal, x = SignatureCluster)) + 
  ylab('# clonal SNVs') + 
  xlab('Signature cluster') +
  geom_boxplot(aes(color = SignatureCluster)) +
  theme_bw() + 
  theme(text = element_text(size = 20), axis.text=element_text(size=13), legend.position = 'none')
dev.off()

pdf(file = "plots/Fig4_B.pdf", width = 6, height = 6)
ggplot(data = statistics[!is.na(statistics$SignatureCluster),], aes(y = nMinor, x = SignatureCluster)) + 
  ylab('# minor SNVs') + 
  xlab('Signature cluster') +
  scale_y_continuous(breaks = seq(0,max(statistics$nMinor), 10)) + 
  geom_boxplot(aes(color = SignatureCluster)) +
  theme_bw() +
  theme(text = element_text(size = 20), axis.text=element_text(size=13), legend.position = 'none')
dev.off()

stat_box_data <- function(y, upper_limit = max(VARIANTdf$VF, na.rm = T) * 1.1) {
  return(
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('m =', length(y), '\n')
    )
  )
}

pdf(file = "plots/Fig4_C.pdf", width = 6, height = 6)
ggplot(data = VARIANTdf[!is.na(VARIANTdf$SignatureCluster),], aes(y = VF, fill = SignatureCluster, x = SignatureCluster)) +
  geom_violin(scale = 'width', draw_quantiles = 0.5)+
  stat_summary(
    fun.data = stat_box_data,
    geom = "text",
    hjust = 0.5,
    vjust = 0.5
  ) + 
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1.05)) + 
  geom_hline(yintercept = 0.2, linetype = 'dashed') +
  xlab("") +
  ylab('Variant Frequency') +
  theme_bw() + 
  theme(text = element_text(size = 20),legend.position="none")
dev.off()



# Compute 

NtoNcount <- data.frame()

for(sc in unique(VARIANTdf$SignatureCluster[!is.na(VARIANTdf$SignatureCluster)])){
  for(orf in reference_Seq$ORFname[reference_Seq$genomicRegion == "CDS"]) {
    seq <- reference_Seq$sequence[reference_Seq$ORFname == orf]
    Atot <- stringi::stri_count(str = seq, regex = "A")
    Gtot <- stringi::stri_count(str = seq, regex = "G")
    Ctot <- stringi::stri_count(str = seq, regex = "C")
    Ttot <- stringi::stri_count(str = seq, regex = "T")
    
    posStart <- reference_Seq$START[reference_Seq$ORFname == orf]
    posEnd <- reference_Seq$END[reference_Seq$ORFname == orf]
    
    CtoT_i = length(unique(VARIANTdf[as.numeric(VARIANTdf$pos) >= posStart 
                                     & as.numeric(VARIANTdf$pos) <= posEnd 
                                     & VARIANTdf$ref == 'C' 
                                     & VARIANTdf$alt == 'T' 
                                     & VARIANTdf$SignatureCluster == sc,'MuID']))
    
    GtoA_i = length(unique(VARIANTdf[as.numeric(VARIANTdf$pos) >= posStart 
                                     & as.numeric(VARIANTdf$pos) <= posEnd 
                                     & VARIANTdf$ref == 'G' 
                                     & VARIANTdf$alt == 'A' 
                                     & VARIANTdf$SignatureCluster == sc,'MuID']))
    
    GtoT_i = length(unique(VARIANTdf[as.numeric(VARIANTdf$pos) >= posStart 
                                     & as.numeric(VARIANTdf$pos) <= posEnd 
                                     & VARIANTdf$ref == 'G' 
                                     & VARIANTdf$alt == 'T' 
                                     & VARIANTdf$SignatureCluster == sc,'MuID']))
    
    CtoA_i = length(unique(VARIANTdf[as.numeric(VARIANTdf$pos) >= posStart 
                                     & as.numeric(VARIANTdf$pos) <= posEnd 
                                     & VARIANTdf$ref == 'C' 
                                     & VARIANTdf$alt == 'A' 
                                     & VARIANTdf$SignatureCluster == sc,'MuID']))
    
    AtoG_i = length(unique(VARIANTdf[as.numeric(VARIANTdf$pos) >= posStart 
                                     & as.numeric(VARIANTdf$pos) <= posEnd 
                                     & VARIANTdf$ref == 'A' 
                                     & VARIANTdf$alt == 'G' 
                                     & VARIANTdf$SignatureCluster == sc,'MuID']))
    
    TtoC_i = length(unique(VARIANTdf[as.numeric(VARIANTdf$pos) >= posStart 
                                     & as.numeric(VARIANTdf$pos) <= posEnd 
                                     & VARIANTdf$ref == 'T' 
                                     & VARIANTdf$alt == 'C' 
                                     & VARIANTdf$SignatureCluster == sc,'MuID']))
    
    
    #Aggiunte
    CtoG_i = length(unique(VARIANTdf[as.numeric(VARIANTdf$pos) >= posStart 
                                     & as.numeric(VARIANTdf$pos) <= posEnd 
                                     & VARIANTdf$ref == 'C' 
                                     & VARIANTdf$alt == 'G' 
                                     & VARIANTdf$SignatureCluster == sc,'MuID']))
    
    GtoC_i = length(unique(VARIANTdf[as.numeric(VARIANTdf$pos) >= posStart 
                                     & as.numeric(VARIANTdf$pos) <= posEnd 
                                     & VARIANTdf$ref == 'G' 
                                     & VARIANTdf$alt == 'C' 
                                     & VARIANTdf$SignatureCluster == sc,'MuID']))
    
    AtoT_i = length(unique(VARIANTdf[as.numeric(VARIANTdf$pos) >= posStart 
                                     & as.numeric(VARIANTdf$pos) <= posEnd 
                                     & VARIANTdf$ref == 'A' 
                                     & VARIANTdf$alt == 'T' 
                                     & VARIANTdf$SignatureCluster == sc,'MuID']))
    
    AtoC_i = length(unique(VARIANTdf[as.numeric(VARIANTdf$pos) >= posStart 
                                     & as.numeric(VARIANTdf$pos) <= posEnd 
                                     & VARIANTdf$ref == 'A' 
                                     & VARIANTdf$alt == 'C' 
                                     & VARIANTdf$SignatureCluster == sc,'MuID']))
    
    TtoA_i = length(unique(VARIANTdf[as.numeric(VARIANTdf$pos) >= posStart 
                                     & as.numeric(VARIANTdf$pos) <= posEnd 
                                     & VARIANTdf$ref == 'T' 
                                     & VARIANTdf$alt == 'A' 
                                     & VARIANTdf$SignatureCluster == sc,'MuID']))
    
    TtoG_i = length(unique(VARIANTdf[as.numeric(VARIANTdf$pos) >= posStart 
                                     & as.numeric(VARIANTdf$pos) <= posEnd 
                                     & VARIANTdf$ref == 'T' 
                                     & VARIANTdf$alt == 'G' 
                                     & VARIANTdf$SignatureCluster == sc,'MuID']))
    
    
    
    
    NtoNcount <- rbind.data.frame(NtoNcount, data.frame(ORF = orf,
                                                        SigClst = sc,
                                                        CtoT = CtoT_i,
                                                        GtoA = GtoA_i,
                                                        GtoT = GtoT_i,
                                                        CtoA = CtoA_i,
                                                        AtoG = AtoG_i,
                                                        TtoC = TtoC_i,
                                                        
                                                        #Aggiunte
                                                        CtoG = CtoG_i,
                                                        GtoC = GtoC_i,
                                                        AtoT = AtoT_i,
                                                        AtoC = AtoC_i,
                                                        TtoA = TtoA_i,
                                                        TtoG = TtoG_i,
                                                        
                                                        A_tot = Atot,
                                                        T_tot = Ttot,
                                                        G_tot = Gtot,
                                                        C_tot = Ctot,
                                                        
                                                        CtoT_freq = CtoT_i / (Ctot + Gtot),
                                                        GtoA_freq = GtoA_i / (Ctot + Gtot),
                                                        GtoT_freq = GtoT_i / (Ctot + Gtot),
                                                        CtoA_freq = CtoA_i / (Ctot + Gtot),
                                                        AtoG_freq = AtoG_i / (Atot + Ttot),
                                                        TtoC_freq = TtoC_i / (Atot + Ttot),
                                                        
                                                        #Aggiunte
                                                        CtoG_freq = CtoG_i / (Ctot + Gtot),
                                                        GtoC_freq = GtoC_i / (Ctot + Gtot),
                                                        AtoT_freq = AtoT_i / (Atot + Ttot),
                                                        AtoC_freq = AtoC_i / (Atot + Ttot),
                                                        TtoA_freq = TtoA_i / (Atot + Ttot),
                                                        TtoG_freq = TtoG_i / (Atot + Ttot)
    )
    )
    
  }
}



NtoNcount_m <- reshape2::melt(NtoNcount[,c("ORF","SigClst","CtoT_freq","GtoA_freq","GtoT_freq","CtoA_freq","AtoG_freq","TtoC_freq",    
                                           "CtoG_freq","GtoC_freq","AtoT_freq","AtoC_freq","TtoA_freq","TtoG_freq")])

NtoNcount_m$ORF <- factor(NtoNcount_m$ORF, levels = c("ORF1ab","S","ORF3a","E","M","ORF6","ORF7a","ORF8","N","ORF10"))
colnames(NtoNcount_m) <- c("ORF","SigClst","Nucleotide_variation", "ratio")
NtoNcount_m$OmoG <- NA
NtoNcount_m$OmoG[NtoNcount_m$Nucleotide_variation== 'CtoT_freq' | NtoNcount_m$Nucleotide_variation== 'GtoA_freq'] <- 'Omo 1'
NtoNcount_m$OmoG[NtoNcount_m$Nucleotide_variation== 'GtoT_freq' | NtoNcount_m$Nucleotide_variation== 'CtoA_freq'] <- 'Omo 2'
NtoNcount_m$OmoG[NtoNcount_m$Nucleotide_variation== 'AtoG_freq' | NtoNcount_m$Nucleotide_variation== 'TtoC_freq'] <- 'Omo 3'

NtoNcount_m$OmoG[NtoNcount_m$Nucleotide_variation== 'CtoG_freq' | NtoNcount_m$Nucleotide_variation== 'GtoC_freq'] <- 'Omo 4'
NtoNcount_m$OmoG[NtoNcount_m$Nucleotide_variation== 'AtoT_freq' | NtoNcount_m$Nucleotide_variation== 'TtoA_freq'] <- 'Omo 5'
NtoNcount_m$OmoG[NtoNcount_m$Nucleotide_variation== 'TtoG_freq' | NtoNcount_m$Nucleotide_variation== 'AtoC_freq'] <- 'Omo 6'

NtoNcount_m$ORF_OMOG <- factor(paste0(NtoNcount_m$ORF,'_',NtoNcount_m$OmoG), 
                               levels = paste0(rep(levels(NtoNcount_m$ORF), each = 6), '_', unique(NtoNcount_m$OmoG)))

breaks <- as.character(NtoNcount_m$ORF_OMOG)
labels <- sapply(breaks, FUN = function(x){strsplit(x,'_')[[1]][1]}, USE.NAMES = F)

NtoNcount_m$Nucleotide_variation <- factor(NtoNcount_m$Nucleotide_variation, levels = c("CtoT_freq","GtoA_freq",
                                                                                        "GtoT_freq","CtoA_freq",
                                                                                        "AtoG_freq","TtoC_freq",
                                                                                        "CtoG_freq","GtoC_freq",
                                                                                        "AtoT_freq","TtoA_freq",
                                                                                        "TtoG_freq","AtoC_freq"))
pdf(file = "plots/FigS3.pdf", width = 19, height = 5)
ggplot(data = NtoNcount_m, mapping = aes(x = ORF_OMOG, y = ratio, fill = Nucleotide_variation))+
  geom_bar(width = 0.8, position = "stack", stat = "identity") +
  facet_wrap(. ~ SigClst, ncol = 4) + 
  scale_x_discrete("ORF", breaks = breaks[endsWith(breaks, suffix = 'Omo 1')], labels = labels[endsWith(breaks, suffix = 'Omo 1')] ) +
  scale_fill_manual(values = c("#e5736b","#e4aaa6","#7cac83","#3cac4e","#6387b8","#b3c1d4", "#a79531","#e8dda3", "#33b5b1","#9de4e1","#ae6a9e","#e0abd4"))+
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.x = element_text(size = 15, angle = 45, hjust = 1)) 
dev.off()


# ######################
# # FIGURE 5 - HEATMAP #
# ######################
 


#####################
# HEATMAP X LINEAGE - Versione 2 #
# X = nLineage
# Y = nPatient
# color = Count
#####################

library('pheatmap')
library('ggplot2')
library('RColorBrewer')
library('dplyr')
library('ggforce')
library('viridis')

inference_step2 <- readRDS("data/inference_step2.rds")

genotypes <- unique(inference_step2$corrected_genotypes)

genotypes <- cbind(genotypes, genotypes[,'28882_G_A|28883_G_C'])
colnames(genotypes)[which(colnames(genotypes) == '28882_G_A|28883_G_C')] <- '28882_G_A'
colnames(genotypes)[ncol(genotypes)] <- '28883_G_C'

nMUTpat_list <- list()

VARIANTdffilt <- VARIANTdf[(!(VARIANTdf$MuID %in% colnames(genotypes)) & 
                             VARIANTdf$VF <= thrCLONAL),]


for(SC in unique(VARIANTdffilt$SignatureCluster[!is.na(VARIANTdffilt$SignatureCluster)])) {
  
  mtx_pat_lineage <- matrix(data = 0, ncol = length(unique(VARIANTdffilt$VERSO)), nrow = length(unique(VARIANTdffilt$Sample)))
  
  nMUTpat_list[[SC]] <- data.frame(MuID = unique(VARIANTdffilt$MuID[VARIANTdffilt$SignatureCluster == SC & 
                                                                  !is.na(VARIANTdffilt$SignatureCluster)]))
  
  nMUTpat_list[[SC]]$nPat <- 0
  
  colnames(mtx_pat_lineage) <- 1:ncol(mtx_pat_lineage)
  rownames(mtx_pat_lineage) <- 1:nrow(mtx_pat_lineage)
  
  for(m in nMUTpat_list[[SC]]$MuID) {
    
    nP <- length(unique(VARIANTdffilt$Sample[VARIANTdffilt$MuID == m & 
                                               VARIANTdffilt$SignatureCluster == SC & 
                                           !is.na(VARIANTdffilt$SignatureCluster)]))
    
    nL <- length(unique(VARIANTdffilt$VERSO[VARIANTdffilt$MuID == m & 
                                                      VARIANTdffilt$SignatureCluster == SC & 
                                                  !is.na(VARIANTdffilt$SignatureCluster)]))
    
    nMUTpat_list[[SC]]$nPat[nMUTpat_list[[SC]]$MuID == m] <- nP
    
    mtx_pat_lineage[nP,nL] <- mtx_pat_lineage[nP,nL] + 1
    
  }

  
  rs <- rowSums(mtx_pat_lineage)
  if(sum(rs>0) > 1) {
    mtx_pat_lineage <- mtx_pat_lineage[rs!=0,, drop = FALSE]
    
    mtx_pat_lineage <- mtx_pat_lineage[nrow(mtx_pat_lineage):1,,drop = FALSE]
    
    mtx_pat_lineage[which(mtx_pat_lineage == 0)] <- NA
    
    
    pdf(paste0("plots/Fig5_C_heatmap_nPatXnLin_", SC, ".pdf"), width = 6, height = 16)
    print(
      pheatmap(log(mtx_pat_lineage, base = 10), 
               color = viridis::plasma(100, direction = 1),
               border_color = "grey60",
               cellwidth = 10,
               cellheight = 10,
               na_col = "#FFFFFF",
               
               scale = "none",
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               
               legend = TRUE,
               
               
               legend_breaks = seq(0,max(log10(mtx_pat_lineage), na.rm = T),1),
               legend_labels = round(10^seq(0,max(log10(mtx_pat_lineage), na.rm = T),1)),
               
               drop_levels = TRUE,
               show_rownames = T,
               show_colnames = T, 
               angle_col = 45,
               
               fontsize = 10,
               
               silent = FALSE,
      )
    )
    dev.off()
  }
  
}


######
# PIECHART private - mono lineage - other
######
library(dplyr)
library(ggforce)

Mut_nP_nL <- data.frame(MuID = setdiff(MUTATIONdf$MuID, colnames(genotypes)))
Mut_nP_nL$nP <- 0
Mut_nP_nL$nL <- 0

for(m in Mut_nP_nL$MuID){
  Mut_nP_nL$nP[Mut_nP_nL$MuID == m] <- length(unique(VARIANTdf$Sample[
    VARIANTdf$VF <= thrCLONAL &
    VARIANTdf$MuID == m
  ]))
  
  Mut_nP_nL$nL[Mut_nP_nL$MuID == m]  <- length(unique(VARIANTdf$VERSO[
    VARIANTdf$VF <= thrCLONAL &
    VARIANTdf$MuID == m
  ]))
  
}
Mut_nP_nL <- Mut_nP_nL[Mut_nP_nL$nP > 0 & Mut_nP_nL$nL > 0,]
Mut_nP_nL$class <- NA
Mut_nP_nL$class[Mut_nP_nL$nL == 1] <- '1 clade'
Mut_nP_nL$class[Mut_nP_nL$nP == 1] <- 'private'
Mut_nP_nL$class[Mut_nP_nL$nL >= 2] <- '>=2 clades'

Mut_nP_nL$class <- factor(Mut_nP_nL$class, levels = c('private', '1 clade', '>=2 clades'))

df <- as.data.frame(table(Mut_nP_nL$class)) %>%
  mutate(end = 2 * pi * cumsum(Freq)/sum(Freq),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1)
  )

pdf("plots/Fig5_E.pdf", width = 10, height = 6)
ggplot(df) +
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                   start = start, end = end, fill = Var1)) +
  geom_text(aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = round((Freq/sum(Freq)*100), digits = 2),
                hjust = hjust, vjust = vjust, angle = 0), size = 4) +
  coord_fixed() +
  scale_x_continuous(limits = c(-1.5, 1.5),  # Adjust so labels are not cut off
                     name = "always minor and mixed", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1.5, 1.5),      # Adjust so labels are not cut off
                     name = "", breaks = NULL, labels = NULL) +
  scale_fill_manual(values = c("#ffcd30","#C23B22","#238B45"))+
  theme_minimal() +
  theme(text = element_text(size = 15))
dev.off()

################################
# VARIATION LANDSCAPE - FIGURE #
################################

pdf(file = 'plots/Fig4_D.pdf', width = 10, height = 3)
ggplot(data = VARIANTdf[VARIANTdf$type == "always minor",],
       mapping = aes(x = nALWmin, fill = Equivalence)) +
  geom_histogram(binwidth = 1, position = 'fill') + 
  #coord_cartesian(ylim = c(0,100)) + 
  #scale_y_continuous(breaks = seq(0,12,by = 2))+
  facet_grid(. ~ SignatureCluster)+
  scale_fill_manual(values = c("#3cac4e","#a79531","#e5736b", "#33b5b1", "#6387b8", "#ae6a9e"))+
  theme_bw()
dev.off()

###########################
# SLIDING WINDOWS - KA/KS #
###########################

SignaturePrevalence <- readRDS('data/SignaturePrevalance.rds')
rownames(SignaturePrevalence) <- c("R","SC#1","SC#2","SC#3")

reference_Seq <- reference_Seq[order(reference_Seq$START),]

CDSseq_ref <- NULL
GENpos_ref <- NULL
for(i in 1:nrow(reference_Seq)) {
  
  if(reference_Seq$genomicRegion[i] == "CDS") {
    
    CDSseq_ref <- c(CDSseq_ref, unlist(strsplit(reference_Seq$sequence[i], '')))
    
    if(reference_Seq$ORFname[i] == 'ORF1ab') {
      GENpos_ref <- c(GENpos_ref, reference_Seq$START[i]:13468, 13468:reference_Seq$END[i])
    } else {
      GENpos_ref <- c(GENpos_ref, reference_Seq$START[i]:reference_Seq$END[i])
    }
  }
}

### Parameters ###
WinSize <- 3*100

SlideSize <- 3*1
#################

dNdS_sign <- list()

for(SC in c('SC#1','SC#2','SC#3')) {
  
  dNdS_sign[[SC]] <- list()
  
  cat('\nCluster: ', SC, '\n')
  
  SignPrevSC <- c(SignaturePrevalence[SC,], SignaturePrevalence[SC,])
  names(SignPrevSC) <- c("C>A","C>G","C>T","T>A","T>C","T>G","G>T","G>C","G>A","A>T","A>G","A>C")
  

  MUTATIONsc <- MUTATIONdf[MUTATIONdf$MuID %in% VARIANTdf$MuID[VARIANTdf$SignatureCluster == SC & 
                                                                 !is.na(VARIANTdf$SignatureCluster)],]
  
  Nstart = 1
  Nend = Nstart + WinSize - 1
  
  pos <- c()
  
  pb <- txtProgressBar(min = 0, max = length(CDSseq_ref), style = 3)
  while(Nend <= length(CDSseq_ref)) {
    
    NNseq_ref <- CDSseq_ref[Nstart:Nend]
    POSref <- GENpos_ref[Nstart:Nend]
    AAref <- seqinr::translate(seq = NNseq_ref, frame = 0)
    
    S_NS_obs <- table(MUTATIONsc$effect[MUTATIONsc$pos %in% POSref])
    
    type_obs <- table(MUTATIONsc$type[MUTATIONsc$pos %in% POSref])
    
    dNdSobs <- as.numeric(S_NS_obs['NS'] / S_NS_obs['S'])
    
    dNdS_sign[[SC]]$nNS <-  c(dNdS_sign[[SC]]$nNS, as.numeric(S_NS_obs['NS']))
    dNdS_sign[[SC]]$nS <-  c(dNdS_sign[[SC]]$nS, as.numeric(S_NS_obs['S']))
    
    dNdS_sign[[SC]]$alw_clonal <- c(dNdS_sign[[SC]]$alw_clonal, as.numeric(type_obs['always clonal']))
    dNdS_sign[[SC]]$alw_minor <- c(dNdS_sign[[SC]]$alw_minor, as.numeric(type_obs['always minor']))
    dNdS_sign[[SC]]$mixed <- c(dNdS_sign[[SC]]$mixed, as.numeric(type_obs['mixed']))
    
    S_sig <- 0

    NS_sig <- 0

    for(i in 1:length(NNseq_ref)){
      
      NNseq_alt <- NNseq_ref 
      for(Nalt in c("A", "T", "C", "G")) {
        
        if(NNseq_ref[i] != Nalt) {
          
          NNseq_alt[i] <- Nalt
          
          AApos <- ceiling(i/3)
          
          CODalt <- NNseq_alt[((AApos*3) -2) : (AApos*3)]
          
          if(AAref[AApos] ==  seqinr::translate(seq = CODalt, frame = 0)) {
            #Synonime Mutation
            
            S_sig <- S_sig + as.numeric(SignPrevSC[paste0(NNseq_ref[i],">",Nalt)])

          } else {
            #Non Mutation
            NS_sig <- NS_sig + as.numeric(SignPrevSC[paste0(NNseq_ref[i],">",Nalt)])
          }
          
        }
        
      }
      
    }
    
    dNdSexp_sig <- NS_sig/S_sig
    

    dNdS_sign[[SC]]$res <- c(dNdS_sign[[SC]]$res, (dNdSobs / dNdSexp_sig))
    pos <- c(pos, POSref[1])

    Nstart <- Nstart + SlideSize
    Nend <- (Nstart + WinSize - 1)
    setTxtProgressBar(pb, Nend)
  }
  dNdS_sign[[SC]]$res[which(dNdS_sign[[SC]]$res == Inf)] <- NA
  
  type_sc <- table(MUTATIONsc$type)
  
  dNdS_sign[[SC]]$alw_clonal_norm <- dNdS_sign[[SC]]$alw_clonal / as.numeric(type_sc['always clonal'])
  dNdS_sign[[SC]]$alw_minor_norm <- dNdS_sign[[SC]]$alw_minor / as.numeric(type_sc['always minor'])
  dNdS_sign[[SC]]$mixed_norm <- dNdS_sign[[SC]]$mixed / as.numeric(type_sc['mixed'])
  
  dNdS_sign[['POS']] <- pos
}

reference_Seq$ORFname <- factor(reference_Seq$ORFname, levels = reference_Seq$ORFname)

dNdSdf <- as.data.frame(dNdS_sign)
dNdSdf <- reshape2::melt(dNdSdf, id = 'POS')

pdf(file = "plots/Fig4_E_dnds_nonSmoot.pdf", width = 18, height = 9)
ggplot(dNdSdf[endsWith(as.character(dNdSdf$variable), "res"),], mapping = aes(y = value, x = POS, color = variable)) + 
  geom_line() +
  facet_grid(variable ~ ., scales = "free_y") +
  scale_color_manual(values = as.character(SignatureColor[1:3])) + 
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_rect(data=reference_Seq[reference_Seq$genomicRegion=="CDS",], 
            mapping=aes(xmin=START, xmax=END, ymin=-0.2, ymax=0, fill=ORFname), 
            inherit.aes = FALSE, color="black", alpha=0.75) +
  theme_bw()
dev.off()

########### smoothing start ##########
dNdS_res_mean <- list()
dNdS_res_mean[['POS']] <- GENpos_ref

### smoothing with moving average
for(SC in c('SC#1','SC#2','SC#3')) {
  dNdS_res_mean[[SC]] <- list()
  for(p in dNdS_res_mean[['POS']]) {
    dNdS_res_mean[[SC]]$res <- c(dNdS_res_mean[[SC]]$res, mean(dNdS_sign[[SC]]$res[dNdS_sign$POS <= p & (dNdS_sign$POS + WinSize) > p], na.rm = T))
    dNdS_res_mean[[SC]]$nNS <- c(dNdS_res_mean[[SC]]$nNS, mean(dNdS_sign[[SC]]$nNS[dNdS_sign$POS <= p & (dNdS_sign$POS + WinSize) > p], na.rm = T))
    dNdS_res_mean[[SC]]$nS <- c(dNdS_res_mean[[SC]]$nS, mean(dNdS_sign[[SC]]$nS[dNdS_sign$POS <= p & (dNdS_sign$POS + WinSize) > p], na.rm = T))
    
    dNdS_res_mean[[SC]]$alw_clonal <- c(dNdS_res_mean[[SC]]$alw_clonal, mean(dNdS_sign[[SC]]$alw_clonal[dNdS_sign$POS <= p & (dNdS_sign$POS + WinSize) > p], na.rm = T))
    dNdS_res_mean[[SC]]$alw_minor <- c(dNdS_res_mean[[SC]]$alw_minor, mean(dNdS_sign[[SC]]$alw_minor[dNdS_sign$POS <= p & (dNdS_sign$POS + WinSize) > p], na.rm = T))
    dNdS_res_mean[[SC]]$mixed <- c(dNdS_res_mean[[SC]]$mixed, mean(dNdS_sign[[SC]]$mixed[dNdS_sign$POS <= p & (dNdS_sign$POS + WinSize) > p], na.rm = T))
    
    dNdS_res_mean[[SC]]$alw_clonal_norm <- c(dNdS_res_mean[[SC]]$alw_clonal_norm, mean(dNdS_sign[[SC]]$alw_clonal_norm[dNdS_sign$POS <= p & (dNdS_sign$POS + WinSize) > p], na.rm = T))
    dNdS_res_mean[[SC]]$alw_minor_norm <- c(dNdS_res_mean[[SC]]$alw_minor_norm, mean(dNdS_sign[[SC]]$alw_minor_norm[dNdS_sign$POS <= p & (dNdS_sign$POS + WinSize) > p], na.rm = T))
    dNdS_res_mean[[SC]]$mixed_norm <- c(dNdS_res_mean[[SC]]$mixed_norm, mean(dNdS_sign[[SC]]$mixed_norm[dNdS_sign$POS <= p & (dNdS_sign$POS + WinSize) > p], na.rm = T))
  }
  dNdS_res_mean[[SC]]$totSNV <- dNdS_res_mean[[SC]]$nNS + dNdS_res_mean[[SC]]$nS 
  cat(paste0(SC, ' done\n'))
}



save(dNdS_sign,
     dNdS_res_mean,
     file = "data/dnds_analyses.RData")

dNdSdf_res_mean <- as.data.frame(dNdS_res_mean)
dNdSdf_res_mean <- reshape2::melt(dNdSdf_res_mean, id = c('POS'))

pdf(file = "plots/Fig4_E_dNdS_smoot.pdf", width = 18, height = 9)
ggplot(dNdSdf_res_mean[endsWith(as.character(dNdSdf_res_mean$variable), "res"), ], mapping = aes(y = value, x = POS, color = variable)) +
  geom_line() +
  facet_grid(variable ~ ., scales = "free_y") +
  scale_color_manual(values = as.character(SignatureColor[1:3])) + 
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_rect(data=reference_Seq, mapping=aes(xmin=START, xmax=END, ymin=-0.2, ymax=0, fill=ORFname), inherit.aes = FALSE, color="black", alpha=0.75) +
  theme_bw()
dev.off()


resDF <- dNdSdf_res_mean[endsWith(as.character(dNdSdf_res_mean$variable), "res"), ]
resDF$variable <- gsub('.res', '',resDF$variable)

totSNVdf <- dNdSdf_res_mean[endsWith(as.character(dNdSdf_res_mean$variable), "totSNV"), ]
totSNVdf$variable <- gsub('.totSNV', '',totSNVdf$variable)

pdf(file = "plots/Fig4_E_ratioS_NS.pdf", width = 18, height = 9)
ggplot() +
  geom_line(data = totSNVdf, mapping = aes(y = value, x = POS), color = "black", linetype = "dotted") + 
  scale_color_manual(values = as.character(SignatureColor[1:3])) + 
  facet_grid(variable ~ ., scales = "free_y") +
  geom_rect(data=reference_Seq, mapping=aes(xmin=START, xmax=END, ymin=-0.2, ymax=0, fill=ORFname), inherit.aes = FALSE, color="black", alpha=0.75) +
  theme_bw()
dev.off()

############## smoothin end ##########
 

#######################
# MANHATTAN-LIKE plot #
#######################

reference_Seq <- reference_Seq[order(reference_Seq$START),]
reference_Seq$ORFname <- factor(reference_Seq$ORFname, levels = reference_Seq$ORFname)
pdf(file = "plots/Fig2_A.pdf", width = 18, height = 5)
ggplot(data = VARIANTdf, mapping = aes(x = as.numeric(pos), y=VF, color = type)) + 
  geom_point(size = 0.5) + 
  geom_rect(data=reference_Seq, mapping=aes(xmin=START, xmax=END, ymin=-0.05, ymax=0, fill=ORFname), inherit.aes = FALSE, color="black", alpha=0.75) +
  ylab("Variant frequency") +
  xlab("Genome position") +
  scale_color_manual(values = c("#ca644f", "#ebaba7", "#d3b3e3"))+
  scale_x_continuous(breaks = seq(0,30000, by = 2000))+
  theme_bw() +
  theme(text = element_text(size = 20))
dev.off()

###############################
# dN/dS su TrueTrans (33 mut) #
###############################

MuID_TT <- MUTATIONdf$MuID[MUTATIONdf$TrueTrans == 'yes']

MuID_TT <- MUTATIONdf$MuID[MUTATIONdf$MuID %in% unique(VARIANTdf$MuID[VARIANTdf$SignatureCluster == 'SC#3' & !is.na(VARIANTdf$SignatureCluster)])]


reference_Seq <- reference_Seq[order(reference_Seq$START),]

CDSseq_ref <- NULL
GENpos_ref <- NULL
for(i in 1:nrow(reference_Seq)) {
  
  if(reference_Seq$genomicRegion[i] == "CDS") {
    
    CDSseq_ref <- c(CDSseq_ref, unlist(strsplit(reference_Seq$sequence[i], '')))
    
    if(reference_Seq$ORFname[i] == 'ORF1ab') {
      GENpos_ref <- c(GENpos_ref, reference_Seq$START[i]:13468, 13468:reference_Seq$END[i])
    } else {
      GENpos_ref <- c(GENpos_ref, reference_Seq$START[i]:reference_Seq$END[i])
    }
  }
}


AAref <- seqinr::translate(seq = CDSseq_ref, frame = 0)

S_NS_obs <- table(MUTATIONdf$effect[MUTATIONdf$MuID %in% MuID_TT])

dNdSobs <- as.numeric(S_NS_obs['NS'] / S_NS_obs['S'])

S_sig <- 0

NS_sig <- 0

pb <- txtProgressBar(min = 0, max = length(CDSseq_ref), style = 3)

for(i in 1:length(CDSseq_ref)){
  
  setTxtProgressBar(pb, i)
  
  CDSseq_alt <- CDSseq_ref 
  for(Nalt in c("A", "T", "C", "G")) {
    
    if(CDSseq_ref[i] != Nalt) {
      
      CDSseq_ref[i] <- Nalt
      
      AApos <- ceiling(i/3)
      
      CODalt <- CDSseq_alt[((AApos*3) -2) : (AApos*3)]
      
      if(AAref[AApos] ==  seqinr::translate(seq = CODalt, frame = 0)) {
        #S Mutation
        
        S_sig <- S_sig + 1
        
      } else {
        #NS Mutation
        NS_sig <- NS_sig + 1
      }
      
    }
    
  }
  
}

dNdSexp_sig <- NS_sig/S_sig
dNdS_res <- dNdSobs / dNdSexp_sig

############################################
# SPEARMAN CORRELATION BETWEEN WEEK AND VF #
############################################

VARIANTdfsub <- VARIANTdf[VARIANTdf$MuID %in% MUTATIONdf$MuID[MUTATIONdf$TrueTrans == 'yes'],]

MUTATIONdf$SpearmanWeekVF_rho <- NA
MUTATIONdf$SpearmanWeekVF_pv <- NA
MUTATIONdf$SpearmanWeekVF_pv_adj <- NA

pb <- txtProgressBar(min = 0, max = length(unique(VARIANTdfsub$MuID)), style = 3)
i = 0
for(m in unique(VARIANTdfsub$MuID)) {
  i = i + 1
  setTxtProgressBar(pb, i)
  VF_m <- VARIANTdfsub$VF[VARIANTdfsub$MuID == m]
  week_m <- as.numeric(VARIANTdfsub$Week[VARIANTdfsub$MuID == m])
  
  #Remove VF with NA week
  VF_m <- VF_m[!is.na(week_m)]
  week_m <- week_m[!is.na(week_m)]
  
  if(length(VF_m) < 3) {
    next
  } else {
    Sp_test <- cor.test(x = VF_m, y = week_m, method = "spearman")
    MUTATIONdf$SpearmanWeekVF_rho[MUTATIONdf$MuID == m] <- Sp_test$estimate
    MUTATIONdf$SpearmanWeekVF_pv[MUTATIONdf$MuID == m] <- Sp_test$p.value
  }
  
}

MUTATIONdf$SpearmanWeekVF_pv_adj[!is.na(MUTATIONdf$SpearmanWeekVF_pv)] <- p.adjust(p = MUTATIONdf$SpearmanWeekVF_pv[!is.na(MUTATIONdf$SpearmanWeekVF_pv)], 
                                                                                   method = "fdr")

############################################
#### VF vs COUNT S-NS-NC
mergedVF<-VARIANTdf[,c("MuID","VF","Sample", "SignatureCluster")]

Sample_S1 <- mergedVF$Sample[mergedVF$SignatureCluster == "SC#1" & !is.na(mergedVF$SignatureCluster)]
Sample_S2 <-  mergedVF$Sample[mergedVF$SignatureCluster == "SC#2" & !is.na(mergedVF$SignatureCluster)]
Sample_S3 <-  mergedVF$Sample[mergedVF$SignatureCluster == "SC#3" & !is.na(mergedVF$SignatureCluster)]

bin <- sort(unique(round(var_freq_matrix_filt[var_freq_matrix_filt>thrDETECTION & var_freq_matrix_filt <= thrCLONAL], digits = 2)))

subMtx_S1 <- var_freq_matrix_filt[Sample_S1,]
subMtx_S2 <- var_freq_matrix_filt[Sample_S2,]
subMtx_S3 <- var_freq_matrix_filt[Sample_S3,]

nS1 <- c()
nS2 <- c()
nS3 <- c()

pb <- txtProgressBar(min = 0, max = length(bin), style = 3)
for(i in 2:length(bin)) {
  nS1 <- c(nS1, sum(subMtx_S1 >= bin[1] & subMtx_S1 <= bin[i]))
  nS2 <- c(nS2, sum(subMtx_S2 >= bin[1] & subMtx_S2  <= bin[i]))
  nS3 <- c(nS3, sum(subMtx_S3 >= bin[1] & subMtx_S3 <= bin[i]))
  setTxtProgressBar(pb, i)
}

nS1 <- nS1 / max(nS1)
nS2 <- nS2 / max(nS2)
nS3 <- nS3 / max(nS3)

dfEffect <- data.frame(bin = rep(bin[1:length(nS1)], 3), freq = c(nS1, nS2,nS3), effect = c(rep("SC#1", length(nS1)) , rep("SC#2", length(nS2)),rep("SC#3", length(nS3))))
dfEffect$effect <- factor(as.character(dfEffect$effect), levels = c("SC#1", "SC#2", "SC#3"))

#########################
# ANNOTATION MrBayes tree
#########################

library('pheatmap')
library('ggplot2')
library('RColorBrewer')
library('dplyr')
library('ggforce')


inference_step2 <- readRDS("data/inference_step2.rds")
assignments <- readRDS("data/assignments.rds")

memoSort <- function(M, sortCol = TRUE) {
  if(sortCol == TRUE) {
    geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
  } else {
    geneOrder <- 1:nrow(M)
  }
  
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]) {
        score <- score + 2^(length(x)-i);
      }
    }
    return(score);
  }
  scores <- apply(M[geneOrder, ], 2, scoreCol);
  sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
  M <- M[geneOrder, sampleOrder]
  
  scores <- apply(M, 1, scoreCol);
  mutOrder <- sort(scores, decreasing=FALSE, index.return=TRUE)$ix;
  
  return(M[mutOrder, ]);
}

genotypes <- unique(inference_step2$corrected_genotypes)
rownames(genotypes) <- inference_step2$C[match(rownames(genotypes), rownames(inference_step2$C))]

genotypes <- cbind(genotypes, genotypes[,'28882_G_A|28883_G_C'])
colnames(genotypes)[which(colnames(genotypes) == '28882_G_A|28883_G_C')] <- '28882_G_A'
colnames(genotypes)[ncol(genotypes)] <- '28883_G_C'


genotypesOrdered <- memoSort(genotypes, sortCol = TRUE)

MutAnn <- MUTATIONdf[match(colnames(genotypesOrdered), MUTATIONdf$MuID), c("pos", "Variation", "genomicRegion", "effect", "AAvar")]

MutAnn$MuID <- paste(MutAnn$pos, MutAnn$Variation, sep = ' ')
MutAnn$AAvar <- sapply(MutAnn$AAvar, FUN = function(x){strsplit(x,'\\.')[[1]][2]})

MutAnn <- MutAnn[,c("MuID", "genomicRegion", "effect", "AAvar")]

colnames(MutAnn) <- c("Variation", "Genomic location", "Effect", "Amino Acid Variations")
MutAnn$`Amino Acid Variations`[is.na(MutAnn$`Amino Acid Variations`)] <- '-'

write.table(MutAnn, "data/mutation_annotations.txt", sep = '\t', col.names = T, row.names = F,  quote = F)


pdf(file = "plots/Fig5_B.pdf", width = 10, height = 10)
pheatmap(t(genotypesOrdered), 
         color = colorRampPalette(c("white","#fff08e"))(100),
         border_color = "black",
         cellwidth = 15,
         cellheight = 15,
         
         scale = "none",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         
         legend = TRUE, 
         
         drop_levels = TRUE,
         show_rownames = T,
         show_colnames = T, 
         angle_col = 90,
         
         fontsize = 10,
         
         silent = FALSE,
)
dev.off()

TableCountMut <- as.data.frame(table(statistics$VERSO))


############ VIOLIN VF AMONG CLADES ############
Mut_nP_nL <- data.frame(MuID = setdiff(MUTATIONdf$MuID, colnames(genotypesOrdered)))
Mut_nP_nL$nP <- 0
Mut_nP_nL$nL <- 0

VARIANTdffilt <- VARIANTdf[VARIANTdf$VF <= thrCLONAL,]

for(m in Mut_nP_nL$MuID){
  Mut_nP_nL$nP[Mut_nP_nL$MuID == m] <- length(unique(VARIANTdffilt$Sample[VARIANTdffilt$MuID == m]))
  
  Mut_nP_nL$nL[Mut_nP_nL$MuID == m]  <- length(unique(VARIANTdffilt$VERSO[VARIANTdffilt$MuID == m]))
}
Mut_nP_nL <- Mut_nP_nL[Mut_nP_nL$nP > 0 & Mut_nP_nL$nL > 0,]

Mut_nP_nL$class <- NA
Mut_nP_nL$class <- paste0(Mut_nP_nL$nL, ' clades')
Mut_nP_nL$class[Mut_nP_nL$nL == 1] <- '1 clade'
Mut_nP_nL$class[Mut_nP_nL$nP == 1] <- 'private'
Mut_nP_nL$class[Mut_nP_nL$nL >= 6] <- '>= 6 clades'

VARIANTdfmerged <- merge(x = VARIANTdffilt[VARIANTdffilt$VF <= thrCLONAL,], y = Mut_nP_nL[,c('MuID', 'class')], by = 'MuID', all.x = F)

VARIANTdfmerged$class <- factor(VARIANTdfmerged$class, levels = c("private","1 clade",paste0(2:5, ' clades'), '>= 6 clades'))
VARIANTdfmerged$type <- factor(VARIANTdfmerged$type, levels = c("always minor","mixed","always clonal"))

stat_box_data <- function(y, upper_limit = max(VARIANTdfmerged$VF[!is.na(VARIANTdfmerged$class)], na.rm = T) * 1.1) {
  return(
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('m =', length(y), '\n')
    )
  )
}
pdf(file = "plots/Fig5_D.pdf", width = 14, height = 4)
ggplot(data = VARIANTdfmerged, aes(x = class, fill = class, y = VF)) +
  geom_violin(scale = "width", draw_quantiles = 0.5)+
  stat_summary(
    fun.data = stat_box_data,
    geom = "text",
    hjust = 0.5,
    vjust = 0.5
  )+ 
  scale_y_continuous(breaks = seq(0.0,1,0.2), limits = c(0,1.05)) + 
  xlab("") +
  ylab('Variant Frequency') +
  theme_bw() + 
  scale_fill_manual(values = c("#ffcd30","#c23b22", rev(colorRampPalette(c("gray90","gray40"))(11))))+
  theme(text = element_text(size = 20), legend.position="none")
dev.off()


for(cl in unique(VARIANTdfmerged$class)) {
  
  nS = length(unique(VARIANTdfmerged$Sample[VARIANTdfmerged$class==cl]))
  nM = length(unique(VARIANTdfmerged$MuID[VARIANTdfmerged$class==cl]))
  
  cat(cl, " | sample = ", nS, " | mut = ", nM, '\n')
}
