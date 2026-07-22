## conda activate r_trajectory
library(UpSetR)
allDBs<-read.csv('../csv/nichesphereDB_pmid.csv')
allDBs$Ligand<-tolower(allDBs$Ligand)
listInput <- list(Collagens = allDBs$Ligand[which(allDBs$category=='Collagens')], 
                  Chemokine = allDBs$Ligand[which(allDBs$category=='Chemokine')], 
                  SecretedFactors = allDBs$Ligand[which(allDBs$category=='SecretedFactors')],
                  Cytokine = allDBs$Ligand[which(allDBs$category=='Cytokine')],
                  ECMaffiliated = allDBs$Ligand[which(allDBs$category=='ECMaffiliated')],
                  ECMglycoprots = allDBs$Ligand[which(allDBs$category=='ECMglycoprots')],
                  ECMregulators = allDBs$Ligand[which(allDBs$category=='ECMregulators')],
                  Proteoglycans = allDBs$Ligand[which(allDBs$category=='Proteoglycans')],
                  GrowthFactor = allDBs$Ligand[which(allDBs$category=='GrowthFactor')],
                  Inhibitory = allDBs$Ligand[which(allDBs$category=='Inhibitory')]
                  #Phagocytosis = allDBs$Ligand[which(allDBs$category=='phagocytosis')],
                  #Apoptosis = allDBs$Ligand[which(allDBs$category=='apoptosis')]
                  )

pdf('../figures/revision/DBsUpset_new.pdf', width=6, height=4)
upset(fromList(listInput), order.by = "freq", nintersects=NA, nsets=10)
dev.off()


