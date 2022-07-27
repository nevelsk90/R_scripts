#### prepare data ####
### load bulkRNAseq of glvl counts
# Wt1het.del

# Nphs2mut.

# adriamycin

### load DE results
ll <- list("/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_Wt1/DEexonIntron_4w_glvl.rda",
     "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_Wt1/DEexonIntron_12w_glvl.rda",
     "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_NPHS2/Nphs2_IntExon_glvl_4wDE.rda",
     "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_NPHS2/Nphs2_IntExon_glvl_8wDE.rda",
     "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_adriamycin/adriamycin.bulkRNAseq_glvl_DE.rda")
DElist <- lapply( ll, readRDS )
names(DElist) <- c("Wt1het.del_4w","Wt1het.del_12w","Nphs2.mut_4w","Nphs2.mut_8w","adriam_4w")
### compare global results
ggname <- na.omit(Reduce(union, lapply(DElist, function(X) X$gName[X$padj<0.1])))
datt <- Reduce(cbind, lapply(DElist, function(X) X$log2FoldChange[match(ggname, X$gName)]))
colnames(datt) <- names(DElist)
# datt[is.na(datt)] <- 0
corrplot::corrplot(cor(datt, use="pairwise.complete.obs"))

### perform GO and pathway annotation


#### buble plot of GO and pathways ####

#### GSEA analysis ####


#### comparison of PPI (STRING) ####

