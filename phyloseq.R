library("phyloseq"); packageVersion("phyloseq")  ##‘1.34.0’
library("ggplot2"); packageVersion("ggplot2") ##‘3.3.3’
theme_set(theme_bw())
library("plyr")
library("cluster")
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames#packageVersion("cluster")
library("knitr")
library("BiocStyle")

set1 <- phyloseq(OTU, TAX, samples)
set1
saveRDS(set1, "set1.rds")

sample_names(set1)
rank_names(set1)
sample_variables(set1)

##filter by family
rank_names(set1)
table(tax_table(set1)[, "Family"], exclude = NULL)
set1 <- subset_taxa(set1, !is.na(Family) & !Family %in% c("", "uncharacterized"))
prevdf = apply(X = otu_table(set1),
               MARGIN = ifelse(taxa_are_rows(set1), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(set1),
                    tax_table(set1))
plyr::ddply(prevdf, "Family", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
filterFamily = c("UBA5755","UBA3700","UBA1390","UBA1381","Paenibacillaceae","CAG-272","CAG-822","CAG-917")
set1f = subset_taxa(set1, !Family %in% filterFamily)
set1f

table(tax_table(set1f)[, "Phylum"], exclude = NULL)
prevdf = apply(X = otu_table(set1f),
               MARGIN = ifelse(taxa_are_rows(set1f), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(set1f),
                    tax_table(set1f))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
filterPhyla = c("Firmicutes_B", "Fusobacteriota")
ps1 = subset_taxa(set1f, !Phylum %in% filterPhyla)
ps1

prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
postscript("barplot/abund_prev.eps",hor=T)
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps1),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, linetype = 2) +  geom_point(size = 2) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
dev.off()

prevalenceThreshold = 0.1 * nsamples(ps1)
prevalenceThreshold
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps1)
ps2

set1n  = transform_sample_counts(ps2, function(x) x / sum(x) )
set1n

library(phyloseqCompanion)
phyloseq2lefse(ps = set1n, covars = "group", file.name = "lefse_data.txt")

ps2.genus <- tax_glom(ps2, taxrank = "Genus")
set1n.genus <- tax_glom(set1n, taxrank = "Genus")
set1n.family <- tax_glom(set1n, taxrank = "Family")

postscript("barplot/Family_filter.eps",hor=T)
plot_bar(set1n, "Family", fill="Phylum", facet_grid=group~.)
dev.off()

postscript("barplot/Order_filter.eps",hor=T)
plot_bar(set1n, "Order", fill="Phylum", facet_grid=group~.)
dev.off()

#p_Proteobacteria
set1n_proteobacteria <- subset_taxa(set1n, Phylum %in% c("Proteobacteria"))
postscript("barplot/group_proteobacteria.eps",hor=T)
plot_bar(set1n_proteobacteria, "Species", fill="Species", facet_grid=group~., title="Proteobacteria")
dev.off()

##p_Bacteroidota
set1n_Bacteroidota <- subset_taxa(set1n, Phylum %in% c("Bacteroidota"))
postscript("barplot/group_Bacteroidota_family.eps",hor=T)
plot_bar(set1n_Bacteroidota, x="Family", fill = "Family", facet_grid = group~., title="Bacteroidota") 
dev.off()

set1n_Bacteroides <- subset_taxa(set1n, Genus %in% c("Bacteroides"))
postscript("barplot/group_Bacteroides-species.eps",hor=T)
plot_bar(set1n_Bacteroides, x="Species", fill = "Genus", facet_grid = group~., title="Bacteroides") 
dev.off()

##p_Actinobacteriota
set1n_Actinobacteriota <- subset_taxa(set1n, Phylum %in% c("Actinobacteriota"))
postscript("barplot/group_Actinobacteriota.eps",hor=T)
plot_bar(set1n_Actinobacteriota, x="Species", fill = "Species", facet_grid = group~., title="Actinobacteriota") 
#  geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
dev.off()

##p_Firmicutes
set1n_Firmicutes <- subset_taxa(set1n, Phylum %in% c("Firmicutes"))
postscript("barplot/group_Firmicutes_family.eps",hor=T)
plot_bar(set1n_Firmicutes, x="Family", fill = "Order", facet_grid = group~., title="Firmicutes") 
#  geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
plot_bar(set1n_Firmicutes, x="Order", fill = "Order", facet_grid = group~., title="Firmicutes")
plot_bar(set1n_Firmicutes, x="Genus", fill = "Genus", facet_grid = group~., title="Firmicutes") #  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
dev.off()

#Lachnospiraceae
set1n_Lachnospiraceae <- subset_taxa(set1n, Family %in% c("Lachnospiraceae"))
postscript("barplot/group_Lachnospiraceae_species.eps",hor=T)
plot_bar(set1n_Lachnospiraceae, x="Species", fill = "Genus", facet_grid = group~., title="Lachnospiraceae") 
dev.off()

#Lactobacillaceae
set1n_Lactobacillaceae <- subset_taxa(set1n, Family %in% c("Lactobacillaceae"))
postscript("barplot/group_Lactobacillales_species.eps",hor=T)
plot_bar(set1n_Lactobacillaceae, x="Species", fill = "Genus", facet_grid = group~., title="Lactobacillaceae") 
dev.off()

#Clostridiaceae
set1n_clostridiaceae <- subset_taxa(set1n, Family %in% c("Clostridiaceae"))
postscript("barplot/group_clostridiaceae_species.eps",hor=T)
plot_bar(set1n_clostridiaceae, x="Species", fill = "Genus", facet_grid = group~., title="clostridiaceae") 
dev.off()

#erysipelotrichaceae
set1n_erysipelotrichaceae <- subset_taxa(set1n, Family %in% c("Erysipelotrichaceae"))
postscript("barplot/group_erysipelotrichaceae_species.eps",hor=T)
plot_bar(set1n_erysipelotrichaceae, x="Species", fill = "Genus", facet_grid = group~., title="clostridiaceae") 
dev.off()

#oscillospiraceae
set1n_oscillospiraceae <- subset_taxa(set1n, Family %in% c("Oscillospiraceae"))
postscript("barplot/group_oscillospiraceae_species.eps",hor=T)
plot_bar(set1n_oscillospiraceae, x="Species", fill = "Genus", facet_grid = group~., title="oscillospiraceae") 
dev.off()

#Peptostreptococcaceae
set1n_peptostreptococcaceae <- subset_taxa(set1n, Family %in% c("Peptostreptococcaceae"))
postscript("barplot/group_peptostreptococcaceae_species.eps",hor=T)
plot_bar(set1n_peptostreptococcaceae, x="Species", fill = "Genus", facet_grid = group~., title="peptostreptococcaceae") 
dev.off()

#ruminococcaceae
set1n_ruminococcaceae <- subset_taxa(set1n, Family %in% c("Ruminococcaceae"))
postscript("barplot/group_ruminococcaceae_species.eps",hor=T)
plot_bar(set1n_ruminococcaceae, x="Species", fill = "Genus", facet_grid = group~., title="ruminococcaceae") 
dev.off()

##p_Verrucomicrobiota
set1n_Verrucomicrobiota <- subset_taxa(set1n, Phylum %in% c("Verrucomicrobiota"))
postscript("barplot/group_Verrucomicrobiota.eps",hor=T)
plot_bar(set1n_Verrucomicrobiota, x="Species", fill = "Species", facet_grid = group~., title="Verrucomicrobiota") 
#  geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")
#plot_bar(set1n_Firmicutes, x="Genus", fill = "Genus", facet_grid = group~., title="Firmicutes") +
#  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
dev.off()

###alpha
alpha_meas = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson")(p <- plot_richness(ps2, "group", "group", measures=alpha_meas))
postscript("alpha/diversity_all.eps",hor=T)
p + geom_boxplot(data=p$data, aes(x=group, y=value, color=NULL), alpha=0.1)
dev.off()

library("ggpubr")
postscript("alpha/diversity_group_significance_chao1.eps",hor=T)
a_my_comparisons <- list( c("Abx", "Abx_Tac"), c("Ctl", "Ctl_Tac"), c("Abx", "Ctl"), c("Abx", "Ctl_Tac"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
plot_richness(set1, x="group", measures="Chao1", color = "group")+
  geom_boxplot() + 
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) +
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args) 
#  facet_wrap(~day)
dev.off()

###beta
set1n.ord <- ordinate(set1n, "NMDS", "bray")
set1.ord <- ordinate(set1, "NMDS", "bray")
p2 = plot_ordination(set1n, set1n.ord, type="samples", color="group", shape="group") 
p2 + geom_polygon(aes(fill=group)) + geom_point(size=5) 
dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist){
  ordi = ordinate(physeq, method=i, distance=dist)
  plot_ordination(physeq, ordi, "samples", color="group")
}, set1n, dist)
names(plist) <- ord_meths
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
p = ggplot(pdataframe, aes(Axis_1, Axis_2, color=group, shape=group, fill=group))
p = p + geom_point(size=4) + geom_polygon()
p = p + facet_wrap(~method, scales="free")
p = p + scale_fill_brewer(type="qual", palette="Set1")
p = p + scale_colour_brewer(type="qual", palette="Set1")
postscript("beta/bray6.eps",hor=T)
p
dev.off()
postscript("beta/bray_CCA.eps",hor=T)
p = plist[[2]] + scale_colour_brewer(type="qual", palette="Set1")
p = p + scale_fill_brewer(type="qual", palette="Set1")
p = p + geom_point(size=5) + geom_polygon(aes(fill=group))
p
dev.off()

set1.ord <- ordinate(set1n, "PCoA", "jsd")
postscript("beta/PCoA_jsd.eps",hor=T)
plot_ordination(set1, set1.ord, type="samples", color="group", 
                shape="group", title="Samples") + geom_point(size=3)
dev.off()

###

