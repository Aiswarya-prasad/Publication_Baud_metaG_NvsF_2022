library(data.table)
library(ggplot2)
library(ape)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(ggrepel)
library(qvalue)
library(vegan)
library(tidyverse)
library(ggnewscale)
library(ggrepel)
library(vegan)
library(ggpubr)

# This function is from https://github.com/Russel88/MicEco/blob/91f8e6f5d67e0bfd018dc3b46da3994b1eadeb46/R/adonis_OmegaSq.R#L10
#' Calculate (partial) Omega-squared (effect-size calculation) for PERMANOVA and add it to the input object
#'
#' @param adonisOutput An adonis object
#' @param partial Should partial omega-squared be calculated (sample size adjusted). Default TRUE
#' @return Original adonis object with the (partial) Omega-squared values added
#' @import vegan
#' @export
adonis_OmegaSq <- function(adonisOutput, partial = TRUE){
    if(!(is(adonisOutput, "adonis") || is(adonisOutput, "anova.cca")))
        stop("Input should be an adonis object")
    if (is(adonisOutput, "anova.cca")) {
        aov_tab <- adonisOutput
        aov_tab$MeanSqs <- aov_tab$SumOfSqs / aov_tab$Df
        aov_tab$MeanSqs[length(aov_tab$Df)] <- NA
    } else {
        aov_tab <- adonisOutput$aov.tab
    }
    heading <- attr(aov_tab, "heading")
    MS_res <- aov_tab[pmatch("Residual", rownames(aov_tab)), "MeanSqs"]
    SS_tot <- aov_tab[rownames(aov_tab) == "Total", "SumsOfSqs"]
    N <- aov_tab[rownames(aov_tab) == "Total", "Df"] + 1
    if(partial){
        omega <- apply(aov_tab, 1, function(x) (x["Df"]*(x["MeanSqs"]-MS_res))/(x["Df"]*x["MeanSqs"]+(N-x["Df"])*MS_res))
        aov_tab$parOmegaSq <- c(omega[1:(length(omega)-2)], NA, NA)
    } else {
        omega <- apply(aov_tab, 1, function(x) (x["SumsOfSqs"]-x["Df"]*MS_res)/(SS_tot+MS_res))
        aov_tab$OmegaSq <- c(omega[1:(length(omega)-2)], NA, NA)
    }
    if (is(adonisOutput, "adonis"))
        cn_order <- c("Df", "SumsOfSqs", "MeanSqs", "F.Model", "R2",
                      if (partial) "parOmegaSq" else "OmegaSq", "Pr(>F)")
    else
        cn_order <- c("Df", "SumOfSqs", "F", if (partial) "parOmegaSq" else "OmegaSq",
                      "Pr(>F)")
    aov_tab <- aov_tab[, cn_order]
    attr(aov_tab, "names") <- cn_order
    attr(aov_tab, "heading") <- heading
    if (is(adonisOutput, "adonis"))
        adonisOutput$aov.tab <- aov_tab
    else
        adonisOutput <- aov_tab
    return(adonisOutput)
}

source("utilities.R")

sample_order <- c('N01','N02','N03','N04','N06','N07','N08','N09','N10','N11','N12','N13','N14','N15','N16',
                  'F01','F02','F03','F04','F06','F07','F08','F09','F10','F11','F12','F13','F14','F15','F16')
phylo_list <- c("api","bapis","bifido","bom","com","firm4","firm5","fper","gilli","lkun","snod")
phylo_order <- c("firm5","bifido","firm4","gilli","snod","api","fper","bapis","com","lkun","bom")
phylo_order_rv <- rev(phylo_order)
treatments <- c("Nurses","Foragers")
#Set the colors.
set1 <- brewer.pal(9,"Set1")
phylo_col <- c(set1, "light grey", "#DFFF00")
phylo_col_rv <- rev(phylo_col)
treat_colors <- c("#D55E00","#56B4E9","#F0E442")
samp_dt <- fread("SamplingGilles2019.csv")
samp_dt <- samp_dt[, c("Sample", "Sample_type", "Location", "Number_guts", "Average_gut_mass")]

### Gut Weight comparisons

ggplot(samp_dt,
       aes(x=factor(Sample_type,levels=treatments),
           y = Average_gut_mass,
           fill=factor(Sample_type,levels=treatments)))+
  geom_boxplot()+
  scale_fill_manual(name="Enrichment", values=c("Nurses"="#D55E00",
                                                "Foragers"="#56B4E9"))+
  geom_point(aes(fill=factor(Sample_type,levels=treatments)),
             shape=22,
             size = 6,
             position=position_dodge2(width=0.6))+
  labs(subtitle=paste0("Wilcoxon sign rank test: \n P = ",
                       wilcox.test(samp_dt[Sample %like% "N", Average_gut_mass],
                                   samp_dt[Sample %like% "F", Average_gut_mass])$p.value),
       fill="Host")+
  ylab("Average gut weight (g)")+
  scale_y_continuous(trans = "log10")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x=element_blank(),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))
ggsave("figures/gutweight_plot.pdf", width = 10, height = 12)
  
### Bacterial loads comparisons

#qPCR standard curve values from Kesnerova et al. (2017)
actin_intercept = 35.0942119870867
actin_slope = -3.2699442250388
UV_intercept = 36.5821936471122
UV_slope = -3.35085896083287

filepath <- paste0("20190813_AllHives.csv")
CT_dt <- fread(filepath, skip="Well Position")
CT_dt <- fread(filepath, skip="Well Position", nrows = CT_dt[,sum(CT != '')])
CT_dt <- merge(CT_dt[, mean(as.numeric(CT), na.rm=TRUE), keyby=c("Sample Name", "Target Name")], 
               CT_dt[, sd(as.numeric(CT), na.rm=TRUE), keyby=c("Sample Name", "Target Name")])
setnames(CT_dt, c("Sample", "Target", "AverageCT", "SD"))
CT_dt <- CT_dt[ !(Sample %like% "[FN]05" ), ]
CT_dt <- CT_dt[ !(Sample %like% "control" ), ]

CT_dt[ Target %like% "Actin", copies := 10^((AverageCT-actin_intercept)/actin_slope)]
CT_dt[ Target %like% "UV", copies := 10^((AverageCT-UV_intercept)/UV_slope)]

CT_dt[ Sample %like% "N", Host := "Nurses" ]
CT_dt[ Sample %like% "F", Host := "Foragers" ]

sampling_dt <- fread("SamplingGilles2019.csv")
CT_dt <- merge.data.table(CT_dt, sampling_dt[,c("Sample", "DNA_yield")], by = "Sample", all.x = TRUE)
CT_dt[, DNA := copies*DNA_yield/10]
CT_dt[, Hive:=gsub("[NF]", "", Sample)]
CT_dt[, copiesperngDNA:=copies/10]
CT_ratio_dt <- merge.data.table(CT_dt[ Target %like% "Actin",],
                                CT_dt[ Target %like% "UV",],
                                by=c("Sample"))
CT_ratio_dt[, ratio := copies.y / copies.x]
med_act_copies <- median(CT_dt[ Target %like% "Actin", copies])
CT_ratio_dt[ , Norm_copies := ratio * med_act_copies ]
CT_ratio_dt <- merge.data.table(CT_ratio_dt, CT_dt[Target %like% "UV", c("Sample", "copiesperngDNA")], by="Sample", all.x = TRUE)

CT_dt[ Target %like% "UV"]
# calculate median of F and N
median_N <- median(CT_dt[ Target %like% "UV" & Sample %like% "N", copiesperngDNA])
median_F <- median(CT_dt[ Target %like% "UV" & Sample %like% "F", copiesperngDNA])
copies_diff <- median_N / median_F
median_F_size <- median(samp_dt[ Sample_type %like% "F", Average_gut_mass])
median_N_size <- median(samp_dt[ Sample_type %like% "N", Average_gut_mass])
size_diff <- median_N_size / median_F_size
copies_diff
size_diff

ggplot(CT_dt[ Target %like% "UV"],
       aes(x=factor(Host,levels=treatments),
           y=copiesperngDNA,
           fill=factor(Host,levels=treatments)))+
  geom_boxplot()+
  scale_fill_manual(name="Host", values=c("Nurses"="#D55E00",
                                                "Foragers"="#56B4E9"))+
  geom_point(aes(fill=factor(Host,levels=treatments),
                 shape=factor(Host,levels=treatments)), 
             size=6,
             position=position_dodge2(width=0.6)) +
  labs(subtitle = paste0("Wilcoxon sign rank test: \n p = ", 
                         wilcox.test(CT_dt[ Target %like% "UV" & Sample %like% "N", copiesperngDNA], 
                                     CT_dt[ Target %like% "UV" & Sample %like% "F", copiesperngDNA],
                                     paired=TRUE)$p.value),
       fill="Host",
       shape="Host")+
  ylab("16S copy number per ng of DNA")+
  scale_y_continuous()+
  scale_shape_manual(values=c(21,24))+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill = "#f2f2f2", colour = NA),
        panel.grid.major=element_line(color="white"),
        panel.grid.minor=element_line(color="white"),
        plot.subtitle=element_text(hjust=0.5),
        axis.title.x=element_blank())
  ggsave("figures/Fig1A-16S_copy_number_per_ng_DNA.pdf", width = 10, height = 12)


ggplot(CT_dt[ Target %like% "Actin"], 
                        aes(x=factor(Host,levels=treatments), 
                            y=copies, 
                            fill=factor(Host,levels=treatments)))+
  geom_boxplot()+
  scale_fill_manual(name="Host", values=c("Nurses"="#D55E00",
                                          "Foragers"="#56B4E9"))+
  geom_point(aes(fill=factor(Host,levels=treatments)), 
             shape=22,
             size=6)+
  geom_line(aes(group=Hive)) +
  labs(subtitle = paste0("Wilcoxon sign rank test: \n P = ", 
                         wilcox.test(CT_dt[ Target %like% "Actin" & Sample %like% "N", copies], 
                                     CT_dt[ Target %like% "Actin" & Sample %like% "F", copies],
                                     paired=TRUE)$p.value),
       fill="Sample type")+
  ylab("Actin copy number")+
  scale_y_continuous()+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill = "#f2f2f2", colour = NA),
        panel.grid.major=element_line(color="white"),
        panel.grid.minor=element_line(color="white"),
        plot.subtitle=element_text(hjust=0.5))
  ggsave("figures/Actin_copy_number.pdf", width = 10, height = 12)

  
uv_DNA_bplot <- ggplot(CT_dt[ Target %like% "UV"], 
                        aes(x=factor(Host,levels=treatments), 
                            y=DNA, 
                            fill=factor(Host,levels=treatments)))+
  geom_boxplot()+
  scale_fill_manual(name="Host", values=c("Nurses"="#D55E00",
                                          "Foragers"="#56B4E9"))+
  geom_point(aes(fill=factor(Host,levels=treatments)), 
             shape=22,
             size=6)+
  geom_line(aes(group=Hive)) +
  labs(subtitle = paste0("Wilcoxon sign rank test: \n P = ", 
                         wilcox.test(CT_dt[ Target %like% "UV" & Sample %like% "N", DNA], 
                                     CT_dt[ Target %like% "UV" & Sample %like% "F", DNA],
                                     paired=TRUE)$p.value),
       fill="Sample type")+
  ylab("Total estimated 16S gene copies in sample")+
  scale_y_continuous()+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))

CT_ratio_dt <- merge.data.table(CT_dt[ Target %like% "Actin",],
                                CT_dt[ Target %like% "UV",],
                                by=c("Sample"))
CT_ratio_dt[, ratio := copies.y / copies.x]
#CT_ratio_dt[, DNAratio := DNA.y / DNA.x]
#CT_ratio_dt <- CT_ratio_dt[, c("Sample", "ratio")]

ggplot(CT_ratio_dt, 
       aes(x=factor(Host.y, levels=treatments),
           y=ratio,
           fill=factor(Host.y, levels=treatments)))+
  geom_boxplot()+
  scale_fill_manual(name="Host", values=c("Nurses"="#D55E00",
                                          "Foragers"="#56B4E9"))+
  geom_point(aes(fill=factor(Host.y, levels=treatments)), 
             shape=22,
             size=6,
             position=position_dodge2(width=0.6)) +
  #geom_line(aes(group = Hive)) +
  labs(subtitle = paste0("Wilcoxon sign rank test: \n P = ",
                         wilcox.test(CT_ratio_dt[Sample %like% "N", ratio],
                                     CT_ratio_dt[Sample %like% "F", ratio])$p.value),
       fill="Sample type")+
  scale_y_continuous(name="16S / actin gene copy number ratio",trans="log10")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x = element_blank(),
        text=element_text(size=22),
        panel.background = element_rect(fill = "#f2f2f2", colour = NA),
        panel.grid.major=element_line(color="white"),
        panel.grid.minor=element_line(color="white"),
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))
    ggsave("figures/16S_actin_ratio.pdf", width = 10, height = 12)

med_act_copies <- median(CT_dt[ Target %like% "Actin", copies])
CT_ratio_dt[ , Norm_copies := ratio * med_act_copies ]
CT_ratio_dt <- merge.data.table(CT_ratio_dt, CT_dt[Target %like% "UV", c("Sample", "copiesperngDNA")], by="Sample", all.x = TRUE)


### Read mapping proportions
#Read data
readnum_dt <- fread("readnumbers.csv")
setnames(readnum_dt, c("Sample","Raw","Bacterial_db","A_mellifera"))
readnum_dt[, Leftover := Raw-Bacterial_db-A_mellifera]
readnum_dt[, prop_bac := Bacterial_db/Raw]
readnum_dt[, prop_amel := A_mellifera/Raw]
reads_abs_stat <- wilcox.test(readnum_dt[Sample %like% "N", prop_bac], 
                              readnum_dt[Sample %like% "F", prop_bac], paired = TRUE)
readnum_mlt <- melt.data.table(readnum_dt, id.vars = "Sample", value.name = "Reads", variable.name = "Mapping")[!(Mapping %like% "prop" | Mapping %like% "Raw"),]

ggplot(readnum_mlt[ Mapping!="Raw",],
       aes(x=factor(Sample,
                    levels=sample_order),
           y=Reads,
           fill=factor(Mapping,
                       levels=c("Leftover","A_mellifera","Bacterial_db"),
                       labels=c("Not mapped","A. mellifera genome","Bacterial database")))) +
  geom_bar(stat="identity")+
  xlab("Sample")+
  ylab("Number of reads")+
  scale_fill_discrete(name="Mapping")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5))
 ggsave("figures/readnumbers.pdf", width = 10, height = 12)

read_prop_mlt <- melt.data.table(readnum_dt, id.vars = "Sample", value.name = "Proportion", variable.name = "Mapping")[Mapping %like% "prop",]

ggplot(read_prop_mlt,
       aes(x=factor(Sample,
                    levels=sample_order),
           y=Proportion,
           fill=factor(Mapping,
                       levels=c("left_prop","prop_amel","prop_bac"),
                       labels=c("Not mapped","A. mellifera genome","Bacterial database")))) +
  geom_bar(stat="identity") +
  xlab("Sample")+
  ylab("Read proportion")+
  scale_fill_discrete(name="Mapping")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5))
      ggsave("figures/readproportions.pdf", width = 10, height = 12)

read_hostbac_ratio <- copy(readnum_dt)
read_hostbac_ratio[, readratio:=Bacterial_db/A_mellifera]
tmp_ct_reads <- merge.data.table(CT_ratio_dt, read_hostbac_ratio[, c("Sample", "readratio")], by = "Sample")
tmp_ct_reads[ Sample %like% "N", Host := "Nurses" ]
tmp_ct_reads[ Sample %like% "F", Host := "Foragers" ]

ggplot(tmp_ct_reads, aes(x=ratio, y=readratio, fill=factor(Host, levels=treatments)))+
  geom_point(aes(fill = factor(Host, levels = treatments)), 
             shape = 22, size = 6, position = position_dodge2(width = 0.6))+
  labs(subtitle = paste0("Pearson correlation coefficient: \n r = ", 
                         cor(tmp_ct_reads[, ratio], tmp_ct_reads[, readratio], method = "pearson")),
       fill="Host")+
  scale_fill_manual(name="Host", values=c("Nurses"="#D55E00",
                                          "Foragers"="#56B4E9"))+
  xlab("16S/actin gene copy number ratio")+
  ylab("bacterial/host read number")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))


### Phylotype relative abundances
#Read data
phylo_all_dt <- data.table()
for (i in phylo_list) {
  tmp_dt <- fread(paste0("community_quantification/", i, "_corecov_coord.txt"), header = TRUE)
  tmp_all_dt <- as.data.table(aggregate(data=tmp_dt, Cov_ter ~ Sample + Cluster, sum))
  tmp_all_dt[, Phylo := i]
  phylo_all_dt <- rbind(phylo_all_dt, tmp_all_dt)
}

#Compute the proportions
cov_sum_dt <- phylo_all_dt[, sum(Cov_ter), by=c("Phylo", "Sample")]
phylo_all_dt <- merge.data.table(phylo_all_dt, cov_sum_dt, by=c("Phylo", "Sample"), all.x = TRUE)
cov_sum_sample_dt <- phylo_all_dt[, sum(Cov_ter), by=Sample]
phylo_all_dt <- merge.data.table(phylo_all_dt, cov_sum_sample_dt, by="Sample", all.x = TRUE)
setnames(phylo_all_dt, c("Sample", "Phylo", "SDP", "Cov_ter", "Phylo_sum", "Sample_sum"))
phylo_all_dt[ , Prop := Cov_ter / Sample_sum ]
phylo_all_dt[ , Prop_SDP := Cov_ter / Phylo_sum]
phylo_all_dt[is.na(Prop_SDP), Prop_SDP:=0]
phylo_all_dt[ , Norm_abun_SDP := setkey(phylo_all_dt, Sample)[ CT_ratio_dt, Prop * copiesperngDNA]]

ggplot(data=phylo_all_dt,
       aes(x=factor(Sample,levels=sample_order),
           y=Prop,
           fill=factor(Phylo,
                       levels=phylo_order_rv,
                       labels=c("Bombella","L. kunkeii","Commensalibacter",
                                "Bartonella","Frischella","Apibacter",
                                "Snodgrasella","Gilliamella","Firm4","Bifidobacterium", "Firm5")))) + 
  geom_bar(stat="identity",colour="white") +
  scale_fill_manual(values=phylo_col_rv, name="Phylotype") +
  xlab("Sample")+
  ylab("Proportion")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title=element_blank(),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5))
ggsave("figures/Fig1B-phylorelabundstackbarplot.pdf", width = 10, height = 12)

#Normalized absolute abundance
phylo_absabun_dt <- as.data.table(aggregate(data=phylo_all_dt, Norm_abun_SDP ~ Sample + Phylo, sum))
setnames(phylo_absabun_dt, c("Sample", "Phylo", "Norm_abun"))
phylo_absabun_dt[ Sample %like% "F[0-9][0-9]", Host:="Foragers"]
phylo_absabun_dt[ Sample %like% "N[0-9][0-9]", Host:="Nurses"]
stats_phylo_absabun <- data.table(phylo_list, 0)
setnames(stats_phylo_absabun, c("Phylo", "p"))
for (i in phylo_list) {
  temp_test <- wilcox.test(phylo_absabun_dt[(Phylo==i & Sample %like% "N[0-9][0-9]"), Norm_abun],
                           phylo_absabun_dt[(Phylo==i & Sample %like% "F[0-9][0-9]"), Norm_abun],
                           paired = TRUE)
  stats_phylo_absabun[ Phylo==i, p:=temp_test$p.value]
}
stats_phylo_absabun[, q:=qvalue(stats_phylo_absabun$p)$qvalues]

phylo_absabun_bplots <- list()
for (i in phylo_list) {
  phylo_absabun_bplots[[i]] <- ggplot(phylo_absabun_dt[ Phylo==i, ],
                                      aes(x=factor(Host,levels=treatments),
                                          y=Norm_abun,
                                          fill=factor(Host,levels=treatments)))+
    geom_boxplot(width=0.5)+
    geom_point(aes(fill = Host),
               shape=22,
               size=4,
               position=position_dodge2(width=0.6)) +
    scale_fill_manual(name="Host", values=c("Nurses"="#D55E00",
                                            "Foragers"="#56B4E9"))+
    ggtitle(i)+
    labs(fill="Host",subtitle=paste0("qvalue = ",stats_phylo_absabun[Phylo==i, q]))+
    scale_y_continuous("Abundance \n (normalized estimated \n 16S gene copies)",
                       trans="sqrt")+
    theme(axis.text.x=element_text(angle=45,hjust=1),
          text=element_text(size=22),
          axis.title.x=element_blank(),
          plot.title=element_text(hjust=0.5), 
          plot.subtitle=element_text(hjust=0.5))
}
pdf(paste0("figures/SuppFig4-phylorabsabundboxplot.pdf"), onefile=TRUE, width=12)
marrangeGrob(phylo_absabun_bplots, nrow = 1, ncol = 1)
dev.off()

### SDP relative abundances and Nurses vs. Foragers normalized abundances boxplots
stats_SDP_absabun <- data.table(unique(phylo_all_dt[, SDP]), 0)
setnames(stats_SDP_absabun, c("SDP", "p"))
for (i in unique(phylo_all_dt[, SDP])) {
  temp_test <- wilcox.test(phylo_all_dt[(SDP==i & Sample %like% "N[0-9][0-9]"), Norm_abun_SDP],
                           phylo_all_dt[(SDP==i & Sample %like% "F[0-9][0-9]"), Norm_abun_SDP],
                           paired = TRUE)
  stats_SDP_absabun[ SDP==i, p:=temp_test$p.value]
}
stats_SDP_absabun[, q:=qvalue(stats_SDP_absabun$p)$qvalues]
phylo_all_dt[ Sample %like% "F[0-9][0-9]", Host:="Foragers"]
phylo_all_dt[ Sample %like% "N[0-9][0-9]", Host:="Nurses"]
sdp_relabun_barplot <- list()
sdp_absabun_boxplot <- list()
for (i in phylo_list) {
  #Relative abundances stuff
  sdp_relabun_barplot[[i]] <- ggplot(data=phylo_all_dt[Phylo==i,],
                                     aes(x=factor(Sample,levels=sample_order),
                                         y=Prop_SDP,
                                         fill=SDP)) + 
    geom_bar(stat="identity", colour="white") +
    scale_fill_manual(values=phylo_col_rv, name = "SDP") + 
    ggtitle(paste0(i," SDP relative abundances")) +
    xlab("Sample")+
    ylab("Proportion")+
    theme(axis.text.x=element_text(angle=45,hjust=1),
          text=element_text(size=22),
          plot.title=element_text(hjust=0.5))
}
pdf(paste0("figures/SuppFig3-pre-sdprelabunstackbarplot.pdf"), onefile=TRUE, width = 12)
marrangeGrob(sdp_relabun_barplot, nrow = 1, ncol = 1)
dev.off()


stats_SDP_absabun <- merge.data.table(stats_SDP_absabun,
                                      phylo_all_dt[Sample=="F01", c("SDP","Phylo")], 
                                      by="SDP", all.x=TRUE )

pdf(paste0("figures/SuppFig4-pre-sdpabsabunboxplot.pdf"), onefile=TRUE, width = 25)
ggplot(phylo_all_dt[!(SDP %in% c("api_apis_dorsa", "api_bombus", "bifido_1_cerana", 
                                 "bifido_bombus","bom_1", "bom_apis_melli","bom_bombus",
                                 "com_drosophila", "com_monarch", "firm5_bombus",
                                 "gilli_4", "gilli_5", "gilli_6", "gilli_apis_andre",
                                 "gilli_apis_dorsa", "gilli_bombus","snod_2", "snod_bombus")),], 
       aes(x=SDP,
           y=Norm_abun_SDP))+
  geom_boxplot(aes(fill=factor(Host, levels=treatments)))+
  geom_point(data=stats_SDP_absabun[q<0.05 & !(SDP %in% c("api_apis_dorsa", "api_bombus", "bifido_1_cerana", 
                                                          "bifido_bombus","bom_1", "bom_apis_melli","bom_bombus",
                                                          "com_drosophila", "com_monarch", "firm5_bombus",
                                                          "gilli_4", "gilli_5", "gilli_6", "gilli_apis_andre",
                                                          "gilli_apis_dorsa", "gilli_bombus","snod_2", "snod_bombus")), ], 
             aes(x=SDP, y=135000), shape="*", size=8, show.legend = FALSE)+
  xlab("SDP")+
  scale_fill_manual(name="Host", values=c("Nurses"="#D55E00",
                                          "Foragers"="#56B4E9"))+
  facet_grid(. ~ Phylo, scales="free", space="free")+
  scale_y_continuous(name="Abundance \n (normalized estimated \n 16S gene copies)", trans="sqrt")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5))
dev.off()

phylo_all_dt[, Hive:=gsub("[FN]","", Sample)]
sdp_fc_dt <- dcast.data.table(phylo_all_dt, Hive + SDP + Phylo ~ Host, value.var = "Norm_abun_SDP", )
sdp_fc_dt[Nurses>0, logtwoFC:=log2(Foragers/Nurses)]

sdp_fc_dt_bis <- copy(sdp_fc_dt)
sdp_fc_dt_bis <- sdp_fc_dt_bis[ !(SDP %in% c("api_apis_dorsa", "api_bombus", "bifido_1_cerana", 
                            "bifido_bombus", "bom_apis_melli","bom_bombus",
                            "com_drosophila", "com_monarch", "firm5_bombus",
                            "gilli_4", "gilli_5", "gilli_6", "gilli_apis_andre",
                            "gilli_apis_dorsa", "gilli_bombus", "snod_bombus")),]

pdf(paste0("figures/Fig1E-fc_sdp_plot.pdf"), onefile=TRUE, width=15)
ggplot(sdp_fc_dt_bis[!(SDP %in% c("api_1", "bom_1", "firm5_7", "snod_2") ),], 
       aes(x = SDP, 
           y = logtwoFC))+
  geom_point(aes(fill=factor(Phylo,
                             levels=phylo_order_rv,
                             labels = c("Bombella", "L. kunkeii", "Commensalibacter",
                                        "Bartonella", "Frischella", "Apibacter",
                                        "Snodgrasella", "Gilliamella", "Firm4", "Bifidobacterium", "Firm5"))), 
             shape = 22, size = 3, position = position_dodge2(width = 0.2)) +
  geom_hline(yintercept = 0)+
  geom_crossbar(data=sdp_fc_dt_bis[!is.na(logtwoFC) & !(SDP %in% c("api_1", "bom_1", "firm5_7", "snod_2")), median(logtwoFC), by=SDP],
                aes(x=SDP,ymin=V1, ymax=V1,y=V1,group=SDP), width = 0.5)+
  geom_point(data=stats_SDP_absabun[q<0.05 & !(SDP %in% c("api_1", "bom_1", "firm5_7", "snod_2")), ], 
             aes(x=SDP, y=3), shape="*", size=8, show.legend = FALSE)+
  scale_fill_manual(values=phylo_col_rv[c(2,3,4,5,7,8,9,10,11)], name = "Phylotype") + 
  xlab("SDP")+
  scale_y_continuous("log2 fold change to nurse sample")+
  theme(axis.text.x=element_text(angle=45,hjust=1), 
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5), 
        plot.subtitle=element_text(hjust=0.5))
dev.off()

### NMDS plots for the relative and absolute compositions of samples, based on phylotypes and based on SDPs, using Bray-Curtis dissimilarity indices.

nmds_plot <- function(data_dt, plot_title, info_df){
  set.seed(123)
  nmds = metaMDS(data_dt, distance = "bray", trymax = 100000)
  nmds_scores = as.data.table(scores(nmds)$sites)
  nmds_scores[, Sample := info_df$Sample]
  nmds_scores[Sample %like% "N[0-9][0-9]", Host := "Nurses" ]
  nmds_scores[Sample %like% "F[0-9][0-9]", Host := "Foragers" ]
  nmds_scores[Sample %like% "N[0-9][0-9]", Host := "Nurses" ]
  nmds_scores[Sample %like% "F[0-9][0-9]", Host := "Foragers" ]
  nmds_plot <- ggplot(nmds_scores, aes(x=NMDS1, y=NMDS2, 
                                     color=factor(Host,levels=treatments),
                                     shape=factor(Host)))+
    geom_point(stat="identity",size=6, alpha = 0.75)+
    scale_color_manual(values=treat_colors)+
    scale_shape_manual(values=c("Foragers" = 17, "Nurses" = 19))+
    labs(x="NMDS1",
         y="NMDS2",
         title=plot_title,
         subtitle=paste0("Stress = ",nmds$stress),
         shape="Host",
         color="Location")+
    theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=22),
          plot.subtitle = element_text(hjust = 0.5)) +
          make_theme(setCol = F)
  return(nmds_plot)
}

nmds_plot_loc <- function(data_dt, plot_title, info_df){
  set.seed(123)
  nmds = metaMDS(data_dt, distance = "bray", trymax = 100000)
  nmds_scores = as.data.table(scores(nmds)$sites)
  nmds_scores[, Sample := info_df$Sample]
  nmds_scores[Sample %like% "N[0-9][0-9]", Host := "Nurses" ]
  nmds_scores[Sample %like% "F[0-9][0-9]", Host := "Foragers" ]
  nmds_scores[Sample %like% "N[0-9][0-9]", Host := "Nurses" ]
  nmds_scores[Sample %like% "F[0-9][0-9]", Host := "Foragers" ]
  nmds_scores[Sample %like% "[N F]0[1-3]", Location := "UNIL" ]
  nmds_scores[Sample %like% "[N F]0[4-7]", Location := "Liebefeld" ]
  nmds_scores[Sample %like% "[N F]0[4-7]", Location := "Liebefeld" ]
  nmds_scores[Sample %like% "[N F]0[8-9]", Location := "Yens" ]
  nmds_scores[Sample %like% "[N F]10", Location := "Yens" ]
  nmds_scores[Sample %like% "[N F]1[1-3]", Location := "Cugy" ]
  nmds_scores[Sample %like% "[N F]1[4-6]", Location := "Vesancy" ]
  nmds_plot <- ggplot(nmds_scores, aes(x=NMDS1, y=NMDS2, 
                                     shape=factor(Host,levels=treatments),
                                     color=factor(Location)))+
    geom_point(stat="identity",size=6, alpha = 0.75)+
    scale_color_manual(values=c("UNIL"="#d7191c",
                                  "Cugy"="#fdae61",
                                  "Yens"="#ffffbf",
                                  "Liebefeld"="#abd9e9",
                                  "Vesancy"="#2c7bb6"))+
    # scale_color_manual(values=c("UNIL"="#CC79A7",
    #                               "Cugy"="#0072B2",
    #                               "Yens"="#009E73",
    #                               "Liebefeld"="#E69F00",
    #                               "Vesancy"="#999999"))+
    scale_shape_manual(values=c("Foragers" = 17, "Nurses" = 19))+
    labs(x="NMDS1",
         y="NMDS2",
         title=plot_title,
         subtitle=paste0("Stress = ",nmds$stress),
         shape="Host",
         color="Location")+
    theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=22),
          plot.subtitle = element_text(hjust = 0.5)) +
          make_theme(setCol = F)
  return(nmds_plot)
}

#Using phylotype proportions
phylo_dt_umlt_prop <- dcast.data.table(phylo_all_dt, Sample ~ Phylo, value.var = "Prop", fun.aggregate = sum)
phylo_prop_mat <- as.matrix(phylo_dt_umlt_prop[, -1])
phylo_dt_umlt_prop[Sample %like% "N[0-9][0-9]", Host := "Nurses" ]
phylo_dt_umlt_prop[Sample %like% "F[0-9][0-9]", Host := "Foragers" ]
phylo_dt_umlt_prop[Sample %like% "[N F]0[1-3]", Location := "UNIL" ]
phylo_dt_umlt_prop[Sample %like% "[N F]0[4-7]", Location := "Liebefeld" ]
phylo_dt_umlt_prop[Sample %like% "[N F]0[4-7]", Location := "Liebefeld" ]
phylo_dt_umlt_prop[Sample %like% "[N F]0[8-9]", Location := "Yens" ]
phylo_dt_umlt_prop[Sample %like% "[N F]10", Location := "Yens" ]
phylo_dt_umlt_prop[Sample %like% "[N F]1[1-3]", Location := "Cugy" ]
phylo_dt_umlt_prop[Sample %like% "[N F]1[4-6]", Location := "Vesancy" ]
phylo_dt_umlt_prop[, Hive:=gsub("[FN]","", Sample)]
phylo_prop_nmdsplot_loc <- nmds_plot_loc(phylo_prop_mat, "Phylotype proportions-based NMDS plot", phylo_dt_umlt_prop)
phylo_prop_nmdsplot <- nmds_plot(phylo_prop_mat, "Phylotype proportions-based NMDS plot", phylo_dt_umlt_prop)
ggsave("figures/Fig1C-phylo_prop_nmdsplot.pdf", width = 10, height = 12)
# save figure phylo_prop_nmdsplot
adonis2_test_phylo <- adonis2(formula = phylo_prop_mat ~ Host + Location, data = phylo_dt_umlt_prop, permutations = 999, by = "margin")
adonis_test_phylo <- adonis_OmegaSq(adonis2_test_phylo)
df_adonis_phylo = as.data.frame(adonis_test_phylo)
write.csv(df_adonis_phylo, "figures/adonis_result_at_Phylo_level.csv", quote = F)

#Using SDP proportions
sdp_dt_umlt_prop <- dcast.data.table(phylo_all_dt, Sample ~ SDP, value.var = "Prop")
sdp_prop_mat <- as.matrix(sdp_dt_umlt_prop[,-1])
sdp_dt_umlt_prop[Sample %like% "N[0-9][0-9]", Host := "Nurses" ]
sdp_dt_umlt_prop[Sample %like% "F[0-9][0-9]", Host := "Foragers" ]
sdp_dt_umlt_prop[Sample %like% "[N F]0[1-3]", Location := "UNIL" ]
sdp_dt_umlt_prop[Sample %like% "[N F]0[4-7]", Location := "Liebefeld" ]
sdp_dt_umlt_prop[Sample %like% "[N F]0[4-7]", Location := "Liebefeld" ]
sdp_dt_umlt_prop[Sample %like% "[N F]0[8-9]", Location := "Yens" ]
sdp_dt_umlt_prop[Sample %like% "[N F]10", Location := "Yens" ]
sdp_dt_umlt_prop[Sample %like% "[N F]1[1-3]", Location := "Cugy" ]
sdp_dt_umlt_prop[Sample %like% "[N F]1[4-6]", Location := "Vesancy" ]
sdp_prop_nmdsplot_loc <- nmds_plot_loc(sdp_prop_mat, "SDP proportions-based NMDS plot", sdp_dt_umlt_prop)
sdp_prop_nmdsplot <- nmds_plot(sdp_prop_mat, "SDP proportions-based NMDS plot", sdp_dt_umlt_prop)
ggsave("figures/Fig1D-sdp_prop_nmdsplot.pdf", width = 10, height = 12)
adonis_OmegaSq(adonis2(formula = sdp_prop_mat ~ Location + Host, data = sdp_dt_umlt_prop, permutations = 999))
write.csv((as.data.frame(adonis_OmegaSq(adonis2(formula = sdp_prop_mat ~ Location + Host, data = phylo_dt_umlt_prop, permutations = 999, by = "margin")))), "figures/adonis_result_at_SDP_level.csv", quote = F)


df_plot_div_firm5 <- cbind(sdp_dt_umlt_prop[,c("firm5_1", "firm5_2", "firm5_3", "firm5_4", "firm5_7", "Host", "Location", "Sample")], "Shannon" = diversity(sdp_prop_mat[,c("firm5_1", "firm5_2", "firm5_3", "firm5_4", "firm5_7")], index = "shannon"))
df_plot_div_gilli <- cbind(sdp_dt_umlt_prop[,c("gilli_1", "gilli_2", "gilli_3", "gilli_4", "gilli_5", "gilli_6", "gilli_apis_andre", "gilli_apis_dorsa", "gilli_bombus", "Host", "Location", "Sample")], "Shannon" = diversity(sdp_prop_mat[,c("gilli_1", "gilli_2", "gilli_3", "gilli_4", "gilli_5", "gilli_6", "gilli_apis_andre", "gilli_apis_dorsa", "gilli_bombus")], index = "shannon"))
df_plot_div_bifido <- cbind(sdp_dt_umlt_prop[,c("bifido_1.1", "bifido_1.2", "bifido_1.3", "bifido_1.4", "bifido_1.5", "bifido_2", "Host", "Location", "Sample")], "Shannon" = diversity(sdp_prop_mat[,c("bifido_1.1", "bifido_1.2", "bifido_1.3", "bifido_1.4", "bifido_1.5", "bifido_2")], index = "shannon"))
df_plot_div_firm4 <- cbind(sdp_dt_umlt_prop[,c("firm4_1", "firm4_2", "Host", "Location", "Sample")], "Shannon" = diversity(sdp_prop_mat[,c("firm4_1", "firm4_2")], index = "shannon"))
# df_plot_div <- cbind(sdp_dt_umlt_prop, "Shannon" = diversity(sdp_prop_mat))
treat_colors_2 <- treat_colors[2:3]
ggplot(data = df_plot_div_firm5, 
               aes(x = Host, y = Shannon, fill = Host)) +
  geom_boxplot(outlier.fill = NA) +
  geom_point(data = df_plot_div_firm5, aes(x = Host, y = Shannon)) +
  scale_fill_manual(values = treat_colors_2) +
  ggtitle("Firm 5") +
  # compare_means(Shannon ~ Host, data = df_plot_div_firm5) +
  stat_compare_means(aes(label = "..p.signif.."),
                     size = 8,
                     method = "wilcox.test") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=24),
        axis.title.x=element_blank(),
        plot.title=element_text(hjust=0.5), 
        plot.subtitle=element_text(hjust=0.5))
        # stat_compare_means()
        ggsave("figures/SuppFig3-preShannon_diversity_firm5.pdf", width = 10, height = 12)

ggplot(data = df_plot_div_firm4, 
               aes(x = Host, y = Shannon, fill = Host)) +
  geom_boxplot(outlier.fill = NA) +
  geom_point(data = df_plot_div_firm4, aes(x = Host, y = Shannon)) +
  scale_fill_manual(values = treat_colors_2) +
  ggtitle("Firm 4") +
  # compare_means(Shannon ~ Host, data = df_plot_div_firm4) +
  stat_compare_means(aes(label = "..p.signif.."),
                     size = 8,
                     method = "wilcox.test") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=24),
        axis.title.x=element_blank(),
        plot.title=element_text(hjust=0.5), 
        plot.subtitle=element_text(hjust=0.5))
        # stat_compare_means()
        ggsave("figures/SuppFig3-preShannon_diversity_firm4.pdf", width = 10, height = 12)

ggplot(data = df_plot_div_gilli, 
               aes(x = Host, y = Shannon, fill = Host)) +
  geom_boxplot(outlier.fill = NA) +
  geom_point(data = df_plot_div_gilli, aes(x = Host, y = Shannon)) +
  scale_fill_manual(values = treat_colors_2) +
  ggtitle("Giliiamella") +
  # compare_means(Shannon ~ Host, data = df_plot_div_gilli) +
  stat_compare_means(aes(label = "..p.signif.."),
                     size = 8,
                     method = "wilcox.test") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=24),
        axis.title.x=element_blank(),
        plot.title=element_text(hjust=0.5), 
        plot.subtitle=element_text(hjust=0.5))
        # stat_compare_means()
        ggsave("figures/SuppFig3-preShannon_diversity_gilli.pdf", width = 10, height = 12)

ggplot(data = df_plot_div_bifido, 
               aes(x = Host, y = Shannon, fill = Host)) +
  geom_boxplot(outlier.fill = NA) +
  geom_point(data = df_plot_div_bifido, aes(x = Host, y = Shannon)) +
  scale_fill_manual(values = treat_colors_2) +
  ggtitle("Bifidobacterium") +
  # compare_means(Shannon ~ Host, data = df_plot_div_bifido) +
  stat_compare_means(aes(label = "..p.signif.."),
                     size = 8,
                     method = "wilcox.test") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=24),
        axis.title.x=element_blank(),
        plot.title=element_text(hjust=0.5), 
        plot.subtitle=element_text(hjust=0.5))
        # stat_compare_means()
        ggsave("figures/SuppFig3-pre-Shannon_diversity_bifido.pdf", width = 10, height = 12)

# Plot results of profiling unmapped reads
df_motus_raw <- read.csv("community_quantification/unmapped_reads_profiling/samples_merged.txt", sep = "\t", skip = 2, stringsAsFactors=FALSE)
df_motus_raw_samples <- df_motus_raw %>%
                          select(!consensus_taxonomy) %>%
                          select(!X.mOTU) %>%
                          select(!NCBI_tax_id) %>%
                          mutate_all(as.numeric) %>%
                          cbind("taxonomy" = df_motus_raw %>% select(X.mOTU))
                          # mutate(sum_ab = rowSums(across(where(is.numeric)))) %>%
                          #   filter(sum_ab > 0) %>%
                          #     select(!sum_ab) %>%
                              # filter(consensus_taxonomy != "unassigned")
df_unassigned <- df_motus_raw_samples %>%
                  filter(X.mOTU == "unassigned")
tax_info <- df_motus_raw %>%
              select(consensus_taxonomy, X.mOTU) %>%
              unique() %>%
                separate(into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), col = consensus_taxonomy, sep = "\\|")
df_motus_raw_samples_t <- as.data.frame(t(df_motus_raw_samples)) %>%
              rownames_to_column("Sample")
colnames(df_motus_raw_samples_t) <- c("Sample", df_motus_raw_samples$X.mOTU)
df_motus <- df_motus_raw_samples_t %>%
              pivot_longer(!Sample, names_to = "motu", values_to = "rel_ab") %>%
                group_by(Sample, motu) %>%
                  mutate(Present = ifelse(rel_ab > 0, 1, 0)) %>%
                    group_by(motu) %>%
                     mutate(Prevalence_num = sum(Present)) %>%
                      mutate(Prevalence = mean(Present)) %>%
                        left_join(tax_info, by = c("motu" = "X.mOTU"))
readnum_dt <- fread("readnumbers.csv")
setnames(readnum_dt, c("Sample","Raw","Bacterial_db","A_mellifera"))
readnum_dt[, Leftover := Raw-Bacterial_db-A_mellifera]
readnum_dt[, prop_bac := Bacterial_db/Raw]
readnum_dt[, prop_amel := A_mellifera/Raw]
reads <- readnum_dt %>%
          mutate(unmapped_perc = round(Leftover/Raw*100)) %>%
          mutate(unmapped_perc = paste0(unmapped_perc, " %")) %>%
          select(Sample, unmapped_perc)
meta_dt <- fread("SamplingGilles2019.csv", header = TRUE) %>%
            left_join(reads)
df_motus_combined <- df_motus %>%
                      left_join(meta_dt) %>%
                      mutate(rel_ab = as.numeric(rel_ab))
Nlist <- c("N01","N02","N03","N04","N06","N07","N08","N09","N10","N11","N12","N13","N14","N15","N16")
Flist <- c("F01","F02","F03","F04","F06","F07","F08","F09","F10","F11","F12","F13","F14","F15","F16")
samplelist <- c(Nlist, Flist)
sample_order <- c(Nlist, Flist)
df_motus_combined %>% filter(rel_ab > 0.01) %>%
  filter(is.na(Genus))
split_genus <- function(g_name){
  final_name <- g_name
  if (!is.na(g_name)) {
    # if (grepl(" ", g_name)) {
    #   # split by space and keep all but the text after the last space
    #   final_name <- strsplit(g_name, " ")[[1]][-length(strsplit(g_name, " ")[[1]])]
    # }
    final_name <- strsplit(g_name, "g__")[[1]][2]
  }
  return(final_name)
}
genera_added <- df_motus_combined %>% ungroup() %>%
                filter(rel_ab > 0.01) %>% 
                mutate(Genus = Vectorize(split_genus)(Genus)) %>% 
                arrange(Sample, rel_ab) %>% pull(Genus) %>% unlist %>% unique()
ggplot(df_motus_combined %>% 
        filter(rel_ab > 0.01) %>%
        mutate(Genus = Vectorize(split_genus)(Genus)),
       aes(x = factor(Sample, samplelist), 
           y = rel_ab,
          fill = factor(Genus, genera_added))) +
                geom_bar(position = "stack", stat = "identity", color = "black", size = 0.2) +
                geom_text(aes(label = unmapped_perc,
                              y = -0.03,
                              group = Sample,
                          ), size = 3, angle = 45) +
                labs(fill = "Genus", y = "mOTUs2 Relative abundance", x = "Sample") +
                facet_wrap(~Sample_type, scales="free") +
                  make_theme(setFill = T, max_colors = length(genera_added),
                             palettefill = "Set1",
                             leg_pos = "bottom",
                             guide_nrow = 10,
                             leg_size = 12,
                             x_angle=45 ,x_vj=1.2, x_hj=1, x_size=12,
                             y_angle=0 ,y_vj=0, y_hj=0, y_size=12
                  ) +
                  theme(axis.title.x = element_text(size = 18),
                      axis.title.y = element_text(size = 18),
                      legend.title = element_text(size = 18),
                      strip.text.x = element_text(size = 18),
                )
  ggsave("figures/SuppFig2-pre-motus2_unmapped.pdf", width = 12, height = 12)
df_to_write <- df_motus_combined %>% 
        filter(rel_ab > 0.01) %>%
        mutate(Genus = Vectorize(split_genus)(Genus))
write.csv(df_to_write, "community_quantification/unmapped_reads_profiling/motus2_profile.csv", quote = F, row.names = FALSE)