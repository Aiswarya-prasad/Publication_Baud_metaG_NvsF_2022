library(tidyverse)
library(data.table)
library(ggplot2)
library(ape)
library(gridExtra)
library(RColorBrewer)
library(ggnewscale)
library(ggrepel)
library(vegan)
library(qvalue)
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

coverage_in_sample <- function(sdp, sample, limit=20){
  # function to determine the ter coverage of a given species (sdp) in a given sample
  df_covs <- read.csv(paste0("snvs/", sdp, "_corecov_coord.txt"), sep = '\t')
  coverage <- df_covs[df_covs$Sample==sample,]$Cov_ter
  return(coverage)
}

source("utilities.R")
figpath <- "figures/"
datapath <- "snvs/"

Nlist <- c("N01","N02","N03","N04","N06","N07","N08","N09","N10","N11","N12","N13","N14","N15","N16")
Flist <- c("F01","F02","F03","F04","F06","F07","F08","F09","F10","F11","F12","F13","F14","F15","F16")
samplelist <- c(Nlist, Flist)
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
snod_1_Nlike_strains <- c("Ga0326500","CPT77","BGH94","BGH95","BGH96","BGH97","BGI00", 
                          "BGI01","BGI08","BGI09","BGI12","BGI13","BHC42","BHC50", 
                          "BHC52","Ga0227306","SALWKB2","H3V01","H3V10","H3V11",
                          "H3U73","H3U82","H3T75","Ga0227304")
snod_1_Flike_strains <- c("BHC47","Ga0227305","H3T79")
sample_order <- c(Nlist, Flist)
phylo_order <- c("firm5","bifido","firm4","gilli","snod","api","fper","bapis","com","lkun", "bom")
phylo_order_rv <- rev(phylo_order)
phylo_list <- c("api","bapis","bifido","bom","com","firm4","firm5","fper","gilli","lkun","snod")
sdp_list <- c("api_1","bapis","bifido_1.1","bifido_1.2",
              "bifido_1.3","bifido_1.4","bifido_1.5","bifido_2","bom_1","com_1",
              "firm4_1","firm4_2","firm5_1","firm5_2","firm5_3","firm5_4","firm5_7",
              "fper_1","gilli_1","gilli_2","gilli_3","gilli_4","gilli_5","gilli_6",
              "lkun","snod_1","snod_2")
sdp_list_restricted <- c("bapis","bifido_1.1","bifido_1.2","bifido_1.3","bifido_1.4",
                         "bifido_1.5","bifido_2","firm4_1","firm4_2","firm5_1",
                         "firm5_2","firm5_3","firm5_4","fper_1","gilli_1",
                         "gilli_2","snod_1")
treat_colors <- c("#D55E00","#56B4E9","#F0E442")
colnames(sdp_dt_umlt_prop)  
samp_dt <- fread("SamplingGilles2019.csv")
samp_dt <- samp_dt[, c("Sample", "Sample_type", "Location", "Number_guts", "Average_gut_mass")]

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

# plots are made for various cutoffs the one where no coverage cutoff is specified is KE's approach
# which is more or less equivalent to the 0.1 cutoff
# the abbreviation _01 corresponds to the 0.01 cutoff, _1 to the 0.1 cutoff and _05 to the 0.05 cutoff
# all those cutoffs also consider the allele only if it is supported by at least 2 reads
# the _0 cutoff is the case where no cutoffs (frequency or coverage) are applied

# Plot the fraction of polymorphic sites per sample for each species at

plot_perc_poly_sites <- function(in_file, out_name, AP = FALSE) {
  frac_snvs_dt <- fread(in_file,h=F)
  setnames(frac_snvs_dt, c("Phylo", "SDP","Host","Colony","Sample", "Nb_snps","Fraction_var"))
  aggr <- as.data.table(aggregate(Fraction_var ~ SDP + Host + Phylo, frac_snvs_dt, median))

  if (AP) {
    frac_snvs_dt[, Fraction_var:=as.numeric(Fraction_var)]
    frac_snvs_dt <- frac_snvs_dt[!(Fraction_var)==0,]
    aggr <- as.data.table(aggregate(Fraction_var ~ SDP + Host + Phylo, frac_snvs_dt, median))
  }

  stats_frac_var <- data.table(sdp_list_restricted, 0)
  setnames(stats_frac_var, c("SDP", "p"))
  for (i in sdp_list_restricted) {
    temp_test <- wilcox.test(frac_snvs_dt[(SDP==i & Sample %like% "N[0-9][0-9]"), Fraction_var],
                            frac_snvs_dt[(SDP==i & Sample %like% "F[0-9][0-9]"), Fraction_var],
                            paired = FALSE)
    stats_frac_var[ SDP==i, p:=temp_test$p.value]
  }
  stats_frac_var[, q:=qvalue(stats_frac_var$p)$qvalues]
  stats_frac_var <- merge.data.table(stats_frac_var,
                                        frac_snvs_dt[, .N, by=c("SDP","Phylo")][, c("SDP", "Phylo")], 
                                        by="SDP", all.x=TRUE )

  frac_snvs_dt <- frac_snvs_dt[!(Host=="g"),]
  write.csv(stats_frac_var, paste0("figures/", out_name ,"_stats_frac_var.csv"), quote=F)
  ggplot(data=frac_snvs_dt[SDP %in% sdp_list_restricted,], 
        aes(x=SDP,
            y=Fraction_var)) +
    geom_jitter(position=position_jitterdodge(), alpha=0.5, size=4,
                aes(color=factor(Host,levels=c("N","F"),labels=treatments))) +
    geom_point(data = aggr[!(Fraction_var)==0 & SDP %in% sdp_list_restricted,], 
              aes(y=Fraction_var, 
                  x=SDP,
                  fill=factor(Host,levels=c("N","F"),labels=treatments),
                  color=factor(Host,levels=c("N","F"),labels=treatments)),
              position=position_jitterdodge(),
              colour="black",
              show.legend = FALSE,
              shape=3,
              alpha=1)+
    geom_point(data=stats_frac_var[q<0.05, ],
              aes(x=SDP, y=8.5), shape="*", size=8, show.legend = FALSE)+
    theme(axis.text.x=element_text(angle=90),
          text=element_text(size = 16), 
          panel.background = element_rect(fill = "#f2f2f2", color = "black"),
          panel.grid.major=element_line(color="white"),
          panel.grid.minor=element_line(color="white")) +
    xlab("Species")+
    scale_y_continuous(name="Polymorphic sites per sample (%)", limits=c(0,11.5), n.breaks=10)+
    facet_grid(. ~ Phylo, scales="free", space="free")+
    scale_color_manual(values=treat_colors, name="Host") #+
    ggsave(paste0(figpath, "Fig2A_", out_name, ".pdf"), width=15, height=10)

    stats_corr_ct_snvs <- data.table(sdp_list_restricted, 0,0,0,0)
  setnames(stats_corr_ct_snvs, c("SDP", "rN","rF","pN","pF"))
  tmp_dt <- copy(frac_snvs_dt)
  tmp_m_dt <- merge.data.table(tmp_dt, CT_dt[Target %like% "UV", c("Sample", "copiesperngDNA")], by="Sample")
  tmp_phylo_all_dt <- copy(phylo_all_dt)
  tmp_phylo_all_dt[SDP=="fper", SDP:="fper_1"]
  tmp_m_dt <- merge.data.table(tmp_m_dt, tmp_phylo_all_dt[, c("Sample","SDP", "Cov_ter", "Prop")], by=c("Sample","SDP"), all.x=TRUE)
  tmp_m_dt[, ratio := Cov_ter/(copiesperngDNA*Prop)]
  tmp_m_dt[Host=="N", Host:="Nurses"]
  tmp_m_dt[Host=="F", Host:="Foragers"]
  for(i in sdp_list_restricted){
    rn <- format(round(cor(tmp_m_dt[SDP==i & Host=="Nurses", ratio],
                          tmp_m_dt[SDP==i & Host=="Nurses", Fraction_var], method = "pearson"),3),nsmall=3)
    stats_corr_ct_snvs[SDP==i, rN:=rn]
    rf <- format(round(cor(tmp_m_dt[SDP==i & Host=="Foragers", ratio],
                          tmp_m_dt[SDP==i & Host=="Foragers", Fraction_var], method = "pearson"),3),nsmall=3)
    stats_corr_ct_snvs[SDP==i, rF:=rf]
    pn <- try(format(round(cor.test(tmp_m_dt[SDP==i & Host=="Nurses", ratio],
                                    tmp_m_dt[SDP==i & Host=="Nurses", Fraction_var], method = "pearson")$p.value,3),
                    nsmall=3),silent=TRUE)
    stats_corr_ct_snvs[SDP==i, pN:=pn]
    pf <- try(format(round(cor.test(tmp_m_dt[SDP==i & Host=="Foragers", ratio],
                                    tmp_m_dt[SDP==i & Host=="Foragers", Fraction_var], method = "pearson")$p.value,3),
                    nsmall=3),silent=TRUE)
    stats_corr_ct_snvs[SDP==i, pF:=pf]
  }
  stats_corr_ct_snvs <- merge.data.table(stats_corr_ct_snvs, tmp_m_dt[, c("SDP","Phylo")], by="SDP", all.x=TRUE)

  ggplot(tmp_m_dt[SDP %in% sdp_list_restricted,], aes(x=as.numeric(ratio), y=Fraction_var))+
    geom_point(aes(fill=factor(Host,levels=treatments),
                  shape=factor(Host,levels=treatments)), 
              size=4)+
    scale_shape_manual(name="Host",values=c(21,24))+
    scale_fill_manual(name="Host", values=treat_colors)+
    geom_text(data=stats_corr_ct_snvs, 
              aes(x = Inf, y = -Inf, label=paste0("r(Nurses)=",rN,"\n",
                                                  "p(Nurses)=",pN,"\n",
                                                  "r(Foragers)=",rF,"\n",
                                                  "p(Foragers)=",pF)), 
              show.legend = FALSE,
              vjust=-0.5,
              hjust=1.1,
              size=4)+
    xlab("Terminus coverage (reads/bp) \n per bacterial load (SDP proportion * 16S copy/ng of DNA)")+
    ylab("Polymorphic sites per sample (%)")+
    facet_wrap(SDP ~ ., scales="free_y", ncol = 2)+
    theme(axis.text.x=element_text(angle = 45, hjust = 1),
          text=element_text(size = 16), 
          panel.background = element_rect(fill = "#f2f2f2", color = "black"),
          panel.grid.major=element_line(color="white"),
          panel.grid.minor=element_line(color="white"))
    ggsave(paste0(figpath, "SuppFig5_", out_name, ".pdf"), width=25, height=20)
}

plot_perc_poly_sites(paste0(datapath, "all_sample_var_host.txt"), "KE_approach")
# plot_perc_poly_sites(paste0(datapath, "all_sample_var_host_filt_0.txt"), "filt_0", AP = TRUE)
plot_perc_poly_sites(paste0(datapath, "all_sample_var_host_filt_01.txt"), "filt_01", AP = TRUE)
plot_perc_poly_sites(paste0(datapath, "_all_sample_var_host_filt_05.txt"), "filt_05", AP = TRUE)
plot_perc_poly_sites(paste0(datapath, "all_sample_var_host_filt_1.txt"), "filt_1", AP = TRUE)


# Plotting cumulative curves for %polymorphic sites per sample detected

plot_cum_curves <- function(in_file_pre, in_file_suf, out_name, AP = FALSE) {
    if (AP) {
        cum_curves_dt <- data.table()
        for (i in sdp_list_restricted) {
        file_path <- paste0(in_file_pre, i, in_file_suf)
        if (file.exists(file_path)) {
            tmp <-  fread(file_path, h=T)
            cum_curves_dt <- rbind(cum_curves_dt, tmp)
        }
        }
        # as numeric nsamples
        setnames(cum_curves_dt, c("Host","curve","Nsamples","percvar","SDP"))
        cum_curves_dt <- cum_curves_dt[ , Nsamples:=as.numeric(Nsamples)]
        cum_curves_dt <- cum_curves_dt[ , percvar:=as.numeric(percvar)]
        ggplot(data=cum_curves_dt[(Host=="F" | Host=="N"),], 
            aes(x=Nsamples, 
                y=percvar,
                colour=factor(Host,
                                levels=c("N", "F"),
                                labels=treatments))) +
        geom_jitter(position=position_dodge(width=0.7)) + 
        geom_smooth(se=FALSE)+
        xlab("Number of bee samples") + 
        ylab("Polymorphic sites (%)") +
        scale_x_continuous(breaks=c(0,15))+
        scale_color_manual(values = treat_colors, name="Host")+
        facet_wrap(. ~ SDP, ncol=6,nrow=3)+
        theme(text=element_text(size = 16), 
                panel.background = element_rect(fill = "#f2f2f2", color = "black"),
                legend.title=element_blank(),
                panel.grid.major=element_line(color="white"),
                panel.grid.minor=element_line(color="white"))
        ggsave(paste0(figpath, "Fig2B_", out_name,".pdf"), width=15, height=10)

    } else {
        cum_curves_dt <- data.table()
        for (i in sdp_list_restricted) {
        tmp <-  fread(paste0(in_file_pre,i, in_file_suf), h=F)
        cum_curves_dt <- rbind(cum_curves_dt, tmp)
        }
        setnames(cum_curves_dt, c("Host","curve","Nsamples","percvar","SDP"))
        ggplot(data=cum_curves_dt[(Host=="F" | Host=="N"),], 
            aes(x=Nsamples, 
                y=percvar,
                colour=factor(Host,
                                levels=c("N", "F"),
                                labels=treatments))) +
        geom_jitter(position=position_dodge(width=0.7)) + 
        geom_smooth(se=FALSE)+
        xlab("Number of bee samples") + 
        ylab("Polymorphic sites (%)") +
        scale_x_continuous(breaks=c(0,15))+
        scale_color_manual(values = treat_colors, name="Host")+
        facet_wrap(. ~ SDP, ncol=6,nrow=3)+
        theme(text=element_text(size = 16), 
                panel.background = element_rect(fill = "#f2f2f2", color = "black"),
                legend.title=element_blank(),
                panel.grid.major=element_line(color="white"),
                panel.grid.minor=element_line(color="white"))
        ggsave(paste0(figpath, "Fig2B_", out_name, ".pdf"), width=15, height=10)
    }
}

plot_cum_curves("snvs/", "_cum_curve.txt", "KE_approach", FALSE)
plot_cum_curves("snvs/", "_cum_curve_filt_0.txt", "filt_0", AP = TRUE)
plot_cum_curves("snvs/", "_cum_curve_filt_0.01.txt", "filt_01", AP = TRUE)
plot_cum_curves("snvs/", "_cum_curve_filt_0.05.txt", "filt_05", AP = TRUE)
plot_cum_curves("snvs/", "_cum_curve_filt_0.1.txt", "filt_1", AP = TRUE)

pcoa_plot<- function(dt_c, graph_title, adonis_info, locations_enabled=TRUE) {
  if (nrow(dt_c[V1 %like% "F[0-9][0-9]",])==0 | nrow(dt_c[V1 %like% "N[0-9][0-9]",]) == 0) {
    return()
  }
  dt_matrix <- as.matrix(dt_c, rownames = "V1")
  mat_dist <- as.dist(dt_matrix)
  res_pcoa <- pcoa(mat_dist)
  ev1 <- res_pcoa$vectors[,1]
  ev2 <- res_pcoa$vectors[,2]
  df_new <- data.frame(cbind(ev1,ev2))
  host_list <- substr(rownames(df_new),1,1)
  df_new <- cbind(df_new,host_list)
  locations <- rownames(df_new)
  hive <- gsub("[N F]", "", locations)
  df_new$hive <- hive
  locations <- gsub("[N F]0[1-3]", "UNIL", locations)
  locations <- gsub("[N F]0[4-7]", "Liebefeld", locations)
  locations <- gsub("[N F]0[8-9]", "Yens", locations)
  locations <- gsub("[N F]10", "Yens", locations)
  locations <- gsub("[N F]1[1-3]", "Cugy", locations)
  locations <- gsub("[N F]1[4-6]", "Vesancy", locations)
  df_new$locations <- locations
  perc_axis <- round(((res_pcoa$values$Relative_eig[c(1,2)])*100), digits=1)
  axis_x_title <- paste0("PCo1 (",perc_axis[1],"%)")
  axis_y_title <- paste0("PCo2 (",perc_axis[2],"%)")
  
  if (!locations_enabled) {
    p <- ggplot(df_new,aes(x=ev1,
                           y=ev2,
                           color=factor(host_list,
                                        levels=c("N","F"),
                                        labels=c("Nurses","Foragers")), 
                           shape=factor(host_list,
                                        levels=c("N","F"),
                                        labels=c("Nurses","Foragers"))))+
      geom_point(stat="identity",size=4)+
      labs(x=as.character(axis_x_title),
           y=as.character(axis_y_title), 
           title=graph_title,
           subtitle = paste0("Permanova\n Host: effect size = ", round(adonis_info[SDP==graph_title,`R2-Host`],4),", p = ", round(adonis_info[SDP==graph_title,`pvalue-Host`],4), 
           "\n Location: effect size = ", round(adonis_info[SDP==graph_title,`R2-Location`],4),", p = ", round(adonis_info[SDP==graph_title,`pvalue-Location`],4)))+
          #  Model: dist ~ Host + Location + Hive\n Permutations = 999
      scale_color_manual(values = treat_colors, name = "Host")+
      scale_shape_discrete(name="Host")+
      theme(legend.title =element_blank(),
            legend.position="none",
            panel.background=element_rect(fill="#f2f2f2"),
            panel.grid.major=element_line(linewidth=0.5,linetype='solid',colour="white"), 
            panel.grid.minor=element_line(linewidth=0.25,linetype='solid',colour="white"),
            axis.title.x=element_text(size=12),
            axis.title.y=element_text(size=12),
            plot.title=element_text(size=15,face="bold"),
            plot.subtitle=element_text(size=8,face="plain"),
            axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8))
  } else {
    p <- ggplot(df_new,aes(x=ev1,
                           y=ev2,
                           colour=factor(locations), 
                           alpha=factor(host_list,
                                        levels=c("N","F"),
                                        labels=c("Nurses","Foragers")),
                           shape=factor(host_list,
                                        levels=c("N","F"),
                                        labels=c("Nurses","Foragers")),
                           size=factor(host_list,
                                        levels=c("N","F"),
                                        labels=c("Nurses","Foragers"))))+
      geom_point(stat="identity", size=4)+
      geom_line(aes(group=factor(hive)), alpha=0.5, size=0.5)+
      scale_size_manual(values=c("Foragers"=2,"Nurses"=3), guide="none")+
      scale_alpha_manual(values = c("Foragers"=1, "Nurses"=0.5), guide='none')+
      scale_color_manual(values=c("UNIL"="#CC79A7",
                                  "Cugy"="#0072B2",
                                  "Yens"="#009E73",
                                  "Liebefeld"="#E69F00",
                                  "Vesancy"="#999999"),
                         name = "Host")+
      labs(x=as.character(axis_x_title),
           y=as.character(axis_y_title), 
           title=graph_title,
           subtitle = paste0("Permanova\n Host: R2 = ", round(adonis_info[SDP==graph_title,`R2-Host`],4),", p = ", round(adonis_info[SDP==graph_title,`pvalue-Host`],4), 
           "\n Location: R2 = ", round(adonis_info[SDP==graph_title,`R2-Location`],4),", p = ", round(adonis_info[SDP==graph_title,`pvalue-Location`],4)))+
      #  Model: dist ~ Host + Location + Hive\n Permutations = 999 - add this to figure legend
      theme(legend.title =element_blank(),
            legend.position="none",
            panel.background=element_rect(fill="#f2f2f2"),
            panel.grid.major=element_line(linewidth=0.5,linetype='solid',colour="white"), 
            panel.grid.minor=element_line(linewidth=0.25,linetype='solid',colour="white"),
            axis.title.x=element_text(size=12),
            axis.title.y=element_text(size=12),
            plot.title=element_text(size=15,face="bold"),
            plot.subtitle=element_text(size=8,face="plain"),
            axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8))
  }
  
  return(p)
}
# writing a new function to use adonis for 
# permonova as per a model rather than one-way with PERMANOVA

get_permanova<- function(dt_c) {
  if (nrow(dt_c[V1 %like% "F[0-9][0-9]",])==0 | nrow(dt_c[V1 %like% "N[0-9][0-9]",]) == 0) {
    return()
  }
  dt_matrix <- as.matrix(dt_c, rownames = "V1")
  mat_dist <- as.dist(dt_matrix)
  locations <- dt_c$V1
  hive <- gsub("[N F]", "", locations)
  locations <- gsub("[N F]0[1-3]", "UNIL", locations)
  locations <- gsub("[N F]0[4-7]", "Liebefeld", locations)
  locations <- gsub("[N F]0[8-9]", "Yens", locations)
  locations <- gsub("[N F]10", "Yens", locations)
  locations <- gsub("[N F]1[1-3]", "Cugy", locations)
  locations <- gsub("[N F]1[4-6]", "Vesancy", locations)
  host_type <- gsub(dt_c$V1, pattern = "[0-9][0-9]", replacement = "")
  data_info <- data.frame(cbind("Sample" = dt_c$V1, "Location" = locations, "Hive" = hive, "Host" = host_type))
  perm <- adonis2(mat_dist ~ Host + Location, data_info, permutations = 999, by = "margin")
  return(perm)
}

snvs_sdp_dt <- list()
snvs_sdp_plots <- list()
permanovastats <- list()
permanova_vals <- data.table()

#Making the plots

# distance matrix was made from snvs/{sdp}_shared_fraction...txt using the script distance_matrix.R

for (i in sdp_list_restricted) {
  snvs_sdp_dt[[i]] <- fread(paste0(datapath, "/distance_matrices/", i, "_dist_matrix_KE.txt"),h=T)
  snvs_sdp_dt[[i]] <- snvs_sdp_dt[[i]][, c("V1", snvs_sdp_dt[[i]][, V1]), with=FALSE]
  dt_send <- copy(snvs_sdp_dt[[i]])
  db_cols <- dt_send[!(V1 %like% "F[0-9][0-9]" | V1 %like% "N[0-9][0-9]"), V1]
  dt_send <- dt_send[V1 %like% "F[0-9][0-9]" | V1 %like% "N[0-9][0-9]", ]
  dt_send[, (db_cols):=NULL]
  permanovastats[[i]] <- get_permanova(dt_send)
  permanova_info <- data.table(i,
                                                     adonis_OmegaSq(permanovastats[[i]])$`Pr(>F)`[[1]], # host pval
                                                     adonis_OmegaSq(permanovastats[[i]])$`Pr(>F)`[[2]], # location pval
                                                     adonis_OmegaSq(permanovastats[[i]])$`parOmegaSq`[[1]], # host R2
                                                     adonis_OmegaSq(permanovastats[[i]])$`parOmegaSq`[[2]] # location R2
                                          )
  permanova_vals <- rbind(permanova_vals, permanova_info
                    )
  setnames(permanova_info, c("SDP","pvalue-Host", "pvalue-Location","R2-Host", "R2-Location"))
  snvs_sdp_plots[[i]] <- pcoa_plot(dt_send, i, permanova_info, locations_enabled=FALSE)
}

pdf(paste0(figpath, "Fig2C_", "KE_approach", ".pdf"), width=25,height=15)
marrangeGrob(snvs_sdp_plots, nrow = 3, ncol = 6, layout_matrix = matrix(1:20,3,6,TRUE))
dev.off()

setnames(permanova_vals, c("SDP","pvalue-Host", "pvalue-Location","R2-Host", "R2-Location"))
permanova_vals[, Type:="Gilles_KE"]
write.table(permanova_vals, file = paste0("figures/Fig2C_permanova_vals_KE_approach.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)


snvs_sdp_dt <- list()
plots_w_loc <- list()
#Making the plots
for (i in sdp_list_restricted) {
  snvs_sdp_dt[[i]] <- fread(paste0(datapath, "/distance_matrices/", i, "_dist_matrix_KE.txt"),h=T)
  snvs_sdp_dt[[i]] <- snvs_sdp_dt[[i]][, c("V1", snvs_sdp_dt[[i]][, V1]), with=FALSE]
  dt_send <- copy(snvs_sdp_dt[[i]])
  db_cols <- dt_send[!(V1 %like% "F[0-9][0-9]" | V1 %like% "N[0-9][0-9]"), V1]
  dt_send <- dt_send[V1 %like% "F[0-9][0-9]" | V1 %like% "N[0-9][0-9]", ]
  dt_send[, (db_cols):=NULL]
  permanovastats[[i]] <- get_permanova(dt_send)
  permanova_info <- data.table(i,
                                                     adonis_OmegaSq(permanovastats[[i]])$`Pr(>F)`[[1]], # host pval
                                                     adonis_OmegaSq(permanovastats[[i]])$`Pr(>F)`[[2]], # location pval
                                                     adonis_OmegaSq(permanovastats[[i]])$`parOmegaSq`[[1]], # host R2
                                                     adonis_OmegaSq(permanovastats[[i]])$`parOmegaSq`[[2]] # location R2
                                          )
  setnames(permanova_info, c("SDP","pvalue-Host", "pvalue-Location","R2-Host", "R2-Location"))
  plots_w_loc[[i]] <- pcoa_plot(dt_send, i, permanova_info, locations_enabled=T)
}

pdf(paste0(figpath, "SuppFig6_", "KE_approach", ".pdf"), width=25,height=15)
marrangeGrob(plots_w_loc, nrow = 3, ncol = 6, layout_matrix = matrix(1:20,3,6,TRUE))
dev.off()

plot_and_make_dist_mat <- function(outname) {
    for (sdp in sdp_list_restricted) {
    outfile_name <- paste0("snvs/distance_matrices/", sdp, "_dist_matrix_", outname,".txt")
    data <- read.table(paste0("snvs/", sdp, "_shared_fraction_", outname,".txt"),h=T)#, check.names=FALSE)
    # calculate jaccard dist as 1- shared fraction
    jaccard_dist <- 1 - data$Shared_fraction
    data_new <- cbind(data, jaccard_dist)
    #Create a distance matrix from the data, and print it out

    samples <- unique(data_new$Sample1)
    # keep only samples for which coverage in sample is > 20 using the funciton defined above
    data.frame(samples) %>%
        mutate(coverage = Vectorize(coverage_in_sample)(sdp, samples)) %>%
        filter(coverage > 20) %>%
        pull(samples) -> samples_to_keep
    samples <- samples_to_keep
    nb_samples <- length(samples)
    sample_index <- c(1:nb_samples)
    matrix_zeros <- rep(0, nb_samples*nb_samples)
    dist_matrix <- matrix(matrix_zeros, nrow=nb_samples, ncol=nb_samples)
    rownames(dist_matrix) <- samples
    colnames(dist_matrix) <- samples

    for (i in sample_index) {
        for (j in sample_index) {
            if (samples[i] == samples[j]) {
            dist_matrix[i,j] = 0
            }
        else {
            data_line <- subset(data_new, data_new$Sample1 == samples[i] & data_new$Sample2 == samples[j])
            if (dim(data_line)[[1]] > 0) {
                dist_matrix[i,j] = data_line$jaccard_dist

                } else {
                dist_matrix[i,j] = dist_matrix[j,i]
                }
            }
        }
    }
    write.table(dist_matrix,file=outfile_name,row.names=TRUE,col.names=NA,quote=FALSE,sep="\t")
    }


    #Reading the data

    snvs_sdp_dt <- list()
    snvs_sdp_plots <- list()
    permanovastats <- list()
    permanova_vals <- data.table()

    #Making the plots
    for (i in sdp_list_restricted) {
    snvs_sdp_dt[[i]] <- fread(paste0("snvs/distance_matrices/", i, "_dist_matrix_", outname,".txt"),h=T)
    snvs_sdp_dt[[i]] <- snvs_sdp_dt[[i]][, c("V1", snvs_sdp_dt[[i]][, V1]), with=FALSE]
    dt_send <- copy(snvs_sdp_dt[[i]])
    db_cols <- dt_send[!(V1 %like% "F[0-9][0-9]" | V1 %like% "N[0-9][0-9]"), V1]
    dt_send <- dt_send[V1 %like% "F[0-9][0-9]" | V1 %like% "N[0-9][0-9]", ]
    dt_send[, (db_cols):=NULL]
    permanovastats[[i]] <- get_permanova(dt_send)
    permanova_info <- data.table(i,
                                                        adonis_OmegaSq(permanovastats[[i]])$`Pr(>F)`[[1]], # host pval
                                                        adonis_OmegaSq(permanovastats[[i]])$`Pr(>F)`[[2]], # location pval
                                                        adonis_OmegaSq(permanovastats[[i]])$`parOmegaSq`[[1]], # host R2
                                                        adonis_OmegaSq(permanovastats[[i]])$`parOmegaSq`[[2]] # location R2
                                            )
    permanova_vals <- rbind(permanova_vals, permanova_info
                        )
    setnames(permanova_info, c("SDP","pvalue-Host", "pvalue-Location","R2-Host", "R2-Location"))
    snvs_sdp_plots[[i]] <- pcoa_plot(dt_send, i, permanova_info, locations_enabled=FALSE)
    }

    pdf(paste0(figpath, "Fig2C_", outname,"_backgroundLight.pdf"), width=25,height=15)
    marrangeGrob(snvs_sdp_plots, nrow = 3, ncol = 6, layout_matrix = matrix(1:20,3,6,TRUE))
    dev.off()

    setnames(permanova_vals, c("SDP","pvalue-Host", "pvalue-Location","R2-Host", "R2-Location"))
    permanova_vals[, Type:="Gilles_KE"]
    write.table(permanova_vals, file = paste0("snvs/Fig2C", outname,"_permanova_vals.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)


    snvs_sdp_dt <- list()
    plots_w_loc <- list()
    #Making the plots
    for (i in sdp_list_restricted) {
    snvs_sdp_dt[[i]] <- fread(paste0("snvs/distance_matrices/", i, "_dist_matrix_", outname,".txt"),h=T)
    snvs_sdp_dt[[i]] <- snvs_sdp_dt[[i]][, c("V1", snvs_sdp_dt[[i]][, V1]), with=FALSE]
    dt_send <- copy(snvs_sdp_dt[[i]])
    db_cols <- dt_send[!(V1 %like% "F[0-9][0-9]" | V1 %like% "N[0-9][0-9]"), V1]
    dt_send <- dt_send[V1 %like% "F[0-9][0-9]" | V1 %like% "N[0-9][0-9]", ]
    dt_send[, (db_cols):=NULL]
    permanovastats[[i]] <- get_permanova(dt_send)
    permanova_info <- data.table(i,
                                                        adonis_OmegaSq(permanovastats[[i]])$`Pr(>F)`[[1]], # host pval
                                                        adonis_OmegaSq(permanovastats[[i]])$`Pr(>F)`[[2]], # location pval
                                                        adonis_OmegaSq(permanovastats[[i]])$`parOmegaSq`[[1]], # host R2
                                                        adonis_OmegaSq(permanovastats[[i]])$`parOmegaSq`[[2]] # location R2
                                            )
    setnames(permanova_info, c("SDP","pvalue-Host", "pvalue-Location","R2-Host", "R2-Location"))
    plots_w_loc[[i]] <- pcoa_plot(dt_send, i, permanova_info, locations_enabled=T)
    }

    pdf(paste0(figpath, "SuppFig6", outname,"_backgroundLight.pdf"), width=25,height=15)
    marrangeGrob(plots_w_loc, nrow = 3, ncol = 6, layout_matrix = matrix(1:20,3,6,TRUE))
    dev.off()
}

plot_and_make_dist_mat("filt_0")
plot_and_make_dist_mat("filt_01")
plot_and_make_dist_mat("filt_05")
plot_and_make_dist_mat("filt_1")
