#!/usr/bin/env Rscript

# --------------------------------------------------------------------------------
# --- Plot basic pop gen statistics
# --------------------------------------------------------------------------------

options(stringsAsFactors = FALSE)

library(ggplot2)

species = c('AG', 'AR')
islands = c('BANDA', 'BUGALA', 'BUKASA', 'BUWAMA', 'ENTEBBE', 'KAZZI', 'KIYINDI', 
            'KOOME', 'MITYANA', 'NSADZI', 'SSERINYA')

# Removed these pops with no individuals:
# BUKASA-LWANABATYRA KOOME-BUGOMBE KOOME-MYENDE SSERINYA-KAFUNA   
pops    = c('BANDA-BANDA', 'BUGALA-BUGOMA', 'BUGALA-LUTOBOKA', 'BUGALA-MWEENA', 
            'BUKASA-NAKIBANGA', 'BUWAMA-BUWAMA', 
            'ENTEBBE-ENTEBBE', 'KAZZI-NABUGABO', 'KIYINDI-KIYINDI',
            'MITYANA-NAAMA', 'NSADZI-KANSAMBWE', 'SSERINYA-BBOSA', 
            'SSERINYA-KASISA')
chrs    = c('2L', '2R', '3L', '3R', 'X')

# Bring in plate info to infer sequencing ID (e.g. ssese42.PE)
plate1 = read.csv("data/ssese_plate1_individuals.csv", header=FALSE)
plate2 = read.csv("data/ssese_plate2_individuals.csv", header=FALSE)
plates = rbind(plate1, plate2)
plates$seq.id = paste0("ssese", 1:nrow(plates), ".PE")
names(plates)[1] = "mosquito.id"

ind.info = read.csv("data/ssese_individual_info.csv")
plates.info = merge(plates, ind.info, by.x="mosquito.id", by.y="mosquito_id")

table(plates.info$island)
table(plates.info$site_name)

inver = read.table("data/inversion_sites.bed")
names(inver) = c("chrom", "start", "end", "name")

plot.basic.stats = function(chr) {

    prefix = paste0("results/chr", chr, ".pass.snp.")

    # --- Inbreeding coefficient ---------------------------------------------------------
    
    for (sp in species) {
    
        het = read.table(paste0(prefix,sp,".het"), header=TRUE)
        het = merge(het, plates.info, by.x="INDV", by.y="seq.id")
        
        p = ggplot(het, aes(factor(island), F)) +
            geom_boxplot() + 
            xlab("Region") + 
            ylab("Inbreeding coefficient, F, by indiv.")
        ggsave(paste0("reports/inbreeding_coef_by_indiv_avg_by_pop_", 
            sp, "_", chr, ".pdf"), plot = p)
        
        p = ggplot(het, aes(color=factor(island), F)) +
            geom_density()
        ggsave(paste0("reports/inbreeding_coef_by_indiv_avg_by_pop_density_", 
            sp, "_", chr, ".pdf"), 
            plot = p, height=5, width=12)
    }
    
    #    # --- HWE ----------------------------------------------------------------------------
    #    
    #    hwe = read.table(paste0(prefix,".hwe"), header=TRUE)
    #    hwe.sig = hwe[hwe$P_HWE < 0.05,]
    #    
    #    p = ggplot(hwe.sig, aes(x=P_HWE)) +
    #        geom_histogram(binwidth=0.001) +
    #        xlab("Hardy Weinberg p-value (p < 0.05)")
    #    ggsave(paste0("reports/hardy_weinberg_pval_hist_", chr, ".pdf"), plot = p)

    #    # --- HWE by population --------------------------------------------------------------
    #    
    #    hwe.pop = hwe[FALSE,]
    #    
    #    for (pop in pops) {
    #        hwe.pop.tmp = read.table(paste0(gsub(".recode", "", prefix), ".", pop, ".hwe"), 
    #                                    header=TRUE)
    #        hwe.pop.tmp$pop = pop
    #        hwe.pop = rbind(hwe.pop, hwe.pop.tmp)
    #    }
    #    
    #    p = ggplot(hwe.pop, aes(P_HWE, color = pop)) +
    #        geom_density() +
    #        xlab("Hardy Weinberg p-value")
    #    ggsave(paste0("reports/hardy_weinberg_pval_density_by_pop_", chr, ".pdf"), plot = p)
    #    
    #    #p = ggplot(hwe.pop, aes(factor(pop), P_HWE)) +
    #    #    geom_boxplot() + 
    #    #    xlab("Region") + 
    #    #    ylab("Hardy Weinberg Test p-value")
    #    #ggsave(paste0("reports/hardy_weinberg_pval_boxplot_by_pop_", chr, ".pdf"), plot = p)
    #
    #    hwe.pop.sig = hwe.pop[hwe.pop$P_HWE < 0.05,]
    #    
    #    # Percentage of sites that are significant
    #    table(hwe.pop.sig$pop) / table(hwe.pop$pop)

    # --- Diversity ----------------------------------------------------------------------
    
    for (sp in species) {
    
        pi = read.table(paste0(prefix, sp, ".sites.pi"), header=TRUE)
        
        p = ggplot(pi, aes(x=PI)) + 
            geom_histogram(binwidth=0.01, fill='white', col='blue') + 
            xlab("Site diversity (pi)")
        ggsave(paste0("reports/diversity_hist_", sp, "_", chr, ".pdf"), plot = p)
        
        pi.win = read.table(paste0(prefix, sp, ".windowed.pi"), header=TRUE)
    
        p = ggplot(pi.win, aes(y=PI, x=BIN_START)) + geom_line() +
            geom_vline(data=subset(inver,chrom==chr), 
                aes(xintercept = start), lty=2, col='red') +
            geom_vline(data=subset(inver,chrom==chr), 
                aes(xintercept = end), lty=2, col='red')
        ggsave(paste0("reports/diversity_windowed_pi_", sp, "_", chr, ".pdf"),
            plot = p, height=5, width=12)
    
    }
    
    # --- Diversity by population --------------------------------------------------------
    
    pi.pop = pi[FALSE,]
    
    for (pop in pops) {
        pi.pop.tmp = read.table(paste0(prefix, pop, ".windowed.pi"), header=TRUE)
        pi.pop.tmp$pop = pop
        pi.pop = rbind(pi.pop, pi.pop.tmp)
    }
    
    ggplot(pi.pop, aes(PI, color = pop)) +
        geom_density()
    
    p = ggplot(pi.pop, aes(pop, PI, color = pop)) + geom_boxplot()
    ggsave(paste0("reports/diversity_windowed_pi_boxplot_", chr, ".pdf"),
            plot = p, height=5, width=12)
    
    #p = ggplot(pi.pop, aes(factor(pop), PI)) +
    #    geom_boxplot() + xlab("Region") + 
    #    ylab("Site diversity (pi)")
    #ggsave(paste0("reports/diversity_boxplot_by_pop_", chr, ".pdf"), plot = p)
    
    # --- Tajima's D ---------------------------------------------------------------------
    
    for (sp in species) {

        taj.d = read.table(paste0(prefix, sp, ".Tajima.D"), header=TRUE)
        
        p = ggplot(taj.d) +
            geom_point(aes(x=BIN_START, y=TajimaD)) +
            xlab("Bin Start") +
            ylab("Tajima's D") + 
            geom_line(aes(x=BIN_START, y=N_SNPS/1000), col='red')
        ggsave(paste0("reports/tajimaD_across_genome_", sp, "_", chr, ".pdf"), plot = p)
        
        # Get top n% of Tajima D
    
        tajd.top5 = taj.d[taj.d$TajimaD > quantile(taj.d$TajimaD,prob=0.95, na.rm=TRUE),]
        
        # Write top 5% of Tajima D SNPs to file
        write.table(tajd.top5, file=paste0(prefix, sp, ".Tajima.D.top5"))
    
    }
    
    # --- Relatedness --------------------------------------------------------------------
    
    for (sp in species) {

        relate1 = read.table(paste0(prefix, sp, ".relatedness"),  header=TRUE)
        relate2 = read.table(paste0(prefix, sp, ".relatedness2"), header=TRUE)
        
        # Add population of individual 1
        relate1 = merge(relate1, plates.info[,c(2,3,4)], by.x="INDV1", by.y="seq.id")
        names(relate1)[4] = "INDV1.island"
        names(relate1)[5] = "INDV1.site"
        relate2 = merge(relate2, plates.info[,c(2,3,4)], by.x="INDV1", by.y="seq.id")
        names(relate2)[8] = "INDV1.island"
        names(relate2)[9] = "INDV1.site"
        
        # And of individual 2
        relate1 = merge(relate1, plates.info[,c(2,3,4)], by.x="INDV2", by.y="seq.id")
        names(relate1)[6] = "INDV2.island"
        names(relate1)[7] = "INDV2.site"
        relate2 = merge(relate2,plates.info[,c(2,3,4)], by.x="INDV2", by.y="seq.id")
        names(relate2)[10] = "INDV2.island"
        names(relate2)[11] = "INDV2.site"
        
        # Boolean, are they from same population
        relate1$is.same.pop = relate1$INDV1.site == relate1$INDV2.site
        relate2$is.same.pop = relate2$INDV1.site == relate2$INDV2.site
        
        relate1.same = relate1[relate1$is.same.pop,]
        relate2.same = relate2[relate2$is.same.pop,]
        
        p = ggplot(relate1.same, aes(RELATEDNESS_AJK, color = INDV1.site)) +
            geom_density() +
            xlab("Relatedness (Unadjusted Ajk Statistic)") +
            ylab("Density")
        ggsave(paste0("reports/relatedness_ajk_density_by_pop_", sp, "_", chr, ".pdf"), 
            plot = p)
        
        p = ggplot(relate2.same, aes(RELATEDNESS_PHI, color = INDV1.site)) +
            geom_density() + 
            xlim(c(-1,0.52)) +
            xlab("Relatedness (phi statistic)") +
            ylab("Density")
        ggsave(paste0("reports/relatedness_phi_density_by_pop_", sp, "_", chr, ".pdf"), 
            plot = p)
        
    }

}

lapply(chrs, plot.basic.stats)
