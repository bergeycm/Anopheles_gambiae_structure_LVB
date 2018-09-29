#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Parse and plot IBDNe results
# ----------------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)

args = commandArgs(trailingOnly=TRUE)
chr = args[1]    # E.g. 3

library(ggplot2)
library(grid)
library(gridExtra)
library(scales)

ibdne.out = list.files(path="results/",
    pattern=paste0("chr", chr, ".pass.snp.flt.*.ne.ne"), full.names=TRUE)

ibdne.res = do.call(rbind, lapply(ibdne.out, function (in.file) {

    tmp = read.table(in.file, header=TRUE)
    tmp$POP = strsplit(in.file, split="\\.")[[1]][5]
    tmp
}))

island.sites   = c("BANDA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")
mainland.sites = c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")

ibdne.res$is.island = FALSE
ibdne.res[ibdne.res$POP %in% island.sites,]$is.island = TRUE

ibdne.res$POP = gsub("BUGALAML", "BUGALA (M)", ibdne.res$POP)
ibdne.res$POP = gsub("BUGALAIS", "BUGALA (I)", ibdne.res$POP)

simple.cap = function(x) {
    s = strsplit(x, " ")[[1]]
    lstr = paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")

    gsub("\\(i\\)", "(I)", gsub("\\(m\\)", "(M)", lstr))
}

ibdne.res$POP = as.vector(sapply(ibdne.res$POP, simple.cap))

# ----------------------------------------------------------------------------------------

p = ggplot(ibdne.res, aes(GEN, log(NE, base=10),
        col=POP)) +
    geom_line() +
    xlim(c(0,5000)) +  ylim(c(0,25)) +
    xlab("Generation") +
    ylab(expression('log'[10]~'(N'[e]*')')) +
    theme_bw() +
    scale_color_discrete(name="Population")

ggsave(p, file=paste0("reports/IBDNe_by_pop.", chr, ".pdf"), height=4, width=10)

# ----------------------------------------------------------------------------------------

p = ggplot(ibdne.res, aes(GEN, log(NE, base=10),
        col=is.island, group=POP)) +
    geom_line() +
    xlab("Generation") +
    ylab(expression('log'[10]~'(N'[e]*')')) +
    theme_bw(base_size=15) +
    scale_color_manual(values=c("#F8766D", "#00BFC4"),
                       name="Site type",
                       breaks=c(TRUE, FALSE),
                       labels=c("Island", "Mainland")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black")) +
    scale_x_continuous(limits = c(0,3000), expand = c(0, 0)) +
    scale_y_continuous(limits=c(0,25), expand = c(0, 0))

ggsave(p, file=paste0("reports/IBDNe_by_site_type.", chr, ".pdf"), height=4, width=10)

# ----------------------------------------------------------------------------------------

ibdne.res$year = ibdne.res$GEN / 11

p = ggplot(ibdne.res, aes(year, NE,
        col=is.island, group=POP)) +
    geom_line() +
    scale_x_continuous(trans='log10',
            #limits = c(1,3000), expand = c(0, 0),
            labels = trans_format('log10',math_format(10^.x))) +
    scale_y_continuous(trans='log10',
            limits=c(1,10^25), expand = c(0, 0),
            labels = trans_format('log10',math_format(10^.x))) +
    annotation_logticks() +
    xlab("Years before present") +
    ylab(expression(N[e])) +
    theme_bw(base_size=15) +
    scale_color_manual(values=c("#F8766D", "#00BFC4"),
                       name="Site type",
                       breaks=c(TRUE, FALSE),
                       labels=c("Island", "Mainland")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"))

ggsave(p, file=paste0("reports/IBDNe_by_site_type.", chr, ".log.pdf"),
    height=4, width=10)

# ----------------------------------------------------------------------------------------

p = ggplot(ibdne.res, aes(GEN, log(NE, base=10),
        col=is.island, group=POP, lty=POP=="Mityana")) +
    geom_line() +
    xlab("Generation") +
    ylab(expression('log'[10]~'(N'[e]*')')) +
    theme_bw(base_size=15) +
    scale_color_manual(values=c("#F8766D", "#00BFC4"),
                       name="Site type",
                       breaks=c(TRUE, FALSE),
                       labels=c("Island", "Mainland")) +
    scale_linetype_manual(values=c(1,2),
                          name="",
                          breaks=c(TRUE, FALSE),
                          labels=c("Wamala", "Other Sites")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black")) +
    scale_x_continuous(limits = c(0,3000), expand = c(0, 0)) +
    scale_y_continuous(limits=c(0,25), expand = c(0, 0))

ggsave(p, file=paste0("reports/IBDNe_by_site_type_Mityana_highlighted.", chr, ".pdf"),
    height=4, width=10)

# ----------------------------------------------------------------------------------------

ibdne.res$year = ibdne.res$GEN / 11

p = ggplot(ibdne.res, aes(year, NE,
        col=is.island, group=POP, lty=POP=="Mityana")) +
    geom_line() +
    scale_x_continuous(trans='log10',
            #limits = c(1,3000), expand = c(0, 0),
            labels = trans_format('log10',math_format(10^.x)),
            sec.axis = sec_axis(~ . * 11,
                                name = "Generations before present",
                                breaks=c(1,5,10,50,100,500,1000))) +
    scale_y_continuous(trans='log10',
            limits=c(1,10^25), expand = c(0, 0),
            labels = trans_format('log10',math_format(10^.x))) +
    annotation_logticks() +
    xlab("Years before present") +
    ylab(expression(N[e])) +
    theme_bw(base_size=15) +
    scale_color_manual(values=c("#F8766D", "#00BFC4"),
                       name="Site type",
                       breaks=c(TRUE, FALSE),
                       labels=c("Island", "Mainland")) +
    scale_linetype_manual(values=c(1,2),
                          name="",
                          breaks=c(TRUE, FALSE),
                          labels=c("Wamala", "Other Sites")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"))

ggsave(p, file=paste0("reports/IBDNe_by_site_type_Mityana_highlighted.", chr, ".log.pdf"),
    height=4, width=10)

# ----------------------------------------------------------------------------------------

p = ggplot(ibdne.res, aes(GEN, log(NE, base=10),
        col=is.island, group=POP, lty=(POP %in% c("Bugala (M)", "Kaazi", "Wamala")))) +
    geom_line() +
    xlim(c(0,3000)) +  ylim(c(0,25)) +
    xlab("Generation") +
    ylab(expression('log'[10]~'(N'[e]*')')) +
    theme_bw() +
    scale_color_manual(values=c("#F8766D", "#00BFC4"),
                       name="Site type",
                       breaks=c(TRUE, FALSE),
                       labels=c("Island", "Mainland")) +
    scale_linetype_manual(values=c(1,2),
                          name="",
                          breaks=c(TRUE, FALSE),
                          labels=c("Wamala", "Other Sites"))

ggsave(p, file=paste0("reports/IBDNe_by_site_type_several_highlighted.", chr, ".pdf"),
    height=4, width=10)
