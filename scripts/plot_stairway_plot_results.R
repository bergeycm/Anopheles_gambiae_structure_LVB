#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Plot stairway... plot results
# ----------------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)
#options(scipen=999)

library(ggplot2)
library(grid)
library(gridExtra)
library(scales)

island.sites   = c("BANDA", "BUGALA", "BUGALAIS",
                   "BUKASA", "KOOME", "NSADZI", "SSERINYA")
mainland.sites = c("BUGALAML", "BUWAMA", "ENTEBBE",
                   "KAZZI", "KIYINDI", "MITYANA")

sites = gsub(".*\\.", "", list.files(path = "data", pattern = "stairway_plot_island"))

all.sw = do.call(rbind, lapply(sites, function (x) {
    this.sw.file = paste0("data/stairway_plot_island.", x, "/plots/Stairway plot island.",
        x, ".final.summary")
    this.sw = read.table(this.sw.file, header=TRUE)
    this.sw$site = x
    this.sw$is.island = c("Mainland", "Island")[1 + this.sw$site %in% island.sites]
    this.sw
}))

# Remove unsplit Bugala
all.sw = all.sw[all.sw$site != "BUGALA",]
all.sw[all.sw$site == "BUGALAML",]$site = "BUGALA (M)"
all.sw[all.sw$site == "BUGALAIS",]$site = "BUGALA (I)"
mainland.sites[mainland.sites == "BUGALAML"] = "BUGALA (M)"
island.sites  [island.sites   == "BUGALAIS"] = "BUGALA (I)"

# Split into two graphs
all.sw.split = split(all.sw, all.sw$is.island)

simple.cap = function(x) {
    s = strsplit(x, " ")[[1]]
    lstr = paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")

    gsub("\\(i\\)", "(I)", gsub("\\(m\\)", "(M)", lstr))
}

all.sw$site = as.vector(sapply(all.sw$site, simple.cap))
island.sites = as.vector(sapply(island.sites, simple.cap))
mainland.sites = as.vector(sapply(mainland.sites, simple.cap))

options(scipen=999)

p = ggplot(all.sw, aes(year, Ne_median, col=is.island, fill=site)) +
        geom_line() +
        #geom_ribbon(aes(ymin=Ne_2.5.,ymax=Ne_97.5.), alpha=0.2, color=NA) +
        facet_grid(is.island ~ .) +
        scale_x_continuous(trans='log10',
            limits = c(1000, 100000),
            labels = trans_format('log10',math_format(10^.x)),
            sec.axis = sec_axis(~ . * 11,
                                name   = "Generations before present",
                                breaks = c(10000,50000,100000,500000),
                                labels = paste(c(10,50,100,500), "thousand"))) +
        #scale_y_continuous(trans='log10',
        #    limits = c(floor(min(all.sw$Ne_2.5)/1000) * 1000,
        #               ceiling(max(all.sw$Ne_97.5)/1000000) * 1000000)) +
        scale_y_continuous(trans='log10',
            limits = c(10e4, 10e6),
            labels = trans_format('log10',math_format(10^.x))) +
        annotation_logticks() +
        xlab("Years before present") +
        ylab(expression(N[e])) +
        scale_fill_discrete(guide=FALSE) +
        scale_color_manual(guide=FALSE, values=c("#00BFC4", "#F8766D")) +
        theme_bw(base_size=15) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(color = "black"),
              axis.text.x = element_text(angle = 0),
              strip.background = element_blank(),
              strip.text = element_blank()) +
        geom_hline(yintercept=0)

ggsave(p, file="reports/stairway.plot.superimposed.pdf", height=6, width=9)
