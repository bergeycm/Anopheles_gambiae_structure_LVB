#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

library(ggplot2)

# ----------------------------------------------------------------------------------------
# --- Parse results of inferring demography from IBS
# ----------------------------------------------------------------------------------------

parse.output = function (site) {

    in.file = paste0("results/IBS_inferred_size_", site, ".txt")

    pop.sizes = system(paste('grep "Population sizes" ', in.file,
        " | tail -n1 | cut -d':' -f2 | tr -d ','"), intern=TRUE)
    pop.sizes = as.numeric(strsplit(pop.sizes, split=' ')[[1]])

    change.times = system(paste('grep "Population size change times" ', in.file,
        " | tail -n1 | cut -d':' -f2 | tr -d ','"), intern=TRUE)
    change.times = as.numeric(strsplit(change.times, split=' ')[[1]])

    ## Duplicate last time, since it has two sizes associated with it
    ## (ancestral and size after first change)
    #change.times = c(change.times, change.times[length(change.times)])

    change.times = change.times[sort(rep(1:length(change.times), 2))]
    pop.sizes = pop.sizes[sort(rep(1:length(pop.sizes), 2))][-c(1,length(pop.sizes) * 2)]

    data.frame(site, change.times, pop.sizes)
}

is.sites = c("BANDA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")

sites = c(is.sites, ml.sites)

res.1pop = do.call(rbind, lapply(sites, parse.1pop.output))

#p = ggplot(res.1pop, aes(change.times, pop.sizes, col=site %in% is.sites, group=site)) + geom_line()
p = ggplot(res.1pop, aes(change.times, pop.sizes, col=site)) + geom_line()
ggsave(p, file="tmp.1pop.res.pdf", width=14, height=7)

# ----------------------------------------------------------------------------------------

# Get all combos of sites
sites.s = sort(sites)
site.combos = cbind(sites.s[combn(1:length(sites.s), 2)[1,]],
                    sites.s[combn(1:length(sites.s), 2)[2,]])

site.combos.str = apply(site.combos, 1, function(x) {paste(x[1], x[2], sep="_vs_")})

res.2pop = do.call(rbind, lapply(site.combos.str, parse.output))

#p = ggplot(res.1pop, aes(change.times, pop.sizes, col=site %in% is.sites, group=site)) + geom_line()
p = ggplot(res.2pop, aes(change.times, pop.sizes, col=site)) + geom_line()
ggsave(p, file="tmp.2pop.res.pdf", width=14, height=7)
