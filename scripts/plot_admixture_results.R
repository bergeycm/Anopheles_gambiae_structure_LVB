#!/usr/bin/Rscript

library(ggplot2)
library(RColorBrewer)

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly = TRUE)
adm.prefix = args[1]
fam.file   = args[2]

# E.g.:
# adm.prefix = "data/ssese_with_ag1000g/adm_subsets/CLUMPAK_withM."
# fam.file = "data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.thinned_for_ADMIXTURE.withM.fam"

# ========================================================================================
# === Analyze and plot ADMIXTURE results
# ========================================================================================

adm.inds = read.table(fam.file)$V1

ind.info = read.table("data/ssese_individual_info_simple_bugala_split.txt")

# Bring in plate info to infer sequencing ID (e.g. ssese42.PE)
plates = read.table("data/sample_to_seq_id_mapping.txt")
names(plates) = c("mosquito.id", "seq.id")

# Reduce individual info to just sequenced animals
ind.info = ind.info[ind.info$V1 %in% plates$mosquito.id,]

ind.info.ag1000g = read.table("data/ag1000g.phase1.ar3/samples.all.txt",
    sep="\t", quote="", header=TRUE)

# ----------------------------------------------------------------------------------------

make.adm.plot = function(k) {

    this.adm.file = paste0(adm.prefix, k, ".Q")

    adm = read.table(this.adm.file)
    names(adm) = paste0("ADM_", 1:as.numeric(k))
    adm = cbind(adm.inds, adm)

    # Don't plot SE if file isn't present
    do.plot.se = FALSE
    adm.se.file = paste0(this.adm.file, "_se")
    if (file.exists(adm.se.file)) {
        adm.se = read.table(adm.se.file)
        names(adm.se) = paste0("SE_", 1:as.numeric(k))
        plot.se = TRUE
        adm = cbind(adm, adm.se)
    }

    adm$pops = do.call(c, lapply(adm$adm.inds, function(x) {
        pop.ss = as.character(ind.info[ind.info$V2 == x, 3])
        pop.ag = paste(
            as.character(ind.info.ag1000g[ind.info.ag1000g$ox_code == x,]$country),
            as.character(ind.info.ag1000g[ind.info.ag1000g$ox_code == x,]$m_s),
            sep="_")
        if (length(pop.ss)) {
            return(pop.ss)
        } else if (length(pop.ag)) {
            return(pop.ag)
        } else if (grepl("AD", x)) {
            return("AD")
        } else {
            return(x)
        }
    }))

    # Remove poorly represented pops
    adm = adm[!(adm$pops %in% c("Guinea_M/S", "Guinea-Bissau_M", "Guinea-Bissau_M/S")),]

    adm$pop.expanded = adm$pops

    adm$pops[adm$pops %in% c("BANDA", "BUKASA", "BUGALA", "BUGALAIS",
        "NSADZI", "SSERINYA")] = "LVB - Ssese Islands"
    adm$pops[adm$pops %in% c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")] =
        "LVB - Mainland"

    # Change Uganda to "Uganda - Ag1000G"
    adm$pops[adm$pops == "Uganda_S"] = "Uganda - Ag1000G"

    ordered.pops = names(table(adm$pops))

    #adm$pops = factor(adm$pops,
    #                  levels = c(grep("LVB", ordered.pops, invert=TRUE, value=TRUE),
    #                             grep("LVB", ordered.pops, value=TRUE)))

    # Roughly by longitude
    if ("Angola_M" %in% ordered.pops) {
        pop.longitude = c("Angola_M", "Burkina Faso_M", "Guinea-Bissau_S", "Guinea_S",
                          "Burkina Faso_S", "Cameroon_S",
                          "Uganda - Ag1000G", "LVB - Mainland", "LVB - Ssese Islands",
                          "Gabon_S", "Kenya_S")
    } else {
        pop.longitude = c("Guinea-Bissau_S", "Guinea_S",
                          "Burkina Faso_S", "Cameroon_S",
                          "Uganda - Ag1000G", "LVB - Mainland", "LVB - Ssese Islands",
                          "Gabon_S", "Kenya_S")
    }
    adm$pops = factor(adm$pops,
        levels = pop.longitude)

    # Re-order to match now that rearranging done
    ordered.pops = levels(adm$pops)

    adm$is.ssese = grepl("LVB", adm$pops)

    ###pal = brewer.pal(n = length(table(adm$pops)) - 2, "Set1")
    ###new.pal = rep(NA, length(ordered.pops))

    #### Fill in non-LVB colors
    ###new.pal[grep ("LVB", ordered.pops, invert=TRUE)] = pal
    #### Set Ssese islands to blue and mainland to red
    ###new.pal[grep ("Ssese Island", ordered.pops)] = "#00BFC4"
    ###new.pal[grep ("Mainland", ordered.pops)]     = "#F8766D"

    num.ind = nrow(adm)

    order.list = append(
                    append(
                        list(adm$pops),
                        lapply(1:10, function (x) { adm[[paste0("ADM_", x)]] })),
                    list(adm$adm.inds)
        )
    # Remove NULL for missing components (greater than current K)
    order.list = order.list[!sapply(order.list, is.null)]

    # Nice trick to "convert" list to ellipsis from http://stackoverflow.com/a/18545884
    new.order = do.call(function(...) order(...), order.list)

    # Figure out where to draw lines between islands
    cats = as.character(adm[new.order,]$pops)
    last.for.pop = which(cats != c(cats[-1], NA))

    pal = brewer.pal(as.numeric(k), "Set3")

    # Write palette to file so maps can use them
    write.table(data.frame(colors=pal), paste0(adm.prefix, k, ".colors.txt"),
        quote=FALSE, sep="\t", row.names=FALSE)

    for (pdf.out in c(paste0(adm.prefix, k, ".pdf"), paste0(adm.prefix, k, ".boot.pdf"))) {

        pdf(file=pdf.out, width=15, height=5)

            par(mar=c(10,2,4,2))

            mp = barplot(t(as.matrix(adm[new.order, grep("ADM",names(adm))])),
                col=pal,
                xlab="", ylab="",
                yaxt = 'n',
                border=NA, xaxt='n',
                space = c(0, 0))

            title(ylab=paste0("k=", k), line=0, cex.lab=1.8)

            if (grepl("boot", pdf.out)) {

                for (this.k in 1:as.numeric(k)) {

                    # Tally up cumulative ancestry explained by k values from 1 to this k
                    anc.explained = rowSums(cbind(rep(0, nrow(adm)),
                                            adm[new.order,2:(1 + this.k)]))

                    if (do.plot.se) {
                        this.SE = adm[new.order,][[paste0("SE_", this.k)]]

                        segments(mp, anc.explained - this.SE * 2,
                                 mp, anc.explained + this.SE * 2, lwd = 0.01)
                    }
                }
            }

            abline(v=mp[last.for.pop] + 0.5 * mp[1], lwd=3)

            cat.midpoints = c(last.for.pop, nrow(adm)) - 0.5 * table(adm$pops)
            axis(1, at=mp[cat.midpoints],
                labels=gsub("_M$", "\n[coluzzii]", gsub("_S$", "", ordered.pops)),
                las=2, cex.axis=1, lwd=0)

        dev.off()
    }

    # ------------------------------------------------------------------------------------
    # --- Make scaled plot for text
    # ------------------------------------------------------------------------------------

    pdf(file=gsub("pdf", "fortext.pdf", pdf.out), width=7, height=4)

        par(mar=c(10,2,4,2))

        mp = barplot(t(as.matrix(adm[new.order, grep("ADM",names(adm))])),
            col=pal,
            xlab="", ylab="",
            yaxt = 'n',
            border=NA, xaxt='n',
            space = c(0, 0),
            cex.axis=1, cex.names=1,
            xaxs="i", yaxs="i")

        abline(v=mp[last.for.pop] + 0.5 * mp[1], lwd=3)

        abline(v=0, lwd=2)
        abline(v=mp[length(mp)] + 0.5, lwd=2)
        abline(h=0, lwd=2)
        abline(h=1, lwd=2)

        cat.midpoints = c(last.for.pop, nrow(adm)) - 0.5 * table(adm$pops)
        axis(1, at=mp[cat.midpoints],
            labels=gsub("_M$", "\n[coluzzii]", gsub("_S$", "", ordered.pops)),
            las=2, cex.axis=0.8, lwd=0)

    dev.off()

    # ------------------------------------------------------------------------------------
    # --- Now make Ssese-only map
    # ------------------------------------------------------------------------------------

    adm.ss = adm[adm$is.ssese,]

    # Fix names
    adm.ss$pop.expanded[adm.ss$pop.expanded == "KAZZI"] = "KAAZI"
    adm.ss$pop.expanded[adm.ss$pop.expanded == "MITYANA"] = "WAMALA"

    ordered.pops = names(table(adm.ss$pop.expanded))

    # Mainland then island, roughly by longitude
    pop.longitude = c("WAMALA", "KAAZI", "BUWAMA", "KIYINDI",
        "BUGALAML", "BUGALAIS", "SSERINYA", "BANDA", "BUKASA", "NSADZI")

    adm.ss$pop.expanded = factor(adm.ss$pop.expanded,
        levels = pop.longitude)

    # Re-order to match now that rearranging done
    ordered.pops = levels(adm.ss$pop.expanded)

    num.ind = nrow(adm.ss)

    order.list = append(
                    append(
                        list(adm.ss$pop.expanded),
                        lapply(1:10, function (x) { adm.ss[[paste0("ADM_", x)]] })),
                    list(adm.ss$adm.inds)
        )
    # Remove NULL for missing components (greater than current K)
    order.list = order.list[!sapply(order.list, is.null)]

    # Nice trick to "convert" list to ellipsis from http://stackoverflow.com/a/18545884
    new.order = do.call(function(...) order(...), order.list)

    # Figure out where to draw lines between islands
    cats = as.character(adm.ss[new.order,]$pop.expanded)
    last.for.pop = which(cats != c(cats[-1], NA))

    pal = brewer.pal(as.numeric(k), "Set3")

    for (pdf.out in c(paste0(adm.prefix, "LVB_only.", k, ".pdf"), paste0(adm.prefix, "LVB_only.", k, ".boot.pdf"))) {

        pdf(file=pdf.out, width=7, height=5)

            par(mar=c(8,2,4,2))

            mp = barplot(t(as.matrix(adm.ss[new.order, grep("ADM",names(adm.ss))])),
                col=pal,
                xlab="", ylab="",
                yaxt = 'n',
                border=NA, xaxt='n',
                space = c(0, 0))

            title(ylab=paste0("k=", k), line=0, cex.lab=1.8)

            if (grepl("boot", pdf.out)) {

                for (this.k in 1:as.numeric(k)) {

                    # Tally up cumulative ancestry explained by k values from 1 to this k
                    anc.explained = rowSums(cbind(rep(0, nrow(adm.ss)),
                                            adm.ss[new.order,2:(1 + this.k)]))

                    if (do.plot.se) {
                        this.SE = adm.ss[new.order,][[paste0("SE_", this.k)]]

                        segments(mp, anc.explained - this.SE * 2,
                                 mp, anc.explained + this.SE * 2, lwd = 0.01)
                    }
                }
            }

            abline(v=mp[last.for.pop] + 0.5 * mp[1], lwd=1)

            # Add thick line between ML and island
            abline(v=mp[last.for.pop[5]] + 0.5 * mp[1], lwd=3)

            simpleCap = function(s) {
                paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
                    sep="", collapse=" ")
            }

            ss.labels = sapply(ordered.pops, simpleCap)
            ss.labels[ss.labels == "Bugalaml"] = "Bugala (M)"
            ss.labels[ss.labels == "Bugalais"] = "Bugala (I)"

            cat.midpoints = c(last.for.pop, nrow(adm.ss)) - 0.5 * table(adm.ss$pop.expanded)
            axis(1, at=mp[cat.midpoints],
                labels=ss.labels,
                las=2, cex.axis=1, lwd=0)

        dev.off()
    }
}

lapply(2:10, make.adm.plot)
