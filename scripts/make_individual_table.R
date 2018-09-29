#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Make individual info table and site info table (with GPS coords)
# ----------------------------------------------------------------------------------------

library(xtable)

depth = read.table("reports/all.pass.snp.flt.idepth", header=TRUE)

mean(depth$MEAN_DEPTH)

ind.info = read.table("data/ssese_individual_info_simple.txt")

names(ind.info) = c("field.ID", "seq.ID", "island", "site", "sex", "species")

ind.d = merge(ind.info, depth, by.x="seq.ID", by.y="INDV")
ind.d$ID.num = as.numeric(gsub("ssese(.*)\\.PE", "\\1", ind.d$seq.ID))
ind.d = ind.d[order(ind.d$ID.num),]

ind.d$seq.ID = gsub("\\.PE", "", gsub("ssese", 'LVB2015-', ind.d$seq.ID))

ind.d = ind.d[,c(1:4,8)]

# ----------------------------------------------------------------------------------------

names(ind.d) = c("ID", "Field ID", "Island", "Site", "Mean depth")

ind.d$Island = paste0(substr(ind.d$Island, 1, 1),
                      tolower(substr(ind.d$Island, 2, nchar(ind.d$Island))))
ind.d$Site   = paste0(substr(ind.d$Site, 1, 1),
                      tolower(substr(ind.d$Site, 2, nchar(ind.d$Site))))

ind.d[["Mean depth"]] = signif(as.numeric(ind.d[["Mean depth"]]), 3)

# Fix spelling
ind.d[ind.d$Island == "Kazzi",  ]$Island = "Kaazi"
ind.d[ind.d$Island == "Mityana",]$Island = "Wamala"

# ----------------------------------------------------------------------------------------

xt = xtable(ind.d, digits=c(0,0,0,0,0,2), align='lll|ll|r',
    caption="List of individuals included in study with mean depth of sequencing coverage.",
    label = "table:ind_info")

out.file = "reports/individual_table.tex"
sink(out.file)

cat("\\documentclass{article}",
    "\\usepackage{graphicx}",
    "\\usepackage{longtable}",
    "\\DeclareGraphicsExtensions{.pdf}",
    "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}",
    "\\usepackage{caption}",
    "\\captionsetup[table]{font=normalsize}",
    "\\begin{document}", sep="\n")

print(xt, type="latex",
    include.rownames=FALSE,
    tabular.environment="longtable", floating = FALSE,
    size="\\fontsize{9pt}{10pt}\\selectfont",
    sanitize.colnames.function = identity,
    sanitize.text.function = identity,
    caption.placement = "top",
    hline.after = c(0,0, nrow(ind.d)))

cat("\\end{document}", sep="\n")

sink()

# ========================================================================================

gps = read.table("data/site_gps_coordinates.txt")

names(gps) = c("Location", "Latitude", "Longitude")

gps = gps[gps$Location != "Entebbe",]

gps$Location = gsub("_", " - ", gps$Location)

site.ct = c(table(ind.d$Island), table(paste(ind.d$Island, ind.d$Site, sep=" - ")))
site.ct.df = data.frame(Location=names(site.ct), "Sample Count"=site.ct)

# Fix spelling
gps[gps$Location == "Kazzi",  ]$Location = "Kaazi"
gps[gps$Location == "Mityana",]$Location = "Wamala"

gps = merge(gps, site.ct.df)
names(gps)[4] = "Sample Count"

xt = xtable(gps, digits=c(0,0,5,5,0), align='ll|rr|c',
    caption="Sampling sites and coordinates.",
    label = "table:site_coords")

out.file = "reports/site_gps_table.tex"
sink(out.file)

cat("\\documentclass{article}",
    "\\usepackage{graphicx}",
    "\\usepackage{longtable}",
    "\\DeclareGraphicsExtensions{.pdf}",
    "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}",
    "\\begin{document}", sep="\n")

print(xt, type="latex",
    include.rownames=FALSE,
    tabular.environment="longtable", floating = FALSE,
    size="\\fontsize{9pt}{10pt}\\selectfont",
    sanitize.colnames.function = identity,
    sanitize.text.function = identity,
    caption.placement = "top",
    hline.after = c(0,0, nrow(gps)))

cat("\\end{document}", sep="\n")

sink()
