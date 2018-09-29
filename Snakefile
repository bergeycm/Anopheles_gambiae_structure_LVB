# snakemake --jobs 8 --cluster "qsub -V -M cbergey@nd.edu -m abe -pe smp={threads}{params.mem}"

# ----------------------------------------------------------------------------------------
# --- Variables
# ----------------------------------------------------------------------------------------

# These are not all islands, despite the variable name. This includes both mainland and
#   island sites.
# Removed the following size zero islands:
# ENTEBBE KOOME
ISLANDS = ["BANDA", "BUGALA", "BUKASA", "BUWAMA", "KAZZI",
           "KIYINDI", "MITYANA", "NSADZI", "SSERINYA"]
ISLANDS_SPLIT_BUGALA = ["BANDA", "BUGALAIS", "BUGALAML", "BUKASA", "BUWAMA", "KAZZI",
           "KIYINDI", "MITYANA", "NSADZI", "SSERINYA"]
GROUPINGS = ["SMALL_ISLANDS", "ISLANDS", "MAINLAND"]
# Removed the following size zero populations:
# BUKASA-LWANABATYRA ENTEBBE-ENTEBBE KOOME-BUGOMBE KOOME-MYENDE SSERINYA-KAFUNA
POPS = ["BANDA-BANDA", "BUGALA-BUGOMA", "BUGALA-LUTOBOKA", "BUGALA-MWEENA",
        "BUKASA-NAKIBANGA", "BUWAMA-BUWAMA", "KAZZI-NABUGABO", "KIYINDI-KIYINDI",
        "MITYANA-NAAMA", "NSADZI-KANSAMBWE", "SSERINYA-BBOSA", "SSERINYA-KASISA"]
CHRS,  = glob_wildcards("data/chr{chr}.pass.snp.vcf")
CHR_THREE  = ["3L", "3R"]
VCF    = expand("data/chr{chrs}.pass.snp.vcf", chrs=CHRS)
VCF_GZ = expand("{vcf}.gz", vcf=VCF)

# --- Site pairs

SITE_PAIRS = []

for i in range(0, len(ISLANDS_SPLIT_BUGALA)):
    for j in range(i+1 ,len(ISLANDS_SPLIT_BUGALA)):
        SITE_PAIRS.append(ISLANDS_SPLIT_BUGALA[i] + '_' + ISLANDS_SPLIT_BUGALA[j])

# --- Site pairs for dadi. Same as above but separated by a period

DADI_SITE_PAIRS = []

for i in range(0, len(ISLANDS_SPLIT_BUGALA)):
    for j in range(i+1 ,len(ISLANDS_SPLIT_BUGALA)):
        DADI_SITE_PAIRS.append(ISLANDS_SPLIT_BUGALA[i] + '.' + ISLANDS_SPLIT_BUGALA[j])

# Function that maps XP-EHH output file to its input haps files
XPEHH_INPUT = {}
XPEHH_INPUT_JP = {}
for c in CHRS:
    for i in range(0, len(ISLANDS_SPLIT_BUGALA)):
        for j in range(i+1 ,len(ISLANDS_SPLIT_BUGALA)):
            ISL1 = ISLANDS_SPLIT_BUGALA[i]
            ISL2 = ISLANDS_SPLIT_BUGALA[j]
            XPEHH_OUT = "results/selscan/xp-ehh." + ISL1 + "-" + ISL2 + \
                "." + c + ".xpehh.out"
            INHAP1 = "data/chr" + c + ".pass.snp.phased." + ISL1 + ".haps"
            INHAP2 = "data/chr" + c + ".pass.snp.phased." + ISL2 + ".haps"
            XPEHH_INPUT[XPEHH_OUT] = [INHAP1, INHAP2]
            INHAP1_JP = "data/xpehh_input_jp/chr" + c + ".pass.snp.phased." + ISL1 + ".jp.haps.transpose"
            INHAP2_JP = "data/xpehh_input_jp/chr" + c + ".pass.snp.phased." + ISL2 + ".jp.haps.transpose"
            XPEHH_OUT_JP = "results/xpehh_jp/xp-ehh." + ISL1 + "-" + ISL2 + \
                "." + c + ".xpehh.out"
            XPEHH_INPUT_JP[XPEHH_OUT_JP] = [INHAP1_JP, INHAP2_JP]

# Function that maps IBS-demog-inference output file to its input popdata files
IBS_INPUT = {}
for i in range(0, len(ISLANDS_SPLIT_BUGALA)):
    for j in range(i+1 ,len(ISLANDS_SPLIT_BUGALA)):
        ISL1 = ISLANDS_SPLIT_BUGALA[i]
        ISL2 = ISLANDS_SPLIT_BUGALA[j]
        IBS_OUT = "results/IBS_inferred_size_" + ISL1 + "_vs_" + ISL2 + \
            ".demographic_history.txt"
        INPOPDATA1 = "results/IBS." + ISL1 + ".popdata"
        INPOPDATA2 = "results/IBS." + ISL2 + ".popdata"
        IBS_INPUT[IBS_OUT] = [INPOPDATA1, INPOPDATA2]

DADI_MODES=['island']
DADI_OUTPUT=[]

for site in ISLANDS_SPLIT_BUGALA:
    DADI_OUTPUT.append('island-' + site)

DADI_OUTPUT_PERIOD = [s.replace('-', '.') for s in DADI_OUTPUT]

# ----------------------------------------------------------------------------------------
# --- Make all
# ----------------------------------------------------------------------------------------

rule all:
    input:
        # make_island_lists:
        expand("data/ssese.samples.is-{islands}.txt", islands=ISLANDS),
        # make_site_lists:
        expand("data/ssese.samples.site-{pops}.txt",  pops=POPS),
        # map_sample_to_seq_id:
        "data/sample_to_seq_id_mapping.txt",
        # make_seqid_lists:
        expand("data/ssese.seqids.is-{islands}.txt", islands=ISLANDS),
        expand("data/ssese.seqids.site-{pops}.txt",  pops=POPS),
        # make_simple_info_list:
        "data/ssese_individual_info_simple.txt",
        # combine_het_inv:
        "data/inversion_heterochromatin.bed",
        # make_het_inv_sets:
        "data/inversion_simple.set",
        "data/heterochromatin.set",
        # combine_chr_vcfs:
        "data/all.pass.snp.vcf.gz",
        # compute_indiv_missingness:
        "reports/all.pass.snp.imiss",
        # find_high_missing_indivs:
        "reports/inds.high.missing.txt",
        # filter_vcf:
        expand("data/chr{chr}.pass.snp.flt.vcf.gz", chr=CHRS),
        # vcf_to_bed:
        expand("data/chr{chr}.pass.snp.flt.bed", chr=CHRS),
        # ss_ag_vcf_to_bed_missingtoref:
        expand("data/ssese_with_ag1000g.missingtoref/chr{chr}.missingtoref.pass.snp.phased.ag1000g.strict.bed",
            chr=CHRS),
        # LD_prune:
        expand("data/chr{chr}.pass.snp.flt.prune.out", chr=CHRS),
        expand("data/chr{chr}.pass.snp.flt.prune.out.tab", chr=CHRS),
        # combine_LD_prune:
        "data/all.pass.snp.flt.prune.out",
        "data/all.pass.snp.flt.prune.out.tab",
        # mds_all:
        "data/all.pass.snp.flt.mds",
        # mds_all_no_inv:
        "data/all.pass.snp.flt.noinv.mds",
        # mds_chr:
        expand("data/chr{chr}.pass.snp.flt.mds", chr=CHRS),
        expand("data/chr{chr}.pass.snp.flt.mds.eigvals", chr=CHRS),
        # plot_mds_all and plot_mds_all_no_inv:
        expand("reports/all.pass.snp.flt{inv_suffix}.ibs_mds_plot.pdf",
            inv_suffix=['', '.noinv']),
        expand("reports/all.pass.snp.flt{inv_suffix}.mds.txt",
            inv_suffix=['', '.noinv']),
        # plot_mds_chr:
        expand("reports/chr{chr}.pass.snp.flt.ibs_mds_plot.pdf", chr=CHRS),
        expand("reports/chr{chr}.pass.snp.flt.mds.txt", chr=CHRS),
        # pca_all_no_inv:
        "data/chr3.pass.snp.flt.noinv.eigenvec",
        "data/chr3.pass.snp.flt.noinv.eigenval",
        # pca_all_no_inv_rare:
        "data/rare_vars.chr3.pass.snp.flt.noinv.eigenvec",
        "data/rare_vars.chr3.pass.snp.flt.noinv.eigenval",
        # plot_pca:
        "reports/chr3.pass.snp.flt.noinv.pca_plot.pdf",
        "reports/rare_vars.chr3.pass.snp.flt.noinv.pca_plot.pdf",
        # plot_pca_bugala:
        "reports/chr3.pass.snp.flt.noinv.pca_plot.bugala.pdf",
        "data/ssese.samples.is-BUGALAML.txt",
        "data/ssese.samples.is-BUGALAIS.txt",
        "data/ssese.seqids.is-BUGALAML.txt",
        "data/ssese.seqids.is-BUGALAIS.txt",
        "data/ssese_individual_info_bugala_split.csv",
        "data/ssese_individual_info_simple_bugala_split.txt",
        # compute_fst_site:
        expand("results/chr{chr}.{site_pair}.weir.fst", chr=CHRS, site_pair=SITE_PAIRS),
        # compute_fst_win:
        expand("results/chr{chr}.{site_pair}.windowed.weir.fst", chr=CHRS, site_pair=SITE_PAIRS),
        # convert_to_eigenstrat:
        "data/all.pass.snp.flt.eigen.eigenstratgeno",
        # compute_fst_hudson:
        "data/all.pass.snp.flt.eigen.fst.se.out",
        # parse_fst_hudson:
        "data/all.pass.snp.flt.eigen.fst.se.out.fst.txt",
        "data/all.pass.snp.flt.eigen.fst.se.out.sd.txt",
        "data/all.pass.snp.flt.eigen.fst.se.out.fstZ.txt",
        # make_fst_table:
        "reports/hudson_fst.tex",
        "reports/fst_heatmap_hudson.pdf",
        # avg_fst:
        "reports/fst_summary.txt",
        # fst_heatmap:
        "reports/fst_heatmap_no_het_or_inv.pdf",
        "reports/fst_heatmap_no_het_2La.pdf",
        "reports/fst_heatmap_no_het_2Rb.pdf",
        "reports/fst_distributions.pdf",
        # bootstrap_for_fst_pval:
        expand("results/bootstrapped_Fst/3L.{site_pair}.real.log",
            site_pair=SITE_PAIRS),
        expand("results/bootstrapped_Fst/3L.{site_pair}.weighted.fst",
            site_pair=SITE_PAIRS),
        #   # compute_fst_pval:
        #   "reports/Fst_pvals.txt",
        # do_mantel:
        "reports/mantel_results.txt",
        "reports/fst_vs_distance.pdf",
        # compute_pi_window:
        expand("results/chr{chr}.{site}.windowed.pi", chr=CHRS, site=ISLANDS_SPLIT_BUGALA),
        # plot_pi_by_site:
        "reports/pi_boxplot.pdf",
        "reports/pi_boxplot.mainland-island.pdf",
        #   # plot_fst_by_pi:
        #   "reports/fst_by_pi.pdf",
        # compute_tajimad:
        expand("results/chr{chr}.{site}.Tajima.D", chr=CHRS, site=ISLANDS_SPLIT_BUGALA),
        # plot_tajimas_d_by_site:
        "reports/tajimas_d_boxplot.pdf",
        "reports/tajimas_d_boxplot.mainland-island.pdf",
        # plot_tajimas_d_across_genome:
        "reports/tajimas_D_by_site.pdf",
        # compute_inbreeding:
        expand("results/all.{site}.het", chr=CHRS, site=ISLANDS_SPLIT_BUGALA),
        # compute_allele_freq:
        expand("results/chr{chr}.{site}.frq", chr=CHRS, site=ISLANDS_SPLIT_BUGALA),
        # sample_allele_freq:
        expand("results/chr{chr}.{site}.sample.frq", chr=CHRS, site=ISLANDS_SPLIT_BUGALA),
        # plot_maf:
        "reports/MAF_boxplot.pdf",
        "reports/MAF_boxplot.mainland-island.pdf",
        "reports/MAF_histogram_SFS.pdf",
        # phase:
        expand("data/chr{chr}.pass.snp.phased.haps", chr=CHRS),
        # find_inv_het_snps_in_haps_files:
        expand("data/chr{chr}.pass.snp.phased.haps.excluded.snps", chr=CHRS),
        # remove_inv_het_snps_from_haps_files:
        expand("data/nohetinv/chr{chr}.pass.snp.phased.nohetinv.haps", chr=CHRS),
        expand("data/nohetinv/chr{chr}.pass.snp.phased.nohetinv.sample", chr=CHRS),
        # make_pop_haps_files:
        expand("data/chr{chr}.pass.snp.phased.{site}.haps", chr=CHRS, site=ISLANDS_SPLIT_BUGALA),
        #   # make_pop_haps_files_isl_ml:
        #   expand("data/chr{chr}.pass.snp.phased.{grouping}.haps", chr=CHRS, grouping=GROUPINGS),
        # make_phased_vcf:
        expand("data/chr{chr}.pass.snp.phased.haps.vcf", chr=CHRS),
        # sample_for_r2:
        expand("results/chr{chr}.passing_SNPs_for_r2.full.recode.vcf", chr=['3L','3R']),
        expand("results/chr{chr}.passing_SNPs_for_r2.sample.txt", chr=['3L','3R']),
        # compute_r2:
        expand("results/chr{chr}.{site}.list.hap.ld", chr=['3L','3R'], site=ISLANDS_SPLIT_BUGALA),
        # plot_r2:
        "reports/r2_boxplot.pdf",
        "reports/r2_boxplot.mainland-island.pdf",
        "reports/r2_decay.pdf",
        # prep_for_jp_xpehh:
        expand("data/xpehh_input_jp/chr{chr}.pass.snp.phased.{site}.jp.haps.transpose",
            chr=CHRS, site=ISLANDS_SPLIT_BUGALA),
        expand("data/xpehh_input_jp/chr{chr}.pass.snp.phased.jp.map",
            chr=CHRS),
        # compute_xpehh:
        expand("results/selscan/xp-ehh.{site_pair}.{chr}.xpehh.out", chr=CHRS, site_pair=map(lambda x: x.replace("_", "-"), SITE_PAIRS)),
        #   # compute_xpehh_jp:
        #   expand("results/xpehh_jp/xp-ehh.{site_pair}.{chr}.xpehh.out", chr=CHRS, site_pair=map(lambda x: x.replace("_", "-"), SITE_PAIRS)),
        # norm_xpehh:
        expand("results/selscan/xp-ehh.{chr}.norm.sentinel.txt", chr=CHRS),
        expand("results/selscan/xp-ehh.{site_pair}.{chr}.xpehh.out.norm", chr=CHRS, site_pair=map(lambda x: x.replace("_", "-"), SITE_PAIRS)),
        # avg_xpehh:
        expand("results/selscan/xp-ehh.{site_pair}.{chr}.xpehh.out.norm.avg", chr=CHRS, site_pair=map(lambda x: x.replace("_", "-"), SITE_PAIRS)),
        # plot_xpehh_fst:
        "reports/fst_xpehh_by_site.pdf",
        # plot_xpehh:
        expand("reports/selscan/xp-ehh.{site_pair}.{chr}.xpehh.out.norm.pdf", chr=CHRS, site_pair=map(lambda x: x.replace("_", "-"), SITE_PAIRS)),
        # rasterize_pdfs:
        expand("reports/xpehh_and_fst_between_populations_{chr}.all.png", chr=CHRS),
        # find_xpehh_peaks:
        "results/xpehh_peaks.txt",
        # compute_ihs:
        expand("results/selscan/iHS.{site}.{chr}.ihs.out", chr=CHRS, site=ISLANDS_SPLIT_BUGALA),
        # compute_css:
        "results/css.all.txt",
        # plot_css:
        expand("reports/css_{chr}.all.pdf", chr=CHRS),
        # norm_ihs:
        expand("results/selscan/iHS.{chr}.norm.sentinel.txt", chr=CHRS),
        expand("results/selscan/iHS.{site}.{chr}.ihs.out.100bins.norm",
            chr=CHRS, site=ISLANDS_SPLIT_BUGALA),
        # plot_ihs:
        expand("reports/ihs_{chr}.all.pdf", chr=CHRS),
        # make_h12_input_all:
        expand("data/chr{chr}.pass.snp.phased.forH12.txt", chr=CHRS),
        # make_h12_input_bysite:
        expand("data/chr{chr}.pass.snp.phased.{site}.forH12.txt",
            chr=CHRS, site=ISLANDS_SPLIT_BUGALA),
        # download_h12_script:
        "scripts/H12_H2H1.py",
        # compute_h12_all:
        expand("results/chr{chr}.pass.snp.phased.H12.out.txt", chr=CHRS),
        # compute_h12_bysite:
        expand("results/chr{chr}.pass.snp.phased.{site}.H12.out.txt",
            chr=CHRS, site=ISLANDS_SPLIT_BUGALA),
        # plot_h12_all:
        "reports/H12.pdf",
        # plot_h12_bysite:
        "reports/H12_by_site.pdf",
        # find_h12_outliers:
        "reports/island-specific-sweeps.tex",
        "reports/site-specific-sweeps.tex",
        "reports/insecticide-sweeps.tex",
        # prep_for_dapc:
        expand("results/dapc_subset_snps.{chr_subset}.raw",
            chr_subset=["chr3L", "chr3R"]),
        # do_dapc:
        expand("results/dapc_subset_snps.{chr_subset}.dapc.scatter.pdf",
            chr_subset=["chr3L", "chr3R"]),
        expand("results/dapc_subset_snps.{chr_subset}.clusters.by.site.pdf",
            chr_subset=["chr3L", "chr3R"]),
        expand("results/dapc_subset_snps.{chr_subset}.pca.colorplot.pdf",
            chr_subset=["chr3L", "chr3R"]),
        expand("results/dapc_subset_snps.{chr_subset}.BIC.pdf",
            chr_subset=["chr3L", "chr3R"]),
        # download_ag1000g:
        expand("data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.pass.biallelic.{chr}.vcf.gz", chr=CHRS),
        expand("data/accessibility/accessibility.{chr}.vcf.gz", chr=CHRS),
        # merge_ag1000g:
        expand("data/ssese_with_ag1000g/ssese_with_ag1000g.{chr}.flt.strict.vcf.gz",
             chr=CHRS),
        # merge_ag1000g_missingtoref:
        expand("data/ssese_with_ag1000g.missingtoref/ssese_with_ag1000g.{chr}.missingtoref.flt.strict.vcf.gz",
             chr=CHRS),
        # ss_ag_vcf_to_bed:
        expand("data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.bed",
            chr=CHR_THREE),
        # ss_ag_LD_prune:
        expand("data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.prune.out",
            chr=CHR_THREE),
        expand("data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.prune.out.tab",
            chr=CHR_THREE),
        # ss_ag_LD_prune_missingtoref:
        expand("data/ssese_with_ag1000g.missingtoref/chr{chr}.missingtoref.pass.snp.phased.ag1000g.strict.prune.out",
            chr=CHR_THREE),
        expand("data/ssese_with_ag1000g.missingtoref/chr{chr}.missingtoref.pass.snp.phased.ag1000g.strict.prune.out.tab",
            chr=CHR_THREE),
        # ss_ag_combine_LD_prune:
        expand("data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.prune.out"),
        expand("data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.prune.out.tab"),
        # ss_ag_combine_LD_prune_missingtoref:
        expand("data/ssese_with_ag1000g.missingtoref/all.missingtoref.pass.snp.phased.ag1000g.strict.prune.out"),
        expand("data/ssese_with_ag1000g.missingtoref/all.missingtoref.pass.snp.phased.ag1000g.strict.prune.out.tab"),
        # ss_ag_find_inds_to_exclude:
        "data/ag1000g.phase1.ar3/SY.individuals.txt",
        "data/ag1000g.phase1.ar3/AD.individuals.txt",
        "data/ag1000g.phase1.ar3/M.individuals.txt",
        "data/ag1000g.phase1.ar3/excluded.individuals.txt",
        "data/ag1000g.phase1.ar3/excluded.individuals.notM.txt",
        # ss_ag_mds_all:
        expand("data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.mds"),
        expand("data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.mds.eigvals"),
        # ss_ag_mds_all_no_inv:
        expand("data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.noinv.mds"),
        expand("data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.noinv.mds.eigvals"),
        # ss_ag_mds_chr:
        expand("data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.mds", chr=CHR_THREE),
        expand("data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.mds.eigvals", chr=CHR_THREE),
        # ss_ag_plot_mds_all and ss_ag_plot_mds_all_no_inv:
        "reports/all.pass.snp.phased.ag1000g.strict.ibs_mds_plot.pdf",
        "reports/all.pass.snp.phased.ag1000g.strict.noinv.ibs_mds_plot.pdf",
        # ss_ag_plot_mds_chr:
        expand("reports/chr{chr}.pass.snp.phased.ag1000g.strict.ibs_mds_plot.pdf", chr=CHR_THREE),
        # ss_ag_pca_all_no_inv:
        "data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.noinv.eigenvec",
        "data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.noinv.eigenval",
        # ss_ag_pca_all_no_inv_withM:
        "data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.noinv.withM.eigenvec",
        "data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.noinv.withM.eigenval",
        # ss_ag_pca_all_no_inv_withM_noKES:
        "data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.noinv.withMnoKES.eigenvec",
        "data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.noinv.withMnoKES.eigenval",
        # ss_ag_pca_all_no_inv_withM_missingtoref:
        "data/ssese_with_ag1000g.missingtoref/chr3.pass.snp.phased.ag1000g.strict.noinv.withM.eigenvec",
        "data/ssese_with_ag1000g.missingtoref/chr3.pass.snp.phased.ag1000g.strict.noinv.withM.eigenval",
        # ss_ag_plot_pca:
        "reports/chr3.pass.snp.phased.ag1000g.strict.noinv.pca_plot.1vs2.pca_plot.pdf",
        "reports/chr3.pass.snp.phased.ag1000g.strict.noinv.withM.pca_plot.1vs2.pca_plot.pdf",
        # ss_ag_plot_pca_missingtoref:
        "reports/chr3.missingtoref.pass.snp.phased.ag1000g.strict.noinv.withM.pca_plot.1vs2.pca_plot.pdf",
        # ss_ag_plot_pca_nokes:
        "reports/chr3.pass.snp.phased.ag1000g.strict.noinv.withMnoKES.pca_plot.1vs2.pca_plot.pdf",
        # make_dadi_input:
        expand("data/dadi/dadi.{mode}.data", mode=DADI_MODES),
        # polarize_dadi:
        expand("data/dadi_polarized/dadi.chr3.{mode}.data", mode=DADI_MODES),
        # prep_for_dadi_bootstrap:
        "data/random_hunks_for_dadi_bootstrap.bed",
        "data/dadi_bootstraps/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa",
        # bootstrap_for_dadi:
        # Replaced with sentinel file: expand("data/dadi_bootstraps/bs_{iter}.{mode}.data", iter=range(0,1000), mode=DADI_MODES),
        expand("data/dadi_bootstraps/bs_X.{mode}.data.marker",
            mode=DADI_MODES),
        # run_dadi:
        expand("results/dadi_polarized/dadi.{mode_subset}.iter{dadi_iter}.out",
            mode_subset=DADI_OUTPUT,
            dadi_iter=range(0,9)),
        expand("data/dadi_polarized/dadi.chr3.{mode_subset}.iter{dadi_iter}.sfs.txt",
            mode_subset=DADI_OUTPUT_PERIOD,
            dadi_iter=range(0,9)),
        # make_dadi_input_ag1000g:
        "data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.pass.biallelic.all.for.dadi.vcf",
        "data/dadi/dadi_pop_list.ag1000g.txt",
        "data/dadi/dadi.ag1000g.data",
        # polarize_dadi_ag1000g:
        "data/dadi_polarized/dadi.ag1000g.data",
        # bootstrap_for_dadi_ag:
        # Replaced with sentinel file: expand("data/dadi_bootstraps/bs_{iter}.ag1000g.data", iter=range(0,1000)),
        "data/dadi_bootstraps/bs_X.ag1000g.data.marker",
        # run_dadi_2pop_intersite:
        expand("results/dadi_polarized/dadi.island.2pop.{site_pair}.iter{dadi_iter}.out",
            site_pair=DADI_SITE_PAIRS,
            dadi_iter=range(0,9)),
        # fix_dadi_results:
        expand("results/dadi_polarized/dadi.{mode_subset}.iter{dadi_iter}.out.fix",
            mode_subset=DADI_OUTPUT,
            dadi_iter=range(0,9)),
        expand("results/dadi_polarized/dadi.island.2pop.{site_pair}.iter{dadi_iter}.out.fix",
            site_pair=DADI_SITE_PAIRS,
            dadi_iter=range(0,9)),
        # parse_dadi:
        "reports/dadi_three_epoch_pop_history.pdf",
        "reports/dadi.island.out.tex",
        "reports/dadi_island_is-ml_differences.txt",
        "reports/dadi_2pop_is-ml_differences.txt",
        "reports/dadi.2pop.island.island.out.tex",
        "reports/dadi.2pop.island.mainland.out.tex",
        "reports/dadi.2pop.mainland.mainland.out.tex",
        "reports/dadi_migration_matrix.txt",
        "reports/dadi_top_migration.txt",
        "reports/dadi_migration_matrix.pdf",
        "reports/dadi_2pop_is-ml_differences.time-size.txt",
        "results/dadi.best_iterations.txt",
        "results/dadi.best_iterations.2pop.txt",
        "results/asymmetrical.migration.txt",
        # copy_dadi_optimization_plots:
        expand("data/dadi_polarized/dadi.chr3.island.{site}.bestmodel.best-iter.pdf",
            site=ISLANDS_SPLIT_BUGALA),
        expand("data/dadi_polarized/dadi.chr3.island.{site_pair}.best-iter.pdf",
            site_pair=map(lambda x: x.replace("_", "-"), SITE_PAIRS)),
        # explore_theta:
        expand("results/dadi.chr3.island.{site}.theta_along_chr3.txt",
            site=ISLANDS_SPLIT_BUGALA),
        # create_thinned_adm_beds:
        expand("data/ssese_with_ag1000g/adm_subsets/chr3.{M_suffix}replicate{rep}.LD.bed",
            M_suffix=['', 'withM.'], rep=range(1,11)),
        # create_thinned_adm_beds_missingtoref:
        expand("data/ssese_with_ag1000g.missingtoref/adm_subsets/chr3.{M_suffix}replicate{rep}.LD.bed",
            M_suffix=['', 'withM.'], rep=range(1,11)),
        # adm_cv
        expand("data/ssese_with_ag1000g/adm_subsets/chr3{M_part}replicate{rep}.LD.iter{iter}.{k}.ADMIXTURE_log{k}.out",
            M_part=['.withM.'], rep=range(1,11), k=range(2,11), iter=range(1,6)),
        expand("data/ssese_with_ag1000g/adm_subsets/chr3{M_part}replicate{rep}.LD.iter{iter}.{k}.P",
            M_part=['.withM.'], rep=range(1,11), k=range(2,11), iter=range(1,6)),
        expand("data/ssese_with_ag1000g/adm_subsets/chr3{M_part}replicate{rep}.LD.iter{iter}.{k}.Q",
            M_part=['.withM.'], rep=range(1,11), k=range(2,11), iter=range(1,6)),
        # plot_adm_cv:
        expand("data/ssese_with_ag1000g/adm_subsets/chr3{M_part}CV.pdf",
            M_part=['.withM.']),
        # make_clumpak_input:
        "data/ssese_with_ag1000g/adm_subsets/chr3.withM.allQs.zip",
        # make_clumpak_pop_file:
        "data/ssese_with_ag1000g/adm_subsets/pop.list.withM.txt",
        # make_bestk_input:
        "data/ssese_with_ag1000g/adm_subsets/bestK.input.withM.log_prob.txt",
        # warn_abt_clumpak_results:
        expand("data/ssese_with_ag1000g/adm_subsets/CLUMPAK_{m_flag}/K{k}/MajorCluster/CLUMPP.files/ClumppIndFile.output",
            m_flag=['withM'], k=range(2,11)),
        # parse_clumpak_results:
        expand("data/ssese_with_ag1000g/adm_subsets/CLUMPAK_{m_flag}.{k}.Q",
            m_flag=['withM'], k=range(2,11)),
        # adm and adm_M:
        expand("data/ssese_with_ag1000g/chr3.withM.pass.snp.phased.ag1000g.strict.{k}.P",
            k=range(2,11)),
        expand("data/ssese_with_ag1000g/chr3.withM.pass.snp.phased.ag1000g.strict.{k}.Q",
            k=range(2,11)),
        # plot_adm_M:
        expand("data/ssese_with_ag1000g/chr3.withM.pass.snp.phased.ag1000g.strict.{k}.pdf",
            k=range(2,11)),
        expand("data/ssese_with_ag1000g/chr3.withM.pass.snp.phased.ag1000g.strict.LVB_only.{k}.pdf",
            k=range(2,11)),
        # find_ROHs:
        "results/ROH/all.pass.snp.flt.noinv.hom",
        "results/ROH/all.pass.snp.flt.noinv.hom.indiv",
        "results/ROH/all.pass.snp.flt.noinv.hom.overlap",
        "results/ROH/all.pass.snp.flt.noinv.hom.summary",
        # plot_ROHs:
        "reports/roh_plot.pdf",
        "reports/roh_boxplot.pdf",
        # identify_IBD:
        "results/IBD/all.pass.snp.flt.noinv.genome",
        # plot_IBD:
        "reports/ibd.pdf",
        # plot_EHH_ROI:
        "results/ehh/chrX.pass.snp.phased.figs/ehh_chrX.X_9Mb_1.pdf",
        "results/ehh/chrX.pass.snp.phased.figs/ehh_chrX.X_4Mb_1.pdf",
        "results/ehh/chr2L.pass.snp.phased.figs/ehh_chr2L.2L_34Mb_1.pdf",
        "results/ehh/chr2R.pass.snp.phased.figs/ehh_chr2R.CYP6P2_1.pdf",
        "results/ehh/chrX.pass.snp.phased.figs/ehh_chrX.CYP9K1_1.pdf",
        # annotate_snps:
        expand("data/chr{chr}.pass.snp.eff.vcf", chr=CHRS),
        # infer_dendrograms_roi:
        expand("results/for_dendrograms/{locus}.{window}.{ploidy}.with-ag1000g.sans-snpeff.haps.{ending}",
            locus = ['chr2L.2L_34Mb', 'chrX.X_9Mb', 'chrX.X_4Mb', 'chrX.X_CYP9K1', 'chr2R.2R_CYP6P2', 'chr3R.3R_GSTE'],
            window = ['10000', '1e+05'],      # Cut '1e+06'
            ploidy = ['haploid'],             # Cut 'diploid'
            ending = ['fancydendro.pdf', 'haplotypes.txt']),
        # plot_map_o_haps:
        "results/ssese_sampling_sites.png",
        "results/ssese_haplotypes_sweep_2L_34Mb.png",
        "results/ssese_haplotypes_sweep_X_9Mb.png",
        expand("data/ssese_with_ag1000g/adm_subsets/CLUMPAK_withM.{k}.Q.map.png",
            k=range(2,11)),
        # make_zoomed_stats_plot:
        expand("reports/all_stats.{locus}.{site}.pdf",
            locus = ['2L_34Mb', 'X_9Mb', 'X_4Mb', 'CYP6P2', 'CYP9K1', 'GSTE'],
            site=ISLANDS_SPLIT_BUGALA),
        # compute_coverage:
        "reports/all.pass.snp.flt.idepth",
        # stairway_plot_site:
        expand("data/stairway_plot_island.{site}/island.{site}.final.summary.pdf",
            site=ISLANDS_SPLIT_BUGALA),
        # plot_stairway_plots:
        "reports/stairway.plot.superimposed.pdf",
        # estimate_Ne:
        "results/genepop_Ne_estimates.txt",
        "results/genepop_bounds.txt",
        # run_ibdseq:
        expand("results/chr{chr}.pass.snp.flt.ibdseq.ibd",
            chr=['3L', '3R']),
        # run_ibdne:
        expand("results/chr{chr}.pass.snp.flt.{site}.ne.ne",
            chr=['3L', '3R'], site=ISLANDS_SPLIT_BUGALA),
        # plot_ibdne:
        expand("reports/IBDNe_by_pop.pdf",
            chr=['3L', '3R']),
        expand("reports/IBDNe_by_site_type.pdf",
            chr=['3L', '3R']),
        expand("reports/IBDNe_by_site_type_Mityana_highlighted.pdf",
            chr=['3L', '3R']),
        # ibs_demog_prep:
        expand("results/IBS.{site}.popdata",
            site=ISLANDS_SPLIT_BUGALA),
        # ibs_demog_1pop:
        expand("results/IBS_inferred_size_{site}.demographic_history.txt",
            site=ISLANDS_SPLIT_BUGALA),
        expand("results/IBS_inferred_size_{site}.demographic_history.txt.pdf",
            site=ISLANDS_SPLIT_BUGALA),
        # ibs_demog_2pop:
        expand("results/IBS_inferred_size_{site_pair_vs}.demographic_history.txt", site_pair_vs=map(lambda x: x.replace("_", "_vs_"), SITE_PAIRS)),
        # count_snps:
        "reports/SNP_count.txt",
        # make_ind_table:
        "reports/individual_table.tex",
        "reports/site_gps_table.tex"

localrules: all

# ----------------------------------------------------------------------------------------
# --- Make lists of individuals
# ----------------------------------------------------------------------------------------

rule make_island_lists:
    input:
        "data/ssese_individual_info.csv"
    output:
        "data/ssese.samples.is-{island}.txt"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "awk 'BEGIN {{FS=\",\"}} {{ if ($9 ~ \"AG\" && $1 ~ \"{wildcards.island}\" && $3 != \"0\" ) print $3 }}' {input} > {output}"

rule make_site_lists:
    input:
        "data/ssese_individual_info.csv"
    output:
        "data/ssese.samples.site-{pop}.txt"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "ISLAND=`echo {wildcards.pop} | cut -d'-' -f1`;"
        "SITE=`echo {wildcards.pop} | cut -d'-' -f2`;"
        "awk -v island=\"$ISLAND\" -v site=\"$SITE\" 'BEGIN {{FS=\",\"}} {{ if ($9 ~ \"AG\" && $1 ~ island && $2 ~ site && $3 != \"0\" ) print $3 }}' {input} > {output}"

# ----------------------------------------------------------------------------------------
# --- Map sample IDs to sequencing IDs
# ----------------------------------------------------------------------------------------

rule map_sample_to_seq_id:
    input:
        "data/ssese_plate1_individuals.csv",
        "data/ssese_plate2_individuals.csv",
    output:
        "data/sample_to_seq_id_mapping.txt"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/map_sample_to_seq_id.R"

# ----------------------------------------------------------------------------------------
# --- Convert from mosquito ID lists to sequencing ID lists
# ----------------------------------------------------------------------------------------

rule make_seqid_lists:
    input:
        "data/sample_to_seq_id_mapping.txt",
        sample_list="data/ssese.samples.{group}.txt"
    output:
        "data/ssese.seqids.{group}.txt"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/convert_to_seq_id.R {input.sample_list}"

# ----------------------------------------------------------------------------------------
# --- Make simple info list
# ----------------------------------------------------------------------------------------

rule make_simple_info_list:
    input:
        "data/sample_to_seq_id_mapping.txt",
        "data/ssese_individual_info.csv"
    output:
        "data/ssese_individual_info_simple.txt"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/make_simple_list.R"

# ========================================================================================
# --- Make sets and combined BED for inversions and heterochromatin regions
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Combine BEDs of inversions and heterochromatic regions
# ----------------------------------------------------------------------------------------

rule combine_het_inv:
    input:
        inv="data/inversion_simple.bed",
        het="data/heterochromatin.bed"
    output:
        "data/inversion_heterochromatin.bed"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "cat {input.inv} {input.het} > {output};"

# ----------------------------------------------------------------------------------------
# --- Make sets for inversions and heterochromatic regions
# ----------------------------------------------------------------------------------------

rule make_het_inv_sets:
    input:
        inv="data/inversion_simple.bed",
        het="data/heterochromatin.bed"
    output:
        out_inv="data/inversion_simple.set",
        out_het="data/heterochromatin.set"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "cat {input.inv} | sed -e \"s/2L/1/\" -e \"s/2R/2/\" -e \"s/3L/3/\" -e \"s/3R/4/\" "
        "    > {output.out_inv};"
        "cat {input.het} | sed -e \"s/2L/1/\" -e \"s/2R/2/\" -e \"s/3L/3/\" -e \"s/3R/4/\" | "
        "    awk 'BEGIN {{ OFS = \"\t\" }}{{ print $1,$2,$3,HET }}' > {output.out_het};"

# ========================================================================================
# --- LD prune steps
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Combine all chromosomal VCF files
# ----------------------------------------------------------------------------------------

rule combine_chr_vcfs:
    input:
        VCF_GZ
    output:
        "data/all.pass.snp.vcf.gz"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load vcftools;"
        "vcf-concat {input} | bgzip -c > {output}"

# ----------------------------------------------------------------------------------------
# --- Compute individual missingness (overall genome)
# ----------------------------------------------------------------------------------------

rule compute_indiv_missingness:
    input:
        "data/all.pass.snp.vcf.gz"
    output:
        "reports/all.pass.snp.imiss"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load vcftools;"
        "vcftools --gzvcf {input} --missing-indv --out reports/all.pass.snp"

# ----------------------------------------------------------------------------------------
# --- Find individuals with high missingness
# ----------------------------------------------------------------------------------------

rule find_high_missing_indivs:
    input:
        "reports/all.pass.snp.imiss"
    output:
        "reports/inds.high.missing.txt"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sh scripts/find_high_missing_inds.sh"

# ----------------------------------------------------------------------------------------
# --- Filter by missingness and remove invariant sites
# ----------------------------------------------------------------------------------------

# --max-missing 0.9 means no more than 10% missingness allowed

rule filter_vcf:
    input:
        invcf="data/chr{chr}.pass.snp.vcf",
        missing="reports/inds.high.missing.txt"
    output:
        "data/chr{chr}.pass.snp.flt.vcf.gz"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load vcftools;"
        "vcftools --vcf {input.invcf} --min-alleles 2 --max-alleles 2 "
        "--max-missing 0.9 --hwe 0.00001 --remove {input.missing} --recode --stdout | "
        "bgzip -c > {output}"

# -------------------------------------------------------------------------------------- #
# --- Convert *FILTERED* VCF to plink's BED format (via temporary PED)
# -------------------------------------------------------------------------------------- #

rule vcf_to_bed:
    input:
        invcf="data/chr{chr}.pass.snp.flt.vcf.gz",
    output:
        "data/chr{chr}.pass.snp.flt.bed"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load plink;"
        "module load vcftools;"
        "OUT_FILE={output};"
        "mkdir -p tmp;"
        "TMP=`pwd`/tmp;"
        "vcftools --gzvcf {input.invcf} --plink --out ${{OUT_FILE/.bed}} --temp $TMP;"
        # Edit the MAP file and change 0 to the right chromosome index:
        #   1   2L
        #   2   2R
        #   3   3L
        #   4   3R
        #   5   X
        "sh scripts/fix_chr_in_ped.sh ${{OUT_FILE/.bed}}.map;"
        "plink --file ${{OUT_FILE/.bed}} --make-bed --out ${{OUT_FILE/.bed}};"
        "rm ${{OUT_FILE/.bed}}.ped"

# -------------------------------------------------------------------------------------- #
# --- Find SNPs in LD to prune out
# -------------------------------------------------------------------------------------- #

# Slide a sliding window of 50 SNPs across each chromosome and calculate r^2 for each
# pair in the window. Remove one random SNP in each pair with r^2 > 0.5.
# Parameters are window size, step, and r^2 threshold.

# To use this list of pruned SNPs to create an LD-pruned dataset, use:
#   plink --exclude data/chrX.pass.snp.flt.prune.out
#   vcftools --exclude-positions data/chrX.pass.snp.flt.prune.out.tab

rule LD_prune:
    input:
        "data/chr{chr}.pass.snp.flt.bed"
    output:
        "data/chr{chr}.pass.snp.flt.prune.out",
        "data/chr{chr}.pass.snp.flt.prune.out.tab"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load plink;"
        "PREFIX=`echo {input} | sed -e \"s/\.bed//\"`;"
        "plink --bfile $PREFIX --indep-pairwise 50 5 0.5 --out $PREFIX;"
        "sed -e \"s/:/\\t/\" {output[0]} > {output[1]};"
        "mkdir -p reports/plink_LD_pruning_part1/;"
        "mv $PREFIX.log reports/plink_LD_pruning_part1/"

# -------------------------------------------------------------------------------------- #
# --- Combine lists of SNPs in LD to prune out
# -------------------------------------------------------------------------------------- #

# To use this list of pruned SNPs to create an LD-pruned dataset, use:
#   plink --exclude data/all.pass.snp.flt.prune.out
#   vcftools --exclude-positions data/all.pass.snp.flt.prune.out.tab

rule combine_LD_prune:
    input:
        expand("data/chr{chr}.pass.snp.flt.prune.out", chr=CHRS),
        expand("data/chr{chr}.pass.snp.flt.prune.out.tab", chr=CHRS)
    output:
        "data/all.pass.snp.flt.prune.out",
        "data/all.pass.snp.flt.prune.out.tab"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "cat data/chr*.pass.snp.flt.prune.out > data/all.pass.snp.flt.prune.out;"
        "cat data/chr*.pass.snp.flt.prune.out.tab > data/all.pass.snp.flt.prune.out.tab;"

# ========================================================================================
# --- MDS
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Compute MDS from IBS values - full genome, LD pruned
# ----------------------------------------------------------------------------------------

rule mds_all:
    input:
        expand("data/chr{chr}.pass.snp.flt.bed", chr=CHRS),
        "data/all.pass.snp.flt.prune.out"
    output:
        "data/all.pass.snp.flt.mds",
        "data/all.pass.snp.flt.mds.eigvals"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load plink;"
        # Temporary files and prefix (needed to not conflict with other jobs when writing
        # temporary BED, BIM, and FAM files)
        "TMP_BED_LIST=`mktemp tmp.chr_bed_list.XXXXX`;"
        "TMP_PREFIX=`mktemp -u data/tmp.all.pass.snp.flt.XXXXX`;"
        # Create list of input BED files (by chr)
        "ls data/chr*.pass.snp.flt.bed | sed -e \"s/\(.*\).bed/\\1.bed \\1.bim \\1.fam/\" > $TMP_BED_LIST;"
        # Compute MDS
        "plink --merge-list $TMP_BED_LIST --exclude data/all.pass.snp.flt.prune.out --cluster --mds-plot 2 eigvals --out $TMP_PREFIX;"
        # Move results
        "mv $TMP_PREFIX.mds data/all.pass.snp.flt.mds;"
        "mv $TMP_PREFIX.mds.eigvals data/all.pass.snp.flt.mds.eigvals;"
        # Clean up
        "rm $TMP_BED_LIST;"
        "rm $TMP_PREFIX.*"

# ----------------------------------------------------------------------------------------
# --- Compute MDS from IBS values - full genome w/out inversions, LD pruned
# ----------------------------------------------------------------------------------------

rule mds_all_no_inv:
    input:
        expand("data/chr{chr}.pass.snp.flt.bed", chr=CHRS),
        "data/all.pass.snp.flt.prune.out",
        "data/inversion_simple.bed"
    output:
        "data/all.pass.snp.flt.noinv.mds",
        "data/all.pass.snp.flt.noinv.mds.eigvals"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load plink;"
        # Temporary files and prefix (needed to not conflict with other jobs when writing
        # temporary BED, BIM, and FAM files)
        "TMP_BED_LIST=`mktemp tmp.chr_bed_list.XXXXX`;"
        "TMP_INV_SIMP=`mktemp tmp.inversion_simple.XXXXX`;"
        "TMP_EXCLUDE=`mktemp tmp.exclude.XXXXX`;"
        "TMP_PREFIX=`mktemp -u data/tmp.all.pass.snp.flt.noinv.XXXXX`;"
        # Create list of input BED files (by chr) and simple inversions BED
        "ls data/chr*.pass.snp.flt.bed | sed -e \"s/\(.*\).bed/\\1.bed \\1.bim \\1.fam/\" > $TMP_BED_LIST;"
        "sed -e \"s/^2L/1/\" -e \"s/^2R/2/\" data/inversion_simple.bed > $TMP_INV_SIMP;"
        # Make SNP set covering the 2La and 2Rb inversions
        "plink --merge-list $TMP_BED_LIST --exclude data/all.pass.snp.flt.prune.out --make-set $TMP_INV_SIMP --write-set --out data/all.pass.snp.flt.mds.inversions;"
        # Compute MDS
        "cat data/all.pass.snp.flt.prune.out data/all.pass.snp.flt.mds.inversions.set > $TMP_EXCLUDE;"
        "plink --merge-list $TMP_BED_LIST --exclude $TMP_EXCLUDE --cluster --mds-plot 2 eigvals --out $TMP_PREFIX --write-snplist;"
        # Move result
        "mv $TMP_PREFIX.mds data/all.pass.snp.flt.noinv.mds;"
        "mv $TMP_PREFIX.mds.eigvals data/all.pass.snp.flt.noinv.mds.eigvals;"
        # Clean up
        "rm $TMP_BED_LIST;"
        "rm $TMP_EXCLUDE;"
        "rm $TMP_INV_SIMP;"
        "rm data/all.pass.snp.flt.mds.inversions.*;"
        "rm $TMP_PREFIX.*;"

# ----------------------------------------------------------------------------------------
# --- Compute MDS from IBS values - by chromosome
# ----------------------------------------------------------------------------------------

rule mds_chr:
    input:
        bed="data/chr{chr}.pass.snp.flt.bed",
        ld_set="data/chr{chr}.pass.snp.flt.prune.out"
    output:
        "data/chr{chr}.pass.snp.flt.mds",
        "data/chr{chr}.pass.snp.flt.mds.eigvals"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load plink;"
        "plink --bfile data/chr{wildcards.chr}.pass.snp.flt --exclude {input.ld_set} --cluster --mds-plot 2 eigvals --out data/chr{wildcards.chr}.pass.snp.flt;"

# ----------------------------------------------------------------------------------------
# --- Plot MDS
# ----------------------------------------------------------------------------------------

rule plot_mds_all:
    input:
        "data/sample_to_seq_id_mapping.txt",
        "data/all.pass.snp.flt.mds",
    output:
        "reports/all.pass.snp.flt.ibs_mds_plot.pdf",
        "reports/all.pass.snp.flt.mds.txt"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_ibs_mds_or_pca.R data/all.pass.snp.flt.mds"

rule plot_mds_all_no_inv:
    input:
        "data/sample_to_seq_id_mapping.txt",
        "data/all.pass.snp.flt.noinv.mds",
    output:
        "reports/all.pass.snp.flt.noinv.ibs_mds_plot.pdf",
        "reports/all.pass.snp.flt.noinv.mds.txt"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_ibs_mds_or_pca.R data/all.pass.snp.flt.noinv.mds"

rule plot_mds_chr:
    input:
        "data/sample_to_seq_id_mapping.txt",
        "data/chr{chr}.pass.snp.flt.mds"
    output:
        "reports/chr{chr}.pass.snp.flt.ibs_mds_plot.pdf",
        "reports/chr{chr}.pass.snp.flt.mds.txt"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_ibs_mds_or_pca.R data/chr{wildcards.chr}.pass.snp.flt.mds"

# ========================================================================================
# --- PCA
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Do PCA - chr3 only, LD pruned, no het or inv - and split Bugala individuals
# ----------------------------------------------------------------------------------------

rule pca_all_no_inv:
    input:
        expand("data/chr{chr}.pass.snp.flt.bed", chr=["3L", "3R"]),
        "data/all.pass.snp.flt.prune.out",
        "data/inversion_simple.bed",
        "data/heterochromatin.bed"
    output:
        eigenvec="data/chr3.pass.snp.flt.noinv.eigenvec",
        eigenval="data/chr3.pass.snp.flt.noinv.eigenval",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        # Temporary files and prefix (needed to not conflict with other jobs when writing
        # temporary BED, BIM, and FAM files)
        "module load plink;"
        "TMP_BED_LIST=`mktemp tmp.chr_bed_list.XXXXX`;"
        "TMP_INV_HET_SIMP=`mktemp tmp.inversion_het_simple.XXXXX`;"
        "TMP_EXCLUDE=`mktemp tmp.exclude.XXXXX`;"
        "TMP_PREFIX=`mktemp -u data/tmp.all.pass.snp.flt.noinv.XXXXX`;"
        # Create list of input BED files (by chr) and simple inversions/het BED
        "ls data/chr*.pass.snp.flt.bed | sed -e \"s/\(.*\).bed/\\1.bed \\1.bim \\1.fam/\" > $TMP_BED_LIST;"
        "cat data/inversion_simple.bed > $TMP_INV_HET_SIMP;"
        "awk 'BEGIN {{OFS=\"\\t\"}} {{ print $1,$2,$3,\"het\" }}' data/heterochromatin.bed >> $TMP_INV_HET_SIMP;"
        "sed -e \"s/^2L/1/\" -e \"s/^2R/2/\" -e \"s/^3L/3/\" -e \"s/^3R/4/\" -iBACKUP $TMP_INV_HET_SIMP;"
        # Make SNP set covering the inversions and het regions
        "plink --merge-list $TMP_BED_LIST --exclude data/all.pass.snp.flt.prune.out --make-set $TMP_INV_HET_SIMP --write-set --out data/all.pass.snp.flt.pca.inversions;"
        # Compute PCA
        "cat data/all.pass.snp.flt.prune.out data/all.pass.snp.flt.pca.inversions.set > $TMP_EXCLUDE;"
        # Retain 100 PCs for scree plot
        "plink --merge-list $TMP_BED_LIST --exclude $TMP_EXCLUDE --maf 0.01 --cluster --pca 100 --out $TMP_PREFIX;"
        # Move result
        "mv $TMP_PREFIX.eigenvec {output.eigenvec}.forScree;"
        "mv $TMP_PREFIX.eigenval {output.eigenval}.forScree;"
        # Retain 18 PCs based on analysis of scree plot
        "plink --merge-list $TMP_BED_LIST --exclude $TMP_EXCLUDE --maf 0.01 --cluster --pca 18 --out $TMP_PREFIX;"
        # Move result
        "mv $TMP_PREFIX.eigenvec {output.eigenvec};"
        "mv $TMP_PREFIX.eigenval {output.eigenval};"
        # Clean up
        "rm $TMP_BED_LIST;"
        "rm $TMP_EXCLUDE;"
        "rm $TMP_INV_HET_SIMP;"
        "rm data/all.pass.snp.flt.pca.inversions.*;"
        "rm $TMP_PREFIX.*;"

# ----------------------------------------------------------------------------------------
# --- Do PCA - chr3 only, LD pruned, no het or inv - rare variants
# ----------------------------------------------------------------------------------------

rule pca_all_no_inv_rare:
    input:
        expand("data/chr{chr}.pass.snp.flt.bed", chr=["3L", "3R"]),
        "data/all.pass.snp.flt.prune.out",
        "data/inversion_simple.bed",
        "data/heterochromatin.bed"
    output:
        eigenvec="data/rare_vars.chr3.pass.snp.flt.noinv.eigenvec",
        eigenval="data/rare_vars.chr3.pass.snp.flt.noinv.eigenval",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        # Temporary files and prefix (needed to not conflict with other jobs when writing
        # temporary BED, BIM, and FAM files)
        "module load plink;"
        "TMP_BED_LIST=`mktemp tmp.chr_bed_list.XXXXX`;"
        "TMP_INV_HET_SIMP=`mktemp tmp.inversion_het_simple.XXXXX`;"
        "TMP_EXCLUDE=`mktemp tmp.exclude.XXXXX`;"
        "TMP_PREFIX=`mktemp -u data/tmp.all.pass.snp.flt.noinv.XXXXX`;"
        # Create list of input BED files (by chr) and simple inversions/het BED
        "ls data/chr*.pass.snp.flt.bed | sed -e \"s/\(.*\).bed/\\1.bed \\1.bim \\1.fam/\" > $TMP_BED_LIST;"
        "cat data/inversion_simple.bed > $TMP_INV_HET_SIMP;"
        "awk 'BEGIN {{OFS=\"\t\"}} {{ print $1,$2,$3,\"het\" }}' data/heterochromatin.bed >> $TMP_INV_HET_SIMP;"
        "sed -e \"s/^2L/1/\" -e \"s/^2R/2/\" -e \"s/^3L/3/\" -e \"s/^3R/4/\" -iBACKUP $TMP_INV_HET_SIMP;"
        # Make SNP set covering the inversions and het regions
        "plink --merge-list $TMP_BED_LIST --exclude data/all.pass.snp.flt.prune.out --make-set $TMP_INV_HET_SIMP --write-set --out data/rare_vars.pass.snp.flt.pca.inversions;"
        # Compute PCA
        "cat data/all.pass.snp.flt.prune.out data/rare_vars.pass.snp.flt.pca.inversions.set > $TMP_EXCLUDE;"
        "plink --merge-list $TMP_BED_LIST --exclude $TMP_EXCLUDE --maf 0.01 --max-maf 0.05 --cluster --pca 10 --out $TMP_PREFIX;"
        # Move result
        "mv $TMP_PREFIX.eigenvec {output.eigenvec};"
        "mv $TMP_PREFIX.eigenval {output.eigenval};"
        # Clean up
        "rm $TMP_BED_LIST;"
        "rm $TMP_EXCLUDE;"
        "rm $TMP_INV_HET_SIMP;"
        "rm data/rare_vars.pass.snp.flt.pca.inversions.*;"
        "rm $TMP_PREFIX.*;"

# ----------------------------------------------------------------------------------------
# --- Plot PCA
# ----------------------------------------------------------------------------------------

rule plot_pca:
    input:
        "data/sample_to_seq_id_mapping.txt",
        "data/{prefix}.eigenvec"
    output:
        "reports/{prefix}.pca_plot.pdf",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_ibs_mds_or_pca.R data/{wildcards.prefix}.eigenvec"

# ----------------------------------------------------------------------------------------
# --- Plot Bugala-focused PCA
# ----------------------------------------------------------------------------------------

rule plot_pca_bugala:
    input:
        "data/sample_to_seq_id_mapping.txt",
        "data/chr3.pass.snp.flt.noinv.eigenvec"
    output:
        "reports/chr3.pass.snp.flt.noinv.pca_plot.bugala.pdf",
        "data/ssese.samples.is-BUGALAML.txt",
        "data/ssese.samples.is-BUGALAIS.txt",
        "data/ssese.seqids.is-BUGALAML.txt",
        "data/ssese.seqids.is-BUGALAIS.txt",
        "data/ssese_individual_info_bugala_split.csv",
        "data/ssese_individual_info_simple_bugala_split.txt",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_ibs_mds_or_pca.R data/chr3.pass.snp.flt.noinv.eigenvec TRUE"

# ========================================================================================
# --- Differentiation
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Compute per-site Fst
# ----------------------------------------------------------------------------------------

rule compute_fst_site:
    input:
        invcf="data/chr{chr}.pass.snp.flt.vcf.gz",
        invhet="data/inversion_heterochromatin.bed"
    output:
        "results/chr{chr}.{site_pair}.weir.fst"
    threads: 1
    params: runtime="8",
            mem=",mem=5gb"
    shell:
        "module load vcftools;"
        "BOTH_POPS=`echo {wildcards.site_pair} | tr '_' ' '`;"
        "OUT_PREFIX=results/chr{wildcards.chr}.{wildcards.site_pair};"
        "sh scripts/compute_Fst.sh {input.invcf} $OUT_PREFIX $BOTH_POPS 1 0"

# ----------------------------------------------------------------------------------------
# --- Compute windowed Fst
# ----------------------------------------------------------------------------------------

rule compute_fst_window:
    input:
        invcf="data/chr{chr}.pass.snp.flt.vcf.gz",
        invhet="data/inversion_heterochromatin.bed"
    output:
        "results/chr{chr}.{site_pair}.windowed.weir.fst"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load vcftools;"
        "BOTH_POPS=`echo {wildcards.site_pair} | tr '_' ' '`;"
        "OUT_PREFIX=results/chr{wildcards.chr}.{wildcards.site_pair};"
        "sh scripts/compute_Fst.sh {input.invcf} $OUT_PREFIX $BOTH_POPS 0 1 10000"

# ----------------------------------------------------------------------------------------
# --- Compute Hudston Fst and its standard error (convert to correct format first)
# ----------------------------------------------------------------------------------------

rule convert_to_eigenstrat:
    input:
        expand("data/chr{chr}.pass.snp.flt.bed", chr=CHRS),
        "data/ssese_individual_info_simple_bugala_split.txt"
    output:
        "data/all.pass.snp.flt.eigen.eigenstratgeno"
    threads: 1
    params: runtime="24",
            mem=",mem=48gb"
    shell:
        "sh scripts/convert_to_eigenstrat.sh"

rule compute_fst_hudson:
    input:
        "data/all.pass.snp.flt.eigen.eigenstratgeno"
    output:
        "data/all.pass.snp.flt.eigen.fst.se.out"
    threads: 12
    params: runtime="24",
            mem=",mem=96gb"
    shell:
        "sh scripts/compute_fst_hudson.sh"

# ----------------------------------------------------------------------------------------
# --- Parse Hudston Fst and its standard error
# ----------------------------------------------------------------------------------------

rule parse_fst_hudson:
    input:
        "data/all.pass.snp.flt.eigen.fst.se.out"
    output:
        "data/all.pass.snp.flt.eigen.fst.se.out.fst.txt",
        "data/all.pass.snp.flt.eigen.fst.se.out.sd.txt",
        "data/all.pass.snp.flt.eigen.fst.se.out.fstZ.txt",
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/parse_hudson_fst_output.sh"

# ----------------------------------------------------------------------------------------
# --- Make table of Hudston Fst and its standard error
# ----------------------------------------------------------------------------------------

rule make_fst_table:
    input:
        "data/all.pass.snp.flt.eigen.fst.se.out.fst.txt",
        "data/all.pass.snp.flt.eigen.fst.se.out.sd.txt",
        "data/all.pass.snp.flt.eigen.fst.se.out.fstZ.txt",
    output:
        "reports/hudson_fst.tex",
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/make_fst_table.R"

# ----------------------------------------------------------------------------------------
# --- Average Fst
# ----------------------------------------------------------------------------------------

rule avg_fst:
    input:
        expand("results/chr{chr}.{site_pair}.weir.fst", chr=CHRS, site_pair=SITE_PAIRS),
        "data/inversion_sites.bed",
        "data/heterochromatin.bed"
    output:
        "reports/fst_summary.txt"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/average_fst.R > {output}"

# ----------------------------------------------------------------------------------------
# --- Plot Fst heatmap and distributions
# ----------------------------------------------------------------------------------------

rule fst_heatmap:
    input:
        "reports/fst_summary.txt",
        "data/inversion_sites.bed",
        "data/heterochromatin.bed"
    output:
        "reports/fst_heatmap_no_het_or_inv.pdf",
        "reports/fst_heatmap_no_het_2La.pdf",
        "reports/fst_heatmap_no_het_2Rb.pdf",
        "reports/fst_distributions.pdf",
        "data/inversion_sites.bed"
    threads: 1
    params: runtime="4",
            mem=",mem=96gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_fst_heatmap.R"

# ----------------------------------------------------------------------------------------
# --- Do bootstrapping to compute Fst p-value
# ----------------------------------------------------------------------------------------

rule bootstrap_for_fst_pval:
    input:
        vcf="data/chr3L.pass.snp.flt.vcf.gz",
        seq_ids=expand("data/ssese.seqids.is-{site}.txt", site=ISLANDS_SPLIT_BUGALA)
    output:
        "results/bootstrapped_Fst/3L.{site_pair}.real.log",
        "results/bootstrapped_Fst/3L.{site_pair}.weighted.fst"
    threads: 1
    params: runtime="12",
            mem=",mem=5gb"
    shell:
        "BOTH_POPS=\"`echo {wildcards.site_pair} | tr '_' ' '`\";"
        "sh scripts/compute_Fst_pvalue_via_bootstrapping.sh {input.vcf} $BOTH_POPS 1000 5000"

# ----------------------------------------------------------------------------------------
# --- Compute Fst p-value
# ----------------------------------------------------------------------------------------

rule compute_fst_pval:
    input:
        expand("results/bootstrapped_Fst/3L.{site_pair}.real.log", site_pair=SITE_PAIRS),
        expand("results/bootstrapped_Fst/3L.{site_pair}.weighted.fst", site_pair=SITE_PAIRS),
    output:
        "reports/Fst_pvals.txt",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/get_all_fst_pvals.R"

# ----------------------------------------------------------------------------------------
# --- Do Mantel test for correlation of Fst and geographic distance
# ----------------------------------------------------------------------------------------

rule do_mantel:
    input:
        "data/all.pass.snp.flt.eigen.fst.se.out.fst.txt",
        "data/site_gps_coordinates.txt",
    output:
        "reports/mantel_results.txt",
        "reports/fst_vs_distance.pdf"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/do_mantel_test_for_distance.R"

# ----------------------------------------------------------------------------------------
# --- Compute windowed pi
# ----------------------------------------------------------------------------------------

rule compute_pi_window:
    input:
        invcf="data/chr{chr}.pass.snp.flt.vcf.gz",
        invhet="data/inversion_heterochromatin.bed"
    output:
        "results/chr{chr}.{site}.windowed.pi"
    threads: 1
    params: runtime="8",
            mem=",mem=5gb"
    shell:
        "module load vcftools;"
        "OUT_PREFIX=results/chr{wildcards.chr}.{wildcards.site};"
        "vcftools --gzvcf {input.invcf} --window-pi 10000 "
        "    --exclude-bed {input.invhet} "
        "    --keep data/ssese.seqids.is-{wildcards.site}.txt "
        "    --out $OUT_PREFIX"

# ----------------------------------------------------------------------------------------
# --- Plot pi boxplots
# ----------------------------------------------------------------------------------------

rule plot_pi_by_site:
    input:
        expand("results/chr{chr}.{site}.windowed.pi", chr=['3L', '3R'], site=ISLANDS_SPLIT_BUGALA),
    output:
        "reports/pi_boxplot.pdf",
        "reports/pi_boxplot.mainland-island.pdf"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_pi_boxplots.R"

# ----------------------------------------------------------------------------------------
# --- Plot Fst against pi
# ----------------------------------------------------------------------------------------

rule plot_fst_by_pi:
    input:
        expand("results/chr{chr}.{site}.windowed.pi", chr=CHRS, site=ISLANDS_SPLIT_BUGALA),
        "data/inversion_sites.bed",
        "data/heterochromatin.bed"
    output:
        "reports/fst_by_pi.pdf"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_fst_against_pi.R"

# ----------------------------------------------------------------------------------------
# --- Compute Tajima's D
# ----------------------------------------------------------------------------------------

rule compute_tajimad:
    input:
        invcf="data/chr{chr}.pass.snp.flt.vcf.gz",
        indlist="data/ssese.seqids.is-{site}.txt"
    output:
        "results/chr{chr}.{site}.Tajima.D"
    threads: 1
    params: runtime="8",
            mem=",mem=5gb"
    shell:
        "module load vcftools;"
        "OUT_PREFIX=results/chr{wildcards.chr}.{wildcards.site};"
        "vcftools --gzvcf {input.invcf} "
        "    --TajimaD 10000 "
        "    --keep {input.indlist} "
        "    --out $OUT_PREFIX"

# ----------------------------------------------------------------------------------------
# --- Plot Tajima's D boxplots
# ----------------------------------------------------------------------------------------

rule plot_tajimas_d_by_site:
    input:
        expand("results/chr{chr}.{site}.Tajima.D", chr=['3L', '3R'], site=ISLANDS_SPLIT_BUGALA),
    output:
        "reports/tajimas_d_boxplot.pdf",
        "reports/tajimas_d_boxplot.mainland-island.pdf"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_tajimas_d.R"

# ----------------------------------------------------------------------------------------
# --- Plot Tajima's D across genome
# ----------------------------------------------------------------------------------------

rule plot_tajimas_d_across_genome:
    input:
        expand("results/chr{chr}.{site}.Tajima.D", chr=['3L', '3R'], site=ISLANDS_SPLIT_BUGALA),
    output:
        "reports/tajimas_D_by_site.pdf"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_tajD_by_site.R"

# ----------------------------------------------------------------------------------------
# --- Compute inbreeding stat (F)
# ----------------------------------------------------------------------------------------

rule compute_inbreeding:
    input:
        invcf="data/all.pass.snp.flt.vcf.gz",
        invhet="data/inversion_heterochromatin.bed",
    output:
        "results/all.{site}.het"
    threads: 1
    params: runtime="8",
            mem=",mem=5gb"
    shell:
        "module load vcftools;"
        "OUT_PREFIX=results/all.{wildcards.site};"
        "vcftools --gzvcf {input.invcf} "
        "    --het "
        "    --exclude-bed {input.invhet} "
        "    --keep data/ssese.seqids.is-{wildcards.site}.txt "
        "    --out $OUT_PREFIX"

# ----------------------------------------------------------------------------------------
# --- Compute allele frequency
# ----------------------------------------------------------------------------------------

rule compute_allele_freq:
    input:
        invcf="data/chr{chr}.pass.snp.flt.vcf.gz",
        invhet="data/inversion_heterochromatin.bed",
    output:
        "results/chr{chr}.{site}.frq"
    threads: 1
    params: runtime="8",
            mem=",mem=5gb"
    shell:
        "module load vcftools;"
        "OUT_PREFIX=results/chr{wildcards.chr}.{wildcards.site};"
        "vcftools --gzvcf {input.invcf} "
        "    --freq "
        "    --exclude-bed {input.invhet} "
        "    --keep data/ssese.seqids.is-{wildcards.site}.txt "
        "    --out $OUT_PREFIX"

# ----------------------------------------------------------------------------------------
# --- Take random sample of allele frequency results for plotting
# ----------------------------------------------------------------------------------------

rule sample_allele_freq:
    input:
        "results/chr{chr}.{site}.frq"
    output:
        "results/chr{chr}.{site}.sample.frq"
    threads: 1
    params: runtime="8",
            mem=",mem=5gb"
    shell:
        # head ends right away with error code 141, so capture it
        "shuf {input} | head -n 10000 > {output} || if [[ $? -eq 141 ]]; then true; else exit $?; fi"

# ----------------------------------------------------------------------------------------
# --- Plot MAF boxplots
# ----------------------------------------------------------------------------------------

rule plot_maf:
    input:
        expand("results/chr{chr}.{site}.sample.frq", chr=['3L', '3R'], site=ISLANDS_SPLIT_BUGALA),
    output:
        "reports/MAF_boxplot.pdf",
        "reports/MAF_boxplot.mainland-island.pdf",
        "reports/MAF_histogram_SFS.pdf"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_maf.R"

# ========================================================================================
# --- Phase with SHAPEIT
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Phase with SHAPEIT
# ----------------------------------------------------------------------------------------

rule phase:
    input:
        "data/chr{chr}.pass.snp.flt.vcf.gz",
    output:
        "data/chr{chr}.pass.snp.phased.haps"
    threads: 12
    params: runtime="48",
            mem=",mem=5gb"
    shell:
        "module load shapeit;"
        "sh scripts/phase_snps_with_shapeit.sh {wildcards.chr}"

# ----------------------------------------------------------------------------------------
# --- Find SNPs in inversions or het regions to exclude
# ----------------------------------------------------------------------------------------

rule find_inv_het_snps_in_haps_files:
    input:
        "data/chr{chr}.pass.snp.phased.haps",
        inv="data/inversion_simple.bed",
        het="data/heterochromatin.bed"
    output:
        "data/chr{chr}.pass.snp.phased.haps.excluded.snps"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load shapeit;"
        "sh scripts/find_inv_het_snps_in_haps_files.sh {wildcards.chr}"

# ----------------------------------------------------------------------------------------
# --- Make version of haps file with no het or inversion regions
# ----------------------------------------------------------------------------------------

rule remove_inv_het_snps_from_haps_files:
    input:
        inhaps="data/chr{chr}.pass.snp.phased.haps",
        badsnps="data/chr{chr}.pass.snp.phased.haps.excluded.snps"
    output:
        "data/nohetinv/chr{chr}.pass.snp.phased.nohetinv.haps",
        "data/nohetinv/chr{chr}.pass.snp.phased.nohetinv.sample"
    threads: 1
    params: runtime="24",
            mem=",mem=5gb"
    shell:
        "module load shapeit;"
        "mkdir -p data/nohetinv;"
        "shapeit -convert "
        "    --input-haps data/chr{wildcards.chr}.pass.snp.phased "
        "    --output-haps data/nohetinv/chr{wildcards.chr}.pass.snp.phased.nohetinv "
        "    --exclude-snp {input.badsnps};"
        "cp data/chr{wildcards.chr}.pass.snp.phased.sample "
        "    data/nohetinv/chr{wildcards.chr}.pass.snp.phased.nohetinv.sample"

# ----------------------------------------------------------------------------------------
# --- Make population haps files
# ----------------------------------------------------------------------------------------

rule make_pop_haps_files:
    input:
        inhap="data/nohetinv/chr{chr}.pass.snp.phased.nohetinv.haps",
        ind_list="data/ssese.seqids.is-{site}.txt"
    output:
        "data/chr{chr}.pass.snp.phased.{site}.haps",
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load shapeit;"
        "sh scripts/make_pop_haps_files.sh {wildcards.chr} {input.ind_list} {wildcards.site};"

# ----------------------------------------------------------------------------------------
# --- Make phased VCF
# ----------------------------------------------------------------------------------------

rule make_phased_vcf:
    input:
        inhap="data/chr{chr}.pass.snp.phased.haps",
    output:
        "data/chr{chr}.pass.snp.phased.haps.vcf",
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load shapeit;"
        "IN_HAPS={input.inhap};"
        "PREFIX=${{IN_HAPS/.haps/}};"
        "shapeit -convert --input-haps $PREFIX --output-vcf {output}"

# ----------------------------------------------------------------------------------------
# --- Subsample to make r^2 input files
# ----------------------------------------------------------------------------------------

rule sample_for_r2:
    input:
        vcf="data/chr{chr}.pass.snp.phased.haps.vcf"
    output:
        "results/chr{chr}.passing_SNPs_for_r2.full.recode.vcf",
        "results/chr{chr}.passing_SNPs_for_r2.sample.txt"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/sample_SNPs_for_r2_work.sh {wildcards.chr}"

# ----------------------------------------------------------------------------------------
# --- Compute r^2
# ----------------------------------------------------------------------------------------

rule compute_r2:
    input:
        "results/chr{chr}.passing_SNPs_for_r2.full.recode.vcf",
        "results/chr{chr}.passing_SNPs_for_r2.sample.txt",
        "data/ssese.seqids.is-{site}.txt",
        "data/inversion_heterochromatin.bed"
    output:
        "results/chr{chr}.{site}.list.hap.ld"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/compute_r2.sh {wildcards.chr} {wildcards.site}"

# ----------------------------------------------------------------------------------------
# --- Plot r^2 boxplots and LD decay
# ----------------------------------------------------------------------------------------

rule plot_r2:
    input:
        expand("results/chr{chr}.{site}.list.hap.ld", chr=['3L', '3R'], site=ISLANDS_SPLIT_BUGALA),
    output:
        "reports/r2_boxplot.pdf",
        "reports/r2_boxplot.mainland-island.pdf",
        "reports/r2_decay.pdf"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_r2_boxplots.R"

# ========================================================================================
# --- Run haplotype-based tests of selection
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Make haps files for JP XP-EHH
# ----------------------------------------------------------------------------------------

rule prep_for_jp_xpehh:
    input:
        expand("data/chr{{chr}}.pass.snp.phased.{site}.haps",
            site=ISLANDS_SPLIT_BUGALA)
    output:
        expand("data/xpehh_input_jp/chr{{chr}}.pass.snp.phased.{site}.jp.haps.transpose",
            site=ISLANDS_SPLIT_BUGALA),
        "data/xpehh_input_jp/chr{chr}.pass.snp.phased.jp.map"
    threads: 1
    params: runtime="24",
            mem=",mem=5gb"
    shell:
        "sh scripts/prep_for_xpehh_JP.sh {wildcards.chr}"

# ----------------------------------------------------------------------------------------
# --- Compute XP-EHH
# ----------------------------------------------------------------------------------------

rule compute_xpehh:
    input:
        lambda wildcards: XPEHH_INPUT["results/selscan/xp-ehh." + wildcards.site_pair + "." + wildcards.chr + ".xpehh.out"]
    output:
        "results/selscan/xp-ehh.{site_pair}.{chr}.xpehh.out"
    threads: 8
    params: runtime="24",
            mem=",mem=5gb"
    shell:
        "BOTH_POPS=`echo {wildcards.site_pair} | tr '-' ' '`;"
        "module load selscan;"
        "sh scripts/compute_XP-EHH.sh {wildcards.chr} $BOTH_POPS"

# ----------------------------------------------------------------------------------------
# --- Normalize XP-EHH
# ----------------------------------------------------------------------------------------

# Sentinel (for this chromosome) depends on all *.xpehh.out files (for this chromosome)
rule norm_xpehh:
    input:
        expand("results/selscan/xp-ehh.{site_pair}.{{chr}}.xpehh.out", site_pair=map(lambda x: x.replace("_", "-"), SITE_PAIRS)),
    output:
        "results/selscan/xp-ehh.{chr}.norm.sentinel.txt",
        expand("results/selscan/xp-ehh.{site_pair}.{{chr}}.xpehh.out.norm", site_pair=map(lambda x: x.replace("_", "-"), SITE_PAIRS)),
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load gcc;"
        "module load selscan;"
        "norm --xpehh --files {input};"
        "touch {output}"

# ----------------------------------------------------------------------------------------
# --- Window average XP-EHH
# ----------------------------------------------------------------------------------------

rule avg_xpehh:
    input:
        "results/selscan/xp-ehh.{site_pair}.{chr}.xpehh.out.norm"
    output:
        "results/selscan/xp-ehh.{site_pair}.{chr}.xpehh.out.norm.avg"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/window_avg_xpehh.R {input}"

# ----------------------------------------------------------------------------------------
# --- Plot XP-EHH and Fst (and CSS)
# ----------------------------------------------------------------------------------------

rule plot_xpehh_fst:
    input:
        expand("results/selscan/xp-ehh.{site_pair}.{chr}.xpehh.out.norm.avg",
            chr=CHRS,
            site_pair=map(lambda x: x.replace("_", "-"), SITE_PAIRS)),
        expand("results/chr{chr}.{site_pair}.windowed.weir.fst",
            chr=CHRS, site_pair=SITE_PAIRS),
        "results/css.all.txt",
        "data/inversion_simple.bed"
    output:
        "reports/fst_xpehh_by_site.pdf"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_stacked_fst_or_xpehh.R"

# ----------------------------------------------------------------------------------------
# --- Plot XP-EHH alone
# ----------------------------------------------------------------------------------------

rule plot_xpehh:
    input:
        "results/selscan/xp-ehh.{site_pair}.{chr}.xpehh.out.norm"
    output:
        "reports/selscan/xp-ehh.{site_pair}.{chr}.xpehh.out.norm.pdf"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_xp-ehh.R {input}"

# ----------------------------------------------------------------------------------------
# --- Rasterize PDFs of XP-EHH and Fst
# ----------------------------------------------------------------------------------------

rule rasterize_pdfs:
    input:
        expand("reports/xpehh_and_fst_between_populations_{chr}.all.pdf", chr=CHRS),
    output:
        expand("reports/xpehh_and_fst_between_populations_{chr}.all.png", chr=CHRS),
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/rasterize_pdfs.sh"

# ----------------------------------------------------------------------------------------
# --- Find peaks in XP-EHH
# ----------------------------------------------------------------------------------------

rule find_xpehh_peaks:
    input:
        expand("results/selscan/xp-ehh.{site_pair}.{chr}.xpehh.out.norm",
            site_pair=map(lambda x: x.replace("_", "-"), SITE_PAIRS), chr=CHRS),
    output:
        "results/xpehh_peaks.txt",
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/find_peaks_in_xp-ehh.R > {output}"

# ----------------------------------------------------------------------------------------
# --- Compute CSS
# ----------------------------------------------------------------------------------------

rule compute_css:
    input:
        expand("results/selscan/xp-ehh.{site_pair}.{chr}.xpehh.out.norm.avg",
            site_pair=map(lambda x: x.replace("_", "-"), SITE_PAIRS), chr=CHRS),
        expand("results/chr{chr}.{site_pair}.windowed.weir.fst",
            site_pair=SITE_PAIRS, chr=CHRS),
    output:
        "results/css.all.txt"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/compute_css.R"

# ----------------------------------------------------------------------------------------
# --- Plot CSS
# ----------------------------------------------------------------------------------------

rule plot_css:
    input:
        "results/css.all.txt",
        "data/inversion_simple.bed"
    output:
        "reports/css_{chr}.all.pdf"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_css.R {wildcards.chr}"

# ----------------------------------------------------------------------------------------
# --- Compute iHS
# ----------------------------------------------------------------------------------------

rule compute_ihs:
    input:
        "data/chr{chr}.pass.snp.phased.{site}.haps"
    output:
        "results/selscan/iHS.{site}.{chr}.ihs.out"
    threads: 8
    params: runtime="8",
            mem=",mem=5gb"
    shell:
        "sh scripts/compute_iHS.sh {wildcards.chr} {wildcards.site}"

# ----------------------------------------------------------------------------------------
# --- Normalize iHS
# ----------------------------------------------------------------------------------------

# Sentinel (for this chromosome) depends on all *.ihs.out files (for this chromosome)
rule norm_ihs:
    input:
        expand("results/selscan/iHS.{site}.{{chr}}.ihs.out", site=ISLANDS_SPLIT_BUGALA),
    output:
        "results/selscan/iHS.{chr}.norm.sentinel.txt",
        expand("results/selscan/iHS.{site}.{{chr}}.ihs.out.100bins.norm",
            site=ISLANDS_SPLIT_BUGALA),
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load gcc;"
        "module load selscan;"
        "norm --ihs --files {input};"
        "touch {output}"

# ----------------------------------------------------------------------------------------
# --- Plot iHS
# ----------------------------------------------------------------------------------------

rule plot_ihs:
    input:
        expand("results/selscan/iHS.{site}.{{chr}}.ihs.out.100bins.norm",
            site=ISLANDS),
        "data/inversion_simple.bed"
    output:
        "reports/ihs_{chr}.all.pdf"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_iHS.R {wildcards.chr}"

# ----------------------------------------------------------------------------------------
# --- Prepare input files for H12 statistic computation - all individuals
# ----------------------------------------------------------------------------------------

rule make_h12_input_all:
    input:
        "data/nohetinv/chr{chr}.pass.snp.phased.nohetinv.haps",
    output:
        "data/chr{chr}.pass.snp.phased.forH12.txt"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "perl scripts/make_h12_input.pl {input} > {output}"

# ----------------------------------------------------------------------------------------
# --- Prepare input files for H12 statistic computation - by population
# ----------------------------------------------------------------------------------------

rule make_h12_input_bysite:
    input:
        "data/chr{chr}.pass.snp.phased.{site}.haps",
    output:
        "data/chr{chr}.pass.snp.phased.{site}.forH12.txt"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "perl scripts/make_h12_input.pl {input} > {output}"

# ----------------------------------------------------------------------------------------
# --- Download H12 scripts
# ----------------------------------------------------------------------------------------

rule download_h12_script:
    output:
        "scripts/H12_H2H1.py"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "wget -O {output} https://raw.githubusercontent.com/ngarud/SelectionHapStats/master/scripts/H12_H2H1.py"

# ----------------------------------------------------------------------------------------
# --- Compute H12 statistic - All
# ----------------------------------------------------------------------------------------

rule compute_h12_all:
    input:
        in_file="data/chr{chr}.pass.snp.phased.forH12.txt",
        script="scripts/H12_H2H1.py"
    output:
        "results/chr{chr}.pass.snp.phased.H12.out.txt",
    threads: 1
    params: runtime="8",
            mem=",mem=5gb"
    shell:
        "module load python/2.7.8;"
        "COL_CT=`awk -F',' '{{print NF; exit}}' {input.in_file}`;"
        "IND_CT=$((COL_CT - 1));"
        "python {input.script} {input.in_file} $IND_CT > {output} || true;"

# ----------------------------------------------------------------------------------------
# --- Compute H12 statistic - by site
# ----------------------------------------------------------------------------------------

rule compute_h12_bysite:
    input:
        in_file="data/chr{chr}.pass.snp.phased.{site}.forH12.txt",
        script="scripts/H12_H2H1.py"
    output:
        "results/chr{chr}.pass.snp.phased.{site}.H12.out.txt",
    threads: 1
    params: runtime="8",
            mem=",mem=5gb"
    shell:
        "module load python/2.7.8;"
        "COL_CT=`awk -F',' '{{print NF; exit}}' {input.in_file}`;"
        "IND_CT=$((COL_CT - 1));"
        "python {input.script} {input.in_file} $IND_CT > {output} || true;"

# ----------------------------------------------------------------------------------------
# --- Plot H12 statistic - All
# ----------------------------------------------------------------------------------------

rule plot_h12_all:
    input:
        expand("results/chr{chr}.pass.snp.phased.H12.out.txt", chr=CHRS),
    output:
        "reports/H12.pdf"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_H12.R"

# ----------------------------------------------------------------------------------------
# --- Plot H12 statistic - by site
# ----------------------------------------------------------------------------------------

rule plot_h12_bysite:
    input:
        expand("results/chr{chr}.pass.snp.phased.H12.out.txt", chr=CHRS)
    output:
        "reports/H12_by_site.pdf"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_H12_bysite.R"

# ----------------------------------------------------------------------------------------
# --- Find H12 outliers
# ----------------------------------------------------------------------------------------

rule find_h12_outliers:
    input:
        expand("results/chr{chr}.pass.snp.phased.H12.out.txt", chr=CHRS)
    output:
        "reports/island-specific-sweeps.tex",
        "reports/site-specific-sweeps.tex",
        "reports/insecticide-sweeps.tex"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        #"module load r/3.4;"
        "Rscript scripts/find_h12_outliers.R"

# ========================================================================================
# --- Explore structure with DAPC
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Prepare to do DAPC
# ----------------------------------------------------------------------------------------

rule prep_for_dapc:
    input:
        "data/inversion_simple.bed",
        "data/heterochromatin.bed",
        "data/all.pass.snp.flt.prune.out",
    output:
        expand("results/dapc_subset_snps.{chr_subset}.raw",
            chr_subset=["chr3L", "chr3R"])
    threads: 1
    params: runtime="48",
            mem=",mem=24gb"
    shell:
        "sh scripts/do_dapc_prep.sh"

# ----------------------------------------------------------------------------------------
# --- Do DAPC
# ----------------------------------------------------------------------------------------

rule do_dapc:
    input:
        "results/dapc_subset_snps.{chr_subset}.raw"
    output:
        "results/dapc_subset_snps.{chr_subset}.dapc.scatter.pdf",
        "results/dapc_subset_snps.{chr_subset}.clusters.by.site.pdf",
        "results/dapc_subset_snps.{chr_subset}.pca.colorplot.pdf",
        "results/dapc_subset_snps.{chr_subset}.BIC.pdf",
    threads: 12
    params: runtime="48",
            mem=",mem=24gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/do_dapc_selected_SNPs.R results/dapc_subset_snps.{wildcards.chr_subset}"

# ========================================================================================
# --- Bring in Ag1000G data
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Download Ag1000G data
# ----------------------------------------------------------------------------------------

rule download_ag1000g:
    output:
        expand("data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.pass.biallelic.{chr}.vcf.gz", chr=CHRS),
        expand("data/accessibility/accessibility.{chr}.vcf.gz", chr=CHRS)
    threads: 1
    params: runtime="12",
            mem=",mem=5gb"
    shell:
        "sh scripts/download_ag1000g_VCFs.sh"

# ----------------------------------------------------------------------------------------
# --- Merge Ssese and Ag1000G VCF files
# ----------------------------------------------------------------------------------------

rule merge_ag1000g:
    input:
        "data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.pass.biallelic.{chr}.vcf.gz",
        "data/accessibility/accessibility.{chr}.vcf.gz",
        "data/chr{chr}.pass.snp.vcf.gz"
    output:
        "data/ssese_with_ag1000g/ssese_with_ag1000g.{chr}.flt.strict.vcf.gz",
    threads: 1
    params: runtime="24",
            mem=",mem=5gb"
    shell:
        "sh scripts/merge_ssese_ag1000g_VCFs.sh {wildcards.chr}"

rule merge_ag1000g_missingtoref:
    input:
        "data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.pass.biallelic.{chr}.vcf.gz",
        "data/accessibility/accessibility.{chr}.vcf.gz",
        "data/chr{chr}.pass.snp.vcf.gz"
    output:
        "data/ssese_with_ag1000g.missingtoref/ssese_with_ag1000g.{chr}.missingtoref.flt.strict.vcf.gz"
    threads: 1
    params: runtime="24",
            mem=",mem=5gb"
    shell:
        "sh scripts/merge_ssese_ag1000g_VCFs.sh {wildcards.chr} .missingtoref"

# -------------------------------------------------------------------------------------- #
# --- Convert Ssese+Ag100G VCF to plink's BED format (via temporary PED)
# -------------------------------------------------------------------------------------- #

rule ss_ag_vcf_to_bed:
    input:
        invcf="data/ssese_with_ag1000g/ssese_with_ag1000g.{chr}.flt.strict.vcf.gz",
    output:
        "data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.bed"
    threads: 1
    params: runtime="24",
            mem=",mem=5gb"
    shell:
        "module load plink;"
        "module load vcftools;"
        "OUT_FILE={output};"
        "mkdir -p tmp;"
        "TMP=`pwd`/tmp;"
        "vcftools --gzvcf {input.invcf} --plink --out ${{OUT_FILE/.bed}} --temp $TMP;"
        # Edit the MAP file and change 0 to the right chromosome index:
        #   1   2L
        #   2   2R
        #   3   3L
        #   4   3R
        #   5   X
        "sh scripts/fix_chr_in_ped.sh ${{OUT_FILE/.bed}}.map;"
        "plink --file ${{OUT_FILE/.bed}} --make-bed --out ${{OUT_FILE/.bed}};"
        "rm ${{OUT_FILE/.bed}}.ped"

rule ss_ag_vcf_to_bed_missingtoref:
    input:
        invcf="data/ssese_with_ag1000g.missingtoref/ssese_with_ag1000g.{chr}.missingtoref.flt.strict.vcf.gz",
    output:
        "data/ssese_with_ag1000g.missingtoref/chr{chr}.missingtoref.pass.snp.phased.ag1000g.strict.bed"
    threads: 1
    params: runtime="24",
            mem=",mem=5gb"
    shell:
        "module load plink;"
        "module load vcftools;"
        "OUT_FILE={output};"
        "mkdir -p tmp;"
        "TMP=`pwd`/tmp;"
        "vcftools --gzvcf {input.invcf} --plink --out ${{OUT_FILE/.bed}} --temp $TMP;"
        # Edit the MAP file and change 0 to the right chromosome index:
        #   1   2L
        #   2   2R
        #   3   3L
        #   4   3R
        #   5   X
        "sh scripts/fix_chr_in_ped.sh ${{OUT_FILE/.bed}}.map;"
        "plink --file ${{OUT_FILE/.bed}} --make-bed --out ${{OUT_FILE/.bed}};"
        "rm ${{OUT_FILE/.bed}}.ped"

# -------------------------------------------------------------------------------------- #
# --- Ssese+Ag100G: Find SNPs in LD to prune out
# -------------------------------------------------------------------------------------- #

# Slide a sliding window of 50 SNPs across each chromosome and calculate r^2 for each
# pair in the window. Remove one random SNP in each pair with r^2 > 0.1.
# Parameters are window size, step, and r^2 threshold.

# To use this list of pruned SNPs to create an LD-pruned dataset, use:
#   plink --exclude data/ssese_with_ag1000g/chrX.pass.snp.phased.ag1000g.prune.out
#   vcftools --exclude-positions data/ssese_with_ag1000g/chrX.pass.snp.phased.ag1000g.prune.out.tab

rule ss_ag_LD_prune:
    input:
        "data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.bed"
    output:
        "data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.prune.out",
        "data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.prune.out.tab"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load plink;"
        "PREFIX=`echo {input} | sed -e \"s/\.bed//\"`;"
        "plink --bfile $PREFIX --indep-pairwise 50 5 0.1 --out $PREFIX;"
        "sed -e \"s/:/\\t/\" {output[0]} > {output[1]};"
        "mkdir -p reports/plink_LD_pruning_part1/;"
        "mv $PREFIX.log reports/plink_LD_pruning_part1/"

rule ss_ag_LD_prune_missingtoref:
    input:
        "data/ssese_with_ag1000g.missingtoref/chr{chr}.missingtoref.pass.snp.phased.ag1000g.strict.bed"
    output:
        "data/ssese_with_ag1000g.missingtoref/chr{chr}.missingtoref.pass.snp.phased.ag1000g.strict.prune.out",
        "data/ssese_with_ag1000g.missingtoref/chr{chr}.missingtoref.pass.snp.phased.ag1000g.strict.prune.out.tab"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load plink;"
        "PREFIX=`echo {input} | sed -e \"s/\.bed//\"`;"
        "plink --bfile $PREFIX --indep-pairwise 50 5 0.1 --out $PREFIX;"
        "sed -e \"s/:/\\t/\" {output[0]} > {output[1]};"
        "mkdir -p reports/plink_LD_pruning_part1/;"
        "mv $PREFIX.log reports/plink_LD_pruning_part1/"

# -------------------------------------------------------------------------------------- #
# --- Ssese+Ag100G: Combine lists of SNPs in LD to prune out
# -------------------------------------------------------------------------------------- #

# To use this list of pruned SNPs to create an LD-pruned dataset, use:
#   plink --exclude data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.prune.out
#   vcftools --exclude-positions data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.prune.out.tab

rule ss_ag_combine_LD_prune:
    input:
        expand("data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.prune.out", chr=CHR_THREE),
        expand("data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.prune.out.tab", chr=CHR_THREE)
    output:
        "data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.prune.out",
        "data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.prune.out.tab"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "cat data/ssese_with_ag1000g/chr*.pass.snp.phased.ag1000g.strict.prune.out > data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.prune.out;"
        "cat data/ssese_with_ag1000g/chr*.pass.snp.phased.ag1000g.strict.prune.out.tab > data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.prune.out.tab;"

rule ss_ag_combine_LD_prune_missingtoref:
    input:
        expand("data/ssese_with_ag1000g.missingtoref/chr{chr}.missingtoref.pass.snp.phased.ag1000g.strict.prune.out", chr=CHR_THREE),
        expand("data/ssese_with_ag1000g.missingtoref/chr{chr}.missingtoref.pass.snp.phased.ag1000g.strict.prune.out.tab", chr=CHR_THREE)
    output:
        "data/ssese_with_ag1000g.missingtoref/all.missingtoref.pass.snp.phased.ag1000g.strict.prune.out",
        "data/ssese_with_ag1000g.missingtoref/all.missingtoref.pass.snp.phased.ag1000g.strict.prune.out.tab"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "cat data/ssese_with_ag1000g.missingtoref/chr*.missingtoref.pass.snp.phased.ag1000g.strict.prune.out > data/ssese_with_ag1000g.missingtoref/all.missingtoref.pass.snp.phased.ag1000g.strict.prune.out;"
        "cat data/ssese_with_ag1000g.missingtoref/chr*.missingtoref.pass.snp.phased.ag1000g.strict.prune.out.tab > data/ssese_with_ag1000g.missingtoref/all.missingtoref.pass.snp.phased.ag1000g.strict.prune.out.tab;"

# ========================================================================================
# --- Do MDS on combined Ssese + Ag1000G dataset
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Find lists of Ag1000G animals to exclude
# ----------------------------------------------------------------------------------------

rule ss_ag_find_inds_to_exclude:
    input:
        expand("data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.1.haplotypes.{chr}.sample", chr=CHR_THREE),
        "data/ag1000g.phase1.ar3/samples.all.txt"
    output:
        inds_sy="data/ag1000g.phase1.ar3/SY.individuals.txt",
        inds_ad="data/ag1000g.phase1.ar3/AD.individuals.txt",
        inds_m="data/ag1000g.phase1.ar3/M.individuals.txt",
        inds_all="data/ag1000g.phase1.ar3/excluded.individuals.txt",
        inds_all_not_M="data/ag1000g.phase1.ar3/excluded.individuals.notM.txt",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "cut -d' ' -f1 data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.1.haplotypes.*.sample | grep 'SY' | sort | uniq > {output.inds_sy};"
        "cut -d' ' -f1 data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.1.haplotypes.*.sample | grep 'AD' | sort | uniq > {output.inds_ad};"
        # Find M individuals from Ag1000G dataset
        "awk 'BEGIN {{ FS = \"\\t\" }} {{ if ( $11 == \"M\") print $2 }}' data/ag1000g.phase1.ar3/samples.all.txt > {output.inds_m};"
        "cat {output.inds_sy} {output.inds_ad} {output.inds_m} > {output.inds_all};"
        # Have another exclusion list that doesn't include coluzzi individuals
        "cat {output.inds_sy} {output.inds_ad} > {output.inds_all_not_M};"

# ----------------------------------------------------------------------------------------
# --- Ssese+Ag100G: Compute MDS from IBS values - full genome, LD pruned
# ----------------------------------------------------------------------------------------

rule ss_ag_mds_all:
    input:
        beds=expand("data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.bed", chr=CHR_THREE),
        prune="data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.prune.out",
        exclude="data/ag1000g.phase1.ar3/excluded.individuals.txt"
    output:
        mds="data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.mds",
        mds_eig="data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.mds.eigvals"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load plink;"
        # Temporary files and prefix (needed to not conflict with other jobs when writing
        # temporary BED, BIM, and FAM files)
        "TMP_BED_LIST=`mktemp tmp.chr_bed_list.XXXXX`;"
        "TMP_PREFIX=`mktemp -u data/tmp.all.pass.snp.flt.XXXXX`;"
        # Create list of input BED files (by chr)
        "ls data/ssese_with_ag1000g/chr*.pass.snp.phased.ag1000g.strict.bed | sed -e \"s/\(.*\).bed/\\1.bed \\1.bim \\1.fam/\" > $TMP_BED_LIST;"
        # Compute MDS
        "plink --merge-list $TMP_BED_LIST --exclude {input.prune} --remove-fam {input.exclude} --cluster --mds-plot 2 eigvals --out $TMP_PREFIX;"
        # Move results
        "mv $TMP_PREFIX.mds {output.mds};"
        "mv $TMP_PREFIX.mds.eigvals {output.mds_eig};"
        # Clean up
        "rm $TMP_BED_LIST;"
        "rm $TMP_PREFIX.*"

# ----------------------------------------------------------------------------------------
# --- Ssese+Ag100G: Compute MDS from IBS values - full genome w/out inversions, LD pruned
# ----------------------------------------------------------------------------------------

rule ss_ag_mds_all_no_inv:
    input:
        beds=expand("data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.bed", chr=CHR_THREE),
        prune="data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.prune.out",
        exclude="data/ag1000g.phase1.ar3/excluded.individuals.txt",
        inv="data/inversion_simple.bed"
    output:
        mds="data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.noinv.mds",
        mds_eig="data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.noinv.mds.eigvals"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load plink;"
        # Temporary files and prefix (needed to not conflict with other jobs when writing
        # temporary BED, BIM, and FAM files)
        "TMP_BED_LIST=`mktemp tmp.chr_bed_list.XXXXX`;"
        "TMP_INV_SIMP=`mktemp tmp.inversion_simple.XXXXX`;"
        "TMP_EXCLUDE=`mktemp tmp.exclude.XXXXX`;"
        "TMP_PREFIX=`mktemp -u data/tmp.all.pass.snp.flt.noinv.XXXXX`;"
        # Create list of input BED files (by chr) and simple inversions BED
        "ls data/ssese_with_ag1000g/chr*.pass.snp.phased.ag1000g.strict.bed | sed -e \"s/\(.*\).bed/\\1.bed \\1.bim \\1.fam/\" > $TMP_BED_LIST;"
        "sed -e \"s/^2L/1/\" -e \"s/^2R/2/\" data/inversion_simple.bed > $TMP_INV_SIMP;"
        # Make SNP set covering the 2La and 2Rb inversions
        "plink --merge-list $TMP_BED_LIST --exclude {input.prune} --make-set $TMP_INV_SIMP --write-set --out data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.mds.inversions;"
        # Compute MDS
        "cat {input.prune} data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.mds.inversions.set > $TMP_EXCLUDE;"
        "plink --merge-list $TMP_BED_LIST --exclude $TMP_EXCLUDE --remove-fam {input.exclude} --cluster --mds-plot 2 eigvals --out $TMP_PREFIX --write-snplist;"
        # Move result
        "mv $TMP_PREFIX.mds {output.mds};"
        "mv $TMP_PREFIX.mds.eigvals {output.mds_eig};"
        # Clean up
        "rm $TMP_BED_LIST;"
        "rm $TMP_EXCLUDE;"
        "rm $TMP_INV_SIMP;"
        "rm data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.mds.inversions.*;"
        "rm $TMP_PREFIX.*;"

# ----------------------------------------------------------------------------------------
# --- Ssese+Ag100G: Compute MDS from IBS values - by chromosome
# ----------------------------------------------------------------------------------------

rule ss_ag_mds_chr:
    input:
        bed="data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.bed",
        ld_set="data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.prune.out",
        exclude="data/ag1000g.phase1.ar3/excluded.individuals.txt"
    output:
        "data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.mds",
        "data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.mds.eigvals"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load plink;"
        "plink --bfile data/ssese_with_ag1000g/chr{wildcards.chr}.pass.snp.phased.ag1000g.strict --exclude {input.ld_set}  --remove-fam {input.exclude} --cluster --mds-plot 2 eigvals --out data/ssese_with_ag1000g/chr{wildcards.chr}.pass.snp.phased.ag1000g.strict;"

# ----------------------------------------------------------------------------------------
# --- Ssese+Ag100G: Plot MDS
# ----------------------------------------------------------------------------------------

rule ss_ag_plot_mds_all:
    input:
        "data/sample_to_seq_id_mapping.txt",
        "data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.mds",
        "data/ssese_individual_info.csv",
        "data/ag1000g.phase1.ar3/samples.all.txt"
    output:
        "reports/all.pass.snp.phased.ag1000g.strict.ibs_mds_plot.pdf",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_ibs_mds_or_pca_with_ag1000g.R data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.mds"

rule ss_ag_plot_mds_all_no_inv:
    input:
        "data/sample_to_seq_id_mapping.txt",
        "data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.noinv.mds",
        "data/ssese_individual_info.csv",
        "data/ag1000g.phase1.ar3/samples.all.txt"
    output:
        "reports/all.pass.snp.phased.ag1000g.strict.noinv.ibs_mds_plot.pdf",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_ibs_mds_or_pca_with_ag1000g.R data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.noinv.mds"

rule ss_ag_plot_mds_chr:
    input:
        "data/sample_to_seq_id_mapping.txt",
        "data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.mds",
        "data/ssese_individual_info.csv",
        "data/ag1000g.phase1.ar3/samples.all.txt"
    output:
        "reports/chr{chr}.pass.snp.phased.ag1000g.strict.ibs_mds_plot.pdf",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_ibs_mds_or_pca_with_ag1000g.R data/ssese_with_ag1000g/chr{wildcards.chr}.pass.snp.phased.ag1000g.strict.mds"

# ========================================================================================
# --- Ssese+Ag100G: Do PCA
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Ssese+Ag100G: Do PCA - chr3 only, LD pruned
# ----------------------------------------------------------------------------------------

rule ss_ag_pca_all_no_inv:
    input:
        beds=expand("data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.bed", chr=CHR_THREE),
        prune="data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.prune.out",
        exclude="data/ag1000g.phase1.ar3/excluded.individuals.txt",
        inv="data/inversion_simple.bed",
        het="data/heterochromatin.bed"
    output:
        eigenvec="data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.noinv.eigenvec",
        eigenval="data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.noinv.eigenval"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load plink;"
        # Temporary files and prefix (needed to not conflict with other jobs when writing
        # temporary BED, BIM, and FAM files)
        "TMP_BED_LIST=`mktemp tmp.chr_bed_list.XXXXX`;"
        "TMP_INV_HET_SIMP=`mktemp tmp.inversion_simple.XXXXX`;"
        "TMP_EXCLUDE=`mktemp tmp.exclude.XXXXX`;"
        "TMP_PREFIX=`mktemp -u data/tmp.all.pass.snp.flt.noinv.XXXXX`;"
        # Create list of input BED files (by chr) and simple inversions BED
        "ls data/ssese_with_ag1000g/chr3*.pass.snp.phased.ag1000g.strict.bed | sed -e \"s/\(.*\).bed/\\1.bed \\1.bim \\1.fam/\" > $TMP_BED_LIST;"
        "cat data/inversion_simple.bed | "
        "    sed -e \"s/^2L/1/\" -e \"s/^2R/2/\" -e \"s/^3L/3/\" -e \"s/^3R/4/\" > $TMP_INV_HET_SIMP;"
        "cat data/heterochromatin.bed |  "
        "    sed -e \"s/^2L/1/\" -e \"s/^2R/2/\" -e \"s/^3L/3/\" -e \"s/^3R/4/\" | "
        "    sed -e \"s/$/\\thet/\" >> $TMP_INV_HET_SIMP;"
        # Make SNP set covering the 2La and 2Rb inversions
        "plink --merge-list $TMP_BED_LIST --exclude {input.prune} --make-set $TMP_INV_HET_SIMP --write-set --out data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.pca.inversions;"
        # Do PCA with 100 PCs retained
        "cat {input.prune} data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.pca.inversions.set > $TMP_EXCLUDE;"
        "plink --merge-list $TMP_BED_LIST --exclude $TMP_EXCLUDE --remove-fam {input.exclude} --maf 0.01 --cluster --pca 100 --out $TMP_PREFIX;"
        # Move result
        "mv $TMP_PREFIX.eigenvec {output.eigenvec}.forScree;"
        "mv $TMP_PREFIX.eigenval {output.eigenval}.forScree;"
        # Do PCA with 15 PCs retained to match Miles et al
        "plink --merge-list $TMP_BED_LIST --exclude $TMP_EXCLUDE --remove-fam {input.exclude} --maf 0.01 --cluster --pca 15 --out $TMP_PREFIX;"
        # Move result
        "mv $TMP_PREFIX.eigenvec {output.eigenvec};"
        "mv $TMP_PREFIX.eigenval {output.eigenval};"
        # Clean up
        "rm $TMP_BED_LIST;"
        "rm $TMP_EXCLUDE;"
        "rm $TMP_INV_HET_SIMP;"
        "rm data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.pca.inversions.*;"
        "rm $TMP_PREFIX.*;"

# ----------------------------------------------------------------------------------------
# --- Ssese+Ag100G: Do PCA - chr3 only, LD pruned - With coluzzi individuals
# ----------------------------------------------------------------------------------------

rule ss_ag_pca_all_no_inv_withM:
    input:
        beds=expand("data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.bed", chr=CHR_THREE),
        prune="data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.prune.out",
        exclude="data/ag1000g.phase1.ar3/excluded.individuals.notM.txt",
        inv="data/inversion_simple.bed",
        het="data/heterochromatin.bed"
    output:
        eigenvec="data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.noinv.withM.eigenvec",
        eigenval="data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.noinv.withM.eigenval"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load plink;"
        # Temporary files and prefix (needed to not conflict with other jobs when writing
        # temporary BED, BIM, and FAM files)
        "TMP_BED_LIST=`mktemp tmp.chr_bed_list.XXXXX`;"
        "TMP_INV_HET_SIMP=`mktemp tmp.inversion_simple.XXXXX`;"
        "TMP_EXCLUDE=`mktemp tmp.exclude.XXXXX`;"
        "TMP_PREFIX=`mktemp -u data/tmp.all.pass.snp.flt.noinv.XXXXX`;"
        # Create list of input BED files (by chr) and simple inversions BED
        "ls data/ssese_with_ag1000g/chr3*.pass.snp.phased.ag1000g.strict.bed | sed -e \"s/\(.*\).bed/\\1.bed \\1.bim \\1.fam/\" > $TMP_BED_LIST;"
        "cat data/inversion_simple.bed | "
        "    sed -e \"s/^2L/1/\" -e \"s/^2R/2/\" -e \"s/^3L/3/\" -e \"s/^3R/4/\" > $TMP_INV_HET_SIMP;"
        "cat data/heterochromatin.bed |  "
        "    sed -e \"s/^2L/1/\" -e \"s/^2R/2/\" -e \"s/^3L/3/\" -e \"s/^3R/4/\" | "
        "    sed -e \"s/$/\\thet/\" >> $TMP_INV_HET_SIMP;"
        # Make SNP set covering the 2La and 2Rb inversions (Named *.withM.* just to avoid conflict with other PCA run)
        "plink --merge-list $TMP_BED_LIST --exclude {input.prune} --make-set $TMP_INV_HET_SIMP --write-set --out data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.pca.withM.inversions;"
        # Do PCA with 100 PCs retained
        "cat {input.prune} data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.pca.withM.inversions.set > $TMP_EXCLUDE;"
        "plink --merge-list $TMP_BED_LIST --exclude $TMP_EXCLUDE --remove-fam {input.exclude} --maf 0.01 --cluster --pca 100 --out $TMP_PREFIX;"
        # Move result
        "mv $TMP_PREFIX.eigenvec {output.eigenvec}.forScree;"
        "mv $TMP_PREFIX.eigenval {output.eigenval}.forScree;"
        # Do PCA with 15 PCs retained to match Miles et al
        "plink --merge-list $TMP_BED_LIST --exclude $TMP_EXCLUDE --remove-fam {input.exclude} --maf 0.01 --cluster --pca 15 --out $TMP_PREFIX;"
        # Move result
        "mv $TMP_PREFIX.eigenvec {output.eigenvec};"
        "mv $TMP_PREFIX.eigenval {output.eigenval};"
        # Clean up
        "rm $TMP_BED_LIST;"
        "rm $TMP_EXCLUDE;"
        "rm $TMP_INV_HET_SIMP;"
        "rm data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.pca.withM.inversions.*;"
        "rm $TMP_PREFIX.*;"

rule ss_ag_pca_all_no_inv_withM_noKES:
    input:
        beds=expand("data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.bed", chr=CHR_THREE),
        prune="data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.prune.out",
        exclude="data/ag1000g.phase1.ar3/excluded.individuals.notM.txt",
        inv="data/inversion_simple.bed",
        het="data/heterochromatin.bed"
    output:
        eigenvec="data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.noinv.withMnoKES.eigenvec",
        eigenval="data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.noinv.withMnoKES.eigenval"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load plink;"
        # Temporary files and prefix (needed to not conflict with other jobs when writing
        # temporary BED, BIM, and FAM files)
        "TMP_BED_LIST=`mktemp tmp.chr_bed_list.XXXXX`;"
        "TMP_INV_HET_SIMP=`mktemp tmp.inversion_simple.XXXXX`;"
        "TMP_EXCLUDE=`mktemp tmp.exclude.XXXXX`;"
        "TMP_PREFIX=`mktemp -u data/tmp.all.pass.snp.flt.noinv.XXXXX`;"
        # Create list of input BED files (by chr) and simple inversions BED
        "ls data/ssese_with_ag1000g/chr3*.pass.snp.phased.ag1000g.strict.bed | sed -e \"s/\(.*\).bed/\\1.bed \\1.bim \\1.fam/\" > $TMP_BED_LIST;"
        "cat data/inversion_simple.bed | "
        "    sed -e \"s/^2L/1/\" -e \"s/^2R/2/\" -e \"s/^3L/3/\" -e \"s/^3R/4/\" > $TMP_INV_HET_SIMP;"
        "cat data/heterochromatin.bed |  "
        "    sed -e \"s/^2L/1/\" -e \"s/^2R/2/\" -e \"s/^3L/3/\" -e \"s/^3R/4/\" | "
        "    sed -e \"s/$/\\thet/\" >> $TMP_INV_HET_SIMP;"
        # Make SNP set covering the 2La and 2Rb inversions (Named *.withMnoKES.* just to avoid conflict with other PCA run)
        "plink --merge-list $TMP_BED_LIST --exclude {input.prune} --make-set $TMP_INV_HET_SIMP --write-set --out data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.pca.withMnoKES.inversions;"
        # Exclude normal excluded individuals, plus KES
        "TMP_TOREMOVE=`mktemp tmp.toremove.XXXXX`;"
        "cat {input.exclude} > $TMP_TOREMOVE;"
        "grep 'Kenya' data/ag1000g.phase1.ar3/samples.all.txt | cut -f2 >> $TMP_TOREMOVE;"
        # Do PCA with 100 PCs retained
        "cat {input.prune} data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.pca.withMnoKES.inversions.set > $TMP_EXCLUDE;"
        "plink --merge-list $TMP_BED_LIST --exclude $TMP_EXCLUDE --remove-fam $TMP_TOREMOVE --maf 0.01 --cluster --pca 100 --out $TMP_PREFIX;"
        # Move result
        "mv $TMP_PREFIX.eigenvec {output.eigenvec}.forScree;"
        "mv $TMP_PREFIX.eigenval {output.eigenval}.forScree;"
        # Do PCA with 15 PCs retained to match Miles et al
        "plink --merge-list $TMP_BED_LIST --exclude $TMP_EXCLUDE --remove-fam $TMP_TOREMOVE --maf 0.01 --cluster --pca 15 --out $TMP_PREFIX;"
        # Move result
        "mv $TMP_PREFIX.eigenvec {output.eigenvec};"
        "mv $TMP_PREFIX.eigenval {output.eigenval};"
        # Clean up
        "rm $TMP_BED_LIST;"
        "rm $TMP_EXCLUDE;"
        "rm $TMP_INV_HET_SIMP;"
        "rm $TMP_TOREMOVE;"
        "rm data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.pca.withMnoKES.inversions.*;"
        "rm $TMP_PREFIX.*;"

rule ss_ag_pca_all_no_inv_withM_missingtoref:
    input:
        beds=expand("data/ssese_with_ag1000g.missingtoref/chr{chr}.missingtoref.pass.snp.phased.ag1000g.strict.bed", chr=CHR_THREE),
        prune="data/ssese_with_ag1000g.missingtoref/all.missingtoref.pass.snp.phased.ag1000g.strict.prune.out",
        exclude="data/ag1000g.phase1.ar3/excluded.individuals.notM.txt",
        inv="data/inversion_simple.bed",
        het="data/heterochromatin.bed"
    output:
        eigenvec="data/ssese_with_ag1000g.missingtoref/chr3.pass.snp.phased.ag1000g.strict.noinv.withM.eigenvec",
        eigenval="data/ssese_with_ag1000g.missingtoref/chr3.pass.snp.phased.ag1000g.strict.noinv.withM.eigenval"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load plink;"
        # Temporary files and prefix (needed to not conflict with other jobs when writing
        # temporary BED, BIM, and FAM files)
        "TMP_BED_LIST=`mktemp tmp.chr_bed_list.XXXXX`;"
        "TMP_INV_HET_SIMP=`mktemp tmp.inversion_simple.XXXXX`;"
        "TMP_EXCLUDE=`mktemp tmp.exclude.XXXXX`;"
        "TMP_PREFIX=`mktemp -u data/tmp.all.pass.snp.flt.noinv.XXXXX`;"
        # Create list of input BED files (by chr) and simple inversions BED
        "ls data/ssese_with_ag1000g.missingtoref/chr3*.pass.snp.phased.ag1000g.strict.bed | sed -e \"s/\(.*\).bed/\\1.bed \\1.bim \\1.fam/\" > $TMP_BED_LIST;"
        "cat data/inversion_simple.bed | "
        "    sed -e \"s/^2L/1/\" -e \"s/^2R/2/\" -e \"s/^3L/3/\" -e \"s/^3R/4/\" > $TMP_INV_HET_SIMP;"
        "cat data/heterochromatin.bed |  "
        "    sed -e \"s/^2L/1/\" -e \"s/^2R/2/\" -e \"s/^3L/3/\" -e \"s/^3R/4/\" | "
        "    sed -e \"s/$/\\thet/\" >> $TMP_INV_HET_SIMP;"
        # Make SNP set covering the 2La and 2Rb inversions (Named *.withM.* just to avoid conflict with other PCA run)
        "plink --merge-list $TMP_BED_LIST --exclude {input.prune} --make-set $TMP_INV_HET_SIMP --write-set --out data/ssese_with_ag1000g.missingtoref/all.pass.snp.phased.ag1000g.strict.pca.withM.inversions;"
        # Do PCA with 100 PCs retained
        "cat {input.prune} data/ssese_with_ag1000g.missingtoref/all.pass.snp.phased.ag1000g.strict.pca.withM.inversions.set > $TMP_EXCLUDE;"
        "plink --merge-list $TMP_BED_LIST --exclude $TMP_EXCLUDE --remove-fam {input.exclude} --maf 0.01 --cluster --pca 100 --out $TMP_PREFIX;"
        # Move result
        "mv $TMP_PREFIX.eigenvec {output.eigenvec}.forScree;"
        "mv $TMP_PREFIX.eigenval {output.eigenval}.forScree;"
        # Do PCA with 15 PCs retained to match Miles et al
        "plink --merge-list $TMP_BED_LIST --exclude $TMP_EXCLUDE --remove-fam {input.exclude} --maf 0.01 --cluster --pca 15 --out $TMP_PREFIX;"
        # Move result
        "mv $TMP_PREFIX.eigenvec {output.eigenvec};"
        "mv $TMP_PREFIX.eigenval {output.eigenval};"
        # Clean up
        "rm $TMP_BED_LIST;"
        "rm $TMP_EXCLUDE;"
        "rm $TMP_INV_HET_SIMP;"
        "rm data/ssese_with_ag1000g.missingtoref/all.pass.snp.phased.ag1000g.strict.pca.withM.inversions.*;"
        "rm $TMP_PREFIX.*;"

# ----------------------------------------------------------------------------------------
# --- Ssese+Ag100G: Plot PCA
# ----------------------------------------------------------------------------------------

rule ss_ag_plot_pca:
    input:
        "data/sample_to_seq_id_mapping.txt",
        "data/ssese_with_ag1000g/{prefix}.eigenvec",
        "data/ssese_individual_info.csv",
        "data/ag1000g.phase1.ar3/samples.all.txt"
    output:
        "reports/{prefix}.pca_plot.1vs2.pca_plot.pdf",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_ibs_mds_or_pca_with_ag1000g.R data/ssese_with_ag1000g/{wildcards.prefix}.eigenvec"

rule ss_ag_plot_pca_missingtoref:
    input:
        "data/sample_to_seq_id_mapping.txt",
        "data/ssese_with_ag1000g.missingtoref/chr3.pass.snp.phased.ag1000g.strict.noinv.withM.eigenvec",
        "data/ssese_individual_info.csv",
        "data/ag1000g.phase1.ar3/samples.all.txt"
    output:
        "reports/chr3.missingtoref.pass.snp.phased.ag1000g.strict.noinv.withM.pca_plot.1vs2.pca_plot.pdf"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_ibs_mds_or_pca_with_ag1000g.R data/ssese_with_ag1000g.missingtoref/chr3.pass.snp.phased.ag1000g.strict.noinv.withM.eigenvec"

rule ss_ag_plot_pca_nokes:
    input:
        "data/sample_to_seq_id_mapping.txt",
        "data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.noinv.withMnoKES.eigenvec",
        "data/ssese_individual_info.csv",
        "data/ag1000g.phase1.ar3/samples.all.txt"
    output:
        "reports/chr3.pass.snp.phased.ag1000g.strict.noinv.withMnoKES.pca_plot.1vs2.pca_plot.pdf"
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_ibs_mds_or_pca_with_ag1000g.R data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.noinv.withMnoKES.eigenvec"

# ========================================================================================
# --- Run dadi
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Make dadi input files (including just chr3)
# ----------------------------------------------------------------------------------------

rule make_dadi_input:
    input:
        "data/ssese_individual_info_simple_bugala_split.txt",
        expand("data/chr{chr}.pass.snp.flt.bed", chr=CHRS),
        "data/inversion_simple.bed",
        "data/heterochromatin.bed",
    output:
        expand("data/dadi/dadi.{mode}.data", mode=DADI_MODES),
        expand("data/dadi/dadi.chr3.{mode}.data", mode=DADI_MODES)
    threads: 8
    params: runtime="12",
            mem=",mem=48gb"
    shell:
        "sh scripts/make_dadi_input.sh;"

# ----------------------------------------------------------------------------------------
# --- Polarize dadi input files using melas and merus genomes - just chr3
# ----------------------------------------------------------------------------------------

rule polarize_dadi:
    input:
        "data/dadi/dadi.chr3.{mode}.data",
        "data/AGC/mela_ref_ug_vqsr_cnvrt_sort.vcf.gz",
        "data/AGC/meru_ref_ug_vqsr_cnvrt_sort.vcf.gz"
    output:
        "data/dadi_polarized/dadi.chr3.{mode}.data"
    threads: 20
    params: runtime="24",
            mem=",mem=5gb"
    shell:
        "sh scripts/add_dadi_outgroup_all.sh data/dadi/dadi.chr3.{wildcards.mode}.data"

# ----------------------------------------------------------------------------------------
# --- Prep for dadi bootstrapping
# ----------------------------------------------------------------------------------------

rule prep_for_dadi_bootstrap:
    input:
        "data/inversion_simple.bed",
        "data/heterochromatin.bed",
    output:
        "data/random_hunks_for_dadi_bootstrap.bed",
        "data/dadi_bootstraps/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa"
    threads: 8
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/prep_for_bootstrapping.sh;"

# ----------------------------------------------------------------------------------------
# --- Do bootstrapping for dadi (chr3 hunks only)
# ----------------------------------------------------------------------------------------

rule bootstrap_for_dadi:
    input:
        expand("data/chr{chr}.pass.snp.flt.bed", chr=CHRS),
        "data/dadi_polarized/dadi.chr3.{mode}.data",
        "data/random_hunks_for_dadi_bootstrap.bed",
    output:
        ###expand("data/dadi_bootstraps/bs_{iter}.{{mode}}.data",
        ###    iter=range(0,1000))
        "data/dadi_bootstraps/bs_X.{mode}.data.marker",
    threads: 8
    params: runtime="60",
            mem=",mem=24gb"
    shell:
        "sh scripts/create_bootstrapped_vcfs.sh {wildcards.mode};"
        "touch {output}"

# ----------------------------------------------------------------------------------------
# --- Run dadi
# ----------------------------------------------------------------------------------------

rule run_dadi:
    input:
        main_data="data/dadi_polarized/dadi.chr3.{mode}.data",
        #bs=expand("data/dadi_bootstraps/bs_{iter}.{mode}.data",
        #    iter=range(0,1000), mode=DADI_MODES)
        bs="data/dadi_bootstraps/bs_X.{mode}.data.marker"
    output:
        main_out="results/dadi_polarized/dadi.{mode}-{subset}.iter{dadi_iter}.out",
        sfs="data/dadi_polarized/dadi.chr3.{mode}.{subset}.iter{dadi_iter}.sfs.txt"
    threads: 1
    params: runtime="4",
            mem=",mem=96gb"
    shell:
        "mkdir -p results/dadi_polarized/;"
        "module load python/2.7.14-anaconda5.0.1;"
        "export PYTHONPATH=/storage/home/cxb585/local_python;"
        "export PYTHONPATH=$PYTHONPATH:/storage/home/cxb585/local_python//lib64/python2.7/site-packages;"
        "export PYTHONPATH=$PYTHONPATH:/storage/home/cxb585/bin/dadi;"
        "if [ '{wildcards.mode}' == 'island' ]; then NUM_IND=`wc -l data/ssese.seqids.is-{wildcards.subset}.txt | cut -d' ' -f 1`; else NUM_IND=50; fi;"
        #"if [ '{wildcards.subset}' == 'BUGALAIS' ]; then NUM_IND=6; fi;"
        "SAMP=`echo \"2 * $NUM_IND - 2\" | bc`;"
        "python scripts/run_dadi_uncert.py -i {input.main_data} -p {wildcards.subset} -s $SAMP -m 1000 -f False -e .iter{wildcards.dadi_iter} > {output.main_out};"
        #"python scripts/run_dadi_uncert.py -i {input.main_data} -p {wildcards.subset} -s $SAMP -m 1000 -f True  -e .iter{wildcards.dadi_iter} > results/dadi_polarized/dadi.{wildcards.mode}-{wildcards.subset}.folded.iter{dadi_iter}.out;"

# ----------------------------------------------------------------------------------------
# --- Ag100G: Make dadi input with ONLY Ag1000G data (not combined)
# ----------------------------------------------------------------------------------------

rule make_dadi_input_ag1000g:
    input:
        expand("data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.pass.biallelic.{chr}.vcf.gz", chr=CHR_THREE),
        "data/inversion_simple.bed",
        "data/heterochromatin.bed",
    output:
        "data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.pass.biallelic.all.for.dadi.vcf",
        "data/dadi/dadi_pop_list.ag1000g.txt",
        "data/dadi/dadi.ag1000g.data"
    threads: 8
    params: runtime="12",
            mem=",mem=48gb"
    shell:
        "sh scripts/prep_to_make_dadi_input_with_ag1000g.sh"

# ----------------------------------------------------------------------------------------
# --- Polarize Ag1000G dadi input files using melas and merus genomes - just chr3
# ----------------------------------------------------------------------------------------

rule polarize_dadi_ag1000g:
    input:
        dadi="data/dadi/dadi.ag1000g.data",
        mela="data/AGC/mela_ref_ug_vqsr_cnvrt_sort.vcf.gz",
        meru="data/AGC/meru_ref_ug_vqsr_cnvrt_sort.vcf.gz"
    output:
        "data/dadi_polarized/dadi.ag1000g.data"
    threads: 20
    params: runtime="24",
            mem=",mem=5gb"
    shell:
        "sh scripts/add_dadi_outgroup_all.sh {input.dadi}"

# ----------------------------------------------------------------------------------------
# --- Ag100G: Do bootstrapping for dadi with Ag1000G - Ag1000G only
# ----------------------------------------------------------------------------------------

rule bootstrap_for_dadi_ag:
    input:
        "data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.pass.biallelic.all.for.dadi.vcf",
        "data/dadi_polarized/dadi.ag1000g.data",
        "data/random_hunks_for_dadi_bootstrap.bed",
    output:
        ###expand("data/dadi_bootstraps/bs_{iter}.ag1000g.data",
        ###    iter=range(0,1000))
        "data/dadi_bootstraps/bs_X.ag1000g.data.marker"
    threads: 8
    params: runtime="60",
            mem=",mem=24gb"
    shell:
        "sh scripts/create_bootstrapped_vcfs_ag1000g.sh;"
        "touch {output}"

# ----------------------------------------------------------------------------------------
# --- Ag100G: Run dadi - Ag1000G only
# ----------------------------------------------------------------------------------------

rule run_dadi_with_ag1000g:
    input:
        dadi="data/dadi_polarized/dadi.ag1000g.data",
        ###bs=expand("data/dadi_bootstraps/bs_{iter}.ag1000g.data",
        ###    iter=range(0,1000))
        bs="data/dadi_bootstraps/bs_X.ag1000g.data.marker"
    output:
        dadi_out="results/dadi_polarized/dadi.ag1000g.{subset}.iter{dadi_iter}.out",
        sfs_out="data/dadi_polarized/dadi.ag1000g.{subset}.iter{dadi_iter}.sfs.txt"
    threads: 8
    params: runtime="24",
            mem=",mem=24gb"
    shell:
        "mkdir -p results/dadi/;"
        "module load python/2.7.14-anaconda5.0.1;"
        "export PYTHONPATH=/storage/home/cxb585/local_python;"
        "export PYTHONPATH=$PYTHONPATH:/storage/home/cxb585/bin/dadi:/usr/global/python/2.7.8/lib/python2.7/site-packages;"
        "python scripts/run_dadi_uncert.py -i {input.dadi} -p {wildcards.subset} -s 104 -m 1000 -f False -e .iter{wildcards.dadi_iter} > {output.dadi_out};"
        #"python scripts/run_dadi_uncert.py -i {input.dadi} -p {wildcards.subset} -s 104 -m 1000 -f True  -e .iter{wildcards.dadi_iter} > results/dadi_polarized/dadi.ag1000g.{wildcards.subset}.iter{dadi_iter}.out;"

# ----------------------------------------------------------------------------------------
# --- Run dadi - 2 pop model between sites
# ----------------------------------------------------------------------------------------

rule run_dadi_2pop_intersite:
    input:
        dadi="data/dadi_polarized/dadi.chr3.island.data",
        bs="data/dadi_bootstraps/bs_X.island.data.marker"
    output:
        "results/dadi_polarized/dadi.island.2pop.{site_pair}.iter{dadi_iter}.out"
    threads: 1
    params: runtime="96",
            mem=",mem=96gb"
    shell:
        "module load python/2.7.14-anaconda5.0.1;"
        "export PYTHONPATH=/storage/home/cxb585/local_python;"
        "export PYTHONPATH=$PYTHONPATH:/storage/home/cxb585/bin/dadi;"
        "export PYTHONPATH=$PYTHONPATH:/usr/global/python/2.7.8/lib/python2.7/site-packages;"
        "SITE_PAIR={wildcards.site_pair};"
        "SITE_A=`echo $SITE_PAIR | cut -d'.' -f 1`;"
        "SITE_B=`echo $SITE_PAIR | cut -d'.' -f 2`;"
        "NUM_IND_A=`wc -l data/ssese.seqids.is-$SITE_A.txt | cut -d' ' -f 1`;"
        "NUM_IND_B=`wc -l data/ssese.seqids.is-$SITE_B.txt | cut -d' ' -f 1`;"
        #"if [ \"$SITE_A\" == 'BUGALAIS' ]; then NUM_IND_A=6; fi;"
        #"if [ \"$SITE_B\" == 'BUGALAIS' ]; then NUM_IND_B=6; fi;"
        "SAMP_A=`echo \"2 * $NUM_IND_A - 6\" | bc`;"
        "SAMP_B=`echo \"2 * $NUM_IND_B - 6\" | bc`;"
        "python scripts/run_dadi_2pop_IM.py -i {input.dadi} -p $SITE_A -q $SITE_B -s $SAMP_A -z $SAMP_B -m 1000 -e .iter{wildcards.dadi_iter} > {output};"

# ----------------------------------------------------------------------------------------
# --- Fix dadi results (1- and 2-pop)
# ----------------------------------------------------------------------------------------

rule fix_dadi_results:
    input:
        expand("results/dadi_polarized/dadi.{mode_subset}.iter{dadi_iter}.out",
            mode_subset=DADI_OUTPUT,
            dadi_iter=range(0,9)),
        expand("results/dadi_polarized/dadi.island.2pop.{site_pair}.iter{dadi_iter}.out",
            site_pair=DADI_SITE_PAIRS,
            dadi_iter=range(0,19))
    output:
        expand("results/dadi_polarized/dadi.{mode_subset}.iter{dadi_iter}.out.fix",
            mode_subset=DADI_OUTPUT,
            dadi_iter=range(0,9)),
        expand("results/dadi_polarized/dadi.island.2pop.{site_pair}.iter{dadi_iter}.out.fix",
            site_pair=DADI_SITE_PAIRS,
            dadi_iter=range(0,19))
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "for file in results/dadi_polarized/*out; do "
        "    grep \"^[id]\" $file > $file.fix;"
        "done;"

# ----------------------------------------------------------------------------------------
# --- Plot and parse dadi results
# ----------------------------------------------------------------------------------------

rule parse_dadi:
    input:
        expand("results/dadi_polarized/dadi.{mode_subset}.iter{dadi_iter}.out.fix",
            mode_subset=DADI_OUTPUT,
            dadi_iter=range(0,9)),
        expand("results/dadi_polarized/dadi.island.2pop.{site_pair}.iter{dadi_iter}.out.fix",
            site_pair=DADI_SITE_PAIRS,
            dadi_iter=range(0,9)),
        "data/inversion_simple.bed"
    output:
        "reports/dadi_three_epoch_pop_history.pdf",
        "reports/dadi.island.out.tex",
        "reports/dadi_island_is-ml_differences.txt",
        "reports/dadi_2pop_is-ml_differences.txt",
        "reports/dadi.2pop.island.island.out.tex",
        "reports/dadi.2pop.island.mainland.out.tex",
        "reports/dadi.2pop.mainland.mainland.out.tex",
        "reports/dadi_migration_matrix.txt",
        "reports/dadi_top_migration.txt",
        "reports/dadi_migration_matrix.pdf",
        "reports/dadi_2pop_is-ml_differences.time-size.txt",
        "results/dadi.best_iterations.txt",
        "results/dadi.best_iterations.2pop.txt",
        "results/asymmetrical.migration.txt",
    threads: 1
    params: runtime="2",
            mem=",mem=5gb"
    shell:
        "Rscript scripts/plot_parse_dadi_results.R"

# ----------------------------------------------------------------------------------------
# --- Copy in best iteration's optimization plots
# ----------------------------------------------------------------------------------------

rule copy_dadi_optimization_plots:
    input:
        "results/dadi.best_iterations.txt",
        "results/dadi.best_iterations.2pop.txt",
    output:
        expand("data/dadi_polarized/dadi.chr3.island.{site}.bestmodel.best-iter.pdf",
            site=ISLANDS_SPLIT_BUGALA),
        expand("data/dadi_polarized/dadi.chr3.island.{site_pair}.best-iter.pdf",
            site_pair=DADI_SITE_PAIRS)
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sh scripts/copy_dadi_optimization_plots.sh"

# ----------------------------------------------------------------------------------------
# --- Explore theta along chromosome 3
# ----------------------------------------------------------------------------------------

rule explore_theta:
    input:
        "results/dadi.best_iterations.txt",
        "data/dadi_polarized/dadi.chr3.island.data",
    output:
        "results/dadi.chr3.island.{site}.theta_along_chr3.txt",
    threads: 1
    params: runtime="24",
            mem=",mem=96gb"
    shell:
        "module load python/2.7.14-anaconda5.0.1; "
        "export PYTHONPATH=/storage/home/cxb585/local_python; "
        "export PYTHONPATH=$PYTHONPATH:/storage/home/cxb585/local_python//lib64/python2.7/site-packages; "
        "export PYTHONPATH=$PYTHONPATH:/storage/home/cxb585/bin/dadi; "
        "POP={wildcards.site}; "
        "ITER=`grep $POP results/dadi.best_iterations.txt | cut -d'.' -f3`; "
        "python scripts/explore_theta_along_chr3.py -i data/dadi_polarized/dadi.chr3.island.data "
        "    -p $POP -s 10 -m 100 -f False -e'.$ITER'"

# ========================================================================================
# --- Run ADMIXTURE
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Create thinned dataset for ADMIXTURE runs (just chr3) (10,000 SNPs only)
# ----------------------------------------------------------------------------------------

rule create_thinned_adm_beds:
    input:
        beds=expand("data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.bed", chr=CHR_THREE),
        exclude="data/ag1000g.phase1.ar3/excluded.individuals.txt",
        exclude_not_M="data/ag1000g.phase1.ar3/excluded.individuals.notM.txt",
    output:
        expand("data/ssese_with_ag1000g/adm_subsets/chr3.replicate{rep}.LD.bed",
            rep=range(1,11)),
        expand("data/ssese_with_ag1000g/adm_subsets/chr3.withM.replicate{rep}.LD.bed",
            rep=range(1,11))
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/generate_admixture_replicate_datasets.sh;"
        "rm data/ssese_with_ag1000g/adm_subsets/*seed*;"
        "rm data/ssese_with_ag1000g/adm_subsets/chr3.*log;"
        "rm data/ssese_with_ag1000g/adm_subsets/chr3.*nosex;"
        "rm data/ssese_with_ag1000g/adm_subsets/chr3*full*;"

rule create_thinned_adm_beds_missingtoref:
    input:
        beds=expand("data/ssese_with_ag1000g.missingtoref/chr{chr}.missingtoref.pass.snp.phased.ag1000g.strict.bed", chr=CHR_THREE),
        exclude="data/ag1000g.phase1.ar3/excluded.individuals.txt",
        exclude_not_M="data/ag1000g.phase1.ar3/excluded.individuals.notM.txt",
    output:
        expand("data/ssese_with_ag1000g.missingtoref/adm_subsets/chr3.replicate{rep}.LD.bed",
            rep=range(1,11)),
        expand("data/ssese_with_ag1000g.missingtoref/adm_subsets/chr3.withM.replicate{rep}.LD.bed",
            rep=range(1,11))
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/generate_admixture_replicate_datasets.sh missingtoref;"
        "rm data/ssese_with_ag1000g.missingtoref/adm_subsets/*seed*;"
        "rm data/ssese_with_ag1000g.missingtoref/adm_subsets/chr3.*log;"
        "rm data/ssese_with_ag1000g.missingtoref/adm_subsets/chr3.*nosex;"
        "rm data/ssese_with_ag1000g.missingtoref/adm_subsets/chr3*full*;"

# ----------------------------------------------------------------------------------------
# --- Run ADMIXTURE (in CV mode) for given k and given rep input BED, and do it 10 times
# ----------------------------------------------------------------------------------------

rule adm_cv:
    input:
        adm_bed="data/ssese_with_ag1000g/adm_subsets/chr3{M_part}replicate{rep}.LD.bed"
    output:
        adm_out="data/ssese_with_ag1000g/adm_subsets/chr3{M_part}replicate{rep}.LD.iter{iter}.{k}.ADMIXTURE_log{k}.out",
        adm_p="data/ssese_with_ag1000g/adm_subsets/chr3{M_part}replicate{rep}.LD.iter{iter}.{k}.P",
        adm_q="data/ssese_with_ag1000g/adm_subsets/chr3{M_part}replicate{rep}.LD.iter{iter}.{k}.Q",
    threads: 8
    params: runtime="2",
            mem=",mem=5gb"
    shell:
        # Make copy of input BED to avoid overwriting other output files
        "IN_BED={input.adm_bed};"
        "BED_PREFIX=${{IN_BED/.bed/}}.iter{wildcards.iter}.{wildcards.k};"
        "cp $IN_BED $BED_PREFIX.bed;"
        "cp ${{IN_BED/.bed/.bim}} $BED_PREFIX.bim;"
        "cp ${{IN_BED/.bed/.fam}} $BED_PREFIX.fam;"
        # Call ADMIXTURE in CV mode
        "sh scripts/run_admixture_cv.sh $BED_PREFIX.bed {wildcards.k} 8;"
        "rm $BED_PREFIX.bed;"
        "rm $BED_PREFIX.bim;"
        "rm $BED_PREFIX.fam;"
        "mv `basename $BED_PREFIX`.{wildcards.k}.P {output.adm_p};"
        "mv `basename $BED_PREFIX`.{wildcards.k}.Q {output.adm_q};"

# ----------------------------------------------------------------------------------------
# --- Plot ADMIXTURE CV error
# ----------------------------------------------------------------------------------------

rule plot_adm_cv:
    input:
        expand("data/ssese_with_ag1000g/adm_subsets/chr3{{M_part}}replicate{rep}.LD.iter{iter}.{k}.ADMIXTURE_log{k}.out",
            rep=range(1,11),
            iter=range(1,6),
            k=range(2,11))
    output:
        "data/ssese_with_ag1000g/adm_subsets/chr3{M_part}CV.pdf"
    threads: 1
    params: runtime="2",
            mem=",mem=5gb"
    shell:
        "PREFIX=data/ssese_with_ag1000g/adm_subsets/chr3{wildcards.M_part};"
        "grep \"^CV error\" ${{PREFIX}}replicate*.LD.iter*.*.ADMIXTURE_log*.out | "
        "cut -d':' -f2-3 | sed -e \"s/CV error (K=//\" -e \"s/): /\t/\" > ${{PREFIX}}CV.txt;"
        "Rscript scripts/plot_admixture_CV.R ${{PREFIX}}CV.txt;"

# ----------------------------------------------------------------------------------------
# --- Make CLUMPAK input
# ----------------------------------------------------------------------------------------

rule make_clumpak_input:
    input:
        expand("data/ssese_with_ag1000g/adm_subsets/chr3{M_part}replicate{rep}.LD.iter{iter}.{k}.Q",
            M_part=['.withM.'], rep=range(1,11), k=range(2,11), iter=range(1,6))
    output:
        #"data/ssese_with_ag1000g/adm_subsets/chr3.allQs.zip",
        "data/ssese_with_ag1000g/adm_subsets/chr3.withM.allQs.zip"
    threads: 8
    params: runtime="2",
            mem=",mem=5gb"
    shell:
        "sh scripts/make_clumpak_input.sh"

# ----------------------------------------------------------------------------------------
# --- Make CLUMPAK population file
# ----------------------------------------------------------------------------------------

rule make_clumpak_pop_file:
    input:
#        "data/ssese_with_ag1000g/adm_subsets/chr3.replicate1.LD.bed",
        "data/ssese_with_ag1000g/adm_subsets/chr3.withM.replicate1.LD.bed",
        "data/ssese_individual_info_simple_bugala_split.txt",
        "data/sample_to_seq_id_mapping.txt",
        "data/ag1000g.phase1.ar3/samples.all.txt"
    output:
#        "data/ssese_with_ag1000g/adm_subsets/pop.list.txt",
        "data/ssese_with_ag1000g/adm_subsets/pop.list.withM.txt"
    threads: 1
    params: runtime="2",
            mem=",mem=5gb"
    shell:
        "Rscript scripts/make_clumpak_pop_list.R"

# ----------------------------------------------------------------------------------------
# --- Make Best K input
# ----------------------------------------------------------------------------------------

rule make_bestk_input:
    input:
        expand("data/ssese_with_ag1000g/adm_subsets/chr3{M_part}replicate{rep}.LD.iter{iter}.{k}.ADMIXTURE_log{k}.out",
            M_part=['.withM.'], rep=range(1,11), k=range(2,11), iter=range(1,6))
    output:
#        "data/ssese_with_ag1000g/adm_subsets/bestK.input.log_prob.txt",
        "data/ssese_with_ag1000g/adm_subsets/bestK.input.withM.log_prob.txt"
    threads: 8
    params: runtime="2",
            mem=",mem=5gb"
    shell:
        "sh scripts/make_bestK_input.sh"

# ----------------------------------------------------------------------------------------
# --- Parse CLUMPAK results
# ----------------------------------------------------------------------------------------

# At this point, CLUMPAK results must be brought in and renamed to lose "=" in names

rule warn_abt_clumpak_results:
    input:
#        "data/ssese_with_ag1000g/adm_subsets/chr3.allQs.zip",
        "data/ssese_with_ag1000g/adm_subsets/chr3.withM.allQs.zip",
#        "data/ssese_with_ag1000g/adm_subsets/pop.list.txt",
        "data/ssese_with_ag1000g/adm_subsets/pop.list.withM.txt"
    output:
        expand("data/ssese_with_ag1000g/adm_subsets/CLUMPAK_{m_flag}/K{k}/MajorCluster/CLUMPP.files/ClumppIndFile.output",
            m_flag=['withM'], k=range(2,11))
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "echo 'Submit ADM to CLUMPAK website now. Be sure to rename to lose equals sign in names.'"

rule parse_clumpak_results:
    input:
        expand("data/ssese_with_ag1000g/adm_subsets/CLUMPAK_{m_flag}/K{k}/MajorCluster/CLUMPP.files/ClumppIndFile.output",
            m_flag=['withM'], k=range(2,11))
    output:
        expand("data/ssese_with_ag1000g/adm_subsets/CLUMPAK_{m_flag}.{k}.Q",
            m_flag=['withM'], k=range(2,11))
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "sh scripts/parse_clumpak_results.sh"

# ----------------------------------------------------------------------------------------
# --- Plot ADMIXTURE results
# ----------------------------------------------------------------------------------------

rule plot_adm_M:
    input:
        adm_bed="data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.thinned_for_ADMIXTURE.withM.bed",
        #adm_p="data/ssese_with_ag1000g/chr3.withM.pass.snp.phased.ag1000g.strict.2.P",
        #adm_q="data/ssese_with_ag1000g/chr3.withM.pass.snp.phased.ag1000g.strict.2.Q",
        adm_q="data/ssese_with_ag1000g/adm_subsets/CLUMPAK_withM.2.Q",
        #all_adm_q=expand("data/ssese_with_ag1000g/chr3.withM.pass.snp.phased.ag1000g.strict.{k}.Q", k=range(2,11))
    output:
        expand("data/ssese_with_ag1000g/chr3.withM.pass.snp.phased.ag1000g.strict.{k}.pdf", k=range(2,11)),
        expand("data/ssese_with_ag1000g/chr3.withM.pass.snp.phased.ag1000g.strict.LVB_only.{k}.pdf", k=range(2,11)),
        expand("data/ssese_with_ag1000g/chr3.withM.pass.snp.phased.ag1000g.strict.{k}.boot.fortext.pdf", k=range(2,11)),
        expand("data/ssese_with_ag1000g/chr3.withM.pass.snp.phased.ag1000g.strict.LVB_only.{k}.boot.fortext.pdf", k=range(2,11)),
        expand("data/ssese_with_ag1000g/adm_subsets/CLUMPAK_withM.LVB_only.{k}.boot.pdf", k=range(2,11)),
    threads: 1
    params: runtime="2",
            mem=",mem=5gb"
    shell:
        "BED={input.adm_bed};"
        "FAM=${{BED/.bed/.fam}};"
        "ADM_Q={input.adm_q};"
        "PREFIX=${{ADM_Q/2.Q/}};"
        "Rscript scripts/plot_admixture_results.R $PREFIX $FAM;"
        "for k in `seq 2 10`; do "
        "    mv $PREFIX$k.pdf data/ssese_with_ag1000g/chr3.withM.pass.snp.phased.ag1000g.strict.$k.pdf;"
        "    mv ${{PREFIX}}LVB_only.$k.pdf data/ssese_with_ag1000g/chr3.withM.pass.snp.phased.ag1000g.strict.LVB_only.$k.pdf;"
        "    mv ${{PREFIX}}$k.boot.fortext.pdf data/ssese_with_ag1000g/chr3.withM.pass.snp.phased.ag1000g.strict.$k.boot.fortext.pdf;"
        "done;"

# ========================================================================================
# --- Look at Runs of Homozygosity
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Identify Runs of Homozygosity
# ----------------------------------------------------------------------------------------

rule find_ROHs:
    input:
        beds=expand("data/chr{chr}.pass.snp.flt.bed", chr=CHRS),
        prune="data/all.pass.snp.flt.prune.out",
        inv="data/inversion_simple.bed",
        het="data/heterochromatin.bed"
    output:
        "results/ROH/all.pass.snp.flt.noinv.hom",
        "results/ROH/all.pass.snp.flt.noinv.hom.indiv",
        "results/ROH/all.pass.snp.flt.noinv.hom.overlap",
        "results/ROH/all.pass.snp.flt.noinv.hom.summary",
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/identify_ROH.sh"

# ----------------------------------------------------------------------------------------
# --- Plot Runs of Homozygosity
# ----------------------------------------------------------------------------------------

rule plot_ROHs:
    input:
        "results/ROH/all.pass.snp.flt.noinv.hom",
        "results/ROH/all.pass.snp.flt.noinv.hom.indiv",
        "data/sample_to_seq_id_mapping.txt",
    output:
        "reports/roh_plot.pdf",
        "reports/roh_boxplot.pdf"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_ROH.R"

# ========================================================================================
# --- Look at shared IBD tracts
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Identify IBD
# ----------------------------------------------------------------------------------------

rule identify_IBD:
    input:
        beds=expand("data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.strict.bed", chr=CHRS),
        prune="data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.prune.out",
        exclude="data/ag1000g.phase1.ar3/excluded.individuals.txt",
        inv="data/inversion_simple.bed",
        het="data/heterochromatin.bed"
    output:
        "results/IBD/all.pass.snp.flt.noinv.genome",
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/find_shared_IBD.sh"

# ----------------------------------------------------------------------------------------
# --- Plot IBD results
# ----------------------------------------------------------------------------------------

rule plot_IBD:
    input:
        "results/IBD/all.pass.snp.flt.noinv.genome",
    output:
        "reports/ibd.pdf"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_ibd_results.R"

# ========================================================================================
# --- Plot EHH decay
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Plot EHH decay for new X and 2L sweeps
# ----------------------------------------------------------------------------------------

rule plot_EHH_ROI:
    input:
        "data/chrX.pass.snp.phased.haps",
        "data/chr2L.pass.snp.phased.haps",
        "data/chr2R.pass.snp.phased.haps",
        "data/ssese_individual_info_simple.txt"
    output:
        "results/ehh/chrX.pass.snp.phased.figs/ehh_chrX.X_9Mb_1.pdf",
        "results/ehh/chrX.pass.snp.phased.figs/ehh_chrX.X_4Mb_1.pdf",
        "results/ehh/chr2L.pass.snp.phased.figs/ehh_chr2L.2L_34Mb_1.pdf",
        "results/ehh/chr2R.pass.snp.phased.figs/ehh_chr2R.CYP6P2_1.pdf",
        "results/ehh/chrX.pass.snp.phased.figs/ehh_chrX.CYP9K1_1.pdf",
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        #"module load r/3.4;"
        # 2L_34Mb
        "CENTER=`grep '2L_near34Mb' results/xpehh_peaks.txt | cut -d':' -f2`;"
        "Rscript scripts/compute_ehh_and_friends.R 2L $CENTER 100000 2L_34Mb FALSE;"
        # X_9Mb
        "CENTER=`grep 'X_near9Mb' results/xpehh_peaks.txt | cut -d':' -f2`;"
        "Rscript scripts/compute_ehh_and_friends.R X  $CENTER 100000 X_9Mb TRUE;"
        # X_4Mb
        "CENTER=4000000;"
        "Rscript scripts/compute_ehh_and_friends.R X  $CENTER 100000 X_4Mb FALSE;"
        # CYP6P2
        "CENTER=28501972;"
        "Rscript scripts/compute_ehh_and_friends.R 2R $CENTER 100000 CYP6P2 FALSE;"
        # CYP9K1
        "CENTER=15241718;"
        "Rscript scripts/compute_ehh_and_friends.R X  $CENTER 100000 CYP9K1 FALSE;"

# ========================================================================================
# --- Annotate SNPs with snpEff
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Run annotation
# ----------------------------------------------------------------------------------------

rule annotate_snps:
    input:
        "data/chr{chr}.pass.snp.vcf"
    output:
        "data/chr{chr}.pass.snp.eff.vcf"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/call_snpEff.sh {wildcards.chr}"

# ========================================================================================
# --- Infer trees (dendrograms) around regions of interest (ROIs)
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Infer trees
# ----------------------------------------------------------------------------------------

rule infer_dendrograms_roi:
    input:
        expand("data/ssese_with_ag1000g/chr{chr}.pass.snp.phased.ag1000g.haps", chr=CHRS),
        expand("data/chr{chr}.pass.snp.phased.haps", chr=CHRS),
        "results/xpehh_peaks.txt"
    output:
        expand("results/for_dendrograms/{locus}.{window}.{ploidy}.with-ag1000g.sans-snpeff.haps.{ending}",
            locus = ['chr2L.2L_34Mb', 'chrX.X_9Mb', 'chrX.X_4Mb', 'chrX.X_CYP9K1', 'chr2R.2R_CYP6P2', 'chr3R.3R_GSTE'],
            window = ['10000', '1e+05'],
            ploidy = ['haploid'],
            ending = ['fancydendro.pdf', 'haplotypes.txt'])
    threads: 12
    params: runtime="24",
            mem=",mem=5gb"
    shell:
        "sh scripts/explore_all_ROI.sh"

# ----------------------------------------------------------------------------------------
# --- Plot sampling map and map of haplotypes (map o' haps)
# ----------------------------------------------------------------------------------------

rule plot_map_o_haps:
    input:
        "data/site_gps_coordinates.txt",
        "data/ssese_individual_info_simple.txt",
        expand("results/for_dendrograms/{locus}.1e+05.haploid.with-ag1000g.sans-snpeff.haps.haplotypes.txt",
            locus = ['chr2L.2L_34Mb', 'chrX.X_9Mb'])
    output:
        "results/ssese_sampling_sites.png",
        "results/ssese_haplotypes_sweep_2L_34Mb.png",
        "results/ssese_haplotypes_sweep_X_9Mb.png",
        expand("data/ssese_with_ag1000g/adm_subsets/CLUMPAK_withM.{k}.Q.map.png",
            k=range(2,11)),
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        #"module load r/3.4;"
        "Rscript scripts/plot_map_of_haplotypes.R"

# ========================================================================================
# --- Plot pop gen stats around regions of interest (ROIs)
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Plot pop gen around putative sweeps
# ----------------------------------------------------------------------------------------

rule make_zoomed_stats_plot:
    input:
        expand("results/chr{chr}.{site_pair}.windowed.weir.fst",
            chr=CHRS, site_pair=SITE_PAIRS),
        expand("results/selscan/xp-ehh.{site_pair}.{chr}.xpehh.out.norm",
            chr=CHRS, site_pair=map(lambda x: x.replace("_", "-"), SITE_PAIRS)),
        expand("results/chr{chr}.pass.snp.phased.{site}.H12.out.txt",
            chr=CHRS, site=ISLANDS_SPLIT_BUGALA),
        expand("results/chr{chr}.{site}.windowed.pi",
            chr=CHRS, site=ISLANDS_SPLIT_BUGALA),
        expand("results/chr{chr}.{site}.Tajima.D",
            chr=CHRS, site=ISLANDS_SPLIT_BUGALA),
        "results/xpehh_peaks.txt"
    output:
        expand("reports/all_stats.{locus}.{{site}}.pdf",
            locus = ['2L_34Mb', 'X_9Mb', 'X_4Mb', 'CYP6P2', 'CYP9K1', 'GSTE'])
    threads: 12
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        #"module load r/3.4;"
        "module load bedtools;"
        # 2L_34Mb
        "CENTER=`grep '2L_near34Mb' results/xpehh_peaks.txt | cut -d':' -f2`;"
        "LEFT=$((CENTER - 1000000));"
        "RIGHT=$((CENTER + 1000000));"
        "Rscript scripts/plot_stats_around_ROI.R 2L $LEFT $RIGHT {wildcards.site} 2L_34Mb;"
        # X_9Mb
        "CENTER=`grep 'X_near9Mb' results/xpehh_peaks.txt | cut -d':' -f2`;"
        "LEFT=$((CENTER - 1000000));"
        "RIGHT=$((CENTER + 1000000));"
        "Rscript scripts/plot_stats_around_ROI.R X $LEFT $RIGHT {wildcards.site} X_9Mb;"
        # X_4Mb
        "CENTER=4000000;"
        "LEFT=3000000;"
        "RIGHT=5000000;"
        "Rscript scripts/plot_stats_around_ROI.R X $LEFT $RIGHT {wildcards.site} X_4Mb;"
        # CYP6P2[AGAP002869]
        "CENTER=28501972;"
        "LEFT=$((CENTER - 1000000));"
        "RIGHT=$((CENTER + 1000000));"
        "Rscript scripts/plot_stats_around_ROI.R 2R $LEFT $RIGHT {wildcards.site} CYP6P2;"
        # CYP9K1
        "CENTER=15241718;"
        "LEFT=$((CENTER - 1000000));"
        "RIGHT=$((CENTER + 1000000));"
        "Rscript scripts/plot_stats_around_ROI.R X $LEFT $RIGHT {wildcards.site} CYP9K1;"
        # GSTE
        "CENTER=28598038;"
        "LEFT=$((CENTER - 1000000));"
        "RIGHT=$((CENTER + 1000000));"
        "Rscript scripts/plot_stats_around_ROI.R 3R $LEFT $RIGHT {wildcards.site} GSTE;"

# ========================================================================================
# --- Do stairway plot
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Do stairway plot for site
# ----------------------------------------------------------------------------------------

rule stairway_plot_site:
    input:
        expand("data/dadi_polarized/dadi.chr3.island.{site}.sfs.txt", site=ISLANDS_SPLIT_BUGALA),
        "data/ssese_individual_info_simple_bugala_split.txt"
    output:
        "data/stairway_plot_island.{site}/island.{site}.final.summary.pdf"
    threads: 12
    params: runtime="48",
            mem=",mem=8gb"
    shell:
        "sh scripts/do_stairway_plot.sh island.{wildcards.site}"

# ----------------------------------------------------------------------------------------
# --- Make plot of stairway... plots
# ----------------------------------------------------------------------------------------

rule plot_stairway_plots:
    input:
        expand("data/stairway_plot_island.{site}/island.{site}.final.summary.pdf",
            site=ISLANDS_SPLIT_BUGALA)
    output:
        "reports/stairway.plot.superimposed.pdf"
    threads: 1
    params: runtime="4",
            mem=",mem=8gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_stairway_plot_results.R"

# ========================================================================================
# --- Estimate Ne
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Estimate Ne with NeEstimator
# ----------------------------------------------------------------------------------------

rule estimate_Ne:
    input:
        "data/chr3L.pass.snp.flt.vcf.gz",
        "data/inversion_heterochromatin.bed",
        expand("data/ssese.seqids.is-{islands}.txt", islands=ISLANDS),
    output:
        "results/genepop_Ne_estimates.txt",
        "results/genepop_bounds.txt"
    threads: 1
    params: runtime="4",
            mem=",mem=8gb"
    shell:
        "sh scripts/run_NeEstimator.sh"

# ----------------------------------------------------------------------------------------
# --- Estimate Ne with IBDNe - first run IBDseq
# ----------------------------------------------------------------------------------------

rule run_ibdseq:
    input:
        "data/chr{chr}.pass.snp.flt.vcf.gz",
        "data/chr{chr}.pass.snp.phased.haps"
    output:
        "results/chr{chr}.pass.snp.flt.ibdseq.ibd"
    threads: 20
    params: runtime="4",
            mem=",mem=8gb"
    shell:
        "module load ibdseq;"
        "sh scripts/call_ibdseq.sh {wildcards.chr}"

# ----------------------------------------------------------------------------------------
# --- Estimate Ne with IBDNe
# ----------------------------------------------------------------------------------------

rule run_ibdne:
    input:
        "results/chr{chr}.pass.snp.flt.ibdseq.ibd",
        "data/ssese_individual_info_simple_bugala_split.txt",
        "data/chr{chr}.pass.snp.phased.haps",
    output:
        "results/chr{chr}.pass.snp.{site}.ne.ne"
    threads: 8
    params: runtime="48",
            mem=",mem=24gb"
    shell:
        "module load ibdne;"
        "module load r/3.4;"
        "Rscript scripts/call_ibdne.R {wildcards.site} {wildcards.chr}"

# ----------------------------------------------------------------------------------------
# --- Plot IBDNe results
# ----------------------------------------------------------------------------------------

rule plot_ibdne:
    input:
        expand("results/chr{{chr}}.pass.snp.{site}.ne.ne",
            site=ISLANDS_SPLIT_BUGALA),
    output:
        "reports/IBDNe_by_pop.{chr}.pdf",
        "reports/IBDNe_by_site_type.{chr}.pdf",
        "reports/IBDNe_by_site_type_Mityana_highlighted.{chr}.pdf",
    threads: 1
    params: runtime="1",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/plot_ibdne_res.R {wildcards.chr}"

# ----------------------------------------------------------------------------------------
# --- Infer demography from IBS (Harris and Nielsen) - Prepare
# ----------------------------------------------------------------------------------------

rule ibs_demog_prep:
    input:
        "data/chr3L.pass.snp.phased.haps.vcf",
        expand("data/ssese.seqids.is-{site}.txt", site=ISLANDS_SPLIT_BUGALA)
    output:
        expand("results/IBS.{site}.popdata", site=ISLANDS_SPLIT_BUGALA)
    threads: 12
    params: runtime="8",
            mem=",mem=5gb"
    shell:
        "sh scripts/infer_demography_from_ibs_prepare.sh"

# ----------------------------------------------------------------------------------------
# --- Infer demography from IBS (Harris and Nielsen) - 1 pop mode
# ----------------------------------------------------------------------------------------

rule ibs_demog_1pop:
    input:
        "results/IBS.{site}.popdata"
    output:
        "results/IBS_inferred_size_{site}.demographic_history.txt",
        "results/IBS_inferred_size_{site}.demographic_history.txt.pdf"
    threads: 1
    params: runtime="48",
            mem=",mem=24gb"
    shell:
        "sh scripts/infer_demography_from_ibs_1pop.sh {wildcards.site}"

# ----------------------------------------------------------------------------------------
# --- Infer demography from IBS (Harris and Nielsen) - 2 pop mode
# ----------------------------------------------------------------------------------------

rule ibs_demog_2pop:
    input:
        lambda wildcards: IBS_INPUT["results/IBS_inferred_size_" + wildcards.site_pair_vs + ".demographic_history.txt"]
    output:
        "results/IBS_inferred_size_{site_pair_vs}.demographic_history.txt"
    threads: 1
    params: runtime="48",
            mem=",mem=48gb"
    shell:
        "BOTH_POPS=`echo {wildcards.site_pair_vs} | sed -e 's/_vs_/ /'`;"
        "sh scripts/infer_demography_from_ibs_2pop.sh $BOTH_POPS"

# ========================================================================================
# --- Make extra tables for paper
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Compute coverage (and other basic stats) by individual
# ----------------------------------------------------------------------------------------

rule compute_coverage:
    input:
        expand("data/chr{chr}.pass.snp.flt.vcf.gz", chr=CHRS),
        "data/inversion_heterochromatin.bed"
    output:
        "reports/all.pass.snp.flt.idepth"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/compute_coverage.sh"

# ----------------------------------------------------------------------------------------
# --- Count total SNPs in dataset
# ----------------------------------------------------------------------------------------

rule count_snps:
    input:
        expand("data/chr{chr}.pass.snp.flt.vcf.gz", chr=CHRS)
    output:
        "reports/SNP_count.txt",
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "sh scripts/count_SNPs.sh"

# ----------------------------------------------------------------------------------------
# --- Make individual info LaTeX table
# ----------------------------------------------------------------------------------------

rule make_ind_table:
    input:
        "reports/all.pass.snp.flt.idepth"
    output:
        "reports/individual_table.tex",
        "reports/site_gps_table.tex"
    threads: 1
    params: runtime="4",
            mem=",mem=5gb"
    shell:
        "module load r/3.4;"
        "Rscript scripts/make_individual_table.R"

# ========================================================================================
# --- Set default for optional mem parameter
# ========================================================================================

default = ""
name = 'mem'
for r in workflow.rules:
    try:
        getattr(r.params, name)
    except AttributeError:
        r.params.append(default)
        r.params.add_name(name)
