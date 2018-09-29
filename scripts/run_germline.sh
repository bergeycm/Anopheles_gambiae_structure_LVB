
# ========================================================================================
# --- Find regions in IBD with GERMLINE
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Convert from phased VCF to BEAGLE format
# ----------------------------------------------------------------------------------------

IN_VCF=data/chr3L.pass.snp.phased.haps.vcf

cat $IN_VCF | java -jar /storage/home/cxb585/bin/vcf2beagle.jar '.' data/for.germline.chr3L.pass.snp.phased

# ls -lhrt data/for.germline.chr3L.pass.snp.phased.*
# -rw-rw---- 1 cxb585 cxb585_collab   22 Jan 29 15:44 data/for.germline.chr3L.pass.snp.phased.int
# -rw-rw---- 1 cxb585 cxb585_collab 127M Jan 29 15:44 data/for.germline.chr3L.pass.snp.phased.markers
# -rw-rw---- 1 cxb585 cxb585_collab 112M Jan 29 15:44 data/for.germline.chr3L.pass.snp.phased.bgl.gz

# ----------------------------------------------------------------------------------------
# --- Convert SHAPEIT gen/sample to PLINK PED
# ----------------------------------------------------------------------------------------

IN_HAPS=data/chr3L.pass.snp.phased.haps

~/work/bin/gtool_v0.7.5/gtool -G --g $IN_HAPS --s ${IN_HAPS/.haps/.sample} \
    --ped data/for.germline.`basename $IN_HAPS`.ped \
    --map data/for.germline.`basename $IN_HAPS`.map

Note that genotypes with a probability below the --threshold value are set as missing. GTOOL codes missing genotypes using "N" characters. You can use Plink to recode them using "0" as follows:

plink --file myPlinkTextData --missing-genotype N --output-missing-genotype 0 --recode --out myPlinkTextData2

# Use helper script to call GERMLINE (takes executable, PED, MAP, output file)
~/work/bin/germline-1-5-1/phasing_pipeline/gline.sh \
    ~/work/bin/germline-1-5-1/germline \
    data/for.germline.`basename $IN_HAPS`.ped \
    data/for.germline.`basename $IN_HAPS`.map \
    results/germline_IBD/germline_IBD_3L.out

# TEST: Try replacing "." with some ID (e.g. 3L_353)
head $IN_HAPS > tmp.haps
emacs $IN_HAPS

~/work/bin/gtool_v0.7.5/gtool -G --g tmp.haps --s ${IN_HAPS/.haps/.sample} \
    --ped data/for.germline.`basename $IN_HAPS`.ped \
    --map data/for.germline.`basename $IN_HAPS`.map

# Fails to make anything, probably because data isn't in right format (GEN)?

# ----------------------------------------------------------------------------------------
                                                                                          
# TEST:

# Add ID to haps files before conversion
sed -e "s/^\([^ ]*\) . \([0-9]*\)/\1 \1_\2 \2/" $IN_HAPS | head -n 100000 > ${IN_HAPS/.haps/.wID.haps}
cp ${IN_HAPS/.haps/.sample} ${IN_HAPS/.haps/.wID.sample}

# To convert haplotypes to the HAP/LEGEND/SAMPLE usually used for reference panels of haplotypes in Impute2:
module load shapeit
shapeit -convert \
        --input-haps ${IN_HAPS/.haps/.wID} \
        --output-ref data/for.germline.`basename $IN_HAPS`.hap data/for.germline.`basename $IN_HAPS`.leg data/for.germline.`basename $IN_HAPS`.sam

# Then go on to Germline   
# Part of germline_test.sh    

GERMLINE_OPTS=data/for.germline.`basename $IN_HAPS`.germline_options.run

echo '2' > $GERMLINE_OPTS
echo "data/for.germline.`basename $IN_HAPS`.leg" >> $GERMLINE_OPTS
echo "data/for.germline.`basename $IN_HAPS`.hap" >> $GERMLINE_OPTS
echo "data/for.germline.`basename $IN_HAPS`.sam" >> $GERMLINE_OPTS
echo "data/for.germline.`basename $IN_HAPS`.germline_output" >> $GERMLINE_OPTS

GERMLINE=~/work/bin/germline-1-5-1/germline
MATCH=3

$GERMLINE \
    -bits 128 \
    -err_hom 4 \
    -err_het 1 \
    -min_m $MATCH \
    -from_snp 3L_232 \
    -to_snp 3L_11248270 < $GERMLINE_OPTS
                                                                                          
# ----------------------------------------------------------------------------------------
# --- Convert to PLINK PED format
# ----------------------------------------------------------------------------------------

module load plink

plink --bfile data/chr3L.pass.snp.flt --recode --out data/for.germline.chr3L.pass.snp.flt

mkdir -p results/germline_IBD

##### Use helper script to call GERMLINE (takes executable, PED, MAP, output file)
####~/work/bin/germline-1-5-1/phasing_pipeline/gline.sh \
####    ~/work/bin/germline-1-5-1/germline \
####    data/for.germline.chr3L.pass.snp.flt.ped \
####    data/for.germline.chr3L.pass.snp.flt.map \
####    results/germline_IBD/germline_IBD_3L.out

#### # Use helper script to phase with BEAGLE and then call GERMLINE
#### # (takes PED, MAP, output file)
#### ~/work/bin/germline-1-5-1/phasing_pipeline/run.sh \
####     data/for.germline.chr3L.pass.snp.flt.ped \
####     data/for.germline.chr3L.pass.snp.flt.map \
####     results/germline_IBD/germline_IBD_3L.out


GERMLINE_OPTS=data/for.germline.chr3L.pass.snp.flt.germline_options.run

echo '1' > $GERMLINE_OPTS
echo "data/for.germline.chr3L.pass.snp.flt.map" >> $GERMLINE_OPTS
echo "data/for.germline.chr3L.pass.snp.flt.ped" >> $GERMLINE_OPTS
echo "data/for.germline.chr3L.pass.snp.flt.germline_output" >> $GERMLINE_OPTS

GERMLINE=~/work/bin/germline-1-5-1/germline
MATCH=1

$GERMLINE \
    -bits 128 \
    -err_hom 4 \
    -err_het 1 \
    -min_m $MATCH \
    -from_snp 3L_232 \
    -to_snp 3L_11248270 < $GERMLINE_OPTS
            

# ----------------------------------------------------------------------------------------

WARNING:SNPs::setmapNucleotideToBinary():0 is not one of variant alleles for SNP 3L:41959401
Please ensure that your input data is bi-allelic and contains NO MISSING data
See the GERMLINE web-site (http://www.cs.columbia.edu/~itsik/Software.htm) for
additional information on imputing missing data.

3L	41959401	.	C	T	94.33	PASS	AC=3;AF=0.013;AN=230;BaseQRankSum=-4.454;DP=1503;Dels=0.00;FS=0.000;HaplotypeScore=1.5419;InbreedingCoeff=-0.0221;MLEAC=3;MLEAF=0.013;MQ=55.17;MQ0=0;MQRankSum=-0.229;QD=2.36;ReadPosRankSum=2.215;SOR=0.678	GT:AD:DP:GQ:PL	0/0:11,1:12:33:0,33,434	0/0:27,3:30:26:0,26,757	0/0:2,0:2:6:0,6,80	0/0:2,0:2:6:0,6,80	0/0:15,0:15:45:0,45,583	0/0:15,0:15:39:0,39,504	0/0:7,0:7:21:0,21,274	0/0:9,0:9:27:0,27,343	0/0:15,0:15:39:0,39,518	0/0:17,0:17:45:0,45,5930/0:15,0:15:42:0,42,565	0/0:15,0:15:45:0,45,585	0/0:15,0:15:39:0,39,504	0/0:4,0:4:12:0,12,157	0/0:16,1:17:36:0,36,473	0/1:11,5:16:97:97,0,367	0/0:22,3:25:17:0,17,653	0/0:13,0:13:36:0,36,473	0/0:7,0:7:21:0,21,263	0/0:19,1:20:45:0,45,599	0/0:16,0:16:48:0,48,626	0/1:10,3:13:32:32,0,361	0/0:14,0:14:39:0,39,519	0/0:13,0:13:33:0,33,450	0/0:19,0:19:54:0,54,701	0/0:12,0:12:30:0,30,395	0/0:16,0:16:39:0,39,479	0/0:10,0:10:27:0,27,359	0/0:11,0:11:33:0,33,414	0/0:10,0:10:27:0,27,353	0/0:6,0:6:15:0,15,199	0/0:11,0:11:30:0,30,3790/0:9,0:9:27:0,27,357	0/0:10,0:10:30:0,30,391	0/0:18,0:18:54:0,54,706	0/0:20,0:20:54:0,54,708	0/0:8,0:8:24:0,24,306	0/0:13,0:13:36:0,36,470	0/0:7,0:7:18:0,18,233	0/0:22,0:22:63:0,63,819	0/0:11,0:11:30:0,30,396	0/0:8,0:8:24:0,24,306	0/0:2,0:2:6:0,6,77	0/0:15,0:15:39:0,39,493	0/0:18,0:18:48:0,48,616	0/1:9,2:11:41:41,0,271	0/0:16,0:16:45:0,45,576	0/0:18,0:18:48:0,48,631	0/0:35,0:35:93:0,93,1225	0/0:15,0:15:39:0,39,496	0/0:11,0:11:30:0,30,387	0/0:5,0:5:15:0,15,192	0/0:15,0:15:39:0,39,496	0/0:15,0:15:42:0,42,523	0/0:11,0:11:30:0,30,399	0/0:16,0:16:45:0,45,587	0/0:13,0:13:36:0,36,469	0/0:9,0:9:21:0,21,265	0/0:8,0:8:18:0,18,236	0/0:14,0:14:36:0,36,465	0/0:7,0:7:21:0,21,257	0/0:14,1:15:7:0,7,499	0/0:22,0:22:60:0,60,769	0/0:16,0:16:48:0,48,598	0/0:12,0:12:30:0,30,392	./.	0/0:16,0:16:45:0,45,556	0/0:27,0:27:63:0,63,852	0/0:19,0:19:57:0,57,7400/0:12,0:12:30:0,30,394	0/0:10,0:10:24:0,24,310	0/0:8,0:8:21:0,21,257	0/0:21,0:21:57:0,57,750	0/0:9,0:9:21:0,21,274	0/0:11,0:11:30:0,30,386	0/0:10,0:10:30:0,30,377	0/0:11,0:12:27:0,27,348	0/0:8,0:8:24:0,24,317	0/0:15,0:15:45:0,45,557	0/0:13,0:13:33:0,33,436	0/0:8,0:8:24:0,24,311	0/0:11,0:11:33:0,33,419	0/0:16,0:16:39:0,39,494	0/0:6,0:6:15:0,15,196	0/0:12,0:12:33:0,33,389	0/0:19,0:19:57:0,57,725	0/0:7,0:7:21:0,21,260	0/0:14,0:14:30:0,30,414	0/0:13,0:13:33:0,33,434	0/0:20,0:20:51:0,51,656	0/0:17,0:17:48:0,48,6230/0:14,1:15:36:0,36,452	0/0:12,0:12:36:0,36,454	0/0:20,0:20:54:0,54,721	0/0:9,0:9:27:0,27,337	0/0:8,0:8:21:0,21,279	0/0:18,0:18:54:0,54,708	0/0:6,0:6:18:0,18,237	0/0:16,0:16:45:0,45,599	0/0:14,0:14:42:0,42,540	0/0:6,0:6:15:0,15,202	0/0:15,0:15:39:0,39,517	0/0:4,0:4:12:0,12,160	0/0:10,0:10:27:0,27,362	0/0:19,0:19:54:0,54,691	0/0:4,0:4:12:0,12,144	0/0:11,0:11:30:0,30,397	0/0:21,0:21:54:0,54,707	0/0:16,0:16:42:0,42,556	0/0:13,1:14:8:0,8,477	0/0:4,0:4:12:0,12,152	0/0:22,0:22:60:0,60,765	0/0:11,0:11:30:0,30,3940/0:1,0:1:3:0,3,40	0/0:21,0:21:60:0,60,765	0/0:3,0:3:6:0,6,85







# ----------------------------------------------------------------------------------------

WARNING:SNPs::setmapNucleotideToBinary():A is not one of variant alleles for SNP 3L:41960264
Please ensure that your input data is bi-allelic and contains NO MISSING data
See the GERMLINE web-site (http://www.cs.columbia.edu/~itsik/Software.htm) for
additional information on imputing missing data.

3L	41960264	.	C	A	523.68	PASS	0/0:29,0:29:84:0,84,1107	0/0:24,0:24:63:0,63,824./.	0/0:4,0:4:12:0,12,157	0/0:22,0:22:54:0,54,731	0/0:5,0:5:15:0,15,197	0/0:21,0:21:57:0,57,770	0/0:16,0:16:45:0,45,589	0/0:10,0:10:24:0,24,327	0/0:17,0:17:51:0,51,6770/0:15,0:15:39:0,39,518	0/0:29,0:29:81:0,81,1072	0/0:17,0:17:48:0,48,645	0/0:18,0:18:45:0,45,608	0/0:18,0:18:48:0,48,650	0/0:23,0:23:60:0,60,814	0/0:22,0:22:51:0,51,6990/0:25,0:25:69:0,69,921	0/0:21,0:21:57:0,57,753	0/0:16,0:16:45:0,45,599	0/0:14,0:14:42:0,42,557	0/0:14,0:14:36:0,36,482	0/0:19,0:19:48:0,48,654	0/0:9,0:9:21:0,21,285	0/0:29,0:29:75:0,75,1011	0/0:19,0:19:42:0,42,566	0/0:18,0:18:48:0,48,650	0/1:9,8:17:99:232,0,309	0/0:22,0:22:66:0,66,860	0/0:22,0:22:63:0,63,845	0/0:29,0:29:75:0,75,1006	0/0:23,0:23:57:0,57,779	0/0:18,0:18:48:0,48,642	0/0:16,0:16:45:0,45,590	0/0:12,0:12:36:0,36,469	0/0:31,0:31:84:0,84,1121	0/0:37,0:37:96:0,96,1263	0/0:28,0:28:81:0,81,1066	0/0:14,0:14:39:0,39,509	0/0:20,0:20:60:0,60,797	0/0:15,0:15:39:0,39,527	0/0:18,0:18:51:0,51,676	0/0:7,0:7:15:0,15,201	0/0:7,0:7:21:0,21,277	0/0:23,0:23:54:0,54,707	0/0:19,0:19:54:0,54,709	0/0:19,0:19:48:0,48,645	0/0:30,0:30:87:0,87,1131	0/0:19,0:19:54:0,54,717	0/0:24,0:24:72:0,72,920	0/0:24,0:24:63:0,63,854	0/0:17,0:17:45:0,45,607	0/0:23,0:23:60:0,60,790	0/0:21,0:21:60:0,60,789	0/0:27,0:27:75:0,75,995	0/0:14,0:14:42:0,42,554	0/0:38,0:38:99:0,105,1397	0/0:29,0:29:78:0,78,1054	0/0:18,0:19:48:0,48,627	0/0:20,0:20:54:0,54,725	0/0:23,0:23:63:0,63,799	0/0:23,0:23:66:0,66,879	0/0:35,0:35:87:0,87,1153	0/0:29,0:29:84:0,84,1116	0/0:23,0:23:66:0,66,8740/0:3,0:3:9:0,9,120	0/0:26,0:26:69:0,69,928	0/0:17,0:17:42:0,42,574	0/0:32,0:32:90:0,90,1171	0/0:25,0:25:69:0,69,895	0/0:16,0:16:48:0,48,637	0/0:16,0:16:45:0,45,5940/1:11,12:23:99:352,0,272	0/0:9,0:9:27:0,27,360	0/0:30,0:30:81:0,81,1083	0/0:21,0:21:57:0,57,761	0/0:20,0:20:57:0,57,762	0/0:12,0:12:33:0,33,436	0/0:20,0:20:54:0,54,727	0/0:26,0:26:72:0,72,952	0/0:17,0:17:48:0,48,636	0/0:22,0:22:57:0,57,751	0/0:17,0:17:48:0,48,624	0/0:16,0:16:42:0,42,548	0/0:32,0:32:87:0,87,1156	0/0:34,0:34:81:0,81,1080	0/0:25,0:25:69:0,69,921	0/0:22,0:22:60:0,60,802	0/0:26,0:26:72:0,72,970	0/0:19,0:19:48:0,48,647	0/0:20,0:20:54:0,54,715	0/0:17,0:17:48:0,48,619	0/0:17,0:17:45:0,45,604	0/0:30,0:30:81:0,81,1075	0/0:15,0:15:45:0,45,600	0/0:14,0:14:42:0,42,560	0/0:29,0:29:81:0,81,1081	0/0:16,0:16:39:0,39,520	0/0:18,0:18:48:0,48,650	0/0:17,0:17:51:0,51,680	0/0:16,0:16:48:0,48,637	0/0:16,0:16:48:0,48,640	0/0:26,0:26:66:0,66,888	0/0:27,0:27:66:0,66,866	0/0:31,0:31:81:0,81,1093	0/0:12,0:12:30:0,30,410	0/0:14,0:14:36:0,36,470	0/0:20,0:20:51:0,51,671	0/0:42,0:42:99:0,111,1484	0/0:36,0:36:84:0,84,1147	0/0:15,0:15:42:0,42,543	0/0:36,0:36:93:0,93,1233	0/0:7,0:7:21:0,21,277	0/0:14,0:14:39:0,39,504	0/0:16,0:16:42:0,42,567	0/0:3,0:3:9:0,9,120

# ----------------------------------------------------------------------------------------

WARNING:SNPs::setmapNucleotideToBinary():T is not one of variant alleles for SNP 3L:41960294
Please ensure that your input data is bi-allelic and contains NO MISSING data
See the GERMLINE web-site (http://www.cs.columbia.edu/~itsik/Software.htm) for
additional information on imputing missing data.
Reading Markers Complete
Match Markers
Matching Markers Complete
Matches completed ... freeing memory