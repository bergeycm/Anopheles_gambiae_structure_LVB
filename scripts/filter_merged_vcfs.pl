#!/usr/bin/perl

# ----------------------------------------------------------------------------------------
# --- Filter merged VCF - Remove inaccessible sites and tally up missingness
# ----------------------------------------------------------------------------------------

use strict;
use warnings;

use Data::Dumper;

my $in_vcf    = shift;   # Merged VCF
my $acc_vcf   = shift;   # Accessibility VCF, e.g. data/accessibility/accessibility.2L.vcf
chomp $acc_vcf;

# ----------------------------------------------------------------------------------------

# --- Open output file to which to write per-dataset missingness info

# Infer chromosome from accessibility file name
my $chr = $acc_vcf;
$chr =~ s:data/accessibility/accessibility\.(.*)\.vcf:$1:;

print STDERR "Processing SNPs for chr$chr...\n";

#my $missingness_file = "data/ssese_with_ag1000g/ssese_with_ag1000g.missingness.$chr.txt";

my $missingness_file = in_vcf;
$missingness_file =~ s/ssese_with_ag1000g\./ssese_with_ag1000g.missingness/;
$missingness_file =~ s/vcf$/txt/;
open (MISS, ">$missingness_file")
    or die "ERROR: Could not open output missingness file [$missingness_file].\n";

print STDERR "Writing file with missingness info to [$missingness_file].\n";

# ----------------------------------------------------------------------------------------

# --- Load in accessibility info

open (ACC, "<$acc_vcf")
    or die "ERROR: Could not open input accessibility VCF [$acc_vcf].\n";

my %acc;

while (my $line = <ACC>) {

    next if ($line =~ /^#/);

    my @acc_info = split /\t/, $line;

    if ($acc_info[1] % 1000000 == 0) {
        print STDERR "Reading line $acc_info[1]...\n";
    }

    if ($acc_info[6] eq "PASS") {
        $acc{$acc_info[1]} = 1;
    }

}

# ----------------------------------------------------------------------------------------

# Counts of various outcomes for each SNP:
my $missing_in_ssese = 0;
my $missing_in_ag1000g_acc = 0;
my $missing_in_ag1000g_inacc = 0;
my $ok_in_both = 0;

open (VCF, "<$in_vcf")
    or die "ERROR: Could not open input VCF [$in_vcf].\n";

my @inds;
my @dataset;
my @ss_indices;
my @ag_indices;

my $line_count = 0;

while (my $line = <VCF>) {

	$line_count++;

    chomp $line;

    if ($line =~ /^#/) {

        if ($line =~ /^#CHROM/) {
            my @header_info = split (/\t/, $line);
            @inds = @header_info[9 .. (scalar (@header_info) - 1)];

            my @dataset = map { (/^ssese/) ? "SS" : "AG" } @inds;

            # Also get array of indices for each dataset
            @ss_indices = grep { $dataset[$_] eq "SS" } 0..$#dataset;
            @ag_indices = grep { $dataset[$_] eq "AG" } 0..$#dataset;

        }

        next;
    }

    my @info = split(/\t/, $line);

    if ($line_count % 1000 == 0) {
        print STDERR "Processing VCF line $line_count, SNP ";
        print STDERR $info[0] . ":" . $info[1] . "...\n";
    }

    # --- Check accessibility info

    my $is_acc = 0;

    if (exists $acc{$info[1]}) {
       $is_acc = 1;
    }

    # --- Parse genotypes to figure out missingness by dataset

    my @genos = @info[9 .. (scalar (@info) - 1)];

    for (@genos) {
        $_ =~ s/:.*//g;
    }

    # --- Find missing in Ssese dataset

    my $ss_total   = 0;
    my $ss_missing = 0;

    for my $ss_i (@ss_indices) {
        $ss_total++;
        if ($genos[$ss_i] eq "./.") {
            $ss_missing++;
        }
    }

    my $ss_missing_perc = $ss_missing / $ss_total;

    # --- Find missing in Ag1000G dataset

    my $ag_total   = 0;
    my $ag_missing = 0;

    for my $ag_i (@ag_indices) {
        $ag_total++;
        if ($genos[$ag_i] eq "./.") {
            $ag_missing++;
        }
    }

    my $ag_missing_perc = $ag_missing / $ag_total;

    # print "SS MISSING: $ss_missing_perc\n";
    # print "AG MISSING: $ag_missing_perc\n\n";
    print MISS join "\t", @info[0 .. 1];
    print MISS "\t";
    print MISS "$ss_missing_perc\t$ag_missing_perc\n";

    # Inaccessible: remove entirely
    next if (! $is_acc);

    # Missing in Ssese
    if ($ag_missing_perc < 1 && $ss_missing_perc == 1) {

        $missing_in_ssese++;

    # Missing in Ag1000G (but Accessible)
    } elsif ($ag_missing_perc == 1 && $ss_missing_perc < 1) {

        $missing_in_ag1000g_acc++;

    # Not completely missing in both
    } else {

        $ok_in_both++;

    }

    # --- Print line
    # Print first columns (before genotypes)
    print join "\t", @info[0 .. 8];
    print "\t";

    for my $geno (@genos) {
        print "$geno\t";
    }

    print "\n";

}

print STDERR "SNP COUNTS:\n";
print STDERR "missing_in_ssese:\t"         . $missing_in_ssese . "\n";
print STDERR "missing_in_ag1000g_acc:\t"   . $missing_in_ag1000g_acc . "\n";
print STDERR "missing_in_ag1000g_inacc:\t" . $missing_in_ag1000g_inacc . "\n";
print STDERR "ok_in_both:\t"               . $ok_in_both . "\n";

close VCF;
close MISS;
close ACC;

exit;
