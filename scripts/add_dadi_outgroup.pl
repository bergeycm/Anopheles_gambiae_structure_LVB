#!/usr/bin/perl

# ----------------------------------------------------------------------------------------
# --- Add outgroup info to dadi file
# ----------------------------------------------------------------------------------------

use strict;
use warnings;

# module load tabix

my $verbose = 0;

my $in_dadi = shift;
chomp $in_dadi;

my $vcf_mela = "data/AGC/mela_ref_ug_vqsr_cnvrt_sort.vcf.gz";
my $vcf_meru = "data/AGC/meru_ref_ug_vqsr_cnvrt_sort.vcf.gz";

open my $fh, '<', $in_dadi
    or die "Cannot open $in_dadi: $!";

print STDERR "Processing input file [$in_dadi]..." . localtime() . "\n";

my $count = 0;

while (my $line = <$fh>) {

    $count++;
    if ($count % 10000 == 0) {
        print STDERR "Processing line $count... " . localtime() . "\n";
    }

    chomp $line;
    my @row = split( "\t", $line );

    my $chr = $row[(scalar @row) - 2];
    my $pos = $row[(scalar @row) - 1];

    print STDERR "Looking up $chr:$pos..." if $verbose;

    my $tabix_res_mela = `tabix $vcf_mela $chr:$pos-$pos | cut -f 4`;
    my $tabix_res_meru = `tabix $vcf_meru $chr:$pos-$pos | cut -f 4`;

    chomp $tabix_res_mela;
    chomp $tabix_res_meru;

    print STDERR "\t[$tabix_res_mela] - [$tabix_res_meru]\t" if $verbose;

    # Break if positions exist but don't agree or if both are missing
    my $match = $tabix_res_mela eq $tabix_res_meru;
    my $mismatch = $tabix_res_mela ne $tabix_res_meru;
    my $both_missing = ($tabix_res_mela eq "") && ($tabix_res_meru eq "");
    my $both_present = ($tabix_res_mela ne "") && ($tabix_res_meru ne "");

    if (($both_present && $mismatch) || $both_missing) {
        print STDERR "SKIPPING\n" if $verbose;
        next;
    }

    my $out_bp;
    if ($tabix_res_mela ne "") {
        $out_bp = $tabix_res_mela;
    } else {
        $out_bp = $tabix_res_meru;
    }

    print STDERR "BASE: [$out_bp]\t" if $verbose;

    my $orig_val = uc(substr($row[1],1,1));

    print STDERR $row[1] . "\t" if $verbose;
    substr($row[1],1,1) = $out_bp;
    print STDERR $row[1] . "\n" if $verbose;

    if ($orig_val ne $out_bp) {
    	print STDERR "Mismatch!\n" if $verbose;
        next;
    }

    print join "\t", @row;
    print "\n";
}

close($fh);

print STDERR "Finished processing input file [$in_dadi]." . localtime() . "\n";

exit;
