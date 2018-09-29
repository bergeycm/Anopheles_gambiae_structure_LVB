#!/usr/bin/perl

use strict;
use warnings;

# ----------------------------------------------------------------------------------------
# --- Make input file for H12 program
# ----------------------------------------------------------------------------------------

my $in_haps = shift;     # E.g. data/nohetinv/chr3L.pass.snp.phased.nohetinv.haps
chomp $in_haps;

open (my $fh, '<', $in_haps)
    or die "Could not open input file '$in_haps' $!";

while (my $row = <$fh>) {
    chomp $row;

    my @info = split ' ', $row;
    my $ref_bp = $info[3];
    my $alt_bp = $info[4];

    print $info[2];

    for (my $i = 5; $i < scalar @info; $i++) {

        if ($info[$i] == "0") {
            print ",", $ref_bp;
        } else {
            print ",", $alt_bp;
        }
    }

    print "\n";
}

exit;
