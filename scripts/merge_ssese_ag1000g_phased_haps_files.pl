#!/usr/bin/perl

# ----------------------------------------------------------------------------------------
# --- Merge in phased Ag1000G data
# ----------------------------------------------------------------------------------------

#    $ head data/chrX.pass.snp.phased.haps | cut -d' ' -f 1-30
#    5 X:7 7 C T 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#    5 X:23 23 A T 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

#    $ head data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.1.haplotypes.X.haps | cut -d' ' -f 1-30
#    0 . 49 A T 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#    0 . 80 G T 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

# Call with:
# mkdir -p data/ssese_with_ag1000g
# perl scripts/merge_ssese_ag1000g_phased_haps_files.pl X > data/ssese_with_ag1000g/chrX.pass.snp.phased.ag1000g.haps

use strict;
use warnings;

my $chr = shift;
chomp $chr;

my $ssese_file = "data/chr${chr}.pass.snp.phased.haps";
my $ag1000g_file = "data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.1.haplotypes.";
$ag1000g_file .= "${chr}.haps";

print STDERR "Processing [$ssese_file] and [$ag1000g_file]...\n";

# ----------------------------------------------------------------------------------------

sub read_file_line {
    my $fh = shift;

    if ($fh and my $line = <$fh>) {
        chomp $line;
        return [ split(/\s/, $line) ];
    }
    return;
}

sub write_combined_line {
    # do something with the 2 values
    my $pair_ss = shift;
    my $pair_ag = shift;

    #print "\tMatch! [$pair_ss] [$pair_ag]\n";

    my @pair_ss = @{ $pair_ss };
    my @pair_ag = @{ $pair_ag };

    # Save first three columns
    my $first_cols = "$pair_ss[0] $pair_ss[1] $pair_ss[2]";

    # Check  if 1.) bases match, 2.) bases switched, 3.) bases completely different.
    my $flip_bases = 0;
    if ($pair_ss[3] eq $pair_ag[3] && $pair_ss[4] eq $pair_ag[4]) {
        #print STDERR "\tTotal match\n";
        print "$first_cols $pair_ss[3] $pair_ss[4]";
        # Print Ag1000G data
        for (my $i = 5; $i < scalar @pair_ag; $i++) {
            print " " . $pair_ag[$i];
        }
        # Print Ssese data
        for (my $i = 5; $i < scalar @pair_ss; $i++) {
            print " " . $pair_ss[$i];
        }
        print "\n";
    } elsif ($pair_ss[3] eq $pair_ag[4] && $pair_ss[4] eq $pair_ag[3]) {
        #print STDERR "\tFlipped basepairs\n";
        $flip_bases = 1;
        print "$first_cols $pair_ss[4] $pair_ss[3]";
        # Print Ag1000G data
        for (my $i = 5; $i < scalar @pair_ag; $i++) {
            print " " . $pair_ag[$i];
        }
        # Print Ssese data FLIPPED
        for (my $i = 5; $i < scalar @pair_ss; $i++) {
            if ($pair_ss[$i] eq "0") {
                print " 1";
            } elsif ($pair_ss[$i] eq "1") {
                print " 0";
            }
        }
        print "\n";
    } else {
        #print STDERR "\tERROR: Bases don't match\n";
    }

    return;
}

# ----------------------------------------------------------------------------------------

# Inspired by http://stackoverflow.com/a/2499339

open(my $f_ss, $ssese_file);
open(my $f_ag, $ag1000g_file);

my $line_ss = read_file_line($f_ss);
my $line_ag = read_file_line($f_ag);

while ($line_ss and $line_ag) {
    #print "SS: " . $line_ss->[2] . "\n";
    #print "AG: " . $line_ag->[2] . "\n";
    if ($line_ss->[2] < $line_ag->[2]) {
        #print STDERR "\tSkipping ahead in Ssese file...\n";
        $line_ss = read_file_line($f_ss);
    } elsif ($line_ag->[2] < $line_ss->[2]) {
        #print STDERR "\tSkipping ahead in Ag1000G file...\n";
        $line_ag = read_file_line($f_ag);
    } else {
        write_combined_line($line_ss, $line_ag);
        $line_ss = read_file_line($f_ss);
        $line_ag = read_file_line($f_ag);
    }

}

close($f_ss);
close($f_ag);
