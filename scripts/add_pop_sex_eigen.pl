#!/usr/bin/perl

use strict;
use warnings;

# ----------------------------------------------------------------------------------------
# --- Add sex and population information to *ind file in EIGENSTRAT analysis (w/ Ag1000G)
# ----------------------------------------------------------------------------------------

my $eigen_ind = shift;
chomp $eigen_ind;

my $ind_fixed = $eigen_ind;
$ind_fixed =~ s/\.ind$/.fixed.ind/;

my $ind_info    = "data/ssese_individual_info_simple_bugala_split.txt";
my $ind_info_ag = "data/ag1000g.phase1.ar3/samples.all.txt";

open (EIG, "<$eigen_ind")
    or die "ERROR: Could not open input individual file. $!\n";

open (FIX, ">$ind_fixed")
    or die "ERROR: Could not open output individual file. $!\n";

while (<EIG>) {

    chomp;
    my @info = split / +/;

    my $this_ind = $info[1];
    $this_ind =~ s/\.variant.*//g;

    # Find info on this individual
    my $grep_pop = "grep \"" . $this_ind . "\" $ind_info | cut -d ' ' -f3";
    my $grep_sex = "grep \"" . $this_ind . "\" $ind_info | cut -d ' ' -f5";

    my $this_pop = `$grep_pop`;
    chomp $this_pop;
    my $this_sex = `$grep_sex`;
    chomp $this_sex;

    if ($this_pop eq "") {

        # Find info on this individual
        $grep_pop = "grep \"" . $this_ind . "\" $ind_info_ag | cut -f6";
        $grep_sex = "grep \"" . $this_ind . "\" $ind_info_ag | cut -f12";

        $this_pop = `$grep_pop`;
        chomp $this_pop;
        $this_sex = `$grep_sex`;
        chomp $this_sex;
    }

    # --- Set pop and sex for non-Ssese animals
    if ($this_ind =~ /^TZ/) {
        $this_pop = "Tanzania_AG";
        $this_sex = "u";
    }
    if ($this_ind =~ /^SRS/) {
        $this_pop = "Tanzania_AR";
        $this_sex = "u";
    }
    if ($this_ind =~ /^CHRISTYI/) {
        $this_pop = "CHRISTYI";
        $this_sex = "u";
    }

    print STDERR "Ind $this_ind is sex $this_sex and population $this_pop\n";

    print FIX ' ' x 12;
    print FIX $this_ind . ' ';
    print FIX uc($this_sex);
    print FIX ' ' x 8;
    print FIX $this_pop . "\n";

}

close EIG;
close FIX;

exit;
