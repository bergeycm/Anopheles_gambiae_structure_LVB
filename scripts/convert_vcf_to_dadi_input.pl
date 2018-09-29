#!/usr/bin/perl

# ========================================================================================
# --- Script by Kun Wang to generate dadi input from VCF file
# --- Downloaded from dadi user group at:
# ---  https://groups.google.com/forum/#!topic/dadi-user/p1WvTKRI9_0
# ---  and slightly modified.
# ========================================================================================

use strict;
use warnings;
use Bio::SeqIO;

my $file      = shift;
my $vcf       = shift;
my $list      = shift;
my $start_SNP = shift || 0;
my $end_SNP   = shift || -1;

# $file is the location of genome file; $list is like this:
# BEGIN
# sample1    population1
# sample2    population1
# sample3    population2
# END
# the sample name should be same as which in vcf file.
# Usage: perl convert_vcf_to_dadi_input.pl <genome file> <vcf file> <list file>

print STDERR "Building list of individuals...\n";
my %list;
open(IN,"< $list") || die "$!";
while (<IN>) {
    chomp;
    next if(/^\#/);
    my @a = split(/\s+/);
    $list{$a[0]} = $a[1];
}
close IN;

print STDERR "Starting to read VCF...\n";
my $SNP_count = 0;
my %vcf;
my %pop;
my %record;
open(IN,"< $vcf");
while (<IN>) {
    chomp;
    next if(/^##/);
    if (/^#/){
        my @a = split(/\s+/);
        for(my $i=9; $i<@a; $i++){
            next if (!exists $list{$a[$i]});
            #die"$a[$i] does not exists in $list\n" if(!exists $list{$a[$i]});
            $record{$i} = $list{$a[$i]};
        }
        next;
    }

    $SNP_count++;
    if ($SNP_count % 10000 == 0) {
        print STDERR "\tProcessing SNP number $SNP_count...\n";
    }

    # Skip ahead if we're outside of user-specified window
    if ($end_SNP != -1) {
        if ($SNP_count < $start_SNP) {
            next;
        } elsif ($SNP_count > $end_SNP) {
            last;
        }
    }

    my @a = split(/\s+/);
    my ($chr,$pos,$ref,$alt) = ($a[0],$a[1],$a[3],$a[4]);
    next if ($alt=~/,/);

    $vcf{$chr}{$pos}{ref} = $ref;
    $vcf{$chr}{$pos}{alt} = $alt;
    #print "DEBUG\t$ref\t$alt\n";
    foreach my $i (keys %record){
        my $indv = $record{$i};
        $pop{$indv} = 1;
        #print "INDV [$indv]\t";
        my $geno = $a[$i];
        #print "GENO [$geno]\t";
        $geno =~ /^(.)\/(.)/;
        my $tempa = $1;
        my $tempb = $2;
        #print "TEMP: [$tempa $tempb]\n";
        my ($a1,$a2) = (0,0);
        if ($tempa eq "." || $tempb eq "."){
            $a1 = 0;
            $a2 = 0;
        } elsif ($tempa + $tempb == 1) {
            $a1 = 1;
            $a2 = 1;
        }
        elsif ($tempa + $tempb == 2) {
            $a1 = 0;
            $a2 = 2;
        }
        elsif ($tempa + $tempb == 0) {
            $a1 = 2;
            $a2 = 0;
        }
        $vcf{$chr}{$pos}{$indv}{a1} += $a1;
        $vcf{$chr}{$pos}{$indv}{a2} += $a2;
    }
    #last;
}
close IN;

print STDERR "Finished processing VCF file.\n";

# ----------------------------------------------------------------------------------------

# If doing whole file, instead of subset, output is just [pop_list].data
# Otherwise it's [pop_list].[start]-[end].data
my $out_file;
if ($end_SNP == -1) {
    $out_file = "$list.data";
} else {
    $out_file = "$list.${start_SNP}-${end_SNP}.data";
}

print STDERR "Writing output file [$out_file]...\n";

open(O,"> $out_file");
my $title = "NAME\tOUT\tAllele1";
foreach my $pop (sort keys %pop){
    $title .= "\t$pop";
}
$title .= "\tAllele2";
foreach my $pop (sort keys %pop){
    $title .= "\t$pop";
}
$title .= "\tGene\tPostion\n";
print O "$title";
my $fa = Bio::SeqIO -> new(-file=>$file, -format=>'fasta');
while (my $seq=$fa -> next_seq){
    my $id=$seq -> id;
    my $seq=$seq -> seq;
    foreach my $pos (sort {$a<=>$b} keys %{$vcf{$id}}){
        my $ref = substr($seq, $pos-2, 3);
        my $line = "$ref\t$ref\t$vcf{$id}{$pos}{ref}";
        foreach my $pop (sort keys %pop){
            my $num = $vcf{$id}{$pos}{$pop}{a1};
            $line .= "\t$num";
        }
        $line .= "\t$vcf{$id}{$pos}{alt}";
        foreach my $pop (sort keys %pop){
            my $num = $vcf{$id}{$pos}{$pop}{a2};
            $line .= "\t$num";
        }
        $line .= "\t$id\t$pos\n";
        print O "$line";
    }
}
close O;
