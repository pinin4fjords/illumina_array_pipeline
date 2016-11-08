#!/usr/bin/perl -w

use strict;

my ( $illuminafile, $probes, $controls ) = @ARGV;

my $cat = undef;

my %lines;

open( IN, $illuminafile );
while ( my $line = <IN> ) {
    if ( $line =~ /^\[([^\]]+)\]/ ) {
        $cat = $1;
        next;
    }
    elsif ( !$cat ) {
        next;
    }
    else {
        $lines{$cat} .= $line;
    }
}
close(IN);

# The main data is 'Probes' in the expression arrays, 'Assay' in Infinium

open (PROBES, ">$probes");
if (defined $lines{"Probes"}){
    print PROBES $lines{"Probes"};
}
elsif(defined $lines{"Assay"}){
    print PROBES $lines{"Assay"};
}
close(PROBES);

open (CONTROLS, ">$controls");
print CONTROLS $lines{"Controls"};
close(CONTROLS);

