#!/usr/bin/env perl

use warnings;
use strict;

while (my $line = <>) {
	if ($line =~ /^#/) {
		next;
	}
	chomp $line;
	my ($score, $chr, $start, $stop, $strand_symbol, $rmsk_class) = split(/\t/, $line);
	print join("\t", $chr, $start, $stop, $rmsk_class, $score, $strand_symbol) . "\n";
}
