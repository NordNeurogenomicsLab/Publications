#!/usr/bin/perl
# Gene_Annotator.pl by Rinaldo Catta-Preta @UCDavis

use strict;
use warnings 'FATAL' => 'all';
use Data::Dumper;

my ($peak_annot_dir, $prefix) = @ARGV;

# Building the gene symbol list
my $filexref = "/share/nordlab/codes/chipseq/kgXref.txt";

my %genesymbols = ();

open (my $ingene, '<', $filexref) or die "Cannot open $filexref\n"; 
print "\n Building gene symbol list from UCSC gene names ...\n\n ";
while (<$ingene>)
{
	chomp; 
	my @xrefgene = split /\t/;
	$genesymbols{$xrefgene[0]}  = $xrefgene[4];
}
close $ingene;


# Procedure to obtain gene symbol from UCSC name
my $mygenesymbol;
sub genesymbol
{
	my $object = $_[0];
	
	if (exists $genesymbols{$object})	{	$mygenesymbol = $genesymbols{$object};	}
	else					{	$mygenesymbol = $object;	}
	
	return $mygenesymbol;
}


# Convert files
my $newfilename;

my @typelist = ("gappedPeak");

foreach my $type (@typelist)
{
	# _peaks.summary files
	my $experiment = 'peaks.summary.' . "$type";
	my $filename = $peak_annot_dir . $prefix . '_' . $experiment;

	if (-f $filename)
	{
		open (my $in, '<', "$filename") or die "Cannot open $filename\n";
		$newfilename = $filename . '.annotated';
		open (my $out, '>', $newfilename) or die "Cannot open $newfilename\n";

		while (<$in>)
		{
			chomp;
			next if /^\s*$/;
			if (/^Chromosome/)
			{
				print $out "$_\tPAP\tPAE\n";
				next;
			}
		
			my @data = split /\t/;
			my @multigene = split /\,/, $data[3] if (defined $data[3]);
			my @newmultigene = ();
		
			$data[4] = genesymbol($data[4]) if (defined $data[4]);
		
			foreach my $gene (@multigene)
			{
				push @newmultigene, genesymbol($gene) if ($gene ne '');
			}
		
			$data[3] = join ',', @newmultigene;
		
			if (defined $data[5])
			{
				if    ($data[5] < 1000)		{	$data[6] = 1; $data[7] = 0;	}
				elsif ($data[5] < 10000 )	{	$data[6] = 0; $data[7] = 0;	}
				elsif ($data[5] > 10000 &&
				       $data[5] < 1000000)	{	$data[6] = 0; $data[7] = 1;	}
				else				{	$data[6] = 0; $data[7] = 0;	}
			}
		
			print $out $data[0];
		
			for (my $i = 1; $i < @data; $i++)
			{
				print $out "\t$data[$i]";
			}
		
			print $out "\n";
		}
		close $in;
		close $out;
	}
	else
	{ print "$filename does not exists\n";}
	print "\n Finished working on $filename\n";


	# _peaks.ndg files
	$experiment = 'peaks.ndg.' . "$type";
	$filename = $peak_annot_dir . $prefix . '_' . $experiment;

	if (-f $filename)
	{
		open (my $in, '<', "$filename") or die "Cannot open $filename\n";
		$newfilename = $filename . '.annotated';
		open (my $out, '>', $newfilename) or die "Cannot open $newfilename\n";

		while (<$in>)
		{
			chomp;
			next if /^\s*$/;
			if (/^Chromosome/)
			{
				my @header = split /\t/;
				my $newheader = join("\t", @header[0 .. 6]);
				$newheader.= "\tPAP\tPAE\t";
				$newheader.= join("\t", @header[7 .. 9]);
				$newheader.= "\tPAP\tPAE";
				print $out "$newheader\n";
				next;
			}
		
			my @data = split /\t/;
		
			$data[5] = genesymbol($data[4]) if (defined $data[4]);
			$data[8] = genesymbol($data[7]) if (defined $data[7]);

			if (defined $data[6] and $data[6] ne "")
			{
				if    ($data[6] < 1000)		{	$data[10] = 1; $data[11] = 0;	}
				elsif ($data[6] < 10000 )	{	$data[10] = 0; $data[11] = 0;	}
				elsif ($data[6] > 10000 &&
				       $data[6] < 1000000)	{	$data[10] = 0; $data[11] = 1;	}
				else				{	$data[10] = 0; $data[11] = 0;	}
			}
		
			if (defined $data[9])
			{
				if    ($data[9] < 1000)		{	$data[12] = 1; $data[13] = 0;	}
				elsif ($data[9] < 10000 )	{	$data[12] = 0; $data[13] = 0;	}
				elsif ($data[9] > 10000 &&
				       $data[9] < 1000000)	{	$data[12] = 0; $data[13] = 1;	}
				else				{	$data[12] = 0; $data[13] = 0;	}
			}
		
			print $out $data[0];

			for (my $i = 1; $i < 7; $i++)
			{
				print $out "\t$data[$i]" if (defined $data[$i]);
				print $out "\t" if (! defined $data[$i]);
			}
		
			if (defined $data[6] and $data[6] ne "")
			{
				print $out "\t$data[10]\t$data[11]" if (defined $data[10] and defined $data[11]);
			}
			else
			{
				print $out "\t\t";
			}
		
			if (defined $data[9])
			{
				for (my $i = 7; $i < 10; $i++)
				{
					print $out "\t $data[$i]";
				}
		
				print $out "\t$data[12]\t$data[13]";
			}
			print $out "\n";
		}
		close $in;
		close $out;
	}
	else
	{ print "$filename does not exists\n";}
	print " Finished working on $filename\n";


	# 	_peaks.overlap files
	$experiment = 'peaks.overlap.' . "$type";
	$filename = $peak_annot_dir . $prefix . '_' . $experiment;

	if (-f $filename)
	{
		open (my $in, '<', "$filename") or die "Cannot open $filename\n";
		$newfilename = $filename . '.annotated';
		open (my $out, '>', $newfilename) or die "Cannot open $newfilename\n";
	
		while (<$in>)
		{
			chomp;
			next if /^\s*$/;
			if (/^Chromosome/)
			{
				print $out "$_\n";
				next;
			}
		
			my @data = split /\t/;
		
			$data[4] = genesymbol($data[3]) if (defined $data[3]);
		
			print $out $data[0];

			for (my $i = 1; $i < @data; $i++)
			{
				print $out "\t$data[$i]";
			}
		
			print $out "\n";
		}
		close $in;
		close $out;
	}
	else
	{ print "$filename does not exists\n";}
	print " Finished working on $filename\n";


	# 	_peaks.tss files
	$experiment = 'peaks.tss.' . "$type";
	$filename = $peak_annot_dir . $prefix . '_' . $experiment;

	if (-f $filename)
	{
		open (my $in, '<', "$filename") or die "Cannot open $filename\n";
		$newfilename = $filename . '.annotated';
		open (my $out, '>', $newfilename) or die "Cannot open $newfilename\n";
	
		while (<$in>)
		{
			chomp;
			next if /^\s*$/;
			if (/^Chromosome/)
			{
				print $out "$_\tPot_Assoc_Prom\tPot_Assoc_Enh\n";
				next;
			}
		
			my @data = split /\t/;
		
			$data[7] = genesymbol($data[6]) if (defined $data[6]);
		
			if    (abs($data[3]) < 1000)		{	$data[9] = 1; $data[10] = 0;	}
			elsif (abs($data[3]) < 10000 )		{	$data[9] = 0; $data[10] = 0;	}
			elsif (abs($data[3]) > 10000 &&
			       abs($data[3]) < 1000000)		{	$data[9] = 0; $data[10] = 1;	}
			else					{	$data[9] = 0; $data[10] = 0;	}
		
			print $out $data[0];

			for (my $i = 1; $i < @data; $i++)
			{
				print $out "\t$data[$i]";
			}
		
			print $out "\n";
		}
		close $in;
		close $out;
	}
	else
	{ print "$filename does not exists\n";}
	print " Finished working on $filename\n";
}

