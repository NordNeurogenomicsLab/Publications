#!/usr/bin/perl
# ChIPseq_wrapper.pl by Rinaldo Catta-Preta @UCDavis
# License terms (GPL v.3) are described in https://github.com/NordNeurogenomicsLab/Publications/blob/master/Nature%20Neuroscience%202017/LICENSE

use strict;
use warnings 'FATAL' => 'all';
use Env;
use POSIX;
use List::MoreUtils qw(uniq);
use Data::Dumper;
 
die "\n*********
usage: $0 <editable_parameters_file>
*********


" unless @ARGV == 1; 
my ($edit_file) = @ARGV;
my $today = `date "+%Y-%m-%d"`;
chomp($today);

# General variables
my ($file_list, $raw_data_dir, $processed_dir, $genome_dir, $organism, $genome_version) = '';
my ($genome, $bed_genome, $chrom_sizes, $genes_gtf, $genes_bed, $blacklist_file) = '';
my ($mem_fastqc, $mem_bwa, $mem_macs, $mem_bedtools, $mem_chipqc, $new_peak_call) = '';
my ($gcpart_fastq, $gcpart_bwa, $gcpart_macs, $gcpart_bed, $outfile, $macs_date) = '';

# FASTQC variables
my ($fastqc_run, $fastqc_re_run) = '';

# BWA variables
my ($bwa_dir, $bwa_realign, $bwa_threads, $bwa_quality, $bwa_trim, $bwa_align_modsize) = '';

# MACS2 variables
my ($macs_dir, $macs_run, $macs_broad_cutoff, $macs_qvalue, $macs_pvalue) = '';
my ($macs_nolambda, $macs_extsize, $macs_dups) = '';

# MACS2 bdgdiff variables
my ($diff_dir, $diff_run, $diff_re_run, $diff_cutoff, $diff_minlen, $diff_maxgap) = '';
my ($diff_depth1, $diff_depth2) = '';

# ChIPQC variables
my ($chipqc_dir, $chipqc_re_run, $chipqc_overlap_p, $chipqc_pkcaller, $chipqc_colorby) = '';
my ($chipqc_facetby) = '';

# BedTools variables
my ($bedtools_run, $bedtools_re_run) = '';

# DeepTools variables
my ($deep_dir, $deeptools_run, $deeptools_re_run, $deep_comparisons, $deep_pseudonumber) = '';
my ($deep_ratio_type, $deep_ign_dup, $deep_normal) = '';

# hiddenDomains variables
my ($hidden_code, $hidden_dir, $hidden_run, $hidden_pvalue, $hidden_ext) = '';

# PeakAnnotator variables
my ($peakannot_code, $peakannot_run, $peakannot_re_run) = '';

# Homer variables
my ($homer_dir, $homer_re_run, $homer_date, $homer_pk_caller, $homer_size_up, $homer_run) = '';
my ($homer_size_down, $homer_bkgrnd, $homer_rna, $homer_known, $homer_norm, $homer_mask) = '';
my ($homer_scoring, $homer_bin) = '';

# ChIPModule variables
my ($chipmod_run, $chipmod_re_run, $chipmod_lambda, $chipmod_pvalue, $chipmod_format) = '';
my ($chipmod_detailed, $chipmod_code, $chipmod_dir) = '';

# Trackhub variables
my ($trackhub_run, $trackhub_proj, $trackhub_bio_afs, $trackhub_name) = '';
my ($trackhub_pcreator, $chipmod_file, $chipmod_PWMs, $chipmod_seq) = '';
my ($trackhub_bio_url, $trackhub_afs_url, $trackhub_autoScale) = '';

my $version = '';
my @version = ();

system("tr '\r' '\n' < $edit_file > nullfile.txt") != -1 or die $!;
system("mv nullfile.txt $edit_file") != -1 or die $!;

open (my $in, '<' , $edit_file) or die "Can't open $edit_file\n";
print "\nNord Lab ChIP-seq analysis pipeline version 1.9.0 running ...\n\n";
my $now = `date`; 
print "... started running on $now\n\n";
while (<$in>)
{
	chomp;
	print "$_\n" unless (/^#[^#]/);
	next unless (/^\>/);
	my @temp = split / /;
	my ($header, $param) = ($temp[1], $temp[2]);
	
	if ($header =~ m/Dataset/)				{ $file_list = $param; } 
	elsif ($header =~ m/Reference_organism/) 		{ $organism = $param; }	
	elsif ($header =~ m/Genome_version/) 			{ $genome_version = $param; }
	elsif ($header =~ m/Reference_Genome/) 			{ $genome = $param; }
	elsif ($header =~ m/Reference_genome_bed_file/) 	{ $bed_genome = $param; }
	elsif ($header =~ m/Reference_genome_chrom_sizes/) 	{ $chrom_sizes = $param; }
	elsif ($header =~ m/Refer_genome_chrom_directory/)	{ $genome_dir = $param; }
	elsif ($header =~ m/Blacklist_file/) 			{ $blacklist_file = $param; }
	elsif ($header =~ m/Genes_GTF_file/) 			{ $genes_gtf = $param; }
	elsif ($header =~ m/Genes_BED_file/) 			{ $genes_bed = $param; }
	elsif ($header =~ m/PeakAnnotator_Directory/) 		{ $peakannot_code = $param; }
	elsif ($header =~ m/Trackhubs_project_code_location/)	{ $trackhub_pcreator = $param; }
	elsif ($header =~ m/ChIPModule_code_location/)		{ $chipmod_code = $param; }
	elsif ($header =~ m/hiddenDomains_code/)		{ $hidden_code = $param; }
	elsif ($header =~ m/Raw_Data_Directory/) 		{ $raw_data_dir = $param; }
	elsif ($header =~ m/Processed_File_Directory/) 		{ $processed_dir = $param; }
	
	elsif ($header =~ m/BWA_Data_Directory/) 		{ $bwa_dir = $param; }
	elsif ($header =~ m/Realignment_Required/)		{ $bwa_realign = $param; }
	elsif ($header =~ m/Trim_read_files/) 			{ $bwa_trim = $param; }
	elsif ($header =~ m/Size_of_each_alignment_module/) 	{ $bwa_align_modsize = $param; }
	elsif ($header =~ m/Quality_cutoff/) 			{ $bwa_quality = $param; }
	elsif ($header =~ m/No_of_threads_used_in_alignment/)	{ $bwa_threads = $param; }
	
	elsif ($header =~ m/New_Peak_Calls/) 			{ $new_peak_call = $param; }
	elsif ($header =~ m/MACS2_Data_Directory/) 		{ $macs_dir = $param; }
	elsif ($header =~ m/Run_MACS2_in_the_pipeline/) 	{ $macs_run = $param; }
	elsif ($header =~ m/Broad_peak_calling_cutoff/) 	{ $macs_broad_cutoff = $param; }	
	elsif ($header =~ m/q-value_cutoff/) 			{ $macs_qvalue = $param; }	
	elsif ($header =~ m/p-value_cutoff/) 			{ $macs_pvalue = $param; }
	elsif ($header =~ m/extsize_parameter/) 		{ $macs_extsize = $param; }
	elsif ($header =~ m/use_background_as_local_lambda/) 	{ $macs_nolambda = $param; }
	elsif ($header =~ m/No_of_dups_to_keep/) 		{ $macs_dups = $param; }	
	
	elsif ($header =~ m/BDGdiff_Data_Directory/)		{ $diff_dir = $param; }
	elsif ($header =~ m/Run_MACS2_bdgdiff_in_the_pipeline/)	{ $diff_run = $param; }
	elsif ($header =~ m/Re-run_bdgdiff_on_called_peaks/)	{ $diff_re_run = $param; }
	elsif ($header =~ m/logLR_cutoff/)			{ $diff_cutoff = $param; }
	elsif ($header =~ m/Min_length_of_diff_region/)		{ $diff_minlen = $param; }
	elsif ($header =~ m/Maximum_gap/)			{ $diff_maxgap = $param; }
	elsif ($header =~ m/Sequ_depth_cond1/)			{ $diff_depth1 = $param; }
	elsif ($header =~ m/Sequ_depth_cond2/)			{ $diff_depth2 = $param; }
	
	elsif ($header =~ m/Run_hiddenDomains_in_the_pipeline/)	{ $hidden_run = $param; }
	elsif ($header =~ m/hiddenDomains_min_post_probab/)	{ $hidden_pvalue = $param; }
	elsif ($header =~ m/hiddenDomains_bin_width/)		{ $hidden_ext = $param; }
	elsif ($header =~ m/hiddenDomains_Data_Directory/)	{ $hidden_dir = $param; }
	
	elsif ($header =~ m/Run_FASTQC_in_the_pipeline/) 	{ $fastqc_run = $param; }
	elsif ($header =~ m/Re-Run_FASTQC/) 			{ $fastqc_re_run = $param; }
	
	elsif ($header =~ m/PeakCaller_file_format/) 		{ $chipqc_pkcaller = $param; }
	elsif ($header =~ m/ChIPQC_Directory/) 			{ $chipqc_dir = $param; }
	elsif ($header =~ m/FacetBy_ChIPQC/) 			{ $chipqc_facetby = $param; }
	elsif ($header =~ m/colourBy/) 				{ $chipqc_colorby = $param; }
	elsif ($header =~ m/No_random_datasets_for_overlap_p/) 	{ $chipqc_overlap_p = $param; }
	elsif ($header =~ m/Re-run_ChIPQC_on_called_peaks/)	{ $chipqc_re_run = $param; }
	
	elsif ($header =~ m/Run_Bedtools_in_the_pipeline/)	{ $bedtools_run = $param; }
	elsif ($header =~ m/Re-run_Bedtools_on_called_peaks/)	{ $bedtools_re_run = $param; }
	elsif ($header =~ m/Run_PeakAnnot_in_the_pipeline/)	{ $peakannot_run = $param; }
	elsif ($header =~ m/Re-run_PeakAnnot_on_called_peaks/)	{ $peakannot_re_run = $param; }
	
	elsif ($header =~ m/Run_DeepTools_in_the_pipeline/)	{ $deeptools_run = $param; }
	elsif ($header =~ m/Re-run_Deeptools/)			{ $deeptools_re_run = $param; }
	elsif ($header =~ m/DeepTools_Data_Directory/)		{ $deep_dir = $param; }
	elsif ($header =~ m/Comparisons_to_make/)		{ $deep_comparisons = $param; }
	elsif ($header =~ m/Pseudocount_value/)			{ $deep_pseudonumber = $param; }
	elsif ($header =~ m/Ratio_type/)			{ $deep_ratio_type = $param; }
	elsif ($header =~ m/ignore_Duplicates/)			{ $deep_ign_dup = $param; }
	elsif ($header =~ m/Normalization_type/)		{ $deep_normal = $param; }
	
	elsif ($header =~ m/Run_Homer_in_the_pipeline/) 	{ $homer_run = $param; }
	elsif ($header =~ m/Homer_Data_Directory/) 		{ $homer_dir = $param; }
	elsif ($header =~ m/Peak_caller/) 			{ $homer_pk_caller = $param; }
	elsif ($header =~ m/Repeat_masked_sequences/) 		{ $homer_mask = $param; }
	elsif ($header =~ m/Region_size_upstream/) 		{ $homer_size_up = $param; }
	elsif ($header =~ m/Region_size_downstream/) 		{ $homer_size_down = $param; }
	elsif ($header =~ m/Remove_background_positions/) 	{ $homer_bkgrnd = $param; }
	elsif ($header =~ m/Look_for_RNA_motiffs/) 		{ $homer_rna = $param; }
	elsif ($header =~ m/File_of_known_motifs(if\ existent)/){ $homer_known = $param; }
	elsif ($header =~ m/Sequence_normalization/) 		{ $homer_norm = $param; }
	elsif ($header =~ m/Bin_size_for_histogram/) 		{ $homer_bin = $param; }
	elsif ($header =~ m/Enrichment_scoring_method/) 	{ $homer_scoring = $param; }
	elsif ($header =~ m/Re-run_Homer_on_called_peaks/) 	{ $homer_re_run = $param; }
	
	elsif ($header =~ m/Date_for_extra_analyses/) 		{ $macs_date = $param; }
	
	elsif ($header =~ m/Amt_memory_chipqc/) 		{ $mem_chipqc = $param; }
	elsif ($header =~ m/Amt_memory_fastqc/) 		{ $mem_fastqc = $param; }
	elsif ($header =~ m/Amt_memory_bwa/) 			{ $mem_bwa = $param; }
	elsif ($header =~ m/Amt_memory_macs/) 			{ $mem_macs = $param; }
	elsif ($header =~ m/Amt_memory_bedtools/) 		{ $mem_bedtools = $param; }
	elsif ($header =~ m/GC_partition_fastqc/) 		{ $gcpart_fastq = $param; }
	elsif ($header =~ m/GC_partition_bwa/) 			{ $gcpart_bwa = $param; }
	elsif ($header =~ m/GC_partition_macs/) 		{ $gcpart_macs = $param; }
	elsif ($header =~ m/GC_partition_bedtools/) 		{ $gcpart_bed = $param; }
} 
close $in;

print "\nFiles to be analyzed:\n";
my %ChIP_seq = ();
my @ChIP_list = ();
my @controls = ();
my @PE_controls = ();
my @ChIP_PE_list = ();
my @ChIP_R2_list = ();
my @header;
my ($index_target, $index_fastq, $index_factor, $index_type, $index_experim);

system("tr '\r' '\n' < $file_list > nullfile.txt") != -1 or die $!;
system("mv nullfile.txt $file_list") != -1 or die $!;

open ($in, '<', "$file_list") or die "Cannot open $file_list\n";

while (<$in>)
{
	chomp;
	next if /^\s*$/;
	next if /^#/;
	next if /^"/;
	if (/Fastq_name/)
	{
		@header = split /\t/;
		($index_target) = grep { $header[$_] eq 'target' } 0..$#header;
		($index_fastq) = grep { $header[$_] eq 'Fastq_name' } 0..$#header;
		($index_factor) = grep { $header[$_] eq 'factor' } 0..$#header;
		($index_experim) = grep { $header[$_] eq 'experiment' } 0..$#header;
	}
	next if /Fastq_name/;
	my @data = split /\t/;
	
	for (my $i = 0; $i < @data; $i++)
	{
		if ($i ne $index_target and $i ne $index_fastq)
		{
			$ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{$header[$i]} = $data[$i];
		}
	}
	
	if ($ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{pair_ended} eq 'no' or
		$ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{pair_ended} eq '')
	{
		push @controls, [$data[$index_experim], $data[$index_fastq]] if 
			($ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{type} eq 'Input' or 
			 $ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{type} eq 'NegCtl' );

		unless (grep {$_ eq $data[$index_fastq]} @ChIP_list)
		{
			push @ChIP_list, $data[$index_fastq];
			print "Single-ended: $data[$index_fastq]\n";
		}
	}
	else
	{
		push @PE_controls, [$data[$index_experim], $data[$index_fastq]] if 
			($ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{type} eq 'Input' or 
			 $ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{type} eq 'NegCtl' );

		if ($ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{otherID} eq 'R1')
		{
			push @ChIP_PE_list, $data[$index_fastq];
		
			$ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{PE_trim} = 
						$data[$index_target] . '_' .
						$ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{tissue} . '_' .
						$ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{condition} . '_' . 
						$ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{replicate} . '_';
		
			$ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{PE_notrim} = 
						$ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{PE_trim};
		
			$ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{PE_trim} .= 
						$ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{type} . "_trimmed_PE.bam";
		
			$ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{PE_notrim} .= 	
						$ChIP_seq{$data[$index_target]}{$data[$index_fastq]}{type} . "_PE.bam";
			print "Pair-ended: $data[$index_fastq]\n";
		}
		else
		{
			push @ChIP_R2_list, $data[$index_fastq];
		}
	}
}
close $in;

my $slurm_array_id = '';
my $slurm_array_id2 = '';
my $slurm_array_id3 = '';
my $slurm_job_id = '';
my $myqueue = '';

sub submit_slurm_job
{
	my @object = @_;
	@object = @object[3..$#object];
	my $myjob = `sbatch --time=1-0 $_[1] $_[0] $_[2] @object`;
	($myjob =~ /^Submitted batch job (\d+)/) or die "Error executing $myjob";
	$myjob = $1;
	my @sbatch = split /\.sbatch/, $_[0];
	print "\nSLURM $sbatch[0] job is no. $myjob\n";
	$version = `grep 'version' $_[0]`;	
	@version = split('version', $version);
	print "$_[0] version is $version[1]\n";
	
	return $myjob;
}

sub submit_slurm_array
{
	my @object = @_;
	@object = @object[3..$#object];
	my $myarray;
	my $exclusion = '';
	if ($_[0] =~ m/MACS_diff/) { $myarray = `sbatch $exclusion --time=7-0 $_[1] $_[0] $_[2] @object`; }
	else  { $myarray = `sbatch $exclusion --time=1-0 $_[1] $_[0] $_[2] @object`; }
	($myarray =~ /^Submitted batch job (\d+)/) or die "Error executing $myarray";
	$myarray = $1;
	my @sbatch = split /\.sbatch/, $_[0];
	print "\nSLURM $sbatch[0] array is no. $myarray\n";
	$version = `grep 'version' $_[0]`;	
	@version = split('version', $version);
	print "$_[0] version is $version[1]\n";
	
	return $myarray;
}

sub update_myqueue
{
	$myqueue = "$slurm_array_id"  if ($slurm_array_id ne '');
	$myqueue = "$slurm_array_id2" if ($slurm_array_id2 ne '');
	$myqueue = "$slurm_array_id3" if ($slurm_array_id3 ne '');
	$myqueue = "$slurm_array_id" . ':' . "$slurm_array_id2" if ($slurm_array_id ne '' and $slurm_array_id2 ne '');
	$myqueue = "$slurm_array_id" . ':' . "$slurm_array_id3" if ($slurm_array_id ne '' and $slurm_array_id3 ne '');
	$myqueue = "$slurm_array_id2" . ':' . "$slurm_array_id3" if ($slurm_array_id2 ne '' and $slurm_array_id3 ne '');
	$myqueue = "$slurm_array_id" . ':' . "$slurm_array_id2" . ':' . "$slurm_array_id3" 
				if ($slurm_array_id ne '' and $slurm_array_id2 ne '' and $slurm_array_id3 ne '');
}

# Create destination directories if they do not exist
system("mkdir -p $processed_dir") != -1 or die $!;

$bwa_dir  = $processed_dir . "bwa/"   unless (defined $bwa_dir);
system("mkdir -p $bwa_dir") != -1 or die $!;

$macs_dir = $processed_dir . "macs2/" unless (defined $macs_dir);
my $code_dir = $processed_dir . "code/";

my $logfile_dir;

my $no_chips = scalar(@ChIP_list);
my $no_PE_chips = scalar(@ChIP_PE_list);
die "ChIP-seq list to be analyzed is empty\n" unless ($no_chips > 0 or $no_PE_chips > 0);
$no_chips = $no_chips - 1;
$no_PE_chips = $no_chips - 1;

### SEQUENCE ALIGNMENT
# Using BWA  

my @ChIP_align_list = ();
my @ChIP_PE_align_list = ();
my @ChIP_fastqc_list = ();

# check whether .bam files from @ChIP_list are present
my $fastqc_dir = $processed_dir . "fastqc/";

foreach my $chip (@ChIP_list)
{
	my @ChIP_bam_name = split(/\.fastq\.gz/, $chip);
	my $ChIP_bam;
	my $ChIP_fastqc = $fastqc_dir . $ChIP_bam_name[0] . '_fastqc.html';
	
	if ($bwa_trim eq 'yes') { $ChIP_bam = $bwa_dir . $ChIP_bam_name[0] . '_trimmed.fq.gz.srt.bam';}
	else			{ $ChIP_bam = $bwa_dir . $chip . '.fastq.gz.srt.bam'; }
	
	if (! -e $ChIP_bam)	{ push(@ChIP_align_list, $chip); }
	else			{ push(@ChIP_align_list, $chip) unless ($bwa_realign ne 'yes'); }
	
	if (! -e $ChIP_fastqc)	{ push(@ChIP_fastqc_list, $chip); }
	else			{ push(@ChIP_fastqc_list, $chip) unless ($fastqc_re_run ne 'yes'); }
}

foreach my $chip (@ChIP_PE_list)
{
	my @ChIP_bam_name = split(/\.fastq\.gz/, $chip);
	my $ChIP_fastqc = $fastqc_dir . $ChIP_bam_name[0] . '_fastqc.html';

	if (! -e $ChIP_fastqc)	{ push(@ChIP_fastqc_list, $chip); }
	else			{ push(@ChIP_fastqc_list, $chip) unless ($fastqc_re_run ne 'yes'); }
}
	
foreach my $chip (@ChIP_R2_list)
{
	my @ChIP_bam_name = split(/\.fastq\.gz/, $chip);
	my $ChIP_bam;
	my $ChIP_fastqc = $fastqc_dir . $ChIP_bam_name[0] . '_fastqc.html';
	
	if (! -e $ChIP_fastqc)	{ push(@ChIP_fastqc_list, $chip); }
	else			{ push(@ChIP_fastqc_list, $chip) unless ($fastqc_re_run ne 'yes'); }
}
	
foreach my $PE_target (keys %ChIP_seq)
{
	my $ChIP_bam;
	
	foreach my $fastq (keys %{ $ChIP_seq{$PE_target} })
	{
		if (grep {$_ eq $fastq} @ChIP_PE_list)
		{
			if ($bwa_trim eq 'yes') { $ChIP_bam = $bwa_dir . $ChIP_seq{$PE_target}{$fastq}{PE_trim};}
			else			{ $ChIP_bam = $bwa_dir . $ChIP_seq{$PE_target}{$fastq}{PE_notrim}; }
			
			if (! -e $ChIP_bam)	{ push(@ChIP_PE_align_list, $fastq); }
			else			{ push(@ChIP_PE_align_list, $fastq) unless ($bwa_realign ne 'yes'); }
		}
	}
}

my $pre_param;
my $pos_param;
my $no_align_chips = scalar(@ChIP_align_list) - 1;
my $no_align_PE_chips = scalar(@ChIP_PE_align_list) - 1;
my $low_threads;

if ($bwa_threads > 2)	{ $low_threads = int($bwa_threads / 2); }
else			{ $low_threads = $bwa_threads; }

my $fastqc_to_run = @ChIP_fastqc_list;
$fastqc_to_run = 0 if ($fastqc_run ne 'yes');

# Quality control on the raw reads
$logfile_dir = $fastqc_dir . "logfiles/";
system("mkdir -p $logfile_dir") != -1 or die $!;

my $logfstqcout = $logfile_dir . "FASTQC_%A_%a.out";

if ($fastqc_to_run > 0)
{
	system("mkdir -p $fastqc_dir") != -1 or die $!; 
	
	my $no_fastqc_chips = scalar(@ChIP_fastqc_list) - 1;
		
	$pre_param = "--array=0-$no_fastqc_chips --partition=$gcpart_fastq --cpus-per-task=$low_threads ";
	$pre_param.= "--mem=$mem_fastqc -o $logfstqcout -e $logfstqcout";
	
	$pos_param = "$low_threads $raw_data_dir $fastqc_dir";
	
	$slurm_array_id = submit_slurm_array("FASTQC.sbatch", $pre_param, $pos_param, @ChIP_fastqc_list);
	update_myqueue;
}
else
{
	print "No FASTQC analysis selected\n";
}

# Single-ended read alignment
if ($no_align_chips >= 0)
{
	$logfile_dir = $bwa_dir . "logfiles/";
	system("mkdir -p $logfile_dir") != -1 or die $!;

	my $logbwaout = $logfile_dir . "BWA_%A_%a.out";
	
	# split @ChIP_align_list into lower-number-of-element lists for storage management
	my $slice = $bwa_align_modsize;
	$slice = $no_align_chips + 1 if ($bwa_align_modsize == 0);
	
	for (my $j = 0; $j <= $no_align_chips; $j+= $slice)
	{
		my $k = $j + $slice - 1;
		$k = $no_align_chips if ($k >= $no_align_chips);
		my $remaining = $slice - 1;
		$remaining = ($k - $j) if ($k == $no_align_chips);
		
		$pre_param = "--array=0-$remaining --partition=$gcpart_bwa --cpus-per-task=$bwa_threads";
		$pre_param.= " --mem=$mem_bwa -o $logbwaout -e $logbwaout";
		$pre_param.= " --dependency=afterok:$myqueue" if ($myqueue ne '');
		
		$pos_param = "$bwa_threads $genome $raw_data_dir $bwa_dir $bwa_quality $bwa_trim $fastqc_dir ";
		
		$slurm_array_id = submit_slurm_array("BWA.sbatch", $pre_param, $pos_param, @ChIP_align_list[$j..$k]);
		update_myqueue;
	}
}
else
{
	print "No BWA alignment of single-ended reads selected\n\n";
}

# Pair-ended read alignment
if ($no_align_PE_chips >= 0)
{
	my @ChIP_PE_align = ();
	
	foreach my $PE_target (keys %ChIP_seq)
	{
		foreach my $PE_R1 (@ChIP_PE_align_list)
		{
			next unless (exists $ChIP_seq{$PE_target}{$PE_R1}{type});
			foreach my $PE_R2 (keys %{ $ChIP_seq{$PE_target} })
			{
				next unless (exists $ChIP_seq{$PE_target}{$PE_R2}{type});
				next if ($PE_R1 eq $PE_R2);

				if ($ChIP_seq{$PE_target}{$PE_R2}{experiment} eq $ChIP_seq{$PE_target}{$PE_R1}{experiment} and
					$ChIP_seq{$PE_target}{$PE_R2}{type}		  eq $ChIP_seq{$PE_target}{$PE_R1}{type} and
					$ChIP_seq{$PE_target}{$PE_R2}{tissue}	  eq $ChIP_seq{$PE_target}{$PE_R1}{tissue} and
					$ChIP_seq{$PE_target}{$PE_R2}{factor}	  eq $ChIP_seq{$PE_target}{$PE_R1}{factor} and
					$ChIP_seq{$PE_target}{$PE_R2}{condition}  eq $ChIP_seq{$PE_target}{$PE_R1}{condition} and
					$ChIP_seq{$PE_target}{$PE_R2}{treatment}  eq $ChIP_seq{$PE_target}{$PE_R1}{treatment} and
					$ChIP_seq{$PE_target}{$PE_R2}{replicate}  eq $ChIP_seq{$PE_target}{$PE_R1}{replicate} and
					$ChIP_seq{$PE_target}{$PE_R2}{pair_ended} eq $ChIP_seq{$PE_target}{$PE_R1}{pair_ended})
				{
					if ($bwa_trim eq 'yes') { $outfile = substr($ChIP_seq{$PE_target}{$PE_R1}{PE_trim}, 0, -14); }
					else				{ $outfile = substr($ChIP_seq{$PE_target}{$PE_R1}{PE_notrim}, 0, -17); }
					
					push @ChIP_PE_align, $PE_R1;
					push @ChIP_PE_align, $PE_R2;
					push @ChIP_PE_align, $outfile;
					print "R1: $PE_R1\nR2: $PE_R2\noutfile: $outfile\n\n";
				}
			}
		}
	}
	
	die "... Not all reads marked as paired-ended are in pairs" if (scalar(@ChIP_PE_align) ne (3 * ($no_align_PE_chips + 1)));
	
	my $logfile_dir = $bwa_dir . "logfiles/";
	system("mkdir -p $logfile_dir") != -1 or die $!;

	my $logbwaout = $logfile_dir . "BWA_PE_%A_%a.out";
	
	$pre_param = "--array=0-$no_align_PE_chips --partition=$gcpart_bwa --cpus-per-task=$bwa_threads";
	$pre_param.= "--mem=$mem_bwa --exclusive -o $logbwaout -e $logbwaout";
	$pre_param.= " --dependency=afterok:$myqueue" if ($myqueue ne '');
	$pos_param = "$bwa_threads $genome $raw_data_dir $bwa_dir $bwa_quality $bwa_trim $fastqc_dir";

	$slurm_array_id3 = submit_slurm_array("BWA_PE.sbatch", $pre_param, $pos_param, @ChIP_PE_align);
	update_myqueue;
}
else { print "No BWA alignment of pair-ended reads selected\n\n"; }

my ($broad_dir, $narrow_dir) = '';
$today = `date "+%Y-%m-%d"`;
chomp($today);
$broad_dir = $macs_dir . 'broad/' . $today . '/';
$narrow_dir = $macs_dir . 'narrow/' . $today . '/';


# Select samples for peak calling
my @macs_list = ();
my @post_fastqc_list = ();
my $controls = scalar(@controls);

foreach my $target (keys %ChIP_seq)
{
	# Create list for single-ended reads
	foreach my $chip (@ChIP_list)
	{
		next unless (exists $ChIP_seq{$target}{$chip}{type});
		next unless ($ChIP_seq{$target}{$chip}{type} eq 'ChIP');
		my $newchip = substr($chip, 0, -9);
		my $trtfile;
		
		if ($controls != 0)
		{
			for my $element ( @controls )
			{
				my $newctrlchip = substr(@{$element}[1], 0, -9);
				my $ctrlfile;
				if (@{$element}[0] eq $ChIP_seq{$target}{$chip}{experiment})
				{
					$outfile = $target . '_' . $ChIP_seq{$target}{$chip}{tissue} . '_' .
							   $ChIP_seq{$target}{$chip}{condition};
					$outfile.= '_' . $ChIP_seq{$target}{$chip}{otherID} if ($ChIP_seq{$target}{$chip}{otherID} ne '');
					$outfile.= '_' . $ChIP_seq{$target}{$chip}{replicate};
			
					if ($bwa_trim eq 'yes')
					{
						$trtfile  = $newchip . "_trimmed.fq.gz.srt.bam";
						$ctrlfile = $newctrlchip . "_trimmed.fq.gz.srt.bam";
						$outfile.=  '_trim-vs.' . $ChIP_seq{$target}{@{$element}[1]}{type};
					}
					else
					{
						$trtfile  = $chip . ".fastq.gz.srt.bam";
						$ctrlfile = @{$element}[1] . ".fastq.gz.srt.bam";
						$outfile.=  '_notrim-vs.' . $ChIP_seq{$target}{@{$element}[1]}{type};
					}				

					$ChIP_seq{$target}{@{$element}[1]}{trim_file} = $ctrlfile;
					$ChIP_seq{$target}{$chip}{trim_file} = $trtfile;
					$trtfile = $bwa_dir . $trtfile;
					$ctrlfile = $bwa_dir . $ctrlfile;
					$ChIP_seq{$target}{$chip}{i_macs} = $narrow_dir . $outfile if ($ChIP_seq{$target}{@{$element}[1]}{type} eq 'Input');
					$ChIP_seq{$target}{$chip}{n_macs} = $narrow_dir . $outfile if ($ChIP_seq{$target}{@{$element}[1]}{type} eq 'NegCtl');
					push @macs_list, $trtfile;
					push @macs_list, $ctrlfile;
					push @macs_list, $outfile;
					push @macs_list, $ChIP_seq{$target}{$chip}{sple_info};
					push @macs_list, $ChIP_seq{$target}{@{$element}[1]}{sple_info};
					my @name = split /\//, $trtfile;
					print "$name[-1] vs. $ChIP_seq{$target}{@{$element}[1]}{type}\n";
				
					push @post_fastqc_list, $trtfile unless (grep {$_ eq $trtfile} @post_fastqc_list);
					push @post_fastqc_list, $ctrlfile unless (grep {$_ eq $ctrlfile} @post_fastqc_list);
				}
			}
		}
		else 
		{
			$outfile = $target . '_' . $ChIP_seq{$target}{$chip}{tissue} . '_' .
						   $ChIP_seq{$target}{$chip}{condition};
			$outfile.= '_' . $ChIP_seq{$target}{$chip}{otherID} if ($ChIP_seq{$target}{$chip}{otherID} ne '');
			$outfile.= '_' . $ChIP_seq{$target}{$chip}{replicate};
				
			if ($bwa_trim eq 'yes')
			{
				$trtfile  = $newchip . "_trimmed.fq.gz.srt.bam";
				$outfile.=  '_trim';
			}
			else
			{
				$trtfile  = $chip . ".fastq.gz.srt.bam";
				$outfile.=  '_notrim';
			}				

			$ChIP_seq{$target}{$chip}{trim_file} = $trtfile;
			$trtfile = $bwa_dir . $trtfile;
			push @macs_list, $trtfile;
			push @macs_list, 'null';
			push @macs_list, $outfile;
			push @macs_list, 'null';
			push @macs_list, 'null';
		}
	}
	
	# Append pair-ended list to the single-ended one
	foreach my $chip (@ChIP_PE_list)
	{
		next unless (exists $ChIP_seq{$target}{$chip}{type});
		next unless ($ChIP_seq{$target}{$chip}{type} eq 'ChIP');
		
		for my $element ( @PE_controls )
		{
			if (@{$element}[0] eq $ChIP_seq{$target}{$chip}{experiment} and $ChIP_seq{$target}{@{$element}[1]}{otherID} eq 'R1')
			{
				my $newchip = substr($ChIP_seq{$target}{$chip}{PE_trim}, 0, -25);
				my $trtfile;
				my $ctrlfile;
				my $outfile;
					
				if ($bwa_trim eq 'yes')
				{
					$trtfile  = $ChIP_seq{$target}{$chip}{PE_trim};
					$ctrlfile = $ChIP_seq{$target}{@{$element}[1]}{PE_trim};
					$outfile = $newchip . '_trim_PE-vs.' . $ChIP_seq{$target}{@{$element}[1]}{type};
				}
				else
				{
					$trtfile  = $ChIP_seq{$target}{$chip}{PE_notrim};
					$ctrlfile = $ChIP_seq{$target}{@{$element}[1]}{PE_notrim};
					$outfile = $newchip . '_notrim_PE-vs.' . $ChIP_seq{$target}{@{$element}[1]}{type};
				}
				
				$ChIP_seq{$target}{@{$element}[1]}{trim_file} = $ctrlfile;
				$ChIP_seq{$target}{$chip}{trim_file} = $trtfile;
				$trtfile = $bwa_dir . $trtfile;
				$ctrlfile = $bwa_dir . $ctrlfile;
				$ChIP_seq{$target}{$chip}{i_macs} = $narrow_dir . $outfile if ($ChIP_seq{$target}{@{$element}[1]}{type} eq 'Input');
				$ChIP_seq{$target}{$chip}{n_macs} = $narrow_dir . $outfile if ($ChIP_seq{$target}{@{$element}[1]}{type} eq 'NegCtl');
				$ChIP_seq{$target}{$chip}{sple_info} = 'anytext' unless (exists $ChIP_seq{$target}{$chip}{sple_info} and $ChIP_seq{$target}{$chip}{sple_info} eq '');
				$ChIP_seq{$target}{@{$element}[1]}{sple_info} = 'anytext' unless (exists $ChIP_seq{$target}{@{$element}[1]}{sple_info} and $ChIP_seq{$target}{@{$element}[1]}{sple_info} eq '');
				push @macs_list, $trtfile;
				push @macs_list, $ctrlfile;
				push @macs_list, $outfile;
				push @macs_list, $ChIP_seq{$target}{$chip}{sple_info};
				push @macs_list, $ChIP_seq{$target}{@{$element}[1]}{sple_info};
				my @name = split /\//, $trtfile;
				print "$name[-1] vs. $ChIP_seq{$target}{@{$element}[1]}{type}\n";
				
				push @post_fastqc_list, $trtfile unless (grep {$_ eq $trtfile} @post_fastqc_list);
				push @post_fastqc_list, $ctrlfile unless (grep {$_ eq $ctrlfile} @post_fastqc_list);
			}
		}
	}
}

# Quality check on samples after trimming and filtering
if ($fastqc_to_run > 0)
{
	my $uplim = @post_fastqc_list - 1;
	$pre_param = "--array=0-$uplim --partition=$gcpart_fastq --cpus-per-task=$low_threads ";
	$pre_param.= "--mem=$mem_fastqc -o $logfstqcout -e $logfstqcout --dependency=afterok:$myqueue";
	$pos_param = "$low_threads $bwa_dir $fastqc_dir";
	submit_slurm_array("FASTQC_post.sbatch", $pre_param, $pos_param, @ChIP_fastqc_list);
}

print "\n";

my $no_macs = ceil(scalar(@macs_list) / 5) - 1;

### PEAK CALLING

$logfile_dir = $macs_dir . "logfiles/";
system("mkdir -p $logfile_dir") != -1 or die $!;

my $logbedout = $logfile_dir . "BedTools_%A_%a.out";

my $peakannot_dir = $processed_dir . "peak_annotation/";

# Using MACS2
if ($macs_run eq 'yes')
{
	if ($myqueue ne '' or $new_peak_call eq 'yes')
	{
		print "\n\nRunning MACS2 on ChIP-seq samples\n";
		
		system("mkdir -p $macs_dir") != -1 or die $!;  
		system("mkdir -p $broad_dir") != -1 or die $!;
		system("mkdir -p $narrow_dir") != -1 or die $!;

		$logfile_dir = $macs_dir . "logfiles/";
		system("mkdir -p $logfile_dir") != -1 or die $!;

		my $logmacsout = $logfile_dir . "MACS2_%A_%a.out";

		$macs_nolambda = '--nolambda' if ($macs_nolambda eq 'yes');
		$macs_pvalue = 0 unless (defined $macs_pvalue);

		$pre_param = "--array=0-$no_macs --partition=$gcpart_macs --mem=$mem_macs --cpus-per-task=4 ";
		$pre_param.= "-o $logmacsout -e $logmacsout";
		$pos_param = "$broad_dir $narrow_dir $organism $macs_broad_cutoff $macs_qvalue $macs_pvalue $macs_nolambda ";
		$pos_param.= "$macs_extsize $macs_dups $mem_macs $logfile_dir";

		$pre_param.= " --dependency=afterok:$myqueue" if ($myqueue ne '');

		$slurm_array_id2 = submit_slurm_array("MACS_exp.sbatch", $pre_param, $pos_param, @macs_list) if
							($controls > 0);
		update_myqueue;
		$pre_param.= " --dependency=afterok:$myqueue";
		$slurm_array_id = submit_slurm_array("MACS_sple.sbatch", $pre_param, $pos_param, @macs_list);
		update_myqueue;
	}
	else { print "No MACS2 peak calling selected\n"; }	
}
else { print "No MACS2 peak calling selected\n"; }

my $new_narrow_dir;
my $new_broad_dir;

# Using deeptools and hiddenDomains
$deep_dir = $processed_dir . "deeptools/" unless (defined $deep_dir);
system("mkdir -p $deep_dir") != -1 or die $!; 
$hidden_dir = $processed_dir . "hiddenDomains/" unless (defined $hidden_dir);
system("mkdir -p $hidden_dir") != -1 or die $!; 

# Preparing the Comparison Lists
my $chosen_sample;
my @compare_list;
my @compare_labels;
my @complete_ChIP;
push @complete_ChIP, @ChIP_list;
# push @complete_ChIP, @ChIP_R2_list;
my @compare_param = split /,/, $deep_comparisons;
@compare_param = sort {$b <=> $a} @compare_param;
print "\ncomparisons:\n", Dumper(@compare_param);

COMPLIST: for my $comptype (@compare_param)
{
	if ($comptype eq 4)
	{
		for (my $i = 0; $i < @complete_ChIP - 1; $i++)
		{
			for (my $j = $i + 1; $j < @complete_ChIP; $j++)
			{
				my $newprim = $complete_ChIP[$i];
				($newprim =~ s/\.fastq\.gz/_trimmed\.fq\.gz\.srt\.bam/) if ($bwa_trim eq 'yes');
				($newprim =~ s/\.fastq\.gz/\.fastq\.gz\.srt\.bam/) if ($bwa_trim ne 'yes');
				push @compare_list, "$bwa_dir$newprim";
				my $newsec = $complete_ChIP[$j];
				($newsec =~ s/\.fastq\.gz/_trimmed\.fq\.gz\.srt\.bam/) if ($bwa_trim eq 'yes');
				($newsec =~ s/\.fastq\.gz/\.fastq\.gz\.srt\.bam/) if ($bwa_trim ne 'yes');
				push @compare_list, "$bwa_dir$newsec";
				my $name = "$newprim" . '_vs_' . "$newsec";
				$_ = $name;
				$name =~ s/_trimmed\.fq\.gz\.srt\.bam//g;
				$name =~ s/\.fastq\.gz\.srt\.bam//g;
				push @compare_list, $name;
			}
		}
		my @prefixes;
		
		for (my $i = 0; $i < @macs_list; $i+= 5)
		{
			my $prefix = $macs_list[$i + 2];
			$prefix =~ s/\-vs\.Input//;
			push @prefixes, $prefix;
		}
		
		for (my $i = 0; $i < @prefixes - 1; $i++)
		{
			for (my $j = $i + 1; $j < @prefixes; $j++)
			{
				if ($prefixes[$i] =~ /KO/)
				{
					my $prefix = $prefixes[$j] . '_vs_' . $prefixes[$i];
					push @compare_labels, $prefix;
				}
				else
				{
					my $prefix = $prefixes[$i] . '_vs_' . $prefixes[$j];
					push @compare_labels, $prefix;
				}
			}
		}
		last COMPLIST;
	}
	elsif ($comptype eq 1)
	{
		for (my $i = 0; $i < @macs_list; $i+= 5)
		{
			push @compare_list, $macs_list[$i];
			push @compare_list, $macs_list[$i + 1];
			push @compare_list, $macs_list[$i + 2];
		}
	}
	elsif ($comptype eq 2)
	{
		for (my $i = 0; $i < @macs_list; $i+= 5)
		{
			if ($macs_list[$i + 2] =~ /WT/)
			{
				for (my $j = 0; $j < @macs_list; $j+= 5)
				{
					if ($macs_list[$j + 2] =~ /KO/)
					{
						push @compare_list, $macs_list[$i];
						push @compare_list, $macs_list[$j];
						my $name = "$macs_list[$i + 2]" . '_vs.' . "$macs_list[$j + 2]";
						$_ = $name;
						$name =~ s/vs\.Input//g;
						$name =~ s/vs\.NegCtl//g;
						push @compare_list, $name;
					}
				}
			}
		}
	}
	elsif ($comptype =~ /^3\([\S*\+]{6}\S*\)$/)
	{
		chop($comptype);
		my @chosen_sample = split /3\(/, $comptype;
		@chosen_sample = split /\+/, $chosen_sample[1];
		
		foreach my $target (keys %ChIP_seq)
		{
			foreach my $fastq (keys %{ $ChIP_seq{$target} })
			{
				if ( ($ChIP_seq{$target}{$fastq}{experiment} eq $chosen_sample[0]) and 
					 ($ChIP_seq{$target}{$fastq}{factor} 	 eq $chosen_sample[1]) and
					 ($ChIP_seq{$target}{$fastq}{type} 	 eq $chosen_sample[2]) and
					 ($ChIP_seq{$target}{$fastq}{tissue} 	 eq $chosen_sample[3]) and
					 ($ChIP_seq{$target}{$fastq}{condition}  eq $chosen_sample[4]) and
					 ($ChIP_seq{$target}{$fastq}{treatment}  eq $chosen_sample[5]) and
					 ($ChIP_seq{$target}{$fastq}{replicate}  eq $chosen_sample[6]) )
				{
					for (my $i = 0; $i < @macs_list; $i+= 5)
					{
						push @compare_list, $macs_list[$i];
						my $newchosen = $fastq;
						($newchosen =~ s/\.fastq\.gz/_trimmed\.fq\.gz\.srt\.bam/) if ($bwa_trim eq 'yes');
						($newchosen =~ s/\.fastq\.gz/\.fastq\.gz\.srt\.bam/) if ($bwa_trim ne 'yes');
						push @compare_list, $newchosen;
						my $name = "$macs_list[$i]" . '_vs.' . "$newchosen";
						$_ = $name;
						$name =~ s/_trimmed\.fq\.gz\.srt\.bam//g;
						$name =~ s/\.fastq\.gz\.srt\.bam//g;
						push @compare_list, $name;
					}
				}
			}
		}
	}
	else { die "\n\nbam Comparison type is not valid\n\n" }
}

my @ratios = split /\,/, $deep_ratio_type;

if ($organism eq 'hs')	{ $genome_version = 'hg' . $genome_version; }
else			{ $genome_version = $organism . $genome_version; }

if ($deeptools_run eq 'yes' or $hidden_run eq 'yes')
{
	$logfile_dir = $hidden_dir . "logfiles/";
	system("mkdir -p $logfile_dir") != -1 or die $!;

	my $logdeept = $logfile_dir . "hiddenDomains_%A_%a.out";
	
	$pre_param = "--array=0-$no_macs --partition=$gcpart_bwa --mem=$mem_macs ";
	$pre_param.= "-o $logdeept -e $logdeept";
	
	if ($deeptools_re_run eq 'yes' or ! defined $myqueue)
	{
		$new_narrow_dir = $macs_dir . 'narrow/' . $macs_date . '/';
		$new_broad_dir = $macs_dir . 'broad/' . $macs_date . '/';
	}
	else
	{
		$pre_param.= " --dependency=afterok:$myqueue";
		
		$new_narrow_dir = $narrow_dir;
		$new_broad_dir = $broad_dir;
	}
	
	$pos_param = "$hidden_ext $hidden_pvalue $chrom_sizes $hidden_code $hidden_dir";
	submit_slurm_array("hiddenDomains.sbatch", $pre_param, $pos_param, @macs_list) if 
			($hidden_run eq 'yes');
	
	# currently, it works only for single-ended reads
	
	if ($deep_ign_dup eq 'yes')	{ $deep_ign_dup = '--ignoreDuplicates' }
	else 				{ $deep_ign_dup = 'no' }
	
	my @normalizations = split /\,/, $deep_normal;
	
	$logfile_dir = $deep_dir . "logfiles/";
	system("mkdir -p $logfile_dir") != -1 or die $!;
	$logdeept = $logfile_dir . "deeptools_%A_%a.out";
	$pre_param = "--array=0-$no_macs --partition=$gcpart_bwa --mem=$mem_macs ";
	$pre_param.= "-o $logdeept -e $logdeept";
	$pre_param.= " --dependency=afterok:$myqueue" if (defined $myqueue and $myqueue ne '');
	
	for my $ratio (@ratios)
	{
		for my $normaliz (@normalizations)
		{
			$pos_param = "$new_broad_dir $new_narrow_dir $deep_dir $bwa_threads $macs_extsize $deep_ign_dup $normaliz $genome_version";
			$slurm_array_id = submit_slurm_array("bamCoverage.sbatch", $pre_param, $pos_param, @macs_list) if 
					($deeptools_run eq 'yes');

			$pos_param = "$deep_pseudonumber $ratio $bwa_threads $normaliz $deep_ign_dup $deep_dir $genome_version $macs_extsize";
			$slurm_array_id = submit_slurm_array("bamCompare.sbatch", $pre_param, $pos_param, @compare_list) if 
					($deeptools_run eq 'yes');
			update_myqueue;
			$pre_param.= " --dependency=afterok:$myqueue" if (defined $myqueue and $myqueue ne '');
		}
	}
}

# Perform Peak Annotation (PeakAnnotation)
if ($peakannot_run eq 'yes')
{
	system("mkdir -p $peakannot_dir") != -1 or die $!; 
	$logfile_dir = $peakannot_dir . "logfiles/";
	system("mkdir -p $logfile_dir") != -1 or die $!;
	
	my $logpeakannot = $logfile_dir . "PeakAnnot_%A_%a.out";

	$pre_param = "--array=0-$no_macs --partition=$gcpart_bed ";
	$pre_param.= "--mem=$mem_bedtools -o $logpeakannot -e $logpeakannot"; 
	
	unless ($peakannot_re_run eq 'yes' or ! defined $myqueue)
	{
		$peakannot_dir = $peakannot_dir . $today . '/';
		system("mkdir -p $peakannot_dir") != -1 or die $!; 

		if ($myqueue =~ /$slurm_array_id/) 	{ $pre_param.= " --dependency=afterok" . ':' . "$myqueue"; }
		else					{ $pre_param.= " --dependency=afterok" . ':' . "$myqueue" . ':' . "$slurm_array_id"; }
		
		$new_narrow_dir = $narrow_dir;
		$new_broad_dir = $broad_dir;
	}
	else
	{
		if ($slurm_array_id ne '') 	{ $pre_param.= " --dependency=afterok:$slurm_array_id"; }
		
		$new_narrow_dir = $macs_dir . 'narrow/' . $macs_date . '/';
		$new_broad_dir = $macs_dir . 'broad/' . $macs_date . '/';
	}
	
	$pos_param = "$peakannot_code $peakannot_dir $new_narrow_dir $genes_bed $chipqc_overlap_p ";
	$pos_param.= "$new_broad_dir $homer_pk_caller";

	if ($peakannot_re_run eq 'yes' or $myqueue ne '')
	{
		$slurm_array_id3 = submit_slurm_array("PeakAnnotator.sbatch", $pre_param, $pos_param, @macs_list);
		update_myqueue;
	}
	else { print "No PeakAnnotator selected\n"; }
}
else { print "No PeakAnnotator selected\n"; }

# Perform QC on the ChIP-seq peak callings
my $out;
my %targets = ();

# >> create list of files to be qc'd (one for Input and one for NegCtl)
my %factors = ();
my $factorfile;
my $number = 0;
my $newqueue = '';
$chipqc_dir = $processed_dir . "chipqc/" unless (defined $chipqc_dir);

sub generate_chipqc_table
{
	if ($_[0] eq 'Input') 	  { $factorfile = $chipqc_dir . $_[1]. "_i.txt"; }
	elsif ($_[0] eq 'NegCtl') { $factorfile = $chipqc_dir . $_[1]. "_n.txt"; };

	open ($out, '>', "$factorfile") or die "Cannot open $factorfile\n";

	print $out "SampleID\tTissue\tFactor\tCondition\tTreatment\tReplicate\tbamReads\t"; 
	print $out "ControlID\tbamControl\tPeaks\tPeakCaller\n";
	
	foreach my $target (keys %ChIP_seq)
	{
		foreach my $chip (keys %{ $ChIP_seq{$target} })
		{
			next unless (exists $ChIP_seq{$target}{$chip}{type});
			
			if ($_[2] ne 'project')
			{
				next if ($ChIP_seq{$target}{$chip}{factor} ne $_[2]);
			}
			
			next unless ($ChIP_seq{$target}{$chip}{type} eq 'ChIP');
			
			unless ($ChIP_seq{$target}{$chip}{pair_ended} eq 'no' or $ChIP_seq{$target}{$chip}{pair_ended} eq '')
			{
				next if ($ChIP_seq{$target}{$chip}{otherID} ne 'R1');
				
				foreach my $ctrlchip (keys %{ $ChIP_seq{$target} })
				{
					next unless (exists $ChIP_seq{$target}{$ctrlchip}{type} and 
								 		$ChIP_seq{$target}{$ctrlchip}{type} eq $_[0]);
					
					next if ($ChIP_seq{$target}{$ctrlchip}{otherID} 	ne 'R1');
					next if ($ChIP_seq{$target}{$ctrlchip}{experiment} 	ne $ChIP_seq{$target}{$chip}{experiment});
					next if ($ChIP_seq{$target}{$ctrlchip}{tissue} 		ne $ChIP_seq{$target}{$chip}{tissue});
					next if ($ChIP_seq{$target}{$ctrlchip}{factor} 		ne $ChIP_seq{$target}{$chip}{factor});
					next if ($ChIP_seq{$target}{$ctrlchip}{condition} 	ne $ChIP_seq{$target}{$chip}{condition});
					next if ($ChIP_seq{$target}{$ctrlchip}{treatment} 	ne $ChIP_seq{$target}{$chip}{treatment});
					next if ($ChIP_seq{$target}{$ctrlchip}{replicate} 	ne $ChIP_seq{$target}{$chip}{replicate});
					
					my $extension;
					$extension = '.no_model_peaks.narrowPeak' if ($chipqc_pkcaller eq 'narrow');
					$extension = '.no_model_peaks.xls' 	  if ($chipqc_pkcaller eq 'macs');
					$extension = '.bed' 			  if ($chipqc_pkcaller eq 'bed');
					$extension = '.txt' 			  if ($chipqc_pkcaller eq 'raw');
					
					my $Peaks;
					$Peaks = "$ChIP_seq{$target}{$chip}{i_macs}" 
							 if ($ChIP_seq{$target}{$ctrlchip}{type} eq 'Input');
					$Peaks = "$ChIP_seq{$target}{$chip}{n_macs}" 
							 if ($ChIP_seq{$target}{$ctrlchip}{type} eq 'NegCtl');
					
					my @temp_peaks = split /\//, $Peaks;
					
					print $out "$ChIP_seq{$target}{$chip}{sampleID}\t";
					print $out "$ChIP_seq{$target}{$chip}{tissue}\t"; 
					print $out "$ChIP_seq{$target}{$chip}{factor}\t";
					print $out "$ChIP_seq{$target}{$chip}{condition}\t"; 
					print $out "$ChIP_seq{$target}{$chip}{treatment}\t";
					print $out "$ChIP_seq{$target}{$chip}{replicate}\t"; 
					print $out "$bwa_dir$ChIP_seq{$target}{$chip}{PE_trim}\t"   if ($bwa_trim eq 'yes');
					print $out "$bwa_dir$ChIP_seq{$target}{$chip}{PE_notrim}\t" if ($bwa_trim ne 'yes');
					print $out "$ChIP_seq{$target}{$ctrlchip}{sampleID}\t";
					print $out "$bwa_dir$ChIP_seq{$target}{$ctrlchip}{PE_trim}\t"   if ($bwa_trim eq 'yes');
					print $out "$bwa_dir$ChIP_seq{$target}{$ctrlchip}{PE_notrim}\t" if ($bwa_trim ne 'yes');
					print $out "$Peaks$extension\t";
					print $out "$chipqc_pkcaller\n";
					$number += 1;
				}
			}
			else
			{
				foreach my $ctrlchip (keys %{ $ChIP_seq{$target} })
				{
					next unless (exists $ChIP_seq{$target}{$ctrlchip}{type} and 
								 		$ChIP_seq{$target}{$ctrlchip}{type} eq $_[0]);

					next if ( (exists $ChIP_seq{$target}{$ctrlchip}{otherID} and 
							   exists $ChIP_seq{$target}{$chip}{otherID}) and
							   ($ChIP_seq{$target}{$ctrlchip}{otherID} 	ne $ChIP_seq{$target}{$chip}{otherID}) );
							   
					next if ($ChIP_seq{$target}{$ctrlchip}{experiment} 	ne $ChIP_seq{$target}{$chip}{experiment});
					next if ($ChIP_seq{$target}{$ctrlchip}{tissue} 		ne $ChIP_seq{$target}{$chip}{tissue});
					next if ($ChIP_seq{$target}{$ctrlchip}{factor} 		ne $ChIP_seq{$target}{$chip}{factor});
					next if ($ChIP_seq{$target}{$ctrlchip}{condition} 	ne $ChIP_seq{$target}{$chip}{condition});
					next if ($ChIP_seq{$target}{$ctrlchip}{treatment} 	ne $ChIP_seq{$target}{$chip}{treatment});
					next if ($ChIP_seq{$target}{$ctrlchip}{replicate} 	ne $ChIP_seq{$target}{$chip}{replicate});
					next if ($ChIP_seq{$target}{$ctrlchip}{pair_ended}  	ne $ChIP_seq{$target}{$chip}{pair_ended});
					
					my $extension;
					$extension = '.no_model_peaks.narrowPeak' if ($chipqc_pkcaller eq 'narrow');
					$extension = '.no_model_peaks.xls' 	  if ($chipqc_pkcaller eq 'macs');
					$extension = '.bed' 			  if ($chipqc_pkcaller eq 'bed');
					$extension = '.txt' 			  if ($chipqc_pkcaller eq 'raw');
					
					my $Peaks;
					$Peaks = "$ChIP_seq{$target}{$chip}{i_macs}" 
							 if ($ChIP_seq{$target}{$ctrlchip}{type} eq 'Input');
					$Peaks = "$ChIP_seq{$target}{$chip}{n_macs}" 
							 if ($ChIP_seq{$target}{$ctrlchip}{type} eq 'NegCtl');
					
					my @temp_peaks = split /\//, $Peaks;
					
					print $out "$ChIP_seq{$target}{$chip}{sampleID}\t";
					print $out "$ChIP_seq{$target}{$chip}{tissue}\t"; 
					print $out "$ChIP_seq{$target}{$chip}{factor}\t";
					print $out "$ChIP_seq{$target}{$chip}{condition}\t"; 
					print $out "$ChIP_seq{$target}{$chip}{treatment}\t";
					print $out "$ChIP_seq{$target}{$chip}{replicate}\t"; 
					print $out "$bwa_dir$ChIP_seq{$target}{$chip}{trim_file}\t" if ($bwa_trim eq 'yes');
					print $out "$bwa_dir$chip\t" 								if ($bwa_trim ne 'yes');
					print $out "$ChIP_seq{$target}{$ctrlchip}{sampleID}\t";
					print $out "$bwa_dir$ChIP_seq{$target}{$ctrlchip}{trim_file}\t" if ($bwa_trim eq 'yes');
					print $out "$bwa_dir$ctrlchip\t" if ($bwa_trim ne 'yes');
					print $out "$Peaks$extension\t";
					print $out "$chipqc_pkcaller\n";
					$number += 1;
				}
			}
		}
	}
	close $out;
	return ($factorfile, $number);
}

$now = `date "+%Y%m%d%H%M"`;
chomp $now;
my $chipqc_project = "ChIPQC_" . $now;
my @factors;
my %seen = ();

for my $target (keys %ChIP_seq)
{
	for my $fastq (keys %{ $ChIP_seq{$target} })
	{
		next unless (exists $ChIP_seq{$target}{$fastq}{factor});
		my $factor = $ChIP_seq{$target}{$fastq}{factor};
		unless (exists $seen{$factor})
		{
			push @factors, $ChIP_seq{$target}{$fastq}{factor};
			$seen{$ChIP_seq{$target}{$fastq}{factor}} = 1;
		}
	}
}

$logfile_dir = $chipqc_dir . "logfiles/";
system("mkdir -p $logfile_dir") != -1 or die $!;
my $logchipqcout = $logfile_dir . "ChIPQC_%J.out";
$number = 12 if ($number > 12);
$number = 4  if ($number < 4);

$pre_param = "--cpus-per-task=$bwa_threads --partition=$gcpart_macs ";
$pre_param.= "--mem=$mem_chipqc -o $logchipqcout -e $logchipqcout";

system("mkdir -p $chipqc_dir") != -1 or die $!; 
$chipqc_dir = $chipqc_dir . $today . '/';
system("mkdir -p $chipqc_dir") != -1 or die $!; 

unless ($chipqc_re_run eq 'yes')
{
	$newqueue = $myqueue;
	$pre_param.= " --dependency=afterok:$newqueue";
	
	$new_narrow_dir = $narrow_dir;
}
else
{
	$newqueue = 'empty';
	$new_narrow_dir = $macs_dir . 'narrow/' . $macs_date . '/';
}

my $chipqc_dir_popped = $chipqc_dir;
chop($chipqc_dir_popped);
$blacklist_file = 'null' unless (defined $blacklist_file);

# Generate ChIPQC report for individual factors
EXIT_IF: {
	if ($chipqc_re_run eq 'yes' or $slurm_array_id2 ne '')
	{
		last EXIT_IF if ($controls == 0);
		foreach my $factor (@factors)
		{
			generate_chipqc_table('Input', $factor, $factor);
			$pos_param = "$factorfile $chipqc_dir_popped $blacklist_file $genome_version $bwa_quality $chipqc_facetby $chipqc_colorby $bwa_threads";
	
			$slurm_job_id = submit_slurm_job("ChIPQC.sbatch", $pre_param, $pos_param, " ");
			print "\nSLURM ChIPQC on $factor (Input) is no. $slurm_job_id";
	
			generate_chipqc_table('NegCtl', $factor, $factor);
			$pos_param = "$factorfile $chipqc_dir_popped $blacklist_file $genome_version $bwa_quality $chipqc_facetby $chipqc_colorby $bwa_threads";
	
			$slurm_job_id = submit_slurm_job("ChIPQC.sbatch", $pre_param, $pos_param, " ");
			print "\nSLURM ChIPQC on $factor (NegCtl) is no. $slurm_job_id\n";
		}

		# Generate ChIPQC report for data set
		if (@factors > 0)    
		{
			generate_chipqc_table('Input', $chipqc_project, 'project');
			$pos_param = "$factorfile $chipqc_dir_popped $blacklist_file $genome_version $bwa_quality $chipqc_facetby $chipqc_colorby $bwa_threads";
		
			$slurm_job_id = submit_slurm_job("ChIPQC.sbatch", $pre_param, $pos_param, " ");
			print "\nSLURM ChIPQC on $chipqc_project (Input) is no. $slurm_job_id";

			generate_chipqc_table('NegCtl', $chipqc_project, 'project');
			$pos_param = "$factorfile $chipqc_dir_popped $blacklist_file $genome_version $bwa_quality $chipqc_facetby $chipqc_colorby $bwa_threads";

			$slurm_job_id = submit_slurm_job("ChIPQC.sbatch", $pre_param, $pos_param, " ");
			print "\nSLURM ChIPQC on $chipqc_project (NegCtl) is no. $slurm_job_id\n";
		}
		update_myqueue;
	}
}

# Motif Analysis
# Using Homer
if ($homer_run eq 'yes')
{
	$homer_dir = $processed_dir . "homer/" unless (defined $homer_dir);
	$logfile_dir = $homer_dir . "logfiles/";
	system("mkdir -p $logfile_dir") != -1 or die $!;

	my $loghomer = $logfile_dir . "Homer_%A_%a.out";
	my $homer_file_type;
	my $h_file_type;

	if ($homer_pk_caller eq 'MACS2') { $homer_file_type = '.no_model_peaks.narrowPeak'; }
	
	if ($homer_mask eq 'yes')	{ $homer_mask = '-mask'; }
	else				{ $homer_mask = 'void'; }
	
	if ($homer_bkgrnd eq 'yes')	{ $homer_bkgrnd = '-bg'; }
	else				{ $homer_bkgrnd = 'void'; }
	
	if ($homer_rna eq 'yes')	{ $homer_rna = '-rna'; }
	else				{ $homer_rna = 'void'; }
	
	if (defined $homer_known)	{ $homer_known = '-mknown ' . "$homer_known"; }
	else				{ $homer_known = 'void'; }
	
	if (defined $homer_norm)	{ $homer_norm = '-' . "$homer_norm"; }
	else				{ $homer_norm = '-gc'; }
	
	$pre_param = "--array=0-$no_macs --partition=$gcpart_bed ";
	$pre_param.= "--cpus-per-task=$bwa_threads --mem=$mem_bwa -o $loghomer -e $loghomer";
	
	if ($homer_re_run ne 'yes')
	{
		system("mkdir -p $homer_dir") != -1 or die $!; 
		$homer_dir = $homer_dir . $today;
		system("mkdir -p $homer_dir") != -1 or die $!; 
		
		$pre_param.= " --dependency=afterok:$myqueue";
		
		$new_narrow_dir = $narrow_dir;
	}
	else
	{
		$new_narrow_dir = $macs_dir . 'narrow/' . $macs_date . '/';
	}
	
	$h_file_type = 'null' unless (defined $h_file_type);
	$pos_param = "$genome $homer_dir $homer_file_type $homer_size_up $homer_size_down ";
	$pos_param.= "$homer_bkgrnd $homer_mask $homer_rna $homer_known $homer_norm ";
	$pos_param.= "$bwa_threads $new_narrow_dir $homer_scoring $h_file_type ";

	if ($homer_re_run eq 'yes' or $myqueue ne '')
	{
		$slurm_array_id = submit_slurm_array("Homer.sbatch", $pre_param, $pos_param, @macs_list);
		update_myqueue;
	
		$pre_param = "--dependency=afterok:$slurm_array_id --array=0-$no_macs --partition=$gcpart_bed ";
		$pre_param.= "--cpus-per-task=$low_threads --mem=$mem_bwa -o $loghomer -e $loghomer";
		$pos_param = "$genome $homer_dir $homer_file_type $homer_size_up $homer_size_down ";
		$pos_param.= "$homer_bin $new_narrow_dir";

		$slurm_array_id = submit_slurm_array("Homer_annot.sbatch", $pre_param, $pos_param, @macs_list);
		update_myqueue;
	}
	else { print "No Homer motif analysis selected\n"; }
}
else { print "No Homer motif analysis selected\n"; }

## Differential peak calling
# Using MACS2 bdgdiff
if ($diff_run eq 'yes')
{
	$diff_dir = $processed_dir . "diffcall/" unless (defined $diff_dir);
	system("mkdir -p $diff_dir") != -1 or die $!; 
	$diff_dir = $diff_dir . $today . '/';
	system("mkdir -p $diff_dir") != -1 or die $!; 
	
	my $logfile_dir = $diff_dir . "logfiles/";
	system("mkdir -p $logfile_dir") != -1 or die $!;

	my $logdiff = $logfile_dir . "bdgdiff_%J.out";
	$pre_param = "--partition=$gcpart_bwa --mem=$mem_bwa -o $logdiff -e $logdiff";

	unless ($diff_re_run eq 'yes')
	{
		$pre_param.= " --dependency=afterok:$myqueue";
		$new_broad_dir = $broad_dir;
		$new_narrow_dir = $narrow_dir;
	}
	else
	{
		$new_broad_dir = $macs_dir . 'broad/' . $macs_date . '/';
		$new_narrow_dir = $macs_dir . 'narrow/' . $macs_date . '/';
	}

	$pos_param = "$new_broad_dir $macs_extsize $diff_dir $mem_bwa $diff_cutoff $diff_minlen";
	$pos_param.= " $diff_maxgap $diff_depth1 $diff_depth2 $logfile_dir $controls $new_narrow_dir";

	$slurm_array_id = submit_slurm_job("MACS_diff.sbatch", $pre_param, $pos_param, @macs_list);
	update_myqueue;
}
