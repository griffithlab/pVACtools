#!/usr/bin/env perl

use strict;
use warnings;

use Test::More;
use File::Temp qw(tempdir);
use Cwd qw(abs_path);
use File::Basename qw(dirname);
use File::Compare qw (compare);

my $base_dir       = File::Spec->catdir(dirname(abs_path($0)), File::Spec->updir());
my $executable_dir = File::Spec->catdir($base_dir, 'bin');
my $test_data_dir  = File::Spec->catdir($base_dir, 'test_data');
my $output_dir     = tempdir(CLEANUP => 1);

my $sample_name             = 'Test';
my $peptide_sequence_length = 21;

my $generate_fasta_key_input_file  = File::Spec->join($test_data_dir, "${sample_name}_${peptide_sequence_length}.fa");
my $generate_fasta_key_output_file = File::Spec->join($output_dir, "${sample_name}_${peptide_sequence_length}.key");
my $generate_fasta_key_executable  = File::Spec->join($executable_dir, 'GenerateFastaKey.pl');
my $generate_fasta_key_command = join(' ',
    "perl $generate_fasta_key_executable",
    "-i $generate_fasta_key_input_file",
    "-o $generate_fasta_key_output_file",
);
ok(!system($generate_fasta_key_command), 'GenerateFastaKey command executes successfully');
ok(
    !compare($generate_fasta_key_output_file, File::Spec->join($test_data_dir, "${sample_name}_${peptide_sequence_length}.key")),
    'GenerateFastaKey output as expected',
);

done_testing;
