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
my $epitope_length          = 9;
my $allele                  = 'HLA-A29:02';

my $parse_output_netmhc_executable  = File::Spec->join($executable_dir, 'ParseOutputNetmhc.pl');
my $parse_output_netmhc_compile_output = `$^X -c $parse_output_netmhc_executable 2>&1`;
like( $parse_output_netmhc_compile_output, qr/syntax OK$/, 'ParseOutputNetmhc script compiles' );

my $parse_output_netmhc_input_file  = File::Spec->join($test_data_dir, "${sample_name}.${allele}.${epitope_length}.netmhc.xls");
my $parse_output_netmhc_key_file    = File::Spec->join($test_data_dir, "${sample_name}_${peptide_sequence_length}.key");
my $parse_output_netmhc_output_file = File::Spec->join($output_dir, "${sample_name}.${allele}.${epitope_length}.netmhc.parsed");
my $parse_output_netmhc_command = join(' ',
    "perl $parse_output_netmhc_executable",
    "-i $parse_output_netmhc_input_file",
    "-k $parse_output_netmhc_key_file",
    "-o $parse_output_netmhc_output_file",
);
ok(!system($parse_output_netmhc_command), 'ParseOutputNetmhc command executes successfully');
ok(
    !compare($parse_output_netmhc_output_file, File::Spec->join($test_data_dir, "${sample_name}.${allele}.${epitope_length}.netmhc.parsed")),
    'ParseOutputNetmhc output as expected',
);

done_testing;
