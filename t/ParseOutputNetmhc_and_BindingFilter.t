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

my $parse_output_netmhc_input_file  = File::Spec->join($test_data_dir, "${sample_name}.${allele}.${epitope_length}.netmhc.xls");
my $parse_output_netmhc_key_file    = File::Spec->join($test_data_dir, "${sample_name}_${peptide_sequence_length}.key");
my $parse_output_netmhc_output_file = File::Spec->join($output_dir, "${sample_name}.${allele}.${epitope_length}.netmhc.parsed");
my $parse_output_netmhc_executable  = File::Spec->join($executable_dir, 'ParseOutputNetmhc.pl');
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

my $binding_filter_input_file  = File::Spec->join($test_data_dir, 'annotated_variants.tsv');
my $binding_filter_output_fof  = File::Spec->join($output_dir, 'Test.fof');
open(my $fh, '>', $binding_filter_output_fof);
print $fh $parse_output_netmhc_output_file;
close $fh;
my $binding_filter_output_file = File::Spec->join($output_dir, "${sample_name}_filtered.xls");
my $binding_filter_executable  = File::Spec->join($executable_dir, 'BindingFilter.pl');
my $binding_filter_command = join(' ',
    "perl $binding_filter_executable",
    "-i $binding_filter_input_file",
    "-f $binding_filter_output_fof",
    "-c 0",
    "-b 500",
    "-o $binding_filter_output_file",
);
ok(!system($binding_filter_command), 'BindingFilter command executes successfully');
ok(
    !compare($binding_filter_output_file, File::Spec->join($test_data_dir, "${sample_name}_filtered.xls")),
    'BindingFilter output as expected',
);

done_testing;
