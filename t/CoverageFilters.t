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


my $coverage_filter_input_file  = File::Spec->join($test_data_dir, 'Test_filtered.readcounts.xls');
my $coverage_filter_output_file = File::Spec->join($output_dir, 'Test_filtered.readcounts.covfilt.xls');
my $coverage_filter_executable  = File::Spec->join($executable_dir, 'CoverageFilters.pl');
my $coverage_filter_command = join(' ',
    "perl $coverage_filter_executable",
    "--input-file $coverage_filter_input_file",
    "--output-file $coverage_filter_output_file",
    "--norm --tdna --trna --expn",
);
ok(!system($coverage_filter_command), 'CoverageFilters command executes successfully');
ok(
    !compare($coverage_filter_output_file, File::Spec->join($test_data_dir, 'Test_filtered.readcounts.covfilt.xls')),
    'CoverageFilters output as expected',
);

done_testing;
