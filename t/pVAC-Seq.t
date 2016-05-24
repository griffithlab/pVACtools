#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use Cwd qw(abs_path);
use File::Basename qw(dirname);
use File::Spec;

my $base_dir        = File::Spec->catdir(dirname(abs_path($0)), File::Spec->updir());
my $pvac_seq_script = File::Spec->join($base_dir, 'pVAC-Seq.pl');

my $compile_output = `$^X -c $pvac_seq_script 2>&1`;
like( $compile_output, qr/syntax OK$/, 'pVAC-Seq script compiles' );

done_testing;
