#!/usr/bin/env perl

#
# The BVBRC App wrapper for Genomad with integrated GenomeAssembly2
# and plasmid and viral (LowVan) annotation.

# https://github.com/apcamargo/genomad
#

use strict;
use warnings;
use Carp::Always;
use Bio::KBase::AppService::AppScript;

# Import the functions we need from our module
use MobileElementDetection qw(run preflight);

# Create the AppScript with references to our imported subroutines
my $script = Bio::KBase::AppService::AppScript->new(\&run, \&preflight);

# Run the script with command line arguments
$script->run(\@ARGV);
