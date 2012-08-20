#!/usr/bin/perl 
use strict;
use Getopt::Std;
use Cwd;

my $home_dir = $ENV{"HOME"};

# Read Program Options i
my %opt = ();
getopt("sfrp", \%opt);
die "usage example:writeDefaultAffy.pl -s GSE1234 -f GSE1234_workflow.Rnw -r GSE1234 -p mouse430a2\n".
    "parameter s is the study ID\n".
    "parameter f is the output file name\n".
    "parameter r is the raw data directory\n".
    "parameter p is the platform "
     unless ($opt{s} && $opt{f} && $opt{r} && $opt{p}); 
my $study = $opt{s};
my $file = $opt{f};
my $rawdata = $opt{r};
my $platform = $opt{p};

open O, ">".$file;

print O  '<<load.data>>='."\n";
print O  "library(metaGEO)\n";
print O  'cel.files <- list.celfiles("'.$rawdata.'/")'."\n";
print O  $study.'.abatch <- ReadAffy(filenames=cel.files)'."\n";
print O  '@'."\n";

print O  '<<model.data>>='."\n";
print O  'pms<- pm('.$study.'.abatch)'."\n";
print O  'pms[pms <= 1] <- 1'."\n";
print O  'data <- log2(pms)'."\n";
print O  'bio.var <- NULL'."\n";
print O  'adj.var <- NULL'."\n";
print O  'int.var <- data.frame(array = factor(1:ncol(data)))'."\n";
print O   $study.'.snm.fit <- snm(data, bio.var,adj.var,int.var,diagnose=FALSE)'."\n"; 
print O  'colnames('.$study.'.snm.fit$norm.dat) <- sampleNames('.$study.'.abatch)'."\n";

print O  'jpeg(file="'.$study.'.snm.jpg")'."\n";
print O  'plot('.$study.'.snm.fit)'."\n";
print O  'dev.off()'."\n";

print O  '@'."\n";
print O  '<<summarize.to.genes>>='."\n";

print O  'data('.$platform.')'."\n";
print O  $study.'.fits <- fit.pset('.$study.'.snm.fit$norm.dat, '.$platform.'.gene2row, "'.$study.'")'."\n";
print O  'save('.$study.'.fits,file="'.$study.'.fits.Rda")'."\n";
print O  '@'."\n";
