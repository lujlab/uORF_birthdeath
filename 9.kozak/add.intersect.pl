#!/usr/bin/perl
#

my $file1 = $ARGV[0];
my $file2 = $ARGV[1];
my $col1 = $ARGV[2]-1;
my $col2 = $ARGV[3]-1;
my $col3 = $ARGV[4]-1;
my $info = $ARGV[5];

open FH1 ,"<$file1";
open FH2 ,"<$file2";

while(<FH1>){
	chomp;
	my @inlines = split(/\t/,$_);
	my $key_A = $inlines[$col1];
	$Appear{$key_A} = 1;
	$hash{$key_A} = $inlines[$col2];
}

close(FH1);

while(<FH2>){
	chomp;
	my @inlines = split(/\t/);
	my $key_B = $inlines[$col3];
	if($Appear{$key_B} == 1)
	{ print $_,"\t",$hash{$key_B},"\n"; }
	else{ print $_,"\t",$info,"\n"; }
}

close(FH2);


