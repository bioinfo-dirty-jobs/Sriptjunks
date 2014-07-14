#!/usr/bin/perl -w
use strict;

open(F, "$ARGV[0]")||die;

open(O1, ">$ARGV[0].1.fq")||die;
open(O2, ">$ARGV[0].2.fq")||die;


my $name1 = "";
my $name2 = "";
my $content1 = "";
my $content2 = "";

while(my $r = <F>)
{
	chomp $r;
	my @v = split(/\t/, $r);

	my ($name, $read) = $v[0] =~ /^(.+)\/([12])/;
	
	if($read == 1)
	{
		$content1 = $r;
		$name1 = $name;
	}
	if($read == 2)
	{
		$content2 = $r;
		$name2 = $name;
	}

	if($read == 2)
	{
		if($name2 eq $name1)
		{
			my @v1 = split(/\t/, $content1);
                        my @v2 = split(/\t/, $content2);

			print O1 "\@$v1[0]\n$v1[1]\n+$v1[0]\n$v1[2]\n";
                        print O2 "\@$v2[0]\n$v2[1]\n+$v2[0]\n$v2[2]\n";
		}
	}
}
