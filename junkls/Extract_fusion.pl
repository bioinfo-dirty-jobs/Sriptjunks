#!/usr/bin/perl -w
use strict;
my $name = $ARGV[0];
my $out = $name;

qx{samtools view  -f 72 $name|cut -f1,10,11  > $out.first.txt};
qx{samtools view  -f 136 $name|cut -f1,10,11 > $out.second.txt};
qx{samtools view  -f 64 unmapped.bam|cut -f1,10,11 >> $out.first.txt};
qx{samtools view  -f 128 unmapped.bam|cut -f1,10,11 >> $out.second.txt};
qx{sort $out.first.txt|uniq >tmp.1};
qx{sort $out.second.txt|uniq >tmp.2};
open(T,"tmp.1");
open(Q,">unito.tmp");
while(my $r = <T>)
{
	chomp $r;
	my @v = split(/\t/,$r);
	print Q "$v[0]\/1\t$v[1]\t$v[2]\n";


}
close Q;
close T;
open(G,">>unito.tmp");
open(TR,"tmp.2");
while(my $r = <TR>)
{
        chomp $r;
        my @v = split(/\t/,$r);
        print G "$v[0]\/2\t$v[1]\t$v[2]\n";


}


close G;
close TR;


qx{cat unito.tmp|sort  > uniti_sort.txt};







