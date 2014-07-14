#!/usr/bin/perl -w
use strict;

open(F, $ARGV[0])||die;
open(R1, ">$ARGV[0]_R1.fq")||die;
open(R2, ">$ARGV[0]_R2.fq")||die;

while(my $r = <F>)
{
        print R1 $r;
        $r = <F>;
        print R1 $r;
        $r = <F>;
        print R1 $r;
        $r = <F>;
        print R1 $r;

        $r = <F>;
        print R2 $r;
        $r = <F>;
        print R2 $r;
        $r = <F>;
        print R2 $r;
        $r = <F>;
        print R2 $r;

}


