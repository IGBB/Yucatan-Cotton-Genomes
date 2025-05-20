use strict;
use warnings;
use autodie;

local $, = "\t";
local $\ = "\n";

my $chrs = {};
while (<>){

    my ($qname, $qlen, undef, undef, $strand,
        $tname, undef, undef, undef, undef, $len, undef) = split("\t", $_, 12);

    next unless ($len >= 250000);

    $chrs->{$qname} ||= {len=>$qlen, matches=>{}};
    $chrs->{$qname}->{matches}->{$tname} ||= { name=>$tname,
                                               total=>0,
                                               negstrand=>0};
    $chrs->{$qname}->{matches}->{$tname}->{total} += $len;

    if($strand eq '-'){
        $chrs->{$qname}->{matches}->{$tname}->{negstrand} += $len;
    }
}

my @allnames = sort {$chrs->{$b}->{len} <=> $chrs->{$a}->{len}} keys(%$chrs);

foreach my $name (@allnames){
    my $len=$chrs->{$name}->{len};

    my @sorted_match=sort {$b->{total} <=> $a->{total}}
                      values(%{$chrs->{$name}->{matches}});
    my $match = shift(@sorted_match);

    my $strand=(($match->{negstrand}/$match->{total}) < 0.5)?'+':'-';
    print $match->{name} , 1, $len, 1, 'W', $name, 1, $len, $strand;
}
