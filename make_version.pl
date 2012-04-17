#!/usr/bin/perl
my $out = shift;
my $ver = <>;
chomp $ver;
my $exists = `grep $ver $out`;
if ($exists !~ /$ver/) {
    open(my $fh, ">$out") || die;
    print $fh "#define VERSION \"$ver\"\n";
}
