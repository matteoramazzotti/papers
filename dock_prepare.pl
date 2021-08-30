#!/usr/bin/perl
#used to prepare protein structure for docking
$clean = 1;
$fname = $ARGV[0];

$oname = $ARGV[0];
$oname =~ s/\.pdb.+/.multi.pdb/;

$name = $ARGV[0];
$name =~ s/.+\///;
$name =~ s/\.pdb.+//;

#goto JUMP;
open(P,">script.pml");
print P "load $fname\n";
print P "symexp sym,$name,all,1\n";
print P "save $oname\n";
close P;

`pymol -c script.pml`;
unlink 'script.pml';

JUMP:
open(IN,$oname) or die "Cannot open $oname";
@chains = qw /A B C D E F G H I J K L M N O P Q R S T U V W X Y Z/;
$old = 'XXX';
$cnt = -1;

while($line = <IN>) {
	if ($clean) {
		next if ($line =~ /^HETATM/);
		next if ($line =~ /^ANISOU/);
		next if ($line =~ /^TER/);
		next if ($line =~ /^CONECT/);
	}
	if ($line =~ /^ATOM/) {
		$k = substr($line,21,1);
		$cnt++ if ($k ne $old);
		$line =~ s/^(.....................)(\w)/$1$chains[$cnt]/;
		$old = $k;
	}
	print $line;
}
close IN;

