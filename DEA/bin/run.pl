use warnings;
use strict;



#=============================parse arguments=============================
my ($in_expression,$in_sampleInfor,$in_ppi,$in_disease,$out_path);
if(scalar@ARGV != 10)
{
    &print_help();
    exit();
}else{
    foreach my $i(0..$#ARGV)
    {
        if ($ARGV[$i] eq '-expr')
        {
            $in_expression = $ARGV[$i+1];
        }elsif($ARGV[$i] eq '-sample')
        {
            $in_sampleInfor = $ARGV[$i+1];
        }elsif($ARGV[$i] eq '-ppi')
        {
            $in_ppi = $ARGV[$i+1];
        }elsif($ARGV[$i] eq '-disease')
        {
            $in_disease = $ARGV[$i+1];
        }elsif($ARGV[$i] eq '-out')
        {
            $out_path = $ARGV[$i+1];
        }
    }
}

#print $out_path;
mkdir $out_path;
#end

#=================generate the Non-DI pair===================================
system "Rscript DEA.R $in_expression $in_sampleInfor $in_disease $out_path";

my %deg_gene;
open(IN,"<",$out_path."deg.txt") or die "not open\n";
while (my $line=<IN>) {
	chomp $line;
	my @line = split(/\t/,$line);
	$deg_gene{$line[0]}=1;
}
close IN;

open(INN,"<",$in_ppi) or die "not open\n";
open(OUT,">",$out_path."Non_DI.txt") or die "not open\n"; 
while (my $line=<INN>) {
	chomp $line;
	my @line = split(/\s+/,$line);
	if (exists $deg_gene{$line[0]} && exists $deg_gene{$line[1]}) {
		print OUT "$line[0]\t$line[1]\n";
	}	
}
close INN;
close OUT;

sub print_help {
    printf "DEA USAGE:
    perl run.pl -expr input_expression_file -sample input_sampleInformation_file -ppi input_ppi_file -disease input_disease_state -out /the/directory/of/output/
    such as:
    perl run.pl -expr ../examples/input/GSE29429_expression.txt -sample ../examples/input/GSE29429_sampleInfor.txt -ppi ../examples/input/PPI.txt -disease HMI -out ../examples/output/
    ==================
    arguments:
    ==================
    -expr       the path of expression matrix, such as ../examples/expression.txt
    -sample     the path of sampleInformation
    -ppi        the path of protein-protein interations
	-disease    the disease state(HMI/MMI/MHCI)
    -out        the path of results\n";
}
