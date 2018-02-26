package Boolean_functions_generator;

use strict;
use warnings;

sub boolean_functions{
	my ($array_lines) = @_;
	my @lines = @{$array_lines};
	my @boolean_functions;
	my %genes=();
	my %genes_regulated=();
	foreach my $line(@lines){
	if ($line =~ /(.+)(\s)(\-[\>\|])(\s)(.+)/){
		$genes{$1}=$1;
		$genes{$5}=$5;
	}
	else{}
}
	my @genes_ordered = sort {$a cmp $b} keys %genes;
	my %genes_translated=();
	#for (my $j=0;$j <@genes_ordered;$j++){
	#	$genes_translated{$genes_ordered[$j]}="x".($j+1);
	#}
	foreach my $line(@lines){
		if ($line =~ /(.+)(\s)(\-[\>\|])(\s)(.+)/){
			
			if (exists $genes_regulated{$5}){$genes_regulated{$5}=$genes_regulated{$5}."\/".$1."\:".$3}
			else{$genes_regulated{$5}=$1."\:".$3}
		}
		else{}
	}
	
	for (my $i=0;$i<@genes_ordered;$i++){
		if (exists $genes_regulated{$genes_ordered[$i]}){
			my @regulators = split /\//, $genes_regulated{$genes_ordered[$i]};
			my @activators=();
			my @inhibitors=();
			foreach my $regulator(@regulators){
				$regulator =~ /(.+)(\:)(\-[\>\|])/;
				if ($3 eq '->'){push @activators, $1;}
				elsif($3 eq '-|'){push @inhibitors, $1;}
				else{}
			}
			########
			my $function_string='';
			for (my $j=0;$j<@inhibitors;$j++){
				if ($j==0){
					$function_string = "\(\!$inhibitors[$j]\)";
				}
				else{
					$function_string = $function_string."\&\(\(\!$inhibitors[$j]\)";
				}
			}
			for (my $j=0;$j<@activators;$j++){
				if (($j==0) and (length $function_string >0)){	
					$function_string = $function_string."\&\($activators[$j]";
				}
				elsif (($j==0) and (length $function_string ==0)){	
					$function_string = "$activators[$j]";
				}
				else {
					$function_string = $function_string."\|\($activators[$j]";
				}

			}
			my $number_ending_parentesis= @inhibitors+@activators-1;
			for (my $j=0;$j<$number_ending_parentesis;$j++){
				$function_string = $function_string."\)";
			}
			$function_string = "\($function_string\)";
			push @boolean_functions, $function_string;
		}	
		else{
			push @boolean_functions,"0";
		}
	}
	my @final_boolean_functions = ();
	for (my $i=0;$i<@boolean_functions;$i++){
		my $index_correction = $i+1;
		push @final_boolean_functions, "$genes_ordered[$i]\t= $boolean_functions[$i]\n"
	}
	
	return(\@final_boolean_functions);
}
1;	