#!/usr/bin/perl -
use strict;
use warnings;

package Scoring_edges;

sub se{
	my ($raw_network_interactions_array,$configurations_selection_array) = @_;
	my @raw_network_interactions = @{$raw_network_interactions_array};
	my @configurations_selection = @{$configurations_selection_array};
	my %edges_frequences=();
	foreach my $raw_network_interaction(@raw_network_interactions){
		$edges_frequences{$raw_network_interaction}=0;
	}
	for (my $i=0;$i<@configurations_selection;$i++){
		my %configuration=%{$configurations_selection[$i]};
		#my %configuration = ();
		for (my $j=0;$j<scalar @raw_network_interactions;$j++){
			if (exists $configuration{$raw_network_interactions[$j]}){
				$edges_frequences{$raw_network_interactions[$j]}=($edges_frequences{$raw_network_interactions[$j]})+(1/scalar @configurations_selection);
			}
			else{}
		}
	}
	foreach my $key(keys %edges_frequences){
		if ($edges_frequences{$key} == '1'){$edges_frequences{$key}='1.00'}
		elsif ($edges_frequences{$key} == '0'){$edges_frequences{$key}='0.00'}
		else{$edges_frequences{$key}=sprintf("%2.2f",$edges_frequences{$key})}
	}
	return(\%edges_frequences);
}
1