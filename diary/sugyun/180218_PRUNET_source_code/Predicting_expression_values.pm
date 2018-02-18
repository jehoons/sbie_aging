#!/usr/bin/perl -
use strict;
use warnings;
use boolean_functions_generator;
use Translator;

package Predicting_expression_values;

sub pev{
	my ($updating_scheme,$updating_sequence_array,$raw_network_interactions_array,$configurations_selection_array,$phenotype_array) = @_;
	my @raw_network_interactions = @{$raw_network_interactions_array};
	my @configurations_selection = @{$configurations_selection_array};
	my @updating_sequence = @{$updating_sequence_array};
	my @phenotype = @{$phenotype_array};
	
	my %expression_data_trough_configurations=();
	my %gene_count_through_configurations=();
	#Filling the hashes with all the genes in the raw network
	foreach my $interaction(@raw_network_interactions){
				if ($interaction =~ /(.+)(\s)(\-[\>\|])(\s)(.+)/){
					$expression_data_trough_configurations{$1}='';
					$gene_count_through_configurations{$1}=0;
					$expression_data_trough_configurations{$5}='';
					$gene_count_through_configurations{$5}=0;
				}else{}
	}
	my %initial_state_phenotype;
	for (my $i=0;$i<@phenotype;$i++){
			if ($phenotype[$i] =~ /(.+)(\s)([01])/){
				$initial_state_phenotype{$1} = $3;
			}
			else{}
	}#Filling the hash with missing genes in the configuration: if the gene expression value is not known it becomes '0'
	foreach my $gene(keys %expression_data_trough_configurations){
		if (exists $initial_state_phenotype{$gene}){}
		else{
			$initial_state_phenotype{$gene}=0;
		}
	}
	
	my @genes_raw_network_ordered = sort {$a cmp $b} keys %initial_state_phenotype;	
	my ($states_phenotype_array,$type_of_attractor_phenotype);
	
#Main for###############################################################################################################################
	for (my $i=0;$i<@configurations_selection;$i++){
		my @network_lines = keys %{$configurations_selection[$i]};
			my %genes=();#Genes in the configuration
			foreach my $line(@network_lines){
				if ($line =~ /(.+)(\s)(\-[\>\|])(\s)(.+)/){
					$genes{$1}=$1;
					$genes{$5}=$5;
				}else{}
			}
			foreach my $configuration_gene(keys %genes){
					$gene_count_through_configurations{$configuration_gene}=$gene_count_through_configurations{$configuration_gene}+1;
			}
			
			my @genes_ordered = sort {$a cmp $b} keys %genes;
			my $initial_state_phenotype_string = '';
			foreach my $gene_ordered(@genes_ordered){
				$initial_state_phenotype_string = $initial_state_phenotype_string.$initial_state_phenotype{$gene_ordered};
			}
			my @boolean_functions_view_content = @{Boolean_functions_generator::boolean_functions(\@network_lines)};
			my $number_of_genes = scalar @boolean_functions_view_content;
			my %genes_encoded = ();
			my @boolean_functions_formated = ();
			
			for (my $i=0;$i<@boolean_functions_view_content;$i++){
				if($boolean_functions_view_content[$i] =~ /(.+)(\t)(\=)(\s)(.+)/){
					my $index_correction=$i+1;
					$genes_encoded{$1} = "x$index_correction";
				}else{}	
			}
			
			for (my $i=0;$i<@boolean_functions_view_content;$i++){
				if($boolean_functions_view_content[$i] =~ /(.+)(\t)(\=)(\s)(.+)/){
					my $index_correction=$i+1;
					my $raw_string=$5;
					for (my $j=0;$j<@genes_ordered;$j++){
						if ($raw_string =~ /(.*)([\!\(])($genes_ordered[$j])([\)\&\|])(.*)/){
							$raw_string = $1.$2.$genes_encoded{$3}.$4.$5;
						}
					}
					####################  Newly added ###########################
					## Check again for substitution 
					my @ar=split"\&",$raw_string;
					for(my $i=0; $i<scalar @ar; $i++){
					    if ($ar[$i] =~ /^(\(\(\!)(.+?)(\).*)$/){  
						if($2 !~ /^x\d+$/) {
						    $ar[$i] = $1.$genes_encoded{$2}.$3;
						}
					    }
					    elsif($ar[$i] =~ /^(\(\()(.+?)(\).*)$/){
						if($2 !~ /^x\d+$/) {
						    $ar[$i] = $1.$genes_encoded{$2}.$3;
						}
					    }
					}
					$raw_string = join("\&",@ar);
					print $raw_string,"\n";
					#############################################################
					$boolean_functions_formated[$i]= "f$index_correction \= $raw_string";
					$boolean_functions_formated[$i]=~ s/\!/\~/g;
					$boolean_functions_formated[$i]=~ s/\&/\*/g;
					$boolean_functions_formated[$i]=~ s/\|/\+/g;
				}else{}
			}
			
			my @polynomial_functions = @{Translator::translator(\@boolean_functions_formated,$number_of_genes)};
			my @polynomial_functions_to_eval = ();
			
			foreach my $polynomial_function(@polynomial_functions){
				$polynomial_function =~ s/x(\d+)/\$x\[$1\]/g;
				$polynomial_function =~ /(.+)(\s\=\s)(.+)(\n)/;
				push @polynomial_functions_to_eval, $3;
			}
			if ($updating_scheme eq 'synchronous'){
				($states_phenotype_array,$type_of_attractor_phenotype) = trajectory_synchronous(\@polynomial_functions_to_eval,$initial_state_phenotype_string,$number_of_genes);
			}elsif($updating_scheme eq 'asynchronous'){
				my @updating_sequence_new = ();
				foreach my $gene_ordered(@updating_sequence){
					if (exists $genes{$gene_ordered}){
						push @updating_sequence_new, $gene_ordered;
					}else{}
				}
				my @updating_sequence_encoded_new=@{encoding_updating_sequence(\@updating_sequence_new)};
				($states_phenotype_array,$type_of_attractor_phenotype) = trajectory_asynchronous(\@polynomial_functions_to_eval,$initial_state_phenotype_string,$number_of_genes,\@updating_sequence_encoded_new);
			}else{}
			if ($type_of_attractor_phenotype eq 'fixed point'){
			
				my @temp_array = @{$states_phenotype_array};
				my @states_phenotype = split //, @temp_array;
				my @predicted_phenotype = split //, (pop @temp_array);
				for(my $i=0;$i<@genes_ordered;$i++){
					$expression_data_trough_configurations{$genes_ordered[$i]}=$expression_data_trough_configurations{$genes_ordered[$i]}.$predicted_phenotype[$i];
				}
			
			}
			else{}
	}
#Main for end###############################################################################################################################
	my %missing_nodes_count = ();
	my %expression_predictions = ();
	foreach my $key (keys %expression_data_trough_configurations){
		$missing_nodes_count{$key}=length($expression_data_trough_configurations{$key});
		if ($missing_nodes_count{$key} eq 0){$expression_predictions{$key}='missing'}
		else{
			$expression_predictions{$key}= (eval(join '+',((split //,$expression_data_trough_configurations{$key}))))/($missing_nodes_count{$key});
		}	
	}
	foreach my $key(keys %expression_predictions){
		if ($expression_predictions{$key} eq 'missing'){}
		else{
			$expression_predictions{$key}= sprintf("%.2f",$expression_predictions{$key});
		}	
	}
	return(\%expression_predictions);
}
sub trajectory_synchronous{
		my ($functions_array,$initial_state_string,$number_of_nodes) =@_;
		my $number_of_possible_values=2;
		my $longest_trajectory = 1000;#####################################################################################################################I have to limit this somehow
		my @functions = @{$functions_array};
		my @y = split //, $initial_state_string;
		my @x = (0,@y);
		my @y_index_correction = (0);
		my %states_visited=();
		my %perturbations = ();
		my @states;
		foreach (0..$longest_trajectory){
				my @y_index_correction = (0);
				my $i = 1;
				
				foreach my $f(@functions){	
					if (exists $perturbations{$i-1}){
						$y[$i] = $perturbations{$i-1};
						$y_index_correction[$i] = $y[$i];
						$i++;
					}else{
						$y[$i] = eval($f) % $number_of_possible_values;
						push @y_index_correction, $y[$i];
						$i++;
					}
				}
				shift @y_index_correction;
				my $new_state = join "", @y_index_correction;
				if (exists $states_visited{$new_state}){
					my $last_visited_state = $states[(scalar @states)-1];
					if ($last_visited_state eq $new_state){
						push @states, $new_state;
						return(\@states,'fixed point');
					}else{
						push @states, $new_state;
						return(\@states,'limit cycle');
					}
				}else{
					$states_visited{$new_state}=$new_state;
				}
				push @states, $new_state;
				#Updating x
				foreach my $index(1..$number_of_nodes){
					$x[$index] = $y[$index];
				}
			}
}

sub trajectory_asynchronous{
		my ($functions_array,$initial_state_string,$number_of_nodes,$updating_sequence_array) =@_;
		my @updating_order_indexes=@{$updating_sequence_array};
		my $number_of_possible_values=2;
		my $longest_trajectory = 1000;
		my @functions = @{$functions_array};
		my @y = split //, $initial_state_string;
		my @x = (0,@y);
		my @y_index_correction = (0);
		my %states_visited=();
		my %perturbations = ();
		my @states;
		foreach (0..$longest_trajectory){
				my @y_index_correction = (0);
				for (my $j=0;$j<@updating_order_indexes;$j++){
					my $a = $updating_order_indexes[$j];
					
					if (exists $perturbations{$a}){
						$y[$a+1] = $perturbations{$a};
						$x[$a+1] = $y[$a+1];
						$y_index_correction[$a+1] = $y[$a+1];
					}else{
						$y[$a+1] = eval($functions[$a]) % $number_of_possible_values;
						$x[$a+1] = $y[$a+1];
						$y_index_correction[$a+1] = $y[$a+1];
					}
				}
				shift (@y_index_correction);
				my $new_state = join "", @y_index_correction;
				if (exists $states_visited{$new_state}){
					my $last_visited_state = $states[(scalar @states)-1];
					if ($last_visited_state eq $new_state){
						push @states, $new_state;
						return(\@states,'fixed point');
					}else{
						push @states, $new_state;
						return(\@states,'limit cycle');
					}
				}else{
					$states_visited{$new_state}=$new_state;
				}
				push @states, $new_state;
				#Updating x
				foreach my $index(1..$number_of_nodes){
					$x[$index] = $y[$index];
				}
			}
}

##Subroutine to asignate a number to each gene of the list of genes that are first actualized 
sub encoding_updating_sequence{
	my ($updating_sequence_view_content_array)=@_;
	my @updating_sequence_view_content = @{$updating_sequence_view_content_array};
	my %aa = ();
	my @updating_sequence_encoded = ();
	my @a = sort {$a cmp $b} @updating_sequence_view_content;
	for (my $i=0;$i<@a;$i++){
			$aa{$a[$i]}=$i;
	}
	for (my $j=0;$j<@updating_sequence_view_content;$j++){
		$updating_sequence_encoded[$j]=$aa{$updating_sequence_view_content[$j]};
	}
	return (\@updating_sequence_encoded);
}
1























