#!/usr/bin/perl -w
#### Ana Rodriguez Sanchez-Archidona, Luxembourg 11/7/2013, Luxembourg Center for Systems Biomedicine ####
#### Optimization algorithm module: EDA performing iterative prune of a literature network using populations of alternative pruned networks
#### that are scored and selected using expression data

use strict;
use warnings;

use boolean_functions_generator;
use Translator;

package Xpred_multivariate;
my @population_of_networks_to_return = ();
my @scores_of_population_of_networks_to_return = ();

my @population_of_networks_to_return_network_1 = ();
my @scores_of_population_of_networks_to_return_network_1 = ();

my @population_of_networks_to_return_network_2 = ();
my @scores_of_population_of_networks_to_return_network_2 = ();

my @results = ();
my @elite = ();
my @elite_1 = ();
my @elite_2 = ();
my %positive_circuits =();

###Optimization algorithm subroutine
sub xpred{
	open FILEHANDLE, ">./optimization_report.txt";
	close FILEHANDLE;
	my ($network_view_content_array,                    #Literature raw network 
		$phenotype_1_content_array,                     #Genes expression data: experiment 1
		$phenotype_2_content_array,                     #Genes expression data: experiment 2
		$fixed_interactions_1_content_array,            #Interactions to be keeped from the 1 steady state
		$fixed_interactions_2_content_array,            #Interactions to be keeped from the 2 steady state
		$updating_sequence_array,                       #List of genes that are actualized first
		$maximum_iterations,                            #Number of iterations in the contextualization process
		$population,                                    #Number of subnetworks in the population
		$selection,                                     #Number of subnetworks that are going to be selected between the population
		$elitism,                                       #Subnetworks to be keeped in the next population of subnetworks due to their high scores
		$network_type,                                  #Fixed (one network) or variable (two networks)
		$optimization,                                  #Multivariate or univariate
		$updating_scheme,			        			#Synchronus or asynchronus (attrators calculation)
		$f2,      				        				#Reference to the tab 3 frame
		$heap,					        				#Reference to the @_[HEAP]	
		$contextualization_b1,		                	#Reference to RUN button
		$contextualization_b2,			        		#Reference to STOP button
		$contextualization_fixed_total_iterations,      #Reference to the contain of Total Iterations entry box ##
		$contextualization_fixed_maximum_matching,$CIRCUITS_HASH) = @_;#Reference to the contain of Maximum Matching entry box ##
	
	#Multivariate sampling probability distribution: including interactions and circuits
	my %CIRCUITS = %{$CIRCUITS_HASH};
	my @network_view_content = @{$network_view_content_array};
	
	if ($optimization eq 'multivariate'){
		my @circuit_names = keys %CIRCUITS;
		
		foreach my $circuit_name(@circuit_names){
			my @circuit_lines= keys %{$CIRCUITS{$circuit_name}};
			$positive_circuits{$circuit_name}=\@circuit_lines;	
		}
	}
	else{}
	my %edges_frequences = ();
	my %positive_circuits_frequences;
	my %positive_circuits_percenteges;
	foreach my $edge(@network_view_content){
			$edges_frequences{$edge} = 0.5;
	}
	##Asignate a number to each gene of the list of genes that are first actualized: subroutine called encoding_updating_sequence
	my @updating_sequence_encoded = @{encoding_updating_sequence($updating_sequence_array)};
	
	###Iterative network prunning
	if ($network_type eq 'fixed'){
		my $best_score = 0;
		for (my $i=0;$i<$maximum_iterations;$i++){
			
			my @scores_of_population_of_networks_to_return_ordered = sort {$b<=>$a} @scores_of_population_of_networks_to_return;
			$best_score = $scores_of_population_of_networks_to_return[0];
			$contextualization_fixed_maximum_matching->configure('-textvariable'=>$best_score);
			$f2->update;#Very Important
			my $stop_button_state = $contextualization_b2->cget('-state');
			if ($stop_button_state eq 'disabled'){
				last
			}
			else{
				$heap->{counter_widget}++;
				$contextualization_fixed_total_iterations->configure('-textvariable'=>$heap->{counter_widget});
				$f2->update;#Very Important
				open FILEHANDLE, ">>./optimization_report.txt";
				print FILEHANDLE "\n\nIteration $i\n";
				close FILEHANDLE;
				my @iteration_results = @{iteration(
								$network_view_content_array,
								$phenotype_1_content_array,
								$phenotype_2_content_array,
								$fixed_interactions_1_content_array,
								$fixed_interactions_2_content_array,
								$updating_sequence_array,
								$maximum_iterations,
								$population,
								$selection,
								$elitism,
								$network_type,
								$optimization,
								$updating_scheme,
								\%edges_frequences,
								\@elite,
								\%positive_circuits_frequences,
								)};
								
				%edges_frequences = %{$iteration_results[0]};
				%positive_circuits_frequences =	%{$iteration_results[1]};			
			}
		}
	#in variable networks
	}elsif ($network_type eq 'variable'){
		my $best_score = 0;
		for (my $i=0;$i<$maximum_iterations;$i++){	
			my @scores_of_population_of_networks_to_return_ordered_network_1 = sort {$b<=>$a} @scores_of_population_of_networks_to_return_network_1;
			my $best_score = $scores_of_population_of_networks_to_return_network_1[0];
			$contextualization_fixed_maximum_matching->configure('-textvariable'=>$best_score); ##
			$f2->update;#Very Important
			my $stop_button_state = $contextualization_b2->cget('-state');
			if ($stop_button_state eq 'disabled'){
				$contextualization_b2->configure(-state=>'normal');
				$contextualization_b1->configure(-state=>'disabled');
				$f2->update;#Very Important
				last
			}else{
				$heap->{counter_widget}++;
				$contextualization_fixed_total_iterations->configure('-textvariable'=>$heap->{counter_widget}); ##
				$f2->update;#Very Important
				open FILEHANDLE, ">>./optimization_report.txt";
				print FILEHANDLE "\n\nIteration $i Network 1\n";
				close FILEHANDLE;
				my @iteration_results = @{iteration_differential_network_1(
											$network_view_content_array,
											$phenotype_1_content_array,
											#$phenotype_2_content_array,
											$fixed_interactions_1_content_array,
											#$fixed_interactions_2_content_array,
											$updating_sequence_array,
											$maximum_iterations,
											$population,
											$selection,
											$elitism,
											$network_type,
											$optimization,
											$updating_scheme,
											\%edges_frequences,
											\@elite,
											\%positive_circuits_frequences,
											)};
				%edges_frequences = %{$iteration_results[0]};
				%positive_circuits_frequences =	%{$iteration_results[1]};
			}										
		}
		push @results,\%edges_frequences;
		foreach my $edge(@network_view_content){
			$edges_frequences{$edge}=0.5;
		}
		@elite=();
		$heap->{counter_widget}=0;
		$contextualization_fixed_total_iterations->configure('-textvariable'=>$heap->{counter_widget}); ##
		$best_score = 0;
		$contextualization_fixed_maximum_matching->configure('-textvariable'=>$best_score); ##
		$f2->update;#Very Important
		for (my $i=0;$i<$maximum_iterations;$i++){	
			my @scores_of_population_of_networks_to_return_ordered_network_2 = sort {$b<=>$a} @scores_of_population_of_networks_to_return_network_2;
			my $best_score = $scores_of_population_of_networks_to_return_network_2[0];
			$contextualization_fixed_maximum_matching->configure('-textvariable'=>$best_score); ##
			$f2->update;#Very Important
			my $stop_button_state = $contextualization_b2->cget('-state');
			if ($stop_button_state eq 'disabled'){
				last
			}else{
				$heap->{counter_widget}++;
				$contextualization_fixed_total_iterations->configure('-textvariable'=>$heap->{counter_widget}); ##
				$f2->update;#Very Important
				open FILEHANDLE, ">>./optimization_report.txt";
				print FILEHANDLE "\n\nIteration $i Network 2\n";
				close FILEHANDLE;
				my @iteration_results = @{iteration_differential_network_2(
											$network_view_content_array,
											#$phenotype_1_content_array,
											$phenotype_2_content_array,
											#$fixed_interactions_1_content_array,
											$fixed_interactions_2_content_array,
											$updating_sequence_array,
											$maximum_iterations,
											$population,
											$selection,
											$elitism,
											$network_type,
											$optimization,
											$updating_scheme,
											\%edges_frequences,
											\@elite,
											\%positive_circuits_frequences)};
				%edges_frequences = %{$iteration_results[0]};
				%positive_circuits_frequences =	%{$iteration_results[1]};							
			}										
		}
		push @results,\%edges_frequences;
	}else{}
	
	###Subroutine for iterative network prunning in fixed networks
	sub iteration{
		my (
			$network_view_content_array,
			$phenotype_1_content_array,
			$phenotype_2_content_array,
			$fixed_interactions_1_content_array,
			$fixed_interactions_2_content_array,
			$updating_sequence_array,
			$maximum_iterations,
			$population,
			$selection,
			$elitism,
			$network_type,
			$optimization,
			$updating_scheme,
			$edges_frequences_hash,
			$elite_array,
			$positive_circuits_frequences_hash) = @_;
	
		my @network_view_content = @{$network_view_content_array};
		my @phenotype_1_content = @{$phenotype_1_content_array};
		my @phenotype_2_content = @{$phenotype_2_content_array};
		my @fixed_interactions_1_content = @{$fixed_interactions_1_content_array};
		my @fixed_interactions_2_content = @{$fixed_interactions_2_content_array};
		my @updating_sequence = @{$updating_sequence_array};
		my %edges_frequences = %{$edges_frequences_hash};
		#my %positive_circuits_percenteges=%{$positive_circuits_frequences_hash};
		my %positive_circuits_percenteges=();
		if ($optimization eq 'multivariate'){%positive_circuits_percenteges=%{$positive_circuits_frequences_hash}}
		else{}
		my @population_of_networks=();
		my @scores_of_population_of_networks=();
		##Selection of the interactions included in the next population of subnetworks
		for (my $i=0;$i<$population;$i++){
			my %configuration=();
			foreach my $interaction(@network_view_content){
				if (random_numbers_generator($edges_frequences{$interaction})){
					$configuration{$interaction}=$interaction;
				} else{}
			}
			
			#Interactions to be keeped from the 1 steady state
			foreach my $fixed_interaction(@fixed_interactions_1_content){
					$configuration{$fixed_interaction}=$fixed_interaction;
			}
			#Selection of the positive circuits included in the next population of subnetworks
			if ($optimization eq 'multivariate'){
				my @positive_circuit_names = keys %positive_circuits_percenteges;
				#circuits addition
				for (my $j=0;$j<@positive_circuit_names;$j++){
					my $circuit_frequency = $positive_circuits_percenteges{$positive_circuit_names[$j]};
					my $decision=random_numbers_generator_circuits($circuit_frequency);
					if ($decision==1){
						my @selected_circuit_interactions = @{$positive_circuits{$positive_circuit_names[$j]}};
						foreach my $selected_circuit_interaction(@selected_circuit_interactions){
							$configuration{$selected_circuit_interaction}=$selected_circuit_interaction;
						}
					}
					elsif ($decision==0){}
					else{}
				}
			}
			else{}
			push @population_of_networks, \%configuration;
		}
		#Keeping some subnetworks in the next population due to their high scores
		foreach my $pop_of_net(@population_of_networks){
			print $pop_of_net,"\n";
		}
		for (my $i=0;$i<@elite;$i++){
			$population_of_networks[$i]=$elite[$i];
		}
		
		my %initial_state_phenotype_1 = ();
		my %initial_state_phenotype_2 = ();
		
			for (my $i=0;$i<@phenotype_1_content;$i++){
				if ($phenotype_1_content[$i] =~ /(.+)(\s)([01])/){
					$initial_state_phenotype_1{$1} = $3;
				}else{}
			}
			for (my $i=0;$i<@phenotype_2_content;$i++){
				if ($phenotype_2_content[$i] =~ /(.+)(\s)([01])/){
					$initial_state_phenotype_2{$1} = $3;
				}else{}
			}
			
		my @genes_raw_network_ordered = sort {$a cmp $b} keys %initial_state_phenotype_1;	
		my $score_phenotype_1 = 0;
		my $score_phenotype_2 = 0;
		
		for (my $i=0;$i<@population_of_networks;$i++){
			$score_phenotype_1 = 0;
			$score_phenotype_2 = 0;
			my @network_lines = keys %{$population_of_networks[$i]};
			my %genes=();
			foreach my $line(@network_lines){
				if ($line =~ /(.+)(\s)(\-[\>\|])(\s)(.+)/){
					$genes{$1}=$1;
					$genes{$5}=$5;
				}else{}
			}
			
			my @genes_ordered = sort {$a cmp $b} keys %genes;
			my $initial_state_phenotype_1_string = '';
			my $initial_state_phenotype_2_string = '';
			
			foreach my $gene_ordered(@genes_ordered){
				if ((exists $initial_state_phenotype_1{$gene_ordered}) and (exists $initial_state_phenotype_1{$gene_ordered})){
					$initial_state_phenotype_1_string = $initial_state_phenotype_1_string.$initial_state_phenotype_1{$gene_ordered};
					$initial_state_phenotype_2_string = $initial_state_phenotype_2_string.$initial_state_phenotype_2{$gene_ordered};
				}
				else{
					$initial_state_phenotype_1_string = $initial_state_phenotype_1_string.'0';
					$initial_state_phenotype_2_string = $initial_state_phenotype_2_string.'0';
				}	
			}
			
			my @boolean_functions_view_content = @{Boolean_functions_generator::boolean_functions(\@network_lines)};
			my $number_of_genes = scalar @boolean_functions_view_content;
			my %genes_encoded = ();
			my @boolean_functions_formated = ();
			
			for (my $i=0;$i<@boolean_functions_view_content;$i++){
				if($boolean_functions_view_content[$i] =~ /(.+)(\t)(\=)(\s)(.+)/){
				    my $index_correction=$i+1;
				    $genes_encoded{$1} = "x$index_correction";
				    #print  $1,"\t",$index_correction,"\n";
				}else{}	
			}
			
			for (my $i=0;$i<@boolean_functions_view_content;$i++){
				if($boolean_functions_view_content[$i] =~ /(.+)(\t)(\=)(\s)(.+)/){
					my $index_correction=$i+1;
					my $raw_string=$5;
					print $raw_string,"\n";
					for (my $j=0;$j<@genes_ordered;$j++){
					    if ($raw_string =~ /(.*)([\!\(])($genes_ordered[$j])([\)\&\|])(.*)/){  # Iteratively substitute
						$raw_string = $1.$2.$genes_encoded{$3}.$4.$5;
					    }
					}
					print $raw_string,"\n";

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
					print $boolean_functions_formated[$i],"\n\n";

				}else{}
			}
			
			my @polynomial_functions = @{Translator::translator(\@boolean_functions_formated,$number_of_genes)};
			my @polynomial_functions_to_eval = ();
			
			foreach my $polynomial_function(@polynomial_functions){
				$polynomial_function =~ s/x(\d+)/\$x\[$1\]/g;
				$polynomial_function =~ /(.+)(\s\=\s)(.+)(\n)/;
				push @polynomial_functions_to_eval, $3;
			}
			
			my ($states_phenotype_1_array,$type_of_attractor_phenotype_1);
			my ($states_phenotype_2_array,$type_of_attractor_phenotype_2);
			
			if ($updating_scheme eq 'synchronous'){
				($states_phenotype_1_array,$type_of_attractor_phenotype_1) = trajectory_synchronous(\@polynomial_functions_to_eval,$initial_state_phenotype_1_string,$number_of_genes);
				($states_phenotype_2_array,$type_of_attractor_phenotype_2) = trajectory_synchronous(\@polynomial_functions_to_eval,$initial_state_phenotype_2_string,$number_of_genes);
			}elsif($updating_scheme eq 'asynchronous'){
				my @updating_sequence_new = ();
				foreach my $gene_ordered(@updating_sequence){
					if (exists $genes{$gene_ordered}){
						push @updating_sequence_new, $gene_ordered;
					}else{}
				}
				my @updating_sequence_encoded_new=@{encoding_updating_sequence(\@updating_sequence_new)};
				($states_phenotype_1_array,$type_of_attractor_phenotype_1) = trajectory_asynchronous(\@polynomial_functions_to_eval,$initial_state_phenotype_1_string,$number_of_genes,\@updating_sequence_encoded_new);
				($states_phenotype_2_array,$type_of_attractor_phenotype_2) = trajectory_asynchronous(\@polynomial_functions_to_eval,$initial_state_phenotype_2_string,$number_of_genes,\@updating_sequence_encoded_new);
			}else{}
			
			my @states_phenotype_1 = @{$states_phenotype_1_array};
			my @states_phenotype_2 = @{$states_phenotype_2_array};
			if ($type_of_attractor_phenotype_1 eq 'fixed point'){
					my @last_state_phenotype_1 = split //,$states_phenotype_1[(scalar @states_phenotype_1)-1];
					for (my $g=0;$g<@genes_ordered;$g++){
						$genes{$genes_ordered[$g]}=$last_state_phenotype_1[$g];
					}
					foreach my $grn(@genes_raw_network_ordered){
						if (exists $genes{$grn}){
							if ($genes{$grn} eq $initial_state_phenotype_1{$grn}){$score_phenotype_1=$score_phenotype_1+1}
						}else{}
					}
			}elsif($type_of_attractor_phenotype_1 eq 'limit cycle'){
					$score_phenotype_1=0;
			}else{}
			
			if ($type_of_attractor_phenotype_2 eq 'fixed point'){
					my @last_state_phenotype_2 = split //,$states_phenotype_2[(scalar @states_phenotype_2)-1];
					for (my $g=0;$g<@genes_ordered;$g++){
						$genes{$genes_ordered[$g]}=$last_state_phenotype_2[$g];
					}
					
					foreach my $grn(@genes_raw_network_ordered){
						if (exists $genes{$grn}){
							if ($genes{$grn} eq $initial_state_phenotype_2{$grn}){$score_phenotype_2=$score_phenotype_2+1}
						}else{}
					}
			}elsif($type_of_attractor_phenotype_2 eq 'limit cycle'){
					$score_phenotype_2=0;
			}else{}
			push @scores_of_population_of_networks, $score_phenotype_1+$score_phenotype_2;
		}
		
		my %scores = ();
		for (my $i=0;$i<@scores_of_population_of_networks;$i++){
			$scores{"configuration_$i"}=$scores_of_population_of_networks[$i];
		}
		
		my @population_of_networks_ordered = sort {$scores{$b} <=> $scores{$a}} keys %scores;
		open FILEHANDLE, ">>./optimization_report.txt";
		print FILEHANDLE "\n";
		foreach my $imprimir(@population_of_networks_ordered){
			print FILEHANDLE "$imprimir\: ".$scores{$imprimir}."\n";
		}	
		close FILEHANDLE;
		foreach my $edge(@network_view_content){
			$edges_frequences{$edge}=0;
		}
		for (my $i=0;$i<$selection;$i++){
			$population_of_networks_ordered[$i] =~ /^(configuration_)(\d+)/;
			for (my $j=0;$j<@network_view_content;$j++){
				if (exists $population_of_networks[$2]{$network_view_content[$j]}){
					$edges_frequences{$network_view_content[$j]}=$edges_frequences{$network_view_content[$j]}+1;
				}else{}
			}
		}
		foreach my $edge_frequence(keys %edges_frequences){
			$edges_frequences{$edge_frequence} = $edges_frequences{$edge_frequence}/$selection;
		}
		for (my $i=0;$i<$elitism;$i++){
			$population_of_networks_ordered[$i] =~ /^(configuration_)(\d+)/;
			$elite[$i] = $population_of_networks[$2];
		}
		my @selected_configurations=();
		for (my $i=0;$i<$selection;$i++){
			$population_of_networks_ordered[$i] =~ /^(configuration_)(\d+)/;
			$selected_configurations[$i]= $population_of_networks[$2];
		}
		%positive_circuits_percenteges = %{computing_circuits_percenteges(\@selected_configurations,$selection)};
		open FILEHANDLE, ">>./optimization_report.txt";
		foreach my $key(keys %edges_frequences){
			print "$key : $edges_frequences{$key} \n";print "$key : $edges_frequences{$key} \n";
			print FILEHANDLE "$key : $edges_frequences{$key} \n";print "$key : $edges_frequences{$key} \n";
		}
		close FILEHANDLE;
		@population_of_networks_to_return= @population_of_networks;
		@scores_of_population_of_networks_to_return=@scores_of_population_of_networks;
		my @iteration_return = (\%edges_frequences,\%positive_circuits_percenteges);
		return(\@iteration_return);
	}
	
	##Subroutine for iterative network prunning in variable networks
	sub iteration_differential_network_1{
		my (
			$network_view_content_array,
			$phenotype_1_content_array,
			#$phenotype_2_content_array,
			$fixed_interactions_1_content_array,
			#$fixed_interactions_2_content_array,
			$updating_sequence_array,
			$maximum_iterations,
			$population,
			$selection,
			$elitism,
			$network_type,
			$optimization,
			$updating_scheme,
			$edges_frequences_hash,
			$elite_array,
			$positive_circuits_frequences_hash) = @_;
		
		my @network_view_content = @{$network_view_content_array};
		my @phenotype_1_content = @{$phenotype_1_content_array};
		my @fixed_interactions_1_content = @{$fixed_interactions_1_content_array};
		#my @fixed_interactions_2_content = @{$fixed_interactions_2_content_array};
		my @updating_sequence = @{$updating_sequence_array};
		my %edges_frequences = %{$edges_frequences_hash};
		my %positive_circuits_percenteges=();
		if ($optimization eq 'multivariate'){%positive_circuits_percenteges=%{$positive_circuits_frequences_hash}}
		else{}
		my @elite = @{$elite_array};
		my @population_of_networks=();
		my @scores_of_population_of_networks=();
		for (my $i=0;$i<$population;$i++){
			my %configuration=();
			foreach my $interaction(@network_view_content){
				if (random_numbers_generator($edges_frequences{$interaction})){
					$configuration{$interaction}=$interaction;
				}else{}
			}
			foreach my $fixed_interaction(@fixed_interactions_1_content){
					$configuration{$fixed_interaction}=$fixed_interaction;
			}
			if ($optimization eq 'multivariate'){
				my @positive_circuit_names = keys %positive_circuits_percenteges;
				#circuits addition
				for (my $j=0;$j<@positive_circuit_names;$j++){
					my $circuit_frequency = $positive_circuits_percenteges{$positive_circuit_names[$j]};
					my $decision=random_numbers_generator_circuits($circuit_frequency);
					if ($decision==1){
						my @selected_circuit_interactions = @{$positive_circuits{$positive_circuit_names[$j]}};
						foreach my $selected_circuit_interaction(@selected_circuit_interactions){
							$configuration{$selected_circuit_interaction}=$selected_circuit_interaction;
						}
					}
					elsif ($decision==0){}
					else{}
				}
			}
			else{}
			push @population_of_networks, \%configuration;
		}
		for (my $i=0;$i<@elite_1;$i++){
			$population_of_networks[$i]=$elite_1[$i];
		}
		my %initial_state_phenotype_1 = ();
		my %initial_state_phenotype_2 = ();

			for (my $i=0;$i<@phenotype_1_content;$i++){
				if ($phenotype_1_content[$i] =~ /(.+)(\s)([01])/){
					$initial_state_phenotype_1{$1} = $3;
				}else{}
			}
		my @genes_raw_network_ordered = sort {$a cmp $b} keys %initial_state_phenotype_1;	
		my $score_phenotype_1 = 0;
		
		for (my $i=0;$i<@population_of_networks;$i++){
			$score_phenotype_1 = 0;
			my @network_lines = keys %{$population_of_networks[$i]};
			my %genes = ();
			foreach my $line(@network_lines){
				if ($line =~ /(.+)(\s)(\-[\>\|])(\s)(.+)/){
					$genes{$1}=$1;
					$genes{$5}=$5;
				}else{}
			}		
			my @genes_ordered = sort {$a cmp $b} keys %genes;
			my $initial_state_phenotype_1_string='';
			foreach my $gene_ordered(@genes_ordered){
				if ((exists $initial_state_phenotype_1{$gene_ordered}) and (exists $initial_state_phenotype_1{$gene_ordered})){
					$initial_state_phenotype_1_string = $initial_state_phenotype_1_string.$initial_state_phenotype_1{$gene_ordered};
					#$initial_state_phenotype_2_string = $initial_state_phenotype_2_string.$initial_state_phenotype_2{$gene_ordered};
				}
				else{
					$initial_state_phenotype_1_string = $initial_state_phenotype_1_string.'0';
					#$initial_state_phenotype_2_string = $initial_state_phenotype_2_string.'0';
				}
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

					###################  Newly added ###########################
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
			my ($states_phenotype_1_array,$type_of_attractor_phenotype_1);
			if ($updating_scheme eq 'synchronous'){
				($states_phenotype_1_array,$type_of_attractor_phenotype_1) = trajectory_synchronous(\@polynomial_functions_to_eval,$initial_state_phenotype_1_string,$number_of_genes);
			}	
			elsif($updating_scheme eq 'asynchronous'){
				my @updating_sequence_new=();
				foreach my $gene_ordered(@updating_sequence){
					if (exists $genes{$gene_ordered}){
						push @updating_sequence_new, $gene_ordered;
					}else{}
				}
				my @updating_sequence_encoded_new=@{encoding_updating_sequence(\@updating_sequence_new)};
				($states_phenotype_1_array,$type_of_attractor_phenotype_1) = trajectory_asynchronous(\@polynomial_functions_to_eval,$initial_state_phenotype_1_string,$number_of_genes,\@updating_sequence_encoded_new);
			}else{}
			my @states_phenotype_1 = @{$states_phenotype_1_array};
			if ($type_of_attractor_phenotype_1 eq 'fixed point'){
					my @last_state_phenotype_1 = split //,$states_phenotype_1[(scalar @states_phenotype_1)-1];
					for (my $g=0;$g<@genes_ordered;$g++){
						$genes{$genes_ordered[$g]}=$last_state_phenotype_1[$g];
					}
					
					foreach my $grn(@genes_raw_network_ordered){
						if (exists $genes{$grn}){
							if ($genes{$grn} eq $initial_state_phenotype_1{$grn}){$score_phenotype_1=$score_phenotype_1+1}
						}else{}
					}
			}elsif($type_of_attractor_phenotype_1 eq 'limit cycle'){
					$score_phenotype_1=0;
			}else{}
			push @scores_of_population_of_networks, $score_phenotype_1;
		}
		my %scores = ();
		for (my $i=0;$i<@scores_of_population_of_networks;$i++){
			$scores{"configuration_$i"}=$scores_of_population_of_networks[$i];
		}
		my @population_of_networks_ordered = sort {$scores{$b} <=> $scores{$a}} keys %scores;
		open FILEHANDLE, ">>./optimization_report.txt";
		print FILEHANDLE "\n";
		foreach my $imprimir(@population_of_networks_ordered){
			print FILEHANDLE "$imprimir\: ".$scores{$imprimir}."\n";
		}	
		close FILEHANDLE;
		
		foreach my $edge(@network_view_content){
			$edges_frequences{$edge}=0;
		}
		for (my $i=0;$i<$selection;$i++){
			$population_of_networks_ordered[$i] =~ /^(configuration_)(\d+)/;
			for (my $j=0;$j<@network_view_content;$j++){
				if (exists $population_of_networks[$2]{$network_view_content[$j]}){
					$edges_frequences{$network_view_content[$j]}=$edges_frequences{$network_view_content[$j]}+1;
				}else{}
			}
		}
		foreach my $edge_frequence(keys %edges_frequences){
			$edges_frequences{$edge_frequence} = $edges_frequences{$edge_frequence}/$selection;
		}
		for (my $i=0;$i<$elitism;$i++){
			$population_of_networks_ordered[$i] =~ /^(configuration_)(\d+)/;
			$elite_1[$i] = $population_of_networks[$2];
		}
		############################################################Circuits frequences
		my @selected_configurations=();
		for (my $i=0;$i<$selection;$i++){
			$population_of_networks_ordered[$i] =~ /^(configuration_)(\d+)/;
			$selected_configurations[$i]= $population_of_networks[$2];
		}
		%positive_circuits_percenteges = %{computing_circuits_percenteges(\@selected_configurations,$selection)};
		############################################################Circuits frequences end
		
		
		
		open FILEHANDLE, ">>./optimization_report.txt";
		print FILEHANDLE "Network 1\n";
		foreach my $key(keys %edges_frequences){
			print "$key : $edges_frequences{$key} \n";print "$key : $edges_frequences{$key} \n";
			print FILEHANDLE "$key : $edges_frequences{$key} \n";print "$key : $edges_frequences{$key} \n";
		}
		close FILEHANDLE;
		@population_of_networks_to_return_network_1= @population_of_networks;
		@scores_of_population_of_networks_to_return_network_1=@scores_of_population_of_networks;
		#return(\%edges_frequences);#Falta asynchronous, multivariate, differential..
		my @iteration_return = (\%edges_frequences,\%positive_circuits_percenteges);
		return(\@iteration_return);
	}
	
	##Subroutine for iterative network prunning in variable networks
	sub iteration_differential_network_2{
				
		my (
			$network_view_content_array,
			$phenotype_1_content_array,
			#$phenotype_2_content_array,
			$fixed_interactions_1_content_array,
			#$fixed_interactions_2_content_array,
			$updating_sequence_array,
			$maximum_iterations,
			$population,
			$selection,
			$elitism,
			$network_type,
			$optimization,
			$updating_scheme,
			$edges_frequences_hash,
			$elite_array,
			$positive_circuits_frequences_hash) = @_;
	
		my @network_view_content = @{$network_view_content_array};
		my @phenotype_1_content = @{$phenotype_1_content_array};
		my @fixed_interactions_1_content = @{$fixed_interactions_1_content_array};
		my @updating_sequence = @{$updating_sequence_array};
		my %edges_frequences = %{$edges_frequences_hash};
		my %positive_circuits_percenteges=();
		if ($optimization eq 'multivariate'){%positive_circuits_percenteges=%{$positive_circuits_frequences_hash}}
		else{my %positive_circuits_percenteges=()}
		my @elite = @{$elite_array};
		
		my @population_of_networks = ();
		my @scores_of_population_of_networks = ();
		for (my $i=0;$i<$population;$i++){
			my %configuration = ();
			foreach my $interaction(@network_view_content){
				if (random_numbers_generator($edges_frequences{$interaction})){
					$configuration{$interaction}=$interaction;
				}else{}
			}
			foreach my $fixed_interaction(@fixed_interactions_1_content){
					$configuration{$fixed_interaction}=$fixed_interaction;
			}
			if ($optimization eq 'multivariate'){
				my @positive_circuit_names = keys %positive_circuits_percenteges;
				#circuits addition
				for (my $j=0;$j<@positive_circuit_names;$j++){
					my $circuit_frequency = $positive_circuits_percenteges{$positive_circuit_names[$j]};
					my $decision=random_numbers_generator_circuits($circuit_frequency);
					if ($decision==1){
						my @selected_circuit_interactions = @{$positive_circuits{$positive_circuit_names[$j]}};
						foreach my $selected_circuit_interaction(@selected_circuit_interactions){
							$configuration{$selected_circuit_interaction}=$selected_circuit_interaction;
						}
					}
					elsif ($decision==0){}
					else{}
				}
			}else{}
			push @population_of_networks, \%configuration;
		}
		for (my $i=0;$i<@elite_2;$i++){
			$population_of_networks[$i]=$elite_2[$i];
		}
		my %initial_state_phenotype_1 = ();
		my %initial_state_phenotype_2 = ();

			for (my $i=0;$i<@phenotype_1_content;$i++){
				if ($phenotype_1_content[$i] =~ /(.+)(\s)([01])/){
					$initial_state_phenotype_1{$1} = $3;
				}else{}
			}
		my @genes_raw_network_ordered = sort {$a cmp $b} keys %initial_state_phenotype_1;	
		my $score_phenotype_1 = 0;
		for (my $i=0;$i<@population_of_networks;$i++){
			$score_phenotype_1=0;
			my @network_lines = keys %{$population_of_networks[$i]};
			my %genes=();
			foreach my $line(@network_lines){
				if ($line =~ /(.+)(\s)(\-[\>\|])(\s)(.+)/){
					$genes{$1}=$1;
					$genes{$5}=$5;
				}else{}
			}		
			my @genes_ordered = sort {$a cmp $b} keys %genes;
			my $initial_state_phenotype_1_string = '';
			#my $initial_state_phenotype_2_string='';

			foreach my $gene_ordered(@genes_ordered){
				if ((exists $initial_state_phenotype_1{$gene_ordered}) and (exists $initial_state_phenotype_1{$gene_ordered})){
						$initial_state_phenotype_1_string = $initial_state_phenotype_1_string.$initial_state_phenotype_1{$gene_ordered};
						#$initial_state_phenotype_2_string = $initial_state_phenotype_2_string.$initial_state_phenotype_2{$gene_ordered};
					}
					else{
						$initial_state_phenotype_1_string = $initial_state_phenotype_1_string.'0';
						#$initial_state_phenotype_2_string = $initial_state_phenotype_2_string.'0';
					}
			}
			
			my @boolean_functions_view_content = @{Boolean_functions_generator::boolean_functions(\@network_lines)};
			my $number_of_genes = scalar @boolean_functions_view_content;
			my %genes_encoded = ();
			my @boolean_functions_formated=();
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
					    #if ($ar[$i] =~ /^(\(\()(.+?)(\).*)$/){  
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
			my @polynomial_functions_to_eval =();
			foreach my $polynomial_function(@polynomial_functions){
				$polynomial_function =~ s/x(\d+)/\$x\[$1\]/g;
				$polynomial_function =~ /(.+)(\s\=\s)(.+)(\n)/;
				push @polynomial_functions_to_eval, $3;
			}
			my ($states_phenotype_1_array,$type_of_attractor_phenotype_1);
			#my ($states_phenotype_2_array,$type_of_attractor_phenotype_2);
			if ($updating_scheme eq 'synchronous'){
				($states_phenotype_1_array,$type_of_attractor_phenotype_1) = trajectory_synchronous(\@polynomial_functions_to_eval,$initial_state_phenotype_1_string,$number_of_genes);
				#($states_phenotype_2_array,$type_of_attractor_phenotype_2) = trajectory_synchronous(\@polynomial_functions_to_eval,$initial_state_phenotype_2_string,$number_of_genes);
			}elsif($updating_scheme eq 'asynchronous'){
				my @updating_sequence_new=();
				foreach my $gene_ordered(@updating_sequence){
					if (exists $genes{$gene_ordered}){
						push @updating_sequence_new, $gene_ordered;
					}else{}
				}
				my @updating_sequence_encoded_new=@{encoding_updating_sequence(\@updating_sequence_new)};
				($states_phenotype_1_array,$type_of_attractor_phenotype_1) = trajectory_asynchronous(\@polynomial_functions_to_eval,$initial_state_phenotype_1_string,$number_of_genes,\@updating_sequence_encoded_new);
				#($states_phenotype_2_array,$type_of_attractor_phenotype_2) = trajectory_asynchronous(\@polynomial_functions_to_eval,$initial_state_phenotype_2_string,$number_of_genes,\@updating_sequence_encoded_new);
			}else{}
			my @states_phenotype_1 = @{$states_phenotype_1_array};
			#my @states_phenotype_2 = @{$states_phenotype_2_array};
			if ($type_of_attractor_phenotype_1 eq 'fixed point'){
					my @last_state_phenotype_1 = split //,$states_phenotype_1[(scalar @states_phenotype_1)-1];
					for (my $g=0;$g<@genes_ordered;$g++){
						$genes{$genes_ordered[$g]}=$last_state_phenotype_1[$g];
					}
					
					foreach my $grn(@genes_raw_network_ordered){
						if (exists $genes{$grn}){
							if ($genes{$grn} eq $initial_state_phenotype_1{$grn}){$score_phenotype_1=$score_phenotype_1+1}
						}else{}
					}
			}elsif($type_of_attractor_phenotype_1 eq 'limit cycle'){
					$score_phenotype_1=0;
			}else{}
			push @scores_of_population_of_networks, $score_phenotype_1;
		}
		my %scores=();
		for (my $i=0;$i<@scores_of_population_of_networks;$i++){
			$scores{"configuration_$i"}=$scores_of_population_of_networks[$i];
		}
		my @population_of_networks_ordered = sort {$scores{$b} <=> $scores{$a}} keys %scores;
		open FILEHANDLE, ">>./optimization_report.txt";
		print FILEHANDLE "\n";
		foreach my $imprimir(@population_of_networks_ordered){
			print FILEHANDLE "$imprimir\: ".$scores{$imprimir}."\n";
		}	
		close FILEHANDLE;
		foreach my $imprimir(@population_of_networks_ordered){
			print "$imprimir\: ".$scores{$imprimir}."\n";
		}
		#my %edges_frequences=();
		foreach my $edge(@network_view_content){
			$edges_frequences{$edge}=0;
		}
		for (my $i=0;$i<$selection;$i++){
			$population_of_networks_ordered[$i] =~ /^(configuration_)(\d+)/;
			for (my $j=0;$j<@network_view_content;$j++){
				if (exists $population_of_networks[$2]{$network_view_content[$j]}){
					$edges_frequences{$network_view_content[$j]}=$edges_frequences{$network_view_content[$j]}+1;
				}else{}
			}
		}
		foreach my $edge_frequence(keys %edges_frequences){
			$edges_frequences{$edge_frequence} = $edges_frequences{$edge_frequence}/$selection;
		}
		for (my $i=0;$i<$elitism;$i++){
			$population_of_networks_ordered[$i] =~ /^(configuration_)(\d+)/;
			$elite_2[$i] = $population_of_networks[$2];
		}
		############################################################Circuits frequences
		my @selected_configurations=();
		for (my $i=0;$i<$selection;$i++){
			$population_of_networks_ordered[$i] =~ /^(configuration_)(\d+)/;
			$selected_configurations[$i]= $population_of_networks[$2];
		}
		%positive_circuits_percenteges = %{computing_circuits_percenteges(\@selected_configurations,$selection)};
		############################################################Circuits frequences end
		open FILEHANDLE, ">>./optimization_report.txt";
		print FILEHANDLE "Network 2\n";
		foreach my $key(keys %edges_frequences){
			print "$key : $edges_frequences{$key} \n";print "$key : $edges_frequences{$key} \n";
			print FILEHANDLE "$key : $edges_frequences{$key} \n";print "$key : $edges_frequences{$key} \n";
		}
		close FILEHANDLE;
		
		@population_of_networks_to_return_network_2= @population_of_networks;
		@scores_of_population_of_networks_to_return_network_2=@scores_of_population_of_networks;
		#return(\%edges_frequences);
		my @iteration_return = (\%edges_frequences,\%positive_circuits_percenteges);
		return(\@iteration_return);
	}
	if ($network_type eq 'variable'){
		return(\@results,\@population_of_networks_to_return_network_1,\@scores_of_population_of_networks_to_return_network_1,\@population_of_networks_to_return_network_2,\@scores_of_population_of_networks_to_return_network_2);
	}elsif($network_type eq 'fixed'){
		return(\%edges_frequences,\@population_of_networks_to_return,\@scores_of_population_of_networks_to_return);#Falta asynchronous, multivariate, differential..
	}else{return 1}	
	exit();
}

###Subroutine to select the interactions included in the next population of subnetworks: Random link generator
sub random_numbers_generator{
	my ($internal) = @_;
	my $random_number = rand;
	my $random_state;
	if ($internal<0.2){$internal=0.2}
	elsif ($internal>0.8){$internal=0.8}
	else{}
	#$internal=0.5;#To make it random
	if ($random_number <= $internal){$random_state = 1}
	else {$random_state = 0}
	return $random_state;
}
sub random_numbers_generator_circuits{
	my ($internal) = @_;
	my $random_number = rand;
	my $random_state;
	if ($internal<0.2){$internal=0.2}
	elsif ($internal>0.8){$internal=0.8}
	else{}
	#$internal=0.5;#To make it random
	if ($random_number <= $internal){$random_state = 1}
	else {$random_state = 0}
	return $random_state;
}

sub trajectory_synchronous{
		my ($functions_array,$initial_state_string,$number_of_nodes) =@_;
		my $number_of_possible_values=2;
		my $longest_trajectory = 1000;#####################################################################################################################I have to limit this somehow
		#print "$initial_state_string es initial state\n";
		my @functions = @{$functions_array};
		#print "these are the functions: @functions";
		my @y = split //, $initial_state_string;
		#print "esto es y @y";
		my @x = (0,@y);
		#print "esto es x @x";
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
		#print "$initial_state_string es initial state\n";
		my @functions = @{$functions_array};
		#print "these are the functions: @functions";
		my @y = split //, $initial_state_string;
		#print "esto es y @y";
		my @x = (0,@y);
		#print "esto es x @x";
		my @y_index_correction = (0);
		my %states_visited=();
		my %perturbations = ();
		my @states;
		foreach (0..$longest_trajectory){
				my @y_index_correction = (0);
				#my $i = 1;
				for (my $j=0;$j<@updating_order_indexes;$j++){
					my $a = $updating_order_indexes[$j];
					
					if (exists $perturbations{$a}){
						$y[$a+1] = $perturbations{$a};
						$x[$a+1] = $y[$a+1];
						$y_index_correction[$a+1] = $y[$a+1];
						#$i++;
					}else{
						$y[$a+1] = eval($functions[$a]) % $number_of_possible_values;
						$x[$a+1] = $y[$a+1];
						$y_index_correction[$a+1] = $y[$a+1];
						#$i++;
					}
				}
				shift (@y_index_correction);
				my $new_state = join "", @y_index_correction;
				if (exists $states_visited{$new_state}){
					my $last_visited_state = $states[(scalar @states)-1];
					if ($last_visited_state eq $new_state){
						#print "This is a fixed point...\n";
						#$fixed_points{$new_state}=$new_state;
						push @states, $new_state;
						return(\@states,'fixed point');
					}else{
						#print "This is a limit cycle...\n";
						push @states, $new_state;
						return(\@states,'limit cycle');
					}
				}else{
					$states_visited{$new_state}=$new_state;
					print "$new_state\n";
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
sub computing_circuits_percenteges{
	my ($array_top_configuration_scores,$selection_number)=@_;
	my @top_configuration_scores= @{$array_top_configuration_scores};
	my %percenteges_circuits=();
	my @positive_circuits_names = keys %positive_circuits;
	for (my $i=0;$i<@positive_circuits_names;$i++){
		my @circuit_lines= @{$positive_circuits{$positive_circuits_names[$i]}};
		my %circuit_lines=();
		foreach my $circuit_line(@circuit_lines){
			$circuit_lines{$circuit_line}=$circuit_line;
		}
		my $circuit_percentege=0;
		foreach my $j(@top_configuration_scores){
			my %buffer=%{$j};
			my @configuration_lines=keys %{$j};
			my %configuration_lines=();
			foreach my $configuration_line(@configuration_lines){
				$configuration_lines{$configuration_line}=$configuration_line;
			}
			#Now I have two hashes: %circuit_lines and %configuration_lines. I have to check if all the elements in %circuit_lines are elements in %configuration_lines
			my $missing_interaction_flag='0';
			foreach my $circuit_or_configuration_line(keys %circuit_lines){
				if (exists $configuration_lines{$circuit_or_configuration_line}){}
				else{$missing_interaction_flag=1}
			}
			if ($missing_interaction_flag eq '0'){$circuit_percentege = $circuit_percentege+1}
			elsif($missing_interaction_flag eq '1'){}
			else{}
		}
		$circuit_percentege = $circuit_percentege/$selection_number;	
		$percenteges_circuits{$positive_circuits_names[$i]}=$circuit_percentege;
	}
	return (\%percenteges_circuits);
}
1
