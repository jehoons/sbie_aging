#!/usr/bin/perl

#### Ana Rodriguez Sanchez-Archidona, Luxembourg 15/7/2013, Luxembourg Center for Systems Biomedicine ####
#### Modul that includes the subroutine for the raw network initial format checking
#### and the initial formating of the literature gene regulatory raw network and experimental gene expression data

use strict;
use warnings;

package  Check_network_format;

#Subroutine for the raw network initial format checking
sub check_network_format{
	my ($network_view_content_array) = @_;
	my @network_view_content = @{$network_view_content_array};
	my $line_index_correction;
	my %nodes = ();
	my %edges = ();
	my @format_errors = ();
	my $number_nodes = 0;
	my $number_edges = 0;
	my $number_activations = 0;
	my $number_inhibitions = 0;
	for (my $i=0;$i<@network_view_content;$i++){
		$line_index_correction = $i+1;
		if ((!exists $edges{$network_view_content[$i]})and($network_view_content[$i] =~ /(.+)(\s)(\-[\>\|])(\s)(.+)(\n*)/)){
			$number_edges++;
			$nodes{$1} = $1;
			$nodes{$5} = $5;
			$edges{$network_view_content[$i]} = $network_view_content[$i];
			if ($3 eq '->'){$number_activations++}
			elsif ($3 eq '-|'){$number_inhibitions++}
			else{}
		} elsif(exists $edges{$network_view_content[$i]}){
			push @format_errors, "Line $line_index_correction repeated:\n$network_view_content[$i]";
		} else{
			push @format_errors, "Line $line_index_correction:\n$network_view_content[$i]\n";
		}
	}
	$number_nodes = scalar keys %nodes;
	my $number_errors = scalar @format_errors;
	return ($number_nodes,$number_edges,$number_activations,$number_inhibitions,$number_errors,\@format_errors);
}

#Subroutine for the expression values initial format checking
sub check_phenotype_format{
	my ($phenotype_array_view_content) = @_;
	my @phenotype_content_array = @{$phenotype_array_view_content};
	my $line_correction;
	my %genes = ();
	my @errors_array = ();
	my $number_genes = 0;

	for (my $i=0;$i<@phenotype_content_array;$i++){
		$line_correction = $i+1;
		if ((!exists $genes{$phenotype_content_array[$i]}) and ($phenotype_content_array[$i] =~ /^(\w+\-*\w*)(\s+)(\w+\-*\w*)/)){
			$number_genes++;
			$genes{$phenotype_content_array[$i]} = $phenotype_content_array[$i];
		} elsif(exists $genes{$phenotype_content_array[$i]}){
			push @errors_array, "Line $line_correction repeated:\n$phenotype_content_array[$i]";
		} else{
			push @errors_array, "Line $line_correction:\n$phenotype_content_array[$i]\n";
		}
	}
	$number_genes = scalar keys %genes;
	my $errors_number = scalar @errors_array;
	return ($number_genes,$errors_number,\@errors_array);
}

#Formating subroutine of the raw network
sub format_checking {
    
    my ($raw_input) = @_;
    
    my @input_raw_network = @{$raw_input};

    #Check if the lines of the file follow a pattern of 3 columns and remove empty and new lines
    my@final_input = ();
    
    foreach my$line (@input_raw_network) {
        if ($line=~/^(\w+\-*\w*)(\s)(\-[\>\|])(\s)(\w+\-*\w*)/) {
            my$new_line= $1.$2.$3.$4.$5."\n";
            $new_line =~ s/\t/ /;
        push (@final_input,$new_line);    
        } 
    }
  
    #Check and remove repeated interations
    my %nonrepeated_file;

    foreach my$line (@final_input) {
    $nonrepeated_file{$line} = $line;
    }
    
    return (keys %nonrepeated_file);  
}

#Formating subroutine for expression values
sub expression_format {
    
    my ($gene_expression)= @_;
    
    my @gene_expression_values= @{$gene_expression}; 
    
    #Check if the lines of the file follow a pattern of 3 columns and remove empty and new lines
    my @final_array = ();
    
    foreach my$line (@gene_expression_values) {
        if ($line=~/^(\w+\-*\w*)(\s+)(\w+\-*\w*)/) {
        push @final_array,$1.$2.$3."\n";    
        } 
    }

    #Check and remove repeated iterations
    my %nonrepeated_expression_values;

    foreach my$gene_expression_line (@final_array) {
    $nonrepeated_expression_values{$gene_expression_line} = $gene_expression_line;
    }

    return (keys %nonrepeated_expression_values); 
}
1