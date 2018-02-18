####################################################################################################################
####                                                                                              			    ####
####            Ganna Androsova, Luxembourg 23/08/2013, Luxembourg Center for Systems Biomedicine      			####
####                                                                                              			    ####
####                              Module for network visualization in XPRED 6.0.                   	    	    ####
####                                                                                             			    ####
#### Input: Gene interactions in format of three columns. First column and third columns are names of genes		####
#### which interact, second column contains type of interaction: '->' for activation and '-|' for inhibition.	####
#### Array should contain only white-space separation between its elements.										####
#### For example: @array_of_genes = (gene_A -> gene_B gene_B -| gene_C).										####
####																											####
#### This module includes: 																						####
#### 	- construction of the gene in form of nodes;															####
#### 	- construction of interactions in form of edges;														####
#### 	- mutual interactions are represented by curved edges;													####
#### 	- one-way interactions are represented by straight edges;												####
#### 	- center and adjust the view for full network visualization on canvas;									####
#### 	- layouts, which are dependent on node amount: Alphabetical Grid, Alphabetical Circular, 				####
####      Node Degree Circular;																					####
#### 	- functions for changing edge width, node color, font, zoom in/out;										####
#### 	- visualization of edge frequencies, which are calculated and parsed in XPRED 6.0;						####
#### 	- highlighting of the current node under the mouse;														####
#### 	- selecting of the specific nodes by clicking on them with pressed <Ctrl>;								####
#### 	- selecting of the specific nodes by lasso, while holding pressed <Ctrl>;								####
#### 	- de-selection of all the nodes by clicking on canvas free space;										####
#### 	- de-selection of specific node by clicking on it with pressed <Ctrl>;									####
#### 	- two modes of node manipulation: all nodes or only selected ones, (depending on the mode, it is 		####
#### 	  enabled to apply layout for all nodes or only selected and move only selected nodes or any of them;	####
#### 	- visualization of UP (red color) and DOWN (green color) regulated genes, UNDEFINED regulation is 		####
####	  represented by white color;																			####
####	- node color gradation depending on statistical frequency value, calculated and parsed in XPRED 6.0.	####
####																											####
#### Contact: ganna.androsova@gmail.com																			####
####################################################################################################################

package Visualize_network;

sub visualize_network{

use strict;
use warnings;
use Tk;
use Tk::BrowseEntry;
use Tk::Entry;
use Tk::Balloon;

#-------------------------------------------#
#                Read in data               #
#-------------------------------------------#
	my ($canv, $elements_array, $fixed_canvas) =  @_;         
	my @elements = @{$elements_array};	
	
	#Splitting of elements into new array
	my @elements_array = ();
	foreach $_ (@elements){
		chomp $_;
		my @element = split / /, $_;
		push @elements_array, @element;
	}
	@elements = @elements_array;
	@elements_array = ();
	
#-------------------------------------------#
#       Separating out unique nodes	        #
#-------------------------------------------#
	#Filtering the elements and creating array containing only nodes (without interactions)	
	my %nodes = ();
	my $index=0;
	$nodes{$elements[$index]} = $elements[$index];
	$nodes{$elements[$#elements]} = $elements[$#elements];
	
	for (my $index=2; $index <= $#elements-1; $index+=2){
		$nodes{$elements[$index]} = $elements[$index];
		$nodes{$elements[$index+1]} = $elements[$index+1];
		$index++;
	}
	
#-------------------------------------------#
#    Sorting nodes in alphabetical order    #
#-------------------------------------------#	
	my @sorted_nodes = ();	
	foreach my $key(keys %nodes){
		push @sorted_nodes, $key;
	} 
	@sorted_nodes = sort {$a cmp $b} @sorted_nodes;
	
#-------------------------------------------#	
#        Place objects on canvas            #
#-------------------------------------------#
	# Saving different types of objects into new arrays (ovals, text, groups of items)
	my %uniq_nodes =();
	#if there were previously any objects - clear them
	$canv->delete('all');
	
	#Set initial font parameters
	my $family = 'Arial Unicode MS'; 
	my $size = '11'; 
	my $weight = 'normal';
	my $layout_choise = 'Alphabetical Grid';
	my $type_of_mode = 'canvas';
	
	my ($elements_array_1, $text_items, $group_items, $uniq_nodes) = create($canv, $family, $size, $weight, @sorted_nodes);
	@elements_array = @{$elements_array_1}; 
	my @text_items = @{$text_items}; 
	my @work_items = @{$group_items};
	%uniq_nodes = %{$uniq_nodes};
	
	#Getting edges from the array
	my @for_edges_construction = create_edge_structure(@elements);

	#Applying default layout and scaling
	&alphabet_grid_layout($canv, \%uniq_nodes, \@for_edges_construction, @work_items);
	$size = &adjust_scale($canv, $fixed_canvas);
	return ($elements_array_1, $group_items, $uniq_nodes, \@for_edges_construction, $family, $size, $weight, $layout_choise, $type_of_mode);
}

#-------------------------------------------#
#      Creating the nodes on the canvas     #
#-------------------------------------------#
sub create{
	my ($canv, $family, $size, $weight, @sorted_nodes) = @_;
	my @elements_array = ();
	my @text_items = ();
	my @group_items = ();
	my %uniq_nodes = ();
	
	#Creating the nodes with default options
	foreach $index(@sorted_nodes) {
		my $oval = $canv->createOval(
									0,0,60,60,
									-tags => [$index, 'oval']);
		push @elements_array, $oval;
		
		#Placing the names of the nodes
		my $text = $canv->createText(
									30, 30, 
									-text =>$index, 
									-font => [-family =>$family, -size => $size, -weight => $weight], 
									-tags=> $index );
		push @text_items, $text;
			
		#Grouping of the nodes with their names, to apply specific action simultaneously
		my $group = $canv->createGroup([0, 0], -members => [$oval, $text], -tags => [$index, 'node']);
			
		#Caving the grouped items into new array
		$uniq_nodes{$index} = $group;
		push @group_items, $group;	
	}
	return (\@elements_array, \@text_items, \@group_items, \%uniq_nodes);
}
	
#-------------------------------------------#
#   Determine the availability of layouts   #
#-------------------------------------------#
sub determine_layouts{
	my $number_of_nodes = shift;
	
	#Depending on amount of nodes, there are avaiable specific layouts
	if ($number_of_nodes < 6){
		@layouts = ('Alphabetical Grid', 'Alphabetic Circular');
	}
	elsif ($number_of_nodes >= 6){
		@layouts = ('Alphabetical Grid', 'Alphabetic Circular', 'Node Degree Circular');
	}	
	return @layouts;
}	

#-------------------------------------------#
#               Apply layouts               #
#-------------------------------------------#	
sub apply_layout{
	my ($canv, $layout_choise, $type_of_mode, $uniq_nodes, $edges, $group_items) = @_;
	my @work_items = @{$group_items};
	my @sorted_groups_by_degree = ();
	
	#Detecting which layout is chosen from the menu
	my $initial = 0;
	my $passed_circ_radius = 0;
	my $x_middle_point = 0; 
	my $y_middle_point = 0;	
	my $size;

	if ($layout_choise eq 'Alphabetical Grid') {
		if ($type_of_mode eq 'canvas'){
			&alphabet_grid_layout($canv, $uniq_nodes, $edges, @work_items);
			$size = &adjust_scale($canv);
		}
		if ($type_of_mode eq 'selected'){
			my @selected = $canv->find(qw/withtag selected/);
			if (@selected == 0){
				&alphabet_grid_layout($canv, $uniq_nodes, $edges, @work_items);
				$size = &adjust_scale($canv);
			}else{
				my @sorted_selected = sort { $a <=> $b } @selected;
				&alphabet_grid_layout($canv, $uniq_nodes, $edges, @sorted_selected);
				$size = &adjust_scale($canv);
			}
		}
	}
	if ($layout_choise eq 'Alphabetic Circular') {
		if ($type_of_mode eq 'canvas'){
			&circular_layout($canv, $uniq_nodes, $edges, $initial, $passed_circ_radius, $x_middle_point, $y_middle_point, @work_items);
			$size = &adjust_scale($canv);
		}	
		if ($type_of_mode eq 'selected'){
			my @selected = $canv->find(qw/withtag selected/);
			if (@selected == 0){
				&circular_layout($canv, $uniq_nodes, $edges, $initial, $passed_circ_radius, $x_middle_point, $y_middle_point, @work_items);
				$size = &adjust_scale($canv);
			}else{
				my @sorted_selected = sort { $a <=> $b } @selected;
				&circular_layout($canv, $uniq_nodes, $edges, $initial, $passed_circ_radius, $x_middle_point, $y_middle_point, @sorted_selected);
				$size = &adjust_scale($canv);
			}
		}
	}
	if ($layout_choise eq 'Node Degree Circular') {
		my ($ref_node_degree, @sorted_groups_by_degree) = node_degree($canv);	
		if ($type_of_mode eq 'canvas'){
			&circular_layout($canv, $uniq_nodes, $edges, $initial, $passed_circ_radius, $x_middle_point, $y_middle_point, @sorted_groups_by_degree);
			$size = &adjust_scale($canv);
		}
		if ($type_of_mode eq 'selected'){
			my @selected = $canv->find(qw/withtag selected/);
			if (@selected == 0){
				&circular_layout($canv, $uniq_nodes, $edges, $initial, $passed_circ_radius, $x_middle_point, $y_middle_point, @sorted_groups_by_degree);
				$size = &adjust_scale($canv);
			}else{
				my @sorted_selected = sort { $a <=> $b } @selected;
				&circular_layout($canv, $uniq_nodes, $edges, $initial, $passed_circ_radius, $x_middle_point, $y_middle_point, @sorted_selected);
				$size = &adjust_scale($canv);
			}
		}
	}	
	return $size;
}

#-------------------------------------------#	
#        Alphabetical grid layout           #
#-------------------------------------------#
sub alphabet_grid_layout {
	my ($canv, $ref_uniq_nodes, $ref_array_edges, @working_items) = @_;
	
	#Setting the default values
	my $x1 = 30; my $y1 = 30;
	my $x2 = 90; my $y2 = 90;
	
	#Calculating the homogeneous distribution of all the nodes
	my @members = $canv->itemcget($working_items[0], -members); 
	my @coord_node = $canv->coords($members[0]);
	my $node_radius = ($coord_node[2]-$coord_node[0])/2;
	my $canv_distribution = ((int(sqrt(abs($#working_items))))*$node_radius*3.5);

	#Specify initial coordinates
	foreach my $group_item(@working_items){
		$canv->coords($group_item, 0,0);
		$canv->move($group_item=>$x1, $y1);
		$x1= $x1+$node_radius*3;
		if ($x1+$node_radius*2 >= $canv_distribution){
			$x1 = 30; $y1 = $y1+$node_radius*3;
		}
	}	
	#Refresh edges
	$canv->addtag('TEMP', withtag => "edge");
	$canv->delete('TEMP');
	&call_for_edges($canv, $ref_uniq_nodes, $ref_array_edges);
}
	
#-------------------------------------------#	
#            Circular layout                #
#-------------------------------------------#	
sub circular_layout{
	my ($canv, $uniq_nodes, $edges, $initial, $passed_circ_radius, $x_middle_point, $y_middle_point, @group_items) = @_;
	
	#Setting the default values
	my $x_circ; my $y_circ;
	my $circ_radius;
	my $Pi = 3.1416;
			
	#Calculate the angle at which all elements can be placed
	my $slice = 2*$Pi/($#group_items+1); 
	
	#Get the radius of the nodes
	my @members = $canv->itemcget($group_items[0], -members); 
	my @coord_node = $canv->coords($members[0]);
	my $node_radius = ($coord_node[2]-$coord_node[0])/2;
			
	#Calculate distance between the nodes in circular layout
	my $distance = $node_radius*3;
		
	if ($initial eq 'degree'){
		#Pass the radius of the circular layout (if defined)
		$circ_radius = $passed_circ_radius;
		my $compare_radius = $distance/sin($slice);
		if ($compare_radius > $circ_radius){
			$circ_radius = $compare_radius;
		}
	}else{
		#Calculate the radius of the circular layout
		$circ_radius = $distance/sin($slice);
				
		#Calculate the middle points of the circular layout
		my $x_middle_point = $circ_radius;
		my $y_middle_point = $x_middle_point;
	}
			
	#Return to initial coordinates
	foreach my $group_item(@group_items){
		$canv->coords($group_item, $x_middle_point,$y_middle_point);
	}		
	for (my $i = 0; $i<=$#group_items; $i++){
		my $angle = $slice*$i;
		my $x = int($circ_radius * cos($angle) + $node_radius + $x_middle_point);
		my $y = int($circ_radius * sin($angle) + $y_middle_point);
		$canv->move($group_items[$i]=>$x, $y);
	}	
	$canv->addtag('TEMP', withtag => "edge");
	$canv->delete('TEMP');
	&call_for_edges($canv, $uniq_nodes, $edges);
	$initial = 0;
}

#-------------------------------------------#	
#         Determine Node Degrees            #
#-------------------------------------------#		
sub node_degree{
	my ($canv) = shift;
	my @sorted_groups_by_degree = ();
	
	#Find all edges and count interactions for each node
	my %node_degree = ();
	my @edges = $canv->find(withtag => "edge");
	foreach my $edge_id (@edges){
		my @edge_tags = $canv->gettags($edge_id);
		foreach my $tag (@edge_tags){
			if ($tag ne 'edge' and $tag ne 'curved'){
				$node_degree{$tag}++;
			}
		}
	}
		
	#Sort node by decreasing degree of interactions
	my @sorted_node_degree = decrease_sorting(%node_degree);
	
	sub decrease_sorting {
		my %h = @_;
		return sort {$h{$b} <=> $h{$a}} keys %h;
	}

	foreach my $sorted_node(@sorted_node_degree){
		my @items_with_tag = $canv->find(withtag =>$sorted_node);
		foreach my $item (@items_with_tag){
			my @items_type = $canv->gettags($item);
			if ( grep/node/,@items_type){
				push @sorted_groups_by_degree, $item;
			}
		}	
	}
	return (\%node_degree, @sorted_groups_by_degree);
}

#-------------------------------------------#
#        Reading in the edge structure      #
#-------------------------------------------#
sub create_edge_structure{ 
	my (@work_items) = @_;
	my @for_edges_construction = ();
	for (my $index = 0; $index <= $#work_items-2; $index+=2){
		push @for_edges_construction, $work_items[$index].$work_items[$index+1].$work_items[$index+2];
		$index++;
	}
	return @for_edges_construction;
}
	
#-----------------------------------------------#
# Passing the edge structure for classification #
#-----------------------------------------------#	
sub call_for_edges{
	my ($canv, $hash_of_nodes, $array_of_edges) = @_;
	my @array_of_edges = @{$array_of_edges};
	my @copy_of_array_of_edges = @array_of_edges;
	my %uniq_nodes = %{$hash_of_nodes};
	my $group_from; my $group_to;
	foreach my $interact (@array_of_edges){
		$interact =~ /(.+)(\-[\>\|])(.+)/;
		if (grep/($3)(\-[\>\|])($1)/,@copy_of_array_of_edges){
			add_edge($interact, 1, $canv, $hash_of_nodes);
		}else{
			add_edge($interact, 0, $canv, $hash_of_nodes);
		}	
	}
	foreach my $node (values %uniq_nodes){
		$canv->raise($node, 'all');
	}
}

#-------------------------------------------#
#         Classification of edges           #
#-------------------------------------------#
sub add_edge{
	my ($interact, $type_of_edge, $canv, $hash_of_edges)=@_;
	
	my %uniq_nodes = %{$hash_of_edges};
	my @arrow_shape = ();
	my $group_from; my $group_to;
	
	$interact =~ /(.+)(\-[\>\|])(.+)/;
	if ($1 eq $3){
		if ($2 eq '-|'){
				@arrow_shape = (1, 1, 7);
		}
		if ($2 eq '->'){
				@arrow_shape = (7, 10, 5);
		}
		foreach my $node1(keys %uniq_nodes){
			if ($node1 eq $1){
				self_loop($canv, $uniq_nodes{$node1}, @arrow_shape);
			}
		}
	} else {
		foreach my $node1(keys %uniq_nodes){
			if ($node1 eq $1){
				$group_from = $uniq_nodes{$node1};
			}
			if ($node1 eq $3){
				$group_to = $uniq_nodes{$node1};
			}
			if ($2 eq '-|'){
				@arrow_shape = (3, 1, 7);
			}
			if ($2 eq '->'){
				@arrow_shape = (7, 10, 5);
			}
		}
		create_edge ($canv, $type_of_edge, $group_from, $group_to, @arrow_shape);
	}
}

#-------------------------------------------#
#    Constructing self-loop interactions    #
#-------------------------------------------#
sub self_loop{
	my ($canv, $group, @arrow_shape) = @_;
	
	my @member = $canv->itemcget($group, -members); 
	my @coord_oval = $canv->coords($member[0]);
	my $radius = ($coord_oval[2]-$coord_oval[0])/2;
	my @coord_text = $canv->coords($member[1]);
	my @tag = $canv->gettags($member[1]);

	my $x1 = $coord_text[0];
	my $y1 = $coord_text[1]-$radius;
	my @first = ($x1, $y1);

	my $x2 = $coord_text[0]-$radius/2;
	my $y2 = $y1 - $radius/2;
	my @second = ($x2, $y2);

	my $x3 = $x2 - $radius;
	my $y3 = $coord_text[1]-$radius;
	my @third = ($x3, $y3);

	my $x4 = $coord_text[0] - $radius - 2;
	my $y4 = $coord_text[1];
	my @fourth = ($x4, $y4);

	my $edge = $canv->createLine(@first, @second, @third, @fourth, -smooth => 1,  -arrow => 'last', -fill=> 'snow4',
				-arrowshape => [@arrow_shape],  -splinesteps => 20, -tags => [@tag, 'edge']);
	&main::set_arrow_width;			
	return $edge;
}
	
#-------------------------------------------#
#     Constructing plane and curved edges   #
#-------------------------------------------#
sub create_edge {
	my ($canv, $type_of_edge, $group_from, $group_to, @arrow_shape) = @_;
	my $x_from; my $y_from;
	my $x_from_shifted; my $y_from_shifted;
	my $x_to; my $y_to;
	my $x_to_shifted; my $y_to_shifted;
	my $x_btw; my $y_btw;
		
	my @member1 = $canv->itemcget($group_from, -members); 
	my @coord_rad1 = $canv->coords($member1[0]);
	my $radius1 = ($coord_rad1[2]-$coord_rad1[0])/2+2;
	my @coord_middle1 = $canv->coords($member1[1]);
	my @tag1 = $canv->gettags($member1[1]);
					
	my @member2 = $canv->itemcget($group_to, -members); 
	my @coord_rad2 = $canv->coords($member2[0]);
	my $radius2 = ($coord_rad2[2]-$coord_rad2[0])/2+2;
	my @coord_middle2 = $canv->coords($member2[1]);
	my @tag2 = $canv->gettags($member2[1]);
	
	my $B_length = abs($coord_middle2[1]-$coord_middle1[1]);
	my $C_length = abs($coord_middle2[0]-$coord_middle1[0]);
	my $A_length = sqrt(($coord_middle2[0]-$coord_middle1[0])**2+($coord_middle2[1]-$coord_middle1[1])**2);
	
	my $cos_angle_L;
	my $cos_angle_B;
	
	if ($A_length != 0){	
		$cos_angle_L = $B_length/$A_length;
		$cos_angle_B = $C_length/$A_length;
	} elsif ($A_length == 0){
		$cos_angle_L = 0;
		$cos_angle_B = 0;
	}
	
	my $B1_length = $radius1*$cos_angle_L-1;
	my $C1_length = $radius1*$cos_angle_B-1;
	my $B2_length = $radius2*$cos_angle_L;
	my $C2_length = $radius2*$cos_angle_B;
	
	# Calulation of arrow position depending from different coordinate posibilities
	if ($coord_middle2[0]>$coord_middle1[0] and $coord_middle1[1]>$coord_middle2[1]){
		$x_from = $coord_middle1[0] + $C1_length;
		$y_from = $coord_middle1[1] - $B1_length;
		$x_to = $coord_middle2[0] - $C2_length;
		$y_to = $coord_middle2[1] + $B2_length;
		$x_btw = ($x_to - $x_from)/2 + $x_from + 15;
		$y_btw = ($y_from - $y_to)/2 + $y_to;
		$x_from_shifted = $x_from + 5;
		$y_from_shifted = $y_from - 2;
		$x_to_shifted = $x_to + 3;
		$y_to_shifted = $y_to + 3;
	}
	elsif ($coord_middle1[0]>$coord_middle2[0] and $coord_middle1[1]>$coord_middle2[1]){
		$x_from = $coord_middle1[0] - $C1_length;
		$y_from = $coord_middle1[1] - $B1_length;
		$x_to = $coord_middle2[0] + $C2_length;
		$y_to = $coord_middle2[1] + $B2_length;
		$x_btw = ($x_from - $x_to)/2 + $x_to + 15;
		$y_btw = ($y_from - $y_to)/2 + $y_to;
		$x_from_shifted = $x_from + 3;
		$y_from_shifted = $y_from - 5;
		$x_to_shifted = $x_to + 3;
		$y_to_shifted = $y_to;
	}
	elsif ($coord_middle1[0]>$coord_middle2[0] and $coord_middle2[1]>$coord_middle1[1]){
		$x_from = $coord_middle1[0] - $C1_length;
		$y_from = $coord_middle1[1] + $B1_length;
		$x_to = $coord_middle2[0] + $C2_length;
		$y_to = $coord_middle2[1] - $B2_length;
		$x_btw = ($x_from - $x_to)/2 + $x_to - 15;
		$y_btw = ($y_to - $y_from)/2 + $y_from;
		$x_from_shifted = $x_from - 5;
		$y_from_shifted = $y_from + 2;
		$x_to_shifted = $x_to - 3;
		$y_to_shifted = $y_to - 3;
	}
	elsif ($coord_middle2[0]>$coord_middle1[0] and $coord_middle2[1]>$coord_middle1[1]){
		$x_from = $coord_middle1[0] + $C1_length;
		$y_from = $coord_middle1[1] + $B1_length;
		$x_to = $coord_middle2[0] - $C2_length;
		$y_to = $coord_middle2[1] - $B2_length;
		$x_btw = ($x_to - $x_from)/2 + $x_from - 15;
		$y_btw = ($y_to - $y_from)/2 + $y_from;
		$x_from_shifted = $x_from;
		$y_from_shifted = $y_from + 5;
		$x_to_shifted = $x_to - 3;
		$y_to_shifted = $y_to - 3;
	}
	elsif ($coord_middle2[0] > $coord_middle1[0] and $coord_middle2[1] eq $coord_middle1[1]){
		$x_from = $coord_middle1[0] + $radius1;
		$y_from = $coord_middle1[1];
		$x_to = $coord_middle2[0] - $radius2;
		$y_to = $coord_middle2[1];
		$x_btw = ($x_to - $x_from)/2 + $x_from;
		$y_btw = $y_to + 15;
		$x_from_shifted = $x_from;
		$y_from_shifted = $y_from + 5;
		$x_to_shifted = $x_to;
		$y_to_shifted = $y_to + 5;
	}
	elsif ($coord_middle2[0] eq $coord_middle1[0] and $coord_middle1[1] > $coord_middle2[1]){
		$x_from = $coord_middle1[0];
		$y_from = $coord_middle1[1] - $radius1;
		$x_to = $coord_middle2[0];
		$y_to = $coord_middle2[1] + $radius2;
		$x_btw = $x_to + 15;
		$y_btw = ($y_from - $y_to)/2 + $y_to;
		$x_from_shifted = $x_from + 5;
		$y_from_shifted = $y_from;
		$x_to_shifted = $x_to + 5;
		$y_to_shifted = $y_to;
	}
	elsif ($coord_middle1[0] > $coord_middle2[0] and $coord_middle2[1] eq $coord_middle1[1]){
		$x_from = $coord_middle1[0] - $radius1;
		$y_from = $coord_middle1[1];
		$x_to = $coord_middle2[0] + $radius2;
		$y_to = $coord_middle2[1];
		$x_btw = ($x_from - $x_to)/2 + $x_to;
		$y_btw = $y_to - 15;
		$x_from_shifted = $x_from;
		$y_from_shifted = $y_from - 5;
		$x_to_shifted = $x_to;
		$y_to_shifted = $y_to - 5;
	}
	elsif ($coord_middle2[0] eq $coord_middle1[0] and $coord_middle2[1] > $coord_middle1[1]){
		$x_from = $coord_middle1[0];
		$y_from = $coord_middle1[1] + $radius1;
		$x_to = $coord_middle2[0];
		$y_to = $coord_middle2[1] - $radius2;
		$x_btw = $x_to - 15;
		$y_btw = ($y_to - $y_from)/2 + $y_from;
		$x_from_shifted = $x_from - 5;
		$y_from_shifted = $y_from;
		$x_to_shifted = $x_to - 5;
		$y_to_shifted = $y_to ;
	}
	
	#Constuction of edge
	if ($type_of_edge){
		($x_from_shifted, $y_from_shifted, $x_btw, $y_btw, $x_to_shifted, $y_to_shifted) = correction_of_curled_arrows($x_to, $y_to, $radius1, $x_btw, $y_btw, $x_to_shifted, $y_to_shifted, $x_from_shifted, $y_from_shifted, \@coord_middle1, \@coord_middle2);
		my $edge = $canv->createLine($x_from_shifted, $y_from_shifted, $x_btw, $y_btw, $x_to_shifted, $y_to_shifted, -arrow=>'last', -smooth=>1, -fill=> 'snow4',
								-arrowshape => [@arrow_shape], -tags => [@tag1, @tag2, 'edge', 'curved']);	
	}elsif ($type_of_edge == 0){	
		my $edge = $canv->createLine($x_from, $y_from, $x_to, $y_to, -arrow=>'last', -fill=> 'snow4',
								-arrowshape => [@arrow_shape], -tags => [@tag1, @tag2, 'edge']);							
	}
	&main::set_arrow_width;	
	return $edge;
}

#-------------------------------------------#
# Coordinate correction of edge end points  #
#-------------------------------------------#	
sub correction_of_curled_arrows{
	my ($x_to, $y_to, $radius1, $x_btw, $y_btw, $x_to_shifted, $y_to_shifted, $x_from_shifted, $y_from_shifted, $array_coord_middle1, $array_coord_middle2) = @_;
	my @coord_middle1 = @{$array_coord_middle1};
	my @coord_middle2 = @{$array_coord_middle2};
	
	# Correction of curled arrow position depending on coordinate positioning
	if ([$y_to > ($coord_middle1[1] - $radius1/3*2) and $y_to < ($coord_middle1[1] - $radius1/3)] 
		or 
		[$y_to > ($coord_middle1[1] + $radius1/3) and $y_to < ($coord_middle1[1] + $radius1/3*2)]){
		if ($coord_middle2[0]>$coord_middle1[0] and $coord_middle1[1]>$coord_middle2[1]){
			$y_btw = $y_btw + 15;
			$x_to_shifted = $x_to_shifted - 2;
			$y_to_shifted = $y_to_shifted + 2;
		}
		elsif ($coord_middle1[0]>$coord_middle2[0] and $coord_middle1[1]>$coord_middle2[1]){
			$y_btw = $y_btw - 15;
			$x_from_shifted = $x_from_shifted - 2;
			$y_to_shifted = $y_to_shifted + 2;
		}
		elsif ($coord_middle1[0]>$coord_middle2[0] and $coord_middle2[1]>$coord_middle1[1]){
			$y_btw = $y_btw - 15;
			$x_to_shifted = $x_to_shifted + 2;
			$y_to_shifted = $y_to_shifted - 2;
		}
		elsif ($coord_middle2[0]>$coord_middle1[0] and $coord_middle2[1]>$coord_middle1[1]){
			$y_btw = $y_btw + 15;
			$x_from_shifted = $x_from_shifted + 2;
			$y_to_shifted = $y_to_shifted + 2;
		}
	}
	elsif ([$y_to > ($coord_middle1[1] - $radius1/3) and $y_to < ($coord_middle1[1] + $radius1/3)] 
		or $y_to eq ($coord_middle1[1] - $radius1/3)
		or $y_to eq ($coord_middle1[1] + $radius1/3)){
		if ($coord_middle2[0]>$coord_middle1[0] and $coord_middle1[1]>$coord_middle2[1]){
			$y_btw = $y_btw + 25;
			$x_to_shifted = $x_to_shifted - 2;
			$y_to_shifted = $y_to_shifted + 4;
		}
		elsif ($coord_middle1[0]>$coord_middle2[0] and $coord_middle1[1]>$coord_middle2[1]){
			$x_btw = $x_btw + 25;
			$x_from_shifted = $x_from_shifted - 2;
			$y_to_shifted = $y_to_shifted + 2;
		}
		elsif ($coord_middle1[0]>$coord_middle2[0] and $coord_middle2[1]>$coord_middle1[1]){
			$y_btw = $y_btw - 25;
			$x_to_shifted = $x_to_shifted + 2;
			$y_to_shifted = $y_to_shifted - 4;
		}
		elsif ($coord_middle2[0]>$coord_middle1[0] and $coord_middle2[1]>$coord_middle1[1]){
			$y_btw = $y_btw + 25;			
			$x_from_shifted = $x_from_shifted + 2;
			$y_to_shifted = $y_to_shifted + 2;
		}
	}
	return ($x_from_shifted, $y_from_shifted, $x_btw, $y_btw, $x_to_shifted, $y_to_shifted);
}

#-------------------------------------------#
#            Changing edge width            #
#-------------------------------------------#
sub apply_arrow_width{
	my ($canv, $arrow_width) = @_;
	my @edges = $canv->find(withtag => 'edge');
	foreach my $edge (@edges){
		$canv->itemconfigure($edge, -width =>$arrow_width);
	}
}

#-------------------------------------------#
#    Visualization of edge frequencies      #
#-------------------------------------------#
sub visualize_edge_frequencies{
	my ($canv,$hash_of_edges_frequencies) = @_;
	my %edges_frequencies = %{$hash_of_edges_frequencies};
	my $x_btw; my $y_btw;
	my @edges = $canv->find(withtag => "edge");
	foreach my $key (keys %edges_frequencies){
		foreach my $edge (@edges){
			my @tags = $canv->gettags($edge);
			$key =~ /(.+)(\s\-[\>\|]\s)(.+)/;
			if ($1 eq $3){
				if (($tags[1] eq 'edge') and ($tags[0] eq $1)){
					my @coords = $canv->coords($edge);
					$x_btw = ($coords[2]-$coords[4])/2+$coords[4]-18;
					$y_btw = ($coords[5]-$coords[3])/2+$coords[3];
					my $text = $canv->createText($x_btw, $y_btw, -text=>$edges_frequencies{$key}, -tags=>'label');
					$canv->lower($text, 'all');	
				}	
			}
			if ($key =~ /^$tags[0]\b(.+)$tags[1]$/){
				my @coords = $canv->coords($edge);
				if (defined($coords[5])){
					my $text = $canv->createText($coords[2], $coords[3], -text=>$edges_frequencies{$key}, -tags=>'label');
					$canv->lower($text, 'all');	
				}else{
					if ($coords[2]>$coords[0] and $coords[1]>$coords[3]){
						$x_btw = ($coords[2] - $coords[0])/2 + $coords[0] + 15;
						$y_btw = ($coords[1] - $coords[3])/2 + $coords[3];
					}
					elsif ($coords[0]>$coords[2] and $coords[1] > $coords[3]){
						$x_btw = ($coords[0] - $coords[2])/2 + $coords[2] + 15;
						$y_btw = ($coords[1] - $coords[3])/2 + $coords[3];
					}
					elsif ($coords[0]>$coords[2] and $coords[3] > $coords[1]){
						$x_btw = ($coords[0] - $coords[2])/2 + $coords[2] - 15;
						$y_btw = ($coords[3] - $coords[1])/2 + $coords[1];
					}
					elsif ($coords[2]>$coords[0] and $coords[3] > $coords[1]){
						$x_btw = ($coords[2] - $coords[0])/2 + $coords[0] - 15;
						$y_btw = ($coords[3] - $coords[1])/2 + $coords[1];
					}
					elsif ($coords[2] > $coords[0] and $coords[3] eq $coords[1]){
						$x_btw = ($coords[2] - $coords[0])/2 + $coords[0];
						$y_btw = $coords[3] + 15;
					}
					elsif ($coords[2] eq $coords[0] and $coords[1] > $coords[3]){
						$x_btw = $coords[2] + 15;
						$y_btw = ($coords[1] - $coords[3])/2 + $coords[3];
					}
					elsif ($coords[0] > $coords[2] and $coords[3] eq $coords[1]){
						$x_btw = ($coords[0] - $coords[2])/2 + $coords[2];
						$y_btw = $coords[3] - 15;
					}
					elsif ($coords[2] eq $coords[0] and $coords[3] > $coords[1]){
						$x_btw = $coords[2] - 15;
						$y_btw = ($coords[3] - $coords[1])/2 + $coords[1];
					}
					my $text = $canv->createText($x_btw, $y_btw, -text=>$edges_frequencies{$key}, -tags=>'label');
					$canv->lower($text, 'all');	
				}
				last;
			}
		}
	}
}	

#-------------------------------------------#
# Fit in canvas size newly created network  #
#-------------------------------------------#
sub adjust_scale{
	my ($canv, $fixed_canvas) = @_;
	my $new_size = 12;
	my @text_width = ();
	my $canv_height;
	my $canv_width;
	my $scale_factor;
	
	#Determine area, occupied by all canvas items, and rescale it up to visible canvas dimensions
	my ($x1, $y1, $x2, $y2) = $canv->bbox("all");
	my $contain_height = $y2-$y1+30;
	
	if (($canv->Height) eq 1){
		$canv_height = $canv->reqheight;
	}else{
		$canv_height = $canv->Height;
	}
	if ($contain_height > $canv_height or $contain_height < $canv_height){
		my $scale_factor = $canv_height/$contain_height;
		$canv->scale('all', 0, 0, $scale_factor, $scale_factor);
	}	
	($x1, $y1, $x2, $y2) = $canv->bbox("all");
	my $contain_width = $x2-$x1+30;
	if (($canv->Width) eq 1){
		$canv_width = $canv->reqwidth;
	}else{
		$canv_width = $canv->Width;
	}
	if ($contain_width > $canv_width){
		my $scale_factor = $canv_width/$contain_width;
		$canv->scale('all', 0, 0, $scale_factor, $scale_factor);
	}
	
	if ($canv->bbox("all")){
		$canv->configure(-scrollregion=>[$canv->bbox("all")]);
		$canv->xview('moveto', 0);
		$canv->yview('moveto', 0);
		$canv->delete('separate');
		($x1, $y1, $x2, $y2) = $canv->bbox("all");
		my $x_pos = ($canv_width-$x2+$x1-30)/2;
		my $rect = $canv->createRectangle(
							($x1-$x_pos),($y1-10),$x1,$y1,
							-outline =>'white',
							-fill => 'white',
							-tags=>'separate');
		$canv->lower($rect, 'all');	
		$canv->configure(-scrollregion=>[$canv->bbox("all")]);
		$canv->xview('moveto', 0);
		$canv->yview('moveto', 0);					
	}
	
	#Fitting the text in the node size, which was changed
	my @nodes = $canv->find(qw/withtag node/);
	my @members = $canv->itemcget($nodes[0], -members);
	my ($x1, $y1, $x2, $y2) = $canv->bbox($members[0]);
	my $oval_width = $x2 - $x1;
	foreach my $_ (@nodes){
		@members = $canv->itemcget($_, -members);
		$canv->itemconfigure($members[1], -font=>[-size=>$new_size]);
		my ($x1_, $y1_, $x2_, $y2_) = $canv->bbox($members[1]);
		my $text_width = $x2_ - $x1_ + 15;
		push @text_width, $text_width;
	}
	@text_width = sort {$b cmp $a} @text_width;
	if (defined($text_width[0]) and [$text_width[0] ne 0]){
		my $scale_factor = $oval_width/($text_width[0]);
		$new_size = $new_size*$scale_factor;
		$new_size = sprintf("%.1f", $new_size);
	}
	return $new_size;
}

#-------------------------------------------#	
#   Mode for moving all or selected nodes   #
#-------------------------------------------#
sub apply_mode{
	my ($canv, $type_of_mode, $node_color, $uniq_nodes, $ref_array_edges, $expression_predictions, $hash_of_edges_frequencies) = @_;
	if ($type_of_mode eq 'canvas'){
		lasso_selecting($canv, $node_color, $expression_predictions, $uniq_nodes);
		moving($canv, $hash_of_edges_frequencies);
	}
	if ($type_of_mode eq 'selected'){
		lasso_selecting($canv, $node_color, $expression_predictions, $uniq_nodes);
		group_moving($canv, $uniq_nodes, $ref_array_edges, $hash_of_edges_frequencies); 
	}
}

#-------------------------------------------#
#   Interactivity (moving of the elements)  #
#-------------------------------------------#
sub moving {
	my ($canv, $hash_of_edges_frequencies) = @_;
	my $oldx = 0; my $oldy = 0;
	my @coords;
	my $compare;
	my %chosen_nodes = ();
	my $group_from; my $group_to;
	my @interactive_edges = ();
	my @nodes = $canv->find(withtag=>'node');
	
	$canv->bind('node' => '<1>' =>
		sub{
			my ($x, $y) = ($Tk::event->x, $Tk::event->y);
			$oldx = $x; $oldy = $y;
			my $id = $canv->find(qw/withtag current/);

			$canv->bind($id => '<B1-Motion>' =>
				sub {
					my ($x1, $y1) = ($Tk::event->x, $Tk::event->y);
					$canv->move($id => $x1 - $oldx, $y1 - $oldy);
					$oldx = $x1; $oldy = $y1;
						
					my @members = $canv->itemcget($id, -members); 
					my @tags = $canv->gettags($members[1]);
					my @edges = $canv->find(withtag => "@tags");
						
					foreach my $_ (@edges){
						my @type = $canv->gettags($_);
						if ( grep/oval/,@type){
							$compare = $type[0];
						}
						if ( grep/edge/, @type){
							
							my @arrow_shape = $canv->itemcget($_, -arrowshape);
							if ($type[1] ne 'edge'){
								$type_of_edge = 0;
								if ($type[3]){
									$type_of_edge = 1;
								}
								if ($type[0] eq $compare){
									$group_from = $id;
									my @to = $canv->find(withtag => "$type[1]");
									foreach my $n(@to){
										my @type1 = $canv->gettags($n);
										if ( grep/node/,@type1){
											$group_to = $n;
										}
									}
								}else{
									$group_to = $id;
									my @from = $canv->find(withtag => "$type[0]");
									foreach my $m(@from){
										my @type2 = $canv->gettags($m);
										if ( grep/node/,@type2){
											$group_from = $m;
										}
									}
								}
							}else{
								my @arrow_shape = $canv->itemcget($_, -arrowshape);
								$canv->addtag('TEMP', withtag => "$_");
								self_loop($canv, $id, @arrow_shape);
								$canv->delete('TEMP');
							}
							$canv->addtag('TEMP', withtag => "$_");
							create_edge($canv, $type_of_edge, $group_from, $group_to, @arrow_shape);
							$canv->configure(-scrollregion => [$canv->bbox("all")]);
						}
						$canv->delete('TEMP');
						
					}
					if (defined($hash_of_edges_frequencies)){
						$canv->delete(withtag => "label");
						visualize_edge_frequencies($canv, $hash_of_edges_frequencies);
					}
					foreach my $_ (@nodes){
						$canv->raise($_, 'all');
					}
					$canv->raise($id, 'all');
				}
			);
	});
	$canv->bind('node' => '<B1-ButtonRelease>' => 
		sub {
			@interactive_edges = ();
		});
}

#-------------------------------------------#
#         Moving only selected nodes        #
#-------------------------------------------#
sub group_moving{
	my ($canv, $ref_uniq_nodes, $ref_array_edges) = @_;
	
	my $oldx = 0; 
	my $oldy = 0;
	my @selected = $canv->find(qw/withtag selected/);
	if (@selected){
		$canv->bind('node' => '<1>' =>
		sub{
			my ($x, $y) = ($Tk::event->x, $Tk::event->y);
			$oldx = $x; $oldy = $y;
			my @ID = $canv->find(qw/withtag current/);
			if (grep/@ID/,@selected){
				$canv->bind(@ID => '<B1-Motion>' =>
					sub {
						my ($x_coord, $y_coord) = ($Tk::event->x, $Tk::event->y);
						my $x_move = $x_coord - $oldx;
						my $y_move = $y_coord - $oldy;
						$oldx = $x_coord; $oldy = $y_coord;
						@selected = $canv->find(qw/withtag selected/);
						foreach my $selected_element (@selected){
							my @initial_coords = $canv->coords($selected_element);
							$canv->coords($selected_element, $initial_coords[0]+$x_move, $initial_coords[1]+$y_move);
							$canv->raise($selected_element, 'all');
						}
						$canv->configure(-scrollregion=>[$canv->bbox("all")]);
						$canv->addtag('TEMP', withtag => "edge");
						$canv->delete('TEMP');
						&call_for_edges($canv, $ref_uniq_nodes, $ref_array_edges);
					});	
		$canv->bind('node' => '<B1-ButtonRelease>' => 
			sub {
				@interactive_edges = ();
				my @nodes = $canv->find(withtag=>'node');
				foreach my $_ (@nodes){
					$canv->raise($_, 'all');
				}
			});			
			}else{moving($canv);}
		});	
	}	
	else{moving($canv);}		
}

#-------------------------------------------#
#         Zoom in and out subroutine        #
#-------------------------------------------#
sub zoom_scale{
	my ($amount, $direction, $canv, $size, $check_index) = @_;
	$amount = ($direction eq 'up') ? $amount : 1/$amount;
	$canv->scale('all',0,0,$amount,$amount);
	my $scaling_factor = $amount;
	if ($direction eq 'up'){
		if ($check_index < 1){
			$check_index++;
		}else{
			$size *=1.19;
			$size = sprintf("%.1f", $size);
		}
	}
	if ($direction eq 'down'){
		$size *=0.83333;
		$size = sprintf("%.1f", $size);
		if ($size > 1){
		}else{
			$size = 1;
			$check_index--;
		}
	}	
	my $width_outline = apply_node_outline($canv, $size);

	#Changing arrow width depending on zooming factor
	if ($size eq 1){$arrow_width = 0.5;}
	elsif ($size > 1 and $size <= 7){$arrow_width = 1;}
	elsif ($size > 7 and $size <= 9){$arrow_width = 1.5;}
	elsif ($size > 9){$arrow_width = $width_outline + 1;}
	
	#Changing the scrolling region
	$canv->configure(-scrollregion=>[$canv->bbox("all")]);
	$canv->xview('moveto', 0);
	$canv->yview('moveto', 0);
	return ($check_index, $size, $arrow_width);
};

#-------------------------------------------#
#     Zoom in and out by mouse wheel        #
#-------------------------------------------#
sub zooming_by_mouse_wheel{
	my ($mw, $canv, $scale_out, $scale_in) = @_;
	#Binding of mouse wheel, when it is in the canvas area	
	$canv->CanvasBind('<Enter>' => [$canv,'Tk::focus']);
	$canv->CanvasBind('<Leave>' => [$mw,'Tk::focus']);
	if ($^O eq 'MSWin32') {
        $canv->CanvasBind('<MouseWheel>' =>
			sub { 
				my $delta = $Tk::event->D;
				if ($delta eq '120'){
					$scale_in->invoke;
				}
				elsif ($delta eq '-120'){
					$scale_out->invoke;
				}
		});
	} else {
		$canv->CanvasBind('<4>' => sub {
            $scale_in->invoke;
        });
        $canv->CanvasBind('<5>' => sub {
            $scale_out->invoke;
        });
    }
}
#-----------------------------------------------------#
# Smoothing effect of zooming by outline modification #
#-----------------------------------------------------#
#Changing outline width of ovals depending on the zooming factor
sub apply_node_outline{
	my ($canv, $size) = @_;
	my $width_outline = 1;
	if ($size eq 1){$width_outline = 0;}
	elsif($size > 1 and $size <= 5){$width_outline = 0.5;}
	elsif($size > 5 and $size <= 18){$width_outline = 1;}
	elsif($size > 18 and $size <= 24){$width_outline = 1.5;}
	elsif($size > 24 and $size <= 29){$width_outline = 2;}
	elsif($size > 29 and $size <= 34){$width_outline = 2.5;}
	elsif($size > 34 and $size <= 40){$width_outline = 3;}
	elsif($size > 40 and $size <= 55){$width_outline = 3.5;}
	elsif($size > 55 and $size <= 69){$width_outline = 4;}
	elsif($size > 69 and $size <= 75){$width_outline = 5;}
	elsif($size > 75){
		my $times = ($size-75)/5-1; 
		$width_outline = 5 + $times;
	}
	my @ovals_width_outline = $canv->find(qw/withtag oval/);
	foreach my $oval (@ovals_width_outline){
		$canv->itemconfigure($oval, -width =>$width_outline);
	}
	return ($width_outline);
}	

#-------------------------------------------#
#               Highlighting                #
#-------------------------------------------#	
sub highlight{
	my $canv = shift;
	my @ids;
	my @members;
	my $initial_color;
	my $lighter_color;
	my $darker_color;
	my $mode_of_newly_selected = 0;

	# When the mouse is over the item, color is lighter
	$canv->bind('node', '<Enter>', 
		sub { 
			my @id = $canv->find(qw/withtag current/);
			@members = $canv->itemcget(@id, -members); 
			$initial_color = $canv->itemcget($members[0], -fill);
			($lighter_color, $darker_color) = parse_color($initial_color, 1);
			$canv->itemconfigure($members[0], -fill => $lighter_color);
		});
			
	#Item is selected by pressing <Control> button and clicking on the node
	$canv->bind('node' => '<Control-1>' =>
		sub{my ($x, $y) = ($Tk::event->x, $Tk::event->y);
			@ids = $canv->find(qw/withtag current/);
			my @members = $canv->itemcget(@ids, -members); 
			my @selected = $canv->find(qw/withtag selected/);
			if (grep /@ids/, @selected){ 
				$canv->dtag(@ids, "selected");
				my ($lighter_color1, $darker_color1) = &parse_color($initial_color, 2);
				$canv->itemconfigure($members[0], -fill => $darker_color1);
				$initial_color = $darker_color1;
			}
			elsif(!grep /@ids/, @selected){
				$canv->addtag("selected", "withtag", @ids);
				$canv->itemconfigure($members[0], -fill => $darker_color);	
				$mode_of_newly_selected = 1;
			}
			&main::mode;
		});

	# When the mouse leaves node - it restored its inital color.
	$canv->bind('node', '<Leave>',
		sub {my @id = $canv->find(qw/withtag current/);
			if ($mode_of_newly_selected){
				$canv->itemconfigure($members[0], -fill => $darker_color); 
			}else{
				$canv->itemconfigure($members[0], -fill => $initial_color); 
			}	
			$mode_of_newly_selected = 0;
		});		
}

#Due to preserving of information for frequency rate for each node, highlighting was excluded, but selection option is kept
sub selection_for_resulting_graph{
	my $canv = shift;
			
	#Item is selected by pressing <Control> button and clicking on the node
	$canv->bind('node' => '<Control-1>' =>
		sub{
			@ids = $canv->find(qw/withtag current/);
			my @members = $canv->itemcget(@ids, -members); 
			my $initial_color = $canv->itemcget($members[0], -fill);
			my @selected = $canv->find(qw/withtag selected/);
			if (grep /@ids/, @selected){ 
				$canv->dtag(@ids, "selected");
				my ($lighter_color1, $darker_color1) = &parse_color($initial_color, 2);
				$canv->itemconfigure($members[0], -fill => $darker_color1);
			}
			elsif(!grep /@ids/, @selected){
				$canv->addtag("selected", "withtag", @ids);
				my ($lighter_color, $darker_color) = &parse_color($initial_color, 1);
				$canv->itemconfigure($members[0], -fill => $darker_color);	
			}
		});
}

#-------------------------------------------#
#         Selection of nodes by lasso       #
#-------------------------------------------#
sub lasso_selecting{
	my ($canv, $current_node_color, $expression_predictions, $uniq_nodes) = @_;
	my $lasso;
	my @lasso_coords = ();	
	my $x;
	my $y;
	my @box;
	
	#If click on canvas with Control button - nothing deselects
	$canv->CanvasBind('<Control-1>' =>
		sub {@box = $canv->bbox("all");});
	
	#If click on canvas in empty area everything is deselected	
	$canv->CanvasBind('<1>' =>
		sub {	
			my @item = $canv->find(qw/withtag current/);
			if (@item){
				my @tags = $canv->gettags($item[0]);
				if (!grep/node/,@tags){	
					$canv->dtag("selected");
					if ($current_node_color ne 1){
						apply_node_color($canv, $current_node_color);
					}else{
						dermine_color_shade($canv, $expression_predictions, $uniq_nodes);
					}
				}
			}elsif(@item == 0){
				$canv->dtag("selected");
				if ($current_node_color ne 1){
					apply_node_color($canv, $current_node_color);
				}else{
					dermine_color_shade($canv, $expression_predictions, $uniq_nodes);
				}	
			}
			&main::mode;
		});		
	$canv->CanvasBind('<Control-Button1-Motion>' =>
		sub {
			my ($x_2, $y_2) = ($Tk::event->x, $Tk::event->y);
			if (@box){
				print "@box\n";
				$x = $x_2+$box[0];
				$y = $y_2+$box[1];
			}elsif(@box == 0){
				$x = $x_2;
				$y = $y_2;
			}
			my @item = $canv->find(qw/withtag current/);
			if (@item){
				my @tags = $canv->gettags($item[0]);
				if (!grep/node/,@tags){	
					$canv->delete(-tags => ['lasso']);
					$lasso = $canv -> createPolygon(0,0,
												-outline =>"green",
												-fill =>"green",
												-stipple => 'gray25',
												-tags => ['lasso']);
					push @lasso_coords, ($x, $y);
					$canv->coords($lasso,@lasso_coords);
				}
			}elsif(@item == 0){
				$canv->delete(-tags => ['lasso']);
				$lasso = $canv -> createPolygon(0,0,
												-outline =>"green",
												-fill =>"green",
												-stipple => 'gray25',
												-tags => ['lasso']);
				push @lasso_coords, ($x, $y);
				$canv->coords($lasso,@lasso_coords);
			}
		});		
	$canv->CanvasBind('<ButtonRelease-1>' =>
		sub {
			if ($lasso) {
				#Determine is the graph is homogeneous or heterogeneous colored
				my ($lighter_color, $darker_color) = parse_color($current_node_color, 1);
				
				# Find all nodes contained in the selected area
				my @lasso_area = $canv->bbox($lasso);
				if (@lasso_area){
					my @candidates = $canv->find('enclosed', @lasso_area);
					
					# Add to selection all enclosed nodes
					foreach my $candidate (@candidates) {
						my @tags = $canv->gettags($candidate);
						if (grep /node/, @tags) {
							$canv->addtag("selected", "withtag", $candidate);
							my @members = $canv->itemcget($candidate, -members);
							if ($current_node_color eq 1){
								my $initial_color = $canv->itemcget($members[0], -fill);
								my ($lighter_color, $darker_color) = parse_color($initial_color, 1);
								$canv->itemconfigure($members[0], -fill => $darker_color);
							}else{
								$canv->itemconfigure($members[0], -fill => $darker_color);
							}
						}
					}
				}	
			}	
			&main::mode;
			$canv->delete(-tags => ['lasso']);
			@lasso_coords = ();
    });
}

#------------------------------------------------------#
# Detection of color while highlighting and selection  #
#------------------------------------------------------#
sub parse_color{
	my ($initial_color, $selected) = @_;
	
	#Converting actual hexadecimal color to RGB format
	my @rgb = split(//, $initial_color);
    my $red = hex($rgb[1].$rgb[2]);
	my $green = hex($rgb[3].$rgb[4]);
	my $blue = hex($rgb[5].$rgb[6]);
	
	#Calculating RGB values for lighter color than actual
	my $red_light = (255 - $red) * 0.5 + $red;
	my $green_light = (255 - $green) * 0.5 + $green;
	my $blue_light = (255 - $blue) * 0.5 + $blue;
	my @rgb_light = (int($red_light), int($green_light), int($blue_light));
	
	#Calculating RGB values for darker color than actual
	my $red_dark = $red * 0.75;
	my $green_dark = $green * 0.75;
	my $blue_dark = $blue * 0.75;
	
	if ($selected eq 2){
		$red_dark = $red / 0.75;
		$green_dark = $green / 0.75;
		$blue_dark = $blue / 0.75;
	}
	my @rgb_dark = (int($red_dark), int($green_dark), int($blue_dark));
	
	#Converting RGB format colors back to hexadecimal
	my $hex_color_light = sprintf("#%02X%02X%02X", @rgb_light);
	my $hex_color_dark = sprintf("#%02X%02X%02X", @rgb_dark);

	return ($hex_color_light, $hex_color_dark);
}

#-------------------------------------------#
#             Apply node color              #
#-------------------------------------------#	
sub apply_node_color{
	my ($canv, $node_color) = @_;
	my @elements_array = $canv->find(qw/withtag oval/);
	my @selected = $canv->find(qw/withtag selected/);
	if (@selected){
	foreach my $existing_node (@elements_array){
			foreach my $selected_group (@selected){
				my @members = $canv->itemcget($selected_group, -members); 
				if ($members[0] eq $existing_node){
					my ($lighter_color, $darker_color) = parse_color($node_color, 1);
					$canv->itemconfigure($existing_node, -fill => $darker_color);
					last;
				} 	
				else{
					$canv->itemconfigure($existing_node, -fill => $node_color);
				}					
			}
		}	
	}elsif(@selected == 0){
		foreach my $existing_node (@elements_array){
			$canv->itemconfigure($existing_node, -fill => $node_color);
		}	
	}
	return $node_color;
}

#-------------------------------------------#
#         Setting color gradation           #
#-------------------------------------------#	
#Depending on the frequency of the expression prediction - there is gradation of color

sub dermine_color_shade{
	my ($canv, $expression_predictions, $uniq_nodes) = @_;
	my $color;
	my @expression_predictions = @{$expression_predictions};
	my %uniq_nodes = %{$uniq_nodes};
	foreach my $prediction (@expression_predictions){
		$prediction =~ /(.+)(: )(.+)/;
		if (grep/UNDEF/, $3){
			$color = '#f5f5f5';
		} else {
			my @phenotype = split "\t", $3;
			if (($phenotype[0] eq 'UP') or ($phenotype[0] eq 'ON') or ($phenotype[0] eq 'INVARIANT ON')){
				$color = '#ff0000';
				if ($phenotype[1] > 0.5 and $phenotype[1] <= 0.6){
					$color = '#ff9999';
				}elsif($phenotype[1] < 0.6 and $phenotype[1] <= 0.7){
					$color = '#ff7f7f';
				}elsif($phenotype[1] < 0.7 and $phenotype[1] <= 0.8){
					$color = '#ff4c4c';
				}elsif($phenotype[1] < 0.8 and $phenotype[1] <= 0.9){
					$color = '#ff3232';
				}elsif($phenotype[1] < 0.9 and $phenotype[1] < 1){
					$color = '#ff0000';
				}elsif($phenotype[1] == 1){
					$color = '#ff0000';
				}
			}
			if (($phenotype[0] eq 'OFF') or ($phenotype[0] eq 'DOWN') or ($phenotype[0] eq 'INVARIANT OFF')){
				$color = '#099400';
				if ($phenotype[1] > 0.5 and $phenotype[1] <= 0.6){
					$color = '#b7fcb2';
				}elsif($phenotype[1] < 0.6 and $phenotype[1] <= 0.7){
					$color = '#6ffa66';
				}elsif($phenotype[1] < 0.7 and $phenotype[1] <= 0.8){
					$color = '#0ff800';
				}elsif($phenotype[1] < 0.8 and $phenotype[1] <= 0.9){
					$color = '#0cc600';
				}elsif($phenotype[1] < 0.9 and $phenotype[1] < 1){
					$color = '#099400';
				}elsif($phenotype[1] == 1){
					$color = '#099400';
				}
			}	
		}
		my @member = $canv->itemcget($uniq_nodes{$1}, -members); 
		$canv->itemconfigure($member[0], -fill => $color);	
		my $color = '';
	}
}

#-------------------------------------------#
#        Changing font configurations       #
#-------------------------------------------#
sub configure_font {
	my ($canv, @font_configs) = @_;
	my @members = $canv->find(withtag=>'node');
	foreach my $member (@members){
		my @member = $canv->itemcget($member, -members); 
		$canv->itemconfigure($member[1], -font=>[@font_configs]);
	}
}

#-------------------------------------------#
#     Pop up information of frequencies     #
#-------------------------------------------#	
#Pop up balloon with frequency of the expression prediction, when the pointer is over the node

sub pop_up_info{
	my ($canv, $contain_canv, $expression_predictions) = @_;
	my @expression_predictions = @{$expression_predictions};
	my $value;
	$canv->bind('node' => '<Enter>' =>
		sub{
			my @id = $canv->find(qw/withtag current/);
			my @members = $canv->itemcget(@id, -members); 
			my @gene_name = $canv->itemcget($members[1], -tags);
			foreach my $prediction (@expression_predictions){
				if (grep /$gene_name[0]\b/,$prediction){
					$prediction =~ /(.+)(\t)(.+)/;
					$value = $3;
				}
			}
			my %message = ();
			$message{$id[0]} = $value;
			my $balloon = $contain_canv->Balloon(-bg => 'white');
			$balloon->attach($contain_canv->Subwidget("canvas"), 
				-initwait => 100, 
				-state => 'balloon', 
				-balloonposition => 'mouse', 
				-msg => \%message);
	});	
	$canv->bind(@id => '<B1-Motion>' =>
			sub{
				$canv->delete($balloon);
				return;
			});	
	$canv->bind(@id => '<1>' =>
			sub{
				$canv->delete($balloon);
				return;
			});	
}

1