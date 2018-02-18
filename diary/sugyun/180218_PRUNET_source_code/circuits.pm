package circuits;
my %encoded_neighbors; 
my %encoded_neighbors_string; 
my %circuits_found=(); 
my %final_circuits=(); 
my $NETWORK={};
my %encoded_genes=(); 
my $stop_button_state = 'normal';
#Subrutines

sub mapping_adjacency_matrix_on_the_graph{
	my ($matrix_hash,$network_lines_array) = @_;
	
	my %matrix = %{$matrix_hash};
	my @network_lines = @{$network_lines_array};
	my @matrix_keys = keys %matrix;
	my %CIRCUITS = ();
	for (my $i=0;$i < @matrix_keys;$i++){
		my %temporal_interactions=();
		my @temporal_set_of_genes = split /:/, $matrix{$matrix_keys[$i]};
		for (my $j=0;$j<@temporal_set_of_genes;$j++){
			foreach my $network_line(@network_lines){
				if ($network_line=~ /(.+)(\s)(\-[\>\|])(\s)(.+)/){
					if (($1 eq $temporal_set_of_genes[$j]) and ($5 eq $temporal_set_of_genes[$j+1])){
						$temporal_interactions{$&}=$&;
					}
					else{}
				}
				else{}
			}
			
		}
		$CIRCUITS{"Circuit_$i"}=\%temporal_interactions;	
	}
	return(\%CIRCUITS);
}

sub circuits{
	my ($network_view_content_array,$circuits_stop_button,$f2) = @_;
	my @network_lines = @{$network_view_content_array};
	
	my %genes_in_circuits = %{main(\@network_lines,$circuits_stop_button,$f2)};
	my %CIRCUITS = ();
	undef @network_lines;
	undef%encoded_neighbors;
	undef%encoded_neighbors_string;
	undef%circuits_found;
	undef%final_circuits;
	undef$NETWORK;
	undef%encoded_genes;
	
	%CIRCUITS = %{mapping_adjacency_matrix_on_the_graph(\%genes_in_circuits,$network_view_content_array)};
	print values %{$CIRCUITS{Circuit_0}};
	SALIR:
	return (\%CIRCUITS);
}
sub main{
	my ($network_lines_array,$circuits_stop_button,$f2) = @_;
	my @network_lines = @{$network_lines_array};

	parser($network_lines_array);
	
	my @genesValues=values %$NETWORK; 
	my $numI=0; 
	foreach my $in(@genesValues)
	{
		my %hash=%$in;
		@val=keys %hash;
		$tot=scalar(@val);
		$numI+=$tot;
	}
	mainCircuits($circuits_stop_button,$f2);
	$refCircuitsFinal=deleteRep();
	%final_circuits=%$refCircuitsFinal;
	my $total=scalar(keys %final_circuits);
	my @borrar = values %final_circuits;
	return (\%final_circuits);
	
}
sub parser
	{
		my ($f)=@_;
		my @network_lines_local = @{$f};
		$cont=0;
		foreach my $line(@network_lines_local)
		{
			@words=split(' ',$line); 
			$value='-|'; 
			if ($words[1] eq '->') {$value='->'}; 
			$NETWORK->{$words[0]}->{$words[2]}=$value; 
			
			$key=$words[0];
			$neig=$words[2];
			if (exists $encoded_genes{$key}==0) 
			{
				$cont++;
				$encoded_genes{"$key"}="$cont";
			}
			if (exists $encoded_genes{$neig}==0) 
			{
				$cont++;
				$encoded_genes{"$neig"}="$cont";
			}
			if (exists $encoded_neighbors_string{$key}) 
			{
				@aux=@{$encoded_neighbors_string{$key}};
				push (@aux,$neig);
				$encoded_neighbors_string{"$key"}=["@aux"];
			}
			else
			{
				$encoded_neighbors_string{"$key"}=["$neig"];
			}
			if (exists $encoded_neighbors{$encoded_genes{$key}}) 
			{
				@aux=@{$encoded_neighbors{$encoded_genes{$key}}};
				push (@aux,$encoded_genes{$neig});
				$encoded_neighbors{"$encoded_genes{$key}"}=["@aux"];
			}
			else
			{
				$encoded_neighbors{"$encoded_genes{$key}"}=["$encoded_genes{$neig}"];
			}
		}

	}
sub mainCircuits{
	my ($circuits_stop_button,$f2) = @_;	
	my @stack=(); 
	my %blocked=();
	$start=2;
	my $s=0;
	my $nodes=scalar(keys %encoded_genes);
	my @path=values %encoded_genes; 
		my @Neig=@{$encoded_neighbors{$start}};
		$together=join('',@Neig);
		@Neig=split(' ',$together);
	while ($s<$nodes){
		$f2->update;#Very Important
		$stop_button_state = $circuits_stop_button->cget('-state');
		if ($stop_button_state eq 'disabled'){
			%final_circuits = ();
			goto SALIR;
		}
		else{}
			
		if (scalar @Neig >0)
		{
			$start=pop @path;
			foreach my $i(@path)
			{
				$blocked[$i]="F";
			}
		$found=findCycles($start,$circuits_stop_button,$f2);
		$s++;
		}
		else {$s=$nodes;}
	}
}
sub findCycles
	{	
		my($v,$circuits_stop_button,$f2)=@_;
		$f="F"; 
		my $selfLoop=0;
		push @stack,$v;
		$blocked[$v]="T";
		my @Neig=@{$encoded_neighbors{$v}};
		my $together=join('',@Neig);
		@Neig=split(' ',$together);
		
		for (my $i=0;$i<@Neig;$i++){
			$f2->update;#Very Important
			$stop_button_state = $circuits_stop_button->cget('-state');
			if ($stop_button_state eq 'disabled'){
				%final_circuits = ();
				goto SALIR;
			}
			else{}
			$w=$Neig[$i];
			my $keyStart=findKey($start,\%encoded_genes); 
			if ($w eq $start)
			{
				my $cycle='';
				my $cycleName='';
				my $key='';
				if(@stack==1){$selfLoop=1;};
				for (my $j=0;$j<@stack;$j++){
				$cycle = $cycle.":".$stack[$j];
				$key=findKey($stack[$j],\%encoded_genes); 
				$cycleName=$cycleName.":".$key;	
				}
			$f="T";
		
			$circuits_found{$cycle.":".$start}=$cycleName.":".$keyStart;
			$contCircuits++;
			$head=$w;
			}
			elsif ($blocked[$w] eq "F")
			{
				if (findCycles($w,$circuits_stop_button,$f2) eq "T"){$f="T";}
			}
		}
		if ($f eq "T" && $selfLoop==0 && $head!=$v){unblock($v);}
		else
		{
			foreach my $w(@Neig)
			{
				$enc=findElementIntoArray($v,$matrixB{$w}); 
				if ($enc==0) 
				{
					if (exists $matrixB{$w}) 
					{
						@aux=@{$matrixB{$w}};
						$enc2=findElementIntoArray($v,\@aux);
						if ($enc2==0)
						{
							push @aux,$v;
							$matrixB{$w}=["@aux"];
						}	
					}
					else {$matrixB{$w}=["$v"];}} 
				$i++;
			}		
		}
		delete $stack[scalar(@stack)-1]; 
		if ($f eq "F"){unblock($v)};
	return $f;		      
	}
sub findKey
	{
		my ($num,$hash)=@_;
		%hash=%$hash;
		@keys=keys %hash;
		@values=values %hash;
		$enc=0;
		$i=0;
		while ($i<@values && $enc==0)
		{
			if ($values[$i]==$num){$enc=1;}
			$i++;
		}
		$found=$keys[$i-1];
		return $found;
	}
sub findKeyString
	{
		my ($num,$hash)=@_;
		%hash=%$hash;
		@keys=keys %hash;
		@values=values %hash;
		$enc=0;
		$i=0;
		while ($i<@values && $enc==0)
		{
			if ($values[$i] eq $num){$enc=1;}
			$i++;
		}
		$found=$keys[$i-1];
		return $found;
	}
sub findElementIntoArray
	{
		my ($elem,$vector)=@_;
		my $i=0;
		$found2=0;
		@vect=@$vector;
		my $together=join('',@vect);
		@vect=split(' ',$together);
		while ($i<@vect && $found2==0)
		{
			if ($elem == $vect[$i]) {$found2=1;}
			$i++;
		}
		return $found2;
	}
sub findElementIntoArrayWithString
	{
		my ($elem,$vector)=@_;
		my $i=0;
		$found2=0;
		@vect=@$vector;
		while ($i<@vect && $found2==0)
		{
			if ($elem eq $vect[$i]) {$found2=1;}
			$i++;
		}
		return $found2;
	}
sub unblock
	{
		my ($u)=@_;
		$blocked[$u]="F";
		@aux=@{$matrixB{$u}};
		$together=join('',@aux);
		@aux=split(' ',$together);
		for (my $w=0;$w<@aux;$w++)
		{
			$elem=$aux[$w];
			delete $aux[$w];
		}
		
	}
sub deleteRep
	{
		my %circuitsHash=();
		 foreach $family (values %circuits_found)
		 {
			$fam="";
			@a=split(':',$family);
			$joinedOrig=join(':',@a); 
			delete $a[scalar(@a-1)];
			@a= sort @a;
			$joined=join(':',@a);	
			$circuitsHash{$joined}="$joinedOrig";	
		 }
	return \%circuitsHash;
	}
1