package Translator;
	my $num_nodes;
	my $fcount = 0;
	my $n = 1;
	my $found = 0;
	my $errString = "";
sub translator{
	my ($array_lines,$number_of_nodes) = @_;
	my @lines = @{$array_lines};
	my @polynomial_functions=();
	$num_nodes = $number_of_nodes;
	$fcount = 0;
	$n = 1;
	$found = 0;
	$errString = "";
foreach my $line(@lines){
    #$line = $_;
	$_=$line;
    #remove newline character
    #chomp($line);
    #remove all spaces
    $line =~ s/\s*//g;
	 $line =~ s/\t*//g;
    #remove repetitions of equals
    $line =~ s/=+/=/g;
    my $count = 0;
    #$count = $line =~ s/\*{2}/\^/g;
    #if the line begins with f or F
    if($line =~ m/^f$n=/i)
    {
		$fcount++;
       my $func = (split(/=/,$line))[1]; # just read the function
       if(length($func) == 0)
       {
         $errString =  "ERROR: Empty function no $n\n";
         $found = 1;
         last;
       }
	   
       if($func =~ m/[^(x)(\d)(\()(\))(\+)(\*)(\~)]/g)
       {print "problematic function $line\n";
         $errString =  "ERROR: Found unacceptable character(s) in function $n\n";
         $found = 1;
         last;
       }
       # check to see if there are equal number of opening and closing paranthesis
       if( tr/\(/\(/ != tr/\)/\)/ )
       {
         $errString = "ERROR: Missing paranthesis in function $n.\n";
         $found = 1;
         last;
       }
       #check to see if the index of x is acceptable
       my $err = 0;
       while( $func =~ m/x(\d+)/g )
       {
         if( ($1 > $num_nodes) || ($1 < 1) )
         {
           $errString = "ERROR: Index of x out of range in function $n.\n ";
           $found = 1;
           $err = 1;
           last;
         }
       }
       #check to see if there was any error in the above while loop
       if($err == 1)
       {
		 $found = 1;
         last;
       }
       #Check to see if function is starting properly
       if($func =~ m/^[\)\*\+\~]/)
       {
         $errString = "ERROR: Incorrect syntax in function $n. Inappropriate char at start of function\n";
         $found = 1;
         last;
       }
       #Check to see if function is ending properly
       if($func =~ m/[^\)\d]$/)
       {
         $errString = "ERROR: Incorrect syntax in function $n. Inappropriate char at end of function\n";
         $found = 1;
         last;
       }
       #check to see if x always has an index
       if($func =~ m/x\D/g)
       {
         $errString = "ERROR: Incorrect syntax of function $n. Check x variable\n";
         $found = 1;
         last;
       }
	   if( ($func =~ m/[\+\*\(][\)\+\*]/g) || ($func =~ m/[\+\*][\~]/g) || ($func =~ m/\)[\(\d\~ x]/g) || ($func =~ m/\d[\(\~ x]/g) )
       {
         $errString = "ERROR: Incorrect syntax of function $n.\n";
         $found = 1;
         last;
       }
       if($found == 0)
       {
		  my @express = ();
		  my @testarr = split( /([\(\)\*\+\~ ])/,$func);
          for(my $i = 0; $i < scalar(@testarr); $i++)
		  {
             if($testarr[$i] ne "" && $testarr[$i] ne " " )
			 {
                push(@express,$testarr[$i]);
             }
           }
          push @polynomial_functions, "f$fcount = ".evaluateInfix(@express)."\n";
       }
    }
    else
    {
      $errString = "ERROR: Incorrect start of function declaration in function $n.\n";
      $found = 1;
      last;
    }
    $n++;
}

if($found == 1)
{
  print "Errors found in the input file. See below for description:\n";
  print "-----------------------------------------------------------------------\n";
  print $errString."\n";
  die("Errors with input file..ending program");
}
if(($n-1) < $num_nodes)
{
	print "Errors found in the input file. See below for description:\n";
    print "-----------------------------------------------------------------------\n";
	print "ERROR: Insufficient number of functions in the input file. Check your number of nodes field\n";
    die("Errors with input file..ending program");
}
return (\@polynomial_functions);
}
sub isOperand
{	
    my ($param) = shift;
    if((!isOperator($param)) && ($param ne "(") && ($param ne ")"))
    {
        return 1;
    }
    else
    {
        return;
    }
}

sub isOperator
{
    my ($param) = shift;
    if(($param eq "+") || ($param eq "*") || ($param eq "~"))
    {
        return 1;
    }
    else
    {
        return;
    }
}


sub InfixToPostfix
{
    my @infixStr = @_;
    my @postfixArr = ();
    my @stack = ();
    for(my $i=0; $i<scalar(@infixStr); $i++)
    {
        if(isOperand($infixStr[$i]))  #if its an operand ( * or + or ~)
        {
            push(@postfixArr,$infixStr[$i]); # add it to the end of the postfixarray or as say top of stack
        }
        if(isOperator($infixStr[$i])) #if its an operator (x1, x2, x3 etc)
        {
            push(@stack,$infixStr[$i]);   # add it to the top of stack
        }
        if($infixStr[$i] eq ")") #this marks the latest closed paranthesis
        {
            push(@postfixArr, pop(@stack)); #remove from stack and push it in postFixArray
        }
    }
    return @postfixArr;
}
sub PostfixEval
{
    my @postfixArr = @_;
    my @stackArr = ();
    for(my $i=0; $i<scalar(@postfixArr); $i++)
    {
        if(isOperand($postfixArr[$i]))
        {
            push(@stackArr,$postfixArr[$i]);
        }
        if(isOperator($postfixArr[$i]))
        {
          my $val = pop(@stackArr);
          if($postfixArr[$i] eq "~")
          {
              $val = "(".$val."+1)";
              push(@stackArr, $val);
          }
          if($postfixArr[$i] eq "+")
          {
             my $val2 = pop(@stackArr);
             $val = "((".$val."+".$val2.")+(".$val."*".$val2."))";
             push(@stackArr, $val);
          }
          if($postfixArr[$i] eq "*")
          {
             my $val2 = pop(@stackArr);
             $val = "(".$val."*".$val2.")";
             push(@stackArr, $val);
          }
        }
    }
    return @stackArr;
}

sub joinArr
{
    my (@who)=@_;
    my $who_len=@who;
    my $retVal = "";
    for(my $i=0; $i<$who_len; $i++)
    {
        $retVal.=$who[$i];
    }
    return $retVal;
}
sub evaluateInfix
{
    my @exp = @_;
    return (joinArr(PostfixEval(InfixToPostfix(@exp))));
}
1;