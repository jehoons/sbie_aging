#!/usr/bin/perl

#### Ana Rodriguez Sanchez-Archidona, Luxembourg 27/7/2013, Luxembourg Center for Systems Biomedicine ####
#### Modul that includes different gene regulatory raw network and experimental gene expression data examples

use strict;
use warnings;

package  Examples;

#Example network and experimental data 1
sub example1{
	my @network_example = (
		'SNAI1 -> ZEB1',
		'SNAI1 -> ZEB2',
		'SNAI1 -| CDH1',
		'SNAI1 -| MIR-200',
		'SNAI1 -| MIR-203',
		'ZEB1 -| CDH1',
		'ZEB1 -| MIR-200',
		'ZEB1 -| MIR-203',
		'ZEB2 -| CDH1',
		'ZEB2 -| MIR-200',
		'ZEB2 -| MIR-203',
		'MIR-200 -| ZEB1',
		'MIR-200 -| ZEB2',
		'MIR-203 -| SNAI1',
		'MIR-203 -| ZEB2');
        
	my @phenotype = (
		"CDH1 DOWN\n",
		"MIR-200 DOWN\n",
		"MIR-203 DOWN\n",
		"SNAI1 UP\n",
		"ZEB1 UP\n",
		"ZEB2 UP\n");
	
    return (\@network_example,\@phenotype);
}

#Example network and experimental data 2
sub example2{
	my @network_example = (
		'FOXP3 -> FOXP3',
		'FOXP3 -| GATA3',
		'FOXP3 -| RORGT',
		'FOXP3 -| TBET',
		'GATA3 -> GATA3',
		'GATA3 -> IL10',
		'GATA3 -> IL4',
		'GATA3 -| RORGT',
		'GATA3 -| STAT4',
		'GATA3 -| TBET',
		'IFNB -> IFNBR',
		'IFNBR -> STAT1',
		'IFNG -> IFNGR',
		'IFNGR -> JAK1',
		'IL10 -> IL10R',
		'IL10R -> STAT3',
		'IL12 -> IL12R',
		'IL12R -> STAT4',
		'IL18 -> IL18R',
		'IL18R -> IRAK',
		'IL2 -> IL2R',
		'IL23 -> IL23R',
		'IL23R -> STAT3',
		'IL2R -> STAT5',
		'IL4 -> IL4R',
		'IL4R -> STAT6',
		'IL6 -> IL6R',
		'IL6R -> JAK3',
		'IRAK -> IFNG',
		'JAK1 -> STAT1',
		'JAK3 -> RORGT',
		'NFAT -> IFNG',
		'RORGT -> IL17',
		'RORGT -> IL6',
		'RORGT -> RORGT',
		'RORGT -| FOXP3',
		'RORGT -| GATA3',
		'RORGT -| TBET',
		'SOCS1 -| IL4R',
		'SOCS1 -| JAK1',
		'STAT1 -> SOCS1',
		'STAT1 -> TBET',
		'STAT1 -| IL4',
		'STAT3 -> RORGT',
		'STAT3 -| IFNG',
		'STAT4 -> IFNG',
		'STAT5 -> FOXP3',
		'STAT6 -> GATA3',
		'STAT6 -| IL12R',
		'STAT6 -| IL18R',
		'TBET -> IFNG',
		'TBET -> SOCS1',
		'TBET -> TBET',
		'TBET -| FOXP3',
		'TBET -| GATA3',
		'TBET -| RORGT',
		'TCR -> NFAT',
		'TGFB -> FOXP3',
		'TGFB -> TGFBR',
		'TGFBR -> RORGT');

	my @phenotype = (
		"FOXP3 DOWN\n",
		"GATA3 UP\n",
		"IFNB DOWN\n",
		"IFNBR DOWN\n",
		"IFNG DOWN\n",
		"IFNGR DOWN\n",
		"IL10 UP\n",
		"IL10R UP\n",
		"IL12 DOWN\n",
		"IL12R DOWN\n",
		"IL17 DOWN\n",
		"IL18 DOWN\n",
		"IL18R DOWN\n",
		"IL2 DOWN\n",
		"IL23 DOWN\n",
		"IL23R DOWN\n",
		"IL2R DOWN\n",
		"IL4 UP\n",
		"IL4R UP\n",
		"IL6 DOWN\n",
		"IL6R DOWN\n",
		"IRAK DOWN\n",
		"JAK1 DOWN\n",
		"JAK3 DOWN\n",
		"NFAT DOWN\n",
		"RORGT DOWN\n",
		"SOCS1 DOWN\n",
		"STAT1 DOWN\n",
		"STAT3 UP\n",
		"STAT4 DOWN\n",
		"STAT5 DOWN\n",
		"STAT6 UP\n",
		"TBET DOWN\n",
		"TCR DOWN\n",
		"TGFB DOWN\n",
		"TGFBR DOWN\n");
	
    return (\@network_example,\@phenotype);
}

#Example network and experimental data 3
sub example3{
	my @network_example = (
		'AR -> HSPB1',
		'AR -> JUN',
		'CHEK1 -| NFKB1',
		'CHEK1 -> RB1',
		'DNAdamage -> CHEK1',
		'HSPB1 -| IKBKB',
		'HSPB1 -> AR',
		'IGF1 -> IRS1',
		'IKBKB -| IRS1',
		'IKBKB -> NFKB1',
		'IKBKB -> RELA',
		'IRS1 -> RPS6KB1',
		'JUN -> IL6',
		'LowNutrition -| RPS6KB1',
		'NFKB1 -> IKBKB',
		'NFKB1 -> IL6',
		'NFKB1 -> JUN',
		'NFKB1 -> RELA',
		'PRKCQ -| IRS1',
		'PRKCQ -> IKBKB',
		'RB1 -> AR',
		'RB1 -> JUN',
		'RELA -| AR',
		'RELA -> JUN',
		'RELA -> NFKB1',
		'RPS6KB1 -| HSPB1',
		'RPS6KB1 -| IRS1',
		'IGF1 -> PRKCQ');

	my @phenotype = (
		"HSPB1 DOWN\n",
		"IGF1 UP\n",
		"IKBKB UP\n",
		"IL6 UP\n",
		"IRS1 UP\n",
		"LowNutrition DOWN\n",
		"NFKB1 UP\n",
		"PRKCQ UP\n",
		"RELA UP\n",
		"RPS6KB1 UP\n");
	
    return (\@network_example,\@phenotype);
}
1
