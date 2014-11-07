#!/usr/bin/perl

use warnings;
use strict;
use File::Slurp;
use LWP::Simple;
use List::Util qw(min);

my $all_uniprots=read_file('./uniprot/all_uniprot.txt');
my @all_uniprots=split("//\n",$all_uniprots);
my @proteins;
my @domains=("Protein kinase");
#foreach(@domains){
#    my $domain=$_;
foreach(@all_uniprots){
    my @line=split(/\.\n/,$_);
    my $protein_info='';
    my($uniprot_id, $protein_name, $short_name, @temporary_domains);
    my @pdbs;
    my $raw_sequence='';
    my($sequence_annotation,$domain_start,$domain_end,$domain_type,$sequence);
    foreach(@line){
	if(/\bAC\s+(\w{6})/){$uniprot_id=$1}
	if(/RecName:\sFull=(.+);/){$protein_name=$1}
	if(/GN\s+Name=([a-zA-Z0-9\-]+);/){$short_name=$1}
	if(/DR\s+PDB; ([a-zA-Z0-9\- ;\/\.=]+)/){push(@pdbs, $1)}
	if(/FT\s+(DOMAIN)\s+(\d+)\s+(\d+)\s+($domains[0])/){
	    $sequence_annotation=$1;
	    $domain_start=$2;
	    $domain_end=$3;
	    my $domain=$4;
	    if($domain=~/([a-zA-Z0-9\- ]+)[0-9]?/){$domain_type=$1}
	}
	my @tail=split(/\n/,$line[-1]);
	foreach(@tail){
	    if(/^\s/){
		#need to remove spaces from the raw sequence lines
		my @sequence_fragment=split(/\s+/,$_); 
		foreach(@sequence_fragment){$raw_sequence.=$_}
	    }
	}
	if($domain_start && $domain_end && $domain_type){
	    my $kinase_sequence=substr $raw_sequence,$domain_start-1,$domain_end-$domain_start+1;
	    $protein_info=$domain_type."; ".$uniprot_id."; ".$protein_name.
		"; ".$short_name."; ".$domain_start."; ".$domain_end.";".
		$kinase_sequence."|";
	    foreach(@pdbs){
		if(/=([0-9]+)\-([0-9]+)/){
		    if($domain_start>$2 || $domain_end<$1){
			next; #the crystal does not contain the domain
		    }
		    $protein_info .= "$_|"}    
	    }
	}
    }
    if($protein_info){
	push(@proteins,$protein_info);
    }
}

#foreach(@proteins){print "$_\n"}

#make_protein_list(\@proteins);
sub make_protein_list{
    my @proteins=@{$_[0]};
    foreach(@proteins){
	my @structures=split(/\|/,$_);
	my @header=split(/;/,shift(@structures)); #get everything before the first |
	my $uniprotID=$header[1];
	my $protein=$header[3];
	my $domain_start=$header[4];
	my $domain_end=$header[5];
	my %missing;
	foreach(@structures){
	    my $pdb_id=substr $_, 0, 4;
	    my $pdb_text=get("http://www.pdb.org/pdb/files/$pdb_id.pdb") or die $!;
	    my @pdb_lines=split(/\n/,$pdb_text);
	    my %chains;
	    foreach(@pdb_lines){
		if(/^DBREF.+$uniprotID/){ #only take chains that contain our protein	
		    my $chainID = substr $_, 12, 1;
		    my $crystal_start = substr $_, 14, 4;
		    my $crystal_end = substr $_, 20, 4;
		    $crystal_start=~s/^\s+|\s+$//g;
		    $crystal_end=~s/^\s+|\s+$//g;
		    $chains{$chainID}="$crystal_start $crystal_end";
		    #print "$pdb_id $chainID $chains{$chainID}\n";
		}
	    }
	    while( my ($chain_id,$value) = each %chains){
		my @seq_info=split(/ /,$value);
		my $crystal_start=$seq_info[0];
		my $crystal_end=$seq_info[1];
		my $nterm_missing_res = (($crystal_start-$domain_start)>0)? 
		    ($crystal_start-$domain_start) : 0;
		my $cterm_missing_res = (($domain_end-$crystal_end)>0)? 
		    ($domain_end-$crystal_end) : 0;
		my $missing_ends=$nterm_missing_res+$cterm_missing_res;
		#print "$pdb_id $chain_id $crystal_start $domain_start $crystal_end $domain_end $missing_ends\n";
		my ($res_idx,%seq_conflict);
		if($pdb_text!~/REMARK 465/){
		    my $rmsd=rmsd("$pdb_id-$chain_id");
		    if($rmsd eq -1){
			print "$protein $pdb_id $chain_id $missing_ends not a kinase domain\n";
			next;
		    }
		    if($rmsd eq 0){
			print "$protein $pdb_id $chain_id $missing_ends align failed\n";
			next;
		    }
		    print "$protein $pdb_id $chain_id $missing_ends $rmsd\n";
		    next; #no need to look for missing residues
		}
		foreach(@pdb_lines){
		    my ($mutant_resname,$chainID,$wt_resname,$mutant_idx,$info);
		    if(/SEQADV.+ENGINEERED MUTATION/){
			$mutant_resname = substr $_, 12, 3;
			$chainID = substr $_, 16, 1;
			$wt_resname = substr $_, 39, 3;
			$mutant_idx = substr $_, 19, 5;
			$mutant_idx=~s/^\s+|\s+$//g; #trim spaces
			if($chainID && $mutant_idx>=$domain_start && $mutant_idx<=$domain_end){
			    $seq_conflict{$chainID}.="$wt_resname-$mutant_idx-$mutant_resname "}
		    }
		    if(/^MODRES/){
			$mutant_resname = substr $_, 12, 3;
			$chainID = substr $_, 16, 1;
			$wt_resname = substr $_, 24, 3;
			$mutant_idx = substr $_, 18, 5;
			$mutant_idx=~s/^\s+|\s+$//g; #trim spaces
			if($chainID && $mutant_idx>=$domain_start && $mutant_idx<=$domain_end){
			    $seq_conflict{$chainID}.="$wt_resname-$mutant_idx-$mutant_resname "}
		    }
		}
		foreach(@pdb_lines){
		    if(/REMARK 465\s+[A-Z]{3}\s+$chain_id\s+([0-9]+)/){
			$res_idx=$1;
			if($res_idx>=$domain_start && $res_idx<=$domain_end){
			    #print "$pdb_id $res_idx $chain_id\n";
			    if($seq_conflict{$chain_id}){
				if(exists $missing{"$protein $pdb_id $chain_id $seq_conflict{$chain_id}"}){
				    $missing{"$protein $pdb_id $chain_id $seq_conflict{$chain_id}"}+=1;
				} else {$missing{"$protein $pdb_id $chain_id $seq_conflict{$chain_id}"}=1+$missing_ends}
			    }
			    if(!%seq_conflict){
				if(exists $missing{"$protein $pdb_id $chain_id"}){
				    $missing{"$protein $pdb_id $chain_id"}+=1;
				} else {$missing{"$protein $pdb_id $chain_id"}=1+$missing_ends}
			    }
			}
			if($res_idx<$domain_start || $res_idx>$domain_end){
			    #print "$pdb_id $res_idx $chain_id\n";
			    if($seq_conflict{$chain_id}){
				if(exists $missing{"$protein $pdb_id $chain_id $seq_conflict{$chain_id}"}){
				    $missing{"$protein $pdb_id $chain_id $seq_conflict{$chain_id}"}+=0;
				} else {$missing{"$protein $pdb_id $chain_id $seq_conflict{$chain_id}"}=$missing_ends}
			    }
			    if(!%seq_conflict){
				if(exists $missing{"$protein $pdb_id $chain_id"}){
				    $missing{"$protein $pdb_id $chain_id"}+=0;
				} else {$missing{"$protein $pdb_id $chain_id"}=$missing_ends}
			    }
			}
		    }
		}
	    }
	}
	while(my ($key,$val) = each %missing){
	    my @key_info=split(/ /,$key);
	    my $pdb_id=$key_info[2];
	    my $chain_id=$key_info[3];
	    my $rmsd=rmsd("$pdb_id-$chain_id");
	    if($rmsd eq -1){
		print "$key $val not a kinase domain\n";
		next;
	    }
	    if($rmsd eq 0){
		print "$key $val align failed\n";
		next;
	    }
	    print "$key $val $rmsd\n";
	}
    }
}

get_protein_chains(\@proteins);
sub get_protein_chains{
    my @proteins=@{$_[0]};
    foreach(@proteins){
	my @structures=split(/\|/,$_);
	my @header=split(/;/,shift(@structures));
	my $uniprotID=$header[1];
	my $protein=$header[3];
	my $domain_start=$header[4];
	my $domain_end=$header[5];
	my $sequence=$header[6];
	my @sequence=split(/(.{67})/,$sequence);
	my %missing;
	{local $|=1;
	 printf "\r%-10s",$protein}
	foreach(@structures){
	    my $pdb_id=substr $_, 0, 4;
	    my $pdb_text=get("http://www.pdb.org/pdb/files/$pdb_id.pdb") or die $!;
	    my @pdb_lines=split(/\n/,$pdb_text);
	    my @chain_list;
	    my %chains;
	    foreach(@pdb_lines){
		if(/^DBREF.+$uniprotID/){ #only take chains that contain our protein	
		    my $chainID = substr $_, 12, 1;
		    my $crystal_start = substr $_, 14, 4;
		    my $crystal_end = substr $_, 20, 4;
		    $crystal_start=~s/^\s+|\s+$//g;
		    $crystal_end=~s/^\s+|\s+$//g;
		    $chains{$chainID}="$crystal_start $crystal_end";
		    #print "$pdb_id $chainID $chains{$chainID}\n";
		}
	    }
	    while( my ($chain_id,$value) = each %chains){
		my @seq_info=split(/ /,$value);
		my $crystal_start=$seq_info[0];
		my $crystal_end=$seq_info[1];
		my $filename="./PDBs/".$pdb_id."-".$chain_id.".pdb";
		open PDBFILE, '>', $filename or die $!;
		print PDBFILE "HEADER    : $pdb_id-$chain_id: $domain_start: $chain_id: $domain_end: $chain_id :$protein : : :\n";
		foreach(@sequence){
		    if(!/[A-Z]/){next}
		    print PDBFILE "REMARK 300 $_\n";
		}
		foreach(@pdb_lines){
		    if(/ATOM [A-Z0-9 ]{16}[$chain_id]\s{0,4}([0-9]{1,5})/){
			my $res_idx=$1;
			if($res_idx>=$domain_start && $res_idx<=$domain_end){
			    print PDBFILE "$_\n";
			}
		    }
		    if(/HETATM[A-Z0-9 ]{15}[$chain_id]\s{0,4}([0-9]{1,5})/){
			my $res_idx=$1;
			if($res_idx>=$domain_start && $res_idx<=$domain_end){
			    print PDBFILE "$_\n";
			}
		    }
		}
		close PDBFILE;
	    }
	}
    }
    print "\n";
}

sub rmsd{
    my $structure=$_[0];
    my $structure_file="./PDBs/$structure.pdb";
    my $active_template='active-egfr-monomer.pdb';
    my $inactive_template='egfr-inactive-monomer.pdb';
    my @actives=qw/2HZI-A 4C3P-A 3NR9-B 4PF4-A 4HNF-B 4C3P-A 2ITW-A 3IW4-A 1LFR-A 3C0H-B 3E7E-A 3KUL-B3BHH-B 2CMW-A 4IDV-A/; 
    my @inactives=qw/2G1T-A 3TV4-B 1LG3-A 4K9Y-A 1QCF-A 3W32-A 3KMW-A 4C57-B 3G33-A 4FGB-A 3TT0-A 4OBO-A/;
    my (@rmsd_act,@rmsd_inact);
    foreach(@actives){
	my $active="./PDBs/$_.pdb";
	my $command= <<EOF;
load $active, active;
load $structure_file, template;
align active, template, cycles=0
EOF
	my $rmsd_active =`pymol -cd "$command"`;
	if($rmsd_active=~/Executive: RMS =\s+([0-9]{1,2}\.[0-9]{3})/){
	    push(@rmsd_act,$1);
	}
    }
    foreach(@inactives){
	my $inactive="./PDBs/$_.pdb";
	my $command= <<EOF;
load $inactive, inactive;
load $structure_file, template;
align inactive, template, cycles=0
EOF
	my $rmsd_inactive =`pymol -cd "$command"`;
	if($rmsd_inactive=~/Executive: RMS =\s+([0-9]{1,2}\.[0-9]{3})/){
	    push(@rmsd_inact,$1);
	}
    }
    my $rmsd_act=min @rmsd_act;
    my $rmsd_inact=min @rmsd_inact;
    if($rmsd_act && $rmsd_inact){
	return "rmsd-act: $rmsd_act, rmsd-inact: $rmsd_inact";
    } else {
	my $died=0;
	open PDBFILE, $structure_file or $died=-1;
	return $died;
    }
}

