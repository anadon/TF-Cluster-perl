#=========================================
# perl MSGA.pl -m matrix -p percentage -q percentage -t times  
# p default is 0.6(0.1~1), which is the least ratio of the authentic connections to the theoretic connections
# q default is 0.75(0~1), which control the importance of the two part in the formula. 
# t defalut is 1, which is the times of Standard Vatiance
#=========================================
use File::Basename;
use Getopt::Std;
use strict;

my %opt=();
getopts("m:p:q:t:",\%opt);
open(GN,"<$opt{m}")or die "can not open genelist file\n";

my ($x,$y,$z); my (@gname,@num); my (%pcon,%scon,%con,%merge,%score,%incre);
my $vv=0.75; 
my $lmd=1;    
my $aaa=0.6;    
#my $vv=$opt{q};
#my $lmd=$opt{t};
#my $aaa=$opt{p};

my $ff=1;
while(<GN>){                       #get all gene_IDs and the non-zero connectivity
    s/\r\n?/\n/;
    chomp();
    my @arr=split;
    push(@gname,$arr[0]);
    my $fm=0;
    foreach my $item(@arr){
        if($item>0 && $ff != $fm){
            push(@num,$item);
        }
	$fm++;
    }
    $ff++;
}
close(GN);

my ($sum,$size,$avg,$std,$sqr);         #count average and standard_deviation of the non_zero connectivity
foreach(@num){
    $sum+=$_;
}
$size=@num;
$avg=$sum/$size;
foreach(@num){
    $sqr+=($_-$avg)**2;
}
$std=($sqr/$size)**0.5;
$sum=$avg+$lmd*$std;
$sum=int($sum);

my $expbasename = basename($opt{m});
my $out_file="Cluster_".$expbasename;
open(OT,">>$out_file");
#open(OT,">>Decomposit-result.txt");
open(MT,"<$opt{m}")or die "can not open the MT file\n";
while(<MT>){                                #get TF pairs and connectivty,record in hash %pcon,  and calculate node's degree saved in %scon 
    s/\r\n?/\n/;
    chomp;
    my @array=split;
    my $degree;
    for(my $i=1;$i<scalar(@array);$i++){
        my $str=$array[0]."\t".$gname[$i-1];
        if($array[$i]>=$sum){
            $pcon{$str}=$array[$i];
            $degree++;
        }
    }
    $scon{$array[0]}=$degree;
}
close(MT);
%con=%pcon;
&get_three();

sub get_three{                       #get three genes with maximal connectivity as the seed of a subnetwork
    my @fthree;
    foreach my $key(keys %pcon){
        my @arr=split(/\t/,$key);
        if($arr[0] ne $arr[1]){
	    foreach my $item(keys %scon){
		if(($item ne $arr[0]) && ($item ne $arr[1])){
		    my $str1=$item."\t".$arr[0];
		    my $str2=$item."\t".$arr[1];
		    if(($con{$str1})&&($con{$str2})){
			if(!@fthree){
			    push(@fthree,$arr[0]);
			    push(@fthree,$arr[1]);
			    push(@fthree,$item);
			}
			else{
			    my ($con1,$con2,$con3,$num);
			    $con1=$fthree[0]."\t".$fthree[1];
			    $con2=$fthree[1]."\t".$fthree[2];
			    $con3=$fthree[2]."\t".$fthree[0];
			    $num=$con{$con1}+$con{$con2}+$con{$con3};
			    if($num<($con{$str1}+$con{$str2}+$con{$key})){
				undef @fthree;
				push(@fthree,$arr[0]);
				push(@fthree,$arr[1]);
				push(@fthree,$item);
			    }
			}
		    }
		}
	    }
        }
    }
    foreach my $it(@fthree){
        delete $scon{$it};
    }
    if(@fthree){
        &get_seed(\@fthree);
    }
    else{
        &merge();
    }
}
my %seed;
my $order=1;
sub get_seed{                       #get a full connected seed 
    my $sd=$_[0];
    my @seed=@$sd;
    my $sumsd;
    my $Wvalue;
    my $mgene;
    my $sseed=@seed;
    foreach my $in(keys %scon){
        my $wv;
        my $uv;
        my $dv;
        my $number;
        foreach my $ce(@seed){
            my $str=$in."\t".$ce;
            if($con{$str}){
                $uv+=$con{$str}**2;
                $dv+=$con{$str};
                $number++;
            }
            else{
                last;
            }
        }
        if($sseed == $number){
            if($dv){
                $wv=$uv/$dv;
            }
            if(!$Wvalue){
                $Wvalue=$wv;
                $mgene=$in;
            }
            else{
                if($wv>$Wvalue){
                    $Wvalue=$wv;
                    $mgene=$in;
                }
            }
        }
    }
    if($mgene){
        push(@seed,$mgene);
	delete $scon{$mgene};
        &get_seed(\@seed);
    }
    else{
        foreach my $str(@seed){
            foreach my $kk(keys %pcon){
                if($kk =~ /$str/){
                    delete $pcon{$kk};
                }
            }
        }
        $seed{$order}=\@seed;
        $order++;
        &get_three();
    }
}

sub increase{                    # join one node into the seed with maximum score
    my $flag=0;
    my $score;
    my $skey;
    my $node;
    foreach my $item(keys %scon){
        foreach my $key(keys %seed){
            my $nub;
            my $wv;
            my $dv;
            my $ar=$seed{$key};
            my @arr=@$ar;
            foreach (@arr){
                my $str=$_."\t".$item;
                if($con{$str}){
                    $nub++;
                    $wv+=$con{$str}**2;
                    $dv+=$con{$str};
                }
            }
	    my $ss=@arr;
	    $ss=$aaa*$ss;
            if($nub>=$ss && $dv){
                $wv=100*(1-$vv)*$nub/@arr+$vv*$wv/$dv;
		if(!$score){
		    $score=$wv;
		    $skey=$key;
		    $node=$item;
		}
		else{
		    if($wv>$score){
			$score=$wv;
			$skey=$key;
			$node=$item;
		    }
		}
	    }
        }
    }

	if($skey){
            my $a=$seed{$skey};
            my @aa=@$a;
            push(@aa,$node);
            $seed{$skey}= \@aa;
            delete $scon{$node};
	    foreach (keys %pcon){
		if($_ =~ /$node/){
		    delete $pcon{$_};
		}
	    }
            $flag=1;
        }
    
    if($flag==0){
        &output();
    }
    else{
        &increase();
    }
}


sub merge{                       # calculate any two seeds' connection scroe, save in %merge
    foreach my $fkey(keys %seed){
        my $far=$seed{$fkey};
        my @farr=@$far;
	my $fsize=@farr;
        foreach my $tkey(keys %seed){
	    if($fkey ne $tkey){
            my $tar=$seed{$tkey};
            my @tarr=@$tar;
	    my $tsize=@tarr;
	    my $tot=$fsize*$tsize;
	    my $vedges=0;
	    my $up;
	    my $down;
	    foreach my $fitem(@farr){
		foreach my $titem(@tarr){
		    my $str=$fitem."\t".$titem;
		    if($con{$str}){
			$up+=$con{$str}**2;
			$down+=$con{$str};
			$vedges++;
		    }
		}
	    }
	    if($tot && $vedges/$tot>=$aaa){
		my $skey=$fkey."\t".$tkey;
		$merge{$skey}=(1-$vv)*100*$vedges/$tot+$vv*$up/$down;
	    }
	    }
	}
    }
    &reduce();
}
sub reduce{             #merge the two seed with maximum score, and delete all the other relations of the two seeds in the %merge 
    my $ta;
    my $fa;
    my $str;
    my $fmerge=0;
    foreach my $fkey(sort{$merge{$b}<=>$merge{$a}} keys %merge){
	if($merge{$fkey}){
	    my @tseeds=split("\t",$fkey);
	    $fa=$tseeds[0];
	    $ta=$tseeds[1];
	    my $far=$seed{$tseeds[0]};
	    my @farr=@$far;
	    my $tar=$seed{$tseeds[1]};
	    my @tarr=@$tar;
	    my @total=(@farr,@tarr);
	    if($fa>$ta){
                $str=$ta."-".$fa;
            }
            else{
                $str=$fa."-".$ta;
            }
	    $seed{$str}=\@total;
	    delete $seed{$ta};
	    delete $seed{$fa};
	    $fmerge=1;
	    foreach my $fkey(keys %merge){
		my @tseed=split("\t",$fkey);
		if($tseed[0] =~ /^$ta$/ ||$tseed[0] =~ /^$fa$/ ||$tseed[1] =~ /^$ta$/ ||$tseed[1] =~ /^$fa$/){
		    delete $merge{$fkey};
		}
	    }
	    last;
	}
    }
    if($fmerge ==1){
	&add($str);
    }
    else{
	&increase();
    }
}


sub add{               # add all the connection scores of the new merge seed into the %merge
    my $str=$_[0];
    my $far=$seed{$str};
    my @farr=@$far;
    my $fsize=@farr;
    foreach my $key(keys %seed){
	if($key ne $str){
	my $tar=$seed{$key};
	my @tarr=@$tar;
	my $tsize=@tarr;
	my $tot=$fsize*$tsize;
	my $vedges=0;
	my ($up,$down);
	foreach my $fitem(@farr){
	    foreach my $titem(@tarr){
		my $str=$fitem."\t".$titem;
		if($con{$str}){
		    $up+=$con{$str}**2;
		    $down+=$con{$str};
		    $vedges++;
		}
	    }
	}
	if($tot && $vedges/$tot>=0.6){
	    my $skey=$str."\t".$key;
	    $merge{$skey}=(1-$vv)*100*$vedges/$tot+$vv*$up/$down;
	}
	}
    }
    &reduce();
}

sub score{                # give all the subnetwork scores 
    foreach my $key(keys %seed){
        my $clu=$seed{$key};
        my @clu=@$clu;
        my $wu;
        my $wd;
        my $degree;
        my $value;
        my $size=@clu;
        my $ss=@clu;
        $size=$size*($size-1);
        foreach my $fitem(@clu){
            foreach my $titem(@clu){
                if($fitem ne $titem){
                    my $str=$fitem."\t".$titem;
                    if($con{$str}){
                        $wu+=$con{$str}**2;
                        $wd+=$con{$str};
                        $degree++;
                    }
                }
            }
        }
        if($size && $wd){
            $value=100*(1-$vv)*$degree/$size+$vv*$wu/$wd;
            $score{$key}=$value;
        }
    }
}

sub output{                  #depend on the score(from big to small) to output all the subnetworks 
    &score();
    my $num=1;
    foreach my $fkey(sort{$score{$b}<=>$score{$a}} keys %score){
        my $bb=$seed{$fkey};
        my @bb=@$bb;
	my $sbb=@bb;
	if($sbb>4){
	print OT "Cluster_$num:\t@bb\t$score{$fkey}\n";
	$num++;
	}
    }
}
