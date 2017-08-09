#=========================================
#   perl SSGA.pl -m matrix -p ratio -t times
#   p default is 0.6(0.1~1), which is the least ratio of the authentic connections to the theoretic connections
#   t defalut is 1, which is the times of Standard Vatiance  
#=========================================
use File::Basename;
use Getopt::Std;
use strict;

my (@gname,@NonzeroCon);
my (%pcon);
my $out=1;
#my $ii=0.6;
my %opt=();
getopts("m:p:t:",\%opt);
my $ff=1;
my $ii=$opt{p};
my $times=$opt{t};
open(GN,"<$opt{m}")or die "can not open the gene matrix";
while(<GN>){                     #get all gene_IDs and the non-zero connectivity
    s/\r\n?/\n/;
    chomp();
    my @array=split;
    push(@gname,$array[0]);
    my $fm=0;
    foreach my $item(@array){
        if($item>0 && $ff != $fm){
            push(@NonzeroCon,$item);
	}
        $fm++;
    }
    $ff++;
}
close(GN);

my ($avg,$std,$sum,$clen,$sqrsum);           #count average and standard_deviation of the non_zero connectivity
    foreach(@NonzeroCon){
	$sum += $_;
}
$clen=@NonzeroCon;
$avg = $sum/$clen;
foreach(@NonzeroCon){
    $sqrsum += ($_ - $avg)**2;
}
$std = ($sqrsum/@NonzeroCon)**0.5;
$sum=$avg+$times*$std;
$sum=int($sum);
open(MT,"<$opt{m}")or "can not open the matrix file\n";
while(<MT>){                          #get TF pairs and connectivty, record in hash %pcon
    s/\r\n?/\n/;
    chomp();
    my @array=split;
    for(my $i=1;$i<scalar(@array);$i++){
        my $str=$array[0]."\t".$gname[$i-1];
        if($array[$i]>=$sum){
            $pcon{$str}=$array[$i];
        }
    }
}
close(MT);

my @cpgname;
my %cppcon;
my $expbasename = basename($opt{m});
my $out_file="Cluster_".$expbasename;
open(OT,">>$out_file");
@cpgname=@gname;       # copy TF pairs and gene_IDs
%cppcon=%pcon;
&get_three();

sub get_three{                  #get three genes with maximal connectivity as the seed of subnetwork
    my $ratio=$ii;
    my @seed;
    foreach my $key(keys %cppcon){
        my @array=split(/\t/,$key);
        if($array[0] ne $array[1]){
            foreach my $item(@cpgname){
                if(($item ne $array[0])&&($item ne $array[1])){
                    my $edge1=$item."\t".$array[0];
                    my $edge2=$item."\t".$array[1];
                    if($cppcon{$edge1} && $cppcon{$edge2}){
                        if(!@seed){
                            push(@seed,$array[0]);
                            push(@seed,$array[1]);
                            push(@seed,$item);
                        }
                        else{
                            my ($SeedCon1,$SeedCon2,$SeedCon3,$SeedConSum);
                            $SeedCon1=$seed[0]."\t".$seed[1];
                            $SeedCon2=$seed[1]."\t".$seed[2];
                            $SeedCon3=$seed[2]."\t".$seed[0];
                            $SeedConSum=$cppcon{$SeedCon1}+$cppcon{$SeedCon2}+$cppcon{$SeedCon3};
                            if($SeedConSum<($cppcon{$key}+$cppcon{$edge1}+$cppcon{$edge2})){
				undef(@seed);
                                push(@seed,$array[0]);
                                push(@seed,$array[1]);
                                push(@seed,$item);
                            }
                        }
                    }
                }
            }
        }
    }
    if(@seed){
        foreach my $item(@seed){
            @cpgname=grep {$_ ne $item}@cpgname;
        }
        &get_cluster(\@seed,$ratio);
    }
}

sub get_cluster{            #grow seed to subnetwork
    my $flag=$_[0];
    my $ratio=$_[1];
    my @flag=@$flag;
    my @nodes;
    my $SizeFlag=@flag;
    my $point;
    my $value;
    foreach my $fitem1(@cpgname){                   #get the node with the maximum $value;
        my ($Degree1,$MulCon1,$SumCon1);
        foreach my $titem1(@flag){
            my $str1=$fitem1."\t".$titem1;
            if($cppcon{$str1}){
                $Degree1++;
                $MulCon1+=$cppcon{$str1}*$cppcon{$str1};
                $SumCon1+=$cppcon{$str1};
            }
        }
        if($Degree1>=$ratio*$SizeFlag && $SumCon1){
            my $Tabel1=$MulCon1/$SumCon1;
            if(!$point){
                $point=$fitem1;
                $value=$Tabel1;
	    }
            else{
                if($Tabel1>$value){
                    $point=$fitem1;
                    $value=$Tabel1;
                }
            }
        }
    }
    @cpgname=grep {$_ ne $point}@cpgname;
    push(@nodes,$point);
    foreach my $fitem2(@cpgname){                  #get all nodes with the maximum $value;
        my ($Degree2,$MulCon2,$SumCon2);
        foreach my $titem2(@flag){
            my $str2=$fitem2."\t".$titem2;
            if($cppcon{$str2}){
                $Degree2++;
                $MulCon2+=$cppcon{$str2}*$cppcon{$str2};
                $SumCon2+=$cppcon{$str2};
            }
        }
	if(($Degree2>=$ratio*$SizeFlag) && $value==($MulCon2/$SumCon2)){
            push(@nodes,$fitem2);
            @cpgname=grep {$_ ne $fitem2}@cpgname;
        }
    }
    if(!$nodes[0]){
	&get_subnetwork(\@flag);                   # no node can join into the seed
    }
    my $SizeNode=@nodes;
    if($nodes[0] && $SizeNode ==1){               # only one node with the maximum $value
	my @array=(@flag,@nodes);
        foreach my $item2(@array){
            @cpgname=grep {$_ ne $item2}@cpgname;
        }
	&get_cluster(\@array,$ratio);
    }
    if($nodes[0] && $SizeNode>1){                 # more than one node with the maximum $value
	my $Kvalue=1000;
        my @MergeNodes;
        foreach my $item4(@nodes){                # get the node with the minimum $Kvalue
            my @cpflag=@flag;
            push(@cpflag,$item4);
	    my @Degree;
            my $mul5=1;
	    foreach my $fitem5(@cpflag){
                my $Degree5;
                foreach my $titem5(@cpflag){
                    if($fitem5 ne $titem5){
                        my $str5=$fitem5."\t".$titem5;
                        if($pcon{$str5}){
                            $Degree5++;
                        }
                    }
                }
                push(@Degree,$Degree5);
		$mul5=$mul5*$Degree5;
	    }
	    my ($SD,$kvalue);
            my ($avg2,$sum2,$clen2,$sqrsum2);
            foreach(@Degree){
                $sum2 += $_;
            }
            $clen2=@Degree;
            $avg2 = $sum2/$clen2;
            foreach(@Degree){
                $sqrsum2 += ($_ - $avg2)**2;
            }
            $SD = ($sqrsum2/@Degree)**0.5;
	    $kvalue=$SD/$mul5;
	    if($Kvalue>$kvalue){
		$Kvalue=$kvalue;
	    }
	}
	foreach my $item6(@nodes){             # get all nodes with the minimum $Kvalue
	    my @cpflag=@flag;
	    push(@cpflag,$item6);
	    my @Degree1;
	    my $mul7=1;
	    foreach my $fitem7(@cpflag){
		my $Degree7;
		foreach my $titem7(@cpflag){
		    if($fitem7 ne $titem7){
			my $str7=$fitem7."\t".$titem7;
			if($pcon{$str7}){
			    $Degree7++;
			}
		    }
		}
		push(@Degree1,$Degree7);
		$mul7=$mul7*$Degree7;
	    }
	    my ($SD,$kvalue);
	    my ($avg3,$sum3,$clen3,$sqrsum3);
	    foreach(@Degree1){
		$sum3 += $_;
	    }
	    $clen3=@Degree1;
	    $avg3 = $sum3/$clen3;
	    foreach(@Degree1){
		$sqrsum3 += ($_ - $avg3)**2;
	    }
	    $SD = ($sqrsum3/@Degree1)**0.5;
	    $kvalue=$SD/$mul7;
	    if($kvalue==$Kvalue){
                push(@MergeNodes,$item6);
                @cpgname=grep {$_ ne $item6}@cpgname;
            }
	}
	my @array=(@flag,@MergeNodes);
	&get_cluster(\@array,$ratio);
    }
}


sub get_subnetwork{
    my $f_subnetwork=$_[0];
    my @subnetwork=@$f_subnetwork;
    my $w_size=@subnetwork;
    $w_size--;
    my $flag=0;

    foreach my $fitem(@subnetwork){                      # delete some node that the ratio is smaller than $ii 
	my $conn=0;
	foreach my $titem(@subnetwork){
	    if($titem ne $fitem){
		my $str=$titem."\t".$fitem;
		if($pcon{$str}){
		    $conn++;
		}
	    }
	}
	if($conn/$w_size <$ii){
	    @subnetwork=grep {$_ ne $fitem}@subnetwork;
	    push(@cpgname,$fitem);
	    $flag=1;
	    last;
	}
    }
    if($flag==1){
	&get_subnetwork(\@subnetwork);
    }

    else{                                               #output 
	foreach my $item(@subnetwork){
	    foreach my $key(keys %cppcon){
		if($key =~/$item/){
		    delete $cppcon{$key};
		}
	    }
	}
	print OT "Cluster_$out:\t@subnetwork\n";
	$out++;
	&get_three();
    }
}

