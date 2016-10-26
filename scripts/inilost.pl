if( ($ARGV[0] eq "") || ($ARGV[1] eq "") || ($ARGV[2] eq "") ){
    print "usage: perl inilost.pl distr lost out\n ";
    exit(0);
}
open(fpr1, $ARGV[0]) || die "Cannot open $ARGV[0] $!\n";
open(fpr2, $ARGV[1]) || die "Cannot open $ARGV[1] $!\n";
open(fpw, ">".$ARGV[2]) || die "Cannot open $ARGV[2] $!\n";

$n=0;
while( <fpr1> ){
    @buf=split ;
    if( ($buf[0] !~ /Phys/) && ($buf[0] !~ /Norm/) ){
        $x0[$n]    =$buf[0]; 
        $px0[$n]   =$buf[1];
        $y0[$n]    =$buf[2];
        $py0[$n]   =$buf[3];
        $z0[$n]    =$buf[4];
        $p0[$n]    =$buf[5];
        $wt[$n]    =$buf[6];
        $n0[$n]    =$buf[7];
        $n=$n+1;
    }
}
$n1=$n;
printf "Number of particles read from distr: %d \n", $n1;
close(fpr1);
#
$n=0;
while( <fpr2> ){
    @buf=split ;
        $xl[$n]    =$buf[0]; 
        $pxl[$n]   =$buf[1];
        $yl[$n]    =$buf[2];
        $pyl[$n]   =$buf[3];
        $zl[$n]    =$buf[4];
        $pl[$n]    =$buf[5];
        $wh[$n]    =$buf[6];
        $whr[$n]   =$buf[7];
        $tl[$n]    =$buf[8];
        $nl[$n]    =$buf[9];
        $n=$n+1;
}
$n2=$n;
printf "Number of particles read from lost: %d \n", $n2;
close(fpr2);
#
printf fpw "# x px y py z pz weight turn why\n";
for($i=0;$i<$n2;$i++){
    $j=0;
    do{
	if( $nl[$i] == $n0[$j] ){
	    printf fpw "%f %f %f %f %f %f %d %d %s\n",
	    $x0[$j],$px0[$j],
            $y0[$j],$py0[$j],
            $z0[$j],$p0[$j],
	    $wt[$j],$tl[$i],$wh[$i];
	}
	$j=$j+1;
    }while($j<$n1);
}
close(fpw);
