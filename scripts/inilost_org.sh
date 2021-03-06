if [ $# == 0 ]; then
echo "This script extracts initial coordinates of lost particles in lifetrac tracking"
echo "by A.Valishev (x2875)"
echo " "
echo "usage: $0 LTR_File"
echo "output: inilost.1.out - initial coordinates of particles lost on the first step"
echo "        inilost.out - initial coordinates of all lost particles in subsequent steps"
echo "initial distribution file can be specified via command DISTR=file_path $0 LTR_File"
echo "initial distribution must be in 'Norm' coordinates"
exit 0;
fi


DISTRDEFAULT='/home/valishev/lifetrac/src/distr_shatilov.txt'
DISTR=${DISTR:-$DISTRDEFAULT}

n0=`grep '(died)' $1 |head -1|cut -f5 -d '='`
grep '%' $1 | grep -v 'part'| grep -v 'norm'|head -n $n0 > tmp.lost1
grep '%' $1 | grep -v 'part'| grep -v 'norm'|tail -n +$(($n0+1)) > tmp.lost
perl ../scripts/inilost.pl $DISTR tmp.lost1 inilost.1.out
perl ../scripts/inilost.pl $DISTR tmp.lost inilost.out
#perl /home/valishev/lifetrac/bin/inilost.pl $DISTR tmp.lost1 inilost.1.out
#perl /home/valishev/lifetrac/bin/inilost.pl $DISTR tmp.lost inilost.out
rm tmp.lost tmp.lost1
