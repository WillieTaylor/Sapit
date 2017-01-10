tcsh home/fsspZ.csh $argv[1] $argv[2] 1
mv zed.dat $argv[1].$argv[2].evd
set fevd = `cat fits.dat | grep -v 9.9000 | awk '{s+=$1 ; t+=$1*$1 ; n++ ; print s/n, n, t/n}' | tail -1`
tcsh home/fsspZ.csh $argv[1] $argv[2] 0
mv zed.dat $argv[1].$argv[2].zed
set fzed = `cat fits.dat | grep -v 9.9000 | awk '{s+=$1 ; t+=$1*$1 ; n++ ; print s/n, n, t/n}' | tail -1`
echo Fits: evd = $fevd, norm = $fzed
awk '{print $1,$2,$3,$4}' $argv[1].$argv[2].evd > tmp
paste tmp $argv[1].$argv[2].zed > $argv[1].$argv[2].dat
grep -v 9.9000 $argv[1].$argv[2].dat > $argv[1].$argv[2].no9
