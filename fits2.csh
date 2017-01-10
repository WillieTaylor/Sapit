# 1 = dir (eg 3chy5nul), 2 = score mode (U = damped sqrt), 3 = EVD (Y/N = 1/0)


#gnuplot> plot [20:180] 'dist.out' using 1:2 notitle with lines, 'dist.out' using 1:3 notitle with lines, 'dist.out' using 1:4 notitle with lines, 'dist.out' using 1:5 notitle with lines

# set name = `echo $argv[1] | sed 's/[A-Z]$//' | sed 's/....$//' | sed 's/[A-Z]/-&/'`
# set code = `echo $argv[1] | sed 's/....//' | sed 's/^[A-Z]//' | sed 's/[A-Z]/-&/'`
# [A-Z] is no longer case sensitive so using...

#set name = `echo $argv[1] | sed 's/[ABCDEFGHIJKLMNOPQRSTUVWXYZ]$//' | sed 's/....$//'`
#set code = `echo $argv[1] | sed 's/....//' | sed 's/^[ABCDEFGHIJKLMNOPQRSTUVWXYZ]//' | sed 's/[ABCDEFGHIJKLMNOPQRSTUVWXYZ]/-&/'`

#echo $name $code
#echo grep $code $name.fssp "| head -1 > name.txt" > temp.csh
#tcsh temp.csh

if ( ! -e $argv[1]/nat1rev2.rat ) exit
if ( ! -e $argv[1]/pdb2pdb1.dat ) exit

cat $argv[1]/[nr]*[12][nr]*[12].rat > temp.dat
wc -l temp.dat
sort -k 2 -nu temp.dat > sort.dat
wc -l sort.dat
grep -v nan sort.dat > $argv[1]/random.dat
rm sort.dat

# home/fit99 $argv[1] T > fit99.sap.log
# home/cdfit $argv[1] T > cdfit.sap.log
# set sap = `grep 'z =' cdfit.sap.log`
# echo SAP score $sap

# home/fit99 $argv[1] R > fit99.rms.log
# home/cdfit $argv[1] R > cdfit.rms.log
# set rms = `grep 'z =' cdfit.rms.log`

main/fits2 $argv[1] $argv[2] $argv[3] > fits2.rms.log
set rms = `grep 'ymax =' fits2.rms.log`
set fit = `grep 'fit =' fits2.rms.log`
if ( $#fit < 1 ) then
	set fits = "-9.90000"
else
	set fits = $fit[12]
endif
echo RMS scores $rms $fits

#echo $rms[3] $rms[6] $rms[9] $rms[12] `cat name.txt` >> zed.dat
echo $rms[3] $rms[6] $rms[9] $rms[12] >> zed.dat
echo ""
grep copt fits2.rms.log >> bcopt.dat
echo $fits >> fits.dat
