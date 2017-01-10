#!/bin/csh

# 1 = prot1, 2 = prot2 (both in /pdb)

# gnuplot>  plot [0:200] 'nat1rev2.rat', 'pdb1pdb2.dat'notitle  with lines, 'pdb1pdb2.dat' notitle, 'nat1rve2.rat', 'nta1rev2.rat',14*(1-exp(-x*x*0.0001)) notitle,'nta1rve2.rat','nta1rev2.rat',28*(1-exp(-x*x*0.0001)) notitle

# gnuplot> plot [10:90][0:15]  'nat1rev2.rat' not, 'pdb1pdb2.dat' with lines not, 'pdb1pdb2.dat' not, 'nat1rve2.rat' not, 'nta1rev2.rat' not, 14*(1-exp(-x*x*0.0001)) not,'nta1rve2.rat' not, 28*(1-exp(-x*x*0.0001)) not, 'rev1nat2.rat' not, 'rev1nta2.rat' not,  'rve1nat2.rat' not, 'rve1nta2.rat' not


set prot1 = `echo $argv[1] | tr "-" " "`
# set prot1 = `echo $STR1 | tr "-" " "`
if ( $#prot1 > 1 ) then
	echo Protein 1 is $prot1[1] chain $prot1[2]
	set str1 = "$prot1[1]$prot1[2]"
	if ( -e pdb/$str1 ) then
		echo found
	else
		echo added
		tcsh util/getpdb.csh $prot1[1] $prot1[2]
	endif
else
	echo Protein 1 is $prot1[1] 
	set str1 = "$prot1[1]"
	if ( -e pdb/$str1 ) then
		echo found
	else
		echo added
		tcsh util/getpdb.csh $prot1[1]
	endif
endif
set prot2 = `echo $argv[2] | tr "-" " "`
# set prot2 = `echo $STR2 | tr "-" " "`
if ( $#prot2 > 1 ) then
	echo Protein 2 is $prot2[1] chain $prot2[2]
	set str2 = "$prot2[1]$prot2[2]"
	if ( -e pdb/$str2 ) then
		echo found
	else
		echo added
		tcsh util/getpdb.csh $prot2[1] $prot2[2]
	endif
else
	echo Protein 2 is $prot2[1] 
	set str2 = "$prot2[1]"
	if ( -e pdb/$str2 ) then
		echo found
	else
		echo added
		tcsh util/getpdb.csh $prot2[1]
	endif
endif
set dir = "$str1$str2"
echo Directory is $dir
# if ( -e $dir ) mv $dir $dir.last
if ( -e $dir ) rm -r $dir
mkdir $dir
if ( -e pdb/$str1 ) cp pdb/$str1 $dir/$str1
if ( -e pdb/$str2 ) cp pdb/$str2 $dir/$str2
cd $dir
ln -s ~/sapit main
ln -s ~/util util

main/revcas $str1 N
mv test.out $str1.nat
main/revcas $str1 R
mv test.out $str1.rev
main/revcas $str2 N
mv test.out $str2.nat
main/revcas $str2 R
mv test.out $str2.rev

main/revcas $str1 n
mv test.out $str1.nta
main/revcas $str1 r
mv test.out $str1.rve
main/revcas $str2 n
mv test.out $str2.nta
main/revcas $str2 r
mv test.out $str2.rve

#exit

set nres = `cat $str1 | grep ATOM | grep ' CA ' | wc -l`
echo $nres residues in str1
@ m = $nres[1] + 1
set nres = `cat $str2 | grep ATOM | grep ' CA ' | wc -l`
echo $nres[1] residues in str2
@ n = $nres[1] + 1

@ CYCS = 5

echo nat1 vs rev2
if ( -e rms.dat) rm rms.dat
tcsh main/sapit2.csh $m $n $str1.nat $str2.rev $CYCS
awk '{print $8, $5, $12, $14, $15, $16, $17}' rms.dat > nat1rev2.rat

	echo nat2 vs rev1
	rm rms.dat
	tcsh main/sapit2.csh $n $m $str2.nat $str1.rev $CYCS
	awk '{print $8, $5, $12, $14, $15, $16, $17}' rms.dat > rev1nat2.rat

echo nta1 vs rev2
rm rms.dat
tcsh main/sapit2.csh $m $n $str1.nta $str2.rev $CYCS
awk '{print $8, $5, $12, $14, $15, $16, $17}' rms.dat > nta1rev2.rat

	echo nta2 vs rev1
	rm rms.dat
	tcsh main/sapit2.csh $n $m $str2.nta $str1.rev $CYCS
	awk '{print $8, $5, $12, $14, $15, $16, $17}' rms.dat > rev1nta2.rat

echo nat1 vs rve2
rm rms.dat
tcsh main/sapit2.csh $m $n $str1.nat $str2.rve $CYCS
awk '{print $8, $5, $12, $14, $15, $16, $17}' rms.dat > nat1rve2.rat

	echo nat2 vs rve1
	rm rms.dat
	tcsh main/sapit2.csh $n $m $str2.nat $str1.rve $CYCS
	awk '{print $8, $5, $12, $14, $15, $16, $17}' rms.dat > rve1nat2.rat

echo nta1 vs rve2
rm rms.dat
tcsh main/sapit2.csh $m $n $str1.nta $str2.rve $CYCS
awk '{print $8, $5, $12, $14, $15, $16, $17}' rms.dat > nta1rve2.rat

	echo nta2 vs rve1
	rm rms.dat
	tcsh main/sapit2.csh $n $m $str2.nta $str1.rve $CYCS
	awk '{print $8, $5, $12, $14, $15, $16, $17}' rms.dat > rve1nta2.rat

cat *.rat > random.dat

# pdb[12] differs from nat[12] only in cyclisation distortion 
#	echo pdb1 vs rev1
#	rm rms.dat
#	tcsh main/sapit1.csh $m $m $str1 $str1.rev $CYCS
#	awk '{print $8, $5, $12, $14, $15, $16, $17}' rms.dat > pdb1rev1.rat

#	echo pdb1 vs rev2
#	rm rms.dat
#	tcsh main/sapit1.csh $m $n $str1 $str2.rev $CYCS
#	awk '{print $8, $5, $12, $14, $15, $16, $17}' rms.dat > pdb1rev2.rat

#	echo pdb2 vs rev2
#	rm rms.dat
#	tcsh main/sapit1.csh $n $n $str2 $str2.rev $CYCS
#	awk '{print $8, $5, $12, $14, $15, $16, $17}' rms.dat > pdb2rev2.rat

#	echo pdb2 vs rev1
#	rm rms.dat
#	tcsh main/sapit1.csh $n $m $str2 $str1.rev $CYCS
#	awk '{print $8, $5, $12, $14, $15, $16, $17}' rms.dat > pdb2rev1.rat

echo pdb1 vs pdb2
rm pdb1pdb2.dat
main/sapit $str1 $str2 > sap.log
cat sap.log | grep ' over best ' | awk '{print $8, $5, $11}' >> pdb1pdb2.dat
cat sap.log | grep ' over all '  | awk '{print $8, $5, $12}' >> pdb1pdb2.dat
tail -1 pdb1pdb2.dat > overall.dat

echo pdb2 vs pdb1
rm pdb2pdb1.dat
main/sapit $str2 $str1 | grep ' over all ' | awk '{print $8, $5, $12, $14, $15, $16, $17}' >> pdb2pdb1.dat

#	echo pdb1 vs nat1
#	rm pdb1nat1.dat
#	main/sapit $str1 $str1.nat | grep ' over best ' | awk '{print $8, $5, $11}' >> pdb1nat1.dat
#	echo pdb2 vs nat2
#	rm pdb2nat2.dat
#	main/sapit $str2 $str2.nat | grep ' over best ' | awk '{print $8, $5, $11}' >> pdb2nat2.dat

if ( $n > $m ) then
	echo "plot [0:$n][0:15]\\" > temp.gnu
else
	echo "plot [0:$m][0:15]\\" > temp.gnu
endif
cat ../plot.gnu >> temp.gnu
# gnuplot temp.gnu
cd ..
rm *.[eo][0-9][0-9]*
