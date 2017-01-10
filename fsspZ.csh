rm zed.dat
rm zeds.dat
rm fits.dat
rm bcopt.dat
foreach str (`grep ' [1-9][a-z]' $argv[1].fssp | awk '{print $2}'`)
	set dir = `echo $argv[1]$str | tr -d '-'`
	echo Running $argv[1] $str from $dir
	if ( -e $dir ) then
		tcsh home/fits2.csh $dir $argv[2] $argv[3]
	else
		echo directory $dir not found
		# exit
	endif
#	set sap = `grep 'z =' cdfit.sap.log`
#	set rms = `grep 'z =' cdfit.rms.log`
#	echo $sap[3] $rms[3] >> zeds.dat
end
