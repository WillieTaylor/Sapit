foreach str (`grep ' [1-9][a-z]' $argv[1].fssp | awk '{print $2}'`)
	echo $str
	set prot1 = `echo $str | tr "-" " "`
	if ( $#prot1 > 1 ) then
        	echo Protein is $prot1[1] chain $prot1[2]
        	set str1 = "$prot1[1]$prot1[2]"
        	if ( -e pdb/$str1 ) then
                	echo found
        	else
                	echo getting pdb
                	tcsh util/getpdb.csh $prot1[1] $prot1[2]
        	endif
	else
		echo Protein is $prot1
        	set str1 = "$prot1"
        	if ( -e pdb/$str1 ) then
                	echo found
        	else
                	echo getting pdb
                	tcsh util/getpdb.csh $prot1
		endif
	endif
end
rm *.[eo][0-9]*
foreach str (`grep ' [1-9][a-z]' $argv[1].fssp | awk '{print $2}'`)
	set pdb = `echo $str | tr -d "-"`
	if ( ! -e pdb/$pdb ) then
		echo PDB $pdb not found
		continue
	endif
	set name = `echo sap$argv[1]$str | tr -d "-"`
	echo Running $argv[1] on $str as $name
	tcsh home/sapit.csh $argv[1] $str 5
#	qsub -N $name -l cput=22:00:00 -v STR1=$argv[1],STR2=$str,CYCS=$argv[2] home/sapit.csh
#	sleep $argv[3]
end
