@ i = 0
cp $argv[3] temp1.out
while (1)
	@ i++
	set left1 = `wc -l temp1.out`
	if ( $left1[1] < $argv[1] ) exit
	head -$argv[1] temp1.out > temp1.pdb
	@ j = 0
	cp $argv[4] temp2.out
	while (1)
		@ j++
		set left2 = `wc -l temp2.out`
		if ( $left2[1] < $argv[2] ) break
		head -$argv[2] temp2.out > temp2.pdb
		main/sapit temp1.pdb temp2.pdb $argv[5] | grep ' Un-weighted RMSd' | tail -1 >> rms.dat
		@ rest = $left2[1] - $argv[2]
		if ( $rest < 10 ) break
		tail -$rest temp2.out > tmp
		mv tmp temp2.out
	end
	@ rest = $left1[1] - $argv[1]
	if ( $rest < 10 ) exit
	tail -$rest temp1.out > tmp
	mv tmp temp1.out
end
