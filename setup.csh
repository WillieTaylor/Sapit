# 1 = new ID, 2 = pdb ID, 3 = chain

mkdir dir$argv[1]
cd dir$argv[1]
ln -s .. home
ln -s ~/util
mkdir pdb
cp ~/Downloads/$argv[2].pdb pdb
cd pdb
cat $argv[2].pdb | cut -c 1-66 | grep ' CA ' > tmp1.pdb
eval "grep  ' $argv[3] ' tmp1.pdb" > tmp2.pdb
grep -v ' CA B' tmp2.pdb | sed 's/ CA A/ CA  /' > tmp1.pdb
cat tmp1.pdb | sed 's/HETATM/ATOM  /' | grep '^ATOM ' | sed 's/MSE/MET/' > $argv[1]
rm tmp[123].pdb
