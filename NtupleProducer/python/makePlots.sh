# SCRIPT MAKES PLOTS FOR JET AND HT/MHT TRIGGER PERFORMANCE FOR FIRST 12 JETS

# Define variables
sig=/eos/user/l/lroberts/P2_Jets/WideCone/working_dir/ttbar/ttbar.root
bkg=/eos/user/l/lroberts/P2_Jets/WideCone/working_dir/singneut/singneut.root
outdir=plots/bigplots
what=doubBinSize

# Make JECs
python3 scripts/makeJecs.py "$sig" -A -o "$what"JECs.root

# Make plots
for i in {1..12}; do
    echo "Making plots for jet $i"
    python3 scripts/jetHtSuite.py "$sig" "$bkg" "$outdir"/jet"$i" -w "$what" -v jet"$i" --jecs "$what"JECs.root
done
python3 scripts/jetHtSuite.py "$sig" "$bkg" "$outdir"/ht -w "$what" -v ht
python3 scripts/jetHtSuite.py "$sig" "$bkg" "$outdir"/mht -w "$what" -v mht