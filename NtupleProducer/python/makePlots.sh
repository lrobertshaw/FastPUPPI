# SCRIPT MAKES PLOTS FOR JET AND HT/MHT TRIGGER PERFORMANCE FOR FIRST 12 JETS

# Define variables
sig=/eos/user/l/lroberts/P2_Jets/CMSSW_14_0_0_pre3/src/FastPUPPI/condor/jobs/HSC_DoubleBinned_HH4b_1723822787/data/hh4b.root
bkg=/eos/user/l/lroberts/P2_Jets/CMSSW_14_0_0_pre3/src/FastPUPPI/condor/jobs/HSC_DoubleBinned_SingNeut_1723822744/data/singneut.root
outdir=plots/HSCDoubBinnedWithJECs_hh4b
what=doubBinSize

# Make JECs
python3 scripts/makeJecs.py "$sig" -A -o "$what"JECs.root

# Make plots
for i in {1..10}; do
    echo "Making plots for jet $i"
    python3 scripts/jetHtSuite.py "$sig" "$bkg" "$outdir"/jet"$i" -w "$what" -v jet"$i" --jecs "$what"JECs.root --gendr 0.4
done
python3 scripts/jetHtSuite.py "$sig" "$bkg" "$outdir"/ht -w "$what" -v ht --raw
python3 scripts/jetHtSuite.py "$sig" "$bkg" "$outdir"/mht -w "$what" -v mht --raw