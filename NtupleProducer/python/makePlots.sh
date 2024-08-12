# SCRIPT MAKES PLOTS FOR JET AND HT/MHT TRIGGER PERFORMANCE FOR FIRST 12 JETS

# Define variables
sig=/eos/user/l/lroberts/P2_Jets/CMSSW_14_0_0_pre3/src/FastPUPPI/condor/jobs/variousMassCutsTTbar_1722943032/data/ttbar.root
bkg=/eos/user/l/lroberts/P2_Jets/CMSSW_14_0_0_pre3/src/FastPUPPI/condor/jobs/variousMassCutsSingNeut_1722812263/data/singneut.root
outdir=plots/MASSCUTSTTBAR
what=jetMassCut

# Make JECs
# python3 scripts/makeJecs.py "$sig" -A -o "$what"JECs.root

# Make plots
for i in {1..1}; do
    echo "Making plots for jet $i"
    python3 scripts/jetHtSuite.py "$sig" "$bkg" "$outdir"/jet"$i" -w "$what" -v jet"$i" --jecs "$what"JECs.root --gendr 0.4 --raw
done
# python3 scripts/jetHtSuite.py "$sig" "$bkg" "$outdir"/ht -w "$what" -v ht --raw
# python3 scripts/jetHtSuite.py "$sig" "$bkg" "$outdir"/mht -w "$what" -v mht --raw