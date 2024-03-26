import ROOT
import sys
import os
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetOptStat(0)
if __name__ == "__main__":
    sys.path.append(os.environ['ANALYSIS_PATH'])

    import Common.Utilities as Utilities

x_bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5 , 6.5, 7.5]
x_bins_vec = Utilities.ListToVector(x_bins, "double")


f1 = ROOT.TFile.Open("Studies/Triggers/files/genkind_distrib_0/gen_kinds_TTTo2L2Nu_tautau_Medium_Medium.root", "READ")
histo1 = f1.Get("tau1_gen_kind")
histo2 = f1.Get("tau2_gen_kind")
hist2D_1 = ROOT.TH2D("tau1_tau2_genKind_VVLoose_0", "tau1_tau2_genKind_VVLoose_0", x_bins_vec.size()-1, x_bins_vec.data(), x_bins_vec.size()-1, x_bins_vec.data())
for x_bin_1 in range(0, histo1.GetNbinsX()+1):
    for x_bin_2 in range(0, histo2.GetNbinsX()+1):
        content_1 = histo1.GetBinContent(x_bin_1) * histo2.GetBinContent(x_bin_2)
        hist2D_1.Fill(histo1.GetBinCenter(x_bin_1), histo2.GetBinCenter(x_bin_2), content_1)

hist2D_1.GetXaxis().SetTitle("#tau_{1} genKind")
hist2D_1.GetYaxis().SetTitle("#tau_{2} genKind")

f2 = ROOT.TFile.Open("Studies/Triggers/files/genkind_distrib_0/gen_kinds_TTTo2L2Nu_tautau_Medium_Medium_VSeLoose.root", "READ")
histo3 = f2.Get("tau1_gen_kind")
histo4 = f2.Get("tau2_gen_kind")
hist2D_2 = ROOT.TH2D("tau1_tau2_genKind_Loose_0", "tau1_tau2_genKind_Loose_0", x_bins_vec.size()-1, x_bins_vec.data(), x_bins_vec.size()-1, x_bins_vec.data())
for x_bin_3 in range(0, histo3.GetNbinsX()+1):
    for x_bin_4 in range(0, histo4.GetNbinsX()+1):
        content_2 = histo3.GetBinContent(x_bin_3) * histo4.GetBinContent(x_bin_4)
        hist2D_2.Fill(histo3.GetBinCenter(x_bin_3), histo4.GetBinCenter(x_bin_4), content_2)

hist2D_2.GetXaxis().SetTitle("#tau_{1} genKind")
hist2D_2.GetYaxis().SetTitle("#tau_{2} genKind")



f3 = ROOT.TFile.Open("Studies/Triggers/files/genkind_distrib_0/gen_kinds_TTTo2L2Nu_tautau_Medium_Medium_VSeVLoose.root", "READ")
histo5 = f3.Get("tau1_gen_kind")
histo6 = f3.Get("tau2_gen_kind")
hist2D_3 = ROOT.TH2D("tau1_tau2_genKind_VLoose_0", "tau1_tau2_genKind_VLoose_0", x_bins_vec.size()-1, x_bins_vec.data(), x_bins_vec.size()-1, x_bins_vec.data())
for x_bin_5 in range(0, histo5.GetNbinsX()+1):
    for x_bin_6 in range(0, histo6.GetNbinsX()+1):
        content_3 = histo5.GetBinContent(x_bin_5) * histo6.GetBinContent(x_bin_6)
        hist2D_3.Fill(histo5.GetBinCenter(x_bin_5), histo6.GetBinCenter(x_bin_6), content_3)

hist2D_3.GetXaxis().SetTitle("#tau_{1} genKind")
hist2D_3.GetYaxis().SetTitle("#tau_{2} genKind")

################################### sum == 1 #####################################

f4 = ROOT.TFile.Open("Studies/Triggers/files/genkind_distrib_1/gen_kinds_TTTo2L2Nu_tautau_Medium_Medium.root", "READ")
histo7 = f4.Get("tau1_gen_kind")
histo8 = f4.Get("tau2_gen_kind")
hist2D_4 = ROOT.TH2D("tau1_tau2_genKind_VVLoose_1", "tau1_tau2_genKind_VVLoose_1", x_bins_vec.size()-1, x_bins_vec.data(), x_bins_vec.size()-1, x_bins_vec.data())
for x_bin_7 in range(0, histo7.GetNbinsX()+1):
    for x_bin_8 in range(0, histo8.GetNbinsX()+1):
        content_4 = histo7.GetBinContent(x_bin_7) * histo8.GetBinContent(x_bin_8)
        hist2D_4.Fill(histo7.GetBinCenter(x_bin_7), histo8.GetBinCenter(x_bin_8), content_4)

hist2D_4.GetXaxis().SetTitle("#tau_{1} genKind")
hist2D_4.GetYaxis().SetTitle("#tau_{2} genKind")

f5 = ROOT.TFile.Open("Studies/Triggers/files/genkind_distrib_1/gen_kinds_TTTo2L2Nu_tautau_Medium_Medium_VSeLoose.root", "READ")
histo9 = f5.Get("tau1_gen_kind")
histo10 = f5.Get("tau2_gen_kind")
hist2D_5 = ROOT.TH2D("tau1_tau2_genKind_Loose_1", "tau1_tau2_genKind_Loose_1", x_bins_vec.size()-1, x_bins_vec.data(), x_bins_vec.size()-1, x_bins_vec.data())
for x_bin_9 in range(0, histo9.GetNbinsX()+1):
    for x_bin_10 in range(0, histo10.GetNbinsX()+1):
        content_5 = histo9.GetBinContent(x_bin_9) * histo10.GetBinContent(x_bin_10)
        hist2D_5.Fill(histo9.GetBinCenter(x_bin_9), histo10.GetBinCenter(x_bin_10), content_5)

hist2D_5.GetXaxis().SetTitle("#tau_{1} genKind")
hist2D_5.GetYaxis().SetTitle("#tau_{2} genKind")


f6 = ROOT.TFile.Open("Studies/Triggers/files/genkind_distrib_1/gen_kinds_TTTo2L2Nu_tautau_Medium_Medium_VSeVLoose.root", "READ")
histo11 = f6.Get("tau1_gen_kind")
histo12 = f6.Get("tau2_gen_kind")
hist2D_6 = ROOT.TH2D("tau1_tau2_genKind_VLoose_1", "tau1_tau2_genKind_VLoose_1", x_bins_vec.size()-1, x_bins_vec.data(), x_bins_vec.size()-1, x_bins_vec.data())
for x_bin_11 in range(0, histo11.GetNbinsX()+1):
    for x_bin_12 in range(0, histo12.GetNbinsX()+1):
        content_6 = histo11.GetBinContent(x_bin_11) * histo12.GetBinContent(x_bin_12)
        hist2D_6.Fill(histo11.GetBinCenter(x_bin_11), histo12.GetBinCenter(x_bin_12), content_6)

hist2D_6.GetXaxis().SetTitle("#tau_{1} genKind")
hist2D_6.GetYaxis().SetTitle("#tau_{2} genKind")


################################### sum == 2 #####################################


f7 = ROOT.TFile.Open("Studies/Triggers/files/genkind_distrib_2/gen_kinds_TTTo2L2Nu_tautau_Medium_Medium.root", "READ")
histo13 = f7.Get("tau1_gen_kind")
histo14 = f7.Get("tau2_gen_kind")
hist2D_7 = ROOT.TH2D("tau1_tau2_genKind_VVLoose_2", "tau1_tau2_genKind_VVLoose_2", x_bins_vec.size()-1, x_bins_vec.data(), x_bins_vec.size()-1, x_bins_vec.data())
for x_bin_13 in range(0, histo13.GetNbinsX()+1):
    for x_bin_14 in range(0, histo14.GetNbinsX()+1):
        content_4 = histo13.GetBinContent(x_bin_13) * histo14.GetBinContent(x_bin_14)
        hist2D_7.Fill(histo13.GetBinCenter(x_bin_13), histo14.GetBinCenter(x_bin_14), content_4)

hist2D_7.GetXaxis().SetTitle("#tau_{1} genKind")
hist2D_7.GetYaxis().SetTitle("#tau_{2} genKind")

f8 = ROOT.TFile.Open("Studies/Triggers/files/genkind_distrib_2/gen_kinds_TTTo2L2Nu_tautau_Medium_Medium_VSeLoose.root", "READ")
histo15 = f8.Get("tau1_gen_kind")
histo16 = f8.Get("tau2_gen_kind")
hist2D_8 = ROOT.TH2D("tau1_tau2_genKind_Loose_2", "tau1_tau2_genKind_Loose_2", x_bins_vec.size()-1, x_bins_vec.data(), x_bins_vec.size()-1, x_bins_vec.data())
for x_bin_15 in range(0, histo15.GetNbinsX()+1):
    for x_bin_16 in range(0, histo16.GetNbinsX()+1):
        content_8 = histo15.GetBinContent(x_bin_15) * histo16.GetBinContent(x_bin_16)
        hist2D_8.Fill(histo15.GetBinCenter(x_bin_15), histo16.GetBinCenter(x_bin_16), content_8)

hist2D_8.GetXaxis().SetTitle("#tau_{1} genKind")
hist2D_8.GetYaxis().SetTitle("#tau_{2} genKind")


f9 = ROOT.TFile.Open("Studies/Triggers/files/genkind_distrib_2/gen_kinds_TTTo2L2Nu_tautau_Medium_Medium_VSeVLoose.root", "READ")
histo17 = f9.Get("tau1_gen_kind")
histo18 = f9.Get("tau2_gen_kind")
hist2D_9 = ROOT.TH2D("tau1_tau2_genKind_VLoose_2", "tau1_tau2_genKind_VLoose_2", x_bins_vec.size()-1, x_bins_vec.data(), x_bins_vec.size()-1, x_bins_vec.data())
for x_bin_17 in range(0, histo17.GetNbinsX()+1):
    for x_bin_18 in range(0, histo18.GetNbinsX()+1):
        content_9 = histo17.GetBinContent(x_bin_17) * histo18.GetBinContent(x_bin_18)
        hist2D_9.Fill(histo17.GetBinCenter(x_bin_17), histo18.GetBinCenter(x_bin_18), content_9)

hist2D_9.GetXaxis().SetTitle("#tau_{1} genKind")
hist2D_9.GetYaxis().SetTitle("#tau_{2} genKind")


################## sum of all histograms ###############

hist2D_all_VVLoose = hist2D_1.Clone("hist2D_all_VVLoose")
hist2D_all_VVLoose.Add(hist2D_4)
hist2D_all_VVLoose.Add(hist2D_7)

hist2D_all_Loose = hist2D_2.Clone("hist2D_all_Loose")
hist2D_all_Loose.Add(hist2D_5)
hist2D_all_Loose.Add(hist2D_8)

hist2D_all_VLoose = hist2D_3.Clone("hist2D_all_VLoose")
hist2D_all_VLoose.Add(hist2D_6)
hist2D_all_VLoose.Add(hist2D_9)


####### printOuts ######
print(f"nEvents for VVLoose are : 0 = {histo1.GetEntries()}, 1 ={histo7.GetEntries()}, 2 = {histo13.GetEntries()} ")
total_VVLoose = histo1.GetEntries() + histo7.GetEntries() + histo13.GetEntries()
print(f"total nEvents for VVLoose are : 0 = {total_VVLoose} ")
print(f"percentage for VVLoose are : 0 = {100*histo1.GetEntries()/total_VVLoose}, 1 ={100* histo7.GetEntries()/total_VVLoose}, 2 = {100*histo13.GetEntries()/total_VVLoose} ")

print()
print(f"nEvents for Loose are : 0 = {histo3.GetEntries()}, 1 ={histo9.GetEntries()}, 2 = {histo15.GetEntries()}")
total_Loose = histo3.GetEntries() + histo9.GetEntries() + histo15.GetEntries()
print(f"total nEvents for Loose are : 0 = {total_Loose} ")
print(f"percentage for Loose are : 0 = {100*histo3.GetEntries()/total_Loose}, 1 ={100* histo9.GetEntries()/total_Loose}, 2 = {100*histo15.GetEntries()/total_Loose}")
print()
total_VLoose = histo5.GetEntries() + histo11.GetEntries() + histo17.GetEntries()
print(f"nEvents for VLoose are : 0 = {histo5.GetEntries()}, 1 ={histo11.GetEntries()}, 2 = {histo17.GetEntries()}")
print(f"total nEvents for VLoose are : 0 = {histo5.GetEntries() + histo11.GetEntries() + histo17.GetEntries()} ")
print(f"percentage for VLoose are : 0 = {100*histo5.GetEntries()/total_VLoose}, 1 ={100* histo11.GetEntries()/total_VLoose}, 2 = {100*histo17.GetEntries()/total_VLoose}")


outFile = ROOT.TFile.Open(f"Studies/Triggers/files/2Ddistrib.root", "RECREATE")
hist2D_1.Write()
hist2D_2.Write()
hist2D_3.Write()

hist2D_4.Write()
hist2D_5.Write()
hist2D_6.Write()

hist2D_7.Write()
hist2D_8.Write()
hist2D_9.Write()

hist2D_all_VVLoose.Write()
hist2D_all_Loose.Write()
hist2D_all_VLoose.Write()
outFile.Close()
