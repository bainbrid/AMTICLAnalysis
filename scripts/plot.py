import ROOT
import sys
from collections import OrderedDict as odict

################################################################################
# I/O

inFileName=None
outFileName=None
if len(sys.argv) == 1:
    inFileName="../TiCLTreeProducer/D49_FineCalo/Reco/pt50_eta21/Debug/RecoTree_Unconv.root"
    outFileName="./output.root"
    plotFileName = "./plots.pdf"
elif len(sys.argv) == 2:
    inFileName = sys.argv[1]
    outFileName="./output.root"
    plotFileName = "./plots.pdf"
elif len(sys.argv) == 3:
    inFileName = sys.argv[1]
    outFileName = sys.argv[2]
    plotFileName = "./plots.pdf"
elif len(sys.argv) == 4:
    inFileName = sys.argv[1]
    outFileName = sys.argv[2]
    plotFileName = sys.argv[3]
else:
    print "USAGE: %s <input root file> <output root file> <output pdf file>"%(sys.argv[0])
    sys.exit(1)
print "Reading from", inFileName, "and writing to", outFileName, " and ", plotFileName

inFile = ROOT.TFile.Open(inFileName ,"READ")
tree = inFile.Get("RecoTree")

################################################################################
# Utility

class Histo :
    def __init__(self,
                 htype, # allowed: "TH1F","TH2F","TProfile"
                 name,title="",
                 xvar="",xbin=0,xmin=0.,xmax=0.,
                 yvar="",ybin=0,ymin=0.,ymax=0.,
                 logx=False,logy=False,logz=False,
                 cut="",
                 options="h",
                 optstat=1110000):
        self.htype=htype
        self.name=name
        self.title=title
        self.xvar=xvar
        self.xbin=xbin
        self.xmin=xmin
        self.xmax=xmax
        self.yvar=yvar
        self.ybin=ybin
        self.ymin=ymin
        self.ymax=ymax
        self.logx=logx
        self.logy=logy
        self.logz=logz
        self.cut=cut
        self.options=options
        self.optstat=optstat

################################################################################
# Configure histograms

config = [
        Histo("TH1F",
              "nTotTS","nTotTS",
              "nTotTS",10,-0.5,9.5,
              ),
        Histo("TH1F",
              "eTotTS","eTotTS",
              "eTotTS",110,0.,110.,
              ),
        Histo("TH1F",
              "eTS1","eTS1",
              "eTS1",110,0.,110.,
              ),
        Histo("TH1F",
              "eTS1_nTotTS_eq1","eTS1",
              "eTS1",110,0.,110.,
              cut="nTotTS==1",
              ),
        Histo("TH1F",
              "eTS1_nTotTS_eq2","eTS1",
              "eTS1",110,0.,110.,
              cut="nTotTS==2",
              ),
        Histo("TH1F",
              "eTS2","eTS2",
              "eTS2",110,0.,110.,
              ),
        Histo("TH1F",
              "eTS2_nTotTS_eq2","eTS2",
              "eTS2",110,0.,110.,
              cut="nTotTS==2",
              ),
        Histo("TH1F",
              "dist_2d_lt10","dist_2d",
              "dist_2d",100,0.,10.,
              ),
        Histo("TH2F",
              "eTotTS_vs_dist_2d_nTotTS_eq1","eTotTS_vs_dist_2d",
              "dist_2d",100,0.,10.,
              "eTotTS",110,0.,110.,
              cut="nTotTS==1",
              options="colz",
              ),
        Histo("TH2F",
              "eTotTS_vs_dist_2d_nTotTS_eq2","eTotTS_vs_dist_2d",
              "dist_2d",100,0.,10.,
              "eTotTS",110,0.,110.,
              cut="nTotTS==2",
              options="colz",
              ),
        Histo("TH2F",
              "nTotTS_vs_dist_2d","nTotTS_vs_dist_2d",
              "dist_2d",100,0.,10.,
              "nTotTS",4,-0.5,3.5,
              options="colz",
              ),
        Histo("TH2F",
              "eTotTS_vs_dist_2d","eTotTS_vs_dist_2d",
              "dist_2d",100,0.,10.,
              "eTotTS",110,0.,110.,
              options="colz",
              ),
        Histo("TProfile",
              "prof_nTotTS_vs_dist_2d","nTotTS_vs_dist_2d",
              "dist_2d",100,0.,10.,
              "nTotTS",ymin=-0.5,ymax=3.5,
              options="",
              ),
#        Histo("TProfile",
#              "prof_nTotTS_vs_dist_3d","nTotTS_vs_dist_3d",
#              "dist_3d",100,0.,10.,
#              "nTotTS",ymin=-0.5,ymax=3.5,
#              options="",
#              ),
        ]

################################################################################
# Create histograms

histos = {}
for his in config:
    if his.htype=="TH1F":
        histos[his.name] = ROOT.TH1D(his.name,his.title,
                                     his.xbin,his.xmin,his.xmax)
        histos[his.name].Sumw2()
    elif his.htype=="TH2F":
        histos[his.name] = ROOT.TH2D(his.name,his.title,
                                     his.xbin,his.xmin,his.xmax,
                                     his.ybin,his.ymin,his.ymax)
        histos[his.name].Sumw2()
    elif his.htype=="TProfile":
        histos[his.name] = ROOT.TProfile(his.name,his.title,
                                         his.xbin,his.xmin,his.xmax,
                                         his.ymin,his.ymax)
        histos[his.name].Sumw2()
    else:
        print "Unknown histo type!",his.htype
        sys.exit(1)

################################################################################
# Event loop

for his in config:
    if his.htype=="TH1F":
        passed = tree.Draw(his.xvar+">>"+his.name,his.cut,"goff")
    elif his.htype=="TH2F":
        passed = tree.Draw(his.yvar+":"+his.xvar+">>"+his.name,his.cut,"goff")
    elif his.htype=="TProfile":
        passed = tree.Draw(his.yvar+":"+his.xvar+">>"+his.name,his.cut,"goff")
    else:
        print "Unknown histo type!",his.htype
        sys.exit(1)

################################################################################
# Cosmetics and write to PDF

canvas = ROOT.TCanvas("canvas")
canvas.cd()
canvas.Print(plotFileName+"[") # Write all plots to same file

#def RearrangeStatBox(hist,pad):
#	StatsBox = hist.FindObject("stats").Clone()
#	StatsBox.SetX1NDC(0.05)
#	StatsBox.SetX2NDC(0.95)
#	StatsBox.SetY1NDC(0.05)
#	StatsBox.SetY2NDC(0.95)
#	pad.cd()
#	StatsBox.Draw()
#	#hist.SetStats(False)

for his in config:
    ROOT.gStyle.SetOptStat(his.optstat)
    histos[his.name].Draw(his.options)
    histos[his.name].SetTitle(his.title if len(his.cut)==0 else his.title+" {"+his.cut+"}")
    canvas.SetLogx(his.logx)
    canvas.SetLogy(his.logy)
    canvas.SetLogz(his.logz)
    canvas.Print(plotFileName)
    canvas.Update()
    histos[his.name].SetDirectory(0) # Do this after TTree.Draw
    #ROOT.gPad.Update()
    #RearrangeStatBox(histos[his.name],canvas)
    #canvas.Update()
    #ROOT.gPad.Update()
    ##stats = histos[his.name].FindObject("stats")
    #stats = histos[his.name].FindObject("stats") 
    #if not stats: stats.__class__ = ROOT.TPaveStats
    #stats.SetX1NDC(0.5)
    #stats.SetX2NDC(0.5)

canvas.Print(plotFileName+"]") # Close plot file 

################################################################################
# Tidy up

inFile.Close()

outHistFile = ROOT.TFile.Open(outFileName ,"RECREATE")
outHistFile.cd()
for name in histos.keys():
    histos[name].Write() 
outHistFile.Close()

################################################################################
# Misc

#for entryNum in range(0,tree.GetEntries()): 
#    tree.GetEntry(entryNum)
#
#    # Histo
#    nTotTS="nTotTS"
#    histos[nTotTS].Fill(getattr(tree,nTotTS))
#
#    # Some filters
#    if getattr(tree,nTotTS) != 2: continue # Check if two tracksters
#
#    # More histos
#    for name,(htype,title,xbin,xmin,xmax) in config.items():
#        if htype=="TH1F":
#            if name!=nTotTS: 
#                histos[name].Fill(getattr(tree,title))
