#!/bin/bash

#loop
cd /nfs/dust/cms/user/karavdia/TstarTstar/102X_v1/Preselection/
for mass in 700 800 900 1000 1100 1200 1300 1400 1500 1600
#for mass in 700
do
    for skip in skip0 skip1 skip2
    #for skip in skip1
    do
	cp RunII_2016_MuonHihjPtId_cutBasedPhotonIDlooseFall17_nonIsoandIsoHLT_tgtgRecoWithModifiedttbar_addMETandSTcuts_dr12_${skip}jets/tgtg/uhh2.AnalysisModuleRunner.MC.MC_TstarTstarToTgluonTgluon_M-${mass}_Run2016v3.root RunII_2016_MuonHihjPtId_cutBasedPhotonIDlooseFall17_nonIsoandIsoHLT_tgtgRecoWithModifiedttbar_addMETandSTcuts_dr12_skipjets_summary/tgtg/uhh2.AnalysisModuleRunner.MC.MC_TstarTstarToTgluonTgluon_M-${mass}_Run2016v3_${skip}.root
    done
done

cd /nfs/dust/cms/user/karavdia/CMSSW_10_2_10/src/UHH2/TstarTstar/macros/