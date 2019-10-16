#!/bin/bash
#loop over samples
for sample in MC\_TstarTstarToTgammaTgamma\_M-700\_Run2016v3 MC\_TstarTstarToTgammaTgamma\_M-800\_Run2016v3 MC\_TstarTstarToTgammaTgamma\_M-900\_Run2016v3 MC\_TstarTstarToTgammaTgamma\_M-1100\_Run2016v3 MC\_TstarTstarToTgammaTgluon\_M-1500\_Run2016v3 MC\_TstarTstarToTgammaTgluon\_M-1600\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-700\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-800\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-900\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-1100\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-1200\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-1300\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-1500\_Run2016v3
do
    #loop over triggers
    for step in SemiLepTTBarMatchGENRECO_triggerSingleLeptonMu SemiLepTTBarMatchGENRECO_triggerSingleJet SemiLepTTBarMatchGENRECO_triggerSingleLeptonEle SemiLepTTBarMatchGENRECO_triggerHT
    do
        #loop over hists
	for hist in Pt_mu Eta_mu AK4_Pt_b AK8_Pt_gluon AK8_Pt_top
	do
	    root -l -b -q TriggerEffPlots.C\(\"uhh2.AnalysisModuleRunner.MC.${sample}.root\",\"${sample}\",\"${step}\",\"${hist}\"\)
	done
    done
done
