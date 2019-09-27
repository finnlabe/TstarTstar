#!/bin/bash
#loop over samples
denomstep=SemiLepTTBarMatchGENRECO_mu
for sample in MC\_TstarTstarToTgammaTgamma\_M-700\_Run2016v3 MC\_TstarTstarToTgammaTgamma\_M-800\_Run2016v3 MC\_TstarTstarToTgammaTgamma\_M-900\_Run2016v3 MC\_TstarTstarToTgammaTgamma\_M-1100\_Run2016v3 MC\_TstarTstarToTgammaTgluon\_M-1500\_Run2016v3 MC\_TstarTstarToTgammaTgluon\_M-1600\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-700\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-800\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-900\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-1100\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-1200\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-1300\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-1500\_Run2016v3
do
    #loop over triggers
    for step in SemiLepTTBarMatchGENRECO_triggerSingleLeptonMu SemiLepTTBarMatchGENRECO_triggerSingleJet_mu SemiLepTTBarMatchGENRECO_triggerHT_mu
    do
        #loop over hists
	for hist in Pt_mu Eta_mu AK4_Pt_b AK8_Pt_gluon AK8_Pt_top AK4_Eta_b AK8_Eta_gluon AK8_Eta_top
	do
	    root -l -b -q TriggerEffPlots.C\(\"uhh2.AnalysisModuleRunner.MC.${sample}.root\",\"${sample}\",\"${step}\",\"${hist}\",\"${denomstep}\"\)
	done
    done
done
