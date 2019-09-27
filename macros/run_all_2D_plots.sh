#!/bin/bash
#loop over samples
for sample in MC\_TstarTstarToTgammaTgamma\_M-700\_Run2016v3 MC\_TstarTstarToTgammaTgamma\_M-800\_Run2016v3 MC\_TstarTstarToTgammaTgamma\_M-900\_Run2016v3 MC\_TstarTstarToTgammaTgamma\_M-1100\_Run2016v3 MC\_TstarTstarToTgammaTgluon\_M-1500\_Run2016v3 MC\_TstarTstarToTgammaTgluon\_M-1600\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-700\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-800\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-900\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-1100\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-1200\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-1300\_Run2016v3 MC\_TstarTstarToTgluonTgluon\_M-1500\_Run2016v3
do
    #loop over selection steps
    for step in After2D SemiLepTTBarMatch
    do
        #loop over hists
	for hist in pt_mu_pt_ak8jet1 pt_ele_pt_ak8jet1 pt_mu_pt_ak4jet1 pt_ele_pt_ak4jet1
	do
	    root -l -b -q Create2Dplots.C\(\"uhh2.AnalysisModuleRunner.MC.${sample}.root\",\"${sample}\",\"${step}\",\"${hist}\"\)
	done
    done
done