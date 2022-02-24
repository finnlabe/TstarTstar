#!/usr/bin/env python3

# this very nice script was adapted from Christopher's version

import os
import sys
import csv
from collections import OrderedDict
import argparse
from itertools import permutations
from termcolor import colored

# to include CrossSectionHelper:
sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/common/UHH2-datasets'))
from CrossSectionHelper import MCSampleValuesHelper
helper = MCSampleValuesHelper()

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/TstarTstar/config'))
from constants import _YEARS


class configContainer:

   '''Container class for configurating XML files'''

   userName = str()
   uhh2Dir = str()
   userMail = str()
   yearVars = dict()
   used_samples = OrderedDict()


   def __init__(self):

      self.userName = os.environ.get('USER')
      self.uhh2Dir = os.environ.get('CMSSW_BASE')+'/src/UHH2/'
      self.userMail = '@'.join([self.userName, 'mail.desy.de']) # Avoid spam due to public code on GitHub

      self.yearVars['yearVersions'] = { # e.g. 'v3' in case  of 2016v3
         'UL16preVFP': '',
         'UL16postVFP': '',
         'UL17': '',
         'UL18': '',
      }

      # Set these values such that there are no more than 2,500 jobs per preselection. This way, you can submit two preselections in parallel to avoid going over 5,000 jobs (current user limit for NAF)
      self.yearVars['preselFileSplit'] = {
         'UL16preVFP': '50',
         'UL16postVFP': '50',
         'UL17': '50',
         'UL18': '50',
      }

      self.yearVars['targetLumis'] = {
         'UL16preVFP': _YEARS['UL16preVFP'].get('lumi_pb'),
         'UL16postVFP': _YEARS['UL16postVFP'].get('lumi_pb'),
         'UL17': _YEARS['UL17'].get('lumi_pb'),
         'UL18': _YEARS['UL18'].get('lumi_pb'),
      }

      self.yearVars['lumiFiles'] = {
         'UL16preVFP': self.uhh2Dir+'common/UHH2-data/UL16preVFP/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_UL16preVFP_normtag.root',
         'UL16postVFP': self.uhh2Dir+'common/UHH2-data/UL16postVFP/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_UL16postVFP_normtag.root',
         'UL17': self.uhh2Dir+'common/UHH2-data/UL17/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON_normtag.root',
         'UL18': self.uhh2Dir+'common/UHH2-data/UL18/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON_normtag.root',
      }

      self.yearVars['pileupFiles'] = {
         'mc': {
            'UL16preVFP': self.uhh2Dir+'common/UHH2-data/UL16preVFP/MyMCPileupHistogram_UL16preVFP.root',
            'UL16postVFP': self.uhh2Dir+'common/UHH2-data/UL16postVFP/MyMCPileupHistogram_UL16postVFP.root',
            'UL17': self.uhh2Dir+'common/UHH2-data/UL17/MyMCPileupHistogram_UL17.root',
            'UL18': self.uhh2Dir+'common/UHH2-data/UL18/MyMCPileupHistogram_UL18.root',
         },
         'data': {
            'UL16preVFP': self.uhh2Dir+'common/UHH2-data/UL16preVFP/MyDataPileupHistogram_UL16preVFP.root',
            'UL16postVFP': self.uhh2Dir+'common/UHH2-data/UL16postVFP/MyDataPileupHistogram_UL16postVFP.root',
            'UL17': self.uhh2Dir+'common/UHH2-data/UL17/MyDataPileupHistogram_UL17.root',
            'UL18': self.uhh2Dir+'common/UHH2-data/UL18/MyDataPileupHistogram_UL18.root',
         },
         'dataUp': {
            'UL16preVFP': self.uhh2Dir+'common/UHH2-data/UL16preVFP/MyMCPileupHistogram_UL16preVFP_72383.root',
            'UL16postVFP': self.uhh2Dir+'common/UHH2-data/UL16postVFP/MyMCPileupHistogram_UL16postVFP_72383.root',
            'UL17': self.uhh2Dir+'common/UHH2-data/UL17/MyDataPileupHistogram_UL17_72383.root',
            'UL18': self.uhh2Dir+'common/UHH2-data/UL18/MyDataPileupHistogram_UL18_72383.root',
         },
         'dataDown': {
            'UL16preVFP': self.uhh2Dir+'common/UHH2-data/UL16preVFP/MyDataPileupHistogram_UL16preVFP_66017.root',
            'UL16postVFP': self.uhh2Dir+'common/UHH2-data/UL16postVFP/MyDataPileupHistogram_UL16postVFP_66017.root',
            'UL17': self.uhh2Dir+'common/UHH2-data/UL17/MyDataPileupHistogram_UL17_66017.root',
            'UL18': self.uhh2Dir+'common/UHH2-data/UL18/MyDataPileupHistogram_UL18_66017.root',
         },
      }

      self.yearVars['HOTVRSFs'] = {
         'UL16preVFP': '',
         'UL16postVFP': '',
         'UL17': self.uhh2Dir+'HOTVR/data/TopTaggingScaleFactors_RunIISummer19UL17_PUPPIv15.root',
         'UL18': self.uhh2Dir+'HOTVR/data/TopTaggingScaleFactors_RunIISummer19UL18_PUPPIv15.root',
      }

      self.additionalBranches = {
         'Presselection': "",
         'Selection': "is_muevt evt_weight is_triggered is_highpt",
         'Analysis': "is_muevt evt_weight neutrino is_btagevent weight_sfmu_id weight_sfmu_id_down weight_sfmu_id_up weight_sfmu_isolation weight_sfmu_isolation_down weight_sfmu_isolation_up weight_sfelec_id weight_sfelec_id_down weight_sfelec_id_up TopTagSF TopTagSF_down TopTagSF_merged_down TopTagSF_merged_up TopTagSF_non_down TopTagSF_non_up TopTagSF_semi_up TopTagSF_semi_down TopTagSF_up",
         'DNN': "is_btagevent is_muevt evt_weight ST_weight DNN_Inputs neutrino TstarTstar_Hyp_gHOTVR TstarTstar_Hyp_gAK4 ST weight_sfmu_id weight_sfmu_id_down weight_sfmu_id_up weight_sfmu_isolation weight_sfmu_isolation_down weight_sfmu_isolation_up weight_sfelec_id weight_sfelec_id_down weight_sfelec_id_up TopTagSF TopTagSF_down TopTagSF_merged_down TopTagSF_merged_up TopTagSF_non_down TopTagSF_non_up TopTagSF_semi_up TopTagSF_semi_down TopTagSF_up"
      }

      self.systematics = list()

   @staticmethod
   def read_database_106X_v2(years: list, groups: list):

      from sample_db import samplesDict
      used_samples = OrderedDict()
      for year in years:
          temp = OrderedDict()
          for group in groups:
             temp[group] = list()
             for k, v in samplesDict.items():
                use_me_year = v.get('years') == None or year in v.get('years', [])

                # block to group config files
                use_me_group = False
                if k.startswith("DATA"):
                    if (group == "DATA"): use_me_group = True
                elif k.startswith("TstarTstar"):
                    if group == "Signal": use_me_group = True
                else:
                    if group == "SM": use_me_group = True
                # end grouping block

                if (use_me_year and use_me_group):

                   sample_entity = sampleEntity((k, v, year,))
                   if not os.path.isfile(sample_entity.xmlPath):
                      # mabe this was split manually?
                      print(colored('XML for sample  '+sample_entity.nickName+' ('+year+')  does not exist. Skipping this sample', 'red'))
                      continue
                   temp[group].append(sample_entity)
          used_samples[year] = temp
      configContainer.used_samples = used_samples

class sampleEntity:

   '''Container to hold information about a data or MC sample, as read from CSV database'''

   def __init__(self, csvRow):

      k, v, year = csvRow
      self.year = year
      self.is_data = k.startswith('DATA_')
      self.nickName = k
      self.n_das = None
      self.n_pnfs = None
      self.pnfs_sum_of_weights = helper.get_nevt(v['db_name'], '13TeV', self.year)
      self.xs_gen = None
      self.xs_theo = None
      self.lumi = 1. if self.is_data else helper.get_lumi(v['db_name'], '13TeV', self.year, kFactor=v.get('kfac', False), Corrections=v.get('corr', False))
      self.xsection = self.pnfs_sum_of_weights / self.lumi
      self.xmlPath = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/common/UHH2-datasets', helper.get_xml(v['db_name'], '13TeV', self.year))


class xmlCreator:

   '''Creates XML files for SFrame'''

   def __init__(self, step: str, year: str, group: str):

      confCon = configContainer()
      self.uhh2Dir = confCon.uhh2Dir
      self.userName = confCon.userName
      self.userMail = confCon.userMail
      self.yearVars = confCon.yearVars
      self.sample_list = confCon.used_samples[year][group]
      self.systematics = confCon.systematics
      self.additionalBranches = confCon.additionalBranches

      if step not in ['Preselection', 'Selection', 'Analysis', 'DNN']:
         sys.exit('Given value of argument "selection" not valid. Abort.')
      self.step = step
      self.is_presel = True if step=='Preselection' else False
      if (step == "Preselection"):
          self.previousFolder = "none"
          self.analysisModule = "TstarTstarPreselectionModule"
      elif (step == "Selection"):
          self.previousFolder = "Preselection"
          self.analysisModule = "TstarTstarSelectionModule"
      elif (step == "Analysis"):
          self.previousFolder = "Selection"
          self.analysisModule = "TstarTstarAnalysisModule"
      elif (step == "DNN"):
          self.previousFolder = "Analysis"
          self.analysisModule = "TstarTstarDNNModule"

      if year not in ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']:
         sys.exit('Given value of argument "year" not valid. Abort.')
      self.year = year
      self.yearVersion = self.yearVars['yearVersions'][year]

      self.outputDirBase = '/nfs/dust/cms/user/flabe/TstarTstar/data/'
      if not os.path.isdir(self.outputDirBase):
         sys.exit('Warning: Make sure to create output directory via "ln -s". Abort.')

      self.xmlFileName = '_'.join(['parsedConfigFile', self.step, self.year, group])+'.xml'
      self.xmlFilePathBase = self.uhh2Dir+'TstarTstar/config/'+'_'.join(['config', self.step, self.year])+'/'
      os.makedirs(self.xmlFilePathBase, exist_ok=True)
      self.xmlFilePath = self.xmlFilePathBase+self.xmlFileName
      self.workdirName = '_'.join(['workdir', self.step, self.year, group])

      self.write_xml_successful = False


   def write_xml(self):

      with open(self.xmlFilePath, 'w') as file:
         file.write('''<?xml version="1.0" encoding="UTF-8"?>\n''')
         file.write('''\n''')
         file.write('''<!DOCTYPE JobConfiguration PUBLIC "" "JobConfig.dtd"[\n''')
         file.write('''\n''')
         file.write('''<!ENTITY TargetLumi "'''+str(self.yearVars['targetLumis'][self.year])+'''">\n''')
         if not self.is_presel:
            file.write('''<!ENTITY INPUTdir "'''+(self.outputDirBase+self.previousFolder+'/'+self.year)+'''">\n''')
            file.write('''<!ENTITY INPUTfilename "uhh2.AnalysisModuleRunner">\n''')
         file.write('''<!ENTITY OUTPUTdir "'''+(self.outputDirBase+self.step+'/'+self.year)+'''">\n''')
         file.write('''<!ENTITY b_Cacheable "False">\n''')
         file.write('''<!ENTITY NEVT "-1">\n''')
         file.write('''<!ENTITY YEARsuffix "_'''+self.year+self.yearVersion+'''">\n''')
         file.write('''<!ENTITY PROOFdir "/nfs/dust/cms/user/'''+self.userName+'''/.proof2">\n''')
         file.write('''\n''')
         for s in self.sample_list:
            if not self.is_presel:
               file.write('''<!ENTITY '''+s.nickName+''' "&INPUTdir;/&INPUTfilename;'''+('.DATA.' if s.is_data else '.MC.')+s.nickName+'''&YEARsuffix;.root">\n''')
            else:
               file.write('''<!ENTITY '''+s.nickName+''' SYSTEM "'''+s.xmlPath+'''">\n''')
         file.write('''\n''')
         file.write(''']>\n''')
         file.write('''\n''')
         file.write('''<!--\n''')
         file.write('''<ConfigParse NEventsBreak="'''+('500000' if not self.is_presel else '0')+'''" FileSplit="'''+('0' if not self.is_presel else self.yearVars['preselFileSplit'][self.year])+'''" AutoResubmit="5"/>\n''')
         file.write('''<ConfigSGE RAM="4" DISK="3" Mail="'''+self.userMail+'''" Notification="as" Workdir="'''+self.workdirName+'''"/>\n''')
         file.write('''-->\n''')
         file.write('''\n''')
         file.write('''<!-- OutputLevel controls which messages are printed; set to VERBOSE or DEBUG for more verbosity, to WARNING or ERROR for less -->\n''')
         file.write('''<JobConfiguration JobName="ExampleCycleJob" OutputLevel="INFO">\n''')
         file.write('''<Library Name="libSUHH2TstarTstar"/>\n''')
         file.write('''<Package Name="SUHH2TstarTstar.par"/>\n''')
         file.write('''<Cycle Name="uhh2::AnalysisModuleRunner" OutputDirectory="&OUTPUTdir;/" PostFix="" TargetLumi="&TargetLumi;">\n''')
         file.write('''\n''')
         for s in self.sample_list:
             if self.is_presel:
                 file.write('''<InputData Lumi="'''+str(s.lumi)+'''" NEventsMax="&NEVT;" Type="'''+('DATA' if s.is_data else 'MC')+'''" Version="'''+s.nickName+'''&YEARsuffix;" Cacheable="&b_Cacheable;"> &'''+s.nickName+'''; <InputTree Name="AnalysisTree"/> <OutputTree Name="AnalysisTree"/> </InputData>\n''')
             else:
                 file.write('''<InputData Lumi="'''+str(s.lumi)+'''" NEventsMax="&NEVT;" Type="'''+('DATA' if s.is_data else 'MC')+'''" Version="'''+s.nickName+'''&YEARsuffix;" Cacheable="&b_Cacheable;"> <In FileName="&'''+s.nickName+''';" Lumi="0.0"/> <InputTree Name="AnalysisTree"/> <OutputTree Name="AnalysisTree"/> </InputData>\n''')
         file.write('''\n''')
         file.write('''<UserConfig>\n''')
         file.write('''\n''')
         file.write('''<Item Name="PrimaryVertexCollection" Value="offlineSlimmedPrimaryVertices"/>\n''')
         file.write('''<Item Name="METName" Value="slimmedMETs"/>\n''')
         file.write('''<Item Name="ElectronCollection" Value="slimmedElectronsUSER"/>\n''')
         file.write('''<Item Name="PhotonCollection" Value="slimmedPhotonsUSER"/>\n''')
         file.write('''<Item Name="MuonCollection" Value="slimmedMuonsUSER"/>\n''')
         file.write('''<Item Name="JetCollection" Value="jetsAk4Puppi"/>\n''')
         file.write('''<Item Name="GenJetCollection" Value="slimmedGenJets"/>\n''')
         file.write('''<Item Name="TopJetCollection" Value="hotvrPuppi"/>\n''')
         file.write('''<Item Name="GenTopJetCollection" Value="hotvrGen"/>\n''')
         file.write('''<Item Name="GenParticleCollection" Value="GenParticles"/>\n''')
         file.write('''<Item Name="GenInfoName" Value="genInfo"/>\n''')
         file.write('''\n''')
         if not (self.additionalBranches[self.step] == ""):
             file.write('''<Item Name="additionalBranches" Value="'''+self.additionalBranches[self.step]+'''"/>\n''')
         file.write('''\n''')
         file.write('''<Item Name="lumi_file" Value="'''+self.yearVars['lumiFiles'][self.year]+'''"/>\n''')
         file.write('''<Item Name="lumihists_lumi_per_bin" Value="500."/>\n''')
         file.write('''\n''')
         file.write('''<Item Name="pileup_directory" Value="'''+self.yearVars['pileupFiles']['mc'][self.year]+'''"/>\n''')
         file.write('''<Item Name="pileup_directory_data" Value="'''+self.yearVars['pileupFiles']['data'][self.year]+'''"/>\n''')
         file.write('''\n''')
         file.write('''<!-- Tell AnalysisModuleRunner NOT to use the MC event weight from SFrame; rather let MCLumiWeight (called via CommonModules) calculate the MC event weight. The MC event weight assigned by MCLumiWeight is InputData.Lumi / Cycle.TargetLumi. -->\n''')
         file.write('''<Item Name="use_sframe_weight" Value="false"/>\n''')
         file.write('''<Item Name="AnalysisModule" Value="'''+self.analysisModule+'''"/>\n''')
         file.write('''<Item Name="uhh2Dir" Value="'''+self.uhh2Dir+'''"/>\n''')
         if(self.step == "Selection"):
             file.write('''\n''')
             file.write('''<!-- scale factor configuration -->\n''')
             file.write('''<Item Name="HOTVRTopTagSFs" Value="'''+self.yearVars['HOTVRSFs'][self.year]+'''"/>\n''')
             file.write('''<Item Name="SF_path" Value="/nfs/dust/cms/user/flabe/TstarTstar/CMSSW_10_2_17/src/UHH2/TstarTstar/factors/" />\n''')
             file.write('''<Item Name="BTagCalibration" Value="/nfs/dust/cms/user/flabe/TstarTstar/CMSSW_10_2_17/src/UHH2/TstarTstar/factors/btag/reshaping_deepJet_106XUL18_v2.csv" />\n''')
             file.write('''<Item Name="NLOCorrections" Value = "/nfs/dust/cms/user/deleokse/RunII_102X_v2/CMSSW_10_2_17/src/UHH2/ZprimeSemiLeptonic/data/" />''')
         file.write('''\n''')
         file.write('''<!-- Switch for debugging of the central AnalysisModule -->\n''')
         file.write('''<Item Name="debug" Value="false"/>\n''')
         file.write('''\n''')
         file.write('''</UserConfig>\n''')
         file.write('''\n''')
         file.write('''</Cycle>\n''')
         file.write('''</JobConfiguration>\n''')

      self.write_xml_successful = True

      print('Created '+self.xmlFilePath)

      return self.xmlFilePath

   def write_bash_scripts(self):

       scriptFilePath_sframe_batch = os.path.join(self.xmlFilePathBase, 'run_all_sframe_batch.sh')
       with open(scriptFilePath_sframe_batch, 'w') as outfile:
          outfile.write('#!/bin/bash\n')
          newline_base = 'sframe_batch.py $1 '
          outfile.write(newline_base+self.xmlFilePath+'\n')
          for systXmlFilePath in self.systXmlFilePaths:
             outfile.write(newline_base+systXmlFilePath+'\n')
       print('Created '+scriptFilePath_sframe_batch)

       scriptFilePath_run_local = os.path.join(self.xmlFilePathBase, 'run_all_local.sh')
       with open(scriptFilePath_run_local, 'w') as outfile:
          outfile.write('#!/bin/bash\n')
          newline_base = 'python run_local.py $1 '
          outfile.write(newline_base+self.workdirName+'\n')
          for systWorkdirName in self.systWorkdirNames:
             outfile.write(newline_base+systWorkdirName+'\n')
       print('Created '+scriptFilePath_run_local)
       # copy_run_local_command = 'cp -n '+os.path.join(self.xmlFilePathBase, '..', 'run_local.py')+' '+os.path.join(self.xmlFilePathBase, '.')
       link_run_local_command = 'ln -s '+os.path.abspath(os.path.join(self.xmlFilePathBase, '..', 'run_local.py'))+' '+os.path.join(self.xmlFilePathBase, 'run_local.py')
       os.system(link_run_local_command)



if __name__=='__main__':

   selections = ['Preselection', 'Selection', 'Analysis', 'DNN']
   years = ['UL16preVFP','UL16postVFP','UL17', 'UL18']

   if not sys.argv[1:]: sys.exit('No arguments provided. Exit.')
   parser = argparse.ArgumentParser()
   parser.add_argument('--all', action='store_true', help='Create XML files for all selections and years.')
   parser.add_argument('--syst', action='store_true', help='Create XML files for systematic uncertainties. Will also create bash scripts with lists of SFrameBatch/run_local.py commands.')
   parser.add_argument('-s', '--selections', choices=selections, nargs='*', default=[])
   parser.add_argument('-y', '--years', choices=years, nargs='*', default=[])
   parser.add_argument('-a', '--auto-complete', action='store_true', help='Auto-complete arguments if not all arguments for selections and years are given.')
   args = parser.parse_args(sys.argv[1:])

   if(args.all == True):
      if(len(args.selections) + len(args.years) != 0):
         sys.exit('Not allowed to use "--all" option jointly with manually given selection or year argument. Exit.')
      args.selections = selections
      args.years = years
      if(args.auto_complete):
         print('Warning: You already specified "--all". Therefore, "--auto-complete" will not have any effect.')
   else:
      if args.auto_complete:
         print('Auto-completing arguments.')
         if not args.selections: args.selections = selections
         if not args.years: args.years = years
      else:
         for p in permutations([args.selections, args.years]):
            if p[0] and not p[1]:
               sys.exit('You specified arguments for at least one of the two options: "--selections" or "--years", but not for both of them. Also, you did not specify "--auto-complete" to compensate for this. Exit.')

   print('Going to create XML files for:')
   print('  Selections: '+', '.join(str(x) for x in args.selections))
   print('  Years: '+', '.join(str(x) for x in args.years))

   groups = ['DATA', 'Signal', 'SM']
   configContainer.read_database_106X_v2(args.years, groups)

   for selection in args.selections:
      for year in args.years:
          for group in groups:
             x = xmlCreator(selection, year, group)
             x.write_xml()
