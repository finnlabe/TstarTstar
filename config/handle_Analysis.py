import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--option', required=True)
parser.add_argument('--doData', action='store_true')
args = parser.parse_args()

years = ["UL16preVFP", "UL16postVFP", "UL17", "UL18"]

for year in years:
    os.system("sframe_batch.py -" + args.option + " config_Analysis_" + year + "/parsedConfigFile_Analysis_" + year + "_SM.xml")
    os.system("sframe_batch.py -" + args.option + "  config_Analysis_" + year + "/parsedConfigFile_Analysis_" + year + "_Signal.xml")
    os.system("sframe_batch.py -" + args.option + "  config_Analysis_" + year + "/parsedConfigFile_Analysis_" + year + "_Signal_Spin32.xml")
    if(args.doData): os.system("sframe_batch.py -" + args.option + "  config_Analysis_" + year + "/parsedConfigFile_Analysis_" + year + "_DATA.xml")
