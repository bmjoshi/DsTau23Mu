import sys, os
from Filelist.filelist_2017_DoubleMuonLowMass_Run2017B_17Nov2017 import inputFiles
os.environ['TERM'] = 'vt100'
server_url="root://cms-xrd-global.cern.ch/"
#sys.argv[1]=dataest
#sys.argv[2]=year
dataset=''
year='2017'
cms_run_cert=''
if (year=='2016'): cms_run_cert=""
if (year=='2017'): cms_run_cert="Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.json"

count=0

command = "echo \"{}\" > "
os.system(command)

for dataFile in inputFiles:
	command=str("edmLumisInFiles.py "+server_url+dataFile+" > tmp_json.json")
	os.system(command)
	os.system("mergeJSON.py tmp_json.json doubleMuonLowMass_runII_2017b.json --output=doubleMuonLowMass_runII_2017b.json")
	print command
	count= count+1
	if (count>3): break

os.system("compareJSON.py --and doubleMuonLowMass_runII_2017b.json"+cms_run_cert+" --output=cert_doubleMuonLowMass_runII_2017b.json")

