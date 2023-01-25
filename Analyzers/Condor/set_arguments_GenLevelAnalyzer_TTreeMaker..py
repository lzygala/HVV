import json
import os
import argparse
import sys
import decimal

inFile = open('../../Data/samples_allYears.json')

samples = json.load(inFile)

outputdir = ''
eospath = '/eos/cms/store/user/lzygala/HVV/GenLevelInfo/'
outputfile = ''
scale = 0.0

years = ["16APV","16","17","18"]

arguments = []

for year in years:
	for i in samples[year]:
		if i['type'] != "signal":
			continue
		outputdir = eospath + year + "/" + i['name']
		if not os.path.exists(outputdir):
			os.makedirs(outputdir)

		name = i['name']
		print(outputdir)
		scale = float(i['xs']) / float(i['nEvts'])
		counter = 0
		sampleType = i['type']
		for j in i['files']:
			outputfile = outputdir + "/" + name + "_" + str(counter) + ".root"	
			if not os.path.exists(outputfile):
				arguments.append("{} {} {}".format(name, outputfile, j))
			counter=counter+1 



print("Njobs: ", len(arguments))
    
with open("arguments_GenLevelAnalyzer_TTreeMaker.txt", "w") as args:
	args.write("\n".join(arguments))

inFile.close()
