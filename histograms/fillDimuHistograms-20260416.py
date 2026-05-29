#!/usr/bin/python3

import subprocess, optparse
from pathlib import Path
import ROOT, math

ETA_MASS = 0.547862
MU_MASS = 0.105658


def parseOptions():
	usage = ('usage: %prog [options]\n' + '%prog -h for help')
	parser = optparse.OptionParser(usage)
	parser.add_option('-i', '--input', dest='INPUT', type='string', help='input file', default='')
	parser.add_option('-o', '--output', dest='OUTPUT', type='string', help='output file')
	parser.add_option('-n', '--njobs', dest='NJOBS', type=int, help='total njobs')
	parser.add_option('-j', '--job', dest='JOB', type=int, help='job index')
	parser.add_option('-l', dest='LIST', help='input file list', default='/afs/cern.ch/user/j/jfriesen/CMSSW_13_0_10/src/Run3DimuonAnalysisTools/Plotting/FillHistogram/muMuGammaTree_ntuples_22_23_24.txt')
	global opt, args
	(opt, args) = parser.parse_args()


def processCmd(cmd, quite=0):
	status, output = subprocess.getstatusoutput(cmd)
	if status != 0 and not quite:
		print('Error in processing command:\n \t[' + cmd + ']')
		print('Output:\n   [' + output + '] \n')
	return output


pt_thresh_bb = 7.36
pt_thresh_ee = 7.06


def map_path_to_category(path, global_index):
	"""Map an input file path to a category number.

	- 2022 -> 1
	- 2023 -> 2
	- 2024 -> 3 or 4 split by parity of global_index (even->3, odd->4)
	- fallback -> 0
	"""
	if path is None:
		return 0
	# Try common capitalized markers first
	if "Run2022" in path:
		return 1
	if "Run2023" in path:
		return 2
	if "Run2024" in path:
		return 3 if (global_index % 2 == 0) else 4

	# fallback to lowercase search
	pl = path.lower()
	if "run2022" in pl:
		return 1
	if "run2023" in pl:
		return 2
	if "run2024" in pl:
		return 3 if (global_index % 2 == 0) else 4

	return 0

def fillHistograms():
	verbose = False
	ROOT.gROOT.SetBatch()

	global opt, args
	parseOptions()

	print('Opening list of input files', opt.LIST)
	file_list = Path(opt.LIST).read_text().splitlines()
	print('\tFound', len(file_list), 'items in list')

	njobs_actual = min(opt.NJOBS, len(file_list))
	file_range = (
		int(((float(opt.JOB) - 1) * len(file_list)) / njobs_actual),
		int(((float(opt.JOB)) * len(file_list)) / njobs_actual)
	)

	job_file_indices = list(range(file_range[0], file_range[1]))

	print(f'\tJob {opt.JOB} of {njobs_actual} has files ({file_range[0]}, {file_range[1]}] with {len(job_file_indices)} files')

	redirector = 'root://cmsxrootd.fnal.gov//'
	redirector = 'root://xrootd.cmsaf.mit.edu:1094//'
	print('\tUsing redirector', redirector)

	itree_name = 'tree/tree'
	itree = ROOT.TChain(itree_name)
	# slice of the global file list that this job will process
	job_files = file_list[file_range[0]:file_range[1]]
	for i, f in enumerate(job_files, start=file_range[0]):
		print('Getting', itree_name, 'from', redirector + f)
		itree.Add(redirector + f)
		print(itree.GetEntries(), 'total entries in TChain')

	print('Creating ' + str(opt.OUTPUT) + str(opt.JOB) + '.root')
	outfile = ROOT.TFile(str(opt.OUTPUT) + str(opt.JOB) + '.root', 'recreate')

	bin_width = 0.0001
	mumu_mass_range = (round(0.0 / bin_width) * bin_width, round(0.6 / bin_width) * bin_width)
	mumu_mass_bins = round((mumu_mass_range[1] - mumu_mass_range[0]) / bin_width)

	histos = {}

	print('Processing events...')
	i_event = 0
	for ev in itree:
		# determine which file in the global list this event came from
		tree_num = itree.GetTreeNumber()
		global_index = file_range[0] + int(tree_num)
		file_path = file_list[global_index]
		cat = map_path_to_category(file_path, global_index)

		if ev.mass > mumu_mass_range[1]:
			continue
		if verbose or i_event % 1000 == 0:
			print('mumu_mass <', mumu_mass_range[1], 'event', i_event)
		i_event += 1

		abs_eta1 = abs(ev.eta1)
		abs_eta2 = abs(ev.eta2)

		is_bb = (abs_eta1 < 1.48) and (abs_eta2 < 1.48)
		is_ee = (abs_eta1 >= 1.48 and abs_eta1 < 3.0) and (abs_eta2 >= 1.48 and abs_eta2 < 3.0)

		if is_bb:
			if ev.pt < pt_thresh_bb:
				continue
			flavor = 'bb'
		elif is_ee:
			if ev.pt < pt_thresh_ee:
				continue
			flavor = 'ee'
		else:
			continue

		if ev.custom_softMvaRun3Value1 > 0.6 and ev.custom_softMvaRun3Value2 > 0.6:
			srcr = '_SR'
		else:
			srcr = '_CR'

		histo_name = f'Cat_{cat}_{flavor}{srcr}'

		if histo_name not in histos:
			histos[histo_name] = ROOT.TH1F(histo_name, histo_name, mumu_mass_bins, mumu_mass_range[0], mumu_mass_range[1])

		histos[histo_name].Fill(ev.mass)

	print('Saving histograms...')
	for pt_key in histos:
		histo = histos[pt_key]
		name = histo.GetName()
		if histo.GetEntries() > 0:
			print('\tSaving', name, 'to', str(opt.OUTPUT) + str(opt.JOB) + '.root', 'with', histo.GetEntries(), 'entries')
			outfile.WriteObject(histo, name)
		else:
			print('\tNOT saving', name)

	outfile.Close()


if __name__ == '__main__':
	fillHistograms()
