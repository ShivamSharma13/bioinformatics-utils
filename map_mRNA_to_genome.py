#!/usr/bin/env python3
"""
NCBI coordinates for a matured mRNA NM_xxx like entry are nothing more than an index.
Something like: 
	1 gtcgttggtg tgttgcgcga ctggccttga gggagagctg gggcctgctc ccggagagat
	61 acggctatgt cgatcgaaat cgaatcttcg gatgtgatcc gccttattat gcagtacttg
	121 aaggagaaca gtttacatcg ggcgttagcc accttgcagg aggagactac tgtgtctctg
On the other hand, popular miRNA binding site databases also use these coordinates instead of Genome Coordinates.
For analysis of these binding sites for our research, we need genomic coordinates instead of the coordinates mentioned here.
This script uses NCBI mRNAs and maps that mRNA's exon coordinates to the genome using the GENCODE annotation GFF3 file. 
"""
from hit_entrez import get_ncbi_data
import os
import sys


gencode_file_path = '../../data/gencode/gencode.v32.chr_patch_hapl_scaff.annotation.gff3'
mirwalk_data_file = '../../data/miRNA-mRNA-targets/miRWalk/hsa_miRWalk_3UTR.txt'

def get_nm_ids_from_mirwalk_db(mirwalk_data_file):
	#Check if a temporary file with the results of this function already exists.
	mirwalk_temp_file_name = mirwalk_data_file.split('/')[-1]
	if mirwalk_data_file.split('/')[-1] in os.listdir():
		#Read this file and return the data.
		return None

	#If the temporary file isn't there, then proces the root file and create one.
	nm_like_ids = {}
	with open(mirwalk_data_file, 'r') as f:
		for idx, line in enumerate(f):
			#Smart print funtion.
			if idx % 1000 == 0:
				print(idx, " lines processed.")
				sys.stdout.write("\033[F")
			
			#Splitting the line by tab characters to obatin just the NM_xxxx ID column values.
			values = line.split('\t')
			nm_like_id = values[1]

			#Checking if the file exists.
			if nm_like_id.startswith('NM_'):
				#Get or increment in nm_like_ids dictionary.
				nm_like_ids[nm_like_id] = nm_like_ids.get(nm_like_id, 1) + 1


	print("Writing the obtained values into a temporary text file.")
	with open(mirwalk_temp_file_name, 'w') as f:
		f.write('\n'.join(list(nm_like_ids.keys())))

	return nm_like_ids

def make_dict_from_gencode_file(gencode_file_path):
	pass

def get_exon_coordinates_from_ncbi(parsed_data):
	try:
		features = parsed_data["GBSeq_feature-table"]
	except KeyError:
		print("Feature list not present, returning...")
		return None
	
	exons = {}
	exon_count = 0

	for feature in features:
		if 'GBFeature_key' not in feature: 
			print("GBFeature_key absent from parsed data. Please check {}, skipping...".format(feature))
			continue
		if feature['GBFeature_key'] == 'exon':
			if "GBFeature_location" not in feature:
				print("GBFeature_location not present in exon, please check {}, skipping...".format(feature_exon))
			
			location_coordinates = feature['GBFeature_location']
			
			ncbi_start, ncbi_stop = location_coordinates.split('..')[0], location_coordinates.split('..')[1]

			#print(ncbi_start, ncbi_stop)
			exons['exon'+str(exon_count)] = {'start': ncbi_start, 'stop': ncbi_stop}
			exon_count+=1

	return exons


if __name__ == "__main__":
	#Get data from miRWalk.
	example_database = "nucleotide"
	example_ids = ["NM_001126116", "NM_000546"]
	
	print("Processing miRWalk file, this might take a while...")
	nm_like_ids = get_nm_ids_from_mirwalk_db(mirwalk_data_file)

	exit()
	print("Processing Gencode file, this might take a while as well...")
	gencode_annotations = make_dict_from_gencode_file(gencode_file_path)

	print(gen)

	for nm_like_id in example_ids:
		ncbi_parsed_data = get_ncbi_data(example_database, example_ids[1])

		exons = get_exon_coordinates_from_ncbi(ncbi_parsed_data)

		print(exons)
		break








