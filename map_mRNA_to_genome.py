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
	mirwalk_temp_file_name = "_tmp_" + mirwalk_data_file.split('/')[-1]
	if mirwalk_temp_file_name in os.listdir():
		#Read this file and return the data.
		print("Found an already existing miRWalk _tmp file, reading that...")
		
		#Not creating a dictionary for id vs occurances, but have it just in case.
		nm_like_ids_list = []
		with open(mirwalk_temp_file_name, 'r') as f:
			for line in f:
				nm_like_ids_list.append(line.split('\t')[0])

		return nm_like_ids_list


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
		for nm_like_id, number_of_occurances in nm_like_ids.items():
			line = nm_like_id+'\t'+str(number_of_occurances)+'\n'
			f.write(line)

	return list(nm_like_ids.keys())

def make_dict_from_gencode_file(gencode_file_path):
	#Output dictionary.
	transcript_to_exonic_coordinates = {}

	#Check if a temporary file with the results of this function already exists.
	gencode_tmp_file_name = '_tmp_' + gencode_file_path.split('/')[-1]
	if gencode_tmp_file_name in os.listdir():
		#Read this file and return the data.
		print("Found an already existing gencode _tmp file, reading that...")

		with open(gencode_tmp_file_name, 'r') as f:
			for line in f:
				values = line.rstrip().split('\t')
				
				transcript_id = values[3]
				#Add stuff to dictionary.
				if transcript_id not in transcript_to_exonic_coordinates:
					transcript_to_exonic_coordinates[transcript_id] = []
				transcript_to_exonic_coordinates[transcript_id].append([*values[:3], values[4]])
				
		return transcript_to_exonic_coordinates
	
	#A dictionary that holds the information vs column index number.
	#Type could either be a gene or exon or CDS etc...
	gencode_indices = {'chr': 0, 'type': 2, 'start': 3, 'stop': 4, 'meta_data': 8}

	with open(gencode_file_path, 'r') as f, open(gencode_tmp_file_name, 'w') as w:
		for idx, line in enumerate(f):
			if line.startswith('#'):
				continue
			if idx % 1000 == 0:
				print(idx, " lines processed.")
				sys.stdout.write("\033[F")
			values = line.split('\t')
			
			#Check if the data seems to be intact.
			if len(values) == 9:

				#Check if type of entry is exon.
				if values[gencode_indices['type']] == 'exon':
					chromosome = values[gencode_indices['chr']]
					chromosome_start = values[gencode_indices['start']]
					chromosome_stop = values[gencode_indices['stop']]

					meta_data = values[gencode_indices['meta_data']]
					
					'''
					meta_data looks like:

					ID=exon:ENST00000456328.2:1;Parent=ENST00000456328.2;gene_id=ENSG00000223972.5;transcript_id=ENST00000456328.2;
					gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;transcript_type=lncRNA;transcript_name=DDX11L1-202;
					exon_number=1;exon_id=ENSE00002234944.1;level=2;transcript_support_level=1;hgnc_id=HGNC:37102;tag=basic;
					havana_gene=OTTHUMG00000000961.2;havana_transcript=OTTHUMT00000362751.1
					'''

					meta_data_of_exon = meta_data.split(';')[0] 
					
					#meta data of exon looks like: ID=exon:ENST00000456328.2:1

					transcript_id, exon_number = meta_data_of_exon.split(':')[1], meta_data_of_exon.split(':')[2]

					#Store information to a dictionary.
					if transcript_id not in transcript_to_exonic_coordinates:
						transcript_to_exonic_coordinates[transcript_id] = []
					transcript_to_exonic_coordinates[transcript_id].append([chromosome, chromosome_start, chromosome_stop, exon_number])

					#Also store it into the temporary file.
					write_line = [chromosome, chromosome_start, chromosome_stop, transcript_id, exon_number]
					w.write('\t'.join(write_line)+'\n')
	

	return transcript_to_exonic_coordinates



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

	print("Processing Gencode file, this might take a while as well...")
	gencode_annotations = make_dict_from_gencode_file(gencode_file_path)

	print(gencode_annotations)

	for nm_like_id in example_ids:
		ncbi_parsed_data = get_ncbi_data(example_database, example_ids[1])

		exons = get_exon_coordinates_from_ncbi(ncbi_parsed_data)

		#print(exons)
		break
