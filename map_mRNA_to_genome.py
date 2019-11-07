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
import re


gencode_file = '../../data/gencode/gencode.v32.chr_patch_hapl_scaff.annotation.gff3'
mirwalk_data_file = '../../data/miRNA-mRNA-targets/miRWalk/hsa_miRWalk_3UTR.txt'
ncbi_annotation_file = '../../data/ncbi/ref_GRCh38.p12_top_level.gff3'

chromosome_identification_table = {}


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
			
			#Check if the a particular entry is intact.
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


def make_dict_from_ncbi_annotation(ncbi_annotation_file_path):
	#A fucntion to parse the NCBI annotation file and yield a dictionary for downstream intersection.
	#Functionality of gencode parser wasn't used becuase the metadata was organized differently for both the files.
	#Because of this, this function and the gencode parser function are very similar; redundancy.
	
	#Output dictionary.
	nm_like_id_to_exomic_coordinates = {}

	ncbi_tmp_annotation_file_name = '_tmp_' + ncbi_annotation_file_path.split('/')[-1]
	if ncbi_tmp_annotation_file_name in os.listdir():
		

		return nm_like_id_to_exomic_coordinates

	ncbi_annotation_indices = {'chr': 0, 'type': 2, 'start': 3, 'stop': 4, 'meta_data': 8}	
	exon_count = None


	#Open the annotation file and create a temporary file to write the parsed results.
	with open(ncbi_annotation_file_path, 'r') as f, open(ncbi_tmp_annotation_file_name, 'w') as w:
		for idx, line in enumerate(f):
			if line.startswith('#'):
				continue
			if idx % 1000 == 0:
				print(idx, " lines processed.")
				sys.stdout.write("\033[F")
			values = line.rstrip('\n').split('\t')
			'''
			An element of values (list) looks like:
			['NC_000001.11', 'BestRefSeq', 'mRNA', '19644173', '19658456', '.', '+', '.', 
			'ID=rna1988;Parent=gene589;Dbxref=GeneID:4681,Genbank:NM_005380.7,HGNC:HGNC:7650,MIM:600613;Name=NM_005380.7;
			gbkey=mRNA;gene=NBL1;product=neuroblastoma 1%2C DAN family BMP antagonist%2C transcript variant 2;transcript_id=NM_005380.7']
			'''		

			#Check if an entry is intact.
			if len(values) == 9:

				#Check if the entry is for an exon.
				if values[ncbi_annotation_indices['type']] == 'exon':
					#NCBI provides an NC_*** like id instead of chr**, so extracting that through another function. 
					chromosome = get_info_for_chromosomes(values[ncbi_annotation_indices['chr']])

					if not chromosome:
						continue

					chromosome_start = values[ncbi_annotation_indices['start']]
					chromosome_stop = values[ncbi_annotation_indices['stop']]

					meta_data = values[ncbi_annotation_indices['meta_data']]

					#Get all the meta data into a list for filtering.
					_transcript_id_entry = meta_data.split(';')[-1]
					transcript_id = _transcript_id_entry.split('=')[-1]

					#Select entries starting with NM_ only.
					if not transcript_id.startswith('NM'):
						continue
					
					#Store the entry into a dictionary.
					#See if this transcript_id doesn't already have an entry in the ..._exomic_coordinate dictionary.
					if transcript_id not in nm_like_id_to_exomic_coordinates:
						nm_like_id_to_exomic_coordinates[transcript_id] = []
					
					exon_count = nm_like_id_to_exomic_coordinates[transcript_id][-1][-1] + 1 if len(nm_like_id_to_exomic_coordinates[transcript_id]) > 0 else 1
					nm_like_id_to_exomic_coordinates[transcript_id].append([chromosome, chromosome_start, chromosome_stop, exon_count])

					#Also store the information into a file.
					write_line = [chromosome, chromosome_start, chromosome_stop, transcript_id, str(exon_count)]
					w.write('\t'.join(write_line)+'\n')

	return nm_like_id_to_exomic_coordinates


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


def get_info_for_chromosomes(nc_like_id):
	#If record is already there then just return that.
	if nc_like_id in chromosome_identification_table:
		return chromosome_identification_table[nc_like_id]

	#Record is not there, let's create it.
	#Fetching data from ncbi.
	#print("Hitting ncbi to get information for {}".format(nc_like_id))
	

	#Special case for X and Y.
	if nc_like_id == 'NC_000023.11':
		chromosome_identification_table[nc_like_id] = 'chrX'
	elif nc_like_id == 'NC_000024.10':
		chromosome_identification_table[nc_like_id] = 'chrY'
	
	ncbi_data = get_ncbi_data("nucleotide", nc_like_id)
	gen_bank_definition = ncbi_data['GBSeq_definition']

	try:
		chromosome_info = re.search("chromosome\ \d+", gen_bank_definition).span()

		#Extract the chromosome number from the definition.
		chromosome_number = gen_bank_definition[chromosome_info[0]:chromosome_info[1]].split(' ')[-1]

		chromosome_identification_table[nc_like_id] = 'chr' + chromosome_number

		return chromosome_identification_table[nc_like_id]
	except AttributeError:
		print("No information on chromosome found for {}".format(nc_like_id))
		chromosome_identification_table[nc_like_id] = False
		return False
	
	return chromosome_identification_table[nc_like_id]


if __name__ == "__main__":
	#Get data from miRWalk.
	example_database = "nucleotide"
	example_ids = ["NM_001126116", "NM_000546"]

	if '_tmp_nc_id_reference.txt' in os.listdir():
		print("Preparing chromosome identification table")
		with open('_tmp_nc_id_reference.txt', 'r') as f:
			for line in f:
				chromosome_identification_table[line.split('\t')[0].rstrip('\n')] = line.split('\t')[1].rstrip('\n')

	print("Processing miRWalk file, this might take a while...")
	#############miRWalk File Map#############
	#mirwalk_nm_like_ids = get_nm_ids_from_mirwalk_db(mirwalk_data_file)


	print("Processing Gencode file, this might take a while as well...")	
	#############Gencode File Map#############
	#gencode_annotations = make_dict_from_gencode_file(gencode_file)

	
	print("Processing NCBI annotation file, this might take a while as well...")
	#############NCBI annotation File Map#############
	#nm_like_id_to_exomic_coordinates = make_dict_from_ncbi_annotation(ncbi_annotation_file)

	#Dump the chromosome reference table.
	'''
	Commented out for precaution.
	if '_tmp_nc_id_reference.txt' not in os.listdir():
		with open('_tmp_nc_id_reference.txt', 'w') as f:
			for key, val in chromosome_identification_table.items():
				f.write(str(key) + '\t' + str(val) + '\n')
	'''
	print("Job done, Bye.")

