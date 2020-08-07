#!/usr/bin/env python3
import argparse
import re

# Your code here
def make_batches_from_fold_value(sequence, fold):
	'''
	This function returns a list of batches based on fold value.
	'''
	if len(sequence) <= fold:
		return [sequence]
	else:
		return [sequence[i:i+fold] for i in range(0, len(sequence), fold)]


def get_file_format_based_on_sequence(entries):
	'''
	Takes the data structure entries and returns a file extension based sequence.
	'''
	nucleotide_regex = re.compile('^[ACGTacgt]+')

	#A, C, and G has been removed from amino acid regex.
	amino_acid_regex = re.compile('^[SPVILNDKQEMHFRYWspvilndkqemhfryw]+')

	#We'll merge the whole sequence and see if there's anything but ACGT, then it's probably an amino acid sequence.
	merged_sequence_for_testing = ''.join([i[1] for i in entries])
	
	regex_result = nucleotide_regex.match(merged_sequence_for_testing)


	if regex_result is not None:
		if regex_result.span()[1] == len(merged_sequence_for_testing):
			file_extension = '.fna'
		else:
			file_extension = '.faa'
	else:
		#Else becuase it's already been matched.
		file_extension = '.faa'

	return file_extension


def embl_to_fasta(embl_file_path, fold):
	entries = []
	sequence = ''
	sequence_id = None
	sequence_flag = False
	nucleotide_regex = re.compile('^[ACGTacgt]+')

	#A, C, and G has been removed from amino acid regex.
	amino_acid_regex = re.compile('^[SPVILNDKQEMHFRYWspvilndkqemhfryw]+')

	ends_with_number_regex = regex = re.compile('[0-9]+')
	with open(embl_file_path, 'r') as f:
		for line in f:
			line = line.rstrip('\n')
			if line == '//':
				entries.append([sequence_id, sequence])
				sequence_flag = False
				sequence = ''
			elif sequence_flag:
				#This means that the coming lines are all part of a sequence.
				if_numbers_at_end = ends_with_number_regex.search(line)
				if if_numbers_at_end is not None:
					cut_start = if_numbers_at_end.span()[0]
					cut_stop = if_numbers_at_end.span()[1]
					line = line[:cut_start].strip()
					line = line.replace(' ', '')
					sequence += line.upper()
			elif line.startswith('ID'):
				id_line_values = line.split(';')
				sequence_id = id_line_values[0].lstrip('ID').lstrip()
			elif line.startswith('DE'):
				sequence_id += " " + line.lstrip('DE').lstrip()
			elif line.startswith('SQ'):
				sequence_flag = True

	file_extension = get_file_format_based_on_sequence(entries)

	file_name = embl_file_path.split('/')[-1]
	if file_name.endswith('.txt'):
		new_file_name = file_name.replace('.txt', file_extension)
	else:
		new_file_name = file_name.split('.')[0] + file_extension
	output_file_path = embl_file_path.rstrip(embl_file_path.split('/')[-1]) + new_file_name

	#Will write to a file now.
	with open(output_file_path, 'w') as w:
		for entry in entries:
			fasta_id = entry[0]
			fasta_sequence = entry[1]
			w.write('>' + 'ENA|' + fasta_id + '\n')
			w.write('\n'.join(make_batches_from_fold_value(fasta_sequence, fold)) + '\n')

	return True


def fastq_to_fasta(fastq_file_path, fold):
	entries = []
	sequence = ''
	sequence_id = None
	nucleotide_regex = re.compile('^[ACGTacgt]+')

	#A, C, and G has been removed from amino acid regex.
	amino_acid_regex = re.compile('^[SPVILNDKQEMHFRYWspvilndkqemhfryw]+')
	with open(fastq_file_path, 'r') as f:
		for line in f:
			line = line.rstrip('\n')
			if line.startswith('@'):
				#New line is starting.
				sequence_id = line.lstrip('@')
				sequence = ''
			elif line.startswith('+'):
				#Sequence ended.
				entries.append([sequence_id, sequence])
			elif nucleotide_regex.match(line) or amino_acid_regex.match(line):
				sequence += line
	
	file_extension = get_file_format_based_on_sequence(entries)

	file_name = fastq_file_path.split('/')[-1]
	if file_name.endswith('.txt'):
		new_file_name = file_name.replace('.txt', file_extension)
	else:
		new_file_name = file_name.split('.')[0] + file_extension
	output_file_path = fastq_file_path.rstrip(fastq_file_path.split('/')[-1]) + new_file_name

	#Will write to a file now.
	with open(output_file_path, 'w') as w:
		for entry in entries:
			fasta_id = entry[0]
			fasta_sequence = entry[1]
			w.write('>' + fasta_id + '\n')
			w.write('\n'.join(make_batches_from_fold_value(fasta_sequence, fold)) + '\n')

	return True


def genbank_to_fasta(gb_file_path, fold):
	entries = []
	sequence = ''
	sequence_id = None
	sequence_flag = False
	nucleotide_regex = re.compile('^[ACGTacgt]+')

	#A, C, and G has been removed from amino acid regex.
	amino_acid_regex = re.compile('^[SPVILNDKQEMHFRYWspvilndkqemhfryw]+')

	ends_with_number_regex = regex = re.compile('[0-9]+')
	with open(gb_file_path, 'r') as f:
		for line in f:
			line = line.rstrip('\n')
			if line == '//':
				entries.append([sequence_id, sequence])
				sequence_flag = False
				sequence = ''
			elif sequence_flag:
				#This means that the coming lines are all part of a sequence.
				if_numbers_at_end = ends_with_number_regex.search(line)
				if if_numbers_at_end is not None:
					cut_start = if_numbers_at_end.span()[0]
					cut_stop = if_numbers_at_end.span()[1]
					line = line[cut_stop:].strip()
					line = line.replace(' ', '')
					sequence += line.upper()
			elif line.startswith('LOCUS'):
				line = line.lstrip('LOCUS').lstrip()
				sequence_id = line.split(' ')[0]
			elif line.startswith('DEFINITION'):
				sequence_id += " " + line.lstrip('DEFINITION').lstrip()
			elif line.startswith('ORIGIN'):
				sequence_flag = True

	file_extension = get_file_format_based_on_sequence(entries)

	file_name = gb_file_path.split('/')[-1]
	if file_name.endswith('.txt'):
		new_file_name = file_name.replace('.txt', file_extension)
	else:
		new_file_name = file_name.split('.')[0] + file_extension
	output_file_path = gb_file_path.rstrip(gb_file_path.split('/')[-1]) + new_file_name
	#Will write to a file now.
	with open(output_file_path, 'w') as w:
		for entry in entries:
			fasta_id = entry[0]
			fasta_sequence = entry[1]
			w.write('>' + 'ENA|' + fasta_id + '\n')
			w.write('\n'.join(make_batches_from_fold_value(fasta_sequence, fold)) + '\n')

	return True


def mega_to_fasta(mega_file_path, fold):
	entries_dict = {}
	entries = []
	interleaved_id_entries = []
	sequence_flag = False
	sequence_id = ''
	sequence = ''
	nucleotide_regex = re.compile('^[ACGTacgt]+')

	#A, C, and G has been removed from amino acid regex.
	amino_acid_regex = re.compile('^[SPVILNDKQEMHFRYWspvilndkqemhfryw]+')

	ends_with_number_regex = regex = re.compile('[0-9]+')

	with open(mega_file_path, 'r') as f:
		for line in f:
			line = line.rstrip('\n')
			if line.startswith('#MEGA') or line.startswith('#mega'):
				continue
			elif line.startswith('#'):
				if ' ' in line:					
					if line.lstrip('#') not in entries_dict:
						entries_dict[line.lstrip('#').split(' ')[0]] = line.lstrip('#').split(' ')[1].upper()
					else:
						entries_dict[line.lstrip('#').split(' ')[0]] += line.lstrip('#').split(' ')[1]/upper()
					interleaved_id_entries.append(line.lstrip('#').split(' ')[0])
				else:
					sequence_id = line.lstrip('#')
					sequence_flag = True
			elif sequence_flag:
				sequence += line.upper()

	if len(entries_dict) > 0:
		#An interleaved file was there.
		for interleaved_id_entry in interleaved_id_entries:
			sequence = entries_dict[interleaved_id_entry]
			entries.append([interleaved_id_entry, sequence])

	else:
		entries.append([sequence_id, sequence])

	file_extension = get_file_format_based_on_sequence(entries)

	file_name = mega_file_path.split('/')[-1]
	if file_name.endswith('.txt'):
		new_file_name = file_name.replace('.txt', file_extension)
	else:
		new_file_name = file_name.split('.')[0] + file_extension
	output_file_path = mega_file_path.rstrip(mega_file_path.split('/')[-1]) + new_file_name

	#Will write to a file now.
	with open(output_file_path, 'w') as w:
		for entry in entries:
			fasta_id = entry[0]
			fasta_sequence = entry[1]
			w.write('>' + fasta_id + '\n')
			w.write('\n'.join(make_batches_from_fold_value(fasta_sequence, fold)) + '\n')

	return True



def get_format_and_pass_to_converters(file_path, fold):
	#Get the format of the input file.
	with open(file_path, 'r') as f:
		for line in f:
			#The loop shall not parse through the whole file.
			#Not using readline() because: maybe the file starts from second line and has nothing in the first line.
			if line.startswith('@'):
				#It is a fastq file.
				fastq_to_fasta(file_path, fold)
				break
			elif line.startswith('ID'):
				#It is an embl file.
				embl_to_fasta(file_path, fold)
				break
			elif line.startswith('#MEGA') or line.startswith('#mega'):
				#It is a Mega file.
				mega_to_fasta(file_path, fold)
				break
			elif line.startswith('LOCUS') or line.startswith('locus'):
				#It is a genbank file.
				genbank_to_fasta(file_path, fold)
				break
			elif line == '\n' or line == '':
				#Maybe the first line is empty?
				continue
			else:
				#Format not recognized.
				#print("Sorry, the converter couldn't recognize the given input file.")
				break

	return True

def get_arguments():
	# Argparse code
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--fold", help="Specifies the line fold, i.e. after how many bases should a new line be inserted.", required=False)
	parser.add_argument("-i", "--input-file", help="Path to a sequence file.", required=True)

	args = vars(parser.parse_args())

	file_path = args['input_file']
	fold = args['fold']

	if fold is not None:
		try:
			fold = int(args['fold'])
		except ValueError:
			#Wrong input given, switching to default.
			fold = 70


	return file_path, fold


if __name__ == "__main__":
	file_path, fold = get_arguments()
	if fold is None:
		fold = 70
	get_format_and_pass_to_converters(file_path, fold)

