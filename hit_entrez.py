#!/usr/bin/env python3

from Bio import SeqIO, Entrez
import time
Entrez.email = 'ssharma454@gatech.edu'

def get_ncbi_data(database, any_id):
	#Hit NCBI [respectfully] to obtain data.
	handle = Entrez.efetch(db=database, id=any_id, retmode="xml")
	record = Entrez.read(handle, "genbank")
	handle.close()
	
	#No more than three in a second.
	time.sleep(0.4)

	#Parsing through the output that NCBI gave.
	if len(record) == 1:
		parsed_data = record[0]
		return parsed_data
	else:
		print("Something wrong with getting the parsed_data.")
		return None

if __name__ == "__main__":
	example_database = "nucleotide"
	example_ids = ["NM_001126116", "NM_000546"]
	print([print(get_ncbi_data(example_database, any_id)) for any_id in example_ids])
