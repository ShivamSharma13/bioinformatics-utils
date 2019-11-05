data_path = "../../../data/mirna_org/"


def	parse_pred_file(temp_file_path):
	with open(data_path + temp_file_path, 'r') as f:
		r = f.read()

	rows = r.split('\n')
	columns = []
	for idx, row in enumerate(rows):
			columns.append(row.split('\t'))
	del rows
	
	dump = []

	for _ in columns[1:]:
		try:
			for idx, col in enumerate(_):
				if idx == 0:
					mirna = col
				if idx == 1:
					gene = col
				if col.startswith('[hg19:'):
					try:
						chromosome, chr_start, chr_stop, strand, if_special = get_chr_start_stop(col)
					except IndexError:
						#print(_)
						continue
				if idx == 3:
					a_score = col
				if idx == 4:
					svr_score = col
			#print(chr_start, chr_stop)
			if if_special:
				dump.append(['chr'+chromosome, chr_start.split('-')[0], chr_start.split('-')[1], 
							mirna + ';' + gene, a_score + ';' + svr_score, strand])
				dump.append(['chr'+chromosome, chr_stop.split('-')[0], chr_stop.split('-')[1],
							mirna + ';' + gene, a_score + ';' + svr_score, strand])
			else:
				dump.append(['chr'+chromosome, chr_start, chr_stop, mirna + ';' + gene, a_score + ';' + svr_score, strand])
			del chromosome, chr_start, chr_stop, strand, mirna, gene, a_score, svr_score
		except UnboundLocalError:
			continue
	
	with open(data_path + "mirna_org_bed.bed", 'w') as f:
		for row in dump:
			f.write('\t'.join(row) + '\n')	

	print('Done')
	return 0


def get_chr_start_stop(hg_like):
	#print(hg_like)
	'''
	return chr7 234234 234234 +/-
	'''
	hg_like = hg_like.strip('[]').replace('hg19:', '')
	cols = hg_like.split(':')
	
	#Checking if multiple coordinates exist for a single entry.
	if_special = False
	if ',' in cols[1]:
		if_special = True
		special_cols = cols[1].split(',')
		co_ods_1 = special_cols[0]
		co_ods_2 = special_cols[1]
		return cols[0], co_ods_1, co_ods_2, cols[2], if_special
	return cols[0], cols[1].split('-')[0], cols[1].split('-')[1], cols[2], if_special



if __name__ == "__main__":
	parse_pred_file('temp.txt')