from Bio import SeqIO
import os, pyBigWig

conf_file = 'data/circos.conf'

window = 200
step_size = 200
genome_fasta = 'blast/U18466.2.fa'
#genome_fasta = 'blast/FR682468.1.fa'
blast_output = 'blast/output1'
#blast_output = 'blast/output2'

circos_path = '/Users/dor/apps/circos-0.69-9'

add_expression = True
expression_files = ['data/RNAseq_S3_plus.bw', 'data/RNAseq_S3_minus.bw', 
                    'data/RNAseq_S5_plus.bw', 'data/RNAseq_S5_minus.bw'] ## bigwigs

add_gff = True
gff_file = 'data/genes.gff3'

chr_name = 'ASFV_BA71V'

def blast():
	## assume one chromosome
	for rec in SeqIO.parse(open(genome_fasta), 'fasta'):
		seq = str(rec.seq)

	genome_db = genome_fasta.split('.fa')[0]
	cmd = 'makeblastdb -dbtype nucl -in %s -out %s' % (genome_fasta, genome_db)
	os.system(cmd)

	os.system('rm '+blast_output)

	output = []
	n = len(seq)
	for i in range(0, n, step_size):
		if i + window < n:
			end = i + window
		else:
			end = n
		query = seq[i : end]
		with open('query', 'w') as f:
			print('>seq_%d-%d'%(i,end), query, sep='\n', file=f)
		cmd = 'blastn -word_size 4 -outfmt 6 -query query -strand both -db %s -out o -evalue 0.001' % genome_db
		os.system(cmd)
		os.system('cat o >> '+blast_output)
		output += open('o').readlines()
	os.system('rm o')
	os.system('rm query')

def overlap( s1, e1, s2, e2 ):
	return max(0, min(e1, e2) - max(s1, s2) + 1)

def get_table(min_hit_length=100, max_evalue=0.00001):
	table = []
	for line in open(blast_output):
		tab = line.split()
		start = int(tab[0].split('seq_')[1].split('-')[0]) + 1
		add_start, add_end = int(tab[6]) - 1, int(tab[7]) - 1
		start = start + add_start
		end = start + add_end
		evalue = float(tab[-2])
		if evalue > max_evalue:
			continue
		hit_length = int(tab[3])
		if hit_length < min_hit_length:
			continue
		hit_start, hit_end = int(tab[-4]), int(tab[-3])
		orientation = '+'
		if hit_start > hit_end:
			orientation = '-'
			hit_start, hit_end = hit_end, hit_start
		if overlap(start, end, hit_start, hit_end):
			continue
		if hit_start < start or hit_end < end:
			continue
		table.append(((start, end), (hit_start, hit_end), orientation, evalue))
	return table

def write_data(bigwig, data_file, window=100):
	bw = pyBigWig.open(bigwig)
	chrom_d = bw.chroms()
	chrom = list(chrom_d.keys())[0]
	chrom_size = chrom_d[chrom]
	with open(data_file, 'w') as f:
		for i in range(0, chrom_size, window):
			j = min(i+window, chrom_size)
			print('ASFV', i, j, bw.stats(chrom, i, j, 'mean')[0], sep=' ', file=f)

def make_circos(table):
	for rec in SeqIO.parse(open(genome_fasta), 'fasta'):
		genome_length = len(rec.seq)

	## make karyotype file
	karyo_file = os.path.join(circos_path, 'data', 'karyotype', 'ASFV.karyotype.txt')
	with open(karyo_file, 'w') as f:
		print('chr - ASFV %s 0 %d vlgrey' % (chr_name, genome_length), file=f)



	## make links file
	links_file = 'data/links.txt'
	evalue_levels = [1e-50, 1e-10]
	colors = ['greys-5-seq-%d' % x for x in range(2, 5)]
	colors.reverse()
	with open(links_file, 'w') as f:
		for ((start, end), (hit_start, hit_end), orientation, evalue) in table:
			for i,evalue_level in enumerate(evalue_levels):
				if evalue <= evalue_level:
					color = colors[i]
					break
			else:
				color = colors[-1]
			if orientation == '-':
				hit_start, hit_end = hit_end, hit_start
			print('ASFV', start, end, 'ASFV', hit_start, hit_end, 'color='+color, sep=' ', file=f)

	## make expression tracks
	if add_expression:
		for expression_file in expression_files:
			data_file = '%s.txt' % expression_file.split('.')[0]
			write_data(expression_file, data_file)

	## highlight genes
	if add_gff:
		genes_file = 'data/genes.txt'
		new_genes_file = 'data/new_genes.txt'
		with open(new_genes_file, 'w') as f2:
			with open(genes_file, 'w') as f1:
				for line in open(gff_file):
					tab = line.split()
					if 'id-' in tab[3]:
						continue
					print('ASFV', int(tab[1])-1, int(tab[2])-1, sep=' ', file=f1)
					if 'pNG' in tab[3]:
						print('ASFV', int(tab[1])-1, int(tab[2])-1, sep=' ', file=f2)

	cmd = '%s/bin/circos -conf %s' % (circos_path, conf_file)
	os.system(cmd)



if __name__ == '__main__':
	#blast()
	table = get_table()
	make_circos(table)






