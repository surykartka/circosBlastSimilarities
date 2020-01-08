from Bio import SeqIO
import os

window = 200
step_size = 200
genome_fasta = 'U18466.2.fa'
blast_output = 'output'

def blast():
	## assume one chromosome
	for rec in SeqIO.parse(open(genome_fasta), 'fasta'):
		seq = str(rec.seq)

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
		cmd = 'blastn -word_size 4 -outfmt 6 -query query -strand both -db U18466.2 -out o -evalue 0.001'
		os.system(cmd)
		os.system('cat o >> '+blast_output)
		output += open('o').readlines()
	os.system('rm o')

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
	print(table)
	print(len(table))

def make_circos(table):


if __name__ == '__main__':
	get_table()

