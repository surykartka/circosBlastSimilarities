from Bio import SeqIO
import os, csv
from collections import OrderedDict

window = 100
step_size = 100

blast_output_dir = 'blast_output'
circos_output_dir = 'circos'

def blast(genome_fasta, blast_db, blast_output):
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
		cmd = 'blastn -word_size 4 -outfmt 6 -query query -strand both -db %s -out o -evalue 0.001' % blast_db
		os.system(cmd)
		os.system('cat o >> '+blast_output)
		output += open('o').readlines()
	os.system('rm o')

def overlap( s1, e1, s2, e2 ):
	return max(0, min(e1, e2) - max(s1, s2) + 1)

e2color = OrderedDict({1e-40:'blues-4-seq-4', 1e-30:'blues-4-seq-3', 
						1e-20:'blues-4-seq-3', 1e-10:'blues-4-seq-2', 
						1e-2:'blues-4-seq-2'})
def get_color(evalue):
	for e in e2color:
		if evalue <= e:
			return e2color[e]
	return e2color[e]

def get_table(blast_output, min_hit_length=50, max_evalue=0.00001, name1='Georgia', name2='BA71V'):
	
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
		table.append(((start, end), (hit_start, hit_end), orientation, evalue))
	
	print(len(table), 'links')

	with open(circos_output_dir + '/links.txt', 'w') as f:
		for ((start, end), (hit_start, hit_end), orientation, evalue) in table:
			c = get_color(evalue)
			print(name1, start, end, name2, hit_start, hit_end, 'color='+c, sep=' ', file=f)

	## highlight genes

	BA71V = {'early':set(), 'late':set()}

	for row in csv.DictReader(open('data/BA71V_expression_table.txt'), delimiter='\t'):
		t = 'early' if 'early' in row['Gene Type'].lower() else 'late'
		BA71V[t].add(row['Gene Name'])

	i = 0
	with open(circos_output_dir+'/MGF.txt', 'w') as f1:
		with open(circos_output_dir+'/genes.txt', 'w') as f2:

			with open(circos_output_dir+'/late_genes.txt', 'w') as f3:
				with open(circos_output_dir+'/early_genes.txt', 'w') as f4:
					with open(circos_output_dir+'/NC_genes.txt', 'w') as f5:

						for line in open('data/U18466.2.gff3'):
							tab = line.split()
							if 'id-U' in line:
								continue
							if 'member of multigene family ' in line:
								f = f1
							else:
								f = f2
							print('BA71V', int(tab[1])-1, int(tab[2])-1, sep=' ', file=f)
							if tab[3] in BA71V['early']:
								f = f4
							elif tab[3] in BA71V['late']:
								f = f3
							else:
								f = f5
							print('BA71V', int(tab[1])-1, int(tab[2])-1, sep=' ', file=f)

						for row in csv.DictReader(open('data/Georgia_expression_table.txt'), delimiter='\t'):
							if 'MGF' in row['Gene']:
								f = f1
							else:
								f = f2
							print('Georgia', row['ORF start'], row['ORF stop'], sep=' ', file=f)
							if row['Expression Stage'] == 'Early':
								f = f4
							elif row['Expression Stage'] == 'Late':
								f = f3
							else:
								f = f5
							print('Georgia', row['ORF start'], row['ORF stop'], sep=' ', file=f)


if __name__ == '__main__':
	#get_table()
	#blast('data/U18466.2.fa', 'data/Georgia', blast_output_dir + '/BA71V_Georgia')
	blast_output = blast_output_dir + '/Georgia_BA71V'
	#blast('data/FR682468.1.fa', 'data/BA71V', blast_output)
	get_table(blast_output)
