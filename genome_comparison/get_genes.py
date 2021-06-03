from Bio import SeqIO, Seq

def overlap( s1, e1, s2, e2 ):
	return max(0, min(e1, e2) - max(s1, s2) + 1)

def get_BA71V():
	for rec in SeqIO.parse(open('data/U18466.2.fa'), 'fasta'):
		genome = str(rec.seq)

	with open('data/BA71V_prot.fa', 'w') as f1:
		with open('data/BA71V_cds.fa', 'w') as f2:
			for line in open('data/U18466.2.gff3'):
				tab = line.split('\t')
				if not 'id' in tab[3]:
					start, end, strand, descr = int(tab[1]), int(tab[2]), tab[5], tab[-1].strip()
					name = tab[3]
					seq = genome[start-1:end]
					if strand == '-':
						seq = Seq.reverse_complement(seq)
					if len(seq)%3!=0 or (seq[-1]!='*' and seq[:3]!='ATG'):
						print(name, seq)
					else:
						prot = Seq.translate(seq).strip('*')
						print('>'+name, prot, sep='\n', file=f1)
					print('>'+name, seq, sep='\n', file=f2)


if __name__ == '__main__':
	get_BA71V()