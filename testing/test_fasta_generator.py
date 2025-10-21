import random

length = 100000  # 100 kb
seq = ''.join(random.choices('ATCG', k=length))
with open('sequence_1kb.fasta', 'w') as f:
    f.write('>synthetic_100kb_sequence\n')
    for i in range(0, len(seq), 70):
        f.write(seq[i:i+70] + '\n')