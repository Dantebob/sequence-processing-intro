import Bio
from Bio.Seq import Seq

def Unique_Dates(date, dates):
    for d in dates:
        if date == d:
            return False
    dates.append(date)
    return True

def Find_VarSites(seq_list):
    var_sites = 0
    pos_vs = []
    standard = seq_list[0]
    for pos in range(len(standard)):#compare sequences
        for seq in seq_list:
            if standard[pos] != seq[pos]:
                var_sites +=1
                break
    return var_sites

def Protein_fasta(dna_sequences):
    count = 0
    protein_sequences = []
    for s in dna_sequences:
        st = Seq(s)
        st = st.translate()
        protein_sequences.append(st)
    with open('Protein.fasta', 'w') as pf:
        while count < len(protein_sequences):
            pf.write(">Sample"+str(count)+"_1/24/2023 4310 gene\n")
            pf.write(str(protein_sequences[count]+"\n"))
            count +=1
    return protein_sequences
with open('ExampleAlignment.fasta', 'r') as ea:
        fasta = ea.readlines()

sequences = []
dates = []
num_of_samples = 0
for line in fasta:
    if line[0] != ">":
        sequences.append(line.strip())
    else:
        num_of_samples += 1
        str_start = line.find('_')
        str_end = line.find(' ')
        Unique_Dates(line[str_start+1 :str_end], dates)
unique_dates = len(dates)

vs_in_dna = Find_VarSites(sequences)
protein_sequences = Protein_fasta(sequences)
vs_in_protein = Find_VarSites(sequences)

with open('log.txt', 'w') as log:
    log.write('Number of samples: '+ str(num_of_samples) + '\nunique dates: ' + str(unique_dates) +
            '\nvariable sites in DNA: ' + str(vs_in_dna) + '\nvariable sites in protein: ' + str(vs_in_protein))
