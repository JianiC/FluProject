from Bio import Entrez, SeqIO
from StringIO import StringIO

Entrez.email = "richard.neher@tuebingen.mpg.de"     # Always tell NCBI who you are

outgroups = {#'H9':'U26830',
             'H9_77':'CY113261.1',
             #'H1N1pdm':'AF455680',
             #'Vic':'CY018813',
             #'Yam':'CY019707'
             }

for virus, genbank_id in outgroups.iteritems(): 
    handle = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="gb")
    seq= SeqIO.read(StringIO(handle.read()), format = 'genbank')
    SeqIO.write(seq, '/Users/yujiazhou/Documents/nextflu/H9_nextflu-master/augur/source-data/'+virus+'_outgroup.gb', format='genbank')

