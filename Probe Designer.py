""" 
Installing biopython: uncomment (#)
"""
#pip install biopython --upgrade
#pip install cython
#pip install primer3
#pip install melt
#pip install pandas
#pip install xlsxwriter
#pip install xlrd
#pip install openpyxl

""" 
Load Packages
"""
import pandas as pd
from Bio import Entrez, SeqIO  #biopython
from Bio.Seq import Seq
import melting
from Bio.Blast import NCBIWWW, NCBIXML
from collections import Counter 
from difflib import SequenceMatcher 
import re
import warnings
warnings.filterwarnings("ignore")


"""
Filename: Probe Designer.py

This code will design H-Probes for PLISH
"""

def main():
    global Final_Table, vb, vbs, species, gene
    
    Entrez.email = "stockmca@stanford.edu"  # Always tell NCBI who you are
    print('   ')
    print('If you are searching for multiple genes, please enter one VB per gene in comma separated lists (not case-sensitive).')
    print('Sftpc, Acta2, Emcn and 8, 7, 6 would make probes for Sftpc in VB8, Acta2 in VB7, Emcn in VB6')
    print('There is a limit of 5 genes at a time')
    
    genes = input("what gene(s) do you want to search? ")
    
    if ',' in genes:
        genes=genes.replace(" ", "")
        genes=genes.split(',') 
    elif ' ' in genes:
        genes=genes.split(' ')
          
    species = input("What species do you want to search? (Mus musculus = M, Homo sapiens = H): ").upper()
    vbs=input("What variable Bridge(s) do you want to use? (1, 2, 3, etc.): ")
     
    if ',' in vbs:
        vbs=vbs.replace(" ", "")
        vbs=vbs.split(',')
    elif ' ' in vbs:
        vbs=vbs.split(' ')
        
    # Run Each Gene for alignment
    
    if type(genes)==list:
        for i in range(len(genes)):
            if species=='M' or species=='Mus musculus':
                gene=genes[i].capitalize()
            elif species=='H' or species=='Homo Sapiens':
                gene=genes[i].upper()
            vb=int(vbs[i])-1# need -1 since index begins at 0
        
            # Individual gene search
            accession_number=gene_search(gene, species)
            fasta_sequence=fasta_download(accession_number)
            Final_Table=probe_search(fasta_sequence)
            Final_Table=H_Probe_Table
            Final_Table['VB']=['VB'+str(vb+1) for i in range(len(Final_Table))]
            for i in range(len(Final_Table)):
                if Final_Table.Target_Genes[i]=='a':
                    Final_Table=Final_Table.drop(i)
            Final_Table.to_excel(gene+'PLISH.xlsx', sheet_name=species+gene)
    
            print("______________________________________________________________________________________________________________________________________")
            print("   ")
            print("You have successfully designed H Probes for " +gene+ " in " +species+ " using variable bridge " +str(vb+1))
            print("   ")
            print("A table of H Probes has been saved to the working directory  --  PLISH_Probes.xlsx.")
            print('   ')     
    elif type(genes) == str:
        if species=='M' or species=='Mus musculus':
            gene=genes.capitalize()
        elif species=='H' or species=='Homo Sapiens':
            gene=genes.upper()
        vb=int(vbs)-1# need -1 since index begins at 0
    
        # Individual gene search
        accession_number=gene_search(gene, species)
        fasta_sequence=fasta_download(accession_number)
        Final_Table=probe_search(fasta_sequence)
        Final_Table=H_Probe_Table
        Final_Table['VB']=['VB'+str(vb+1) for i in range(len(Final_Table))]
        for i in range(len(Final_Table)):
            if Final_Table.Target_Genes[i]=='a':
                Final_Table=Final_Table.drop(i)
        Final_Table.to_excel(gene+'PLISH.xlsx', sheet_name=species+gene)

        print("______________________________________________________________________________________________________________________________________")
        print("   ")
        print("You have successfully designed H Probes for " +gene+ " in " +species+ " using variable bridge " +str(vb+1))
        print("   ")
        print("A table of H Probes has been saved to the working directory  --  PLISH_Probes.xlsx.")
        print('   ')     
        

"""
Instructions for searching for gene accession numbers

This will output a list of strings titles acession_numbers for eithera human 
or mouse gene

The user will first input the gene name and then tell the program
whether they want to search for human or mouse accession numbers
"""

def gene_search(gene, species):    
    if species == ('M' or 'Mus Musculus'):
        # Searching ENTREZ, Indicating that you want to search mus musculus 
        # gene acession numbers
        handle = Entrez.esearch(db="nucleotide", term=str('Mus musculus[Orgn] AND ' + gene +'[Gene]'), idtype="acc")
        record = Entrez.read(handle)
        accession_numbers = record["IdList"]
    elif species == ('H' or 'Homo sapiens'):
        # Searching ENTREZ, Indicating that you want to search homo sapiens 
        # gene acession numbers
        handle = Entrez.esearch(db="nucleotide", term=str('Homo sapiens[Orgn] AND ' + gene +'[Gene]'), idtype="acc")
        record = Entrez.read(handle)
        accession_numbers = record["IdList"]
    # Can add more elif commands for other species to search, copy above code, 
    # and change species name to be able to complete the search
    accession_numbers=[seq_ID for seq_ID in accession_numbers if 'NM' in seq_ID]
    if len(accession_numbers)==0:
        accession_numbers=input("No NM accession numbers found, please enter accession number:")    
    elif len(accession_numbers)>=2:
        accession_numbers=accession_numbers[0]
    
    print('   ')
    print("The accession number used for " +gene+ " in "+species+ " is: " + str(accession_numbers))
    
    return accession_numbers
    
"""
Instructions for downloading the FASTA sequence

This will output a string of letters (ATCG) where primers can be drawn

The user will not need to input any information and will run automatically from
gene selected in gene search

"""
def fasta_download(accession_number):
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return str(record.seq) 

""" 
Instructions to search for H-probes

This will output a list of H-probes

The user will not have to enter any information
"""
def probe_search(fasta_sequence):
    global H_Probes
    global H_Probe_Table

    H_Probe_Table=h_probe_design(fasta_sequence)
    
    print('   ')
    print(str(len(H_Probe_Table)) + ' H_Probes were identified to be scored and aligned.')
    print('   ')
    
    H_Probe_Table=append_variable_bridge(vb)
    H_Probe_Table['Score']=score_h_probes()    
    probe_number=len(H_Probe_Table)
    i=1
    for i in range(probe_number):
        if H_Probe_Table.Score[i]>5:
            H_Probe_Table=H_Probe_Table.drop(i)
    H_Probe_Table = H_Probe_Table.reset_index(drop=True)
    print("______________________________________________________________________________________________________________________________________")
    print('   ')
    print('Alignment Running')
    print('   ')
    print("______________________________________________________________________________________________________________________________________")
    H_Probe_Table=blast_h_probes(species)
    return H_Probe_Table
    
""" 
Functions used in probe_search:

h_probe_design(fasta_sequence) # input = fasta sequence
    Makes a list of all L/R H-Probes for the length of the gene
    
blast_h_probes() no inputs necessary
    blasts all 40bp sequences that will bind to mRNAs within the cell - this will 
    add a column of gene targets to the H_Probe_Table
    
append_variable_bridge # input=vb (variable bridge)
    adds VB+CC to each H_Probe - generates sequences to order
    
score_h_probes() no inputs necessary
    will create a score based on TM of each H_Probe, GC content, and ability to 
    form a hairpin

"""
def h_probe_design(fasta_sequence):
    global H_Probes
    fw_probes=[]
    rv_probes=[]
    count=0
    probe_number=len(fasta_sequence)-41
    for i in range(probe_number):
        fw_probe=str(fasta_sequence[i:(i+20)])
        rv_probe=str(fasta_sequence[(i+20):(i+40)])
        fw_probe_1=Seq(fw_probe)
        rv_probe_1=Seq(rv_probe)
        if (fw_probe_1.endswith('T') and rv_probe_1.startswith('A')) or (fw_probe_1.endswith('A') and rv_probe_1.startswith('G')):
            fw_probe_2=fw_probe_1.reverse_complement()
            rv_probe_2=rv_probe_1.reverse_complement()
            if count==0:
                fw_probes=['delete', str(fw_probe_2)]
                rv_probes=['delete', str(rv_probe_2)]
                count+=1
                fw_probes2=fw_probes.pop(0)
                rv_probes2=rv_probes.pop(0)
            else:
                fw_probes.append(str(fw_probe_2)) 
                rv_probes.append(str(rv_probe_2))  
    probes=pd.DataFrame({'Left_Probe':fw_probes,'Right_Probe':rv_probes})
    probes['Location']=range(len(probes))
    ends=probes["Right_Probe"].str.endswith(("CT","AT"))
    H_Probes=probes[ends]
    H_Probes = H_Probes.reset_index(drop=True)
   
    return H_Probes

def blast_h_probes(species):  
    H_Probe_Table["Blast_Sequence"]=H_Probe_Table.Right_Probe+H_Probe_Table.Left_Probe
    result_handle=[]
    probe_number=len(H_Probe_Table)
    target_genes=['a' for i in range(len(H_Probe_Table))]
    count=1
    
    if species=='M':
        species='Mus musculus'
        taxid='txid10090[orgn]'
    elif species=='H':
        species='Homo sapiens'
        taxid='txid9606[orgn]'
        
    for i in range(probe_number):    
        search=H_Probe_Table.Blast_Sequence[i].lower()
        result_handle = NCBIWWW.qblast("blastn",'nt', search, short_query=True, entrez_query=taxid)
        blast_record = NCBIXML.read(result_handle)
        E_VALUE_THRESH = 0.05
        count=1
        print('   ')
        print(gene+ ': ' +search)
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                count+=1
                if (hsp.expect < E_VALUE_THRESH) and (species in alignment.title) and (('ref|NM_' in alignment.title) or ('ref|XM_' in alignment.title)) and str(hsp.strand) == "('Plus', 'Minus')": 
                    a=str(alignment.title)
                    ot_gene= a[a.find('(')+1:a.find('),')]
                    if '(' in ot_gene:
                        loc=ot_gene.find('(')+1
                        ot_gene=ot_gene[loc:]
                    H_Probe_Table.Location[i]=hsp.sbjct_start
                    if ((ot_gene == gene) or (gene in alignment.title)) and (target_genes[i]=='a'):
                        target_genes[i]=ot_gene
                    if ot_gene != gene and  hsp.positives > 20 and (ot_gene not in target_genes[i]):
                        target_genes[i]+=', '
                        target_genes[i]+= str(ot_gene)   
                        print('****Alignment****')
                        print('sequence:', alignment.title)
                        print('length:', alignment.length)
                        print('e value:', hsp.expect)
                        print('Start NT: ' +str(hsp.sbjct_start)+ ' | Gaps: ' +str(hsp.gaps)+ ' | Identical Residues: '+str(hsp.positives)+ '/40 | Strands: ' + str(hsp.strand))
                        print(hsp.query[0:75])
                        print(hsp.match[0:75])
                        print(hsp.sbjct[0:75])
                                  
    H_Probe_Table["Target_Genes"]=target_genes
    for i in range(len(H_Probe_Table)):
        if H_Probe_Table.Target_Genes[i] == 0:
            H_Probe_Table.drop(i)
    return H_Probe_Table
  
def append_variable_bridge(vb):
    cc_l='AGGTCAGGAAT'
    cc_r='ATAGCCAGGTT'
    vbs=variable_bridges()
    cc_lvb=cc_l+vbs.VB_L[vb]
    rvb_cc=vbs.VB_R[vb]+cc_r
    H_Probe_Table['Left_Probe_VB'] =[cc_lvb + s for s in H_Probe_Table.Left_Probe]
    H_Probe_Table['Right_Probe_VB'] =[s + rvb_cc for s in H_Probe_Table.Left_Probe]
    return H_Probe_Table

def score_h_probes():
    TM_Table=tm_probe_calculator()
    GC_Table=gc_probe_calculator()
    tm_score=score_tm(TM_Table)
    gc_score=score_gc(GC_Table)
    hairpin_score=hairpin_calculator()
    dimer_score=dimer_calculator()     
    score=[a + b + c + d for a, b, c, d in zip(tm_score,gc_score,hairpin_score,dimer_score)]
    return score    
     
""" 
Calculators for TM, GC, Hairpin and Dimer Formation

Commands to calculate a 'score' for each H_Probe design with the attached 
variable bridge
"""

def tm_probe_calculator():
    left_probes=H_Probe_Table.Left_Probe.tolist()
    right_probes=H_Probe_Table.Right_Probe.tolist()
    left_probes_mT=[]
    right_probes_mT=[]    
    probe_number=len(H_Probe_Table)
    i=1
    for i in range(probe_number):    
        left_probe_mT=melting.temp(left_probes[i])
        right_probe_mT=melting.temp(right_probes[i])
        left_probes_mT.append(left_probe_mT)
        right_probes_mT.append(right_probe_mT)
    return pd.DataFrame({'Left_Probe_mT':left_probes_mT,'Right_Probe_mT':right_probes_mT})
def score_tm(TM_Table):
    left_probes_mT=TM_Table.Left_Probe_mT.tolist()
    right_probes_mT=TM_Table.Right_Probe_mT.tolist()
    # Score probe pairs on 45<TM<65:
    scoretm=[]
    for x in range(len(left_probes_mT)):
        if 45<x<65:
           scoretm.append(0)
        else:
            scoretm.append(1)
    for x in range(len(right_probes_mT)):
        if 45>x>65:
            scoretm[x]+=1
    return scoretm


def gc_probe_calculator():
    left_probes=H_Probe_Table.Left_Probe_VB.tolist()
    right_probes=H_Probe_Table.Right_Probe_VB.tolist()
    left_probes_gc=[]
    right_probes_gc=[]    
    probe_number=len(H_Probe_Table)
    i=1
    for i in range(probe_number):
        l_gc_counts=Counter(left_probes[i])
        left_probes_gc.append((l_gc_counts["G"]+l_gc_counts["C"])/20)
        r_gc_counts=Counter(right_probes[i])
        right_probes_gc.append((r_gc_counts["G"]+r_gc_counts["C"])/20)
    return pd.DataFrame({'Left_Probe_GC':left_probes_gc,'Right_Probe_GC':right_probes_gc})
def score_gc(GC_Table):
    left_probes_gc=GC_Table.Left_Probe_GC.tolist()
    right_probes_gc=GC_Table.Right_Probe_GC.tolist()
    # Score probe pairs on .2<GC<.8:
    scoregc=[]
    for x in range(len(left_probes_gc)):
        if .2<left_probes_gc[x]<.8:
            scoregc.append(0)
        else:
            scoregc.append(1)
    for x in range(len(right_probes_gc)):
        if .2>right_probes_gc[x]>.8:
            scoregc[x]=+1
    return scoregc  

def hairpin_calculator():
    left_probes=H_Probe_Table.Left_Probe_VB.tolist()
    right_probes=H_Probe_Table.Right_Probe_VB.tolist()
    left_probes_hairpin=[]
    right_probes_hairpin=[]    
    probe_number=len(H_Probe_Table)
    i=1 
    for i in range(probe_number):     
        str1=left_probes[i]
        str1=Seq(str1)
        str2=str1.reverse_complement()
        seqMatch = SequenceMatcher(None,str1,str2) 
        match1 = seqMatch.find_longest_match(0, len(str1), 0, len(str2))
        match1=str(str1[match1.a:match1.a+match1.size])
        matchrv=str(Seq(match1).reverse_complement())
        str1=str(str1)
        for match in re.finditer(match1,str1):
            left=match.end()
        for match in re.finditer(matchrv,str1):
            right=match.start()
        if len(match1)>=4 and (right-left>4):
            left_probes_hairpin.append(str(match1))
        else:
            left_probes_hairpin.append('No')
        str1=right_probes[i]
        str1=Seq(str1)
        str2=str1.reverse_complement()
        seqMatch = SequenceMatcher(None,str1,str2) 
        match1 = seqMatch.find_longest_match(0, len(str1), 0, len(str2))
        match1=str(str1[match1.a:match1.a+match1.size])
        matchrv=str(Seq(match1).reverse_complement())
        str1=str(str1)
        for match in re.finditer(match1,str1):
            left=match.end()
        for match in re.finditer(matchrv,str1):
            right=match.start()
        if int(len(match1))>=4 and (int(right-left)>4):
            right_probes_hairpin.append(match1)
        else:
            right_probes_hairpin.append('No')
            
    scorehp=[]      
    for x in range(len(left_probes_hairpin)):
        if left_probes_hairpin[x] == 'No':
            scorehp.append(0)
        else:
            scorehp.append(5)
    for x in range(len(right_probes_hairpin)):
        if  right_probes_hairpin[x] != 'No':
            scorehp[x]+=5
        return scorehp
def dimer_calculator():
    left_probes=H_Probe_Table.Left_Probe_VB.tolist()
    right_probes=H_Probe_Table.Right_Probe_VB.tolist()
    left_probes_dimer=[]
    right_probes_dimer=[]    
    probe_number=len(H_Probe_Table)
    i=1 
    for i in range(probe_number):     
        str1=left_probes[i]
        str1=Seq(str1)
        str2=str1.reverse_complement()
        seqMatch = SequenceMatcher(None,str1,str2) 
        match1 = seqMatch.find_longest_match(0, len(str1), 0, len(str2))
        match1=str(str1[match1.a:match1.a+match1.size])
        matchrv=str(Seq(match1).reverse_complement())
        str1=str(str1)
        for match in re.finditer(match1,str1):
            left=match.end()
        for match in re.finditer(matchrv,str1):
            right=match.start()
        if len(match1)>=2 and (right-left>0):
            left_probes_dimer.append(match1)
        else:
            left_probes_dimer.append('No')
        str1=right_probes[i]
        str1=Seq(str1)
        str2=str1.reverse_complement()
        seqMatch = SequenceMatcher(None,str1,str2) 
        match1 = seqMatch.find_longest_match(0, len(str1), 0, len(str2))
        match1=str(str1[match1.a:match1.a+match1.size])
        matchrv=str(Seq(match1).reverse_complement())
        str1=str(str1)
        for match in re.finditer(match1,str1):
            left=match.end()
        for match in re.finditer(matchrv,str1):
            right=match.start()
        if len(match1)>=2 and (right-left>0):
            left_probes_dimer.append(match1)
        else:
            right_probes_dimer.append('No')
        
    scorehp=[]      
    for x in range(len(left_probes_dimer)):
        if left_probes_dimer[x] == 'No':
            scorehp.append(0)
        else:
            scorehp.append(5)
    for x in range(len(right_probes_dimer)):
        if  right_probes_dimer[x] != 'No':
            scorehp[x]+=5
        return scorehp

def energy_calculator(l_two):
    dg=0
    for y in range(len(l_two)):
        if l_two[y]=='CT' or l_two[y]=='GA':
            dg+=1.6
        elif l_two[y]=='GC' or l_two[y]=='GG'or l_two[y]=='CC':                    
            dg+=3.1
        elif l_two[y]=='CG':
            dg+=3.6
        elif l_two[y]=='CA' or l_two[y]=='TC' or l_two[y]=='TG' or l_two[y]=='AC' or l_two[y]=='AG':
            dg+=2.0
        elif l_two[y]=='GT':
            dg+=1.3
        elif l_two[y]=='TA':
            dg+=1.0
        elif l_two[y]=='AT':
            dg+=1.5
        elif l_two[y]=='AA' or l_two[y]=='TT':
            dg+=1.9
        elif l_two[y]=='A' or l_two[y]=='T':
            dg+=1.5
        elif l_two[y]=='C' or l_two[y]=='G':
            dg+=3
        return dg
        
'''
Other Commands
'''
# List of VBs
def variable_bridges():
    vb1_l='ACTTACGTCGTTATGG'
    vb1_r='TTATAGGTCGAGTAGT'
    vb2_l='ACTTAGCTATTGATGG'
    vb2_r='TTCTGTGTAGACGACT'
    vb3_l='ACCAGGTTGTAATGG'
    vb3_r='TTAGTATGATGACAAT'
    vb4_l='ACGACTACGAGTAGG'
    vb4_r='TTCTTGTGGTGTAAGT'
    vb5_l='ACGTGGAGTGACCCGG'
    vb5_r='TTACCCCTTGTGAGAT'
    vb6_l='AAGACATGCTGACAGG'
    vb6_r='TTCACACAATATAGAT'
    vb7_l='ATTCATATTGGTACGG'
    vb7_r='TCCAGCATGGAAGATT'
    vb8_l='AGGACAACAAGGTCGG'
    vb8_r='TCAGCGCTAATCACAT'
    vb9_l='ACACTGGGCACGGAGG'
    vb9_r='TCATGAAGAAAAGAAT'
    vb10_l='AGTAAACAACCCATGG'
    vb10_r='TGCGAGAGCCGGAATT'
    vb11_l='AATCGGAAGGCAGTGG'
    vb11_r='TTGAGGCTGCTTATCT'
    vb12_l='AGTCTCAGCTCCGAGG'
    vb12_r='TGGATTTACTTCCAGT'
    vb13_l='ATAGGCCGAAAAAAGG'
    vb13_r='TTGTATACTATTAAAT'
    vb14_l='AGTGGCTTGTTATGGG'
    vb14_r='TGAATGGACTAGCACT'
    vb15_l='ACGGTGTTCGTTTTGG'
    vb15_r='TCCGCTGTTCTACACT'
    vb16_l='AGACTACTCAAACTGG'
    vb16_r='TGATTACACCTTCAGT'
    vb17_l='ATGTCAATTGATGGGG'
    vb17_r='TTAAGAGAACTTTCTT'
    vb18_l='AAATTTTAGCCAATGG'
    vb18_r='TAGCATATCTTCTATT'
    vb19_l='ATGACTTTTCCGCGGG'
    vb19_r='TCTTTTTGTAACAATT'
    vb20_l='ACTTCGCGAAACCAGG'
    vb20_r='TATAGTAACGAGGACT'
    vb21_l='AAAGGGCCCCTGCCGG'
    vb21_r='TGCATCCGCACGATTT'
    vb22_l='AAGCGACAAGTTTGGG'
    vb22_r='TAGGAGCGACGACGAT'
    vb23_l='ATGTCAGGCGAGGCGG'
    vb23_r='TTTATACACGAGACTT'
    vb24_l='ACTGTACATACGGGGG'
    vb24_r='TCCTAACCGTGCGCCT'
    vb_l=[vb1_l,vb2_l,vb3_l,vb4_l,vb5_l,vb6_l,vb7_l,vb8_l,vb9_l,vb10_l,vb11_l,
          vb12_l,vb13_l,vb14_l,vb15_l,vb16_l,vb17_l,vb18_l,vb19_l,vb20_l,vb21_l,
          vb22_l,vb23_l,vb24_l]
    vb_r=[vb1_r,vb2_r,vb3_r,vb4_r,vb5_r,vb6_r,vb7_r,vb8_r,vb9_r,vb10_r,vb11_r,
          vb12_r,vb13_r,vb14_r,vb15_r,vb16_r,vb17_r,vb18_r,vb19_r,vb20_r,vb21_r,
          vb22_r,vb23_r,vb24_r]
    return pd.DataFrame({'VB_L':vb_l,'VB_R':vb_r})    
  

# This provided line is a required at the end of a Python file to call main()
if __name__ =='__main__':
    main()