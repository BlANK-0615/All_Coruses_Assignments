#!/usr/bin/python3
# Script of ICA2 for BPSM
# Author B224911
# Version 10, 20th Nov 2022

# import modules
import os
import subprocess
import shutil
import pandas as pd
import numpy as np
import re

# Function1 control the input of users
# replace the 'space' in the species and protein names 
# with '_', for later proccess.
def input_judge() :
    while True :
        print('\n##################################################################\
        \nStep2: choose the species and protein you want to analysis\n\
        ####################################################\n')
        spec = '_'.join((input("What's the species you want to analysis?\n")).split())
        pro = '_'.join((input("What's the protein you want to analysis?\n")).split())
        if spec == "" or pro == "" :
            print('You enter the invaid species or protien, please enter again.')
            continue
        else :
            return  spec, pro

# Function2 count the number of sequences of research.
def count_function(filename) :
    subprocess.call(f"esearch -db protein -query '{protein}[Protein Name] AND {species}[Organism]' |\
    efetch -format fasta > {filename}.fasta",shell=True)
    count = subprocess.run(f"grep -o '>' {filename}.fasta | wc -l", shell=True, capture_output=True)
    num_seq = int(count.stdout.decode())
    return num_seq

# Function3 count the number of species
def species_count(filename) :
    gene_list = []
    my_file = open(f"{filename}.fasta")
    for eachline in my_file :
        if eachline.startswith(">") :
            split1 = eachline.split("[")[1]
            split2 = split1.split("]")[0]  # extract the spieces by '.split()'
            gene_list.append(split2)
            # print(gene_list)       
    my_file.close()   
    nodup_list = list(set(gene_list))
    num_spec = len(nodup_list)
    return num_spec

# Function4 Extract information of fasta data
# and calculate the sequence length and store it in the list
def extract_fasta(filename) :
    acc_list_in = []
    spec_list_in = []
    my_file_contents = open(f'{filename}.fasta').read()
    acc_list_in = re.findall(r'>([^\s]*)', my_file_contents)  # extract accession number
    spec_list_in = re.findall(r'\[(.*)]',my_file_contents)  # extract spieces 
    
    seq_list_in = []
    seq_len_list_in = []
    line_list = []
    one_seq = []
    with open(f'{filename}.fasta') as my_file :
        for eachline in my_file:
            if eachline.startswith('>') :
                if line_list != [] :
                    one_seq = ''.join(line_list)
                    seq_list_in.append(one_seq)
                    length = len(one_seq)
                    seq_len_list_in.append(length)
                    line_list = []
            else :
                if not eachline.startswith('>') :    
                    line_list.append(eachline.replace('\n','')) 
    one_seq = ''.join(line_list)
    seq_list_in.append(one_seq)  # extract sequence
    length = len(one_seq)
    seq_len_list_in.append(length)  # extract the length of sequence
    return acc_list_in, spec_list_in,seq_list_in, seq_len_list_in

# Function5 runing clustalo
def clustalo(filename) :
    print('\n########################\nStep5:runing clustalo\n########################\n')
    subprocess.call('mkdir clustalo_dic', shell=True)  # make a dictionary for clustalo 
    print('\n#######################\n runing multi-alignment\n#######################\n')
    subprocess.call(f'clustalo -i {filename}.fasta -o ./clustalo_dic/clustalo_aglin.txt', shell=True)  # run multi-align
    print('\n#########################\n generating distance_matrix\n#########################\n')
    subprocess.call(f'clustalo -i {filename}.fasta \
    --distmat-out=./clustalo_dic/distance_matrix.ma --full', shell=True)  # generate matrix
    print('\n######################\n generating plotcon\n######################\n\
    For plotcon,please enter the number of winsize(default is 4)\n')
    subprocess.call(f'plotcon {filename}.fasta -graph pdf', shell=True)  # generate plot of similarity
    subprocess.call('mv *.pdf ./clustalo_dic', shell=True)  # move to clustalo dic
    return

# Function6：runing optional package and generate report
def option(package) :
    print('\n######################\n\
    Step6.runing patmatmotifs,garnier,pepstats\n######################\n')
    subprocess.call("touch medium.fasta", shell=True) 
    subprocess.call(f"mkdir {package}_dic", shell=True)  # make a dictionary for optional package
    for i in range(len(acc_list)) :
        medium_file = open("medium.fasta","w")  # make a medium file to convey the sequence for emboss analysis
        medium_file.write(f">{acc_list[i]}\n{seq_list[i]}")  # rewrite the medium file
        medium_file.close()
        subprocess.call(f"{package} -sequence medium.fasta \
        -outfile ./{package}_dic/{acc_list[i]}.{package}", shell=True)  # run optional packages 
        if package == 'pepstats' :  
            subprocess.call("touch pepstats_report.txt", shell=True)
            with open (f"./pepstats_dic/{acc_list[i]}.pepstats","r")as pepstats_count :
                        report_file = open("pepstats_report.txt","a")
                        line_num = 0
                        for eachline in pepstats_count :
                            if line_num < 8 :  # extract the first eight lines of content
                                if line_num == 0 :
                                    report_file.write(f"\n###{eachline}###")
                                else :
                                    report_file.write(f"{eachline}")  # iteratively write to file
                                line_num += 1
                            else :
                                break
                        report_file.close()

        if package == 'patmatmotifs' :
            subprocess.call("touch patmatmotifs_report.txt", shell=True)        
            with open (f"./patmatmotifs_dic/{acc_list[i]}.patmatmotifs","r")as patmatmotifs_count:  # genetate report
                patmatmotifs_str=''
                report_file=open("patmatmotifs_report.txt","a")
                for eachline in patmatmotifs_count :
                    Hitcount=re.search(r"HitCount: (.*)", eachline)
                    Motif=re.search(r"Motif = (.*)", eachline)
                    if Hitcount :
                        count_str = Hitcount.group(1)  # extract HitCount
                    if Motif :
                        patmatmotifs_str += f',{Motif.group(1)}'  # extract motifs
                report_file.write(f"{acc_list[i]}: HitCount:{count_str}, Motif:{patmatmotifs_str}\n")  # iteratively write to file
                report_file.close() 

    if package == 'garnier' :
            garnier_list=os.listdir('garnier_dic')
            subprocess.call("touch garnier_report.txt", shell=True)
            for i in range(len(garnier_list)) :
                with open (f"./garnier_dic/{garnier_list[i]}","r")as garnier_count :
                            report_file=open("garnier_report.txt","a")
                            for eachline in garnier_count :
                                Hitcount=re.search(r"HitCount: (.*)", eachline)
                                residue=re.search(r"Residue totals: (.*)", eachline)
                                per=re.search(r"percent: (.*)", eachline)
                                if Hitcount :
                                    count_str=Hitcount.group(1)  # extract motifs
                                if residue :
                                    residue_str=f'{residue.group(1)}'  # extract residues
                                if per :
                                    per_str=f'{per.group(1)}'  # extract residues percent
                            report_file.write(f"{garnier_list[i]}: HitCount:{count_str}, Residue totals:{residue_str}, percent: {per_str}\n")  # iteratively write to file
                            report_file.close()
    return


# step1. set work directory
workspace = input("Step1:Please enter the name of your workspace?\n")
os.mkdir(workspace) 
os.chdir(workspace)

# step2:
# get protein names and species name from users
# constrain the number of sequences
# call function
species, protein = input_judge()
filename = species+'__'+protein  # define the parameter of function
count_seq = count_function(filename)
# judge the number of sequence
while count_seq > 1000 :
    print("\n################################################\n\
    The number of sequences from search is \
    larger than 1000.\nThe program could not run.\nPlease enter more specific species or protein names\n\
    ####################################################\n")
    choice1 = input("Do you want to continue?  y/n\n")
    if choice1 == 'y' :
        species, protein = input_judge()
        count_seq = count_function(species, protein)    
    if choice1 == 'n' : 
        print('Thank you for your using.\n')
        exit()  # quit program 
    else :
        print('you dont enter invaild value.\n')
        continue

# step3:
# count the spices and ask the user option
count_spec = species_count(filename)  # count the species
while True :
    print('\n####################################################\n\
    Step3:your sequences contain '+str(count_spec)+' species.\n\
    ####################################################\n')
    choice2 = input("Do you want to continue?  y/n\n")
    if choice2 == 'y' :
        print('Contine processing')
        break
    if choice2 == 'n' :
        print('Thank you for your using.\n')
        exit()
    else :
        print('you dont enter invaild value.\n')
        continue

# step4:
#Extract sequence information and build dictionaries and csv filefor further analysis later
acc_list, spec_list, seq_list, seq_len_list = extract_fasta(filename)
fasta_inform = {'Accession number':acc_list, 'Spieces':spec_list, 'Sequence':seq_list, 'Length':seq_len_list}
df = pd.DataFrame(fasta_inform)
df.to_csv("fasta_inform.csv",sep=",",header=True)
print("\n####################################################\n\
Step4:There is a csv file about fasta information in your workspace,\n\
you could look at that,and decide whether to limit the minimum length of \n\
the sequence for subsequent analysis.\n\
####################################################\n")

# step5:
# Perform sequence filter, multiple sequence comparison.
# determine the level of sequence conservation, and plot.
# (1)judge the choice3 input and run pullseq
while True :
    choice3 = input("Do you want to limit the minimum length of the sequence? y/n\n")
    if choice3 == 'y' :
        subprocess.call(f'cp /localdisk/data/BPSM/ICA2/pullseq ./', shell=True)
        mini = int(input('Please enter your minimum length length of sequence\n'))
        subprocess.call(f'./pullseq -i {filename}.fasta -m {mini} > filter.fasta', shell=True)
        my_fliter_contents = open('filter.fasta').read()
        fliter_list = re.findall(r'>([^\s]*)', my_fliter_contents)  # fliter list for further process
        break
    elif choice3 == 'n' :
        print('Contine processing')
        break
    elif choice3 != 'y' and choice3 != 'n' :
        continue

# (2)clustalo-align, distance_matrix, plotcon
if choice3 == 'y' and len(acc_list) != len(fliter_list) :  # compare if the pullseq change the number of sequences
    clustalo('filter')
else :
    clustalo(filename)

# step6:
# Scan protein sequence(s) with motifs from the PROSITE database
# By patmatmotifs, and generate a brief report of the results
# Run optional choice
# (1)patmatmotifs
if choice3 == 'y' and len(acc_list) != len(fliter_list) :  # compare if the pullseq change the number of sequences
    acc_list, schoipec_list, seq_list, seq_len_list = extract_fasta('filter')  # change the list for exacted sequences
    option('patmatmotifs')
else :
    option('patmatmotifs')
# (2)option1-----garnier
option('garnier')
# (3)option2-----pepstats
option('pepstats')

