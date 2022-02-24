import sys
import Bio
from Bio import SeqIO
import random

def merge_cds_to_gene(gff_file,gene_file):
    cds_dict = {}
    with open(gff_file, 'r') as f:
        for line in f:
            line = line.strip()
            #line don't start with ##
            if not line.startswith('#'):
                #get line that have column 3 is CDS and column 6 is +
                if line.split()[2] == 'CDS' :
                    #get chromosome name
                    chr_name = line.split()[0]
                    #get cds name
                    # cds_name = line.split()[8].split(';')[0].split('=')[1]
                    #get cds start and end
                    cds_start = int(line.split()[3])
                    cds_end = int(line.split()[4])
                    strand = line.split()[6]
                    #get gene name
                    # cds_name = line.split()[8].split(';')[0].split('=')[1]
                    gene_name = line.split()[8].split(';')[2].split('=')[1]
                    #get number of cds in cds name that after "CDS"
                    # cds_num = int(cds_name.split('CDS')[1])
                    #add to dict by use gene_name and chr_name as key
                    # cds_name = cds_name + '|' + chr_name
                    gene_name = gene_name + '|' + chr_name
                    # cds_dict[cds_name] = [cds_start, cds_end, strand]
                    if gene_name in cds_dict:
                        #if cds_start is less than cds_start in dict, set cds_start to cds_start in dict
                        if cds_start < cds_dict[gene_name][0]:
                            cds_dict[gene_name][0] = cds_start
                        #if cds_end is greater than cds_end in dict, set cds_end to cds_end in dict
                        if cds_end > cds_dict[gene_name][1]:
                            cds_dict[gene_name][1] = cds_end
                    else:
                        cds_dict[gene_name] = [cds_start, cds_end, strand]
                
    f.close()
    #write to file
    with open(gene_file, 'w') as f:
        for gene in cds_dict:
            gene_name = gene.split('|')[0]
            chr_name = gene.split('|')[1]
            f.write(chr_name + '\t' + str(cds_dict[gene][0]) + '\t' + str(cds_dict[gene][1]) + '\t' + gene_name + '\t' + cds_dict[gene][2] + '\n')
    f.close()
    
#create bed file from cds file
def gff_to_bed(cds_file,bed_file):
    with open(cds_file, 'r') as f:
        for line in f:
            line = line.strip()
            #get gene name
            start = line.split('\t')[1]
            end = line.split('\t')[2]
            gene_name = line.split('\t')[3]
            #check strand
            gene_length = (int(end) - int(start))+1 #/3 ok
 
            with open(bed_file, 'a') as f:
                f.write(gene_name + '\t' + '0' + '\t' + str(gene_length-1) + '\t' + "cds" + '\t' + str(gene_length-1) + '\t' + '+' + '\n') #use strand + only
    
    f.close()

def extract(fasta_file,coordinate_file,extract_file):
    #read fasta file
    fasta_dict = SeqIO.index(fasta_file, "fasta")
    #read coordinate file
    with open(coordinate_file, 'r') as f:
        for line in f:
            line = line.strip()
            chr_name = line.split()[0]
            #get start and end
            start = int(line.split()[1])
            end = int(line.split()[2])
            name = line.split()[3]
            strand = line.split()[4]
            start = start -1 
            #get sequence
            if strand == '+':
                seq = str(fasta_dict[chr_name].seq[start:end])
            elif strand == '-':
                seq = str(fasta_dict[chr_name].seq[start:end])
                #reverse complement sequence of seq
                seq = seq.translate(str.maketrans("ATGC", "TACG"))
                seq = seq[::-1]
                
            #write to file
            with open(extract_file, 'a') as f:
                f.write(">" + name + '\n' + seq + '\n')
    f.close()
            
#get distance of closest gene from gff with start position of chromosome
def get_distance_from_gff(cds_file, chromosomeCEC, startCEC, strandCEC):
    #new list for genes
    genes = []
    # distance = 10000000
    with open(cds_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.split('\t')[0] == chromosomeCEC and line.split('\t')[4] == strandCEC:
                gene_start = int(line.split('\t')[1])
                gene_end = int(line.split('\t')[2])
                gene_name = line.split('\t')[3]
                startCEC = int(startCEC)
                # print(start,strand,gene_name,gene_start,gene_end)
                if startCEC >= gene_start and startCEC <= gene_end: #check if startCEC is in gene
                    if strandCEC == '+':
                        distance = startCEC - gene_start 
                        if genes == []:
                        #add gene name and distance to list
                            genes = [gene_name, distance]
                        elif distance < genes[1]:
                            genes = [gene_name, distance]
                    elif strandCEC == '-':
                        distance = gene_end - startCEC
                        # value in genes list
                        if genes == []:
                            #add gene name and distance to list
                            genes = [gene_name, distance]
                        elif distance < genes[1]:
                            genes = [gene_name, distance]  
                
    f.close()
    return genes
  
def checkCEC_inCDS(gff_file, cec_file):
    #check position of cec file in cds file
    cdsdict = {}
    with open(gff_file, 'r') as f:
        for line in f:
            line = line.strip()
            #line don't start with ##
            if not line.startswith('#'):
                #get line that have column 3 is CDS
                if line.split()[2] == 'CDS' :
                    #get chromosome name
                    chr_name = line.split()[0]
                    #get cds start and end
                    cds_start = int(line.split()[3])
                    cds_end = int(line.split()[4])
                    strand = line.split()[6]
                    #get gene name
                    cds_name = line.split()[8].split(';')[0].split('=')[1]
                    #add to dict
                    cdsdict[chr_name, cds_start, cds_end, strand] = cds_name
    f.close()
    #open cec file
    with open(cec_file, 'r') as f:
        #skip header
        next(f)
        for line in f:
            line = line.strip()
            #split line with comma   
            name = line.split(',')[0]   
            chromosome = line.split(',')[1]
            startCEC = int(line.split(',')[2])
            strand = line.split(',')[3]
            #check if cec is in cds
            #if startCEC is between cds start and end
            for cds in cdsdict:
                if cds[0] == chromosome and startCEC >= cds[1] and startCEC <= cds[2] and strand == cds[3]:
                    #write to file                 
                    with open(cec_file + '_inCDS.csv', 'a') as f:
                        f.write(name + ',' + chromosome + ',' + str(startCEC) + ',' + strand + ',' + cdsdict[cds] + '\n')
    f.close()

#read file with chromosome and start position
def read_file(file,cds_file,distance_file):
    cec_dict = {}
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            #split line with comma    
            name = line.split(',')[0]  
            chromosome = line.split(',')[1]
            startCEC = int(line.split(',')[2])
            strand = line.split(',')[3]
            #keep in dict to check duplicate
            if name not in cec_dict:
                cec_dict[name] = 1
                genesdistance = get_distance_from_gff(cds_file, chromosome, startCEC, strand)
                # #write gene name and distance to file
                if genesdistance != []:
                    # print(line,genesdistance[0], genesdistance[1])
                    with open(distance_file, 'a') as d:
                        d.write(line + ',' + str(genesdistance[1]) + ',' + genesdistance[0] + '\n')
                
    f.close()
    d.close()

def dis2csv(distancefile,csvfile):
    with open(distancefile,'r') as dis, open(csvfile,'w') as csv:
        csv.write('id,strand,rname,start,stop,copy_number,sequence,cigar,mdz,number_of_mappings,query_length,alignment_length\n')
        #auto increment id
        id = 0
        for line in dis:
            id += 1
            line = line.strip()
            #split line with comma
            line = line.split(',')
            rname = line[6]
            start = int(line[5])
            stop = start +80
            #random copy number from 1 to 10
            copy_number = 1
            #create random sequence with random length
            query_length = random.randint(20,5000)
            alignment_length = query_length
            sequence = ''.join(random.choice(['A','T','G','C']) for i in range(alignment_length))
            #random cigar
            cigar = str(query_length) + "M"
            mdz = ""
            number_of_mappings = 1
            #write to csv file
            csv.write(str(id) + ',' + "1" + ',' + rname + ',' + str(start) + ',' + str(stop) + ',' + str(copy_number) + ',' + sequence + ',' + cigar + ',' + mdz + ',' + str(number_of_mappings) + ',' + str(query_length) + ',' + str(alignment_length) + '\n')
    dis.close()
    csv.close()


#get dict from merge_cds_to_gene #input:gff file and print dict to output file
#merge_cds_to_gene("PlasmoDB-55_Pfalciparum3D7.gff","PlasmoDB-55_Pfalciparum3D7_cds.txt")

#create bed file from cds file    #input file, output file
#gff_to_bed("PlasmoDB-55_Pfalciparum3D7_cds.txt","PlasmoDB-55_Pfalciparum3D7_gene.bed")

#Extract CDS sequence from genome.fasta plasmo_cds.txt => extracted_seq.fasta
#extract("PlasmoDB-55_Pfalciparum3D7_Genome.fasta","PlasmoDB-55_Pfalciparum3D7_cds.txt","extracted_seq_plasmo55.fasta")
#position start,cds,output file
#checkCEC_inCDS("PlasmoDB-55_Pfalciparum3D7.gff", "CEC2_positions.csv")
read_file("CEC2_positions.csv_inCDS.csv","PlasmoDB-55_Pfalciparum3D7_cds.txt","distance_from_position.txt")
#Create csv to insert to SQLite database distance_from_position.txt => transcriptSQL.csv
dis2csv('distance_from_position.txt','transcriptSQL.csv')