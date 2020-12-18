library_valid = True
try:
    from ftplib import FTP
    from pathlib import Path
    import os, sys, os.path
    import urllib.request
    from ete3 import NCBITaxa
    import Bio
    import time
    import sys
    import random
except:
    print("Missing required model(s), make sure you have installed:")
    print("ftplib")
    print("pathlib")
    print("urllib")
    print("ete3")
    print("Bio")

if library_valid:
    def findLength(taxid):
        url = 'https://api.ncbi.nlm.nih.gov/genome/v0/expected_genome_size?species_taxid=' + taxid
        urllib.request.urlretrieve(url, "API_temp")

        file = open("API_temp",'r')
        lines = file.readlines()
        file.close

        for line in lines:
            if "expected_ungapped_length" in line:
                trim1 = line.split('>')
                trim2 = trim1[1].split('<')
        return (int(trim2[0]))
            
    def checkNdownload(ftp, filename, os, asmbly_search, err):
        try: #in case the "file" is a directory, open it and create a folder in local and open it
            ftp.cwd(filename)
            if not os.path.exists(filename):
                os.mkdir(filename)
            os.chdir(filename)
            filename_list = ftp.nlst()
            for element in filename_list:
                ftp, err = checkNdownload(ftp, element, os, asmbly_search, err)
            os.chdir("../")
            ftp.cwd("../")
        except: #in case the "file" is actually a file
            #print("Working FTP directory:", ftp.pwd() + '/' + filename)
            if asmbly_search == 'a':
                print("downloading", filename, "in", os.getcwd())
                local_filename = os.path.join(filename)
                file = open(local_filename, 'wb')
                ftp.retrbinary('RETR '+ filename, file.write)
                file.close()
            elif asmbly_search == 's':
                if "genomic.fna.gz" in filename or "gbff.gz" in filename:
                    print("downloading", filename, "in", os.getcwd())
                    local_filename = os.path.join(filename)
                    try:
                        file = open(local_filename, 'wb')
                        ftp.retrbinary('RETR '+ filename, file.write)
                        file.close()
                    except:
                        print(filename + "cannot be downloaded")
                        err = True
                else:
                    pass
            else:
                pass


        return(ftp, err)

    need_input = True
    while need_input:
        user_input = input("Please enter a valid working directory: ")
        if os.path.isdir(user_input):
            path = user_input
            need_input = False
        else:
            print("Invalid directory! Example directory: /mnt/c/Users/username/genome_data")
    os.chdir(path)

    print("Loading NCBI database")
    ncbi = NCBITaxa()
    print("Loading complete")
    
    filename = "assembly_summary_genbank.txt"
    with open(filename,'r') as file:#read in files
        lines = file.readlines()
    
    genome_count = 0
    genome_count_b = 0
    genome_count_a = 0
    genome_count_e = 0
    genome_count_v = 0
    full_count = 0
    partial_count = 0
    sample_list = []
    genome_b = []
    genome_a = []
    genome_e = []
    sample_list_b = []
    sample_list_a = []
    sample_list_e = []
    taxid_list = [['a','b']]

        
      
    waiting = ['1'] * 8
    waiting[0] = "Processing |"
    waiting[1] = "Processing /"
    waiting[2] = "Processing -"
    waiting[3] = "Processing \\"
    waiting[4] = "Processing |"
    waiting[5] = "Processing /"
    waiting[6] = "Processing -"
    waiting[7] = "Processing \\"
    n_wait = 0
    passed = 0
    total = len(lines)
    #Read the data table
    print("Importing table")
    for line in lines:
        each_line = line.split("\t") #split the cells to read
        
        percent = genome_count / total * 100
        print ("Percentage: {:.4f}%".format(percent), end= '\r') #print dynamic progress report
        
        genome_rep = "None" #initialize genome replication status
        try:#attempt to acquire genome replication status
            genome_rep = each_line[13] 
        except:#if it is not a regular line, pass
            pass
        if genome_rep == "Full": #increment count
            full_count += 1
        if genome_rep == "Partial":#increment count
            partial_count += 1
            
        try:#investigating taxonomy lineage
            taxid = each_line[6] #get taxid
            lineage = ncbi.get_lineage(taxid) #get full lineage
            lineage = int(lineage[2]) #focus on the domain taxid
            if lineage == 2: #bacteria
                genome_count_b += 1
                genome_b += [taxid]
            elif lineage == 2157: #archaea
                genome_count_a += 1
                genome_a += [taxid]
            elif lineage == 2759: #eukaryota
                genome_count_e += 1
                genome_e += [taxid]
            else:
                genome_count_v += 1
            genome_count += 1
        except:
            pass
    
    print ("Percentage: {:.4f}%".format(100.0000)) #dynamic progress report end
    print("Table imported") #progress report
    print("===========================Result===========================") #print results
    print("There are", genome_count, "genomes in the databank, including: ")
    print(full_count, "fully replicated genomes and ", partial_count,"partially replicated genomes")
    print("Among", genome_count, "genomes,")
    print(genome_count_b, "are bacterial")
    print(genome_count_a, "are archaeal")
    print(genome_count_e, "are eukaryotic")
    print(genome_count_v, "are viral")
    print("============================================================") 
    #Ask user to process averaging
    need_input = True
    sample = False
    while need_input:
        user_input = input("Do you want to sample average genome length? (y/n): ")
        if user_input == "y":
            sample = True
            need_input = False
        if user_input == "n":
            need_input = False

    #start averaging if user demands
    if sample == True:
        sample_total = 0
        sample_total_b = 0
        sample_total_a = 0
        sample_total_e = 0
        try:
            sample_size_n = eval(input("Please enter the number of sample you need for averaging total genome size (default 500): "))
            sample_size = eval(input("Please enter the number of sample you need for averaging domain genome size (default 100): "))
        except:
            sample_size_n = 500
            sample_size = 100

        print("Averaging total genome length")
        #sample average length for overall genome

        for i in range (sample_size_n):
            percent = (i+1)/sample_size_n * 100
            print ("Percentage: {:.4f}%".format(percent), end= '\r')
            valid_length = False
            while valid_length == False: #the program will try to sample a valid length
                temp = random.randint(0, len(lines) -1)
                this_line = lines[temp].split('\t')
                taxid = this_line[6]
                try:#try to get the length from NCBI but the genome does not necessary have a length
                    length = findLength(taxid)
                    sample_total += length
                    valid_length = True
                except:
                    pass
        print ("Percentage: {:.4f}%".format(100.0000))        
        print("Averaging bacteria genome length")
        #Sample average length for bacteria
        for i in range (sample_size):
            percent = (i+1)/sample_size * 100
            print ("Percentage: {:.4f}%".format(percent), end= '\r')
            valid_length = False
            while valid_length == False: #the program will try to sample a valid length
                temp = random.randint(0, len(genome_b) -1)
                try:#try to get the length from NCBI but the genome does not necessary have a length
                    length = findLength(genome_b[temp])
                    sample_total_b += length
                    valid_length = True
                except:
                    pass

        print ("Percentage: {:.4f}%".format(100.0000))
        print("Averaging archaea genome length")
        #Sample average length for archaea
        for i in range (sample_size):
            percent = (i+1)/sample_size * 100
            print ("Percentage: {:.4f}%".format(percent), end= '\r')
            valid_length = False
            while valid_length == False: #the program will try to sample a valid length
                temp = random.randint(0, len(genome_a) -1)
                try:#try to get the length from NCBI but the genome does not necessary have a length
                    length = findLength(genome_a[temp])
                    sample_total_a += length
                    valid_length = True
                except:
                    pass
                
        print("Percentage: {:.4f}%".format(100.0000))
        print("Averaging eukaryota genome length")

        #Sample average length for eukaryota
        for i in range (sample_size):
            percent = (i+1)/sample_size * 100
            print ("Percentage: {:.4f}%".format(percent), end= '\r')
            valid_length = False #the program will try to sample a valid length
            while valid_length == False:
                temp = random.randint(0, len(genome_e) -1)
                try: #try to get the length from NCBI but the genome does not necessary have a length
                    length = findLength(genome_e[temp])
                    sample_total_e += length
                    valid_length = True
                except:
                    pass
        print ("Percentage: {:.4f}%".format(100.0000))

        #printing result
        print("===========================Result===========================")
        print("Average genome length: ", sample_total/sample_size)
        print("Average bacteria genome length: ", sample_total_b/sample_size)
        print("Average archaea genome length: ", sample_total_a/sample_size)
        print("Average eukaryota genome length: ", sample_total_e/sample_size)

    print("============================================================")
    #This is the genome downloader part
    #Ask whether user wants to download genome
    need_input = True
    need_download = False
    while need_input:
        user_input = input("Do you want to download genomes? (y/n): ")
        if user_input == "y":
            need_download = True
            need_input = False
        if user_input == "n":
            need_download = False
            need_input = False

    while need_download:
        #initialize search parameters
        taxid_search = 'x'
        genus_search = 'x'
        ref_search = 'x'
        n_search = 0
        asmbly_search = 'x'
        download_dir = 'x'

        #ask user whether use taxid or genus to search organisms
        need_input = True
        while need_input:
            user_input = input("Do you want to search with taxid or genus? (t/g/n): ")
            if user_input == 't' or user_input == 'g' or user_input == 'n':
                need_input = False
        #get the query
        if user_input == 't':
            taxid_search = input("Please enter the taxid you want to search: ")
        if user_input == 'g':
            genus_search = input("Please enter the genus you want to search: ")

        #ask user whether to limit to reference genome
        need_input = True
        while need_input:
            user_input = input("Do you want to lmited your search to reference genome, reference genome + representative genome, or no requirement? (r/h/n): ")
            if user_input == 'r' or user_input == 'h' or user_input == 'n':
                need_input = False
        #get the query
        ref_search = user_input

        #ask user how much genome to download
        need_input = True
        while need_input:
            user_input = input("How many genomes do you want to download: ")
            try:
                n_search = int(user_input)
                need_input = False
            except:
                pass
        #ask user whether to download all the data or only genomic assembly data
        need_input = True
        while need_input:
            user_input = input("Do you want to download all genome files (a) or only assembly files (s)? (a/s): ")
            if user_input == 'a' or user_input == 's':
                asmbly_search = user_input
                need_input = False

        #ask user to enter a deposite directory
        need_input = True
        while need_input:
            user_input = input("Please enter an existing directory for download, press enter to use current directory: ")
            if user_input == '':
                need_input = False
                download_dir = path
            elif os.path.isdir(user_input):
                need_input = False
                download_dir = user_input
            else:
                print("Invalid directory! Example directory: /mnt/c/Users/username/genome_data")

        print("Search parameter acquired, downloading...")
        error_genome = [] #initialize error genome

        os.chdir(download_dir)

        n_line = 0
        n_downloaded = 0
        for line in lines:
            met = True #initialize condition
            each_line = line.split("\t") #split the cells to read
            lineage = [] #initialize lineage
            if n_line == 0: #skip header
                met = False 
            if taxid_search != 'x':
                try:#investigating taxonomy lineage
                    lineage = ncbi.get_lineage(each_line[6]) #get full lineage
                    met = False
                    for taxid in lineage:
                        if int(taxid) == int(taxid_search):
                            met = True
                except:
                    met = False
            if met and genus_search != 'x':
                if genus_search not in each_line[7]:
                    met = False
            if met and ref_search != 'x':
                if ref_search == 'h':
                    if each_line[4] == "na":
                        met = False
                if ref_search == 'r':
                    if each_line[4] != "reference genome":
                        met = False
            #if all requirements are met, begin downloading
            if met:
                ftp = FTP('ftp.ncbi.nlm.nih.gov')
                ftp.login()
                print("==============Downloading genome", n_downloaded + 1, "==============")
                ftp_link = each_line[19]
                ftp_link = ftp_link[27:] #trim the ncbi home directory out of the link
                ftp_link_list = ftp_link.split('/')
                for each_level in ftp_link_list:
                    if not os.path.exists(each_level):
                        os.mkdir(each_level)
                    os.chdir(each_level)
                ftp.cwd(ftp_link)
                err = False
                filenames = ftp.nlst() # get filenames within the directory
                for filename in filenames:
                    ftp, err = checkNdownload(ftp, filename, os, asmbly_search, err)

                if err == True:
                    error_genome += each_line[0]
                
                for i in range (len(ftp_link_list)):
                    os.chdir("../")
                    ftp.cwd("../")
                ftp.quit()
                n_downloaded += 1
                if n_downloaded == n_search:
                    break
            n_line += 1

        if len(error_genome) == 0:
            print("Download successfuly and complete!")
        else:
            print("Download complete, but error had occured in genome: ")
            for genome in error_genome:
                print(genome)
        need_input = True
        user_input = input("Do you want to download other genomes? (y/n): ")
        while need_input:
            if user_input == "y":
                need_download = True
                need_input = False
            if user_input == "n":
                need_download = False
                need_input = False
                

    


        
        
        


