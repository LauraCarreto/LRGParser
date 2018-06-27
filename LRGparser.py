"""
LRGparser.py

Created: 14 December 2016

Developed and tested on python versions 3.5.2 and 2.7.8, in a Linux/Unix environment.
This script will run with python versions 3.5 and 2.7.

@authors: Laura Carreto, Rosie Coates-Brown

usage: python LRGparser.py -g [LRG file name] -d [True/False] -a [True/False] -s [file/url] -p [path to working directory]

required parameters:
-g, --gene    [name of LRG file without .xml suffix]

optional parameters:
-h, --help shows this message and quits
-d, --diff, --difference = [True/False] triggers or supresses the output of [LRG]_diffs.csv
-a, --annot, --annotations = [True/False] triggers or suppresses [LRG]_annotation.csv
-s, --source = [url/file] default is from file. Adding -s url will trigger LRGparser.py to grab the xml from http://ftp.ebi.ac.uk/pub/databases/lrgex 
-p, --path =  defaults to cwd if no path is supplied. Path to your locally downloaded LRG files and list_LRGs_GRCh38.txt file. This can be downloaded from http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_GRCh38.txt

NB if you are using locally saved LRG files and list_LRGs_GRCh38.txt these must be in the same location.

output:
[LRG]_t1.bed: a tab separated bed file containing the chromosome number, exon start position, exon end position  
[LRG]_diffs.csv: a csv file containing the differences between 37 and 38; -d False will suppress this.
[LRG]_annotation.csv: a csv file of gene information including synonyms, lsdb, long gene name; -i False will suppress this.

"""

try:
    import xml.etree.ElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import sys, os, csv, getopt


def check_status(in_opt, refpath, version, genein):
    """
    reads in the list_LRGs_GRCh38.txt file from file or URL and checks the status of the LRG file is public. 
    LRGparser.py doesn't deal with pending files for data integrity
    
    parameters: in_opt (str from command line), version (str), (genein (str from command line))
    returns: status (str)
    """

    gene = genein

    # create dictionary of LRG = status
    LRG_status = {}
    
    # open file with LRG status and add key=value pairs to LRG_status dictionary
    #a flow control is required to deal with the difference in modules required for python 2.7 and python 3.5    
    if in_opt == "url":
        print ("Using list_LRGs_GRCh38.txt from: http://ftp.ebi.ac.uk/pub/databases/lrgex/")
        
        if version == "3.5":
            from urllib.request import urlretrieve
            status_filename, headers = urlretrieve('http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_GRCh38.txt') 
            
            with open(status_filename, 'r') as status_file:
                for line in status_file:
                    split_line = line.split()
                        
                    if split_line[2] == 'modified:':             
                        last_modified = split_line[3]
                        next
                    elif split_line[1] == 'LRG_ID': # line with headers
                        next
                    else:
                        LRG_status[split_line[0]] = split_line[2]
                    
                    
        elif version == "2.7":
            import urllib2
            status_filename = urllib2.urlopen('http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_GRCh38.txt')
        
            for line in status_filename:
                split_line = line.split()
                        
                if split_line[2] == 'modified:':             
                    last_modified = split_line[3]
                    next
                elif split_line[1] == 'LRG_ID': # line with headers
                    next
                else:
                    LRG_status[split_line[0]] = split_line[2]
               

    elif in_opt == "file":
        status_filename = refpath+'/list_LRGs_GRCh38.txt'
        print ('Using status file:'+status_filename)
        try:
            with open(status_filename, 'r') as status_file:
                for line in status_file:
                    split_line = line.split()
                        
                    if split_line[2] == 'modified:':             
                        last_modified = split_line[3]
                        next
                    elif split_line[1] == 'LRG_ID': # line with headers
                        next
                    else:
                        LRG_status[split_line[0]] = split_line[2]
        except:
            print("Couldn't open file... check if you have a trailing / in your reference filepath; remove it, and try again.")
            usage()
            sys.exit(2)
                    
                
    # inform user if LRG status is public or pending; exit parser if the latter occurs
    status = LRG_status.get(gene, None)
    if status is None:
        print ("No such LRG file exists")
        exit(0)
        
    print ("LRG status for " + gene + " : " + status)
    print ("(last status update: " + last_modified + ")")
    
    if status == 'pending':
        print ("Sorry! LRGparser does not process LRGs with 'pending' status.\n")
        print ("To check the last status update, please go to:")
        print ("http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_GRCh38.txt"+"\n")
        exit (0)

    return status

def read_file(genein, in_opt, refpath, version):
    """
    read in the LRG.xml file
    test: does the file exist

    parameters: genein(string returned from command line arguments)
    returns: root(ElementTree root node), gene(user input string variable)

    """
    gene = genein
    
    file_name = gene+'.xml'
    #file_path = '/home/swc/Desktop/LRGParser/Assignment/'
    #file_path = '/Users/rosiecoates/Documents/Clinical_bioinformatics_MSc/programming/assignment/'
    file_path = refpath+'/'
    full_path = file_path+file_name
    
    #get LRG from reference file location
    if in_opt == "file":
        # check if LRG file exists
        try:
            tree = ET.parse(full_path)
        except:
            print("Couldn't open file... check you have supplied an LRG file name without the .xml extension. Check if you have a trailing / in your reference filepath; remove it, and try again.")
            usage()
            sys.exit(2)
    
        tree = ET.ElementTree(file=full_path)
        
    #Get LRG file form URL (again, flow control required due to different modules required to open urls for py2 and py3)   
    elif in_opt == "url":
        print ('Requesting url: '+ 'http://ftp.ebi.ac.uk/pub/databases/lrgex/'+ file_name)
        if version == "2.7":
            import urllib2
            try: # parse xml and create root object
                tree = ET.ElementTree(file=urllib2.urlopen('http://ftp.ebi.ac.uk/pub/databases/lrgex/'+ file_name))
        
            except:
                print("Couldn't open URL with LRG file... check internet provision and that you have supplied a valid LRG gene name.")
                usage()
                sys.exit(2)
                
        elif version == "3.5":
            from urllib.request import urlopen
            try: # parse xml and create root object
                with urlopen('http://ftp.ebi.ac.uk/pub/databases/lrgex/'+ file_name) as response:    
                    tree = ET.ElementTree(file=response)
            except:
                print("Couldn't open URL with LRG file... check internet provision and that you have supplied a valid LRG gene name.")
                usage()
                sys.exit(2)
    
    root = tree.getroot()
    
    return root, gene, file_name
    
def write_csv(mylist, myfilename, mode):
    """
    Output to CSV file from a list
    
    Parameters: mylist (list), myfilename(string), mode (string); the mode options used in the code are: 'a'= append, 'w'=write
    
    """
    with open(myfilename, mode) as csvfile:
        out = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
        out.writerow(mylist)
    csvfile.close()
        
    return


def bed_file(root, gene):
    """
    produces a bed file containing the chromosome, genomic start position and genomic end position

    test: is start position bigger than end position?
    test: is the calculate chromosome position between the given chromosome start and end?

    parameters: root(ElementTree root node), gene(user input string variable)

    returns: exon_ranges(dict)
    """
    # initialise dictionary of exon ranges; format: exon = list(exon_start, exon_end)    
    exon_ranges={}
    strand=''
    
    for mapping in root.findall("./updatable_annotation/annotation_set[@type='lrg']/mapping[@type='main_assembly']"):
        # for each branch, get start and end coordinates in the reference genome and chromosome number 
    
        chr = mapping.attrib['other_name'] # chromosome
        ref_start = int(mapping.attrib['other_start']) # start (converted to integer)
        ref_end = int(mapping.attrib['other_end']) # end (converted to integer)
        
        # account for + or - strand when converting exon lrg coordinates to genomic coordinates
        for mapping_span in mapping:  
            strand = mapping_span.attrib['strand']
            print (gene+' is coded in strand '+ strand)
        
        if strand == '1': # CONFIRM OFFSET
            offset = (ref_start)     # offset genomic start 
            offset_int = int(offset) # convert offset start to integer
            
        elif strand == '-1': # CONFIRM OFFSET 
            offset = (ref_end)        # offset end 
            offset_int = int(offset)  # convert offset start to integer

    for id in root.findall("./fixed_annotation/id"):
        id_tag = id.text

    for transcript in root.findall("./fixed_annotation/transcript"):
        # for each transcript, get exon coordinates into bed file
        trans_num =  (transcript.attrib['name']) # transcript number
        ref_name = gene + "_" + trans_num + ".bed" # define bed file name
        
        bedfile = open(ref_name, 'w') # open bed file
        
        for exon in transcript:
            #for each exon, get start and end coordinates as integers; offset lrg coordinates to get genomic coordinates 
            if exon.tag == 'exon':
                exon_number =  (exon.attrib['label'])

                for coordinates in exon:
                    if (coordinates.attrib['coord_system']) == id_tag:
                        start=(coordinates.attrib['start'])
                        start_int = int(start)
                        end = coordinates.attrib['end']
                        end_int = int(end)
                        
                        exon_ranges[exon_number]=[start_int,end_int]

                        if strand == '1': # genomic coordinates calculations
                            gen_start = str(start_int + offset_int)
                            gen_start_int = int(gen_start)
                            gen_end = str(end_int + offset_int + 1) 
                            gen_end_int = int(gen_end)                     
                             
                        if strand == '-1': # genomic coordinates calculations
                            gen_start = str(offset_int - end_int - 1)
                            gen_start_int = int(gen_start)
                            gen_end = str(offset_int - start_int)
                            gen_end_int = int(gen_end)
                        
                        # assert that exon genomic coordinates are within genomic start and end coordinates
                        assert (gen_start_int >= ref_start or gen_start_int < ref_end), "ERROR: calculated genomic position of exon start is not within given genomic coordinates"
                        assert (gen_end_int > ref_start or gen_end_int <= ref_end), "ERROR: calculated genomic position of exon end is not within given genomic coordinates"
                        # assert that genomic coordinates for exon start are below genomic coordinates for exon end
                        assert (gen_start_int < gen_end_int), "ERROR: calculated genomic start position is above genomic end position"
                        
                        bed_list = [chr, gen_start, gen_end]
                        bedfile.write("\t".join(bed_list))
                        bedfile.write("\n")
                            
                   

    return exon_ranges, strand


def get_diffs(exon_ranges, gene, root):
    """
    produces a csv file of the differences in the gene between build 37 and build 38

    parameters: exon_ranges(dict), gene(string), root(elementTree root node)


    """
    # initialise data structures
    diffexons={} # dictionary linking differences to exon number; value defaults to 'intronic' if position do not map within exon coordinates
    lrgstartlist=[] # list of differences' start position
    
    # define file name for differences file and create list of column headers to write to file    
    diff_file = gene + "_diffs.csv"
    diff_headers = ["position", "type", "lrg_start", "lrg_end", "other_start", "other_end", "LRG_seq", "other_seq"]

    # Check if differences exist with respect to assembly GRCh38.p7
    diff_count = 0
    for diffs in root.findall("./updatable_annotation/annotation_set[@type='lrg']/mapping[@type='main_assembly']/mapping_span/diff"):
        diff_count += 1
    
    if (diff_count == 0):
        print ("No differences exist with respect to assembly GRCh38.p7, therefore no differences file was generated.\n")
        return
    else:
        pass

    # write headers into csv file
    write_csv(diff_headers, diff_file, 'w') # mode 'w' truncates any file with same name in the directory

    for mapping in root.findall("./updatable_annotation/annotation_set[@type='lrg']/mapping[@type='main_assembly']"):
        # get reference assembly id and print to console
        ref_assem = (mapping.attrib['coord_system'])
        print ("Reference assembly= "+ref_assem +" so differences are with respect to this build.\n")
        
        # for each difference in xml, get attributes into list
        for span in mapping:
            for diff in span:
                lrg_start = int(diff.attrib['lrg_start'])
                lrgstartlist.append(lrg_start) #for each difference, append to list of start positions
                
                typeattrib = diff.attrib['type']
                lrg_start_str = (diff.attrib['lrg_start'])
                lrg_end = (diff.attrib['lrg_end'])
                other_start = (diff.attrib['other_start'])
                other_end = (diff.attrib['other_end'])
                LRG_seq = diff.attrib['lrg_sequence']
                other_seq = diff.attrib['other_sequence']

                for pos in lrgstartlist:
                    # check if start position of difference is within exon start and end coordinates
                    # using exon = (exon_start, exon_end) dictionary; values are a list;
                    for key, value in exon_ranges.items(): # for key, value pair in exon_ranges dictionary                     
                        # for exon in exon_ranges dictionary
                        # compare pos to exon_start and exon_end coordinates in all exons                       
                        if pos >= value[0] and pos <= value[1]: # value [0] refers to exon_start, value[1] refers to exon_end
                        # if pos within exon range,
                        # add pos = exon to diffexons dictionary
                            diffexons[pos] = 'exon '+str(key)
                            
                        else:
                            next

                    # if pos not in exons, it is not in diffexons yet
                    # must be added to diffexons as 'intronic'
                    if pos not in diffexons.keys():
                        diffexons[pos] = 'intronic'

                    for position, location in diffexons.items():
                        # for position in diffexons, location is either the exon number or 'intronic'
                        if position == lrg_start:
                            # set up list with attributes for k difference
                            diff_list = [location, typeattrib, lrg_start_str, lrg_end, other_start, other_end, LRG_seq, other_seq]
                            # write diff_list content to line in csv file                          
                            write_csv(diff_list, diff_file, 'a') #mode 'a' to append to existing file diff_file
    return diff_list

def get_annotations(gene, root):    
   """
    Outputs gene annotations, including overlapping genes and respective gene name synonyms, to CSV file
    
    parameters: gene(string), root(elementTree root node)
    
   """
   #initialise file with headers
   annot_file = gene+"_annotation.csv"
   annot_headers = ['NCBI_ID','HGNC_symbol', 'LRG_start','LRG_end','Strand','Description','Synonyms' ]
   write_csv (annot_headers, annot_file, 'w')#write mode ('w') truncates file with same name in directory to avoid appending to an old file

   # loop through LRG file to get annotations into annotation_list
   # one list for each overlapping gene (if present)
   for gene in root.findall("./updatable_annotation/annotation_set[@type='ncbi']/features/gene"):
        annotation_list=[]
        
        if gene.attrib.get('source')=='NCBI-Gene':
            NCBI = gene.attrib['accession']
            annotation_list.append(NCBI)
        else: #in case no NCBI accession
            annotation_list.append('')
            
        for branch in gene:
            if branch.tag == 'symbol':
                    
                if branch.attrib.get('source')=='HGNC':
                    HGNC = branch.attrib['name']
                    annotation_list.append(HGNC)
                    
                    #create a synonyms list
                    synonym_list=[]
                    for synonym in branch:
                        synonym_list.append(synonym.text)
                        
                    synonym_string = ', '.join(synonym_list) #convert synonym_list into string to append later into annotation_list as one element
                   
                    
                    
                else:#in case no HGNC name
                    annotation_list.append('')
                    synonym_string = ''

                    
            if branch.tag == 'coordinates':
                LRG_start = branch.attrib['start']
                LRG_end = branch.attrib['end']
                Strand = branch.attrib['strand']
                                
                annotation_list.append(LRG_start)
                annotation_list.append(LRG_end)
                annotation_list.append(Strand)
                                
            if branch.tag == 'long_name':
                ln = branch.text
                annotation_list.append(ln)
                                    
                annotation_list.append(synonym_string)
                                    
                write_csv (annotation_list, annot_file, 'a') #mode 'a' to append to annot_file
                
        return ln

def versiontest():
    versionbool = ''
    version =".".join(map(str, sys.version_info[:2]))
    if version in ["3.5", "2.7"]:
        versionbool = 0
        print ("Goodnews! you are using a compatible version of python: " + version)
    else:
        veresionbool = 1
        print ("Friendly warning: LRGparser has not been tested on this python version: " + version)

    return versionbool, version

    
    
def main():
    """
    Parses and handles command line arguments. Calls the other functions in the script.
    Functions produce a bedfile, an optional annotation file and an optional information file
    """
    bool, version = versiontest()
    
    #parses the command line arguments to check that all flags passed are valid, exits if not
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hg:d:a:s:p:', ['help', 'gene=', 'difference=', 'info='])
    except getopt.GetoptError as err:
        print (err)
        usage()
        sys.exit(2)
    #defines all parameters to allow the possibility of optional arguments
    refpath = os.getcwd()
    genein = ''
    diff = "false"
    info = "false"
    #sets input file source to reference file as default
    in_opt = "file"
    for opt, arg in opts:
        if opt == '-h':
            print (__doc__)
            sys.exit(2)
        if opt == '--help':
            print (__doc__)
            sys.exit(2)
            
        elif opt in ('-g', '--gene'):
            genein = arg
        elif opt in ('-d', '--difference', '-diff'):
            diff = arg
        elif opt in ('-a', '--annotation', '--annot'):
            info = arg
        elif opt in ('-s', '--source'):
            in_opt = arg
        elif opt in ('-p', '--path'):
            refpath = arg
            print ('Using reference file path:'+refpath)
        else:
            usage()
            sys.exit(2)
    #checks if an LRG file name has been passed, exits if not
    if genein == '':
        print ('Please supply LRG file name without extension')
        usage()
        sys.exit(2)
    
    status = check_status(in_opt, refpath, version, genein)
      
    # read xml; function returns root object and variable with gene name
    root, gene, filename = read_file(genein, in_opt, refpath, version)
    # create bed file and return dictionary of exon ranges
    exon_ranges, strand = bed_file(root, gene)

    # create csv file with sequence differences
    if diff == "True":
        diff_list = get_diffs(exon_ranges, gene, root)
        print ("-d = ", diff, " therefore differences file requested. The file may not have been produced, if there are no differences to report.")
    else:
        print ("-d = ", diff, " therefore no differences file produced.")

    # create csv file with annotations for overlaping genes and respective synonyms    
    if info == "True":
        last_ln = get_annotations(gene, root)
        print ("-a = ", info, " therefore annotation file produced.")
    else:
        print ("-a = ",  info, "therefore no annotation file produced.")
        
    
def usage():
    """
    helpful hints about the usage of the script
    """

    print ("usage:")
    print ("python LRGparser.py -g [LRG file name] -d [True/False] -a [True/False] -s [file/url] -p [path to working directory]")
    print ("NB no trailing / is required when providing a file path to locally saved LRG files.")    
#runs script if it is run as a script from the command line as opposed to as a function
if __name__ == "__main__":
    main()





