#parse through harrington_clinical_data and print to the screen loop through each line

import sys
print(sys.executable)


import argparse
import os
import subprocess
import shutil
import pysam 


class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""

    def __init__(self, filePath, headerSymbols=['@', '+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...

        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'r')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols

    def __iter__(self):
        return self

    def __next__(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1  ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else:
                elemList.append(None)

        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4, \
            "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
            Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (
                self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]), \
            "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
            Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (
            self._hdSyms[0], self._currentLineNumber)
        assert elemList[2].startswith(self._hdSyms[1]), \
            "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
            Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (
            self._hdSyms[1], self._currentLineNumber)
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (
            self._currentLineNumber)

        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)

#### MY CODE BELOW ######

def get_clinical_data(file_path):

    #go through the harrington_clinical_data.txt file and return a dictionary
    clinical_data = {}
    with open(file_path, 'r') as file:
        #grab headers
        headers = file.readline().strip().split('\t')
        for line in file:
            values = line.strip().split('\t')
            data = dict(zip(headers, values))
            #the code name the key- pull out the value 
            code_name = data.pop("Barcode")
            clinical_data[code_name] = data
    return clinical_data

def trim_ends(read):
    # Initialize a variable to store the index of the start of the degraded end
    start_index = None
    # Find the index of the first occurrence of at least two consecutive 'D' or 'F' from the end
    for i in range(len(read)):
        if i > 0  and read[i] in ['D', 'F'] and read[i - 1] in ['D', 'F']:
            start_index = i - 1
            break
    # If degraded end found, trim the read
    if start_index is not None:
        trimmed_read = read[:start_index]
    else:
        trimmed_read = read
    return trimmed_read

def align_fastq_to_reference(reference_file, fastq_file, output_folder):
    # Check if the reference file is indexed, if not, index it
    # if not os.path.exists(reference_file + ".bwt"):
    #     subprocess.run(["bwa", "index", reference_file])


    # Extract the name of the fastq file without the extension
    fastq_name = os.path.splitext(os.path.basename(fastq_file))[0]

    # Define the name for the output SAM file
    sam_file = "./"+output_folder+"/"+fastq_name + ".sam"

    # Perform alignment using BWA mem
    # subprocess.run(["bwa", "mem", reference_file, fastq_file], stdout=open(sam_file, "w"))
    with open(sam_file, 'w') as output_file:
        subprocess.run(["bwa", 'mem', reference_file, fastq_file], stdout=output_file)

def build_sam_files(fastq_directory,output_folder):
    # Path to the reference file
    reference_file = "dgorgon_reference.fa"

    # Iterate over each FASTQ file in the directory
    for fastq_file in os.listdir(fastq_directory):
        if fastq_file.endswith(".fastq"):
            # Construct the full path to the FASTQ file
            fastq_path = os.path.join(fastq_directory, fastq_file)

            # Align the FASTQ file to the reference sequence
            align_fastq_to_reference(reference_file,fastq_path,output_folder)

def parse_fastqs(folder_name, clinical_data,fastqParsArgsString):

    # print(clinical_data["TATGG"])
    fastqfile = ParseFastQ(fastqParsArgsString)


    #use the "barcode" column 3 to pasre through the hawkins_pooled_sequence loop through each line
    #1a remove the first five nucleotides
    #1b remove reads that have DD, FF, FD, or DF at the end of the read
    #output 50 fastq files {name}_trimmed.fastq (ie trim_trimmed.Fastq)

    #Grab one fastq sequence at a time and crunch it.
    for fastq_obj in fastqfile:
        header = fastq_obj[0]

        sequence = fastq_obj[1]

        #Get the Barcode from the first 5 characters
        barcode = sequence[:5]
        if barcode in clinical_data:

           #Get the s equence from the rest of the string (eating the barcode)
            sequence = sequence[5:]
            # This is the separator
            separator = fastq_obj[2]

            # This is the quality score
            qualScore = trim_ends(fastq_obj[3])
            
            #Get the s equence from the rest of the string (eating the barcode) 
            sequence = sequence[:len(qualScore)]

            qualScore = qualScore[:]
            # sequence = sequence[:len(qualScore)]
            #save the sequence
            string_to_write = header + "\n" + sequence + "\n" + separator + "\n" + qualScore + "\n"

            file_name = "./" + folder_name + "/" + clinical_data[barcode]["Name"] + "_trimmed.fastq"
            with open(file_name, 'a') as f:
                f.write(string_to_write)

def make_folder(folder_name):
    delete_folder(folder_name)

    os.makedirs(folder_name)

def delete_folder(folder_name):
    if os.path.exists(folder_name):
        # Folder exists, delete it
        shutil.rmtree(folder_name)

def convert_to_bam(sam_file, output_folder):

    # Extract the name of the fastq file without the extension
    name = os.path.splitext(os.path.basename(sam_file))[0]

    # Define the name for the output SAM file
    bam_file = "./"+output_folder+"/"+name + ".bam"

    print(bam_file)
    print("./"+sam_file)

    
    # sam_file = "./sam_files/Abbey_trimmed.sam"
    # bam_file = "./bam_files/test.bam"
    with open(bam_file, 'wb') as output_file:
        subprocess.run(f"samtools view -bS {sam_file} > {bam_file}", shell=True)
    # exit(0)

    # Maybe for the server Save for now
    # with open(sam_file, 'wb') as output_file:
    #     subprocess.run(["samtools", "view", "-bS", sam_file], stdout=output_file)
    # exit(0)


def convert_to_bam_sorted(input_bam, output_folder):

     # Extract the name of the fastq file without the extension
    name = os.path.splitext(os.path.basename(input_bam))[0]
    name = name[:len(name)-8]
    # Define the name for the output SAM file
    # output_bam = "./"+output_folder+"/"+ name + ".bam"
    output_sorted_bam = os.path.join(output_folder, name + ".sorted.bam")
    name = os.path.splitext(os.path.basename(output_sorted_bam))[0]

    # with open(output_sorted_bam, 'wb') as output_file:
    subprocess.run(["samtools", "sort", "-m", "100M", "-o", output_sorted_bam, input_bam])

    sorted_bam_name = output_sorted_bam.split('.')[0]  # Extracting the base name without extension

    # Command to index the sorted BAM file
    index_command = ["samtools", "index", f"{sorted_bam_name}.sorted.bam"]

def build_bam_files(sam_folder,output_folder):

    for fastq_file in os.listdir(sam_folder):
        if fastq_file.endswith(".sam"):
            # Construct the full path to the FASTQ file
            fastq_path = os.path.join(sam_folder, fastq_file)

            # Align the FASTQ file to the reference sequence
            convert_to_bam(fastq_path,output_folder)

def build_sorted_bam_files(sam_folder, output_folder):

    for fastq_file in os.listdir(sam_folder):
        if fastq_file.endswith(".bam"):
            # Construct the full path to the FASTQ file
            fastq_path = os.path.join(sam_folder, fastq_file)

            # Align the FASTQ file to the reference sequence
            convert_to_bam_sorted(fastq_path,output_folder)

def build_pileup_bam_files(sam_folder):
    data=[]
    for fastq_file in os.listdir(sam_folder):
        # if fastq_file.endswith(".sam"):
        # Construct the full path to the FASTQ file
        fastq_path = os.path.join(sam_folder, fastq_file)

        # Align the FASTQ file to the reference sequence
        data.append(pileup(fastq_path))
    return data



def pileup(samfile_name):
    positions = []
    #test file, replaced with the sorted.bam you are using. Make sure it is indexed! (Use samtools index yourbam.sorted.bam)
    samfile = pysam.AlignmentFile(samfile_name,"rb")

    #Since our reference only has a single sequence, we're going to pile up ALL of the reads. Usually you would do it in a specific region (such as chromosome 1, position 1023 to 1050 for example)
    for pileupcolumn in samfile.pileup():
        # print ("coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
        #use a dictionary to count up the bases at each position
        ntdict = {}
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # You can uncomment the below line to see what is happening in the pileup. 
                # print('\tbase in read %s = %s' % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]))
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                ########## ADD ADDITIONAL CODE HERE #############

                if base in ntdict:
                    ntdict[base] += 1
                else:
                    ntdict[base] = 1
                # Populate the ntdict with the counts of each base 
                # This dictionary will hold all of the base read counts per nucletoide per position.
                # Use the dictionary to calculate the frequency of each site, and report it if if the frequency is NOT  100% / 0%. 
                #############################################

        positions.append(ntdict)
    samfile.close()

    temp = {}
    temp["file"] = samfile_name[5:len(samfile_name)-11]
    temp["pos"] = positions
    return temp

def get_color_from_clinical_data(clinical_data, name):
    for barcode, data in clinical_data.items():
        if data["Name"] == name:
            return data["Color"]

def generate_report(reportData,clinical_data,file_name):
    with open(file_name, 'w') as f:
        for report in reportData:
            pos = 0
            for item in report["pos"]:
                if len(item) > 1:
                    # Get percentage
                    baseCount = 0
                    baseSum = 0
                    minBase = ""
                    for base, count in item.items():
                        baseSum += count
                        if baseCount == 0 or count < baseCount:
                            minBase = base
                            baseCount = count
                    percentage = baseCount/baseSum*100
                    color = get_color_from_clinical_data(clinical_data,report["file"])

                    f.write("Sample " + report["file"] + " had a " + color + " mold, " + str(baseSum) + " reads, and had " + str(int(percentage)) + "% of the reads at position " + str(pos) + " had the mutation " + minBase + ".\n\n")
                    # print("Position " + str(pos))
                    # print(item)
                pos = pos + 1
def main():

    fastq_folder_name = "fastqs"
    sam_folder_name = "sam_files"
    bam_folder_name = "bam_files"
    sorted_bam_folder_name = "bams"
    fastqParsArgsString = "./hawkins_pooled_sequences.fastq"

    harrington_clinical_data_file = './harrington_clinical_data.txt'
    clinical_data = get_clinical_data(harrington_clinical_data_file)


    make_folder(fastq_folder_name)
    parse_fastqs(fastq_folder_name,clinical_data, fastqParsArgsString)
    
    make_folder(sam_folder_name)
    build_sam_files(fastq_folder_name, sam_folder_name)

    make_folder(bam_folder_name)
    build_bam_files(sam_folder_name, bam_folder_name)

    make_folder(sorted_bam_folder_name)
    build_sorted_bam_files(bam_folder_name, sorted_bam_folder_name)

    delete_folder(bam_folder_name)
    delete_folder(sam_folder_name)

    reportData = build_pileup_bam_files(sorted_bam_folder_name)

    generate_report(reportData,clinical_data,"report.txt")

main()
    


