#! /usr/bin/env python

import zipfile
import re
import matplotlib.pyplot as plt
import os
import glob
import subprocess
import sys
import numpy as np
import PrettyTable from prettytable
import wcwidth
from datetime import datetime
import exml.etree.ElementTree as et
import json
import xmltodict

### Please ignore these comments (these were system specific bugs, but leaving here in case any errors arise and it's helpful)

# Add my modules dir to the system path so python finds it
# I don't like this, but I don't have permissions to download these modules where pythonpath is set to search for modules within the cluster
#sys.path.append('/home/sam/bin/python_modules/prettytable/src/prettytable')
#sys.path.append('/home/sam/bin/python_modules/wcwidth/wcwidth')
#sys.path.append('/home/sam/bin/python_modules/xmltodict')

# import the modules I need -- NOTE: these must be in this order due to reference to parent/local directories in wcwdith script. I edited these to just refer to the modules without relation to one another
# if this present problems, simply import prettytable and wcwidth above as you normally would
#import table_wide
#import table_zero
#import unicode_versions
#from prettytable import PrettyTable
#import wcwidth
#import xmltodict

# MODULE ERROR - PRINT PATH WHERE MODULES ARE LOOKED FOR
#print("sys.path:\n" + "\n".join(sys.path))
#exit()

###

'''
Author: Sam Czerski
Generation of Report following Microbiome Analysis Pipeline (CCS & DEMUX)
Must be run from the run_id subdirectory in the SCRATCH drive with conventionally named subfolders (eg. 1_A01, 2_B01, 3_C01/ccs/outputs/, 4_D01/demux_no_peek/outputs/)
Last updated: 03/18/2022
'''


def check_smrt_cells_and_ccs_files():
    print("Confirming Existence of CCS Files...")

    #Creating list of possible smrt cell directories to search for
    smrt_cells = ["1_A01", "2_B01", "3_C01", "4_D01"]
    smrt_cells_in_working_dir = []

    #First, see what smrt cells are present and add them to a new list to work from
    for cell in smrt_cells:
        is_dir_here = os.path.exists("%s" % cell)
        if is_dir_here == True:
            smrt_cells_in_working_dir.append(cell)
            print("%s found" % cell)
        else:
            print("%s not found" % cell)

    #Now, we have a list of the smrt cells within the working directory.
    #Next, we will confirm the existence of files with the CCS data we are interested in.

    for present_cell in smrt_cells_in_working_dir:
        os.chdir(f'{present_cell}/ccs/outputs')

    #if ccs.report.zip must be unzipped, try unzipping it.
        try:
            with zipfile.ZipFile("ccs.report.csv.zip") as z:
                z.extractall()
                print("\nExtracted Output CCS Statistics File Successfully.")
        #if ccs.report.csv.zip is already unzipped, try to open ccs_statistics.csv file (unzipped from report.csv.zip)
        except:
            try:
                with open("ccs.report.csv"):
                    print("Located CCS Summary Output File Successfully.")
            except:
                #if this file is not found, stop.
                raise FileNotFoundError("ccs.report.csv.zip not found.")

        finally:
            is_passes_tsv_here = os.path.exists("ccs_passes.tsv")
            if is_passes_tsv_here == True:
                print("ccs_passes.tsv File Found.\n")
                print("All %s CCS Files Found.\n" % present_cell)
                os.chdir("..")
                os.chdir("..")
                os.chdir("..")
                
            else:
                raise FileNotFoundError("ccs_passes.tsv not found.")

          # better way to do this...? seems a little risky, but should be fine as long as files are named conventionally


    return smrt_cells_in_working_dir

def get_ccs_summary(smrt_cells):

    # create new output statistics report file
    with open("ccs_demux_summary_report.txt", "a") as fout:
        # write date and time to the top of the output file
        now = datetime.now()
        # reformat
        dt_string = now.strftime("%Y_%m_%d %H:%M:%S")
        fout.write(str("Date and Time: "))
        fout.write(str(dt_string))
        fout.write(str("\n"))
        # Begin Report - write title
        fout.write(str('\nCCS Sub-Report\n'))

        #initialize list for comparing read counts with polymerase reads
        total_ccs_counts_all_cells = []
        #initialize variable for keeping track of num of iterations to keep list accurate
        iteration = 0
        
        # report stats for all cells
        for cell in smrt_cells:
            fout.write(str(f'\nSMRT CELL: {cell}\n'))
            # initializing variables
            read_length = []
            count = 0
            read_sum = 0
            pattern = re.compile(r',(\d+),')
        

            # change into directory with existing CCS files
            os.chdir(f'{cell}/ccs/outputs/')
    
            # obtain/generate summary data and statistics
            with open("ccs.report.csv") as f:
                # skip the header line
                next(f)

                for line in f:
                    count += 1
                    qlen = pattern.findall(line)
                    res = [int(i) for i in qlen]
                    read_length.extend(res)

                for num in read_length:
                    read_sum += num

                # Append count to list  returneded with totals for all cells
                if iteration == 0 and cell == "2_B01":
                    total_ccs_counts_all_cells.append(0)
                    total_ccs_counts_all_cells.append(count)
                    iteration += 1

                elif iteration == 0 and cell == "3_C01":
                    total_ccs_counts_all_cells.append(0)
                    total_ccs_counts_all_cells.append(0)
                    total_ccs_counts_all_cells.append(count)
                    iteration += 1

                elif iteration == 0 and cell == "4_D01":
                    total_ccs_counts_all_cells.append(0)
                    total_ccs_counts_all_cells.append(0)
                    total_ccs_counts_all_cells.append(0)
                    total_ccs_counts_all_cells.append(count)
                    iteration += 1

                else:
                    total_ccs_counts_all_cells.append(count)


                # Summary Statistics

                #fout.write("Total number of reads:\n")
                #fout.write(str(count))
                #fout.write("\nMin. Read Length:\n")
                min_read_length = min(read_length)
                #fout.write(str(min_read_length))
                max_read_length = max(read_length)
                #fout.write("\nMax. Read Length:\n")
                #fout.write(str(max_read_length))
                avg_read_length = round(read_sum / count)
                
                #fout.write("\nAverage Read Length:\n")
                #fout.write(str(round(avg_read_length)))
                fout.write("\nHistogram of Average Read Length Saved to png file\n")
                plt.hist(np.log(read_length))
                plt.xlim([np.log(min_read_length),np.log(max_read_length)])
                X_ticks = np.arange(np.log(min_read_length),np.log(max_read_length))
                plt.xticks(X_ticks)
                plt.title(f'{cell}: Post-CCS Log-Scaled Read Length')
                plt.xlabel("Log-Scaled Lengths of Reads")
                plt.ylabel("Frequency")
                plt.savefig(f'{cell}_avg_read_length_hist.png')
                plt.clf()
                
                # Make table of Summary Stats
                summary_stats_headers = ["Avg. Read Length", "Min. Read Length", "Max. Read Length", "Total Num. of Reads"]
                summary_stats_table = PrettyTable(summary_stats_headers)
                summary_stats_table.add_row([avg_read_length, min_read_length, max_read_length, count])
                # Write table to output file
                fout.write(str(summary_stats_table))

            # CCS Passes Stats

            # Initiate variables
            num_passes = []
            passes_sum = 0
            pcount = 0
            pattern = re.compile(r'\t(\d+)')

            # Open passes file and compute average
            with open("ccs_passes.tsv") as f:

                for line in f:
                    pcount += 1
                    passes = pattern.findall(line)
                    res = [int(i) for i in passes]
                    num_passes.extend(res)

                for num in num_passes:
                    passes_sum += num

                # Summary Statistics

                #fout.write("\nAverage Number of CCS Passes:\n")
                avg_num_ccs_passes = round(passes_sum / pcount)
                min_num_ccs_passes = min(num_passes)
                max_num_ccs_passes = max(num_passes)
                #fout.write(str(round(avg_num_ccs_passes)))
                fout.write("\nHistogram of CCS Passes Saved to png file\n")
                plt.hist(num_passes)
                plt.xlim([0,100])
                X_ticks = np.arange(0,100,5)
                plt.xticks(X_ticks)
                plt.xlabel("Number of CCS Passes")
                plt.ylabel("Frequency")
                plt.title(f'{cell}: CCS Passes')
                plt.savefig(f'{cell}_ccs_passes_hist.png')
                plt.clf()

                # Make CCS Passes Table
                passes_headers = ["Avg. CCS Passes", "Min. CCS Passes", "Max. CCS Passes"]
                passes_table = PrettyTable(passes_headers)
                passes_table.add_row([avg_num_ccs_passes, min_num_ccs_passes, max_num_ccs_passes])

                # Write table to output file
                fout.write(str(passes_table))
                fout.write(str('\n'))

                # Return to beginning directory

            os.chdir("..")
            os.chdir("..")
            os.chdir("..")
            
#again, i dont really like this, but look later for a better/more secure way. this could easily get messed up if the file convention is wrong. I could save a pwd command in a variable to get the run directory and use that.... I should do that.
            
            # Lastly, return the total number of post-ccs reads (count) to compare with polymerase reads
        return total_ccs_counts_all_cells

    print("\nCCS Sub-Report Complete.\n")

def check_demux_files(smrt_cells):
    print("Confirming Existence of DEMUX files...")
    for cell in smrt_cells:
        os.chdir(f'{cell}/')

        barcode_to_sample_csv = glob.glob('*.csv')
        #is the barcode to sample name csv in the smrt cell subdirectory?
        try:
            #glob returns a list, so just take the first element.
            with open(barcode_to_sample_csv[0]):
                print("Located Barcode Input File Successfully.")
        except:
            #if this file is not found, stop.
            raise FileNotFoundError("barcode to sample csv not found. Make sure this csv is in SMRT cell subdirectory")

        finally:
            try:
                is_summary_csv_here = os.path.exists("demux_no_peek/outputs/barcode_ccs_summary.csv")
                if is_summary_csv_here == True:
                    print("DEMUX Summary File Found.")
                    print("\nAll %s DEMUX Files Found.\n" % cell)
                    os.chdir("..")
                    
            except:
                raise FileNotFoundError("DEMUX Summary File Not Found.")

def get_demux_summary(smrt_cells, total_post_ccs_reads_list):
    #Reminder: Now we are back in the run_id directory with smrt cells as subdirectories in this folder.
    
    # Adding to the summary report file created above
    with open("ccs_demux_summary_report.txt", "a") as fout:
        # want to report stats for all cells
        fout.write(str('\n'))
        fout.write(str('\nDEMUX Sub-Report\n'))
        for cell in smrt_cells:
            fout.write(str(f'\nSMRT CELL: {cell}\n'))

            # cd into directory with existing CCS files
            os.chdir(f'{cell}/')

            # obtain/generate summary data and statistics
            barcode_to_sample_csv = glob.glob('*.csv')

            # again, glob returns a list, so take first element 
            with open(barcode_to_sample_csv[0]) as f:
                # skip the header line
                next(f)
                # get list of expected barcodes and write to output file
                expected_barcodes = [line.split(",")[0] for line in f]
                fout.write("\nList of Input Barcodes:\n")
                for barcode in expected_barcodes:
                    fout.write(str(f'{barcode}\n'))

            # get into DEMUX output directory
            os.chdir("demux_no_peek/outputs/")
            # open report file
            with open("barcode_ccs_summary.csv", 'r') as f:
                # skip header line
                next(f)
                # get list of demultiplexed barcodes and write to output file.
                demultiplexed_barcodes = [line.split(",")[1] for line in f]
                fout.write("\nList of Demultiplexed Barcodes:\n")
                # record the number of barcodes demultiplexed for table entries
                num_of_barcodes = 0
                for barcode in demultiplexed_barcodes:
                    num_of_barcodes +=1
                    fout.write(str(f'{barcode}\n'))

            # Add list of barcodes expected to find but did not find
                # initiate list
                barcodes_not_found = []
                # compare lists of expected and found barcodes
                for input_barcode in expected_barcodes:
                     if input_barcode not in demultiplexed_barcodes:
                         barcodes_not_found.append(input_barcode)
                # write to output file
                fout.write("\nList of Expected Barcodes Not Found:\n")
                for barcode in barcodes_not_found:
                    fout.write(str(f'{barcode}\n'))

                # notify if all expected barcodes were found
                if len(barcodes_not_found) == 0:
                    fout.write(str("**All Expected Barcodes Found**\n"))

                # Summary Stats

            # Polymerase Reads
            with open("barcode_ccs_summary.csv", 'r') as f:
                # Initiated variables    
                count = 0
                total_pr = 0
                l = []
                polymerase_reads = []

                # Get num from line in file
                for line in f:
                    l.append(line)
                for line in l:
                    polymerase_reads.append(line.split(",")[3])
                    
                # Get rid of header
                polymerase_reads.pop(0)

                # get total of polymerase reads
                for num in polymerase_reads:
                    num = int(num)
                    total_pr += int(num)
                    count += 1

                #fout.write("\nTotal Number of Polymerase Reads:\n")
                #fout.write(str(total_pr))

                # compute average
                polymerase_reads_avg = int(total_pr / count)


                # while I'm here, I'm also going to compute the average for the ccs reads
                ccs_reads_avg = []
                for i in total_post_ccs_reads_list:
                    reads_avg = int(i / count)
                    ccs_reads_avg.append(reads_avg)
                    
                #fout.write("\nAverage Number of Polymerase Reads:\n") #Confirm: Number or Length??
                #fout.write(str(round(polymerase_reads_avg)))
                # compute min and max for axes
                max_polymerase_reads = max(polymerase_reads)
                min_polymerase_reads = min(polymerase_reads)

                # Make bar chart
                fout.write("\nBar Chart of Polymerase Reads Saved to png file\n")
                plt.clf()
                polymerase_reads_f = [float(i) for i in polymerase_reads]
                plt.bar(np.arange(len(demultiplexed_barcodes)), polymerase_reads_f, align='center', alpha=0.5)
                plt.xticks(np.arange(len(demultiplexed_barcodes)), demultiplexed_barcodes, rotation=90)
                plt.xlabel("Barcodes")
                plt.ylabel("Polymerase Reads")
                plt.title(f'{cell}: Bar Chart of Polymerase Reads')
                plt.tight_layout()
                plt.savefig(f'{cell}_polymerase_reads_bar_chart.png')
                plt.clf()

                # Mean Barcode Quality
            with open("barcode_ccs_summary.csv", "r") as f:

                # Initiate Variables
                barcode_quality_scores = [] 
                total_score = 0
                count = 0
                bqs_lines = []

                # Get barcode quality scores from lines in file
                for line in f:
                    bqs_lines.append(line)
                for line in bqs_lines:
                    barcode_quality_scores.append(line.split(",")[4])

                # Remove header (First element)
                barcode_quality_scores.pop(0)

                # Compute summary stats
                min_barcode_score = min(barcode_quality_scores)
                max_barcode_score = max(barcode_quality_scores)

                for score in barcode_quality_scores:
                    total_score += int(score)
                    count += 1

                avg_barcode_quality_score = total_score / count

                fout.write("\nMin Barcode Quality Score:\n")
                fout.write(str(min_barcode_score))
                fout.write("\nMax Barcode Quality Score:\n")
                fout.write(str(max_barcode_score))
                fout.write("\nAverage Barcode Quality Score:\n")
                fout.write(str(round(avg_barcode_quality_score)))
                fout.write("\n")
                fout.write("\n")


                # Create a table structure for organizing output data
                with open("barcode_ccs_summary.csv", 'r') as f:
                    # initiate variables
                    all_lines = []
                    sample_names = []
                    demux_data = []
                    table_headers = ["Barcode", "Sample_Name", "Polymerase_Reads"] #,"Avg_Read_Length"]

                    # Get sample names from lines in file
                    for line in f:
                        all_lines.append(line)
                    for line in all_lines:
                        sample_names.append(line.split(",")[0])
                    
                    #remove header from list
                    sample_names.pop(0)

                    # for each barcode/sample, use list nesting to group associated data
                    # floats cannot be accessed by indexing, so make those values ints
                    for i in range(0, num_of_barcodes):
                        demux_data.append([demultiplexed_barcodes[i], sample_names[i], polymerase_reads[i]]) #, polymerase_reads_avg])
                        
                    # initiate table with headers
                    demux_table = PrettyTable(table_headers)

                    # populate table
                    for data in demux_data:
                        demux_table.add_row(data)
                        
                    # write to output file
                    fout.write(str(demux_table))
                    fout.write(str('\n'))
                    # Print this as well so you can get a quick understanding
                    #print(str(demux_table))

                    # Report the average difference in ccs reads to polymerase reads
                    #print(total_post_ccs_reads_list)
                    
                    if cell == "1_A01":
                        if ccs_reads_avg[0] < polymerase_reads_avg:
                            avg_ccs_polymerase_perc_increase = round( ( (polymerase_reads_avg - ccs_reads_avg[0]) / ccs_reads_avg[0] ) * 100 )
                            fout.write(str(f'\nOn average, the difference in CCS reads to Polymerase reads is {avg_ccs_polymerase_perc_increase}%.\n'))
                        elif ccs_reads_avg[0] > polymerase_reads_avg:
                            avg_ccs_polymerase_perc_decrease = round( ( (ccs_reads_avg[0] - polymerase_reads_avg) / ccs_reads_avg[0] ) * 100 )
                            fout.write(str(f'\nOn average, the difference in CCS reads to Polymerase reads is {avg_ccs_polymerase_perc_decrease}%.\n'))
                            
                    elif cell == "2_B01":
                        if ccs_reads_avg[1] < polymerase_reads_avg:
                            avg_ccs_polymerase_perc_increase = round( ( (polymerase_reads_avg - ccs_reads_avg[1]) / ccs_reads_avg[1] ) * 100 )
                            fout.write(str(f'\nOn average, the difference in CSS reads to Polymerase reads is {avg_ccs_polymerase_perc_increase}%.\n'))
                        elif ccs_reads_avg[1] > polymerase_reads_avg:
                            avg_ccs_polymerase_perc_decrease = round( ( (ccs_reads_avg[1] - polymerase_reads_avg) / ccs_reads_avg[1] ) * 100 )
                            fout.write(str(f'\nOn average, the difference in CCS reads to Polymerase reads is {avg_ccs_polymerase_perc_decrease}%.\n'))
                            
                    elif cell == "3_C01":
                        if ccs_reads_avg[2] < polymerase_reads_avg:
                            avg_ccs_polymerase_perc_increase = round( ( (polymerase_reads_avg - ccs_reads_avg[2]) / ccs_reads_avg[2] ) * 100 )
                            fout.write(str(f'\nOn average, the difference in CCS reads to Polymerase reads is {avg_ccs_polymerase_perc_increase}%.\n'))
                        elif ccs_reads_avg[2] > polymerase_reads_avg:
                            avg_ccs_polymerase_perc_decrease = round( ( (ccs_reads_avg[2] - polymerase_reads_avg) / ccs_reads_avg[2] ) * 100 )
                            fout.write(str(f'\On average, the difference in CCS reads to Polymerase reads is {avg_ccs_polymerase_perc_decrease}%.\n'))
                            
                    elif cell == "4_D01":
                        if ccs_reads_avg[3] < polymerase_reads_avg:
                            avg_ccs_polymerase_perc_increase = round( ( (polymerase_reads_avg - ccs_reads_avg[3]) / ccs_reads_avg[3] ) * 100 )
                            fout.write(str(f'\nOn average, the difference in CCS reads to Polymerase reads is {avg_ccs_polymerase_perc_increase}%.\n'))
                        elif ccs_reads_avg[3] > polymerase_reads_avg:
                            avg_ccs_polymerase_perc_decrease = round( ( (ccs_reads_avg[3] - polymerase_reads_avg) / ccs_reads_avg[3] ) * 100 )
                            fout.write(str(f'\nOn average, the difference in CCS reads to Polymerase reads is {avg_ccs_polymerase_perc_decrease}%.\n')) 
                    
                # Return to original run directory for accessing other cells
                os.chdir("..")
                os.chdir("..")
                os.chdir("..") 


    print("DEMUX Sub-Report Complete.")


def consolidate_files(smrt_cells):
    # Move all files to the output directory
    for cell in smrt_cells:
        # Create target variables
        run_dir = subprocess.run(["pwd"], stdout=subprocess.PIPE).stdout.decode('utf-8')
        # New line is attached to the end of run_dir to make my life harder, so remove that before appending output dir
        run_dir_corrected = run_dir.strip()
        target_dir = run_dir_corrected + "/ccs_demux_report_output"
        
        # Move CCS files
        os.chdir(f'{cell}/ccs/outputs/')
        os.rename(f'{cell}_avg_read_length_hist.png', f'{target_dir}/{cell}_avg_read_length_hist.png')
        #os.rename("ccs_passes_histogram.txt", "%s/ccs_passes_histogram.txt" % target_dir)
        os.rename(f'{cell}_ccs_passes_hist.png', f'{target_dir}/{cell}_ccs_passes_hist.png')
        
        # Move Demux files
        os.chdir("..")
        os.chdir("..")
        os.chdir("demux_no_peek/outputs/")

        os.rename(f'{cell}_polymerase_reads_bar_chart.png', f'{target_dir}/{cell}_polymerase_reads_bar_chart.png')

        # Go back to run directory
        os.chdir("..")
        os.chdir("..")
        os.chdir("..")

        # Move the json files to output dir
        os.rename(f'{cell}_data.json', f'{target_dir}/{cell}_data.json')

    # Finally, move the summary text file to the output dir
    os.rename("ccs_demux_summary_report.txt", "%s/ccs_demux_summary_report.txt" % target_dir)

    
def get_productivity_values(smrt_cells):
    print("\nObtaining Productivity Values...\n")
    for cell in smrt_cells:
        # use pwd to obtain the run id, which should be the directory above our current directory
        run_dir = subprocess.run(["pwd"], stdout=subprocess.PIPE) .stdout.decode('utf-8')
        # because run_ids have the same length and there is a naming convention, I can use str splicing to take the last X characters to make up the run ID.
        run_id = run_dir[-23:]
        run_id = run_id.rstrip()
        # now I need to get tge path to the sts.xml file associated with the run post sequencing
        pacbio_sequel_data_dir = "/data/pacbio/sequel/userdata/data_root/"
        # full path to file str cat
        data_dir = pacbio_sequel_data_dir + run_id + "/" + cell
        # let's create a softlink to the file so we can easily get back to the current
        os.system(f'ln -s {data_dir}/*.sts.xml')
        # assign file to variable - glob outputs list
        target_file = glob.glob('*.sts.xml')
        
        # XML Parsing is proving to be complicated... I think the sts.xml file must be compressed or not in standard xml format, so parsing/iterating through the file returns nothing. Going to try to convert to json file and pull info that way.

        # Convert XML to json file
        json_file = convert_xml_to_json(target_file[0])
        # Extract Data
        with open(json_file) as j_file:
            data = json.load(j_file)
            # Productivity Values
            p_zero = data['PipeStats']['ProdDist']['ns:BinCounts']['ns:BinCount'][0]
            p_one = data['PipeStats']['ProdDist']['ns:BinCounts']['ns:BinCount'][1]
            p_two = data['PipeStats']['ProdDist']['ns:BinCounts']['ns:BinCount'][2]
        
            # Make a table in output summary report file
            p_zero_header = u'P\u2080'
            p_one_header = u'P\u2081'
            p_two_header = u'P\u2082'

            p_table_headers = [p_zero_header, p_one_header, p_two_header]

            p_table = PrettyTable(p_table_headers)
            p_table.add_row([p_zero, p_one, p_two])
            
            with open("ccs_demux_summary_report.txt", "a") as fout:
                # Write to output file
                fout.write(str("\n"))
                fout.write(str(f'\nProductivity Values {cell}:\n'))
                fout.write(str(p_table))
                
                # remove the soft link of the .sts.xml file so we can repeat with next cell(s)
                os.system(f'rm *.sts.xml')
        # rename file for clarity
        os.rename('data.json', f'{cell}_data.json')
                    

def convert_xml_to_json(xml_file):
    # read xml file and convert to dict object
    with open(xml_file, 'r') as xml_file:
        data_dict = xmltodict.parse(xml_file.read())
        xml_file.close()
        
        # convert to json object
        json_data = json.dumps(data_dict, indent=2)

        # write json data to output json file
        with open("data.json", "w") as json_file:
            json_file.write(json_data)
            json_file.close()

        return "data.json"
    

def main():
    #First, Make Output Directory
    print("Creating Output Directory...\n")
    os.system("mkdir ccs_demux_report_output")
    #Second, Make Sure CCS Files Exist
    print("Generating CCS Report Summary...\n")
    smrt_cells = check_smrt_cells_and_ccs_files()
    #Third, Make CSS Sub-Report
    total_post_ccs_reads = get_ccs_summary(smrt_cells)
    #Fourth, Make Sure DEMUX Files Exist
    print("Generating DEMUX Report Summary...\n")
    check_demux_files(smrt_cells)
    #Fifth, Make DEMUX Sub-Report
    get_demux_summary(smrt_cells, total_post_ccs_reads)
    #Sixth, obtain productivity values 
    get_productivity_values(smrt_cells)
    #Seventh, consolidate files to output dir
    consolidate_files(smrt_cells)
    # Fin
    print("\nCCS_DEMUX_Report Successfully Completed.\n")
    print("Please See ccs_demux_report_output/ For All Output Files.\n")

if __name__=="__main__":
	main()

