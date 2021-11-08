# microbiome_analysis_pipeline_report
Generates a report following successful completing Circular Consensus Sequencing, Demultiplexing (Microbiome Analysis Pipeline), and preparation for the MCSMRT Microbiome Classifier Pipeline using Pacbio sequencing data.

Please note the standardardized file-naming convention that must be used in order to run this program.

Organization is as follows: 
project_dir/SMRTcell_subdir/CCS_or_DEMUX_subdir/outputs

example: 
m54257_12345_20211108/1_A01/ccs/outputs/
another example: 
m54257_67890_202111/2_B01/demux_no_peek/outputs/

It may be possible to edit file names and directory structure as per your own needs, but this is the convention we use.

Required Input Files:
CCS Sub-report:
1) ccs_passes.tsv --> 
a tsv file with recorded number of passess associated with each read (required for MCSMRT pipeline, script can be found on MCSMRT repository)
2) ccs_report.zip OR ccs_statistics.csv --> 
a csv file output from pbcromwell pb_ccs 

DEMUX Sub-report:
1) barcode_to_sample.csv --> 
a csv file with a translation of barcode ID and bio sample ID. (barcode,sample)
2) barcode_ccs_summary.csv --> 
a csv file output from pbcromwell pb_demux_ccs

Output:
1) ccs_demux_report_output/ (directory)
a) ccs_demux_summary_report.txt
b) SMRTcell_avg_read_length_hist.png
c) SMRTcell_ccs_passes_hist.png
d) SMRTcell_polymerase_reads_hist.png
where SMRTcell represents a prefix of the cell-associated data

ccs_demux_summary_report.txt contains useful information/reports and plots of statistics including of the number of reads, avg read length, avg number of ccs passes, polymerase reads, lists of input barcodes, demultiplexed barcodes, barcodes not found during demultiplexing, and relational tables of bio sample id, associated barcode, polymerase reads, and avg polymerase reads associated with each SMRTcell.
