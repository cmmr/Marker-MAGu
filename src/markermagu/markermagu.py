#!/usr/bin/env python

import argparse
import sys, os
import subprocess

__version__='0.3.0'

## function to allow boolean (true/false) input arguments
def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
    
pathname = os.path.dirname(__file__)

markermagu_script_path = os.path.abspath(pathname)      
print(markermagu_script_path) 

## this is the main entry point function
## users run `markermagu` on the command line with the arguments for argparse
## these arguments are formatted, dependencies are checked
## then arguments are fed to Marker-MAGu_mapper.sh
def markermagu():

    parser = argparse.ArgumentParser(description='Marker-MAGu is a read mapping pipeline which uses marker genes to detect and measure bacteria, phages, archaea, and microeukaryotes. Version ' + str(__version__))

    required_args = parser.add_argument_group(' REQUIRED ARGUMENTS for Marker-MAGu ')

    required_args.add_argument("-r", "--reads", nargs="+",
                            dest="READS", required=True, 
                            help='read file(s) in .fastq format. You can specify more than one separated by a space')
    required_args.add_argument("-s", "--sample", 
                            dest="SAMPLE", type=str, required=True, 
                            help='Sample name. No space characters, please.')
    required_args.add_argument("-t", "--cpu", 
                            dest="CPU", type=int, required=True, 
                            help='Example: 32 -- Number of CPUs available for Marker-MAGu. ')
    required_args.add_argument("-o", "--output_dir", 
                            dest="OUTPUT_DIR", type=str, required=True, 
                            help='Output directory name. Will be created if it does not exist. Can be shared with other samples. No space characters, please. ')

    

    optional_args = parser.add_argument_group(' OPTIONAL ARGUMENTS for Marker-MAGu.')

    optional_args.add_argument('--version', action='version', version=str(__version__))
    optional_args.add_argument('-q', "--qual", dest="QUAL", type=str2bool, default='False',
                            help='True or False. Remove low-quality reads with fastp?')
    optional_args.add_argument('-f', "--filter_seqs", dest="FILTER_SEQS", type=str2bool, default='False',
                            help='True or False. Remove reads aligning to sequences at filter_seqs/filter_seqs.fna ?')
    optional_args.add_argument("--filter_dir", 
                            dest="FILTER_DIR", type=str, default='default',
                            help='path to directory of sequences to filter. If not set, Marker-MAGu looks for environmental variable MARKERMAGU_FILTER. Then, if this variable is unset, it this is unset, DB path is assumed to be ' + markermagu_script_path.replace("src/markermagu", "filter_seqs"))
    optional_args.add_argument("--temp", 
                            dest="TEMP_DIR", type=str, default='default',
                            help='path of temporary directory. Default is {OUTPUT_DIR}/{SAMPLE}_temp/')
    optional_args.add_argument("--keep", 
                            dest="KEEP", type=str2bool, default='False',
                            help='True of False. Keep the intermediate files, located in the temporary directory? These can add up, so it is not recommended if space is a concern.')
    optional_args.add_argument("--db", 
                            dest="DB", type=str, default='default',
                            help='DB path. If not set, Marker-MAGu looks for environmental variable MARKERMAGU_DB. Then, if this variable is unset, it this is unset, DB path is assumed to be ' + markermagu_script_path.replace("src", "DBs/v1.0"))

    args = parser.parse_args()

    ## joins read files together when provided with spaces in between
    READS = ' '.join(map(str,args.READS))

    ## DB path check/change
    if args.DB == "default" and os.getenv('MARKERMAGU_DB') != None:
        args.DB = os.getenv('MARKERMAGU_DB')
    elif args.DB == "default":
        args.DB = markermagu_script_path.replace("src", "DBs/v1.0")

    ## filter seq path check/change
    if args.FILTER_DIR == "default" and os.getenv('MARKERMAGU_FILTER') != None:
        args.FILTER_DIR = os.getenv('MARKERMAGU_FILTER')
    elif args.FILTER_DIR == "default":
        args.FILTER_DIR = markermagu_script_path.replace("src/markermargu", "filter_seqs")


    ## check if R script with library check returns good exit code
    completedProc = subprocess.run(['Rscript', str(markermagu_script_path) + '/check_R_libraries1.R'])

    print(completedProc.returncode)
    if completedProc.returncode != 0 :
        print ("some required R packages are not found. Required:")
        print ("dplyr, data.table, stringr")
        print ("Did you activate the conda environment?")
        print ("see yml. Exiting")
        quit()


    def is_tool(name):
        """Check whether `name` is on PATH."""
        from distutils.spawn import find_executable
        return find_executable(name) is not None

    if not is_tool("coverm") :
        print ("coverm is not found. Exiting.")
        quit()
    if not is_tool("minimap2") :
        print ("minimap2 is not found. Exiting.")
        quit()
    if not is_tool("samtools") :
        print ("samtools is not found. Exiting.")
        quit()
    if not is_tool("seqkit") :
        print ("seqkit is not found. Exiting.")
        quit()
    if not is_tool("fastp") :
        print ("fastp is not found. Exiting.")
        quit()
        
    ## this actually calls the main mapper bash script 
    ## and provides all the arguments taken in this file
    subprocess.call(['bash', str(markermagu_script_path) + '/Marker-MAGu_mapper.sh', 
                str(READS), str(args.SAMPLE), str(args.CPU), str(args.OUTPUT_DIR), 
                str(args.QUAL), str(args.FILTER_SEQS), str(args.FILTER_DIR), str(args.TEMP_DIR), 
                str(args.KEEP), str(args.DB), str(__version__), str(markermagu_script_path)])

