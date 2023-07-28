from tqdm.contrib.telegram import tqdm
from pathos.pools import ProcessPool
from gdbr.version import get_version

import subprocess
import argparse
import datetime
import p_tqdm
import shutil
import vcfpy
import json
import glob
import os

def get_telegram_data(file_loc):
    if file_loc is None:
        return None
    
    if os.path.exists(file_loc):
        with open(file_loc, 'r') as f:
            data = json.load(f)
            return str(data['token']), str(data['chat_id'])
    else:
        return None


def get_proper_thread(suggest_num_cpus, num_cpus, num_works=-1):
    suggest_num_cpus = max(1, suggest_num_cpus)

    if num_cpus <= suggest_num_cpus:
        hard_num_cpus = num_cpus
        loop_num_cpus = 1
    elif num_cpus // num_works >= suggest_num_cpus:
        hard_num_cpus = num_cpus // num_works
        loop_num_cpus = num_works
    else:
        cpu_usage_list = [num_cpus % i for i in range(suggest_num_cpus, suggest_num_cpus + 3)]
        hard_num_cpus = suggest_num_cpus + cpu_usage_list.index(min(cpu_usage_list))
        loop_num_cpus = num_cpus // hard_num_cpus

    return hard_num_cpus, loop_num_cpus


def p_map(f, it, num_cpus=1, pbar=True, telegram_token_loc='telegram.json', desc=''):
    if pbar:
        telegram_data = get_telegram_data(telegram_token_loc)
        if telegram_data is None:
            return p_tqdm.p_map(f, it, num_cpus=num_cpus, desc=desc)
        else:
            # only telegram taskbar; silent stdout
            return p_tqdm.p_map(f, it, num_cpus=num_cpus, tqdm=tqdm, token=telegram_data[0], chat_id=telegram_data[1], desc=desc, file=open(os.devnull, 'w'))
    else:
        pool = ProcessPool(num_cpus)
        return pool.map(f, it)


def logprint(s):
    print(f'[{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] ' + s)


def gdbr_parser():
    parser = argparse.ArgumentParser(prog='gdbr', description='GDBr is tool that identify Double-strand Break Repair (DSBR) using genome and variant calling.')

    parser.add_argument('-v', '--version',
                        action='version',
                        version=f'%(prog)s v{get_version()}')
    
    subparsers = parser.add_subparsers(help='modes', dest='command')
    subparsers.required = True

    parser_anl = subparsers.add_parser('analysis', help='find DSBR by corrected variant calling', add_help=False)
    parser_pre = subparsers.add_parser('preprocess', help='preprocess the genome by scaffolding and have same chromosome name with reference', add_help=False)
    parser_cor = subparsers.add_parser('correct', help='correct the variant calling using BLAST', add_help=False)
    parser_ant = subparsers.add_parser('annotate', help='find DSBR by corrected variant calling', add_help=False)
    
    # analysis_main
    # analysis required arguments
    re_parser_anl = parser_anl.add_argument_group('required argument')
    re_parser_anl.add_argument('-r', '--reference',
                               type=os.path.abspath,
                               help='reference sequence location',
                               required=True)

    re_parser_anl.add_argument('-q', '--query',
                               type=os.path.abspath,
                               nargs='+',
                               help='raw query sequence locations',
                               required=True)
    
    re_parser_anl.add_argument('-s', '--species',
                               type=str,
                               help="Species of data for RepeatMasker, All unique clade names occurring in this database (http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html) can be used.",
                               required=True)
    
    # analysis functional optional argument
    fo_parser_anl = parser_anl.add_argument_group('functional optional arguments')
    fo_parser_anl.add_argument('-h', '--help',
                               action='help',
                               help='show this help message and exit')
    
    fo_parser_anl.add_argument('-o', '--output',
                               type=os.path.abspath,
                               default='gdbr_output',
                               help='preprocessed query save location')
        
    fo_parser_anl.add_argument('-t', '--threads',
                               type=int,
                               default=1,
                               help='number of threads')
    
    fo_parser_anl.add_argument('--workdir',
                               type=os.path.abspath,
                               default='gdbr_data',
                               help='program work directory')
    
    fo_parser_anl.add_argument('--silent',
                               action='store_true',
                               help='turn off progress bar')
    
    fo_parser_anl.add_argument('--telegram_data_loc',
                               type=os.path.abspath,
                               default=None,
                               help='for telegram progress bar suppport, this file must have "token" and "chat_id" with .json format')

    # analysis preprocess optional arguments
    po_parser_anl = parser_anl.add_argument_group('preprocess optional arguments')
    po_parser_anl.add_argument('--pre_min_sv_size',
                               type=int,
                               default=50,
                               help='minimum variant size at preprocess')
    
    po_parser_anl.add_argument('--low_memory',
                               action='store_true',
                               help='turn off query multiprocessing and reduce memory usage')
    
    po_parser_anl.add_argument('--trust_workdir',
                               action='store_true',
                               help='When you run the preprocess again, speed up the preprocess with data left in the work directory.')
    
    # analysis correct optional arguments
    co_parser_anl = parser_anl.add_argument_group('correct optional arguments')
    co_parser_anl.add_argument('--sv_find_len',
                               type=int,
                               default=2000,
                               help='sequence length to correct variant calling')
    
    co_parser_anl.add_argument('--repeat_find_len',
                               type=int,
                               default=50,
                               help='sequence length to identity repeat nearby variant')
    
    co_parser_anl.add_argument('--cor_min_sv_size',
                               type=int,
                               default=50,
                               help='minimum variant size at correct')
    

    # analysis annotate optional arguments
    ao_parser_anl = parser_anl.add_argument_group('annotate optional arguments')
    ao_parser_anl.add_argument('--hom_find_len',
                               type=int,
                               default=2000,
                               help='sequence length to find DSBR')

    ao_parser_anl.add_argument('--temp_indel_find_len',
                               type=int,
                               default=100,
                               help='sequence length to find template insertion')
    
    ao_parser_anl.add_argument('--near_gap_find_len',
                               type=int,
                               default=5,
                               help='search length to remove gaps at template insertion')
    
    ao_parser_anl.add_argument('--temp_gap_baseline',
                               type=int,
                               default=3,
                               help='maximum number of gaps allowed at template insertion')
    
    ao_parser_anl.add_argument('--near_sv_kb_baseline',
                               type=float,
                               default=100.0,
                               help='minimum length (kb) to be recognized as different variant')
    
    ao_parser_anl.add_argument('--diff_locus_hom_baseline',
                               type=int,
                               default=3,
                               help='minimum homology length to find different locus DSBR')
    
    # preprocess_main
    # preprocess required arguments
    re_parser_pre = parser_pre.add_argument_group('required argument')
    re_parser_pre.add_argument('-r', '--reference',
                               type=os.path.abspath,
                               help='reference sequence location',
                               required=True)

    re_parser_pre.add_argument('-q', '--query',
                               type=os.path.abspath,
                               nargs='+',
                               help='raw query sequence locations',
                               required=True)

    # preprocess optional arguments
    op_parser_pre = parser_pre.add_argument_group('optional arguments')
    op_parser_pre.add_argument('-h', '--help',
                               action='help',
                               help='show this help message and exit')
    
    op_parser_pre.add_argument('-o', '--preprocess_save',
                               type=os.path.abspath,
                               default='gdbr_prepro',
                               help='preprocessed query save location')
        
    op_parser_pre.add_argument('-t', '--threads',
                               type=int,
                               default=1,
                               help='number of threads')
    
    op_parser_pre.add_argument('--min_sv_size',
                               type=int,
                               default=50,
                               help='minimum variant size')
    
    op_parser_pre.add_argument('--low_memory',
                               action='store_true',
                               help='turn off query multiprocessing and reduce memory usage')
    
    op_parser_pre.add_argument('--trust_workdir',
                               action='store_true',
                               help='When you run the preprocess again, speed up the preprocess with data left in the work directory.')
    
    op_parser_pre.add_argument('--workdir',
                               type=os.path.abspath,
                               default='gdbr_data',
                               help='program work directory')
    
    op_parser_pre.add_argument('--silent',
                               action='store_true',
                               help='turn off progress bar')
    
    op_parser_pre.add_argument('--telegram_data_loc',
                               type=os.path.abspath,
                               default=None,
                               help='for telegram progress bar suppport, this file must have "token" and "chat_id" with .json format')

    # correct_main
    # correct required arguments
    re_parser_cor = parser_cor.add_argument_group('required arguments')
    re_parser_cor.add_argument('-r', '--reference',
                               type=os.path.abspath,
                               help='reference sequence location',
                               required=True)
    
    re_parser_cor.add_argument('-q', '--query',
                               type=os.path.abspath,
                               nargs='+',
                               help='preprocessed query sequence locations (same choromosome name with reference)',
                               required=True)
    
    re_parser_cor.add_argument('-v', '--vcf',
                               type=os.path.abspath,
                               nargs='+',
                               help='variant calling file location (only SVIM-asm supported)',
                               required=True)

    re_parser_cor.add_argument('-s', '--species',
                               type=str,
                               help="Species of data for RepeatMasker, All unique clade names occurring in this database (http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html) can be used.",
                               required=True)

    op_parser_cor = parser_cor.add_argument_group('optional arguments')
    op_parser_cor.add_argument('-h', '--help',
                               action='help',
                               help='show this help message and exit')
    
    op_parser_cor.add_argument('-o', '--sv_save',
                               type=os.path.abspath,
                               default='gdbr_sv',
                               help='corrected variant CSV save location')
    
    op_parser_cor.add_argument('-t', '--threads',
                               type=int,
                               default=1,
                               help='number of threads')
    
    op_parser_cor.add_argument('--sv_find_len',
                               type=int,
                               default=2000,
                               help='sequence length to correct variant calling')
    
    op_parser_cor.add_argument('--repeat_find_len',
                               type=int,
                               default=50,
                               help='sequence length to identity repeat nearby variant')

    op_parser_cor.add_argument('--min_sv_size',
                               type=int,
                               default=50,
                               help='minimum variant size (you must use this option when not using variant calling file from GDBr )')
    
    op_parser_cor.add_argument('--workdir',
                               type=os.path.abspath,
                               default='gdbr_data',
                               help='program work directory')
    
    op_parser_cor.add_argument('--silent',
                               action='store_true',
                               help='turn off progress bar')
    
    op_parser_cor.add_argument('--telegram_data_loc',
                               type=os.path.abspath,
                               default=None,
                               help='for telegram progress bar suppport, this file must have "token" and "chat_id" with .json format')

    # annotate main
    # annotate required arguments
    re_parser_ant = parser_ant.add_argument_group('required arguments')
    re_parser_ant.add_argument('-r', '--reference',
                               type=os.path.abspath,
                               help='reference sequence location',
                               required=True)

    re_parser_ant.add_argument('-q', '--query',
                               type=os.path.abspath,
                               nargs='+',
                               help='preprocessed query sequence locations (same choromosome name with reference)',
                               required=True)
    
    re_parser_ant.add_argument('-v', '--sv_csv',
                               type=os.path.abspath,
                               nargs='+',
                               help='corrected variant CSV file',
                               required=True)

    op_parser_ant = parser_ant.add_argument_group('optional arguments')
    op_parser_ant.add_argument('-h', '--help',
                               action='help',
                               help='show this help message and exit')
    
    op_parser_ant.add_argument('-o', '--result_save',
                               type=os.path.abspath,
                               default='gdbr_result',
                               help='DSBR annotate save location')
    
    op_parser_ant.add_argument('-t', '--threads',
                               type=int,
                               default=1,
                               help='number of threads')
    
    op_parser_ant.add_argument('--hom_find_len',
                               type=int,
                               default=2000,
                               help='sequence length to find DSBR')

    op_parser_ant.add_argument('--temp_indel_find_len',
                               type=int,
                               default=100,
                               help='sequence length to find template insertion')
    
    op_parser_ant.add_argument('--near_gap_find_len',
                               type=int,
                               default=5,
                               help='search length to remove gaps at template insertion')
    
    op_parser_ant.add_argument('--temp_gap_baseline',
                               type=int,
                               default=3,
                               help='maximum number of gaps allowed at template insertion')
    
    op_parser_ant.add_argument('--near_sv_kb_baseline',
                               type=float,
                               default=100.0,
                               help='minimum length (kb) to be recognized as different variant')
    
    op_parser_ant.add_argument('--diff_locus_hom_baseline',
                               type=int,
                               default=3,
                               help='minimum homology length to find different locus DSBR')
    
    op_parser_ant.add_argument('--workdir',
                               type=os.path.abspath,
                               default='gdbr_data',
                               help='program work directory')
    
    op_parser_ant.add_argument('--silent',
                               action='store_true',
                               help='turn off progress bar')
    
    op_parser_ant.add_argument('--telegram_data_loc',
                               type=os.path.abspath,
                               default=None,
                               help='for telegram progress bar suppport, this file must have "token" and "chat_id" with .json format.')

    return parser

def check_proper_species(spec, workdir):
    os.makedirs(workdir, exist_ok=True)
    test_seq_loc = os.path.join(workdir, 'test_rpm.fa')
    with open(test_seq_loc, 'w') as f:
        f.write('>test\nATGC')
    
    rpm_result = subprocess.run(['RepeatMasker', '-spec', spec, '-nopost', test_seq_loc], capture_output=True, text=True, cwd=workdir)

    # remove RepeatMasker output
    for file in glob.glob(test_seq_loc + '*'):
        if os.path.isfile(file):
            os.remove(file)

    if rpm_result.stderr != '':
        # remove RepeatMasker Temp folder when failed
        for folder in glob.glob(os.path.join(workdir, 'RM_*')):
            if os.path.isdir(folder):
                os.rmdir(folder)
        raise Exception(rpm_result.stderr)


def check_file_exist(file_loc_data, file_type_list):
    for file_loc_list, file_type in zip(file_loc_data, file_type_list):
        for file_loc in file_loc_list:
            if not os.path.isfile(file_loc):
                raise Exception(f'{file_type} file : {file_loc} missing')


def check_unique_basename(file_loc_list):
    file_basename_list = list(map(os.path.basename, file_loc_list))

    if len(file_basename_list) != len(set(file_basename_list)):
        raise Exception('Query basename must be diffrent')


def check_variant_caller(vcf_loc_list):
    for vcf_loc in vcf_loc_list:
        find = False
        with vcfpy.Reader.from_path(vcf_loc) as record:
            for header in filter(lambda t: type(t) == vcfpy.HeaderLine, record.header.lines):
                # svim-asm caller
                if header.key == 'source' and 'SVIM-asm' in header._value:
                    find = True
                    break

            if not find:
                raise Exception(f"Can't recongnize VCF file source at header : {record.path}")
        

def get_min_sv_size(vcf_loc):
    with vcfpy.Reader.from_path(vcf_loc) as record:
        for header in filter(lambda t: type(t) == vcfpy.HeaderLine, record.header.lines):
            if header.key == 'gdbr_min_sv_size':
                return int(header._value)

    return None


def clean_workdir(workdir):
    for item in os.listdir(workdir):
        if os.path.isdir(os.path.join(workdir, item)) and (item == 'db' or item.isdigit()):
            shutil.rmtree(os.path.join(workdir, item))

    if len(os.listdir(workdir)) == 0:
        os.rmdir(workdir)

def remove_gdbr_postfix(basename):
    if '.GDBr.' in basename:
        basename_list = basename.split('.')
        gdbr_index = basename_list.index('GDBr')
        return '.'.join(basename_list[:gdbr_index])
    else:
        return basename