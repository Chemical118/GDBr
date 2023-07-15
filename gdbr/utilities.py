from concurrent.futures import ProcessPoolExecutor
from tqdm.contrib.concurrent import process_map
from tqdm.contrib.telegram import tqdm
from gdbr.version import get_version

import argparse
import datetime
import p_tqdm
import json
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
        cpu_usage_list = [num_cpus % i for i in range(suggest_num_cpus - 1, suggest_num_cpus + 2)]
        hard_num_cpus = suggest_num_cpus - 1 + cpu_usage_list.index(min(cpu_usage_list))
        loop_num_cpus = num_cpus // hard_num_cpus

    return hard_num_cpus, loop_num_cpus


def p_map(f, it, num_cpus=1, pbar=True, telegram_token_loc='telegram.json', desc=''):
    if pbar:
        telegram_data = get_telegram_data(telegram_token_loc)
        if telegram_data is None:
            return process_map(f, it, max_workers=num_cpus, desc=desc, chunksize=1)
        else:
            # only telegram taskbar; silent stdout
            return process_map(f, it, max_workers=num_cpus, tqdm_class=tqdm, token=telegram_data[0], chat_id=telegram_data[1], desc=desc, file=open(os.devnull, 'w'), chunksize=1)
    else:
        executor = ProcessPoolExecutor(max_workers=num_cpus)
        return executor.map(f, it, chunksize=1)


def logprint(s):
    print(f'[{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] ' + s)


def gdbr_parser():
    parser = argparse.ArgumentParser(prog='gdbr', description='GDBr is tool that identify Double-strand Break Repair (DSBR) using genome and variant calling.')

    parser.add_argument('-v', '--version',
                        action='version',
                        version=f'%(prog)s v{get_version()}')
    
    subparsers = parser.add_subparsers(help='modes', dest='command')
    subparsers.required = True

    parser_pre = subparsers.add_parser('preprocess', help='preprocess the genome by scaffolding and have same chromosome name with reference', add_help=False)
    parser_cor = subparsers.add_parser('correct', help='correct the variant calling using BLAST', add_help=False)
    parser_anl = subparsers.add_parser('analysis', help='find DSBR by corrected variant calling', add_help=False)
    
    # preprocess_main
    # preprocess required arguments
    re_parser_pre = parser_pre.add_argument_group('required argument')
    re_parser_pre.add_argument('-r', '--reference',
                               type=os.path.abspath,
                               help='reference sequence location')

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
    
    op_parser_pre.add_argument('--low_memory',
                               action='store_true',
                               help='turn off query multiprocessing and reduce memory usage')
    
    op_parser_pre.add_argument('--workdir',
                               type=os.path.abspath,
                               default='gdbr_data',
                               help='program work directory')
    
    op_parser_pre.add_argument('--min_sv_size',
                               type=int,
                               default=50,
                               help='minimum variant size')
    
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
                               help='reference sequence location')
    
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

    op_parser_cor = parser_cor.add_argument_group('optional arguments')
    op_parser_cor.add_argument('-h', '--help',
                               action='help',
                               help='show this help message and exit')
    
    op_parser_cor.add_argument('-o', '--sv_save',
                               type=os.path.abspath,
                               default='gdbr_sv',
                               help='corrected variant CSV save location (default : same folder with variant calling file)')
    
    op_parser_cor.add_argument('-t', '--threads',
                               type=int,
                               default=1,
                               help='number of threads')
    
    op_parser_cor.add_argument('--sv_find_len',
                               type=int,
                               default=2000,
                               help='sequence length to correct variant calling')

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

    # analysis main
    # analysis required arguments
    re_parser_anl = parser_anl.add_argument_group('required arguments')
    re_parser_anl.add_argument('-r', '--reference',
                               type=os.path.abspath,
                               help='reference sequence location')

    re_parser_anl.add_argument('-q', '--query',
                               type=os.path.abspath,
                               nargs='+',
                               help='preprocessed query sequence locations (same choromosome name with reference)',
                               required=True)
    
    re_parser_anl.add_argument('-v', '--sv_csv',
                               type=os.path.abspath,
                               nargs='+',
                               help='corrected variant CSV file',
                               required=True)

    op_parser_anl = parser_anl.add_argument_group('optional arguments')
    op_parser_anl.add_argument('-h', '--help',
                               action='help',
                               help='show this help message and exit')
    
    op_parser_anl.add_argument('-o', '--dsbr_save',
                               type=os.path.abspath,
                               default='gdbr_dsbr',
                               help='DSBR analysis save location')
    
    op_parser_anl.add_argument('-t', '--threads',
                               type=int,
                               default=1,
                               help='number of threads')
    
    op_parser_anl.add_argument('--hom_find_len',
                               type=int,
                               default=2000,
                               help='sequence length to find DSBR')

    op_parser_anl.add_argument('--temp_indel_find_len',
                               type=int,
                               default=100,
                               help='sequence length to find template insertion')
    
    op_parser_anl.add_argument('--near_gap_find_len',
                               type=int,
                               default=5,
                               help='search length to remove gaps at template insertion')
    
    op_parser_anl.add_argument('--temp_gap_baseline',
                               type=int,
                               default=2,
                               help='maximum number of gaps allowed at template insertion')
    
    op_parser_anl.add_argument('--near_sv_kb_baseline',
                               type=float,
                               default=100.0,
                               help='minimum length (kb) to be recognized as different variant')
    
    op_parser_anl.add_argument('--diff_locus_hom_baseline',
                               type=int,
                               default=3,
                               help='minimum homology length to find different locus DSBR')
    
    op_parser_anl.add_argument('--workdir',
                               type=os.path.abspath,
                               default='gdbr_data',
                               help='program work directory')
    
    op_parser_anl.add_argument('--silent',
                               action='store_true',
                               help='turn off progress bar')
    
    op_parser_anl.add_argument('--telegram_data_loc',
                               type=os.path.abspath,
                               default=None,
                               help='for telegram progress bar suppport, this file must have "token" and "chat_id" with .json format')

    return parser