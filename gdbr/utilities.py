from tqdm.contrib.telegram import tqdm
from pathos.pools import ProcessPool
from gdbr.version import get_version
from pyfaidx import Fasta

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

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
    
    fo_parser_anl.add_argument('--min_sv_size',
                               type=int,
                               help='minimum variant size for preprocess and correct step')
    
    fo_parser_anl.add_argument('--benchmark',
                               type=float,
                               nargs='?',
                               const=None,
                               metavar='BENCHMARK_INTERVAL',
                               default=argparse.SUPPRESS,
                               help='turn on benchmark mode, if you specify an argument you can adjust the interval of memory benchmark')

    fo_parser_anl.add_argument('--workdir',
                               type=os.path.abspath,
                               default='gdbr_data',
                               help='program work directory')
    
    fo_parser_anl.add_argument('--overwrite_output',
                               action='store_true',
                               help='overwrite the output directory')

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
                               help='minimum variant size at preprocess step')
    
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
                               help='minimum variant size at correct step')
    

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
    
    op_parser_pre.add_argument('--benchmark',
                               type=float,
                               nargs='?',
                               const=None,
                               metavar='BENCHMARK_INTERVAL',
                               default=argparse.SUPPRESS,
                               help='turn on benchmark mode, if you specify an argument you can adjust the interval of memory benchmark')
    
    op_parser_pre.add_argument('--workdir',
                               type=os.path.abspath,
                               default='gdbr_data',
                               help='program work directory')
    
    op_parser_pre.add_argument('--overwrite_output',
                               action='store_true',
                               help='overwrite the output directory')

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
                               help='minimum variant size (you must use this option when not using variant calling file from GDBr )')
    
    op_parser_cor.add_argument('--benchmark',
                               type=float,
                               nargs='?',
                               const=None,
                               metavar='BENCHMARK_INTERVAL',
                               default=argparse.SUPPRESS,
                               help='turn on benchmark mode, if you specify an argument you can adjust the interval of memory benchmark')

    op_parser_cor.add_argument('--workdir',
                               type=os.path.abspath,
                               default='gdbr_data',
                               help='program work directory')

    op_parser_cor.add_argument('--overwrite_output',
                               action='store_true',
                               help='overwrite the output directory')

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

    op_parser_ant.add_argument('--benchmark',
                               type=float,
                               nargs='?',
                               const=None,
                               metavar='BENCHMARK_INTERVAL',
                               default=argparse.SUPPRESS,
                               help='turn on benchmark mode, if you specify an argument you can adjust the interval of memory benchmark')

    op_parser_ant.add_argument('--workdir',
                               type=os.path.abspath,
                               default='gdbr_data',
                               help='program work directory')

    op_parser_ant.add_argument('--overwrite_output',
                               action='store_true',
                               help='overwrite the output directory')

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


def check_fai(fasta_loc):
    return os.path.isfile(fasta_loc + '.fai')


def check_query_chrom(qry_loc, ref_chr_list):
    if not check_fai(qry_loc):
        subprocess.run(['samtools', 'faidx', qry_loc], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    qry_seq = Fasta(qry_loc, build_index=False)
    qry_chr_list = list(set(ref_chr_list) & set(qry_seq.records.keys()))
    if qry_chr_list == []:
        raise Exception(f'Wrong query {qry_loc} : Query must have at least one of the same chromosomes as the reference')


def get_min_sv_size(vcf_loc):
    with vcfpy.Reader.from_path(vcf_loc) as record:
        for header in filter(lambda t: type(t) == vcfpy.HeaderLine, record.header.lines):
            if header.key == 'gdbr_min_sv_size':
                return int(header._value)

    return None


def clean_workdir(workdir):
    for item in os.listdir(workdir):
        if os.path.isdir(os.path.join(workdir, item)) and item[-5:] == '_gdbr' and (item[:-5] == 'db' or item[:-5].isdecimal()):
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


def safe_makedirs(targetdir, overwrite_output=False):
    if not os.path.exists(targetdir) or len(os.listdir(targetdir)) == 0 or overwrite_output:
        os.makedirs(targetdir, exist_ok=True)
    else:
        index = 1
        while True:
            new_workdir = f'{targetdir}.{index}'
            if not os.path.exists(new_workdir):
                os.rename(targetdir, new_workdir)
                os.makedirs(targetdir)
                break
            index += 1


def check_reference(ref_loc):
    if not check_fai(ref_loc):
        subprocess.run(['samtools', 'faidx', ref_loc], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    ref_seq = Fasta(ref_loc, build_index=False)
    ref_chr_list = list(map(lambda t: t[0], filter(lambda t: len(t[1]) > 1e6, ref_seq.records.items())))

    if ref_chr_list == []:
        raise Exception(f'Wrong reference {ref_loc} : Reference must have a chromosome larger than 1 Mbp')


def save_fig(fig, savedir, name):
    fig.savefig(os.path.join(savedir, name + '.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(savedir, name + '.pdf'), bbox_inches='tight')


def saferound(v, d=0, md=5):
    for i in range(d, max(md, d) + 1):
        r = round(v, i)
        if r != 0:
            return r if i else int(r)
    return r


def draw_result(savedir, pre_type_cnt, cor_type_cnt, del_type_cnt, ins_type_cnt, sub_type_cnt,
                indel_hom_cnt, temp_ins_hom_cnt, diff_locus_dsbr_hom_cnt, tot_sv_len, query_num):
    # Default setting
    legend_fontsize = 8.5
    width = 0.17

    sns.set_palette('colorblind')
    cb_hex_list = sns.color_palette('colorblind').as_hex()

    code_palette_data = [('DEL', 0),
                        ('INS', 1),
                        ('SUB', 2),
                        ('BND', 6),
                        ('INV', 4),
                        ('REPEAT', 7),
                        ('EXCEPT', 6),
                        ('HOM', 0),
                        ('HOM_GT_SV_90', 1),
                        ('NO_HOM', 6),
                        ('SUB_HOM_DUP', 0),
                        ('SUB_HOM_GT_SV_90', 1),
                        ('DIFF_LOCUS_DSBR', 2),
                        ('SUB_UNIQUE_NO_HOM', 3),
                        ('SUB_REPEAT', 8),
                        ('SUB_NOT_SPECIFIED', 7),
                        ('REPEAT:TRF_FIRST', 7),
                        ('REPEAT:BLAST', 5),
                        ('REPEAT:TRF_LAST', 4),
                        ('REPEAT:RPM', 8),
                        ('SV_COR_ERR', 3),
                        ('SV_SIZE_FILTER', 9)]
    
    code_palette_dict = dict((s, cb_hex_list[p]) for s, p in code_palette_data)

    # Preprocess classification
    target_data = [(k, pre_type_cnt[k] / query_num if k in pre_type_cnt else 0) for k in ['DEL', 'INS', 'BND', 'INV']]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 9.6))
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax1, width=width)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_title('Preprocess variant classification average count')
    ax1.set_xlabel('Avg. count')
    ax1.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val)} / {round(tot_sv_len / query_num)})'for tar, val in target_data])
    target_data = [(k, pre_type_cnt[k] / tot_sv_len * 100 if k in pre_type_cnt else 0) for k in ['DEL', 'INS', 'BND', 'INV']]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax2, width=width)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.set_title('Preprocess variant classification frequency')
    ax2.set_xlabel('Frequency (%)')
    ax2.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val, 1)}%)'for tar, val in target_data])
    save_fig(fig, savedir, 'preprocess_classification')

    # Correct classification
    # verbose figure
    fig, ax = plt.subplots(2, 2, figsize=(6.4 * 2, 9.6))
    cor_type_cnt['EXCEPT'] = sum(cor_type_cnt.get(i, 0) for i in ['FP_SV', 'ERR_POS', 'UNSUP_ID', 'NO_TARGET_CHROM'])
    target_list = ['DEL', 'INS', 'SUB', 'REPEAT:TRF_FIRST', 'REPEAT:BLAST', 'REPEAT:TRF_LAST', 'REPEAT:RPM', 'SV_COR_ERR', 'SV_SIZE_FILTER', 'EXCEPT']
    target_data = [(k, cor_type_cnt[k] / query_num if k in cor_type_cnt else 0) for k in target_list]
    ax1, ax2 = ax[0, 1], ax[1, 1]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax1, width=width)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_title('Corrected variant classification average count')
    ax1.set_xlabel('Avg. count')
    ax1.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val)} / {round(tot_sv_len / query_num)})'for tar, val in target_data])
    target_data = [(k, cor_type_cnt[k] / tot_sv_len * 100 if k in cor_type_cnt else 0) for k in target_list]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax2, width=width)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.set_title('Corrected variant classification frequency')
    ax2.set_xlabel('Frequency (%)')
    ax2.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val, 1)}%)'for tar, val in target_data])

    # silent figure
    cor_type_cnt['REPEAT'] = sum(cor_type_cnt.get(i, 0) for i in ['REPEAT:TRF_FIRST', 'REPEAT:BLAST', 'REPEAT:TRF_LAST', 'REPEAT:RPM'])
    cor_type_cnt['EXCEPT'] = sum(cor_type_cnt.get(i, 0) for i in ['FP_SV', 'ERR_POS', 'UNSUP_ID', 'SV_COR_ERR', 'SV_SIZE_FILTER', 'NO_TARGET_CHROM'])
    target_data = [(k, cor_type_cnt[k] / query_num if k in cor_type_cnt else 0) for k in ['DEL', 'INS', 'SUB', 'REPEAT', 'EXCEPT']]
    ax1, ax2 = ax[0, 0], ax[1, 0]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax1, width=width)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_title('Corrected variant classification average count')
    ax1.set_xlabel('Avg. count')
    ax1.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val)} / {round(tot_sv_len / query_num)})'for tar, val in target_data])
    target_data = [(k, cor_type_cnt[k] / tot_sv_len * 100 if k in cor_type_cnt else 0) for k in ['DEL', 'INS', 'SUB', 'REPEAT', 'EXCEPT']]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax2, width=width)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.set_title('Corrected variant classification frequency')
    ax2.set_xlabel('Frequency (%)')
    ax2.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val, 1)}%)'for tar, val in target_data])
    save_fig(fig, savedir, 'correct_classification')
    
    # Deletion, insertion, indel classification
    cnt = del_type_cnt
    target_list = ['HOM', 'NO_HOM', 'HOM_GT_SV_90']
    target_data = [(k, cnt[k] / query_num if k in cnt else 0) for k in target_list]
    fig, ax = plt.subplots(2, 3, figsize=(6.4 * 3, 9.6))

    ax1, ax2 = ax[0, 0], ax[1, 0]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax1, width=width)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_title('Deletion variant DSBR estimation average count')
    ax1.set_xlabel('Avg. count')
    ax1.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val)} / {round(sum(cnt.values()) / query_num)})'for tar, val in target_data])
    target_data = [(k, cnt[k] / sum(cnt.values()) * 100 if k in cnt else 0) for k in target_list]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax2, width=width)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.set_title('Deletion variant DSBR estimation frequency')
    ax2.set_xlabel('Frequency (%)')
    ax2.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val, 1)}%)'for tar, val in target_data])

    ax1, ax2 = ax[0, 1], ax[1, 1]
    cnt = ins_type_cnt
    target_list = ['HOM', 'NO_HOM', 'HOM_GT_SV_90']
    target_data = [(k, cnt[k] / query_num if k in cnt else 0) for k in target_list]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax1, width=width)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_title('Insertion variant DSBR estimation average count')
    ax1.set_xlabel('Avg. count')
    ax1.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val)} / {round(sum(cnt.values()) / query_num)})'for tar, val in target_data])
    target_data = [(k, cnt[k] / sum(cnt.values()) * 100 if k in cnt else 0) for k in target_list]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax2, width=width)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.set_title('Insertion variant DSBR estimation frequency')
    ax2.set_xlabel('Frequency (%)')
    ax2.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val, 1)}%)'for tar, val in target_data])

    ax1, ax2 = ax[0, 2], ax[1, 2]
    cnt = del_type_cnt + ins_type_cnt
    target_list = ['HOM', 'NO_HOM', 'HOM_GT_SV_90']
    target_data = [(k, cnt[k] / query_num if k in cnt else 0) for k in target_list]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax1, width=width)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_title('Indel variant DSBR estimation average count')
    ax1.set_xlabel('Avg. count')
    ax1.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val)} / {round(sum(cnt.values()) / query_num)})'for tar, val in target_data])
    target_data = [(k, cnt[k] / sum(cnt.values()) * 100 if k in cnt else 0) for k in target_list]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax2, width=width)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.set_title('Indel variant DSBR estimation frequency')
    ax2.set_xlabel('Frequency (%)')
    ax2.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val, 1)}%)'for tar, val in target_data])
    save_fig(fig, savedir, 'result_indel_classification')

    # Substitution classification
    cnt = sub_type_cnt
    target_list = ['SUB_HOM_DUP', 'SUB_HOM_GT_SV_90', 'SUB_UNIQUE_NO_HOM', 'DIFF_LOCUS_DSBR', 'SUB_REPEAT', 'SUB_NOT_SPECIFIED']
    target_data = [(k, cnt[k] / query_num if k in cnt else 0) for k in target_list]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 9.6))
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax1, width=width)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_title('Substitution variant DSBR estimation average count')
    ax1.set_xlabel('Avg. count')
    ax1.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val)} / {round(sum(cnt.values()) / query_num)})'for tar, val in target_data])

    target_data = [(k, cnt[k] / sum(cnt.values()) * 100 if k in cnt else 0) for k in target_list]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax2, width=width)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.set_title('Substitution variant DSBR estimation frequency')
    ax2.set_xlabel('Frequency (%)')
    ax2.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val, 1)}%)'for tar, val in target_data])
    save_fig(fig, savedir, 'result_sub_classification')

    # Indel homology distribution
    if indel_hom_cnt:
        fig, ax_list = plt.subplots(3, 1, figsize=(6.4, 14.4))
        for tar_range, ax in zip([(1, 200), (1, 30), (15, 200)], ax_list):
            s = sns.histplot(x=indel_hom_cnt.keys(), weights=indel_hom_cnt.values(), binrange=tar_range, binwidth=1, element='step', color=cb_hex_list[0], alpha=1, ax=ax)
            s.set(xlabel='(micro)homology (bp)', ylabel='Variant count')
        ax_list[0].set_title('Indel homology distribution')
        save_fig(fig, savedir, 'result_indel_hom_distribution')

    # Templated insertion respective homology distribution
    if temp_ins_hom_cnt:
        fig, ax_list = plt.subplots(2, 1, figsize=(6.4, 9.6))
        for tar_range, ax in zip([(1, 200), (1, 30)], ax_list):
            s = sns.histplot(x=temp_ins_hom_cnt.keys(), weights=temp_ins_hom_cnt.values(), binrange=tar_range, binwidth=1, element='step', color=cb_hex_list[0], alpha=1, ax=ax)
            s.set(xlabel='(micro)homology (bp)', ylabel='Variant count')
        ax_list[0].set_title('Templated insertion respective homology distribution')
        save_fig(fig, savedir, 'result_temp_ins_hom_distribution')

    # Different locus DSBR respective homology distribution
    if diff_locus_dsbr_hom_cnt:
        fig, ax_list = plt.subplots(2, 1, figsize=(6.4, 9.6))
        for tar_range, ax in zip([(1, 2000), (1, 200)], ax_list):
            s = sns.histplot(x=diff_locus_dsbr_hom_cnt.keys(), weights=diff_locus_dsbr_hom_cnt.values(), binrange=tar_range, binwidth=1, element='step', color=cb_hex_list[0], alpha=1, ax=ax)
            s.set(xlabel='(micro)homology (bp)', ylabel='Variant count')
        ax_list[0].set_title('Different locus DSBR respective homology distribution')
        save_fig(fig, savedir, 'result_diff_locus_dsbr_hom_distribution')