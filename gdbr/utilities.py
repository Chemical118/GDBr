from statsmodels.nonparametric.kernel_regression import KernelReg
from matplotlib.ticker import MaxNLocator
from pathos.pools import ProcessPool
from gdbr.version import get_version
from pyfaidx import Fasta

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

import subprocess
import argparse
import datetime
import shutil
import vcfpy
import math
import glob
import os

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


def p_map(f, it, num_cpus=1, hard=True):
    pool = ProcessPool(num_cpus)
    if hard:
        return pool.map(f, it, chunksize=1)
    else:
        return pool.map(f, it)


def logprint(s):
    print(f'[{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] ' + s, flush=True)


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

    ao_parser_anl.add_argument('--diff_locus_dsbr_analysis',
                               action='store_true',
                               help='turn on different locus DSBR analysis (Warning : analysis can give false positives due to partial homology on the sex chromosomes)')

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

    op_parser_ant.add_argument('--diff_locus_dsbr_analysis',
                               action='store_true',
                               help='turn on different locus DSBR analysis (Warning : analysis can give false positives due to partial homology on the sex chromosomes)')

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


def extract_number(s):
    i = len(s) - 1
    while i >= 0:
        if not s[i].isdigit():
            break
        i -= 1
    
    if i < len(s) - 1:
        return int(s[i + 1:])
    else:
        return None


def get_chrom_order(s):
    n = extract_number(s)
    return n is None, n, s


def draw_result(savedir, pre_type_cnt, cor_type_cnt, del_type_cnt, ins_type_cnt, sub_type_cnt,
                indel_hom_cnt, temp_ins_hom_cnt, temp_ins_seq_loc_cnt, temp_ins_seq_len_cnt, temp_ins_del_len_cnt, diff_locus_dsbr_hom_cnt, tot_sv_len, diff_locus_dsbr_analysis, ref_seq, merge_bed_df):
    legend_fontsize = 8.5
    width = 0.17
    linewidth = 1

    sns.set_palette('colorblind')
    cb_hex_list = sns.color_palette('colorblind').as_hex()

    code_palette_data = [('DEL', 0),
                        ('INS', 1),
                        ('SUB', 2),
                        ('BND', 6),
                        ('INV', 4),
                        ('REPEAT', 7),
                        ('EXCEPT', 6),
                        ('ETC_EXCEPT', 6),
                        ('HOM', 0),
                        ('HOM_GT_SV_90', 6),
                        ('NO_HOM', 1),
                        ('TEMP_INS', 0),
                        ('TEMP_INS_HOM_GT_SV_90', 6),
                        ('DIFF_LOCUS_DSBR', 2),
                        ('SUB_UNIQUE_NO_HOM', 3),
                        ('SUB_REPEAT', 8),
                        ('SUB_NOT_SPECIFIED', 1),
                        ('REPEAT:TRF_FIRST', 7),
                        ('REPEAT:BLAST', 5),
                        ('REPEAT:TRF_LAST', 4),
                        ('REPEAT:RPM', 8),
                        ('SV_COR_ERR', 3),
                        ('SV_SIZE_FILTER', 9),
                        ('REF', 0),
                        ('QRY', 1)]
    
    code_palette_dict = dict((s, cb_hex_list[p]) for s, p in code_palette_data)

    # Preprocess classification
    target_data = [(k, pre_type_cnt[k] if k in pre_type_cnt else 0) for k in ['DEL', 'INS', 'BND', 'INV']]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 9.6))
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax1, width=width)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_title('Preprocess variant classification count')
    ax1.set_xlabel('Count')
    ax1.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({val} / {tot_sv_len})'for tar, val in target_data])
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
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
    cor_type_cnt['ETC_EXCEPT'] = sum(cor_type_cnt.get(i, 0) for i in ['FP_SV', 'ERR_POS', 'UNSUP_ID', 'NO_TARGET_CHROM'])
    target_list = ['DEL', 'INS', 'SUB', 'REPEAT:TRF_FIRST', 'REPEAT:BLAST', 'REPEAT:TRF_LAST', 'REPEAT:RPM', 'SV_COR_ERR', 'SV_SIZE_FILTER', 'ETC_EXCEPT']
    target_data = [(k, cor_type_cnt[k] if k in cor_type_cnt else 0) for k in target_list]
    ax1, ax2 = ax[0, 1], ax[1, 1]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax1, width=width)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_title('Corrected variant classification count')
    ax1.set_xlabel('Count')
    ax1.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), prop={'size': legend_fontsize}, labels=[f'{tar} ({val} / {tot_sv_len})'for tar, val in target_data])
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    target_data = [(k, cor_type_cnt[k] / tot_sv_len * 100 if k in cor_type_cnt else 0) for k in target_list]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax2, width=width)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.set_title('Corrected variant classification frequency')
    ax2.set_xlabel('Frequency (%)')
    ax2.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val, 1)}%)'for tar, val in target_data])

    # silent figure
    cor_type_cnt['REPEAT'] = sum(cor_type_cnt.get(i, 0) for i in ['REPEAT:TRF_FIRST', 'REPEAT:BLAST', 'REPEAT:TRF_LAST', 'REPEAT:RPM'])
    cor_type_cnt['EXCEPT'] = sum(cor_type_cnt.get(i, 0) for i in ['FP_SV', 'ERR_POS', 'UNSUP_ID', 'SV_COR_ERR', 'SV_SIZE_FILTER', 'NO_TARGET_CHROM'])
    target_data = [(k, cor_type_cnt[k] if k in cor_type_cnt else 0) for k in ['DEL', 'INS', 'SUB', 'REPEAT', 'EXCEPT']]
    ax1, ax2 = ax[0, 0], ax[1, 0]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax1, width=width)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_title('Corrected variant classification count')
    ax1.set_xlabel('Count')
    ax1.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({val} / {tot_sv_len})'for tar, val in target_data])
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
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
    target_data = [(k, cnt[k] if k in cnt else 0) for k in target_list]
    fig, ax = plt.subplots(2, 3, figsize=(6.4 * 3, 9.6))

    ax1, ax2 = ax[0, 0], ax[1, 0]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax1, width=width)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_title('Deletion variant DSBR estimation count')
    ax1.set_xlabel('Count')
    ax1.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({val} / {sum(cnt.values())})' for tar, val in target_data])
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    target_data = [(k, cnt[k] / sum(cnt.values()) * 100 if k in cnt else 0) for k in target_list]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax2, width=width)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.set_title('Deletion variant DSBR estimation frequency')
    ax2.set_xlabel('Frequency (%)')
    ax2.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val, 1)}%)'for tar, val in target_data])

    ax1, ax2 = ax[0, 1], ax[1, 1]
    cnt = ins_type_cnt
    target_list = ['HOM', 'NO_HOM', 'HOM_GT_SV_90']
    target_data = [(k, cnt[k] if k in cnt else 0) for k in target_list]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax1, width=width)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_title('Insertion variant DSBR estimation count')
    ax1.set_xlabel('Count')
    ax1.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({val} / {sum(cnt.values())})' for tar, val in target_data])
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    target_data = [(k, cnt[k] / sum(cnt.values()) * 100 if k in cnt else 0) for k in target_list]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax2, width=width)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.set_title('Insertion variant DSBR estimation frequency')
    ax2.set_xlabel('Frequency (%)')
    ax2.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val, 1)}%)'for tar, val in target_data])

    ax1, ax2 = ax[0, 2], ax[1, 2]
    cnt = del_type_cnt + ins_type_cnt
    target_list = ['HOM', 'NO_HOM', 'HOM_GT_SV_90']
    target_data = [(k, cnt[k] if k in cnt else 0) for k in target_list]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax1, width=width)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_title('Indel variant DSBR estimation count')
    ax1.set_xlabel('Count')
    ax1.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({val} / {sum(cnt.values())})' for tar, val in target_data])
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    target_data = [(k, cnt[k] / sum(cnt.values()) * 100 if k in cnt else 0) for k in target_list]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax2, width=width)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.set_title('Indel variant DSBR estimation frequency')
    ax2.set_xlabel('Frequency (%)')
    ax2.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val, 1)}%)'for tar, val in target_data])
    save_fig(fig, savedir, 'result_indel_classification')

    # Substitution classification
    cnt = sub_type_cnt
    target_list = ['TEMP_INS', 'DIFF_LOCUS_DSBR', 'SUB_UNIQUE_NO_HOM', 'SUB_NOT_SPECIFIED', 'SUB_REPEAT', 'TEMP_INS_HOM_GT_SV_90'] if diff_locus_dsbr_analysis else ['TEMP_INS', 'SUB_NOT_SPECIFIED', 'TEMP_INS_HOM_GT_SV_90']
    target_data = [(k, cnt[k] if k in cnt else 0) for k in target_list]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 9.6))
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax1, width=width)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_title('Substitution variant DSBR estimation count')
    ax1.set_xlabel('Count')
    ax1.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({val} / {sum(cnt.values())})' for tar, val in target_data])
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    target_data = [(k, cnt[k] / sum(cnt.values()) * 100 if k in cnt else 0) for k in target_list]
    pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax2, width=width)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.set_title('Substitution variant DSBR estimation frequency')
    ax2.set_xlabel('Frequency (%)')
    ax2.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val, 1)}%)'for tar, val in target_data])
    save_fig(fig, savedir, 'result_sub_classification')

    # Indel homology distribution
    if indel_hom_cnt:
        indel_hom_cnt_data = np.array(list(filter(lambda t: 10 <= t[0] <= 60, indel_hom_cnt.items())))

        x, y = indel_hom_cnt_data[:, 0], indel_hom_cnt_data[:, 1]
        kernel_reg = KernelReg(y, x, var_type='c', bw=[2])

        pre_gap = 0.1
        pre_x = np.arange(np.min(x), np.max(x) + 1, pre_gap)
        pre_y = kernel_reg.fit(pre_x)[0]

        slp_y = np.diff(pre_y)
        slp_x = pre_x[:-1] + pre_gap / 2

        ans_x, pos_x, ans_range = None, None, 1
        for i, (sx, sy) in enumerate(zip(slp_x, slp_y)):
            if pos_x is None:
                if sy > 0:
                    pos_x = sx
            else:
                if sy < 0 or i == np.size(slp_x) - 1:
                    if sx - pos_x > ans_range:
                        ans_x = pos_x
                        ans_range = sx - pos_x
                    pos_x = None

        a_ej_baseline = None if ans_x is None else math.trunc(ans_x)

        fig, ax_list = plt.subplots(3, 1, figsize=(6.4, 14.4))
        for tar_range, ax in zip([(1, 200), (1, 30)], ax_list):
            s = sns.histplot(x=indel_hom_cnt.keys(), weights=indel_hom_cnt.values(), binrange=tar_range, binwidth=1, element='step', alpha=1, ax=ax)
            s.set(xlabel='micro/homology (bp)', ylabel='Variant count')

        ax = ax_list[2]
        s = sns.histplot(x=indel_hom_cnt.keys(), weights=indel_hom_cnt.values(), binrange=(15, 200), binwidth=1, element='step', alpha=1, ax=ax)
        s.set(xlabel='micro/homology (bp)', ylabel='Variant count')
        if a_ej_baseline is not None:
            ax.text(a_ej_baseline - 2, ax.get_ylim()[1] * 0.99, 'TMEJ',
                    color='red',
                    horizontalalignment='right',
                    verticalalignment='top')

            ax.text(a_ej_baseline + 2, ax.get_ylim()[1] * 0.99, 'SSA',
                    color='red',
                    horizontalalignment='left',
                    verticalalignment='top')

            ax.text(a_ej_baseline, ax.get_ylim()[1], f'{a_ej_baseline}bp',
                    color='red',
                    horizontalalignment='center',
                    verticalalignment='bottom')

            ax.axvline(x=a_ej_baseline, color='red', linestyle='dashed')

        ax_list[0].set_title('Indel homology distribution')
        save_fig(fig, savedir, 'result_indel_hom_distribution')

    # Templated insertion respective homology distribution
    if temp_ins_hom_cnt:
        fig, ax_list = plt.subplots(2, 1, figsize=(6.4, 9.6))
        for tar_range, ax in zip([(1, 200), (1, 30)], ax_list):
            s = sns.histplot(x=temp_ins_hom_cnt.keys(), weights=temp_ins_hom_cnt.values(), binrange=tar_range, binwidth=1, element='step', color=cb_hex_list[0], alpha=1, ax=ax)
            s.set(xlabel='micro/homology (bp)', ylabel='Variant count')
        ax_list[0].set_title('Templated insertion respective homology distribution')
        save_fig(fig, savedir, 'result_temp_ins_hom_distribution')

        fig, ax_list = plt.subplots(2, 2, figsize=(6.4 * 2, 9.6))

        ax1, ax2 = ax_list[0, 0], ax_list[1, 0]
        cnt = temp_ins_seq_loc_cnt
        target_list = ['REF', 'QRY']
        target_data = [(k, temp_ins_seq_loc_cnt[k] if k in cnt else 0) for k in target_list]
        pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax1, width=width)
        ax1.axes.get_yaxis().set_visible(False)
        ax1.set_title('Templated insertion location count')
        ax1.set_xlabel('Count')
        ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax1.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({val} / {sum(cnt.values())})'for tar, (_, val) in zip(['Reference', 'Query'], target_data)])
        target_data = [(k, cnt[k] / sum(cnt.values()) * 100 if k in cnt else 0) for k in target_list]
        pd.DataFrame(dict(target_data), index=['']).plot.barh(color=code_palette_dict, stacked=True, ax=ax2, width=width)
        ax2.axes.get_yaxis().set_visible(False)
        ax2.set_title('Templated insertion location frequency')
        ax2.set_xlabel('Frequency (%)')
        ax2.legend(loc=2, prop={'size': legend_fontsize}, labels=[f'{tar} ({saferound(val, 1)}%)'for tar, (_, val) in zip(['Reference', 'Query'], target_data)])


        ax1, ax2 = ax_list[0, 1], ax_list[1, 1]
        x = list(temp_ins_seq_len_cnt.keys()) + list(temp_ins_del_len_cnt.keys())
        weights = list(temp_ins_seq_len_cnt.values()) + list(temp_ins_del_len_cnt.values())
        hue = ['Templated insertion'] * len(temp_ins_seq_len_cnt) + ['Deletion'] * len(temp_ins_del_len_cnt)
        sns.histplot(x=x, weights=weights, hue=hue, binrange=(1, 50), binwidth=1, element='step', alpha=1, fill=False, ax=ax1, linewidth=linewidth)
        sns.move_legend(ax1, loc=None, fontsize=legend_fontsize)
        ax1.set_title('Templated insertion - deletion length distribution')
        ax1.set_xlabel('length (bp)')

        sns.histplot(x=x, weights=weights, hue=hue, binrange=(1, 200), binwidth=1, element='step', alpha=1, fill=False, ax=ax2, linewidth=linewidth)
        sns.move_legend(ax2, loc=None, fontsize=legend_fontsize)
        ax2.set_xlabel('length (bp)')
        save_fig(fig, savedir, 'result_temp_ins_sequence_data')

    # Different locus DSBR respective homology distribution
    if diff_locus_dsbr_hom_cnt:
        fig, ax_list = plt.subplots(2, 1, figsize=(6.4, 9.6))
        for tar_range, ax in zip([(1, 2000), (1, 200)], ax_list):
            s = sns.histplot(x=diff_locus_dsbr_hom_cnt.keys(), weights=diff_locus_dsbr_hom_cnt.values(), binrange=tar_range, binwidth=1, element='step', color=cb_hex_list[0], alpha=1, ax=ax)
            s.set(xlabel='micro/homology (bp)', ylabel='Variant count')
        ax_list[0].set_title('Different locus DSBR respective homology distribution')
        save_fig(fig, savedir, 'result_diff_locus_dsbr_hom_distribution')
    
    # Draw chromosome distribution
    merge_bed_df.columns = ['chrom', 'st', 'nd', 'gdbr_type', 'hom_l', 'hom_r', 'merge_id']
    merge_bed_df['sv_loc'] = [round(i) for i in (merge_bed_df['st'] + merge_bed_df['nd']) / 2]
    merge_bed_df['Homology type'] = ['Homology' if i in {'HOM', 'TEMP_INS'} else 'No homology' for i in merge_bed_df['gdbr_type']]
    if a_ej_baseline is not None:
        dsb_repair_type_list = []
        for gdbr_type, hom in zip(merge_bed_df['gdbr_type'], merge_bed_df['hom_l']):
            if gdbr_type == 'TEMP_INS':
                dsb_repair_type_list.append('TMEJ')
            elif gdbr_type == 'HOM':
                if hom <= a_ej_baseline:
                    dsb_repair_type_list.append('TMEJ')
                else:
                    dsb_repair_type_list.append('SSA')
            else:
                dsb_repair_type_list.append('No homology')

        merge_bed_df['DSB repair type'] = dsb_repair_type_list

    chr_list = sorted(set(merge_bed_df['chrom']), key=get_chrom_order)
    chr_res_df_data = [(merge_bed_df.query(f'chrom == "{chrom}"'), chrom, len(ref_seq[chrom])) for chrom in chr_list]

    binwidth = round(sum(len(ref_seq[chrom]) for chrom in chr_list) / 5000)
    nrows, ncols = len(chr_res_df_data), 2

    fig, ax_array = plt.subplots(nrows=nrows + 1, ncols=ncols, figsize=(12 * ncols, 0.02 + 4 * nrows), gridspec_kw={"height_ratios" : [0.02] + [4] * nrows}, layout="constrained")

    fig.suptitle('Homology - no homology chromosome distribution', fontweight='bold')

    ax = ax_array[0, 0]
    ax.axis('off')
    ax.set_title('Variant count histogram')

    ax = ax_array[0, 1]
    ax.axis('off')
    ax.set_title('Varaiant normalized count histogram')

    hue, hue_order = 'Homology type', ['Homology', 'No homology']
    for i, (tdf, chrom, chrom_len) in enumerate(chr_res_df_data):
        ax = ax_array[i + 1, 0]
        sns.histplot(data=tdf, x='sv_loc', hue=hue, binwidth=binwidth, ax=ax, element="step", fill=False, hue_order=hue_order, linewidth=linewidth, alpha=1)
        ax.set(title=chrom, xlim=(0, len(ref_seq[chrom])), xlabel='Variant location (bp)', ylabel='Variant count')
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        if i == 0:
            sns.move_legend(ax, loc='lower right', bbox_to_anchor=(1, 1))
        else:
            ax.get_legend().remove()
        
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        ax = ax_array[i + 1, 1]
        sns.histplot(data=tdf, x='sv_loc', hue=hue, binwidth=binwidth, ax=ax, element="step", fill=False, hue_order=hue_order, linewidth=linewidth, alpha=1)
        ax.set(title=chrom, xlim=(0, len(ref_seq[chrom])), ylim=(0, 1.05), xlabel='Variant location (bp)', ylabel='Varaiant normalized count')
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        if i == 0:
            sns.move_legend(ax, loc='lower right', bbox_to_anchor=(1, 1))
        else:
            ax.get_legend().remove()

        for line2d in ax.lines:
            y = line2d.get_ydata()
            line2d.set_ydata(y / np.max(y))

    save_fig(fig, savedir, 'result_chrom_hom_no_hom_distribution')

    if a_ej_baseline is not None:
        fig, ax_array = plt.subplots(nrows=nrows + 1, ncols=ncols, figsize=(12 * ncols, 0.02 + 4 * nrows), gridspec_kw={"height_ratios" : [0.02] + [4] * nrows}, layout="constrained")
        fig.suptitle('TMEJ - no homology chromosome distribution', fontweight='bold')

        ax = ax_array[0, 0]
        ax.axis('off')
        ax.set_title('Variant count histogram')

        ax = ax_array[0, 1]
        ax.axis('off')
        ax.set_title('Varaiant normalized count histogram')

        hue, hue_order = 'DSB repair type', ['TMEJ', 'No homology']
        for i, (tdf, chrom, chrom_len) in enumerate(chr_res_df_data):
            ax = ax_array[i + 1, 0]
            sns.histplot(data=tdf, x='sv_loc', hue=hue, binwidth=binwidth, ax=ax, element="step", fill=False, hue_order=hue_order, linewidth=linewidth, alpha=1)
            ax.set(title=chrom, xlim=(0, chrom_len), xlabel='Variant location (bp)', ylabel='Variant count')
            ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
            if i == 0:
                sns.move_legend(ax, loc='lower right', bbox_to_anchor=(1, 1))
            else:
                ax.get_legend().remove()
            
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))

            ax = ax_array[i + 1, 1]
            sns.histplot(data=tdf, x='sv_loc', hue=hue, binwidth=binwidth, ax=ax, element="step", fill=False, hue_order=hue_order, linewidth=linewidth, alpha=1)
            ax.set(title=chrom, xlim=(0, chrom_len), ylim=(0, 1.05), xlabel='Variant location (bp)', ylabel='Varaiant normalized count')
            ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
            if i == 0:
                sns.move_legend(ax, loc='lower right', bbox_to_anchor=(1, 1))
            else:
                ax.get_legend().remove()

            for line2d in ax.lines:
                y = line2d.get_ydata()
                line2d.set_ydata(y / np.max(y))

        save_fig(fig, savedir, 'result_chrom_tmej_no_hom_distribution')

        fig, ax_array = plt.subplots(nrows=nrows + 1, ncols=ncols, figsize=(12 * ncols, 0.02 + 4 * nrows), gridspec_kw={"height_ratios" : [0.02] + [4] * nrows}, layout="constrained")
        fig.suptitle('TMEJ - no homology - SSA chromosome distribution', fontweight='bold')

        ax = ax_array[0, 0]
        ax.axis('off')
        ax.set_title('Variant count histogram')

        ax = ax_array[0, 1]
        ax.axis('off')
        ax.set_title('Varaiant normalized count histogram')

        hue, hue_order = 'DSB repair type', ['TMEJ', 'No homology', 'SSA']
        for i, (tdf, chrom, chrom_len) in enumerate(chr_res_df_data):
            ax = ax_array[i + 1, 0]
            sns.histplot(data=tdf, x='sv_loc', hue=hue, binwidth=binwidth, ax=ax, element="step", fill=False, hue_order=hue_order, linewidth=linewidth, alpha=1)
            ax.set(title=chrom, xlim=(0, chrom_len), xlabel='Variant location (bp)', ylabel='Variant count')
            ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
            if i == 0:
                sns.move_legend(ax, loc='lower right', bbox_to_anchor=(1, 1))
            else:
                ax.get_legend().remove()
            
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))

            ax = ax_array[i + 1, 1]
            sns.histplot(data=tdf, x='sv_loc', hue=hue, binwidth=binwidth, ax=ax, element="step", fill=False, hue_order=hue_order, linewidth=linewidth, alpha=1)
            ax.set(title=chrom, xlim=(0, chrom_len), ylim=(0, 1.05), xlabel='Variant location (bp)', ylabel='Varaiant normalized count')
            ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
            if i == 0:
                sns.move_legend(ax, loc='lower right', bbox_to_anchor=(1, 1))
            else:
                ax.get_legend().remove()

            for line2d in ax.lines:
                y = line2d.get_ydata()
                line2d.set_ydata(y / np.max(y))

        save_fig(fig, savedir, 'result_chrom_tmej_no_hom_ssa_distribution')

    return a_ej_baseline