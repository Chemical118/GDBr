from gdbr.utilities import p_map, logprint, check_file_exist, check_unique_basename, check_variant_caller, get_min_sv_size, remove_gdbr_postfix, safe_makedirs, check_query_chrom
from functools import partial
from pyfaidx import Fasta

import numpy as np

import subprocess
import shutil
import vcfpy
import csv
import os

#                [0]       [1]        [2]        [3]        [4]      [5]      [6]      [7]       [8]
blastn_fmt = ['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'bitscore']
trf_fmt = ['2', '5', '7', '80', '10', '50', '500', '-ngs', '-h']

def get_blast_result(start_temp_seq, end_temp_seq, sv_find_len, sv_id, chrom, dbdir, qryworkdir):
    start_temp_seq_loc = os.path.join(qryworkdir, f'ref_start_{sv_id}.fasta')
    end_temp_seq_loc = os.path.join(qryworkdir, f'ref_end_{sv_id}.fasta')

    with open(start_temp_seq_loc, 'w') as f:
        f.write('>' + start_temp_seq.fancy_name + '\n')
        f.write(str(start_temp_seq))

    with open(end_temp_seq_loc, 'w') as f:
        f.write('>' + end_temp_seq.fancy_name + '\n')
        f.write(str(end_temp_seq))

    start_result = subprocess.run(['blastn', '-strand', 'plus', '-db', os.path.join(dbdir, chrom), '-query', start_temp_seq_loc, '-outfmt', '10 ' + ' '.join(blastn_fmt)], capture_output=True, text=True)
    end_result = subprocess.run(['blastn', '-strand', 'plus', '-db', os.path.join(dbdir, chrom), '-query', end_temp_seq_loc, '-outfmt', '10 ' + ' '.join(blastn_fmt)], capture_output=True, text=True)

    start_result_list =  [] if start_result.stdout == '' else map(lambda t: list(map(eval, t.split(','))), start_result.stdout[:-1].split('\n'))
    end_result_list = [] if end_result.stdout == '' else map(lambda t: list(map(eval, t.split(','))), end_result.stdout[:-1].split('\n'))

    os.remove(start_temp_seq_loc)
    os.remove(end_temp_seq_loc)

    start_filter_result = list(filter(lambda t: t[1] > sv_find_len * 0.9, start_result_list))
    end_filter_result = list(filter(lambda t: t[1] > sv_find_len * 0.9, end_result_list))

    return start_filter_result, end_filter_result


def get_blast_single_result(ref_temp_seq, qry_temp_seq, sv_find_len, sv_id, qryworkdir, filter_func):
    ref_temp_seq_loc = os.path.join(qryworkdir, f'sub_ref_{sv_id}.fasta')
    qry_temp_seq_loc = os.path.join(qryworkdir, f'sub_qry_{sv_id}.fasta')

    with open(ref_temp_seq_loc, 'w') as f:
        f.write('>' + ref_temp_seq.fancy_name + '\n')
        f.write(str(ref_temp_seq))

    with open(qry_temp_seq_loc, 'w') as f:
        f.write('>' + qry_temp_seq.fancy_name + '\n')
        f.write(str(qry_temp_seq))

    sub_result = subprocess.run(['blastn', '-subject', qry_temp_seq_loc, '-query', ref_temp_seq_loc, '-strand', 'plus', '-outfmt', '10 ' + ' '.join(blastn_fmt)], capture_output=True, text=True)
    sub_result_list =  [] if sub_result.stdout == '' else map(lambda t: list(map(eval, t.split(','))), sub_result.stdout[:-1].split('\n'))
    sub_filter_result = sorted(filter(lambda t: t[1] > sv_find_len * 0.90 and filter_func(t), sub_result_list), key=lambda t: t[8], reverse=True)

    os.remove(ref_temp_seq_loc)
    os.remove(qry_temp_seq_loc)

    return sub_filter_result


def get_corrected_location(ref_start, ref_end, sv_find_len, start_filter_result, end_filter_result):
    ref_start -= sv_find_len - start_filter_result[0][5]
    ref_end += end_filter_result[0][4] - 1
    ref_len = ref_end - ref_start + 1

    qry_start = start_filter_result[0][7] + 1
    qry_end = end_filter_result[0][6] - 1
    qry_len = qry_end - qry_start + 1
    return ref_start, ref_end, ref_len, qry_start, qry_end, qry_len


def get_correctrd_location_by_idx(ref_temp_seq, qry_temp_seq, ref_start, ref_end, qry_start, qry_end):
    base_ref_start = ref_start
    base_qry_start = qry_start

    while ref_start < ref_end + 1 and qry_start < qry_end and ref_temp_seq[ref_start - base_ref_start] == qry_temp_seq[qry_start - base_qry_start]:
        ref_start += 1
        qry_start += 1

    while ref_start - 1 < ref_end and qry_start - 1 < qry_end and ref_temp_seq[ref_end - base_ref_start] == qry_temp_seq[qry_end - base_qry_start]:
        ref_end -= 1
        qry_end -= 1

    ref_len = ref_end - ref_start + 1
    qry_len = qry_end - qry_start + 1
    return ref_start, ref_end, ref_len, qry_start, qry_end, qry_len


def trf_repeat_check_first(ref_seq, ref_start, ref_end, sv_id, chrom, qryworkdir, sv_seq, repeat_find_len=50):
    trf_temp_seq_loc = os.path.join(qryworkdir, f'trf_first_{sv_id}.fasta')
    start_temp_seq = ref_seq[chrom][ref_start - repeat_find_len - 1:ref_start - 1]
    end_temp_seq = ref_seq[chrom][ref_end:ref_end + repeat_find_len]

    with open(trf_temp_seq_loc, 'w') as f:
        f.write(f'>{sv_id}.ref\n')
        f.write(str(start_temp_seq) + sv_seq + str(end_temp_seq))
    
    trf_result = subprocess.run(['trf', trf_temp_seq_loc] + trf_fmt, capture_output=True, text=True)

    os.remove(trf_temp_seq_loc)
    return trf_result.stdout != ''

def trf_repeat_check_last(ref_seq, qry_seq, ref_start, ref_end, ref_len, qry_start, qry_end, qry_len, sv_id, chrom, qryworkdir, repeat_find_len=50):
    trf_temp_seq_loc = os.path.join(qryworkdir, f'trf_last_{sv_id}.fasta')

    with open(trf_temp_seq_loc, 'w') as f:
        if ref_len > 0:
            ref_temp_seq = ref_seq[chrom][ref_start - repeat_find_len - 1:ref_end + repeat_find_len]
            f.write(f'>{sv_id}.ref\n')
            f.write(str(ref_temp_seq) + '\n\n')

        if qry_len > 0:
            qry_temp_seq = qry_seq[chrom][qry_start - repeat_find_len - 1:qry_end + repeat_find_len]
            f.write(f'>{sv_id}.qry\n')
            f.write(str(qry_temp_seq) + '\n\n')
    
    trf_result = subprocess.run(['trf', trf_temp_seq_loc] + trf_fmt, capture_output=True, text=True)
    
    os.remove(trf_temp_seq_loc)
    return trf_result.stdout != ''


def get_real_sv(record_data, ref_loc, qry_loc, dbdir, qryworkdir, sv_find_len=2000, repeat_find_len=50, min_sv_size=50):
    sv_id, cal_id, record, is_tar = record_data

    if not is_tar:
        return [sv_id, cal_id, 'NO_TARGET_CHROM']
    
    if cal_id not in {'DEL', 'INS'}:
        return [sv_id, cal_id, 'UNSUP_ID']

    ref_seq = Fasta(ref_loc, build_index=False)
    qry_seq = Fasta(qry_loc, build_index=False)

    ref_start, ref_end, sv_len, chrom = record.POS + 1, record.INFO['END'], abs(record.INFO['SVLEN']), record.CHROM

    call_sv_seq = record.REF[1:] if cal_id == 'DEL' else record.ALT[0].value[1:]
    if trf_repeat_check_first(ref_seq, ref_start, ref_end, sv_id, chrom, qryworkdir, call_sv_seq, repeat_find_len=repeat_find_len):
        return [sv_id, cal_id, 'REPEAT:TRF_FIRST']

    start_temp_seq = ref_seq[chrom][ref_start - sv_find_len - 1:ref_start - 1]
    end_temp_seq = ref_seq[chrom][ref_end:ref_end + sv_find_len]

    start_filter_result, end_filter_result = get_blast_result(start_temp_seq, end_temp_seq, sv_find_len, sv_id, chrom, dbdir, qryworkdir)
    
    if len(start_filter_result) != 1 or len(end_filter_result) != 1:
        return [sv_id, cal_id, 'REPEAT:BLAST']

    ref_start, ref_end, ref_len, qry_start, qry_end, qry_len = get_corrected_location(ref_start, ref_end, sv_find_len, start_filter_result, end_filter_result)

    # remove false substitution
    call_err_try = 3
    while call_err_try > 0 and qry_len < 0 and start_filter_result[0][6] < end_filter_result[0][6] and start_filter_result[0][7] < end_filter_result[0][7]:
        ref_end -= qry_len
        qry_end -= qry_len
        
        start_temp_seq = ref_seq[chrom][ref_start - sv_find_len - 1:ref_start - 1]
        end_temp_seq = ref_seq[chrom][ref_end:ref_end + sv_find_len]

        start_filter_result, end_filter_result = get_blast_result(start_temp_seq, end_temp_seq, sv_find_len, sv_id, chrom, dbdir, qryworkdir)
        
        if len(start_filter_result) != 1 or len(end_filter_result) != 1:
            return [sv_id, cal_id, 'REPEAT:BLAST']
        
        ref_start, ref_end, ref_len, qry_start, qry_end, qry_len = get_corrected_location(ref_start, ref_end, sv_find_len, start_filter_result, end_filter_result)

        call_err_try -= 1
    
    # substitution blast start test
    if ref_len > 0 and qry_len > 0:
        check_len = min(ref_len, qry_len)

        start_ref_check_len = check_len + (min(abs(ref_len - qry_len), 10) if ref_len > check_len else 0)
        start_qry_check_len = check_len + (min(abs(ref_len - qry_len), 10) if qry_len > check_len else 0)

        start_ref_temp_seq = ref_seq[chrom][ref_start - sv_find_len - 1:ref_start - 1 + start_ref_check_len]
        start_qry_temp_seq = qry_seq[chrom][qry_start - sv_find_len - 1:qry_start - 1 + start_qry_check_len]

        start_sub_filter_result = get_blast_single_result(start_ref_temp_seq, start_qry_temp_seq, sv_find_len, sv_id, qryworkdir, lambda t: t[5] >= sv_find_len and t[7] >= sv_find_len)

        if len(start_sub_filter_result) > 0:
            ref_start += start_sub_filter_result[0][5] - sv_find_len
            qry_start += start_sub_filter_result[0][7] - sv_find_len

            ref_len = ref_end - ref_start + 1
            qry_len = qry_end - qry_start + 1

    # substitution blast end test
    if ref_len > 0 and qry_len > 0:
        check_len = min(ref_len, qry_len)

        end_ref_check_len = check_len + (min(abs(ref_len - qry_len), 10) if ref_len > check_len else 0)
        end_qry_check_len = check_len + (min(abs(ref_len - qry_len), 10) if qry_len > check_len else 0)
        
        end_ref_temp_seq = ref_seq[chrom][ref_end - end_ref_check_len:ref_end + sv_find_len]
        end_qry_temp_seq = qry_seq[chrom][qry_end - end_qry_check_len:qry_end + sv_find_len]

        end_sub_filter_result = get_blast_single_result(end_ref_temp_seq, end_qry_temp_seq, sv_find_len, sv_id, qryworkdir, lambda t: t[4] <= end_ref_check_len + 1 and t[6] <= end_qry_check_len + 1)

        if len(end_sub_filter_result) > 0:
            ref_end -= end_ref_check_len + 1 - end_sub_filter_result[0][4]
            qry_end -= end_qry_check_len + 1 - end_sub_filter_result[0][6]

            ref_len = ref_end - ref_start + 1
            qry_len = qry_end - qry_start + 1

    # substitution index test
    if ref_len > 0 and qry_len > 0:
        ref_temp_seq = str(ref_seq[chrom][ref_start - 1:ref_end]).upper()
        qry_temp_seq = str(qry_seq[chrom][qry_start - 1:qry_end]).upper()

        ref_start, ref_end, ref_len, qry_start, qry_end, qry_len = get_correctrd_location_by_idx(ref_temp_seq, qry_temp_seq, ref_start, ref_end, qry_start, qry_end)
    
    real_sv_len = max(ref_len, qry_len)
    if qry_len < 0:
        return [sv_id, cal_id, 'ERR_POS']

    elif real_sv_len < min_sv_size:
        sv_name = 'SV_SIZE_FILTER'

    elif real_sv_len > sv_len * 3:
        return [sv_id, cal_id, 'SV_COR_ERR']
    
    elif trf_repeat_check_last(ref_seq, qry_seq, ref_start, ref_end, ref_len, qry_start, qry_end, qry_len, sv_id, chrom, qryworkdir, repeat_find_len=repeat_find_len):
        return [sv_id, cal_id, 'REPEAT:TRF_LAST']
    
    elif ref_len == 0 and qry_len == 0:
        return [sv_id, cal_id, 'FP_SV']

    elif ref_len > 0 and qry_len > 0:
        sv_name = 'SUB'

    elif ref_len == 0:
        sv_name = 'INS'
    
    else:
        sv_name = 'DEL'

    return [sv_id, cal_id, sv_name, chrom, ref_start, ref_end, qry_start, qry_end]  


def makeblastdb_from_location(chr_name, seq_loc, dbdir):
    seq = Fasta(seq_loc, build_index=False)
    with open(os.path.join(dbdir, chr_name + '.fasta'), 'w') as f:
        tseq = seq[chr_name]
        f.write('>' + chr_name + '\n')
        f.write(str(tseq))
    
    subprocess.run(['makeblastdb', '-in', os.path.join(dbdir, chr_name + '.fasta'), '-input_type', 'fasta', '-dbtype', 'nucl', '-out', os.path.join(dbdir, chr_name)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.remove(os.path.join(dbdir, chr_name + '.fasta'))


def get_proper_correct_order(record_data):
    _, cal_id, record, is_tar = record_data

    if not is_tar or cal_id not in {'DEL', 'INS'}:
        return 0
    else:
        return -abs(record.INFO['SVLEN'])


def get_vcf_tot_data(vcf_loc, ref_chr_list):
    with vcfpy.Reader.from_path(vcf_loc) as tot_record:
        vcf_tot_data = [(index, record.ID[0].split('.')[1:-1][0], record, record.CHROM in ref_chr_list) for index, record in enumerate(tot_record)]
    
    return vcf_tot_data


def prepare_repeatmasker_multiprocessing(ref_seq, qry_seq, sv_list, qryworkdir, num_cpus, repeat_find_len):
    rpmworkdir_list = []
    for i, target_ind in enumerate(np.array_split(np.arange(len(sv_list)), num_cpus)):
        if np.size(target_ind) > 0:
            rpmworkdir = os.path.join(qryworkdir, f'{i}_rpm')
            rpmworkdir_list.append(rpmworkdir)
            os.makedirs(rpmworkdir, exist_ok=True)

            rpm_temp_seq_loc = os.path.join(rpmworkdir, 'rpm.fa')
            with open(rpm_temp_seq_loc, 'w') as f:
                for sv_ind in target_ind:
                    sv = sv_list[sv_ind]
                    if sv[2] in {'DEL', 'INS', 'SUB'}:
                        sv_id, _, _, chrom, ref_start, ref_end, qry_start, qry_end = sv

                        ref_len = ref_end - ref_start + 1
                        qry_len = qry_end - qry_start + 1

                        if ref_len > 0:
                            temp_seq = ref_seq[chrom][ref_start - repeat_find_len - 1:ref_end + repeat_find_len]
                            f.write(f'>{sv_id}.ref\n')
                            f.write(str(temp_seq) + '\n\n')

                        if qry_len > 0:
                            temp_seq = qry_seq[chrom][qry_start - repeat_find_len - 1:qry_end + repeat_find_len]
                            f.write(f'>{sv_id}.qry\n')
                            f.write(str(temp_seq) + '\n\n')
    
    return rpmworkdir_list



def get_reapeatmasker_index_list(rpmworkdir, species):
    rpm_temp_seq_loc = os.path.join(rpmworkdir, 'rpm.fa')
    subprocess.run(['RepeatMasker', '-spec', species, '-nopost', rpm_temp_seq_loc], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=rpmworkdir)
        
    rpm_output_loc = rpm_temp_seq_loc + '.cat' if os.path.isfile(rpm_temp_seq_loc + '.cat') else rpm_temp_seq_loc + '.cat.gz'
    subprocess.run(['ProcessRepeats', '-spec', species, '-gff', rpm_output_loc], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=rpmworkdir)
    
    rpm_ind_set = set()
    if os.path.isfile(rpm_temp_seq_loc + '.out.gff'):
        with open(rpm_temp_seq_loc + '.out.gff', 'r') as f:
            tf = csv.reader(f, delimiter='\t')
            rpm_list = [l for l in tf]

        for rpm in rpm_list:
            if rpm[0][0] != '#':
                rpm_ind_set.add(int(rpm[0].split('.')[0]))
    
    shutil.rmtree(rpmworkdir)
    return list(rpm_ind_set)


def get_group_index(tar_list, ind, baseline=10000, reverse=False):
    ans_ind = ind
    walk = -1 if reverse else 1
    while ans_ind + walk >= 0 and ans_ind + walk < len(tar_list) and abs(tar_list[ans_ind + walk] - tar_list[ans_ind]) < baseline:
        ans_ind += walk

    return ans_ind


def get_correct_error_index_list(sv_list):
    err_baseline = 10000
    err_ind_list = []
    real_sv_list = sorted(filter(lambda t: t[2] in {'DEL', 'INS', 'SUB'}, sv_list), key=lambda t: (t[3], t[4]))

    tmp_chrom = None
    tmp_chrom_list = []

    real_chrom_sv_data = []

    for sv in real_sv_list:
        chrom = sv[3]

        if chrom != tmp_chrom:
            if len(tmp_chrom_list) > 0:
                real_chrom_sv_data.append(tmp_chrom_list)
            tmp_chrom, tmp_chrom_list = chrom, [sv]
        else:
            tmp_chrom_list.append(sv)

    real_chrom_sv_data.append(tmp_chrom_list)

    for real_chrom_sv_list in real_chrom_sv_data:
        if len(real_chrom_sv_list) > 1:
            real_chrom_qry_st = [i[6] for i in real_chrom_sv_list]

            qry_diff_data = np.diff(real_chrom_qry_st)
            qry_dup_ind_data = sorted(filter(lambda t: t[1] == 0, enumerate(qry_diff_data)), key=lambda t: -t[0])

            for ind, _ in qry_dup_ind_data:
                sv = real_chrom_sv_list.pop(ind)
                err_ind_list.append(sv[0])

            real_chrom_qry_st = [sv[6] for sv in real_chrom_sv_list]

            qry_diff_data = np.diff(real_chrom_qry_st)
            qry_err_ind_data = list(filter(lambda t: t[1] < -err_baseline, enumerate(qry_diff_data)))

            if len(qry_err_ind_data) > 0:
                real_chrom_ref_qry_diff = [sv[6] - sv[4] for sv in real_chrom_sv_list]
                for ind, _ in qry_err_ind_data:
                    err_lft_group = get_group_index(real_chrom_ref_qry_diff, ind, reverse=True, baseline=err_baseline)
                    err_rht_group = get_group_index(real_chrom_ref_qry_diff, ind + 1, baseline=err_baseline)

                    if err_lft_group == 0 or len(real_chrom_sv_list) - 1 == err_rht_group:
                        if err_lft_group == 0:
                            err_ind_list.extend([sv[0] for sv in real_chrom_sv_list[err_lft_group:ind + 1]])
                        if len(real_chrom_sv_list) - 1 == err_rht_group:
                            err_ind_list.extend([sv[0] for sv in real_chrom_sv_list[ind + 1:err_rht_group + 1]])
                        
                    else:
                        if abs(real_chrom_ref_qry_diff[err_lft_group - 1] - real_chrom_ref_qry_diff[ind + 1]) < abs(real_chrom_ref_qry_diff[err_rht_group + 1] - real_chrom_ref_qry_diff[ind]):
                            err_ind_list.extend([sv[0] for sv in real_chrom_sv_list[err_lft_group:ind + 1]])
                        else:
                            err_ind_list.extend([sv[0] for sv in real_chrom_sv_list[ind + 1:err_rht_group + 1]])

    return err_ind_list


def correct_main(ref_loc, qry_loc_list, vcf_loc_list, species, sv_find_len=2000, repeat_find_len=50, workdir='data', sv_save='sv', min_sv_size=None,
                 num_cpus=1, overwrite_output=False, trust_query=False):
    check_file_exist([[ref_loc], qry_loc_list, vcf_loc_list], ['Reference', 'Query', 'Variant'])
    check_unique_basename(qry_loc_list)
    check_variant_caller(vcf_loc_list)
    
    if min_sv_size is None:
        min_sv_size = get_min_sv_size(vcf_loc_list[0])

        if min_sv_size is None:
            raise Exception('VCF file is not from GDBr, Please add --min_sv_size option')

    if len(qry_loc_list) != len(vcf_loc_list):
        raise Exception('The number of query and variant must be same')
    
    os.makedirs(workdir, exist_ok=True)
    safe_makedirs(sv_save, overwrite_output)

    # read .fasta file
    ref_seq = Fasta(ref_loc, build_index=False)

    # get 1Mbp chromosome
    ref_chr_list = list(map(lambda t: t[0], filter(lambda t: len(t[1]) > 1e6, ref_seq.records.items())))
    if not trust_query:
        p_map(partial(check_query_chrom, ref_chr_list=ref_chr_list), qry_loc_list, num_cpus=num_cpus)

    logprint(f'Task correct start : {len(vcf_loc_list)} query detected')
    for qry_ind, (qry_loc, vcf_loc) in enumerate(zip(qry_loc_list, vcf_loc_list)):
        qry_seq = Fasta(qry_loc, build_index=False)
        qry_chr_list = list(set(ref_chr_list) & set(qry_seq.records.keys()))
        
        qryworkdir = os.path.join(workdir, str(qry_ind) + '_gdbr')
        dbdir = os.path.join(qryworkdir, 'db_gdbr')
        
        os.makedirs(dbdir, exist_ok=True)

        # split query .fasta file and makeblastdb per chromosome
        p_map(partial(makeblastdb_from_location, seq_loc=qry_loc, dbdir=dbdir), qry_chr_list, num_cpus=num_cpus)

        vcf_tot_data = sorted(get_vcf_tot_data(vcf_loc, ref_chr_list), key=get_proper_correct_order)

        sv_list = p_map(partial(get_real_sv, sv_find_len=sv_find_len, repeat_find_len=repeat_find_len, min_sv_size=min_sv_size, ref_loc=ref_loc, qry_loc=qry_loc, dbdir=dbdir, qryworkdir=qryworkdir), vcf_tot_data, num_cpus=num_cpus)

        sv_list = sorted(sv_list, key=lambda t: t[0])

        rpmworkdir_list = prepare_repeatmasker_multiprocessing(ref_seq, qry_seq, sv_list, qryworkdir, num_cpus, repeat_find_len)
        rpm_ind_list_data = p_map(partial(get_reapeatmasker_index_list, species=species), rpmworkdir_list, num_cpus=num_cpus)
        for rpm_ind_list in rpm_ind_list_data:
            for rpm_ind in rpm_ind_list:
                sv_list[rpm_ind][2] = 'REPEAT:RPM'

        call_err_try = 3
        err_ind_list = [-1]
        while call_err_try > 0 and err_ind_list != []:
            call_err_try -= 1
            err_ind_list = get_correct_error_index_list(sv_list)
            for err_ind in err_ind_list:
                sv_list[err_ind][2] = 'SV_COR_ERR'

        qry_basename = os.path.basename(qry_loc)
        with open(os.path.join(sv_save, remove_gdbr_postfix(qry_basename)) + '.GDBr.correct.csv', 'w') as f:
            cf = csv.writer(f)
            cf.writerow(('ID', 'TYPE', 'REAL_TYPE', 'CHR', 'REF_START', 'REF_END', 'QRY_START', 'QRY_END'))
            cf.writerows(sv_list)

        logprint(f'{qry_ind + 1}/{len(qry_loc_list)} : {os.path.basename(vcf_loc)} correct complete')
