from gdbr.utilities import p_map, logprint, check_file_exist, check_unique_basename, check_variant_caller, get_min_sv_size, remove_gdbr_postfix, safe_makedirs, check_query_chrom
from functools import partial
from pyfaidx import Fasta

import subprocess
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


def correct_main(ref_loc, qry_loc_list, vcf_loc_list, species, sv_find_len=2000, repeat_find_len=50, workdir='data', sv_save='sv', min_sv_size=None,
                 num_cpus=1, pbar=True, telegram_token_loc='telegram.json', overwrite_output=False, trust_query=False):
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
    safe_makedirs(sv_save)

    # read .fasta file
    ref_seq = Fasta(ref_loc, build_index=False)

    # get 1Mbp chromosome
    ref_chr_list = list(map(lambda t: t[0], filter(lambda t: len(t[1]) > 1e6, ref_seq.records.items())))
    if not trust_query:
        p_map(partial(check_query_chrom, ref_chr_list=ref_chr_list), qry_loc_list, num_cpus=num_cpus, pbar=False)

    logprint(f'Task correct start : {len(vcf_loc_list)} query detected')
    output_data = []
    for qry_ind, (qry_loc, vcf_loc) in enumerate(zip(qry_loc_list, vcf_loc_list)):
        qry_seq = Fasta(qry_loc, build_index=False)
        qry_chr_list = list(set(ref_chr_list) & set(qry_seq.records.keys()))
        
        qryworkdir = os.path.join(workdir, str(qry_ind) + '_gdbr')
        dbdir = os.path.join(qryworkdir, 'db_gdbr')
        
        os.makedirs(dbdir, exist_ok=True)

        # split query .fasta file and makeblastdb per chromosome
        p_map(partial(makeblastdb_from_location, seq_loc=qry_loc, dbdir=dbdir), qry_chr_list, num_cpus=num_cpus, pbar=False)

        vcf_tot_data = sorted(get_vcf_tot_data(vcf_loc, ref_chr_list), key=get_proper_correct_order)

        sv_list = p_map(partial(get_real_sv, sv_find_len=sv_find_len, repeat_find_len=repeat_find_len, min_sv_size=min_sv_size, ref_loc=ref_loc, qry_loc=qry_loc, dbdir=dbdir, qryworkdir=qryworkdir), vcf_tot_data, 
                        num_cpus=num_cpus, pbar=pbar, telegram_token_loc=telegram_token_loc, desc=f'COR {qry_ind + 1}/{len(qry_loc_list)}')

        sv_list = sorted(sv_list, key=lambda t: t[0])

        rpm_temp_seq_loc = os.path.join(qryworkdir, str(qry_ind)) + '.fasta'
        with open(rpm_temp_seq_loc, 'w') as f:
            for sv in sv_list:
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
        
        subprocess.run(['RepeatMasker', '-spec', species, '-nopost', '-pa', str(num_cpus), rpm_temp_seq_loc], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=workdir)
        
        rpm_output_loc, output_postfix = (rpm_temp_seq_loc + '.cat', '.cat') if os.path.isfile(rpm_temp_seq_loc + '.cat') else (rpm_temp_seq_loc + '.cat.gz', '.cat.gz')
        subprocess.run(['ProcessRepeats', '-spec', species, '-gff', rpm_output_loc], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=workdir)
        
        if os.path.isfile(rpm_temp_seq_loc + '.out.gff'):
            with open(rpm_temp_seq_loc + '.out.gff', 'r') as f:
                tf = csv.reader(f, delimiter='\t')
                rpm_list = [l for l in tf]

            for rpm in rpm_list:
                if rpm[0][0] != '#':
                    sv_list[int(rpm[0].split('.')[0])][2] = 'REPEAT:RPM'
            
        for postfix in ['', output_postfix, '.out', '.out.gff', '.tbl']:
            if os.path.isfile(rpm_temp_seq_loc + postfix):
                os.remove(rpm_temp_seq_loc + postfix)

        qry_basename = os.path.basename(qry_loc)
        with open(os.path.join(sv_save, remove_gdbr_postfix(qry_basename)) + '.GDBr.correct.csv', 'w') as f:
            cf = csv.writer(f)
            cf.writerow(('ID', 'TYPE', 'REAL_TYPE', 'CHR', 'REF_START', 'REF_END', 'QRY_START', 'QRY_END'))
            cf.writerows(sv_list)

        logprint(f'{qry_ind + 1}/{len(qry_loc_list)} : {os.path.basename(vcf_loc)} correct complete')
