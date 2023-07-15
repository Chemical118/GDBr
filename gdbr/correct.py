from gdbr.utilities import p_map, logprint
from functools import partial
from pyfaidx import Fasta

import subprocess
import vcfpy
import csv
import os

#                [0]       [1]        [2]        [3]        [4]      [5]      [6]      [7]       [8]
blastn_fmt = ['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'bitscore']


def get_blast_result(start_temp_seq, end_temp_seq, sv_find_len, sv_id, chrom, dbdir, qryworkdir):
    start_temp_seq_loc = os.path.join(qryworkdir, f'ref_start_{sv_id}.fasta')
    end_temp_seq_loc = os.path.join(qryworkdir, f'ref_end_{sv_id}.fasta')

    with open(start_temp_seq_loc, 'w') as f:
        f.write('>' + start_temp_seq.fancy_name + '\n')
        f.write(str(start_temp_seq))

    with open(end_temp_seq_loc, 'w') as f:
        f.write('>' + end_temp_seq.fancy_name + '\n')
        f.write(str(end_temp_seq))

    start_result = subprocess.run(['blastn', '-db', os.path.join(dbdir, chrom), '-query', start_temp_seq_loc, '-strand', 'plus', '-outfmt', '10 ' + ' '.join(blastn_fmt)], capture_output=True, text=True)
    end_result = subprocess.run(['blastn', '-db', os.path.join(dbdir, chrom), '-query', end_temp_seq_loc, '-strand', 'plus', '-outfmt', '10 ' + ' '.join(blastn_fmt)], capture_output=True, text=True)

    start_result_list =  [] if start_result.stdout == '' else map(lambda t: list(map(eval, t.split(','))), start_result.stdout[:-1].split('\n'))
    end_result_list = [] if end_result.stdout == '' else map(lambda t: list(map(eval, t.split(','))), end_result.stdout[:-1].split('\n'))

    os.remove(start_temp_seq_loc)
    os.remove(end_temp_seq_loc)

    start_filter_result = list(filter(lambda t: t[1] > sv_find_len * 0.90, start_result_list))
    end_filter_result = list(filter(lambda t: t[1] > sv_find_len * 0.90, end_result_list))

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


def get_real_sv(record, ref_loc, qry_loc, dbdir, qryworkdir, sv_find_len=2000):
    tid = record.ID[0].split('.')[1:-1][0]
    
    if tid not in {'DEL', 'INS'}:
        return 'UNSUP_ID',

    ref_seq = Fasta(ref_loc, build_index=False)
    qry_seq = Fasta(qry_loc, build_index=False)

    ref_start, ref_end, ref_len, sv_id, chrom = record.POS + 1, record.INFO['END'], abs(record.INFO['SVLEN']), record.ID[0], record.CHROM
    
    start_temp_seq = ref_seq[chrom][ref_start - sv_find_len - 1:ref_start - 1]
    end_temp_seq = ref_seq[chrom][ref_end:ref_end + sv_find_len]

    start_filter_result, end_filter_result = get_blast_result(start_temp_seq, end_temp_seq, sv_find_len, sv_id, chrom, dbdir, qryworkdir)
    
    if len(start_filter_result) != 1 or len(end_filter_result) != 1:
        return f'FND_IDX:({len(start_filter_result)}, {len(end_filter_result)})',

    ref_start, ref_end, ref_len, qry_start, qry_end, qry_len = get_corrected_location(ref_start, ref_end, sv_find_len, start_filter_result, end_filter_result)

    if qry_len > 1e6 or ref_len > 1e6:
        sv_find_len *= 2

        start_temp_seq = ref_seq[chrom][ref_start - sv_find_len - 1:ref_start - 1]
        end_temp_seq = ref_seq[chrom][ref_end:ref_end + sv_find_len]

        start_filter_result, end_filter_result = get_blast_result(start_temp_seq, end_temp_seq, sv_find_len, sv_id, chrom, dbdir, qryworkdir)
    
        if len(start_filter_result) != 1 or len(end_filter_result) != 1:
            return f'FND_IDX:({len(start_filter_result)}, {len(end_filter_result)})',

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
            return f'FND_IDX:({len(start_filter_result)}, {len(end_filter_result)})',

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
        
    if qry_len < 0:
        sv_name = 'ERR_POS'

    elif ref_len == 0 and qry_len == 0:
        sv_name = 'FP_SV'

    elif ref_len > 0 and qry_len > 0:
        sv_name = 'SUB'

    elif ref_len == 0:
        sv_name = 'INS'
    
    else:
        sv_name = 'DEL'

    return sv_name, chrom, ref_start, ref_end, qry_start, qry_end    


def makeblastdb_from_location(chr_name, seq_loc, dbdir):
    seq = Fasta(seq_loc, build_index=False)
    with open(os.path.join(dbdir, chr_name + '.fasta'), 'w') as f:
        tseq = seq[chr_name]
        f.write('>' + chr_name + '\n')
        f.write(str(tseq))
    
    subprocess.run(['makeblastdb', '-in', os.path.join(dbdir, chr_name + '.fasta'), '-input_type', 'fasta', '-dbtype', 'nucl', '-out', os.path.join(dbdir, chr_name)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.remove(os.path.join(dbdir, chr_name + '.fasta'))


def correct_main(ref_loc, qry_loc_list, vcf_loc_list, sv_find_len=2000, workdir='data', sv_save='sv', file=True, num_cpus=1, pbar=True, telegram_token_loc='telegram.json'):
    qry_basename_list = list(map(os.path.basename, qry_loc_list))

    if len(qry_basename_list) != len(set(qry_basename_list)):
        raise Exception('Query basename must be diffrent')
    
    if len(qry_loc_list) != len(vcf_loc_list):
        raise Exception('The number of query and variant must be same')
    
    os.makedirs(workdir, exist_ok=True)
    os.makedirs(sv_save, exist_ok=True)

    # read .fasta file
    ref_seq = Fasta(ref_loc, build_index=True)

    # get 1Mbp chromosome
    ref_chr_list = list(map(lambda t: t[0], filter(lambda t: len(t[1]) > 1e6, ref_seq.records.items())))

    logprint(f'Task start : {len(vcf_loc_list)} query detected')
    output_data = []
    for qry_ind, (qry_loc, vcf_loc) in enumerate(zip(qry_loc_list, vcf_loc_list)):
        qry_seq = Fasta(qry_loc, build_index=True)

        # If there no same chromosone in query and reference
        qry_chr_list = list(set(ref_chr_list) & set(qry_seq.records.keys()))
        if qry_chr_list == []:
            raise Exception('Chromosone name must same')
        
        qryworkdir = os.path.join(workdir, str(qry_ind))
        dbdir = os.path.join(qryworkdir, 'db')
        
        os.makedirs(dbdir, exist_ok=True)

        # split query .fasta file and makeblastdb per chromosome
        p_map(partial(makeblastdb_from_location, seq_loc=qry_loc, dbdir=dbdir), qry_chr_list, num_cpus=num_cpus, pbar=False)

        with vcfpy.Reader.from_path(vcf_loc) as tot_record:
            vcf_tot_data = [record for record in tot_record if record.CHROM in ref_chr_list]
        
        sv_list = p_map(partial(get_real_sv, sv_find_len=sv_find_len, ref_loc=ref_loc, qry_loc=qry_loc, dbdir=dbdir, qryworkdir=qryworkdir), vcf_tot_data, 
                        num_cpus=num_cpus, pbar=pbar, telegram_token_loc=telegram_token_loc, desc=f'COR {qry_ind + 1}/{len(qry_loc_list)}')
        
        if file:
            qry_basename = os.path.basename(qry_loc)
            with open(os.path.join(sv_save, qry_basename) + '.COR.csv', 'w') as f:
                cf = csv.writer(f)
                cf.writerow(('ID', 'SV_TYPE', 'CHR', 'REF_START', 'REF_END', 'QRY_START', 'QRY_END'))
                cf.writerows([[i] + list(v) for i, v in enumerate(sv_list)])
        else:
            output_data.append([[i] + list(v) for i, v in enumerate(sv_list)])

        logprint(f'{qry_ind + 1}/{len(qry_loc_list)} : {os.path.basename(qry_loc)} correction complete')

    if not file:
        return output_data