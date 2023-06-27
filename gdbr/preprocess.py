from functools import partial
from pyfaidx import Fasta
from gdbr.my_tqdm import p_map

import subprocess
import shutil
import vcfpy
import csv
import os

blastn_fmt = ['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'bitscore']

ref_loc = ''
qry_loc = ''
vcf_loc = ''
dbdir = ''
qryworkdir = ''

def get_blast_result(start_temp_seq, end_temp_seq, sv_find_len, sv_id, chrom):
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

def get_blast_single_result(ref_temp_seq, qry_temp_seq, sv_find_len, sv_id, filter_func):
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

def get_real_sv(record, sv_find_len=2000):
    tid = record.ID[0].split('.')[1:-1][0]
    
    if tid not in {'DEL', 'INS'}:
        return 'UNSUP_ID',

    ref_seq = Fasta(ref_loc, build_index=False)
    qry_seq = Fasta(qry_loc, build_index=False)

    ref_start, ref_end, ref_len, sv_id, chrom = record.POS + 1, record.INFO['END'], abs(record.INFO['SVLEN']), record.ID[0], record.CHROM
    
    start_temp_seq = ref_seq[chrom][ref_start - sv_find_len - 1:ref_start - 1]
    end_temp_seq = ref_seq[chrom][ref_end:ref_end + sv_find_len]

    start_filter_result, end_filter_result = get_blast_result(start_temp_seq, end_temp_seq, sv_find_len, sv_id, chrom)
    
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

        start_filter_result, end_filter_result = get_blast_result(start_temp_seq, end_temp_seq, sv_find_len, sv_id, chrom)
        
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

        start_sub_filter_result = get_blast_single_result(start_ref_temp_seq, start_qry_temp_seq, sv_find_len, sv_id, lambda t: t[5] >= sv_find_len and t[7] >= sv_find_len)

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

        end_sub_filter_result = get_blast_single_result(end_ref_temp_seq, end_qry_temp_seq, sv_find_len, sv_id, lambda t: t[4] <= end_ref_check_len + 1 and t[6] <= end_qry_check_len + 1)

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

def preprocess_main(ref_loc_, qry_loc_list, vcf_loc_list, workdir='data', force=False, file=True, **pbar_arg):
    # working directory setting
    global ref_loc, qry_loc, dbdir, qryworkdir

    ref_loc = ref_loc_
    force = True

    if len(qry_loc_list) != len(vcf_loc_list):
        raise Exception('The number of query and variant must be same')

    if not os.path.isdir(workdir):
        try:
            os.mkdir(workdir)
        except OSError as error:
            raise Exception(str(error))

    # read .fasta file
    ref_seq = Fasta(ref_loc, build_index=True)

    # get 1Mbp chromosome
    ref_chr_list = list(map(lambda t: t[0], filter(lambda t: len(t[1]) > 1e6, ref_seq.records.items())))

    output_data = []
    for qry_ind, (qry_loc, vcf_loc) in enumerate(zip(qry_loc_list, vcf_loc_list)):
        qry_seq = Fasta(qry_loc, build_index=True)

        # check all ref chromosome in all qry
        if set(ref_chr_list) > set(qry_seq.records.keys()):
            raise Exception('Chromone Error')
        
        qryworkdir = os.path.join(workdir, str(qry_ind))
        dbdir = os.path.join(qryworkdir, 'db')

        os.mkdir(qryworkdir)
        os.mkdir(dbdir)

        # split query .fasta file and makeblastdb per chromosome
        for chr_name in ref_chr_list:
            with open(os.path.join(qryworkdir, chr_name + '.fasta'), 'w') as f:
                tseq = qry_seq[chr_name]
                f.write('>' + chr_name + '\n')
                f.write(str(tseq))
            
            subprocess.run(['makeblastdb', '-in', os.path.join(qryworkdir, chr_name + '.fasta'), '-input_type', 'fasta', '-dbtype', 'nucl', '-out', os.path.join(dbdir, chr_name)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            os.remove(os.path.join(qryworkdir, chr_name + '.fasta'))

        vcf_tot_data = []
        record = vcfpy.Reader.from_path(vcf_loc)
        for r in record:
            if r.CHROM in ref_chr_list:
                vcf_tot_data.append(r)
        record.close()

        sv_find_len = 2000
        sv_list = p_map(partial(get_real_sv, sv_find_len=sv_find_len), vcf_tot_data, **pbar_arg)
        
        if file:
            with open(f'sv_call_{qry_ind}.csv', 'w') as f:
                cf = csv.writer(f)
                cf.writerow(('ID', 'SV_TYPE', 'CHROM', 'REF_START', 'REF_END', 'QRY_START', 'QRY_END'))
                cf.writerows([[i] + list(v) for i, v in enumerate(sv_list)])
        else:
            output_data.append([[i] + list(v) for i, v in enumerate(sv_list)])

    if not file:
        return output_data