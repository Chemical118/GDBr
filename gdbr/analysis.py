from Bio.Blast import NCBIXML
from functools import partial
from pyfaidx import Fasta
from gdbr.utilities import p_map, logprint
from gdbr.correct import makeblastdb_from_location

import multiprocessing as mp

import subprocess
import itertools
import csv
import os

#                [0]       [1]        [2]        [3]        [4]      [5]      [6]      [7]       [8]        [9]
blastn_fmt = ['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'bitscore', 'sseqid']


def get_blastn_ref_chr_list(dbdir):
    return [i.split('.')[0] for i in os.listdir(dbdir) if i.split('.')[1] == 'ndb']


def blast_output_to_list(blast_fmt_output):
    blast_list = blast_fmt_output.split(',')
    return [float(blast_list[0])] + list(map(int, blast_list[1:8])) + [float(blast_list[8]), blast_list[9]]


def get_sv_list(sv_loc):
    with open(sv_loc, 'r') as f:
        cf = csv.reader(f)
        sv_list = [[int(i[0]), i[1], i[2], int(i[3]), int(i[4]), int(i[5]), int(i[6])] if len(i) == 7 else [int(i[0]), i[1]] for i in [l for l in cf][1:]]
    return sv_list


def get_one_way_homology(ref_part_seq, qry_part_seq, ref_hom_find_len, qry_hom_find_len, sv_id, qryworkdir):
    blast_filter_result = []
    ref_hom = ref_bitscore = 0
    temp_ref_hom_find_len = 3

    ref_part_seq_loc = os.path.join(qryworkdir, f'ref_{sv_id}.fasta')
    temp_qry_part_seq_loc = os.path.join(qryworkdir, f'qry_{sv_id}.fasta')

    with open(ref_part_seq_loc, 'w') as f:
        f.write('>' + ref_part_seq.fancy_name + '\n')
        f.write(str(ref_part_seq))

    first_flag = True
    while ref_hom == temp_ref_hom_find_len or first_flag:
        first_flag = False

        temp_ref_hom_find_len *= 2 if ref_hom == temp_ref_hom_find_len else 1
        temp_qry_part_seq = qry_part_seq[:qry_hom_find_len + temp_ref_hom_find_len]

        with open(temp_qry_part_seq_loc, 'w') as f:
            f.write('>' + temp_qry_part_seq.fancy_name + '\n')
            f.write(str(temp_qry_part_seq))
        
        blast_result = subprocess.run(['blastn', '-subject', ref_part_seq_loc, '-query', temp_qry_part_seq_loc, '-strand', 'plus', '-outfmt', '10 ' + ' '.join(blastn_fmt)], capture_output=True, text=True)
        blast_result_list =  [] if blast_result.stdout == '' else map(blast_output_to_list, blast_result.stdout[:-1].split('\n'))
        blast_filter_result = sorted(filter(lambda t: t[1] > ref_hom_find_len * 0.90 and t[6] < 0.2 * qry_hom_find_len and t[7] > qry_hom_find_len and t[4] < 0.2 * ref_hom_find_len and t[5] > ref_hom_find_len and t[8] > ref_bitscore, blast_result_list), key=lambda t: t[8])

        if len(blast_filter_result) == 0:
            break
        
        blast_filter_result = blast_filter_result[0]
        ref_hom = blast_filter_result[7] - ref_hom_find_len
        ref_bitscore = blast_filter_result[8]

    os.remove(ref_part_seq_loc)
    os.remove(temp_qry_part_seq_loc)
    
    ref_hom_end = ref_hom_find_len - 1 if blast_filter_result == [] else blast_filter_result[7] - 1 
    qry_hom_end = qry_hom_find_len - 1 if blast_filter_result == [] else blast_filter_result[5] - 1

    ref_part_seq_len = len(ref_part_seq)
    qry_part_seq_len = len(qry_part_seq)

    # index test
    while ref_hom_end < ref_part_seq_len - 1 and qry_hom_end < qry_part_seq_len - 1 and str(ref_part_seq[ref_hom_end + 1]).upper() == str(qry_part_seq[qry_hom_end + 1]).upper():
        ref_hom += 1
        ref_hom_end += 1
        qry_hom_end += 1

    return ref_hom


def get_one_way_templated_insertion(ref_seq, qry_seq, qryworkdir, temp_indel_find_len, near_gap_find_len, sv_id, gap_baseline, chrom, ref_start, ref_end, ref_len, qry_start, qry_end, qry_len):
    ref_subject_temp_seq = ref_seq[chrom][ref_start - 1 - temp_indel_find_len:ref_end + temp_indel_find_len]

    ref_left_hom = ref_right_hom = ref_bitscore = 0
    ref_left_temp_hom_find_len = ref_right_temp_hom_find_len = 3

    ref_subject_temp_seq_loc = os.path.join(qryworkdir, f'ref_sub_{sv_id}.fasta')
    qry_query_temp_seq_loc = os.path.join(qryworkdir, f'qry_qry_{sv_id}.fasta')
    ref_blast_result_loc = os.path.join(qryworkdir, f'ref_blast_{sv_id}.xml')

    with open(ref_subject_temp_seq_loc, 'w') as f:
        f.write('>' + ref_subject_temp_seq.fancy_name + '\n')
        f.write(str(ref_subject_temp_seq))

    first_flag = True
    while ref_left_hom == ref_left_temp_hom_find_len or ref_right_hom == ref_right_temp_hom_find_len or first_flag:
        first_flag = False

        # each side homology length
        ref_left_temp_hom_find_len *= 2 if ref_left_hom == ref_left_temp_hom_find_len else 1
        ref_right_temp_hom_find_len *= 2 if ref_right_hom == ref_right_temp_hom_find_len else 1
        
        qry_query_temp_seq = qry_seq[chrom][qry_start - 1 - ref_left_temp_hom_find_len:qry_end + ref_right_temp_hom_find_len]
        with open(qry_query_temp_seq_loc, 'w') as f:
            f.write('>' + qry_query_temp_seq.fancy_name + '\n')
            f.write(str(qry_query_temp_seq))
        
        with open(ref_blast_result_loc, 'w') as f:
            subprocess.run(['blastn', '-task', 'blastn-short', '-subject', ref_subject_temp_seq_loc, '-query', qry_query_temp_seq_loc, '-strand', 'plus', '-outfmt', '5'], stdout=f, stderr=subprocess.DEVNULL)
        
        with open(ref_blast_result_loc, 'r') as f:
            ref_blast_result = NCBIXML.read(f).alignments

        os.remove(ref_blast_result_loc)
        
        ref_blast_result_list = [] if ref_blast_result == [] else ref_blast_result[0].hsps
        # must align at both side
        ref_blast_filter_result = sorted(filter(lambda t: t.query_start < ref_left_temp_hom_find_len + 1 and ref_left_temp_hom_find_len + qry_len < t.query_end and t.gaps < gap_baseline and t.bits > ref_bitscore, ref_blast_result_list), key=lambda t: t.bits, reverse=True)

        if len(ref_blast_filter_result) == 0:
            break
        
        for blast_result in ref_blast_filter_result:
            # find real location without gap
            sub_real_idx_list = [i for i, v in enumerate(blast_result.sbjct) if v != '-']

            sub_sv_loc_st = temp_indel_find_len + 1 - blast_result.sbjct_start
            sub_sv_loc_nd = temp_indel_find_len + ref_len - blast_result.sbjct_start
            
            # out of blast location
            sub_sv_real_loc_st = sub_real_idx_list[sub_sv_loc_st] if 0 <= sub_sv_loc_st < len(sub_real_idx_list) else sub_sv_loc_st
            sub_sv_real_loc_nd = sub_real_idx_list[sub_sv_loc_nd] if 0 <= sub_sv_loc_nd < len(sub_real_idx_list) else sub_sv_loc_nd

            qry_real_idx_list = [i for i, v in enumerate(blast_result.query) if v != '-']

            qry_sv_real_loc_st = qry_real_idx_list[ref_left_temp_hom_find_len + 1 - blast_result.query_start]
            qry_sv_real_loc_nd = qry_real_idx_list[ref_left_temp_hom_find_len + qry_len - blast_result.query_start]

            # select smaller zone also, exclude subject is out of match zone
            if sub_sv_real_loc_st < 0 or sub_sv_real_loc_nd >= len(blast_result.match) or sub_sv_real_loc_nd - sub_sv_real_loc_st >= qry_sv_real_loc_nd - qry_sv_real_loc_st:
                short_sv_real_loc_st = qry_sv_real_loc_st
                short_sv_real_loc_nd = qry_sv_real_loc_nd

            else:
                short_sv_real_loc_st = sub_sv_real_loc_st
                short_sv_real_loc_nd = sub_sv_real_loc_nd

            # selected zone must all match
            match_sv_seq = blast_result.match[short_sv_real_loc_st:short_sv_real_loc_nd + 1]
            if ' ' in match_sv_seq:
                continue
            
            # smaller zone and nearby near_gap_find_len must have no gap
            # COMBO gap search
            """
            TAAAA
            T-AAA <=> TA-AA TAA-A TAAA-

            at this situation

            i = 4 (at for loop)
            combo_st = 1
            combo_tar = 'A'
            combo_isgap = True
            """

            combo_st = 0
            combo_tar = ''
            combo_isgap = False

            for i in range(len(blast_result.query)):
                if combo_tar == '':
                    combo_isgap = min(blast_result.query[i], blast_result.sbjct[i]) == '-'
                    if blast_result.query[i] == blast_result.sbjct[i] or combo_isgap:
                        combo_st = i
                        combo_tar = max(blast_result.query[i], blast_result.sbjct[i])

                if combo_tar != '':
                    if i + 1 > len(blast_result.query) - 1 or (min(blast_result.query[i + 1], blast_result.sbjct[i + 1]) == '-' and max(blast_result.query[i + 1], blast_result.sbjct[i + 1]) != combo_tar)\
                            or (min(blast_result.query[i + 1], blast_result.sbjct[i + 1]) != '-' and (blast_result.query[i + 1] != blast_result.sbjct[i + 1] or blast_result.query[i + 1] != combo_tar)):
                        
                        if combo_isgap and (short_sv_real_loc_st - near_gap_find_len <= combo_st <= short_sv_real_loc_nd + near_gap_find_len or short_sv_real_loc_st - near_gap_find_len <= i <= short_sv_real_loc_nd + near_gap_find_len):
                            break

                        combo_tar = ''
                        combo_isgap = False

                    else:
                        combo_isgap |= min(blast_result.query[i + 1], blast_result.sbjct[i + 1]) == '-'

            if combo_isgap:
                continue

            ref_left_hom = ref_left_temp_hom_find_len - blast_result.query_start + 1
            ref_right_hom = blast_result.query_end - qry_len - ref_left_temp_hom_find_len
            ref_bitscore = blast_result.bits

            ref_sub_st = ref_start - 1 - temp_indel_find_len + blast_result.sbjct_start
            ref_sub_nd = ref_start - 1 - temp_indel_find_len + blast_result.sbjct_end
            ref_qry_st = qry_start - 1 - ref_left_temp_hom_find_len + blast_result.query_start
            ref_qry_nd = qry_start - 1 - ref_left_temp_hom_find_len + blast_result.query_end
            break
    
    # index search
    if ref_left_hom > 0:
        # left side
        while str(ref_seq[chrom][ref_sub_st - 2]).upper() == str(qry_seq[chrom][ref_qry_st - 2]).upper():
            ref_sub_st -= 1
            ref_qry_st -= 1
            ref_left_hom += 1

        # right side
        while str(ref_seq[chrom][ref_sub_nd]).upper() == str(qry_seq[chrom][ref_qry_nd]).upper():
            ref_sub_nd += 1
            ref_qry_nd += 1
            ref_right_hom += 1

    os.remove(ref_subject_temp_seq_loc)
    os.remove(qry_query_temp_seq_loc)

    return ref_left_hom, ref_right_hom


def is_away_from_locus(ref_start, ref_end, tar_start, tar_end, near_seq_kb_baseline):
    dis1 = ref_start - tar_end
    dis2 = ref_end - tar_start

    # SV and target is overlap
    if dis1 * dis2 <= 0:
        return False
    else:
        dis = min(abs(dis1), abs(dis2))
        return dis > near_seq_kb_baseline * 1e3


def _multiprocessing_find_insertion_chrom(chr_name, near_seq_kb_baseline, refdbdir, len_tar_seq, tar_seq_loc, chrom, ref_start, ref_end):
    with lock:
        if global_num_blast_result.value > 1:
            return []
    
    blast_result = subprocess.run(['blastn', '-db', os.path.join(refdbdir, chr_name), '-query', tar_seq_loc, '-strand', 'plus', '-outfmt', '10 ' + ' '.join(blastn_fmt)], capture_output=True, text=True)
    blast_result_list =  [] if blast_result.stdout == '' else map(blast_output_to_list, blast_result.stdout[:-1].split('\n'))
    blast_filter_result = list(filter(lambda t: t[1] > len_tar_seq * 0.95, blast_result_list)) if chr_name != chrom else \
                            list(filter(lambda t: t[1] > len_tar_seq * 0.95 and is_away_from_locus(ref_start, ref_end, t[6], t[7], near_seq_kb_baseline), blast_result_list))
    
    with lock:
        global_num_blast_result.value += len(blast_filter_result)
        if global_num_blast_result.value > 1:
            return []

    if len(blast_filter_result) == 1:
        return blast_filter_result[0][9], blast_filter_result[0][6], blast_filter_result[0][7]
    else:
        return []


def _multiprocessing_init(_num_blast_result, _lock):
    global global_num_blast_result, lock
    global_num_blast_result = _num_blast_result
    lock = _lock


def get_non_homology_insertion_loaction(tar_seq, near_seq_kb_baseline, refdbdir, sv_id, chrom, qryworkdir, ref_start, ref_end, ref_chr_list, num_cpus):
    len_tar_seq = len(tar_seq)
    tar_seq_loc = os.path.join(qryworkdir, f'diff_locus_{sv_id}.fasta')
    with open(tar_seq_loc, 'w') as f:
        f.write('>' + tar_seq.fancy_name + '\n')
        f.write(str(tar_seq))

    lock = mp.Lock()
    num_blast_result = mp.Value('i', 0)
    ins_chrom = False
    ins_start = ins_end = -1

    with mp.Pool(initializer=_multiprocessing_init, initargs=(num_blast_result, lock), processes=num_cpus) as p:
        blast_total_result = p.map(partial(_multiprocessing_find_insertion_chrom, near_seq_kb_baseline=near_seq_kb_baseline, refdbdir=refdbdir,
                                           len_tar_seq=len_tar_seq, tar_seq_loc=tar_seq_loc, chrom=chrom, ref_start=ref_start, ref_end=ref_end), ref_chr_list)
    
    os.remove(tar_seq_loc)

    if num_blast_result.value != 1:
        ins_chrom = num_blast_result.value
    else:
        ins_chrom, ins_start, ins_end = tuple(itertools.chain.from_iterable(blast_total_result))

    return ins_chrom, ins_start, ins_end


def get_homology_hard(sv_data, ref_loc, qry_loc, refdbdir, qryworkdir, ref_chr_list, num_cpus, hom_find_len=2000, near_seq_kb_baseline=100.0, diff_locus_hom_baseline=3):
    sv_id, tid, chrom, ref_start, ref_end, qry_start, qry_end = sv_data

    dsb_repair_type = 'NO_DSBR'
    dsb_repair_type_code = 0

    ref_seq = Fasta(ref_loc, build_index=False)
    qry_seq = Fasta(qry_loc, build_index=False)

    ref_len = ref_end - ref_start + 1
    qry_len = qry_end - qry_start + 1

    ref_ins_chrom, ref_ins_start, ref_ins_end = get_non_homology_insertion_loaction(ref_seq[chrom][ref_start - 1:ref_end], near_seq_kb_baseline,
                                                                                    refdbdir, sv_id, chrom, qryworkdir, ref_start, ref_end, ref_chr_list, num_cpus)
        
    qry_ins_chrom, qry_ins_start, qry_ins_end = get_non_homology_insertion_loaction(qry_seq[chrom][qry_start - 1:qry_end], near_seq_kb_baseline,
                                                                                    refdbdir, sv_id, chrom, qryworkdir, ref_start, ref_end, ref_chr_list, num_cpus)
    is_find_ins_ref = isinstance(ref_ins_chrom, str)
    is_find_ins_qry = isinstance(qry_ins_chrom, str)

    if is_find_ins_ref or is_find_ins_qry:
        ref_ins_lef_hom = ref_ins_rht_hom = qry_ins_lef_hom = qry_ins_rht_hom = -1

        if is_find_ins_ref:
            ref_ins_len = ref_ins_end - ref_ins_start + 1                    
            ref_ins_lef_hom = get_one_way_homology(ref_seq[chrom][max(0, ref_start - 1 - hom_find_len):ref_end].reverse,
                                                   ref_seq[ref_ins_chrom][ref_ins_start - 1 - hom_find_len:ref_ins_end].reverse,
                                                   ref_len, ref_ins_len, sv_id, qryworkdir)
            
            ref_ins_rht_hom = get_one_way_homology(ref_seq[chrom][ref_start - 1:ref_end + hom_find_len],
                                                   ref_seq[ref_ins_chrom][ref_ins_start - 1:ref_ins_end + hom_find_len],
                                                   ref_len, ref_ins_len, sv_id, qryworkdir)

        if is_find_ins_qry:
            qry_ins_len = qry_ins_end - qry_ins_start + 1 
            qry_ins_lef_hom = get_one_way_homology(qry_seq[chrom][max(0, qry_start - 1 - hom_find_len):qry_end].reverse,
                                                   ref_seq[ref_ins_chrom][qry_ins_start - 1:qry_ins_end + hom_find_len].reverse,
                                                   qry_len, qry_ins_len, sv_id, qryworkdir)
            
            qry_ins_rht_hom = get_one_way_homology(qry_seq[chrom][qry_start - 1:qry_end + hom_find_len],
                                                   ref_seq[ref_ins_chrom][qry_ins_start - 1:qry_ins_end + hom_find_len],
                                                   qry_len, qry_ins_len, sv_id, qryworkdir)
        
        dsb_repair_type_code = 3
        ins_chrom, ins_start, ins_end, lef_hom, rht_hom = (ref_ins_chrom, ref_ins_start, ref_ins_end, ref_ins_lef_hom, ref_ins_rht_hom) \
                                                            if ref_ins_lef_hom + ref_ins_rht_hom + (ref_len > qry_len) - 0.5 > qry_ins_lef_hom + qry_ins_rht_hom else \
                                                            (qry_ins_chrom, qry_ins_start, qry_ins_end, qry_ins_lef_hom, qry_ins_rht_hom)
        if lef_hom + rht_hom < diff_locus_hom_baseline:
            dsb_repair_type = 'SUB_UNIQUE'
        else:
            dsb_repair_type = 'DIFF_LOCUS_DSBR'

    else:
        ins_chrom = ref_ins_chrom if ref_len > qry_len else qry_ins_chrom
        
        if ins_chrom:
            dsb_repair_type = 'SUB_REPEAT'
        else:
            dsb_repair_type = 'SUB_UNKNOWN'


    if dsb_repair_type_code == 3:
        return sv_id, tid, dsb_repair_type, lef_hom, rht_hom, ins_chrom, ins_start, ins_end

    else:
        return sv_id, tid, dsb_repair_type


def get_homology(sv_data, ref_loc, qry_loc, qryworkdir, hom_find_len=2000, temp_indel_find_len=100, near_gap_find_len=5, user_gap_baseline=3):
    if sv_data[1] not in {'DEL', 'INS', 'SUB'}:
        return sv_data[0], 'UNSUP_ID'

    sv_id, tid, chrom, ref_start, ref_end, qry_start, qry_end = sv_data

    dsb_repair_type = 'NO_DSBR'
    dsb_repair_type_code = 0

    ref_seq = Fasta(ref_loc, build_index=False)
    qry_seq = Fasta(qry_loc, build_index=False)

    ref_len = ref_end - ref_start + 1
    qry_len = qry_end - qry_start + 1

    if tid == 'SUB':
        sub_gap_baseline = min(user_gap_baseline, ref_len, qry_len)
        ref_lef_hom, ref_rht_hom = get_one_way_templated_insertion(ref_seq, qry_seq, qryworkdir,
                                                                   temp_indel_find_len, near_gap_find_len, sv_id, sub_gap_baseline, chrom, 
                                                                   ref_start, ref_end, ref_len,
                                                                   qry_start, qry_end, qry_len)
        
        qry_lef_hom, qry_rht_hom = get_one_way_templated_insertion(qry_seq, ref_seq, qryworkdir,
                                                                   temp_indel_find_len, near_gap_find_len, sv_id, sub_gap_baseline, chrom,
                                                                   qry_start, qry_end, qry_len,
                                                                   ref_start, ref_end, ref_len)
                                                                   
        
        if ref_lef_hom > 0 or qry_lef_hom > 0:
            dsb_repair_type = 'TEMP_INS'
            dsb_repair_type_code = 2
            lef_hom, rht_hom = (ref_lef_hom, ref_rht_hom) if ref_lef_hom + ref_rht_hom >= qry_lef_hom + qry_rht_hom else (qry_lef_hom, qry_rht_hom)

        else:
            dsb_repair_type = 'SUB_HARD'

    else:
        lef_hom = rht_hom = 0
        if tid == 'DEL':
            lef_hom = get_one_way_homology(ref_seq[chrom][max(0, ref_start - 1 - hom_find_len):ref_end],
                                           qry_seq[chrom][max(0, qry_start - 1 - hom_find_len):qry_end + ref_len],
                                           hom_find_len, hom_find_len, sv_id, qryworkdir)
            
            rht_hom = get_one_way_homology(ref_seq[chrom][ref_start - 1:ref_end + hom_find_len].reverse,
                                           qry_seq[chrom][max(0, qry_start - 1 - ref_len):qry_end + hom_find_len].reverse,
                                           hom_find_len, hom_find_len, sv_id, qryworkdir)

        else:
            lef_hom = get_one_way_homology(qry_seq[chrom][max(0, qry_start - 1 - hom_find_len):qry_end],
                                           ref_seq[chrom][max(0, ref_start - 1 - hom_find_len):ref_end + qry_len],
                                           hom_find_len, hom_find_len, sv_id, qryworkdir)
            
            rht_hom = get_one_way_homology(qry_seq[chrom][qry_start - 1:qry_end + hom_find_len].reverse,
                                           ref_seq[chrom][max(0, ref_start - 1 - qry_len):ref_end + hom_find_len].reverse,
                                           hom_find_len, hom_find_len, sv_id, qryworkdir)

        hom = lef_hom + rht_hom

        if hom > (ref_len if tid == 'DEL' else qry_len):
            dsb_repair_type = 'HOM_LEN'
        elif hom > 0:
            dsb_repair_type_code = 1
            dsb_repair_type = 'HOM'

    if dsb_repair_type_code == 1:
        return sv_id, tid, dsb_repair_type, hom

    elif dsb_repair_type_code == 2:
        return sv_id, tid, dsb_repair_type, lef_hom, rht_hom

    else:
        return sv_id, tid, dsb_repair_type


def analysis_main(ref_loc, qry_loc_list, sv_loc_list, hom_find_len=2000, temp_indel_find_len=100, near_gap_find_len=5, user_gap_baseline=3, near_seq_kb_baseline=100.0, diff_locus_hom_baseline=3, workdir='data', dsbr_save='dsbr', file=True, num_cpus=1, pbar=True, telegram_token_loc='telegram.json'):
    qry_basename_list = list(map(os.path.basename, qry_loc_list))

    if len(qry_basename_list) != len(set(qry_basename_list)):
        raise Exception('Query basename must be diffrent')
    
    # read .fasta file
    ref_seq = Fasta(ref_loc, build_index=True)
    
    if len(qry_loc_list) != len(sv_loc_list):
        raise Exception('The number of query and variant must be same')

    # get 1Mbp chromosome
    ref_chr_list = list(map(lambda t: t[0], filter(lambda t: len(t[1]) > 1e6, ref_seq.records.items())))

    refdbdir = os.path.join(workdir, 'db')
    os.makedirs(refdbdir, exist_ok=True)
    os.makedirs(dsbr_save, exist_ok=True)
    
    # split reference .fasta file and makeblastdb per chromosome
    p_map(partial(makeblastdb_from_location, seq_loc=ref_loc, dbdir=refdbdir), ref_chr_list, num_cpus=num_cpus, pbar=False)

    # select cpu proper usage
    suggest_num_cpus = max(1, len(ref_chr_list) // 2)

    if num_cpus <= suggest_num_cpus:
        hard_num_cpus = num_cpus
        loop_num_cpus = 1
    else:
        cpu_usage_list = [num_cpus % i for i in range(suggest_num_cpus - 1, suggest_num_cpus + 2)]
        hard_num_cpus = suggest_num_cpus - 1 + cpu_usage_list.index(min(cpu_usage_list))
        loop_num_cpus = num_cpus // hard_num_cpus

    logprint(f'Task start : {len(sv_loc_list)} SV detected')
    for qry_ind, (qry_loc, sv_loc) in enumerate(zip(qry_loc_list, sv_loc_list)):
        qry_seq = Fasta(qry_loc, build_index=True)

        # If there no same chromosone in query and reference
        qry_chr_list = list(set(ref_chr_list) & set(qry_seq.records.keys()))
        if qry_chr_list == []:
            raise Exception('Chromosone name must same')
        
        qryworkdir = os.path.join(workdir, str(qry_ind))
        os.makedirs(qryworkdir, exist_ok=True)

        sv_list = get_sv_list(sv_loc)
        hom_list = p_map(partial(get_homology, ref_loc=ref_loc, qry_loc=qry_loc, qryworkdir=qryworkdir,
                                 hom_find_len=hom_find_len, temp_indel_find_len=temp_indel_find_len,
                                 near_gap_find_len=near_gap_find_len, user_gap_baseline=user_gap_baseline),
                                 sv_list, pbar=pbar, num_cpus=num_cpus, telegram_token_loc=telegram_token_loc, desc=f'HOM {qry_ind + 1}/{len(qry_loc_list)}')
        
        # check SV to analysis hard mode
        hard_sv_list = []
        for hom in hom_list:
            if hom[1] == 'SUB' and hom[2] == 'SUB_HARD':
                hard_sv_list.append(sv_list[hom[0]])

        hard_hom_list = p_map(partial(get_homology_hard, ref_loc=ref_loc, qry_loc=qry_loc, qryworkdir=qryworkdir, refdbdir=refdbdir,
                                      ref_chr_list=ref_chr_list, hom_find_len=hom_find_len, near_seq_kb_baseline=near_seq_kb_baseline,
                                      diff_locus_hom_baseline=diff_locus_hom_baseline, num_cpus=hard_num_cpus),
                                      hard_sv_list, pbar=False, num_cpus=loop_num_cpus)
        
        for hom in hard_hom_list:
            hom_list[hom[0]] = hom

        output_data = []
        if file:
            qry_basename = os.path.basename(qry_loc)
            with open(os.path.join(dsbr_save, qry_basename) + '.DSBR.tsv', 'w') as f:
                tf = csv.writer(f, delimiter='\t')
                tf.writerow(('ID', 'SV_TYPE', 'CHR', 'REF_START', 'REF_END', 'QRY_START', 'QRY_END', 'SV_TYPE', 'REPAIR_TYPE', 'HOM(L)', 'HOM(R)', 'DSBR_CHR', 'DSBR_START', 'DSBR_END'))
                tf.writerows([sv + list(hom)[1:] if len(sv) == 7 else sv for sv, hom in zip(sv_list, hom_list)])
        else:
            output_data.append(hom_list)
        
        logprint(f'{qry_ind + 1}/{len(qry_loc_list)} : {os.path.basename(qry_loc)} analysis complete')
    
    if not file:
        return output_data