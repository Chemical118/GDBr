from Bio.Blast import NCBIXML
from functools import partial
from pyfaidx import Fasta
from my_tqdm import p_map

import subprocess
import csv
import os


def get_sv_list(qry_ind):
    with open(f'sv_call_{qry_ind}.csv', 'r') as f:
        cf = csv.reader(f)
        sv_list = [[int(i[0]), i[1], i[2], int(i[3]), int(i[4]), int(i[5]), int(i[6])] if len(i) == 7 else [int(i[0]), i[1]] for i in [l for l in cf][1:]]
    return sv_list


def get_homology(sv_data, hom_find_len=2000, temp_indel_find_len=100, near_gap_find_len=5):
    if sv_data[1] not in {'DEL', 'INS', 'SUB'}:
        return sv_data[0], 'UNSUP_ID'

    sv_id, tid, chrom, ref_start, ref_end, qry_start, qry_end = sv_data

    blastn_fmt = ['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'bitscore']
    dsb_repair_type = 'NO_HOM'

    ref_seq = Fasta(ref_loc, build_index=False)
    qry_seq = Fasta(qry_loc, build_index=False)

    ref_len = ref_end - ref_start + 1
    qry_len = qry_end - qry_start + 1

    if tid == 'SUB':
        gap_baseline = min(3, ref_len, qry_len)

        # -word_size 7 at -task blastn-short
        ref_left_temp_hom_find_len = ref_right_temp_hom_find_len = qry_left_temp_hom_find_len = qry_right_temp_hom_find_len = 3

        ref_left_hom = ref_right_hom = qry_left_hom = qry_right_hom = ref_bitscore = qry_bitscore = 0

        ref_subject_temp_seq = ref_seq[chrom][ref_start - 1 - temp_indel_find_len:ref_end + temp_indel_find_len]
        qry_subject_temp_seq = qry_seq[chrom][qry_start - 1 - temp_indel_find_len:qry_end + temp_indel_find_len]

        ref_subject_temp_seq_loc = os.path.join(qryworkdir, f'ref_sub_{sv_id}.fasta')
        qry_query_temp_seq_loc = os.path.join(qryworkdir, f'qry_qry_{sv_id}.fasta')

        qry_subject_temp_seq_loc = os.path.join(qryworkdir, f'qry_sub_{sv_id}.fasta')
        ref_query_temp_seq_loc = os.path.join(qryworkdir, f'ref_qry_{sv_id}.fasta')

        ref_blast_result_loc = os.path.join(qryworkdir, f'ref_blast_{sv_id}.xml')
        qry_blast_result_loc = os.path.join(qryworkdir, f'qry_blast_{sv_id}.xml')

        with open(ref_subject_temp_seq_loc, 'w') as f:
            f.write('>' + ref_subject_temp_seq.fancy_name + '\n')
            f.write(str(ref_subject_temp_seq))

        with open(qry_subject_temp_seq_loc, 'w') as f:
            f.write('>' + qry_subject_temp_seq.fancy_name + '\n')
            f.write(str(qry_subject_temp_seq))

        first_flag = True
        while ref_left_hom == ref_left_temp_hom_find_len or ref_right_hom == ref_right_temp_hom_find_len or first_flag:
            first_flag = False

            ref_left_temp_hom_find_len *= 2 if ref_left_hom == ref_left_temp_hom_find_len else 1
            ref_right_temp_hom_find_len *= 2 if ref_right_hom == ref_right_temp_hom_find_len else 1

            qry_query_temp_seq = qry_seq[chrom][qry_start - 1 - ref_left_temp_hom_find_len:qry_end + ref_right_temp_hom_find_len]
            with open(qry_query_temp_seq_loc, 'w') as f:
                f.write('>' + qry_query_temp_seq.fancy_name + '\n')
                f.write(str(qry_query_temp_seq))

            subprocess.run(['blastn', '-task', 'blastn-short', '-subject', ref_subject_temp_seq_loc, '-query', qry_query_temp_seq_loc, '-strand', 'plus', '-outfmt', '5', '-out', ref_blast_result_loc], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            with open(ref_blast_result_loc, 'r') as f:
                ref_blast_result = NCBIXML.read(f).alignments

            os.remove(ref_blast_result_loc)
            
            ref_blast_result_list = [] if ref_blast_result == [] else ref_blast_result[0].hsps
            ref_blast_filter_result = sorted(filter(lambda t: t.query_start < ref_left_temp_hom_find_len + 1 and ref_left_temp_hom_find_len + qry_len < t.query_end and t.gaps < gap_baseline and t.bits > ref_bitscore, ref_blast_result_list), key=lambda t: t.bits, reverse=True)

            if len(ref_blast_filter_result) == 0:
                break
            
            num = -1
            for blast_result in ref_blast_filter_result:
                num += 1
                sub_real_idx_list = [i for i, v in enumerate(blast_result.sbjct) if v != '-']

                sub_sv_loc_st = temp_indel_find_len + 1 - blast_result.sbjct_start
                sub_sv_loc_nd = temp_indel_find_len + ref_len - blast_result.sbjct_start

                sub_sv_real_loc_st = sub_real_idx_list[sub_sv_loc_st] if 0 <= sub_sv_loc_st < len(sub_real_idx_list) else sub_sv_loc_st
                sub_sv_real_loc_nd = sub_real_idx_list[sub_sv_loc_nd] if 0 <= sub_sv_loc_nd < len(sub_real_idx_list) else sub_sv_loc_nd

                qry_real_idx_list = [i for i, v in enumerate(blast_result.query) if v != '-']

                qry_sv_real_loc_st = qry_real_idx_list[ref_left_temp_hom_find_len + 1 - blast_result.query_start]
                qry_sv_real_loc_nd = qry_real_idx_list[ref_left_temp_hom_find_len + qry_len - blast_result.query_start]

                if sub_sv_real_loc_nd - sub_sv_real_loc_st > qry_sv_real_loc_nd - qry_sv_real_loc_st:
                    short_sv_real_loc_st = qry_sv_real_loc_st
                    short_sv_real_loc_nd = qry_sv_real_loc_nd

                else:
                    short_sv_real_loc_st = sub_sv_real_loc_st
                    short_sv_real_loc_nd = sub_sv_real_loc_nd

                match_sv_seq = blast_result.match[short_sv_real_loc_st:short_sv_real_loc_nd + 1]

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

                if combo_isgap or ' ' in match_sv_seq:
                    continue

                ref_left_hom = ref_left_temp_hom_find_len - blast_result.query_start + 1
                ref_right_hom = blast_result.query_end - qry_len - ref_left_temp_hom_find_len
                ref_bitscore = blast_result.bits

                ref_sub_st = ref_start - 1 - temp_indel_find_len + blast_result.sbjct_start
                ref_sub_nd = ref_start - 1 - temp_indel_find_len + blast_result.sbjct_end
                ref_qry_st = qry_start - 1 - ref_left_temp_hom_find_len + blast_result.query_start
                ref_qry_nd = qry_start - 1 - ref_left_temp_hom_find_len + blast_result.query_end
                break

        if ref_left_hom > 0:
            while str(ref_seq[chrom][ref_sub_st - 2]).upper() == str(qry_seq[chrom][ref_qry_st - 2]).upper():
                ref_sub_st -= 1
                ref_qry_st -= 1
                ref_left_hom += 1


            while str(ref_seq[chrom][ref_sub_nd]).upper() == str(qry_seq[chrom][ref_qry_nd]).upper():
                ref_sub_nd += 1
                ref_qry_nd += 1
                ref_right_hom += 1

        os.remove(ref_subject_temp_seq_loc)
        os.remove(qry_query_temp_seq_loc)

        first_flag = True
        while qry_left_hom == qry_left_temp_hom_find_len or qry_right_hom == qry_right_temp_hom_find_len or first_flag:
            first_flag = False
            
            qry_left_temp_hom_find_len *= 2 if qry_left_hom == qry_left_temp_hom_find_len else 1
            qry_right_temp_hom_find_len *= 2 if qry_right_hom == qry_right_temp_hom_find_len else 1

            ref_query_temp_seq = ref_seq[chrom][ref_start - 1 - qry_left_temp_hom_find_len:ref_end + qry_right_temp_hom_find_len]
            with open(ref_query_temp_seq_loc, 'w') as f:
                f.write('>' + ref_query_temp_seq.fancy_name + '\n')
                f.write(str(ref_query_temp_seq))

            subprocess.run(['blastn', '-task', 'blastn-short', '-subject', qry_subject_temp_seq_loc, '-query', ref_query_temp_seq_loc, '-strand', 'plus', '-outfmt', '5', '-out', qry_blast_result_loc], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            with open(qry_blast_result_loc, 'r') as f:
                qry_blast_result = NCBIXML.read(f).alignments
                
            os.remove(qry_blast_result_loc)
            
            qry_blast_result_list = [] if qry_blast_result == [] else qry_blast_result[0].hsps
            qry_blast_filter_result = sorted(filter(lambda t: t.query_start < qry_left_temp_hom_find_len + 1 and qry_left_temp_hom_find_len + ref_len < t.query_end and t.gaps < gap_baseline and t.bits > qry_bitscore, qry_blast_result_list), key=lambda t: t.bits, reverse=True)

            if len(qry_blast_filter_result) == 0:
                break
            
            num = -1
            for blast_result in qry_blast_filter_result:
                num += 1
                sub_real_idx_list = [i for i, v in enumerate(blast_result.sbjct) if v != '-']

                sub_sv_loc_st = temp_indel_find_len + 1 - blast_result.sbjct_start
                sub_sv_loc_nd = temp_indel_find_len + qry_len - blast_result.sbjct_start

                sub_sv_real_loc_st = sub_real_idx_list[sub_sv_loc_st] if 0 <= sub_sv_loc_st < len(sub_real_idx_list) else sub_sv_loc_st
                sub_sv_real_loc_nd = sub_real_idx_list[sub_sv_loc_nd] if 0 <= sub_sv_loc_nd < len(sub_real_idx_list) else sub_sv_loc_nd

                qry_real_idx_list = [i for i, v in enumerate(blast_result.query) if v != '-']

                qry_sv_real_loc_st = qry_real_idx_list[qry_left_temp_hom_find_len + 1 - blast_result.query_start]
                qry_sv_real_loc_nd = qry_real_idx_list[qry_left_temp_hom_find_len + ref_len - blast_result.query_start]

                if sub_sv_real_loc_nd - sub_sv_real_loc_st > qry_sv_real_loc_nd - qry_sv_real_loc_st:
                    short_sv_real_loc_st = qry_sv_real_loc_st
                    short_sv_real_loc_nd = qry_sv_real_loc_nd

                else:
                    short_sv_real_loc_st = sub_sv_real_loc_st
                    short_sv_real_loc_nd = sub_sv_real_loc_nd

                match_sv_seq = blast_result.match[short_sv_real_loc_st:short_sv_real_loc_nd + 1]

                combo_st = 0
                combo_tar = ''
                combo_isgap = False

                for i in range(blast_result.align_length):
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

                if combo_isgap or ' ' in match_sv_seq:
                    continue

                qry_left_hom = qry_left_temp_hom_find_len - blast_result.query_start + 1
                qry_right_hom = blast_result.query_end - ref_len - qry_left_temp_hom_find_len
                qry_bitscore = blast_result.bits
                
                qry_sub_st = qry_start - 1 - temp_indel_find_len + blast_result.sbjct_start
                qry_sub_nd = qry_start - 1 - temp_indel_find_len + blast_result.sbjct_end
                qry_qry_st = ref_start - 1 - qry_left_temp_hom_find_len + blast_result.query_start
                qry_qry_nd = ref_start - 1 - qry_left_temp_hom_find_len + blast_result.query_end
                break

        if qry_left_hom > 0:
            while str(qry_seq[chrom][qry_sub_st - 2]).upper() == str(ref_seq[chrom][qry_qry_st - 2]).upper():
                qry_sub_st -= 1
                qry_qry_st -= 1
                qry_left_hom += 1


            while str(qry_seq[chrom][qry_sub_nd]).upper() == str(ref_seq[chrom][qry_qry_nd]).upper():
                qry_sub_nd += 1
                qry_qry_nd += 1
                qry_right_hom += 1
        
        os.remove(qry_subject_temp_seq_loc)
        os.remove(ref_query_temp_seq_loc)

        if ref_left_hom > 0 or qry_left_hom > 0:
            dsb_repair_type = 'TEMP_INS'

            if ref_left_hom + ref_right_hom >= qry_left_hom + qry_right_hom:
                left_hom = ref_left_hom
                right_hom = ref_right_hom

            else:
                left_hom = qry_left_hom
                right_hom = qry_right_hom

        else:
            dsb_repair_type = 'SUB'

    else:
        if tid == 'DEL':
            start_temp_seq = ref_seq[chrom][ref_start - 1 - hom_find_len:ref_end]
            end_temp_seq = ref_seq[chrom][ref_start - 1:ref_end + hom_find_len]
            subject_temp_seq = qry_seq[chrom][qry_start - 1 - hom_find_len:qry_end + hom_find_len]

        else:
            start_temp_seq = qry_seq[chrom][qry_start - 1 - hom_find_len:qry_end]
            end_temp_seq = qry_seq[chrom][qry_start - 1:qry_end + hom_find_len]
            subject_temp_seq = ref_seq[chrom][ref_start - 1 - hom_find_len:ref_end + hom_find_len]

        start_temp_seq_loc = os.path.join(qryworkdir, f'start_{sv_id}.fasta')
        end_temp_seq_loc = os.path.join(qryworkdir, f'end_{sv_id}.fasta')
        subject_temp_seq_loc = os.path.join(qryworkdir, f'sub_{sv_id}.fasta')

        with open(start_temp_seq_loc, 'w') as f:
            f.write('>' + start_temp_seq.fancy_name + '\n')
            f.write(str(start_temp_seq))

        with open(end_temp_seq_loc, 'w') as f:
            f.write('>' + end_temp_seq.fancy_name + '\n')
            f.write(str(end_temp_seq))

        with open(subject_temp_seq_loc, 'w') as f:
            f.write('>' + subject_temp_seq.fancy_name + '\n')
            f.write(str(subject_temp_seq))

        start_result = subprocess.run(['blastn', '-subject', subject_temp_seq_loc, '-query', start_temp_seq_loc, '-strand', 'plus', '-outfmt', '10 ' + ' '.join(blastn_fmt)], capture_output=True, text=True)
        end_result = subprocess.run(['blastn', '-subject', subject_temp_seq_loc, '-query', end_temp_seq_loc, '-strand', 'plus', '-outfmt', '10 ' + ' '.join(blastn_fmt)], capture_output=True, text=True)

        os.remove(start_temp_seq_loc)
        os.remove(end_temp_seq_loc)
        os.remove(subject_temp_seq_loc)

        start_result_list =  [] if start_result.stdout == '' else map(lambda t: list(map(eval, t.split(','))), start_result.stdout[:-1].split('\n'))
        end_result_list = [] if end_result.stdout == '' else map(lambda t: list(map(eval, t.split(','))), end_result.stdout[:-1].split('\n'))

        start_filter_result = list(filter(lambda t: t[1] > hom_find_len * 0.90, start_result_list))
        end_filter_result = list(filter(lambda t: t[1] > hom_find_len * 0.90, end_result_list))

        if len(start_filter_result) != 1 or len(end_filter_result) != 1:
            dsb_repair_type = f'FND_HOM:({len(start_filter_result)}, {len(end_filter_result)})'

        else:
            start_filter_result = start_filter_result[0]
            end_filter_result = end_filter_result[0]

            if start_filter_result[6] > end_filter_result[7]:
                dsb_repair_type = 'BST_LOC'
            else:
                hom = start_filter_result[7] - end_filter_result[6] + 1

                if hom > (ref_len if tid == 'DEL' else qry_len):
                    dsb_repair_type = 'HOM_LEN'
                elif hom > 0:
                    dsb_repair_type = 'HOM'


    if dsb_repair_type == 'HOM':
        return sv_id, tid, dsb_repair_type, hom

    elif dsb_repair_type == 'TEMP_INS':
        return sv_id, tid, dsb_repair_type, left_hom, right_hom

    elif dsb_repair_type == 'NH_DSBR':
        return sv_id, tid, dsb_repair_type, left_hom, right_hom

    else:
        return sv_id, tid, dsb_repair_type

# read .fasta file
ref_loc = 'refseq/chm13v2.0.fa'
ref_seq = Fasta(ref_loc, build_index=True)

workdir = 'data'

# query must have a chromosome name (ex by ragtag)
qry_loc_list = ['qryseq/NA21309.maternal.fa', 'qryseq/NA21309.paternal.fa']

# get 1Mbp chromosome
ref_chr_list = list(map(lambda t: t[0], filter(lambda t: len(t[1]) > 1e6, ref_seq.records.items())))

for qry_ind, qry_loc in enumerate(qry_loc_list):
    qry_seq = Fasta(qry_loc, build_index=True)

    # check all ref chromosome in all qry
    if set(ref_chr_list) > set(qry_seq.records.keys()):
        raise Exception('Chromone Error')
    
    qryworkdir = os.path.join(workdir, str(qry_ind))
    dbdir = os.path.join(qryworkdir, 'db')

    sv_list = get_sv_list(qry_ind)

    hom_find_len, temp_indel_find_len = 2000, 100
    hom_list = p_map(partial(get_homology, hom_find_len=hom_find_len, temp_indel_find_len=temp_indel_find_len),
                     sv_list, num_cpus=16)
    
    with open(f'hum_hom_{qry_ind}.csv', 'w') as f:
        cf = csv.writer(f)
        cf.writerow(('ID', 'SV_TYPE', 'REPAIR_TYPE', 'HOM', 'HOM2'))
        cf.writerows(hom_list)