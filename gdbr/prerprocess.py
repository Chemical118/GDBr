from gdbr.utilities import p_map, logprint, get_proper_thread, check_file_exist, check_unique_basename, safe_makedirs
from functools import partial

import subprocess
import inspect
import shutil
import os


def preprocess_pipeline(qry_loc, ref_loc, log_save, qry_save, var_save, workdir, min_sv_size, num_cpus, trust_workdir):
    qry_basename = os.path.basename(qry_loc)
    ragtag_save = os.path.join(workdir, 'ragtag_gdbr', qry_basename)
    svim_asm_save = os.path.join(workdir, 'svim_asm_gdbr', qry_basename)

    if not trust_workdir and os.path.isdir(ragtag_save):
        shutil.rmtree(ragtag_save)
    os.makedirs(ragtag_save, exist_ok=True)

    if not trust_workdir and os.path.isdir(svim_asm_save):
        shutil.rmtree(svim_asm_save)
    os.makedirs(svim_asm_save, exist_ok=True)

    module_loc = os.path.dirname(inspect.getfile(inspect.currentframe()))
    with open(os.path.join(log_save, qry_basename + '.log'), 'w') as f:
        pipeline_result = subprocess.run([f'bash -e {os.path.join(module_loc, "script", "pipeline.sh")} {qry_loc} {ref_loc} {qry_save} {var_save} {ragtag_save} {svim_asm_save} {min_sv_size} {num_cpus} {1 if trust_workdir else 0}'], stdout=f, stderr=f, shell=True)

    if pipeline_result.returncode != 0:
        raise Exception(f'{qry_loc} preprocess pipeline failed, please check log in {os.path.join(log_save, qry_basename + ".log")}')
    
    with open(os.path.join(var_save, qry_basename + '.GDBr.preprocess.vcf'), 'r') as f:
        vcf_str_list = f.readlines()
    
    vcf_str_list[0] += f'##gdbr_min_sv_size={min_sv_size}\n'
    with open(os.path.join(var_save, qry_basename + '.GDBr.preprocess.vcf'), 'w') as f:
        f.writelines(vcf_str_list)


def preprocess_main(ref_loc, qry_loc_list, preprocess_save='prepro', workdir='data', min_sv_size=50, num_cpus=1, low_memory=False, trust_workdir=False, pbar=True, telegram_token_loc='telegram.json', overwrite_output=False):
    check_file_exist([[ref_loc], qry_loc_list], ['Reference', 'Raw query'])
    check_unique_basename(qry_loc_list)
    
    safe_makedirs(preprocess_save, overwrite_output)
    
    log_save = os.path.join(preprocess_save, 'log')
    qry_save = os.path.join(preprocess_save, 'query')
    var_save = os.path.join(preprocess_save, 'vcf')

    if not trust_workdir:
        safe_makedirs(workdir)
    
    os.makedirs(qry_save, exist_ok=True)
    os.makedirs(var_save, exist_ok=True)
    os.makedirs(log_save, exist_ok=True)
    
    # select cpu proper usage
    if low_memory:
        preprocess_num_cpus, loop_num_cpus = num_cpus, 1
    else:
        preprocess_num_cpus, loop_num_cpus = get_proper_thread(1 if trust_workdir else 4, num_cpus, len(qry_loc_list))
    
    logprint(f'Task preprocess start : {len(qry_loc_list)} query detected')

    p_map(partial(preprocess_pipeline, ref_loc=ref_loc, log_save=log_save, qry_save=qry_save, var_save=var_save, 
                  workdir=workdir, min_sv_size=min_sv_size, num_cpus=preprocess_num_cpus, trust_workdir=trust_workdir),
                  qry_loc_list, pbar=pbar, num_cpus=loop_num_cpus, telegram_token_loc=telegram_token_loc, desc='PRE')
    
    logprint(f'{"Trust preprocess" if trust_workdir else "Preprocess"} complete')
