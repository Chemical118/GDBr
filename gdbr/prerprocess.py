from gdbr.utilities import p_map, logprint, get_proper_thread
from functools import partial

import subprocess
import inspect
import shutil
import os


def preprocess_pipeline(qry_loc, ref_loc, log_save, qry_save, var_save, workdir, min_sv_size, num_cpus):
    qry_basename = os.path.basename(qry_loc)
    ragtag_save = os.path.join(workdir, 'ragtag', qry_basename)
    svim_asm_save = os.path.join(workdir, 'svim_asm', qry_basename)

    if os.path.isdir(ragtag_save):
        shutil.rmtree(ragtag_save)
    os.makedirs(ragtag_save, exist_ok=True)

    if os.path.isdir(svim_asm_save):
        shutil.rmtree(svim_asm_save)
    os.makedirs(svim_asm_save, exist_ok=True)

    module_loc = os.path.dirname(inspect.getfile(inspect.currentframe()))
    with open(os.path.join(log_save, qry_basename + '.log'), 'w') as f:
        subprocess.run([f'bash {os.path.join(module_loc, "script", "pipeline.sh")} {qry_loc} {ref_loc} {qry_save} {var_save} {ragtag_save} {svim_asm_save} {min_sv_size} {num_cpus}'], stdout=f, stderr=f, shell=True)


def preprocess_main(ref_loc, qry_loc_list, qry_save='querys', var_save='vcfs', workdir='data', min_sv_size=50, num_cpus=1, low_memory=False, pbar=True, telegram_token_loc='telegram.json'):
    qry_basename_list = list(map(os.path.basename, qry_loc_list))

    if len(qry_basename_list) != len(set(qry_basename_list)):
        raise Exception('Query basename must be diffrent')
    
    log_save = os.path.join(workdir, 'logs')

    os.makedirs(workdir, exist_ok=True)
    os.makedirs(qry_save, exist_ok=True)
    os.makedirs(var_save, exist_ok=True)
    os.makedirs(log_save, exist_ok=True)

    # select cpu proper usage
    if low_memory:
        preprocess_num_cpus, loop_num_cpus = num_cpus, 1
    else:
        preprocess_num_cpus, loop_num_cpus = get_proper_thread(4, num_cpus, len(qry_loc_list))
    
    logprint(f'Task start : {len(qry_loc_list)} query detected')
    p_map(partial(preprocess_pipeline, ref_loc=ref_loc, log_save=log_save, qry_save=qry_save, var_save=var_save, 
                  workdir=workdir, min_sv_size=min_sv_size, num_cpus=preprocess_num_cpus),
                  qry_loc_list, pbar=pbar, num_cpus=loop_num_cpus, telegram_token_loc=telegram_token_loc, desc='PRE')