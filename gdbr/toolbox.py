from pathos.multiprocessing import ProcessPool
from tqdm.contrib.telegram import tqdm

import p_tqdm
import json
import os


def get_telegram_data(file_loc):
    if os.path.exists(file_loc):
        with open(file_loc, 'r') as f:
            data = json.load(f)
            return data['token'], data['chat_id']
    else:
        return None


def p_map(f, it, num_cpus=1, pbar=True, telegram_token_loc='telegram.json', desc=''):
    if pbar:
        telegram_data = get_telegram_data(telegram_token_loc)
        if telegram_data is None:
            return p_tqdm.p_map(f, it, num_cpus=num_cpus)
        else:
            return p_tqdm.p_map(f, it, num_cpus=num_cpus, tqdm=tqdm, token=telegram_data[0], chat_id=telegram_data[1], desc=desc, file=open(os.devnull, 'w'))
    else:
        pool = ProcessPool(num_cpus)
        return pool.map(f, it)
