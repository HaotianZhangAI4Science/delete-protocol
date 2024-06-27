import subprocess
from glob import glob 
import os.path as osp
import time
import argparse


query_sdfs = glob('./KRAS/kras_active/*') #contains ['0.sdf','1.sdf',...,'7.sdf']
num_return = 100
library = 'chemdiv' #['chemdiv','specs']

for query_sdf in query_sdfs:

    start_time = time.time()
    query_name = osp.basename(query_sdf).split('.')[0]
    
    if library == 'specs':
        command = f'python search_similarity.py --query_sdf {query_sdf} \
            --molfname ./library/Specs_ExAcD_Aug_2020.sdf --fpdbfname ./library/Specs_ExAcD_Aug_2020.fpbin\
            --saved_sdf ./KRAS/kras_active_sim_search/{query_name}_specs_sim_{num_return}.sdf --num_return {num_return}'
    elif library == 'chemdiv'
        command = f'python search_similarity.py --query_sdf {query_sdf} \
            --molfname ./library/chemdiv.sdf --fpdbfname ./library/chemdiv.fpbin\
            --saved_sdf ./KRAS/kras_active_sim_search/{query_name}_chemdiv_sim_{num_return}.sdf --num_return {num_return}'
    else:
        print('Do not have the specified library')
        
    print(command)
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    if result.returncode == 0:
        print('executed successfully.')
        print('Output:')
        print(result.stdout)
        print('consumed time: ',time.time()-start_time)
    else:
        print('execution failed.')
        print('Error:')
        print(result.stderr)
