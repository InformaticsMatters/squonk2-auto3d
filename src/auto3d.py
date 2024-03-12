#!/usr/bin/env python
"""Template echo utility.
"""
import argparse
import os
import multiprocessing as mp
import time
import logging
import tempfile
import mmap
from io import StringIO
import torch
import shutil
import math

from typing import Union, List
from pathlib import Path

# Import Data Manager DmLog utility.
# Messages emitted using this result in Task Events.
# from dm_job_utilities.dm_log import DmLog

from rdkit import Chem
import rdkit_utils

from Auto3D.auto3D import options, smiles2mols, main
from dm_job_utilities.dm_log import DmLog
import gc

# logger = logging.getLogger(__name__)

SDF_SEP = '$$$$'
CHUNK_SIZE = 4096


path = Path(__file__).parent.absolute()
# smiles = ["CCNCC", "O=C(C1=CC=CO1)N2CCNCC2"]
# args = options(k=1, use_gpu=False)
# mols = smiles2mols(smiles, args)

                
def read_input(input_file_path):
    DmLog.emit_event(Path(".").absolute())
    # temp_dir = tempfile.mkdtemp(dir='data')
    temp_dir = tempfile.mkdtemp(dir=path)
    # temp_dir = tempfile.mkdtemp()
    print('tempdir', temp_dir)
    result = []
    count = 0
    ext = Path(input_file_path).suffix
    if ext in ('.smi', '.smiles'):
        sep = '\n'
    elif ext in ('.sdf', '.sd'):
        sep = SDF_SEP
    else:
        DmLog.emit_event(f'Unknown file type: {ext}')
        raise TypeError
    
    with open(input_file_path, 'r') as input_file:
        count = 0
        buff = StringIO()
        for line in input_file:
            buff.write(line)
            if sep in line:
                count += 1
                fname = str(Path(temp_dir).joinpath(f'input_{str(count)}{ext}'))
                with open(fname, 'w') as output_file:
                    output_file.write(buff.getvalue())
                buff = StringIO()
                result.append(fname)

    return result                

def concat_output(output_file_path, opt_result):
    with open(output_file_path,'ab') as wfd,  open(opt_result,'rb') as fd:
        shutil.copyfileobj(fd, wfd)    

def diagnostics():
    while True:
        # print('proc', proc)
        # print('procdir', dir(proc))
        # print(proc.__dict__)
        time.sleep(10)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Auto3D Package")

    # Input parameters
    # parser.add_argument("--path", type=str, help="A input.smi containing SMILES and IDs")
    parser.add_argument("path", type=str, help="A path of .smi or .SDF file to store all molecules and IDs.")
    parser.add_argument("--output", default="result.sdf", help="The output file")    
    parser.add_argument(
        "--k", default=1, type=int,  help="Outputs the top-k structures for each SMILES"
    )
    parser.add_argument(
        "--window",
        action="store_true",
        help="Outputs the structures whose energies are within x (kcal/mol) from the lowest energy conformer. Only one of --k and --window need to be specified.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="When True, save all meta data while running",
    )
    parser.add_argument(
        "--job_name", type=str, default="", help="A folder name to save all meta data"
    )
    parser.add_argument(
        "--enumerate_tautomer",
        action="store_true",
        help="When True, enumerate tautomers for the input",
    )
    parser.add_argument(
        "--tauto_engine",
        type=str,
        default="rdkit",
        choices=["rdkit", "oechem"],
        help="Programs to enumerate tautomers",
    )
    parser.add_argument(
        "--pKaNorm",
        action="store_true",
        default=True,
        help="When True, the ionization state of each tautomer will be assigned to a predominant state at ~7.4 (Only works when tauto_engine=’oechem’)",
    )
    parser.add_argument(
        "--isomer_engine",
        type=str,
        default="rdkit",
        choices=["rdkit", "omega"],
        help="The program for generating 3D isomers for each SMILES",
    )
    parser.add_argument(
        "--enumerate_isomer",
        action="store_true",
        default=True,
        help="When True, cis/trans and r/s isomers are enumerated",
    )
    parser.add_argument(
        "--mode_oe",
        type=str,
        default="classic",
        choices=["classic", "macrocycle", "dense", "pose", "rocs", "fast_rocs"],
        help="The mode that omega program will take",
    )
    parser.add_argument(
        "--mpi_np",
        type=int,
        default=4,
        help="Number of CPU cores for the isomer generation engine",
    )
    parser.add_argument(
        "--max_confs",
        type=int,
        default=None,
        help="Maximum number of isomers for each SMILES",
    )
    parser.add_argument(
        "--use_gpu",
        action="store_true",
        default=False,
        help="If True, the program will use GPU when available",
    )
    parser.add_argument("--gpu_idx", type=int, nargs="+", default=0, help="GPU index")
    parser.add_argument(
        "--capacity",
        type=int,
        default=42,
        help="Number of SMILES that the model will handle for 1 G memory",
    )
    parser.add_argument(
        "--optimizing_engine",
        type=str,
        default="AIMNET",
        choices=["ANI2x", "ANI2xt", "AIMNET"],
        help="Choose either ANI2x, ANI2xt, or AIMNET for energy calculation and geometry optimization",
    )
    parser.add_argument(
        "--patience",
        type=int,
        default=1000,
        help="If the force does not decrease for a continuous patience steps, the conformer will drop out of the optimization loop",
    )
    parser.add_argument(
        "--opt_steps",
        type=int,
        default=5000,
        help="Maximum optimization steps for each structure",
    )
    parser.add_argument(
        "--convergence_threshold",
        type=float,
        default=0.003,
        help="Optimization is considered as converged if maximum force is below this threshold",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.3,
        help="If the RMSD between two conformers are within threshold, they are considered as duplicates",
    )
    parser.add_argument(
        "--memory",
        type=int,
        default=None,
        help="The RAM size assigned to Auto3D (unit GB)",
    )
    parser.add_argument(
        "--batchsize_atoms",
        type=int,
        default=1024,
        help="The number of atoms in 1 optimization batch for 1GB",
    )

    args = parser.parse_args()

    opts = options(
        args.path,
        k=args.k,
        window=args.window,
        verbose=args.verbose,
        job_name=args.job_name,
        enumerate_tautomer=args.enumerate_tautomer,
        tauto_engine=args.tauto_engine,
        pKaNorm=args.pKaNorm,
        isomer_engine=args.isomer_engine,
        enumerate_isomer=args.enumerate_isomer,
        mode_oe=args.mode_oe,
        mpi_np=args.mpi_np,
        max_confs=args.max_confs,
        use_gpu=args.use_gpu,
        gpu_idx=args.gpu_idx,
        capacity=args.capacity,
        optimizing_engine=args.optimizing_engine,
        patience=args.patience,
        opt_steps=args.opt_steps,
        convergence_threshold=args.convergence_threshold,
        threshold=args.threshold,
        memory=args.memory,
        batchsize_atoms=args.batchsize_atoms,
    )

    num_outputs = 0
    count = -1

    files = read_input(args.path)
    print(files)
    # print(asdfg)


    for f in files:
        print(f)
        count += 1
        opts['path'] = str(f)
        
        try:
            result = main(opts)
            concat_output(args.output, result)
        except RuntimeError as exc:
            print('\n\n\nCuda error?\n\n\n', exc)
        except OSError as exc:
            # this seems to happen with incomplete opt, not enough steps
            print('\n\n\nos error\n\n\n', exc)
        except EOFError as exc:
            # idk, it thows it every time
            print('\n\n\neof error\n\n\n', exc)
        finally:
            torch.cuda.empty_cache()
            gc.collect()
            # this is not enough, BUT it doesn't happen with the
            # molecule alone. garbage collect/kill the thread
            # somehow?

    
    # diag_proc = mp.Process(target=diagnostics)
    # diag_proc.start()
    

    os.chmod(args.output, 0o664)    
    
