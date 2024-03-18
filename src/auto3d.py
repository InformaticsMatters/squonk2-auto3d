#!/usr/bin/env python
"""Template echo utility.
"""
import argparse
import os
import tempfile
import shutil

from pathlib import Path

# Import Data Manager DmLog utility.
# Messages emitted using this result in Task Events.
# from dm_job_utilities.dm_log import DmLog

import rdkit_utils

from Auto3D.auto3D import options, main
from dm_job_utilities.dm_log import DmLog

# logger = logging.getLogger(__name__)

SDF_SEP = '$$$$'
CHUNK_SIZE = 4096


path = Path(__file__).parent.absolute()

# for local development, data directory exists in repo but not copied to container
try:
    path = path.parent.joinpath('data')
except NameError:
    pass


def process_input(input_filename, write_header, delimiter, read_header, id_column, sdf_read_records):
    result = []
    temp_dir = tempfile.mkdtemp(dir=path)
    reader = rdkit_utils.create_reader(
        input_filename,
        delimiter=delimiter,
        read_header=read_header,
        id_column=id_column,
        sdf_read_records=sdf_read_records,
    )

    extra_field_names = reader.get_extra_field_names()
    calc_prop_names = []

    ext = Path(input_filename).suffix

    count = 0
    while True:
        count += 1
        try:
            mol, smi, mol_id, props = reader.read()
        except TypeError as ex:
            DmLog.emit_event(f"{ex}")
            continue
        except StopIteration as ex:
            # end of file
            break    

        output_filename = str(Path(temp_dir).joinpath(f'input_{str(count)}{ext}'))
        result.append(output_filename)
        writer = rdkit_utils.create_writer(
            output_filename,
            delimiter=delimiter,
        )

        if write_header:
            headers = rdkit_utils.generate_header_values(extra_field_names, len(props), calc_prop_names)
            writer.write_header(headers)

        writer.write(
            smiles=smi,
            mol=mol,
            mol_id=mol_id,
            existing_props=props,  # only used in SmilesWriter
        )

        writer.close()
        
    reader.close()
    return result

      
        
def concat_output(output_file_path, opt_result):
    with open(opt_result,'rb') as opt_file:
        count = opt_file.read().count(b'$$$$')
        with open(output_file_path,'ab') as result_file:
            shutil.copyfileobj(opt_file, result_file)
    return count



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Auto3D Package")

    # Input parameters
    parser.add_argument("path", type=str, help="A path of .smi or .SDF file to store all molecules and IDs.")
    parser.add_argument("--output", default="result.sdf", help="The output file")    
    parser.add_argument(
        "--k", default=1, type=int,  help="Outputs the top-k structures for each SMILES"
    )
    parser.add_argument(
        "--window",
        type=float,
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
        default=40,
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
    parser.add_argument("-d", "--delimiter", default="\t", help="Delimiter when using SMILES")
    parser.add_argument(
        "--id-column",
        help="Column for name field (zero based integer for .smi, text for SDF)",
    )
    parser.add_argument(
        "--read-header",
        action="store_true",
        help="Read a header line with the field names when reading .smi or .txt",
    )
    parser.add_argument(
        "--write-header",
        action="store_true",
        help="Write a header line when writing .smi or .txt",
    )
    parser.add_argument(
        "--sdf-read-records",
        default=100,
        type=int,
        help="Read this many SDF records to determine field names",
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

    files = process_input(
        input_filename=args.path,
        write_header=args.write_header,
        delimiter=args.delimiter,
        read_header=args.read_header,
        id_column=args.id_column,
        sdf_read_records=args.sdf_read_records,        
    )

    DmLog.emit_event('Input files prepared')

    for f in files:
        DmLog.emit_event(f'Processing {f}')
        opts['path'] = str(f)
        
        try:
            result = main(opts)
            num_outputs += concat_output(args.output, result)
        except RuntimeError as exc:
            DmLog.emit_event(f"Error processing {f}, no output produced: {exc}")
        except OSError as exc:
            # this SEEMS to happen with incomplete optimisation, not enough steps
            DmLog.emit_event(f"Error processing {f}, no output produced: {exc}")
        except EOFError as exc:
            DmLog.emit_event(f"Error processing {f}, no output produced: {exc}")

    

    if Path(args.output).is_file():
        os.chmod(args.output, 0o664)

    DmLog.emit_event(num_outputs, "outputs from", len(files), "molecules")
    DmLog.emit_cost(num_outputs)        
    
