# Job: auto3d

This describes how to run the `auto3d` job in the `im-auto3d` collection.

## What the job does

This job optimises 3D geometry and generates low-energy conformers of the input 2D or 3D structure.

For more on the Auto3D package see the [Auto3D documentation](https://auto3d.readthedocs.io/en/latest/index.html).


## Implementation details

* Job implementation: [auto3d.py](/src/auto3d.py)
* Job definition: `auto3d` in [jobs.yaml](../../jobs.yaml)


## How to run the job

### Inputs

* **Molecules to optimise**: input molecules in sdf or SMILES format.

### Options
* **Output file**: name of the returned file. Default: `result.sdf`.
* **k**: Outputs the top-k structures for each input. Default `1`.
* **Window**: Outputs the structures whose energies are within x (kcal/mol) from the lowest energy conformer. Only one of --k and --window need to be specified.
* **Enumerate tautomers**: When True, enumerate tautomers for the input. Default `false`
* **Enumerate isomers**: When True, cis/trans and r/s isomers are enumerated. Default `true`.
* **Max conformers**: Maximum number of isomers for each input.
* **Optimizing engine**: Choose either ANI2x, ANI2xt, or AIMNET for energy calculation and geometry optimization. Default `AIMNET`.
* **Patience**: If the force does not decrease for a continuous patience steps, the conformer will drop out of the optimization loop. Default `1000`.
* **Optimisation steps**: Maximum optimization steps for each structure. Default `5000`.
* **Convergence threshold**: Optimization is considered as converged if maximum force is below this threshold. Default `0.003`.
* **Threshold**: If the RMSD between two conformers are within threshold, they are considered as duplicates. Default `0.3`.
* **Delimiter**: Delimiter used in SMILES file. Default: tab. Ignored in sdf.
* **ReadHeader**: Read header from the SMILES input file. Default: False. Ignored in sdf.


## Limitations

1. Currently, only `rdkit` is supported for both tautomer- and isomer generation.

## Related topics

* [Virtual screening](https://github.com/InformaticsMatters/virtual-screening)
