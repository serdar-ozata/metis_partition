### Compiling
```bash
mkdir bin && cd bin
cmake ..
make
```
### CMake Options

The following CMake options are available, all defaulting to `~/local`:

- `PARMETIS_PATH`: Path to the ParMETIS library
- `METIS_PATH`: Path to the METIS library
- `GKLIB_PATH`: Path to the GKlib library

### Running
Usage: `bin/metis_partition <input_file> <num_partitions> <output_directory>`

Example:
```bash
# Assuming you are in the project root directory
mpirun -np 4  bin/metis_partition /path/to/your/datasets/pattern1.mtx 64 /your/output/directory
```
You can download pattern1 from the [SuiteSparse Matrix Collection](https://sparse.tamu.edu/Andrianov/pattern1).



