# idxrange.py

* incrementally computes range of fields in a .idx volume

## Usage

```
$ python idxrange.py --help
usage: idxrange.py [-h] [-f FIELDS [FIELDS ...]] [--min MINVAL] [--max MAXVAL]
                   [--maxmem MAXMEM] [--timestep TIMESTEP] [--profile]
                   idxpath

Incrementally computes ranges of fields within an IDX volume.

positional arguments:
  idxpath               path/url of IDX volume

optional arguments:
  -h, --help            show this help message and exit
  -f FIELDS [FIELDS ...], --fields FIELDS [FIELDS ...]
                        list of fields for which to compute min/max (default:
                        all fields)
  --min MINVAL          minimum value to include (lower values will be
                        ignored). NOTE: if different min/max required for each
                        fields, please run this utility separated for each
                        field
  --max MAXVAL          maxiumu value to include (higher values will be
                        ignored). NOTE: if different min/max required for each
                        fields, please run this utility separated for each
                        field
  --maxmem MAXMEM       maximum amount of memory (in MB) to be used by this
                        utility
  --timestep TIMESTEP   timestep for which to compute range
  --profile             profile execution time of this utility
```
