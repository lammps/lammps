This folder contains files and tools for creating, modifying, and validating
files in JSON format. This is work in progress while we are adding JSON
support.

# JSON file format validation.

## JSON-Schema files

We provide schema files for the file formats that LAMMPS supports following
the specifications available on [JSON-Schema](https://json-schema.org) webpage.
The following files are currently available.

- `molecule-schema.json`  Schema file for JSON-format molecule files.
- `dump-molecules-schema.json`  Schema file for the 'dump' format in the 'molecules' style.
- `color-schema.json` Schema file for the colors and lights file used by 'dump image/movie' and in LAMMPS-GUI
- `force-style-test-schema.json` Schema file for the force-style regression test reference files in `unittest/force-styles/`

These files provide a concise description of the hierarchy and supported fields
in JSON file formats.  Thus they provide a detailed documentation and can also
be used for validating JSON files.

## Validation of JSON files

There are multiple tools for JSON file validation available.  Here are instructions
for how to use a tool called `check-jsonschema` which is available via
[PyPi](https://pypi.org/).

``` bash
# Installation into a virtual environment.
# Once installed only the activation should be needed
python3 -m venv validate-json
source validate-json/bin/activate
pip install --upgrade pip
pip install check-jsonschema

# Validation of two molecule files "rxn1.json" and "twomols.json" with "molecule-schema.json"
check-jsonschema --schemafile molecule-schema.json rxn1.json twomols.json
```

If the files are conforming there should be the output:
```
ok -- validation done
```
Otherwise details about the non-conforming fields are given.

## Validation of YAML files

The same `check-jsonschema` tool can validate YAML files, since YAML is a
superset of JSON.  The force-style regression test reference files in
`unittest/force-styles/` are described by `force-style-test-schema.json` and can
be validated the same way as JSON files:

``` bash
# Validate all force-style test reference files against their schema.
check-jsonschema --schemafile force-style-test-schema.json \
    ../../unittest/force-styles/tests/*.yaml
```

A conforming set of files again produces the `ok -- validation done` output;
otherwise the offending file and field are listed.

These reference files are written by the `YamlWriter` class
(`unittest/force-styles/yaml_writer.cpp`) in a regular YAML format that standard
YAML loaders parse without special handling: block scalars use the plain `|`
literal style and empty values are written as `""`, so the default string tag is
left implicit (no `! ` tag indicator).  Older reference files that still use the
non-specific `! ` tag are read identically and remain valid.

# JSON file format normalization

There are extensions to the strict JSON format that allow for comments
or ignore additional (dangling) commas. The ``reformat-json.cpp`` tool
will read JSON files in relaxed format, but write it out in strict format.
It is also possible to change the level of indentation from -1 (all data
one long line) to any positive integer value.  The original file will be
backed up (.bak added to file name) and then overwritten.

Manual compilation (it will be automatically included in the CMake build
if building tools is requested during CMake configuration):

```bash
g++ -I <path/to/lammps/src> -o reformat-json reformat-json.cpp
```

-------

updated by Axel Kohlmeyer, 2025-05-23
