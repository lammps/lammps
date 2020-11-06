README for MLIAPPY Example

This example runs the Ta06 example from the MLIAP example, but using the python coupling.

To run this, first run convert_mliap_Ta06A.py, which will convert the Ta06 potential into a pytorch model.
It will be saved as "Ta06A.mliap.pytorch.model.pkl".

It will also copy the "torchlink.py" file into the current working directory. torchlink.py contains
class definitions suitable for wrapping an arbitrary energy model MLIAPPY.

From that point you can run the example lmp -in in.mliap.snap.Ta06A -echo both