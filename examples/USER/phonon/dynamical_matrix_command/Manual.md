# dynamical_matrix command

## Syntax

```
dynamical_matrix group-ID style args keyword value ...
```

* group-ID = ID of group of atoms to displace
* style = *regular* or *eskm*
```
*regular* args = gamma
  gamma = finite difference displacement length
*eskm* args = gamma
  gamma = finite difference displacement length
```
* zero or more keyword/value pairs may be appended
* keyword = *file* or *binary*
```
*file* value = output_file
  output_file = name of file to dump the dynamical matrix into
*binary* values = *yes* or *no* or *gzip*
```

## Examples

```
dynamical_matrix 1 regular 0.000001
dynamical_matrix 1 eskm 0.000001
dynamical_matrix 3 regular 0.00004 file dynmat.dat
dynamical_matrix 5 eskm 0.00000001 file dynamical.dat binary yes
```

## Description

Calculate the dynamical matrix of the selected group.

## Restrictions

None

## Related commands

None

## Default

The option defaults are file = "dynmat.dyn", binary = no  
