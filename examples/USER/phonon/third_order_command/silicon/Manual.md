# third_order command

## Syntax

```
third_order group-ID style args keyword value ...
```

* group-ID = ID of group of atoms to displace
* style = *regular* or *ballistico*
```
*regular* args = gamma
  gamma = finite difference displacement length
*ballistico* args = gamma
  gamma = finite difference displacement length
```
* zero or more keyword/value pairs may be appended
* keyword = *file* or *binary*
```
*file* value = output_file
  output_file = name of file to dump the dynamical matrix into
*binary* values = *no* or *gzip*
```

## Examples

```
third_order 1 regular 0.000001
third_order 1 ballistico 0.000001
third_order 3 regular 0.00004 file third_order.dat
third_order 5 ballistico 0.00000001 file third_order.dat binary gzip
```

## Description

Calculate the finite difference third order tensor of the selected group.

## Restrictions

None

## Related commands

None

## Default

The option defaults are file = "third_order.dat", binary = no  
