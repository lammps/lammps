Output structured data from LAMMPS
##################################

LAMMPS can output structured data with the :doc:`print <print>` and :doc:`fix
print <fix_print>` command.  This gives you flexibility since you can build
custom data formats that contain system properties, thermo data, and variables
values. This output can be directed to the screen and/or to a file for post
processing.

Writing the current system state, thermo data, variable values
==============================================================

Use the :doc:`print <print>` command to output the current system state, which
can include system properties, thermo data and variable values.

YAML
----

.. code-block:: LAMMPS

   print """---
   timestep: $(step)
   pe: $(pe)
   ke: $(ke)""" file current_state.yaml screen no

.. code-block:: yaml
   :caption: current_state.yaml

   ---
   timestep: 250
   pe: -4.7774327356321810711
   ke: 2.4962152903997174569

JSON
----

.. code-block:: LAMMPS

   print """{
     "timestep": $(step),
     "pe": $(pe),
     "ke": $(ke)
   }""" file current_state.json screen no

.. code-block:: JSON
   :caption: current_state.json

   {
     "timestep": 250,
     "pe": -4.7774327356321810711,
     "ke": 2.4962152903997174569
   }


Writing continuous data during a simulation
===========================================

The :doc:`fix print <fix_print>` command allows you to output an arbitrary string at defined times during a simulation run.

YAML
----

.. code-block:: LAMMPS

   fix extra all print 50 """
   - timestep: $(step)
     pe: $(pe)
     ke: $(ke)""" file output.yaml screen no

.. code-block:: yaml
   :caption: output.yaml

   # Fix print output for fix extra
   - timestep: 0
     pe: -6.77336805325924729
     ke: 4.4988750000000026219

   - timestep: 50
     pe: -4.8082494418323200591
     ke: 2.5257981827119797558

   - timestep: 100
     pe: -4.7875608875581505686
     ke: 2.5062598821985102582

   - timestep: 150
     pe: -4.7471033686005483787
     ke: 2.466095925545450207

   - timestep: 200
     pe: -4.7509052858544134068
     ke: 2.4701136792591693592

   - timestep: 250
     pe: -4.7774327356321810711
     ke: 2.4962152903997174569

Post-processing of YAML files can be easily be done with Python and other
scripting languages. In case of Python the `yaml` package allows you to load the
data files and obtain a list of dictionaries.

.. code-block:: python

   import yaml

   with open("output.yaml") as f:
      data = yaml.load(f, Loader=yaml.FullLoader)

   print(data)

.. code-block::

   [{'timestep': 0, 'pe': -6.773368053259247, 'ke': 4.498875000000003}, {'timestep': 50, 'pe': -4.80824944183232, 'ke': 2.5257981827119798}, {'timestep': 100, 'pe': -4.787560887558151, 'ke': 2.5062598821985103}, {'timestep': 150, 'pe': -4.747103368600548, 'ke': 2.46609592554545}, {'timestep': 200, 'pe': -4.750905285854413, 'ke': 2.4701136792591694}, {'timestep': 250, 'pe': -4.777432735632181, 'ke': 2.4962152903997175}]

Line Delimited JSON (LD-JSON)
-----------------------------

The JSON format itself is very strict when it comes to delimiters. For continuous
output/streaming data it is beneficial use the *line delimited JSON* format.
Each line represents one JSON object.

.. code-block:: LAMMPS

   fix extra all print 50 """{"timestep": $(step), "pe": $(pe), "ke": $(ke)}""" title "" file output.json screen no

.. code-block:: json
   :caption: output.json

   {"timestep": 0, "pe": -6.77336805325924729, "ke": 4.4988750000000026219}
   {"timestep": 50, "pe": -4.8082494418323200591, "ke": 2.5257981827119797558}
   {"timestep": 100, "pe": -4.7875608875581505686, "ke": 2.5062598821985102582}
   {"timestep": 150, "pe": -4.7471033686005483787, "ke": 2.466095925545450207}
   {"timestep": 200, "pe": -4.7509052858544134068, "ke": 2.4701136792591693592}
   {"timestep": 250, "pe": -4.7774327356321810711, "ke": 2.4962152903997174569}

One simple way to load this data into a Python script is to use the `pandas`
package. It can directly load these files into a data frame:

.. code-block:: python

   import pandas as pd

   data = pd.read_json('output.json', lines=True)
   print(data)

.. code-block:: bash

      timestep        pe        ke
   0         0 -6.773368  4.498875
   1        50 -4.808249  2.525798
   2       100 -4.787561  2.506260
   3       150 -4.747103  2.466096
   4       200 -4.750905  2.470114
   5       250 -4.777433  2.496215
