#   In this example, we construct a large molecule ("Polymer")
#   from multiple smaller molecular subunits ("Monomer").
#   The "Monomer" molecule contains two atoms (type "CA", and "R")
#
#             "Monomer"
#
#                @R
#                 |
#                @CA
#
#   Eventually, we will connect multiple "Monomer" molecular subunits
#   together to form a polymer, as shown below:
#
#                     @R          @R
#                     |           |
#                   _@CA_       _@CA_
#       ...  -.@CA-'     `-@CA-'     `  ...
#               |          |
#               @R         @R
