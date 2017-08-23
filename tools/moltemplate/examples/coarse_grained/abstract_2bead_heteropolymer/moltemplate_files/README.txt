#   In this example, we define two types of molecules: "H" and "P",
#   both containing two atoms, whose ids (names) are "ca" and "r",
#   and whose atom-types vary.
#
#             "H" molecules:          "P" molecules:
#            ("hydrophobic")          ("polar"/"hydrophilic")
#
#                @HR                     @PR
#                 |                       |
#                @CA                     @CA
#
#   Eventually, we will connect multiple "H" and "P" molecules
#   together to form a polymer, as shown below:
#
#                    @HR         @HR
#                     |           |
#                   _@CA_       _@CA_
#       ...  -.@CA-'     `-@CA-'     `  ...
#               |          |
#              @PR         @PR
#
#   The "H" and "P" molecules both share the same type of
#   backbone atom ("CA"), but have their own custom "r"
#   sidechain atoms with different properties:
#   The "HR" atoms belonging to "H" molecules are attracted to each other.
#   The "PR" atoms in "P" molecules are not.
