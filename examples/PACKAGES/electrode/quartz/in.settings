pair_style	hybrid/overlay morse 9.0 coul/long 12.0    #
kspace_style	pppm/electrode 1.0e-6
# -------------------------------------------------------- #
# ------------------------- BOX -------------------------- #
# -------------------------------------------------------- #
read_data	data.quartz
#
pair_coeff      2 2 morse 0.29560 1.717429552 3.4103 # Si-Si
pair_coeff      1 1 morse 0.53630 1.375868904 3.7835  # O-O
pair_coeff      1 2 morse 45.9970 2.725476839 1.6148 # Si-O
pair_coeff	* * coul/long
# -------------------------------------------------------- #
# ------------------------ GROUPS ------------------------ #
# -------------------------------------------------------- #
group all type	1 2
group si type	2
group o type	1
# -------------------------------------------------------- #
# ------------------ compute charges  -------------------- #
# -------------------------------------------------------- #	
variable q atom q
compute qall all reduce sum v_q
compute qsi si reduce sum v_q
compute qo o reduce sum v_q
#
variable charge equal c_qall
variable charge_qsi equal c_qsi/count(si)
variable charge_qo equal c_qo/count(o)
# -------------------------------------------------------- #
