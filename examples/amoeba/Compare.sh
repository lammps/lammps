# compare test runs to original runs

# dimer

kdiff log.water_dimer.amoeba.1 log.water_dimer.amoeba.1.test
kdiff dump.water_dimer.amoeba.1 dump.water_dimer.amoeba.1.test

kdiff log.water_dimer.amoeba.4 log.water_dimer.amoeba.4.test
kdiff dump.water_dimer.amoeba.4 dump.water_dimer.amoeba.4.test

kdiff log.water_dimer.hippo.1 log.water_dimer.hippo.1.test
kdiff dump.water_dimer.hippo.1 dump.water_dimer.hippo.1.test

kdiff log.water_dimer.hippo.4 log.water_dimer.hippo.4.test
kdiff dump.water_dimer.hippo.4 dump.water_dimer.hippo.4.test

# hexamer

kdiff log.water_hexamer.amoeba.1 log.water_hexamer.amoeba.1.test
kdiff dump.water_hexamer.amoeba.1 dump.water_hexamer.amoeba.1.test

kdiff log.water_hexamer.amoeba.4 log.water_hexamer.amoeba.4.test
kdiff dump.water_hexamer.amoeba.4 dump.water_hexamer.amoeba.4.test

kdiff log.water_hexamer.hippo.1 log.water_hexamer.hippo.1.test
kdiff dump.water_hexamer.hippo.1 dump.water_hexamer.hippo.1.test

kdiff log.water_hexamer.hippo.4 log.water_hexamer.hippo.4.test
kdiff dump.water_hexamer.hippo.4 dump.water_dimer.hippo.4.test

# water box

kdiff log.water_box.amoeba.1 log.water_box.amoeba.1.test
kdiff dump.water_box.amoeba.1 dump.water_box.amoeba.1.test

kdiff log.water_box.amoeba.32 log.water_box.amoeba.32.test
kdiff dump.water_box.amoeba.32 dump.water_box.amoeba.32.test

kdiff log.water_box.hippo.1 log.water_box.hippo.1.test
kdiff dump.water_box.hippo.1 dump.water_box.hippo.1.test

kdiff log.water_box.hippo.32 log.water_box.hippo.32.test
kdiff dump.water_box.hippo.32 dump.water_box.hippo.32.test

# ubiquitin

kdiff log.ubi.1 log.ubi.1.test
kdiff dump.ubi.1 dump.ubi.1.test

kdiff log.ubi.32 log.ubi.32.test
kdiff dump.ubi.32 dump.ubi.32.test

