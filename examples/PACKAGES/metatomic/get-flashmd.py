# /// script
# dependencies = [
#   "flashmd == 0.2.*",
# ]
# ///

# Run this script to get the a pre-trained flashmd model for `fix metatomic`
#
# Usage: python get-flashmd.py in a virtual environment with
# https://github.com/lab-cosmo/flashmd installed or `uv run get-flashmd.py`

import flashmd


energy_model, flashmd_model = flashmd.get_pretrained("pet-omatpes", time_step=16)

energy_model.save("flashmd-energy-model.pt")
flashmd_model.save("flashmd-dynamics-model.pt")
