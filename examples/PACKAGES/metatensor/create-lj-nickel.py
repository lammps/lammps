# https://github.com/Luthaf/metatensor-lj-test/
import metatensor_lj_test


model = metatensor_lj_test.lennard_jones_model(
    atomic_type=28,
    cutoff=6.5,
    sigma=1.5808,
    epsilon=0.1729,
    length_unit="Angstrom",
    energy_unit="eV",
    with_extension=False,
)

model.export("nickel-lj.pt")


model = metatensor_lj_test.lennard_jones_model(
    atomic_type=28,
    cutoff=6.5,
    sigma=1.5808,
    epsilon=0.1729,
    length_unit="Angstrom",
    energy_unit="eV",
    with_extension=True,
)
model.export("nickel-lj-extensions.pt", collect_extensions="collected-extensions/")
