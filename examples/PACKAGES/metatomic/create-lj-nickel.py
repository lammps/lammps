try:
    import metatomic_lj_test
except ImportError as e:
    raise ImportError(
        "could not import metatomic_lj_test, please install it with "
        "`pip install git+https://github.com/metatensor/lj-test/`"
    ) from e


model = metatomic_lj_test.lennard_jones_model(
    atomic_type=28,
    cutoff=6.5,
    sigma=1.5808,
    epsilon=0.1729,
    length_unit="Angstrom",
    energy_unit="eV",
    with_extension=False,
)

model.save("nickel-lj.pt")
print("created 'nickel-lj.pt' model")


model = metatomic_lj_test.lennard_jones_model(
    atomic_type=28,
    cutoff=6.5,
    sigma=1.5808,
    epsilon=0.1729,
    length_unit="Angstrom",
    energy_unit="eV",
    with_extension=True,
)
model.save("nickel-lj-extensions.pt", collect_extensions="collected-extensions/")
print("created 'nickel-lj-extensions.pt' model")
