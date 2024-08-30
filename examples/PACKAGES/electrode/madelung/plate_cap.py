#!/usr/bin/env python3

import time

import numpy as np
from scipy.special import erf

SQRT2 = np.sqrt(2)
SQRTPI_INV = 1 / np.sqrt(np.pi)
COULOMB = 332.06371  #  Coulomb constant in Lammps 'real' units
QE2F = 23.060549
NKTV2P = 68568.415  # pressure in 'real' units
LENGTH = 10000  # convergence parameter
LZ = 20


def lattice(length):
    """indices combinations of images in one quadrant"""
    x = np.arange(length)  # range(length)
    y = np.arange(1, length)
    return np.array(np.meshgrid(x, y)).T.reshape(-1, 2)


def a_element(r, eta):
    """Coulomb contribution of two Gaussians"""
    return erf(eta * r) / r


def b_element(r, q, eta):
    """Coulomb contribution of a Gaussian with a point charge"""
    return q * erf(eta * r) / r


def force_gauss(r, qq, eta):
    etar = eta * r
    return (qq / np.square(r)) * (
        erf(etar) - 2 * etar * SQRTPI_INV * np.exp(-np.square(etar))
    )


def force_point(r, qq):
    return qq / np.square(r)


def force_component(dx, d, qq, eta=None):
    if eta:
        return np.sum(dx / d * force_gauss(d, qq, eta))
    else:
        return np.sum(dx / d * force_point(d, qq))


time_start = time.perf_counter()
a = 1  # nearest neighbor distance i.e. lattice constant / sqrt(2)
x_elec = [-2, 2]
x_elyt = [-1, 1]
q_elyt = [0.5, -0.5]
distance_plates = x_elec[1] - x_elec[0]  # distance between plates
v = np.array([-0.5, 0.5]) * (QE2F / COULOMB)

# distances to images within electrode and to opposite electrode
distances = a * np.linalg.norm(lattice(LENGTH), axis=1)
opposite_distances = np.sqrt(np.square(distances) + distance_plates**2)
image_distances = []
for x in x_elec:
    image_distances.append([])
    for y in x_elyt:
        image_distances[-1].append(np.sqrt(np.square(distances) + np.abs(y - x) ** 2))
image_elyt_distances = [[None for _ in range(len(x_elyt))] for _ in range(len(x_elyt))]
for i, (xi, qi) in enumerate(zip(x_elyt, q_elyt)):
    for j, (xj, qj) in list(enumerate(zip(x_elyt, q_elyt)))[i + 1 :]:
        image_elyt_distances[i][j] = np.sqrt(
            np.square(distances) + np.abs(xj - xi) ** 2
        )

for name, eta_elec in [("", [2.0, 2.0]), ("_eta_mix", [0.5, 3.0])]:
    # for name, eta_elec in [("", [2.0, 2.0])]:
    eta_mix = np.prod(eta_elec) / np.sqrt(np.sum(np.square(eta_elec)))
    # self interaction and within original box
    A_11 = np.sqrt(2 / np.pi) * eta_elec[0]
    A_22 = np.sqrt(2 / np.pi) * eta_elec[1]
    A_12 = erf(eta_mix * distance_plates) / distance_plates

    # interaction with periodic images
    A_11 += 4 * np.sum(a_element(distances, eta_elec[0] / SQRT2))
    A_22 += 4 * np.sum(a_element(distances, eta_elec[1] / SQRT2))
    A_12 += 4 * np.sum(a_element(opposite_distances, eta_mix))
    A = np.array([[A_11, A_12], [A_12, A_22]])
    inv = np.linalg.inv(A)
    e = np.array([1, 1])
    inv -= np.matmul(inv, np.matmul(np.outer(e, e), inv)) / np.dot(e, np.dot(inv, e))

    # electrode-electrolyte interaction
    b = []
    for i, (x, eta) in enumerate(zip(x_elec, eta_elec)):
        bi = 0
        for j, (y, q) in enumerate(zip(x_elyt, q_elyt)):
            bi += b_element(np.abs(y - x), q, eta)
            bi += 4 * np.sum(b_element(image_distances[i][j], q, eta))
        b.append(bi)
    b = np.array(b)

    # electrolyte-electrolyte energy
    elyt_11 = 4 * np.sum(1 / distances)
    distance_elyt = x_elyt[1] - x_elyt[0]
    elyt_12 = 1 / distance_elyt + 4 * np.sum(1 / image_elyt_distances[0][1])
    elyt = np.array([[elyt_11, elyt_12], [elyt_12, elyt_11]])
    energy_elyt = 0.5 * np.dot(q_elyt, np.dot(elyt, q_elyt))

    # electrode charges and energy
    q = np.dot(inv, v - b)
    energy = COULOMB * (0.5 * np.dot(q, np.dot(A, q)) + np.dot(b, q) + energy_elyt)

    # forces in out-of-plane direction
    f_elec = np.zeros(len(x_elec))
    f_elyt = np.zeros(len(x_elyt))
    # electrode-electrode
    dx = x_elec[1] - x_elec[0]
    fij_box = force_component(dx, np.abs(dx), q[0] * q[1], eta_mix)
    fij_img = 4 * force_component(dx, opposite_distances, q[0] * q[1], eta_mix)
    f_elec[0] -= fij_box + fij_img
    f_elec[1] += fij_box + fij_img
    # electrode-electrolyte
    for i, (xi, qi, etai) in enumerate(zip(x_elec, q, eta_elec)):
        for j, (xj, qj) in enumerate(zip(x_elyt, q_elyt)):
            dx = xj - xi
            fij_box = force_component(dx, np.abs(dx), qi * qj, etai)
            fij_img = 4 * force_component(dx, image_distances[i][j], qi * qj, etai)
            f_elec[i] -= fij_box + fij_img
            f_elyt[j] += fij_box + fij_img
    # electrolyte-electrolyte
    for i, (xi, qi) in enumerate(zip(x_elyt, q_elyt)):
        for j, (xj, qj) in list(enumerate(zip(x_elyt, q_elyt)))[i + 1 :]:
            dx = xj - xi
            fij_box = force_component(dx, np.abs(dx), qi * qj)
            fij_img = 4 * force_component(dx, image_elyt_distances[i][j], qi * qj)
            f_elyt[i] -= fij_img + fij_box
            f_elyt[j] += fij_img + fij_box
    # force units
    assert np.abs(np.sum(f_elec) + np.sum(f_elyt)) < 1e-8
    f_elec *= COULOMB
    f_elyt *= COULOMB

    # Virial
    volume = a**2 * LZ
    virial = 0.0
    for x, f in [(x_elec, f_elec), (x_elyt, f_elyt)]:
        virial += np.dot(x, f)
    pressure = NKTV2P * virial / volume

    with open(f"plate_cap{name}.csv", "w") as f:
        f.write(
            "length, energy / kcal/mol, q1 / e, q2 / e, inv11 / A, inv12 / A"
            + ", b1 / e/A, b2 / e/A, felec1 / kcal/mol/A, felec2 / kcal/mol/A"
            + ", felyt1 / kcal/mol/A, felyt2 / kcal/mol/A, press\n"
        )
        f.write(
            ", ".join(
                [
                    str(LENGTH),
                    f"{energy:.8f}",
                    f"{q[0]:.10f}",
                    f"{q[1]:.10f}",
                    f"{inv[0, 0]:.10f}",
                    f"{inv[0, 1]:.10f}",
                    f"{b[0]:.8f}",
                    f"{b[1]:.8f}",
                    f"{f_elec[0]:.5f}",
                    f"{f_elec[1]:.5f}",
                    f"{f_elyt[0]:.5f}",
                    f"{f_elyt[1]:.5f}",
                    f"{pressure:.2f}",
                ]
            )
            + "\n"
        )
time_end = time.perf_counter()
print(f"{time_end - time_start:0.4f} seconds")
