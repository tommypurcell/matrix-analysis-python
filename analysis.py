import numpy as np

def estiff(A, Izz, Iyy, J, E, v, L):
    # Rows of the element stiffness matrix
    r1 = [A / L, 0, 0, 0, 0, 0, -A / L, 0, 0, 0, 0, 0]
    r2 = [0, (12 * Izz) / (L**3), 0, 0, 0, (6 * Izz) / (L**2), 0, (-12 * Izz) / (L**3), 0, 0, 0, (6 * Izz) / (L**2)]
    r3 = [0, 0, (12 * Iyy) / (L**3), 0, (-6 * Iyy) / (L**2), 0, 0, 0, (-12 * Iyy) / (L**3), 0, (-6 * Iyy) / (L**2), 0]
    r4 = [0, 0, 0, J / (2 * (1 + v) * L), 0, 0, 0, 0, 0, -J / (2 * (1 + v) * L), 0, 0]
    r5 = [0, 0, -(6 * Iyy) / (L**2), 0, (4 * Iyy) / L, 0, 0, 0, (6 * Iyy) / (L**2), 0, (2 * Iyy) / L, 0]
    r6 = [0, (6 * Izz) / (L**2), 0, 0, 0, (4 * Izz) / L, 0, -(6 * Izz) / (L**2), 0, 0, 0, (2 * Izz) / L]
    r7 = [-A / L, 0, 0, 0, 0, 0, A / L, 0, 0, 0, 0, 0]
    r8 = [0, -(12 * Izz) / (L**3), 0, 0, 0, -(6 * Izz) / (L**2), 0, (12 * Izz) / (L**3), 0, 0, 0, -(6 * Izz) / (L**2)]
    r9 = [0, 0, -(12 * Iyy) / (L**3), 0, (6 * Iyy) / (L**2), 0, 0, 0, (12 * Iyy) / (L**3), 0, (6 * Iyy) / (L**2), 0]
    r10 = [0, 0, 0, -J / (2 * (1 + v) * L), 0, 0, 0, 0, 0, J / (2 * (1 + v) * L), 0, 0]
    r11 = [0, 0, -(6 * Iyy) / (L**2), 0, (2 * Iyy) / L, 0, 0, 0, (6 * Iyy) / (L**2), 0, (4 * Iyy) / L, 0]
    r12 = [0, (6 * Izz) / (L**2), 0, 0, 0, (2 * Izz) / L, 0, -(6 * Izz) / (L**2), 0, 0, 0, (4 * Izz) / L]

    # Element stiffness matrix
    elk = E * np.array([r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12])

    return elk


def etran(beta_ang, xaxis):
    # Normalize xaxis
    xaxis = np.array(xaxis) / np.linalg.norm(xaxis)

    # Determine ey based on the value of xaxis
    if np.array_equal(xaxis, [0, 1, 0]):
        ey = np.array([-1, 0, 0])
    elif np.array_equal(xaxis, [0, -1, 0]):
        ey = np.array([1, 0, 0])
    else:
        ey = np.array([0, 1, 0])

    # Calculate ez_double_prime
    ez_double_prime = np.cross(xaxis, ey) / np.linalg.norm(np.cross(xaxis, ey))

    # Calculate ey_double_prime
    ey_double_prime = np.cross(ez_double_prime, xaxis)

    # Rotation matrix R1
    R1 = np.array([xaxis, ey_double_prime, ez_double_prime])

    # Rotation matrix R2 for beta_ang
    R2 = np.array([
        [1, 0, 0],
        [0, np.cos(np.radians(beta_ang)), np.sin(np.radians(beta_ang))],
        [0, -np.sin(np.radians(beta_ang)), np.cos(np.radians(beta_ang))]
    ])

    # Small gamma matrix
    small_gamma = np.dot(R2, R1)

    # Zero matrix
    zero_matrix = np.zeros((3, 3))

    # Full gamma matrix
    gamma = np.block([
        [small_gamma, zero_matrix, zero_matrix, zero_matrix],
        [zero_matrix, small_gamma, zero_matrix, zero_matrix],
        [zero_matrix, zero_matrix, small_gamma, zero_matrix],
        [zero_matrix, zero_matrix, zero_matrix, small_gamma]
    ])

    return gamma
