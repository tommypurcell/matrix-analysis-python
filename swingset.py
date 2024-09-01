import numpy as np
from plot import plotstructure
from analysis import estiff, etran

# Configuration
np.set_printoptions(precision=4, suppress=True, linewidth=150)

# Input Information
nnodes = 7
ndof = nnodes * 6

# Node coordinates in mm
coord = 1000 * np.array([[-1, 0, 0], [1, 0, 0], [0, 0, 2.5], [0, 1.5, 2.5], [0, 3, 2.5], [-1, 3, 0], [1, 3, 0]])

# External nodal forces
concen = np.array([
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, -5.5, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0]
])

# Fixity conditions (0 = fixed, NaN = free)
fixity = np.full((nnodes, 6), np.nan)
fixity[0, :3] = 0  # Node 1 fixed in UX, UY, UZ
fixity[1, :3] = 0  # Node 2 fixed in UX, UY, UZ
fixity[5, :3] = 0  # Node 6 fixed in UX, UY, UZ
fixity[6, :3] = 0  # Node 7 fixed in UX, UY, UZ

# Number of elements
nele = 6

# End nodes of each element
ends = np.array([[1, 3], [2, 3], [3, 4], [4, 5], [6, 5], [7, 5]])

# Member properties
to_vector = np.ones(nele)
A = 1430 * to_vector
I = 1.26e6 * to_vector
J = 2.52e6 * to_vector
E = 200 * to_vector
v = 0.3 * to_vector
beta_ang = np.zeros(nele)

# DOFs for each node
node_id = np.zeros((nnodes, 6), dtype=int)
for i in range(nnodes):
    node_id[i, :] = i * 6 + np.arange(6)

# DOFs for each member
mem_id = np.zeros((nele, 12), dtype=int)
for i in range(nele):
    mem_id[i, :6] = node_id[ends[i, 0] - 1, :]
    mem_id[i, 6:] = node_id[ends[i, 1] - 1, :]

# Reshape fixity matrix and find fixed/free DOFs
D = fixity.T.flatten(order='F')
fixed_dof = np.where(D == 0)[0]
free_dof = np.where(np.isnan(D))[0]

# External forces vector
P_total = concen.T.flatten(order='F')

# Calculate member lengths and xprime axis coordinates
L = np.zeros(nele)
x_axis = np.zeros((nele, 3))
for i in range(nele):
    x1, y1, z1 = coord[ends[i, 0] - 1]
    x2, y2, z2 = coord[ends[i, 1] - 1]
    L[i] = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
    x_axis[i, :] = (x2 - x1, y2 - y1, z2 - z1) / L[i]

# Create local stiffness matrices
k_stack_local = np.zeros((12, 12, nele))
for i in range(nele):
    k_stack_local[:, :, i] = estiff(A[i], I[i], I[i], J[i], E[i], v[i], L[i])

# Transformation matrices
k_trans_stack = np.zeros((12, 12, nele))
for i in range(nele):
    k_trans_stack[:, :, i] = etran(beta_ang[i], x_axis[i, :])

# Global stiffness matrices
k_stack_global = np.zeros((12, 12, nele))
for i in range(nele):
    k_stack_global[:, :, i] = k_trans_stack[:, :, i].T @ k_stack_local[:, :, i] @ k_trans_stack[:, :, i]

# Assemble global stiffness matrix
k_total = np.zeros((ndof, ndof))
for i in range(nele):
    k_total[np.ix_(mem_id[i, :], mem_id[i, :])] += k_stack_global[:, :, i]

# Partition global stiffness matrix
kff = k_total[np.ix_(free_dof, free_dof)]
ksf = k_total[np.ix_(fixed_dof, free_dof)]

# Solve for free displacements
p_free = P_total[free_dof]
d_free = np.linalg.solve(kff, p_free)

# Insert free displacements into total displacements vector
d_total = np.zeros(ndof)
d_total[free_dof] = d_free

# Reshape displacements into node format
DEFL = d_total.reshape((nnodes, 6))

# Calculate reaction forces
react = ksf @ d_free

# Update total forces with reactions
P_total[fixed_dof] = react

# Calculate internal forces in each member
ELE_FOR = np.zeros((nele, 12))
for i in range(nele):
    d_global = d_total[mem_id[i, :]]
    d_local = k_trans_stack[:, :, i] @ d_global
    f_local = k_stack_local[:, :, i] @ d_local
    ELE_FOR[i, :] = f_local

# Print Reaction Forces for Base Nodes
print('\nReaction Forces for Base Nodes:')
for i, node in enumerate([1, 2, 6, 7]):
    print(f'Node {node}\tFX = {react[3*i]:.4f} kN\tFY = {react[3*i+1]:.4f} kN\tFZ = {react[3*i+2]:.4f} kN')

# Print Node Deflections with better alignment
print('Node Deflections:')
print('Node\t  UX (mm)   \t  UY (mm)   \t  UZ (mm)   \t  RX (rad)  \t  RY (rad)  \t  RZ (rad)')
for i in range(len(DEFL)):
    print(f'{i+1:<4}\t{DEFL[i, 0]:>10.4f}\t{DEFL[i, 1]:>10.4f}\t{DEFL[i, 2]:>10.4f}\t{DEFL[i, 3]:>10.4f}\t{DEFL[i, 4]:>10.4f}\t{DEFL[i, 5]:>10.4f}')

# Plot the structure
plotstructure(coord, ends)
