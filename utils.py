import numpy as np
import math
import matplotlib.pyplot as plt
from tqdm import tqdm, trange
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh


########### type 1: random phi(s), type 0: phi =0
'''Generate a random spin lattice with n x n sites.'''
def generate_lattice(n, seed=1, type = 1):
    np.random.seed(seed)  
    u = np.random.uniform(-1, 1, size=(n, n))  
    if type: phi = np.random.uniform(0, 2 * np.pi, size=(n, n)) 
    else: phi = np.zeros((n, n))
    theta = np.arccos(u)  
    spin_lattice = np.stack([theta, phi], axis=0)
    return spin_lattice


#############
'''Calculates the Hamiltonian matrix for the given spin lattice.'''
def calculate_hamiltonian(spin_lattice, n, t_0=1.0):
    num_sites = n * n
    H = np.zeros((num_sites, num_sites))
    
    theta = spin_lattice[0, :, :]
    phi = spin_lattice[1, :, :]

    def H_element(i, j, ni, nj):
        current_index = i * n + j
        neighbor_index = ni * n + nj
        
        cos_theta_ij = np.sqrt(0.5 * (1 + np.sin(theta[i, j]) * np.sin(theta[ni, nj]) 
        * np.cos(phi[i, j] - phi[ni, nj]) + np.cos(theta[i, j]) * np.cos(theta[ni, nj])))

        t_ij = t_0 * cos_theta_ij
        
        H[current_index, neighbor_index] = -t_ij
        H[neighbor_index, current_index] = -t_ij 

    for i in range(n):
        for j in range(n):
            H_element(i, j, (i + 1) % n, j)
            H_element(i, j, i, (j + 1) % n) # Right
    
    return H


################
'''Computes the eigenvalues and eigenvectors of the Hamiltonian matrix.'''
def compute_eigenvalues(H, k=None):
    H_sparse = csr_matrix(H)
    if k is None or k >= H_sparse.shape[0]:
        eigenvalues, eigenvectors = np.linalg.eigh(H_sparse.toarray())
    else:
        eigenvalues, eigenvectors = eigsh(H_sparse, k=k, which='SA')
    
    sorted_indices = np.argsort(eigenvalues)
    sorted_eigenvalues = eigenvalues[sorted_indices]
    sorted_eigenvectors = eigenvectors[:, sorted_indices]
    
    return sorted_eigenvalues, sorted_eigenvectors


##############
'''Metropolis algorithm for the double exchange model.'''
def metropolis(spin_lattice, n, f, b, t_0=1.0, iter=1000, type = 1):

    H = calculate_hamiltonian(spin_lattice, n, t_0)
    E_ini, _ = compute_eigenvalues(H)
    E_ini = E_ini[:f].sum()
    E_list = []

    for i in trange(iter):
        theta, phi = spin_lattice[0], spin_lattice[1]
        x, y = np.random.randint(0, n, size=2)
        if i%20 == 0:  E_list.append(E_ini)

        # delta_theta = np.random.uniform(-epsilon, epsilon)
        # if type: delta_phi = np.random.uniform(-epsilon, epsilon)
        # else: delta_phi = 0
        # theta_new = (theta[x, y] + delta_theta) % np.pi
        # phi_new = (phi[x, y] + delta_phi) % (2 * np.pi)

        u = np.random.uniform(-1, 1)
        if type: phi_new = np.random.uniform(0, 2 * np.pi)
        else: phi_new = 0
        theta_new = np.arccos(u)

        spin_lattice_new = spin_lattice.copy()
        spin_lattice_new[0, x, y], spin_lattice_new[1, x, y] = theta_new, phi_new    

        H_new = calculate_hamiltonian(spin_lattice_new, n, t_0)
        E_new, _ = compute_eigenvalues(H_new)
        E_new = E_new[:f].sum()

        if E_new < E_ini or np.random.rand() < np.exp(-b * (E_new - E_ini)):
            # if i%(iter/100) == 0 : print(np.exp(-b * (E_new - E_ini)))
            spin_lattice = spin_lattice_new
            E_ini = E_new

    return spin_lattice, E_list


####################
'''Calculates the magnetic moment of the localized spins.'''
def calculate_mag(spin_lattice):

    theta = spin_lattice[0]  # Polar angle θ
    phi = spin_lattice[1]    # Azimuthal angle φ

    mx = np.sum(np.sin(theta) * np.cos(phi))
    my = np.sum(np.sin(theta) * np.sin(phi))
    mz = np.sum(np.cos(theta))

    return math.sqrt(mx**2 + my**2 + mz**2)


####################
'''Calculates the magnetic moment of the hopping electrons.'''
def quantum_mag(spin_lattice, E_vec):
    theta, phi = spin_lattice[0], spin_lattice[1]
    E_vecr = E_vec.reshape(theta.shape)
    # print(E_vecr)
    mx = np.mean(np.sin(theta) * np.cos(phi) * E_vecr)
    my = np.mean(np.sin(theta) * np.sin(phi) * E_vecr)
    mz = np.mean(np.cos(theta) * E_vecr)
    return math.sqrt(mx**2 + my**2 + mz**2)

# def quantum_mag(spin_lattice, E_vec):
#     n = spin_lattice.shape[1]
#     theta = spin_lattice[0]
#     return max(f * np.sum(np.sin(theta)) / (2*n*n), f * np.sum(np.cos(theta)) / (2*n*n))


###############################
'''Calculates the Spin-Spin correlation parameter.'''
def calculate_dq0(spin_lattice, n):
    theta = spin_lattice[0]
    phi = spin_lattice[1]

    Sx = np.sin(theta) * np.cos(phi)
    Sy = np.sin(theta) * np.sin(phi)
    Sz = np.cos(theta)
    
    D_q = 0.0

    for i in range(n):
        for j in range(n):
            for ni in range(n):
                for nj in range(n):
                    Si_dot_Sj = ( Sx[i, j] * Sx[ni, nj] + Sy[i, j] * Sy[ni, nj] + Sz[i, j] * Sz[ni, nj])
                    D_q += Si_dot_Sj
    
    return float(D_q/(n*n))


epsilon = 1.0 #constant