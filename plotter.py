import numpy as np
import matplotlib.pyplot as plt


def plot_classical(spin_lattice, plot_size=(6,3)):
    theta = spin_lattice[0]
    phi = spin_lattice[1]

    Sx = np.sin(theta) * np.cos(phi)  # x-component of the spin
    Sy = np.sin(theta) * np.sin(phi)  # y-component of the spin
    Sz = np.cos(theta)               # z-component of the spin

    fig, axes = plt.subplots(1, 3, figsize=(3 * plot_size[0], plot_size[1]))
    components = [(Sx, "Sx (Spin x-component)"), 
                  (Sy, "Sy (Spin y-component)"), 
                  (Sz, "Sz (Spin z-component)")]

    for ax, (component, title) in zip(axes, components):
        im = ax.imshow(component, cmap='gray', origin='upper', extent=[0, theta.shape[1], 0, theta.shape[0]])
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("x (Lattice)")
        ax.set_ylabel("y (Lattice)")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    plt.tight_layout()
    plt.show()

############
def plot_magnetization(Temp, mmvT_list, f, plot_type):
    plt.figure(figsize=(6, 4))  # Set the figure size
    plt.plot(Temp, mmvT_list, color='blue', linestyle='-', linewidth=2, marker='o', markersize=4, markerfacecolor='red')
    plt.xlabel('Temperature (T)', fontsize=12)
    if plot_type == 'Classical': plt.ylabel('Localized Spin Magnetization (m)', fontsize=12)
    else: plt.ylabel('Magnetization (m)', fontsize=12)
    # plt.title(f'{plot_type} Magnetization vs Temperature (n={f:.2f})', fontsize=14, fontweight='bold')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()


#########################3
def calculate_dq(spin_lattice, q_vector, n):
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
                    r_ij = np.array([ni - i, nj - j]) 
                    dot_q_r = np.dot(q_vector, r_ij)
                    Si_dot_Sj = ( Sx[i, j] * Sx[ni, nj] + Sy[i, j] * Sy[ni, nj] + Sz[i, j] * Sz[ni, nj])
                    D_q += Si_dot_Sj * np.exp(1j * dot_q_r)
    
    return D_q


############################
def plot_dq_vs_q(spin_lattice, n, trans = np.abs):
    theta = spin_lattice[0]
    phi = spin_lattice[1]

    q_vals = np.arange(np.min(np.fft.fftfreq(n, d=1/n) * 2 * np.pi / n), np.max(np.fft.fftfreq(n, d=1/n) * 2 * np.pi / n), 0.1)
    # print(np.fft.fftfreq(n, d=1/n), q_vals)
    q_vectors = [(qx, qy) for qx in q_vals for qy in q_vals]
    # print(len(q_vectors))
    
    dq_values = []
    q_magnitudes = []
    
    for i, (qx, qy) in enumerate(q_vectors):
        q_vector = np.array([qx, qy])
        dq = calculate_dq(spin_lattice, q_vector, n)
        q_magnitudes.append(np.linalg.norm(q_vector))
        # print(i, dq, qx, qy)
        dq_values.append(trans(dq))  
 
    q_vectors = np.array(q_vectors)  
    dq_values = np.array(dq_values)
    q_x = q_vectors[:, 0]
    q_y = q_vectors[:, 1]

    max_idx = np.argmax(dq_values)
    print(f"Max Dq: {dq_values[max_idx]} at q = {q_vectors[max_idx]}, {max_idx}")

    qx_unique = np.unique(q_x)
    qy_unique = np.unique(q_y)

    # Create a 2D array for Dq
    Dq_grid = np.zeros((len(qx_unique), len(qy_unique)))

    for i, q in enumerate(q_vectors):
        qx_idx = np.where(qx_unique == q[0])[0][0]
        qy_idx = np.where(qy_unique == q[1])[0][0]
        Dq_grid[qx_idx, qy_idx] = dq_values[i]

    # Plot the heatmap
    plt.figure(figsize=(8, 6))
    plt.imshow(Dq_grid, extent=[qx_unique.min(), qx_unique.max(), qy_unique.min(), qy_unique.max()],
            origin='lower', aspect='auto', cmap='viridis')
    plt.colorbar(label='D(q)')
    plt.xlabel('q_x')
    plt.ylabel('q_y')
    # plt.title('Heatmap of D(q) vs q_x, q_y')
    plt.show()