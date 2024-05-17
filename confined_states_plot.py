import numpy as np
import functions

coordinates = functions.read_matrix_from_file("coordinates.txt")
hoppings = functions.read_matrix_from_file("hop.txt")

N = len(coordinates)
N_c = 2
wfnums = [13 + i for i in range(N_c)]
wfs = np.zeros((N, N_c))
psis = np.zeros((N_c, N_c))
for i in range(N_c):
    wf = functions.read_matrix_from_file(f"wave_functions/AB{N}/WF_{wfnums[i]}.txt")
    wf = wf.reshape(-1)
    wfs[:, i] = wf

    psi = functions.read_matrix_from_file(f"confined_states/AB{N}/psis_{i + 1}.txt")
    psi = psi.reshape(-1)
    psis[:, i] = psi

for i in range(N_c):
    wf = np.zeros(N)
    for j in range(N_c):
        wf += psis[j, i] * wfs[:, j]
    vmax = abs(wf).max()
    wf = wf / vmax
    output = f"figures/AB{N}_confined_state_{i + 1}.svg"
    functions.make_plot_WF(coordinates, hoppings, wf, -1, 1, "seismic", 500, output) #s: AB33; 500, AB153; 100 