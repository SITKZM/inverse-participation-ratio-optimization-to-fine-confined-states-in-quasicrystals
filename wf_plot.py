import functions

coordinates = functions.read_matrix_from_file("coordinates.txt")
hoppings = functions.read_matrix_from_file("hop.txt")

N = len(coordinates)
wfnums = [65 + i for i in range(6)]
for wfnum in wfnums:
    wf = functions.read_matrix_from_file(f"wave_functions/WF_{wfnum}.txt")
    wf = wf.reshape(-1)
    output = f"figures/AB{N}_WF_{wfnum}.svg"

    vmax = abs(wf).max()
    wf = wf / vmax
    functions.make_plot_WF(coordinates, hoppings, wf, -1, 1, "seismic", 100, output) #s: AB33; 500, AB153; 100 