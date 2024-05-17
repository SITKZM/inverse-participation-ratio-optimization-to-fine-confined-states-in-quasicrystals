import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.size"] = 36
plt.rcParams["svg.fonttype"] = "none"

def read_matrix_from_file(filename):
    # ファイルを開いて行ごとに読み込み、各行をリストに格納する
    with open(filename, 'r') as file:
        lines = file.readlines()

    # 各行をスペースで分割し、数値部分をfloatに変換する
    matrix_data = []
    for line in lines:
        row = [float(num) for num in line.split()]
        matrix_data.append(row)

    # Numpyのarrayに変換して返す
    return np.array(matrix_data)

def make_plot_WF(coordinates, hoppings, data, data_min, data_max, cmaptype, size, name):
    fig, ax = plt.subplots(figsize=(13,10))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)

    for hopping_index in hoppings:
        i = int(hopping_index[0])
        j = int(hopping_index[1])
        ax.plot([coordinates[i, 0], coordinates[j, 0]],
                [coordinates[i, 1], coordinates[j, 1]], c='gray')

    N = len(coordinates)
    x = [coordinates[i, 0] for i in range(N)]
    y = [coordinates[i, 1] for i in range(N)]
    c = [data[i] for i in range(N)]
    mappable = ax.scatter(x, y, c=c, vmin=data_min, vmax=data_max, s=size, cmap=cmaptype, zorder=3) # sの値: AB478=100, AB1393=50, AB2786=30, AB8119=10

    cbar = fig.colorbar(mappable, ax=ax)
    cbar.formatter.set_useMathText(True)
    cbar.formatter.set_powerlimits((0, 0))

    fig.savefig(name)