import re
import numpy as np
import xlrd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


name = "Hybri_Xla.stage9.go" # name of data file
types = "BP" # KEGG BP MF CC
fs = [10, 3] # size of the picture
grid = [1, 4] # layout

config = {
    'figure.dpi': 300,
    'axes.linewidth': 0.5,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'axes.unicode_minus': False,
    "font.family": 'serif',
    "font.serif": ['Arial'],
    "font.size": 6,
}


def get_rgb(s, mi, ma, cm):
    if ma != mi:
        return cm[int(255 * (s - mi) / (ma - mi))]
    else:
        return cm[0]


with PdfPages(f'{name}_{types}.pdf') as pdf:
    plt.rcParams.update(config)
    fig = plt.figure(figsize=fs)
    gs0 = gridspec.GridSpec(1, 1, figure=fig)
    gs00 = gridspec.GridSpecFromSubplotSpec(*grid, subplot_spec=gs0[0], wspace=0.5)
    ax = []
    scatter = []
    inp = dict()
    data = xlrd.open_workbook(name + '.xlsx')
    for name in data._sheet_names:
        result = []
        if True:
            inp[name] = data.sheet_by_name(name)
            if types != "KEGG":
                for n in range(1, len(inp[name].col_values(0))):
                    line = inp[name].row_values(n)
                    if line[0] == types:
                        result.append([int(line[3].split("/")[0])/int(line[3].split("/")[1]), line[2], float(line[9]), float(line[6])])
                if len(result) > 0:
                    print("ok")
                    array = sorted(result, key=lambda x: x[0])[-20:]
                    ax.append(fig.add_subplot(gs00[0, 1 + 2 * len(ax)]))
                    scatter.append(ax[-1].scatter([0 for x in range(0, len(array))], [0 for x in range(0, len(array))],
                                                  s=[x[2] * 0.5 for x in array], c=[x[3] for x in array],
                                                  vmin=min([x[3] for x in array]), vmax=max([x[3] for x in array]),
                                                  zorder=0))
                    ax[-1].grid(linestyle="--", linewidth=0.4, zorder=0)
                    ax[-1].tick_params(width=0.5, zorder=4)
                    for i, dot in enumerate(array):
                        ax[-1].scatter(dot[0], i, s=dot[2] * 0.5, c=dot[3], vmin=min([x[3] for x in array]),
                                       vmax=max([x[3] for x in array]), zorder=10)

                    ax[-1].legend(scatter[-1].legend_elements(prop='sizes', num=3)[0],
                                  [str(2 * int(re.search("[0-9]+", x).group(0))) for x in
                                   scatter[-1].legend_elements(prop='sizes', num=3)[1]], loc='lower right',
                                  title='Count')

                    plt.colorbar(scatter[-1], label='P.adjust', ax=ax[-1])
                    ax[-1].set_xlim(
                        min([x[0] for x in array]) - (max([x[0] for x in array]) - min([x[0] for x in array])) * 0.1,
                        max([x[0] for x in array]) + (max([x[0] for x in array]) - min([x[0] for x in array])) * 0.1)
                    ax[-1].set_yticks(range(0, len(array)), labels=[x[1] for x in array])
                    ax[-1].set_ylim(-1, len(array))
                    ax[-1].yaxis.tick_left()
                    ax[-1].set_xlabel("Gene Ratio")
                    ax[-1].set_title(" ".join(f'{name} {types}'.split("_")), fontsize=6)
            elif types == "KEGG":
                for n in range(1, len(inp[name].col_values(0))):
                    line = inp[name].row_values(n)
                    result.append([int(line[2].split("/")[0]) / int(line[2].split("/")[1]), line[1], float(line[8]), float(line[5])])
                if len(result) > 0:
                    print("ok")
                    array = sorted(result, reverse=True, key=lambda x: x[3])[-20:]
                    ax.append(fig.add_subplot(gs00[0, 1 + 2 * len(ax)]))
                    ax[-1].grid(linestyle="--", linewidth=0.4, zorder=0)
                    ax[-1].tick_params(width=0.5, zorder=4)
                    scatter.append(ax[-1].scatter([0 for x in range(0, len(array))], [0 for x in range(0, len(array))],
                                           s=[0 for x in array], c=[x[3] for x in array],
                                           vmin=min([x[3] for x in array]), vmax=max([x[3] for x in array]), zorder=0))
                    ax[-1].barh(np.arange(len([x[1] for x in array])), [x[2] for x in array],
                                    color=[get_rgb(x[3], min([y[3] for y in array]), max([y[3] for y in array]),
                                                   scatter[-1].cmap.colors)
                                           for x in array], align='center', zorder=2)
                    plt.colorbar(scatter[-1], label='P.adjust', ax=ax[-1])
                    ax[-1].set_yticks(range(0, len(array)), labels=[x[1] for x in array])
                    ax[-1].set_ylim(-1, len(array))
                    ax[-1].yaxis.tick_left()
                    ax[-1].set_xlabel("Gene Count")
                    ax[-1].set_title(" ".join(f'{name} {types}'.split("_")), fontsize=6)
    pdf.savefig(fig)
    plt.close()
