import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np

colors = ['#629563', '#044E75', '#CA8F04', '#645D0D',
          '#43BC52', '#B25E89', '#2E3790', '#F118BE',
          '#50974E', '#3273D6', '#0AF24B', '#A3F159',
          '#933835', '#CEB134', '#226BCF', '#856218']


def plot_cells_positions(data, to_plot='all'):
    groups = data.groupby('cell_type')

    plt.figure(figsize=(10, 10))
    cm = plt.cm.get_cmap('tab10')
    print(len(data['cell_type'].unique()))
    for (name, group), col in zip(groups, colors):
        if to_plot == 'all' or name in to_plot:
            plt.scatter(group['x'], group['y'], marker="o", s=5, label=name, c=col)
    plt.legend()
    plt.show()



