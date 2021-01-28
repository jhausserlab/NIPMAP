import numpy as np


def tumor_misfit_error(original_sites, approx_sites):
    denom = np.sum((original_sites - np.mean(original_sites, axis=0))**2)
    print(denom)
    error = np.sum((np.subtract(approx_sites, original_sites))**2) / denom
    return error


def arch2color(arch):
    d = {0: 'red', 1: 'blue', 2: 'green'}
    return d[arch]


def alfa2rgb(v):
    rgb_vect = (v * 255).astype('uint8')
    return rgb_vect


def alfa2color(v):
    d = {0: [255, 0, 0, 0], 1: [0, 255, 0, 0], 2: [0, 0, 255, 0], 3: [0, 0, 0, 255]}
    return np.array(d[np.argmax(v)])


def color_mapper(color_vector, alfa):
    v = np.array([alfa[i] * color_vector[:, i] for i in range(len(alfa))]).T
    rgb_code = np.sum(v, axis=1)
    return rgb_code


def scale(x, min, max, new_min=0, new_max=1):
    return (new_max-new_min) / (max - min) * (x - max) + max
