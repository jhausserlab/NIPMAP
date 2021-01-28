from src.CellAbundance import CellAbundance
def site_generation_and_count(x, y, t, x_span, y_span, radius, grid_point):        
    site = []
    x_center = (x_span[grid_point[0]+1] + x_span[grid_point[0]]) / 2
    y_center = (y_span[grid_point[1]+1] + y_span[grid_point[1]]) / 2
    site = np.array([(x[c_idx], y[c_idx], t[c_idx]) for c_idx in range(len(x)) if CellAbundance.is_in_cirle(x[c_idx], y[c_idx], x_center, y_center, radius)])
    return site, grid_point

CELL_TYPES = ['Kreatin-positive tumor', 'Treg', 'CD3-T', 'Neutrophils', 'Tumor', 'B', 
              'Macrophages', 'Mesenchymal-like', 'Other Immune', 'CD8-T', 'CD4-T', 
              'Undefined', 'Mono/Neu', 'DC/Mono', 'Endothelial', 'DC', 'NK']
x = np.random.randint(800, size=1000)
y = np.random.randint(800, size=1000)
x_span = np.arange(0, 825, 25)
y_span = np.arange(0, 825, 25)
t = np.random.choice(CELL_TYPES, 1000)
radius = 100
mesh = np.array(np.meshgrid(range(len(x_span)-1), range(len(y_span)-1)))
grid_point = list(mesh.T.reshape(-1, 2))[200]
s, gp = site_generation_and_count(x, y, t, x_span, y_span, radius, grid_point)
