def group_map(id, tumor_yn):
    labels = {0: [
        'Tumor', 'Treg', 'CD4-T', 'CD8-T', 'CD3-T', 'NK', 'B', 'Neutrophils',
        'Macrophages', 'DC', 'DC/Mono', 'Mono/Neu', 'Other Immune'],
        1: ['Undefined', 'Immune', 'Endothelial', 'Mesenchymal-like',
                  'Tumor', 'Kreatin-positive tumor']}
    id = id-1 if tumor_yn == 1 else id
    return labels[tumor_yn][id]


class Cell:
    def __init__(self, x, y, cell_label, cell_type):
        self.x = x
        self.y = y
        self.label = cell_label
        self.cell_type = cell_type
