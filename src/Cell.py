# TODO: to delete?
#
# def group_map2(id, tumor_yn):
#     labels = {0: [
#         'Tumor', 'Treg', 'CD4-T', 'CD8-T', 'CD3-T', 'NK', 'B', 'Neutrophils',
#         'Macrophages', 'DC', 'DC/Mono', 'Mono/Neu', 'Other Immune'],
#         1: ['Undefined', 'Immune', 'Endothelial', 'Mesenchymal-like',
#                   'Tumor', 'Kreatin-positive tumor']}
#     id = id-1 if tumor_yn == 1 else id
#     return labels[tumor_yn][id]


def group_map(cell_id, immune_id):
    cells = {1: 'Unidentified', 2: 'Immune', 3: 'Endothelial', 4:'Mesenchymal-like',
             5: 'Tumor', 6: 'Kreatin-positive tumor'}
    immune_cells = {1: 'Tregs', 2: 'CD4-T', 3: 'CD8-T', 4: 'CD3-T', 5: 'NK',
                    6: 'B', 7: 'Neutrophils', 8: 'Macrophages', 9: 'DC',
                    10: 'DC / Mono', 11: 'Mono / Neu', 12: 'Other immune'}

    return cells[cell_id] if cell_id != 2 else immune_cells[immune_id]


class Cell:
    def __init__(self, x, y, cell_label, cell_type):
        self.x = x
        self.y = y
        self.label = cell_label
        self.cell_type = cell_type
