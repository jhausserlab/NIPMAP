import sys
sys.path.insert(1, '/srv/mfs/hausserlab/fabio/data_analysis/src/utils')
from archetypes import ArchetypalAnalysis
import pandas as pd
import numpy as np
import json

data1 = json.loads(sys.argv[1])
evs = json.loads(sys.argv[2])
means = json.loads(sys.argv[3])
nbArchs = sys.argv[4]
outputPath1 = sys.argv[5]
## FIXME:A print command in the archetypal analysis ==> remove it !! 
# Reading the file of PCA from Wagner dataset
#pcData = pd.read_csv("/scratch/anissa.el/ImmuneStates/wagner_analysis/data/3PCA_Wagner.csv")
pc3dWagner =np.array(data1, dtype = "float64") #np.array(pcData, dtype = "float64")
AA = ArchetypalAnalysis(n_archetypes = nbArchs, 
                          tolerance = 0.001, 
                          max_iter = 200, 
                          random_state = 0, 
                          C = 0.0001, 
                          initialize = 'random',
                          redundancy_try = 30)
AA.fit_transform(pc3dWagner)
#archetypes_wagner = AA.alfa
archsCellsSpace = np.dot(AA.archetypes.T,evs) + np.array(means,dtype="float64") 
pd.DataFrame(AA.archetypes).to_csv(outputPath1+"shuffled_data_archetypes.csv")
#pd.DataFrame(archetypes_wagner).to_csv("/scratch/anissa.el/ImmuneStates/wagner_analysis/outputs/Archetypes_PCA_Wagner.csv")
