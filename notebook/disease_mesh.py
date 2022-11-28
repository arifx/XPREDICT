import numpy as np
import itertools
import math
import numpy as np
import pandas as pd
import os
import shutil

disease_mesh = 'data/input/omim-disease-mesh.csv'
mesh_ann = {}
allmeshterm = []
with open(disease_mesh) as meshfile:
    next(meshfile)
    for line in meshfile:
        line = line.strip().split(',')
        if len(line) != 2: continue
        di = line[0]
        mesh = line[1].split(',')

        #### new code for adding an array of values
        try:
          if mesh_ann[di]:
            a = mesh_ann[di]
            a.append(mesh)
            mesh_ann[di] = a
        except:
            mesh_ann[di] = [mesh]

        #mesh_ann[di] = mesh
        allmeshterm.extend(mesh)

vocabulary = list(set(allmeshterm))
len(vocabulary)

co_mat = np.zeros((len(mesh_ann),len(vocabulary)))


commonDiseases = mesh_ann.keys()
mesh2id= { di:i for i,di in enumerate(mesh_ann.keys())}
# fill in the co-occurrence matrix
for key in mesh_ann:
    annotations = mesh_ann[key]
    #### new code to work with an array of values
    for annotation_array in annotations:
        for item in annotation_array:
            print(item)
            col_index = vocabulary.index(item)
            co_mat[mesh2id[key],col_index] = 1

def cosine_similarity(a,b):
    return  np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b))

values = []
# calculate cosine similarity between diseases using mesh annotation vector
for comb in itertools.combinations(commonDiseases,2) :
    disease1 = comb[0]
    disease2 = comb[1]
    sim = cosine_similarity(co_mat[mesh2id[disease1],:], co_mat[mesh2id[disease2],:])
    values.append([disease1, disease2, sim])

disease_pheno_df = pd.DataFrame(values, columns =['Disease1','Disease2','PHENO-SIM'])
disease_pheno_df.to_csv('data/features/diseases-pheno-sim.csv')

