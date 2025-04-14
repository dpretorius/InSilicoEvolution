import math
import numpy as np
from tensorflow.keras.models import load_model, Model
from tensorflow.keras import backend as K
from tensorflow.keras.layers import Average
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
import tensorflow as tf
from Bio.PDB import *
from math import ceil
from typing import Callable, Union
from tensorflow.keras import backend as K

################ CLASSES #########################

class Matrix:

    def __init__(self, coords, indices):
        self.coords=coords
        self.indices=indices-indices[0]

    #Generates a distance matrix
    # It takes three arguments:
    #l - list of indices,
    #m - list of atom coordinates,
    #n - max x-dimension value of the matrix
    #It creates an array of zeros of appropriate dimensions and subsequently populates it with distance values.



    def make_matrix(self):

        matrix = np.zeros((self.indices[-1]+1, self.indices[-1]+1), np.float)
        for  x in range(0,len(self.coords)):
            for y in range(0, len(self.coords)):
                matrix[self.indices[x], self.indices[y]] = np.linalg.norm(self.coords[x]-self.coords[y])
        #matrix=np.reshape(matrix, (len(matrix),maxdim,1))
        #matrix=matrix.astype('float32')

        #plt.matshow(matrix)
        return matrix

class Protein:

    def __init__(self, path, chain,filetype):
        self.path = path
        self.chain = chain
        self.filetype=filetype
    #global - structure and pdb id
    def set_structure(self):
        if self.filetype=='cif':
            parser = MMCIFParser()
        elif self.filetype=='pdb':
            parser=PDBParser()
        structure=parser.get_structure("sds",self.path)
        self.structure=structure
        #self.pdb_id=self.structure.name

    #more local - chains, chain_ids and first residues
    def set_chains(self):
        self.chain_ids=[x.get_id() for x in self.structure[0]]
        self.first_res=self.structure[0][self.chain].get_unpacked_list()[0].get_id()[1]

    def __str__(self):

        return ("Protein_ID " + str(self.pdb_id) + " with chains " + str(self.chain_ids))

    #obtain arrays of atom coordinates and indices based on specified atom type
    def set_atom_lists(self,atom):

        self.coords=np.asarray([residue[str(atom)].get_coord()
                                            for residue in self.structure[0][self.chain] if len(residue)>=4 and str(atom) in residue])

        self.indices=np.asarray([residue.get_id()[1] for residue in self.structure[0][self.chain] if len(residue)>=4 and str(atom) in residue])

    def set_matrices(self):

        distance_matrix1=Matrix(self.coords, self.indices)
        self.atom_matrices=distance_matrix1.make_matrix()
    def output_matrices(self, path2):

        [np.savetxt(str(path2) + "/dmca_matrix_" + str(self.pdb_id) + "_" + str(self.chain_ids[i]) + "_clean_gemmi.csv", self.atom_matrices.get_matrix()[i], delimiter=',') for i in range(len(self.chain_ids))]

    def set_up(self):
        self.set_structure()
        self.set_chains()
        self.set_atom_lists('CA')
        self.set_matrices()

    def return_matrix_and_first_res(self):
        first_res_chain=self.first_res
        distance_matrix_chain=self.atom_matrices

        return first_res_chain, distance_matrix_chain, self.chain_ids


class Prediction:

    def __init__(self, first_res_chain, distance_matrix_chain,loaded_model,sub_matrix_size, threshold, maxlength):
        self.first_res_chain = first_res_chain
        self.distance_matrix_chain=distance_matrix_chain
        self.loaded_model=loaded_model
        self.sub_matrix_size=sub_matrix_size
        self.threshold=threshold
        self.maxlength=maxlength
        self.alpha_predictions=[]
        self.alphabeta_predictions=[]
        self.beta_predictions=[]

    def process_matrix_for_prediction(self):

        def split(array, nrows, ncols):

            r, h = array.shape
            return(array.reshape(h//nrows, nrows, -1, ncols).swapaxes(1, 2).reshape(-1, nrows, ncols, 1))
        def diagonal_split(matrix,sub_matrix_size):
            half_sub_matrix=int(sub_matrix_size/2)
            slices=[]
            for i in range(0,len(matrix)-half_sub_matrix,half_sub_matrix):
                slices.append(matrix[i:i+sub_matrix_size,i:i+sub_matrix_size])
            return np.asarray(slices)
        # obtain value to pad with
        dm_div = len(self.distance_matrix_chain)/self.sub_matrix_size
        pad_value = int(((ceil(dm_div) - dm_div)*self.sub_matrix_size))
        half_pad=int(pad_value/2)
        lpad=half_pad+int(self.sub_matrix_size/4)
        rpad=pad_value-half_pad+int(self.sub_matrix_size/4)

        # pad distance matrix
        distance_matrix_padded = np.pad(self.distance_matrix_chain, ((lpad, rpad), (lpad, rpad)), 'constant', constant_values=(0, 0))
        distance_matrix_padded=distance_matrix_padded/np.max(distance_matrix_padded)
        #plt.matshow(distance_matrix_padded)


        # split up distance matrix into diagonal slices
        distance_matrix_reshaped = diagonal_split(distance_matrix_padded, self.sub_matrix_size)
        distance_matrix_reshaped=np.expand_dims(distance_matrix_reshaped,axis=-1)

        self.distance_matrix_reshaped=distance_matrix_reshaped
        self.pad_value=lpad

    def make_prediction(self):
        total_predictions=np.array([])
        self.process_matrix_for_prediction()

        predictions=self.loaded_model.predict(self.distance_matrix_reshaped)
        #predictions=predictions>self.threshold

        predictions_argmaxed=np.argmax(predictions, axis=3)
        prediction_matrix=np.zeros((predictions_argmaxed.shape[0]*64,predictions_argmaxed.shape[0]*64))

        for i in range(predictions_argmaxed.shape[0]):
            prediction_matrix[i*64:i*64+64,i*64:i*64+64]=predictions_argmaxed[i]

        #plt.matshow(prediction_matrix)

        self.prediction_matrix=prediction_matrix
        for i in range(self.pad_value-int(self.sub_matrix_size/4),len(self.prediction_matrix)):

            for j in range(self.pad_value-int(self.sub_matrix_size/4),len(self.prediction_matrix)):
                if i==j and i<=self.maxlength+self.pad_value-int(self.sub_matrix_size/4) and j<=self.maxlength+self.pad_value-int(self.sub_matrix_size/4):

                    total_predictions=np.append(total_predictions,self.prediction_matrix[i,j])

        #total_predictions[total_predictions>3]=0
        return total_predictions, predictions, self.pad_value

    def get_final_onehot(self):
        distance_matrix_reshaped=self.process_matrix_for_prediction()
        total_predictions=np.array([])

        predictions=self.loaded_model.predict(self.distance_matrix_reshaped)
        prediction_matrix=np.zeros((predictions.shape[0]*64,predictions.shape[0]*64,6))

        for i in range(predictions.shape[0]):
            prediction_matrix[i*64:i*64+64,i*64:i*64+64]=predictions[i]

        one_hot_preds=[]
        for i in range(self.pad_value-int(self.sub_matrix_size/4),len(prediction_matrix)):
            for j in range(self.pad_value-int(self.sub_matrix_size/4),len(prediction_matrix)):
                if i==j and i<=self.maxlength+self.pad_value-int(self.sub_matrix_size/4) and j<=self.maxlength+self.pad_value-int(self.sub_matrix_size/4):
                    one_hot_preds.append(prediction_matrix[i,j])

        return [i[0:4] for i in one_hot_preds]

    def get_solenoid_predictions(self):

        total_preds, predictions, self.pad_value=self.make_prediction()
        solenoid_list=[]
        beta_list=[]
        alphabeta_list=[]
        alpha_list=[]
        propeller_list=[]
        beta_sandwich_list=[]


        for i in range(len(total_preds)):

            if total_preds[i]==1:
                beta_list.append(i+self.first_res_chain)
            elif total_preds[i]==2:
                alphabeta_list.append(i+self.first_res_chain)
            elif total_preds[i]==3:
                alpha_list.append(i+self.first_res_chain)
            elif total_preds[i]==4:
                propeller_list.append(i+self.first_res_chain)
            elif total_preds[i]==5:
                beta_sandwich_list.append(i+self.first_res_chain)

        return np.unique(beta_list), np.unique(alphabeta_list), np.unique(alpha_list), np.unique(propeller_list), np.unique(beta_sandwich_list)

def get_final_scores(structure_path,chain_id='A'):
    loaded_model=load_model('./Conv2D_dmca_multiclass_norm_filesplit_s128_b64_e100_split0.2_1645094358_diagonal_checkpoint.h5')
    protein1=Protein(structure_path,chain=chain_id, filetype='pdb')
    protein1.set_up()
    first_residue,distance_matrix,chain_list=protein1.return_matrix_and_first_res()
    maxlength=protein1.indices[-1]-first_residue
    prediction1=Prediction(first_res_chain=first_residue, distance_matrix_chain=distance_matrix,loaded_model=loaded_model,sub_matrix_size=128, threshold=0.90, maxlength=maxlength)
    final_scores=prediction1.get_final_onehot()
    return final_scores
