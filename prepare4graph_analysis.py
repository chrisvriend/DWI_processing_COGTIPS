#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to prepare the matrices derived from MRtrix3 anatomical constrainted tractography
(i.e. MRTRIX_construct_connectome.sh) for graph analysis using the Brain Connectivity toolbox in Matlab.

@author: Chris Vriend & Bernardo Maciel - 2020
Amsterdam UMC | location VUmc - dpt. Anatomy & Neurosciences | Psychiatry

"""

########################
# IMPORT MODULES       #
########################

import os
import sys
import pandas as pd
import numpy as np
from scipy.io import savemat
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns


########################
# INPUT VARIABLES #
########################
workdir='/mnt/anw-archive/NP/yvdwerf/archive_COGTIPS/analysis/DWI/DWI_TBSS/connectomes/BNA'
outputdir=os.path.join(workdir,'output2')
yeofile=os.path.join(workdir,'YEO_BNA_subnetwork_labels.txt')

#BNA nodes to exclude (e.g. because outside FOV or no tracts stemming from these region)
EXnodes=[117,118]


#########################
if not os.path.exists(outputdir):
        os.makedirs(outputdir)

os.chdir(workdir)
# list comprehension | find files in work directory that start with sub
subjfiles=[subjfiles for subjfiles in os.listdir(workdir)
      if (subjfiles.endswith('connectome.BNA.csv') )];

# sort files and names
subjfiles.sort()


##########################
#   FUNCTIONS            #
##########################
def symmetrize(mat):
        """
        Transforms an upper triangular matrix in a symmetric matrix by copying the values from position i,j to position j,i

        Parameters:
        mat (numpy.array): an upper triangular (square) matrix

        Returns:
        numpy.array: symmetric matrix
        """

        # Variables to store matrix shape
        num_rows = mat.shape[0]
        num_cols = mat.shape[1]

        # Copy upper triangle to lower triangle
        for i in range(num_rows):
                for j in range(i, num_cols):
                        mat[j, i] = mat[i, j]

        return mat





def modYEOandMatrix(yeofile,matrixfile,nodes):
    """
    Function to reshuffle the subnetwork 'key' file according to the nodes we
    want to remove. and do the same for the SC matrix.
    It writes the output file in the same directory.

    Parameters
    ----------
    yeofile : string
        A string with the filepath to the list of nodes.
    matrixfile : string
        A string with the filepath to the SC matrix
    nodes : list
        The list of nodes to be removed.

    Returns
    -------
    None.

    """

    filebase=(os.path.splitext(yeofile))[0]
    basematrix=(os.path.splitext(matrixfile))[0]
    yeo=pd.read_csv(yeofile, delimiter="\t",names=['BNA','YEO'])
    mat = np.genfromtxt(matrixfile, skip_header = 0, delimiter = ' ')

    # symmetrize matrix + convert to dataFrame
    PDmatrix=pd.DataFrame(symmetrize(mat))

    if yeo.shape[0] != PDmatrix.shape[0]:
       sys.exit('YEO label file and connectivity matrix do not have the same length - before exclusion of nodes')

    for node in nodes:
		# Get names of indexes for
        indexNames = yeo[yeo['BNA'] == node ].index
#        Delete these row indexes from dataFrame
        yeo.drop(indexNames , inplace=True)
        # delete row
        PDmatrix.drop(indexNames,axis=0,inplace=True)
        # delete column
        PDmatrix.drop(indexNames,axis=1,inplace=True)

    if yeo.shape[0] != PDmatrix.shape[0]:
       sys.exit('YEO label file and connectivity matrix do not have the same length - after exclusion of nodes')

    # reset index to be consistent with numpy SC matrix later on
    yeo.reset_index(inplace=True,drop=True)

    oddidx=yeo.index[yeo.BNA.isin([num for num in pd.Series(yeo.BNA) if num % 2 == 1])].tolist()
    evenidx=yeo.index[yeo.BNA.isin([num for num in pd.Series(yeo.BNA) if num % 2 == 0])].tolist()


    # reshuffle legend file with odd and even using BNA label as 'key'
    oddYEO=yeo[yeo.BNA.isin([num for num in pd.Series(yeo.BNA) if num % 2 == 1])]
    evenYEO=yeo[yeo.BNA.isin([num for num in pd.Series(yeo.BNA) if num % 2 == 0])]

        # concatenate odd and even, reset index
    yeo2=pd.concat([oddYEO,evenYEO],ignore_index=True)

        #start index at 1 to be consistent with matlab
    yeo2.index +=1

    if not os.path.exists(os.path.join(outputdir,(filebase + '_mod2sample.txt'))):

        # write to txt w/ index
        # yeo2.to_csv(filebase + '_mod2sample.csv',header=False,sep='\t')

        # write to txt w/o index (only needs to be saved once)
        yeo2.to_csv(os.path.join(outputdir,(filebase + '_mod2sample.txt')),index=False,header=False,sep='\t')

    ###########################
    # reshuffle SC matrix
    ###########################
    # dataframe to numpy matrix
    mat=PDmatrix.to_numpy()
    # select rows
    odd = mat[oddidx,:]
    even = mat[evenidx,:]
    mat=np.concatenate((odd,even),axis=0)
    del(odd,even)
    # select columns
    odd = mat[:,oddidx]
    even = mat[:,evenidx]
    # concatenate
    matrix_reshuffled=np.concatenate((odd,even),axis=1)
    # save to txt file
    np.savetxt(os.path.join(outputdir,(basematrix + '_reshuffled.txt')),matrix_reshuffled)
    # save to mat file
    mdic={"connmatrix":matrix_reshuffled,"label":"none"}
    savemat(os.path.join(outputdir,(basematrix + "_reshuffled.mat")),mdic)

    ###########################
    # make heatmap
    ###########################
    # below code for creating the heatmap is derived from: https://github.com/multinetlab-amsterdam/network_TDA_tutorial


    #Creating a DataFrame with will have the brain areas as rows and column names for easier plotting of heatmap
    matrixdiagNaN = matrix_reshuffled.copy()
    np.fill_diagonal(matrixdiagNaN,np.nan)
    # compute  mean
    mu=np.nanmean(matrixdiagNaN)
    # compute nonzero std
    sig=np.nanstd(matrixdiagNaN)

    # z-transform matrix
    matrixdiagNaN2=(matrixdiagNaN-mu) / sig

    temp=matrixdiagNaN2 - np.nanmin(matrixdiagNaN2)
    matrixdiagNaN3=np.divide(temp,np.nanmax(temp))


    figmatrix = pd.DataFrame(matrixdiagNaN3)
    figmatrix.columns=pd.Series(yeo2.BNA)
    figmatrix.index=pd.Series(yeo2.BNA)
    #figmatrix = figmatrix.sort_index(0).sort_index(1)


    #This mask variable gives you the possibility to plot only half of the correlation matrix.
    #If you want the full version, alter parameter  [ mask=mask --> half matrix, mask=None --> full matrix

    mask = np.zeros_like(figmatrix.values, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True

    plt.figure(figsize = (10, 10))
    _ = sns.heatmap(figmatrix, cmap='coolwarm', cbar=True, square=False, mask=None)
    #save heatmap
    plt.savefig(os.path.join(outputdir,(basematrix + '_heatmap_matrix.png')), format='png')
    # clear filgure
    plt.clf()
    plt.close()


##########################
#   PROCESS THE DATA     #
##########################

if __name__ == '__main__':


    for x in tqdm(range(len(subjfiles))):
        subj=subjfiles[x]
        modYEOandMatrix(yeofile,subj,EXnodes)
