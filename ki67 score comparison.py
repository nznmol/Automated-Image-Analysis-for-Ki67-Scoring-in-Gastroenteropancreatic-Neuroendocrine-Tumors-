import random
#import xlrd
import math
import statistics
import openpyxl
from openpyxl import Workbook
import numpy as np
from operator import add
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from scipy.cluster.vq import vq
from scipy import optimize
import glob
import os
import re
import pandas as pd
from pandas import read_excel
import pickle
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pandas import ExcelWriter
import operator
from scipy import stats
from sklearn.preprocessing import RobustScaler

##############################################################################################
#reading coordinate of total number of detected tumor cells of manual count

ImageJ_dic = {}
ImageJ_count_dic = {}
ImageJ_coordinates = []
tile_nrs = []
list_j_max_x = []
list_j_max_y = []
list_j_min_x = []
list_j_min_y = []

j_counter = 0
for file in list(glob.glob('coordinate/ImageJ/ImageJ*.xlsx')):
    dfJ = pd.read_excel(file, usecols = [5, 6], engine='openpyxl')
    dfJ['Slide'] = os.path.basename(file)
    dfJ['Slide'] = dfJ['Slide'].map(lambda x: x.lstrip('ImageJ').rstrip('.xlsx'))
    tile_nr = dfJ['Slide'].iloc[1]
    ImageJ_y_coordinates = dfJ['Y'].tolist()
    ImageJ_x_coordinates = dfJ['X'].tolist()
    ImageJ_x_y_coordinates = list(zip(ImageJ_x_coordinates, ImageJ_y_coordinates))
    ImageJ_coordinates.append(ImageJ_x_y_coordinates)
#imageJ_dic retuns a dictionary: key is the tile number, values are the coordinate of all tumor cells
    ImageJ_dic [tile_nr] = ImageJ_x_y_coordinates
    tile_nrs.append(tile_nr)
    list_j_max_x.append(dfJ['X'].max())
    list_j_max_y.append(dfJ['Y'].max())
    list_j_min_x.append(dfJ['X'].min())
    list_j_min_y.append(dfJ['Y'].min())
    counter = dfJ["X"].count()
#ImageJ_count_dic returns a dictionary: key:tile number, value: total number of tumor cells in the specified tile
    ImageJ_count_dic [tile_nr] = counter

total_j_max_y = max(list_j_max_y)  
total_j_max_x = max(list_j_max_x)
total_j_min_y = max(list_j_min_y)  
total_j_min_x = max(list_j_min_x)

#make the data float instead of string
imagej_x_coordinates_int = list (map(float, ImageJ_x_coordinates))
imagej_y_coordinates_int = list (map(float, ImageJ_y_coordinates))



##############################################################################################

#reading coordinate of positive tumor cells of manual count

ImageJ_positive_dic = {}
ImageJ_positive_coordinates = []
positive_tile_nrs = []

for file in list(glob.glob('coordinate/ImageJ/PositiveCell/*.xlsx')):
    dfJ_positive = pd.read_excel(file, usecols = [5, 6], engine='openpyxl')
    dfJ_positive['Slide'] = os.path.basename(file)
    dfJ_positive['Slide'] = dfJ_positive['Slide'].map(lambda x: x.rstrip('.xlsx'))
    tile_nr = dfJ_positive['Slide'].iloc[0]
    ImageJ_positive_y_coordinates = dfJ_positive['Y'].tolist()
    ImageJ_positive_x_coordinates = dfJ_positive['X'].tolist()
    ImageJ_positive_x_y_coordinates = list(zip(ImageJ_positive_x_coordinates, ImageJ_positive_y_coordinates))
    ImageJ_positive_coordinates.append(ImageJ_x_y_coordinates)
#ImageJ_positive_dic returns a dictionary : key:tile number, values: positive cells coordinate
    ImageJ_positive_dic [tile_nr] = ImageJ_positive_x_y_coordinates
    positive_tile_nrs.append(tile_nr)
#include the tile numbers that are not exisitng, becasue they didn't have a positive cell
    ImageJ_tile_no_positive = list(set(tile_nrs) - set(positive_tile_nrs))
    ImageJ_no_positive_dic = {key: [] for key in ImageJ_tile_no_positive}
    ImageJ_positive_dic.update(ImageJ_no_positive_dic) 

#############################################################################################

#reading coordinate of total number of detected cells and positive tumor cells by Aiforia 
#aiforia read all tumor cells and positive tumor cells at once

Aiforia_dic = {}
Aiforia_positive_dic = {}
Aiforia_coordinates = []
Aiforia_coordinates_positive = []
list_a_max_x = []
list_a_max_y = []

counter = 0
for file in list(glob.glob('coordinate/Aiforia/IA_details__*.xlsx')):
    dfA = pd.read_excel(file, usecols = [0, 9, 11, 21, 22], engine='openpyxl')
    dfA.columns = ['Slide','Result type','Cell Type','X','Y']
    dfA = dfA[dfA['Result type'].str.contains('Object')]
    counter = dfA['Cell Type'].str.count("positive").sum()
    dfA['Slide'] = dfA['Slide'].str.replace(r"\(.*\)", "")
    dfA['Slide'] = dfA['Slide'].replace(r'_','-', regex=True)
    tile_nr = dfA['Slide'].iloc[1]
    Aiforia_y_coordinates = [y for y in dfA['Y'].tolist() if str(y) != 'nan']
    Aiforia_x_coordinates = [x for x in dfA['X'].tolist() if str(x) != 'nan']
    cell_type = [i for i in dfA['Cell Type'].tolist() if str(i) != 'nan']
    Aiforia_x_y_coordinates = list(zip(Aiforia_x_coordinates, Aiforia_y_coordinates))
    Aiforia_coordinates.append(Aiforia_x_y_coordinates)
#Aiforia_dic returns a list of tuple of tuples: [(coordinates of all tumor cells in a tile),(type:positivr or negative)]
    Aiforia_dic [tile_nr] = list(zip(Aiforia_x_y_coordinates, cell_type))
#Aiforia_positive_dic returns key:tile number, value:number of positive cells
    Aiforia_positive_dic [tile_nr] = counter
    list_a_max_x.append(dfA['X'].max())
    list_a_max_y.append(dfA['Y'].max())

total_a_max_x = max(list_a_max_x)
total_a_max_y = max(list_a_max_y)


############################################################################################

#reading coordinate of total number of detected cells ImageScope (with Nuclear algorithm)

ImageScope_dic = {}
ImageScope_count_dic = {}
ImageScope_coordinates = []
tile_nrs = []
list_scope_max_x = []
list_scope_max_y = []
list_scope_min_x = []
list_scope_min_y = []

j_counter = 0
for file in list(glob.glob('coordinate/ImageScope/ImageScope*.xlsx')):
    df_scope = pd.read_excel(file, usecols = [5, 6], engine='openpyxl')
    df_scope['Slide'] = os.path.basename(file)
    df_scope['Slide'] = df_scope['Slide'].map(lambda x: x.lstrip('ImageScope').rstrip('.xlsx'))
    tile_nr = df_scope['Slide'].iloc[1]
    ImageScope_y_coordinates = df_scope['Y'].tolist()
    ImageScope_x_coordinates = df_scope['X'].tolist()
    ImageScope_x_y_coordinates = list(zip(ImageScope_x_coordinates, ImageScope_y_coordinates))
    ImageScope_coordinates.append(ImageScope_x_y_coordinates)
#ImageScope_dic returns: key is the tile number, values is a list of tupuls (coordinate)
    ImageScope_dic [tile_nr] = ImageScope_x_y_coordinates
    tile_nrs.append(tile_nr)
    list_scope_max_x.append(df_scope['X'].max())
    list_scope_max_y.append(df_scope['Y'].max())
    list_scope_min_x.append(df_scope['X'].min())
    list_scope_min_y.append(df_scope['Y'].min())
    counter = df_scope["X"].count()
#ImageScope_count_dic returns: key is the tile numbers, values are the number of totall cell count
    ImageScope_count_dic [tile_nr] = counter
   
total_scope_max_y = max(list_scope_max_y)  
total_scope_max_x = max(list_scope_max_x)
total_scope_min_y = max(list_scope_min_y)  
total_scope_min_x = max(list_scope_min_x)

#make the data float instead of string
ImageScope_x_coordinates_int = list (map(float, ImageScope_x_coordinates))
ImageScope_y_coordinates_int = list (map(float, ImageScope_y_coordinates))

#################################################################################################

#reading coordinate of positive tumor cells ImageScope (Nuclear algorithm)

ImageScope_positive_dic = {}
ImageScope_positive_coordinates = []
positive_tile_nrs = []

for file in list(glob.glob('coordinate/ImageScope/PositiveCell/*.xlsx')):
    df_scope_positive = pd.read_excel(file, usecols = [5, 6], engine='openpyxl')
    df_scope_positive['Slide'] = os.path.basename(file)
    df_scope_positive['Slide'] = df_scope_positive['Slide'].map(lambda x: x.lstrip('ImageScope-Positive').rstrip('.xlsx'))
    tile_nr = df_scope_positive['Slide'].iloc[0]
    ImageScope_positive_y_coordinates = df_scope_positive['Y'].tolist()
    ImageScope_positive_x_coordinates = df_scope_positive['X'].tolist()
    ImageScope_positive_x_y_coordinates = list(zip(ImageScope_positive_x_coordinates, ImageScope_positive_y_coordinates))
    ImageScope_positive_coordinates.append(ImageScope_x_y_coordinates)
#ImageScope_positive_dic retuns: key is the tile number, values is a list of tupuls (coordinate)
    ImageScope_positive_dic [tile_nr] = ImageScope_positive_x_y_coordinates
    positive_tile_nrs.append(tile_nr)
#include the tiles that are not exisitng, becasue they didn't have a positive cell
    ImageScope_tile_no_positive = list(set(tile_nrs) - set(positive_tile_nrs))
    ImageScope_no_positive_dic = {key: [] for key in ImageScope_tile_no_positive}
    ImageScope_positive_dic.update(ImageScope_no_positive_dic) 

########################################################################################################
#Functions

def scaling_coordinate(a_coord, j_coord):
    # TODO make a_max_x, a_max_y, j_max_x, j_max_y available as arguments to the function
    # and pass them to the function wheneveryou call them
    
    cell_type = [coord[1] for coord in a_coord]
    a_coordinate = [tuple(j for j in i[0]) for i in a_coord] 
    a_x_coord = [coord[0] for coord in a_coordinate]
    a_y_coord = [coord[1] for coord in a_coordinate]
    # TODO it's unfortunate to calculate the maximum every time - move the calculation
    # further (outside of this function!), and do it only once
    a_y_coordinates_ranged = [y*(total_j_max_y) / (total_a_max_y) for y in a_y_coord]
    a_x_coordinates_ranged = [x*(total_j_max_x) / (total_a_max_x) for x in a_x_coord]
    
    a_coord_ranged = list(zip(a_x_coordinates_ranged, a_y_coordinates_ranged))
    a_coord_ranged_celltype = list(zip(a_coord_ranged, cell_type))
    return a_coord_ranged, a_coord_ranged_celltype


def compare_coordinates(a_coord, j_coord):
    distance_square = []
    distance = cdist(a_coord, j_coord, 'euclidean')
    max_distance = np.max(distance)
    if len(j_coord) < len(a_coord):
        zeros = np.zeros((len(a_coord), len(a_coord)-len(j_coord)))
        zeros = zeros+1000*(max_distance+1)
        distance_square = np.append(distance, zeros, 1)
    elif len(j_coord) > len(a_coord):
        zeros = np.zeros((len(j_coord)-len(a_coord), len(j_coord)))
        zeros = zeros+1000*(max_distance+1)
        distance_square = np.append(distance, zeros, 0)
    else:
        distance_square = distance

    distance_square_array = np.array(distance_square)
    row_ind, col_ind = linear_sum_assignment(distance_square_array)
    assignment = list(zip(row_ind, col_ind))

    filtered_assignment = [ (row,col) for row,col in assignment
       if distance_square_array[row,col] <= max_distance ]

    return max_distance, distance_square_array, assignment

def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() 
                            if len(locs)>1)


def find_minimum_index(list_of_costs):
        return np.argmin(list_of_costs)

    
    
def scale_coordinates_dic(Aiforia_dic, ImageJ_dic, tile_nrs):
    a_coord_ranged = {}
    a_positive_ranged = {}
    for tile in tile_nrs:
        Aiforia_coord = Aiforia_dic[tile]
        ImageJ_coord = ImageJ_dic[tile]
        a_coord_ranged[tile], a_coord_ranged_celltype = scaling_coordinate(Aiforia_coord, ImageJ_coord)
        a_positive_ranged[tile] = [i[0] for i in a_coord_ranged_celltype if 'positive' in i[1]]
    return a_coord_ranged, a_positive_ranged


        
def best_match(aiforia_dic, imagej_dic, tile_nrs):
    unique_optimal_matches_dic = {}
    for tile in tile_nrs:
        #print(aiforia_dic[tile])
        aiforia = aiforia_dic[tile]
        imagej = imagej_dic[tile]
        if len(aiforia) != 0 and len(imagej) != 0:
            code, dist = vq (aiforia, imagej)
            # tuple of (n, [a,b]) where n is from imagej and was assigned to both a and b from aiforia [a,b] with a,b as above
            duplicate = []                
            duplicate_index = []
            for i in sorted(list_duplicates(code)):
                duplicate.append(i)
                duplicate_index.append(i[1]) 
            #tuple (c(a),c(b)) where c(a) is the cost of the assingment for a, and same for b
            duplicate_dist = [[dist[number] for number in group] for group in duplicate_index]

            best_imagej_match = {}
            for (imagej_mapping, aiforia_points) in duplicate:
                best_imagej_match[imagej_mapping] = aiforia_points[np.argmin([dist[a] for a in aiforia_points])]
                
            #list of tuples (aiforia, imagej, cost)
            unique_optimal_matches_dic[tile] = [(aiforia, imagej, dist[aiforia]) for (aiforia, imagej) 
                in enumerate(code) if imagej not in best_imagej_match or
                aiforia == best_imagej_match[imagej]]
        else:
            unique_optimal_matches_dic[tile] = []
    return unique_optimal_matches_dic



def visualizing_assignment(aiforia_dic, imagej_dic, unique_optimal_matches_dic, 
                           MATCH_THRESHOLD, TARGET_FOLDER, label_a="aiforia", label_b="imagej"):
    #TARGET_FOLDER = "vq-assignment-visualization/PositiveCell/imagescope-imagej"
    #cell annotation size in aiforia was 6um
    MATCH_THRESHOLD = 30
    
    for tile in tile_nrs:
        aiforia = aiforia_dic[tile]
        imagej = imagej_dic[tile]
        unique_optimal_matches = unique_optimal_matches_dic[tile]
        tileImagePath = os.path.join("80tiles", tile + ".tif")
        # make sure we dont get trouble with tiles with underscores
        if not os.path.exists(tileImagePath):
            tileImagePath = tileImagePath.replace('-','_')
        #plot the jpg image as the background
        im = mpimg.imread(tileImagePath)
        plt.figure(figsize=(15,15))
        plt.imshow(im)

        #ploting the coordinates from Aiforia and ImageJ
        plt.plot([c[0] for c in aiforia], [c[1] for c in aiforia], 'bo', markersize = 12, label=label_a)
        plt.plot([c[0] for c in imagej], [c[1] for c in imagej], 'rs', markersize = 7, label=label_b)
        plt.legend()

        #connection the Aiforia coordinate with the best assigned math
        for (aifo,imaj,cost) in unique_optimal_matches:
            if cost < MATCH_THRESHOLD:
                plt.plot([ aiforia[aifo][0], imagej[imaj][0]], 
                    [aiforia[aifo][1], imagej[imaj][1]],'k-')
        targetFile = os.path.splitext(os.path.basename(tileImagePath))[0] + ".png"
        targetPath = os.path.join(TARGET_FOLDER, targetFile)
        plt.savefig(targetPath)
        plt.close('all')
    return 


def max_x(l):
    return max(map(lambda c: c[0], l))

def max_y(l):
    return max(map(lambda c: c[1], l))

def min_x(l):
    return min(map(lambda c: c[0], l))

def min_y(l):
    return min(map(lambda c: c[1], l))

def scale(data, data_max, other_max):
    return [d*other_max / data_max for d in data]

def transform_array(array, transform_vector): 
    """ transform_array is a simple helper function to apply transform_vector 
    to the array or list of points """ 
    return np.array(list(map(lambda point: (point + transform_vector[:2]) * 
        transform_vector[2:], array))) 
 
def vq_match(one_list, other_list): 
    """ simple vq match between two lists, returning a list of the form 
    (i,j,c(i,j)) where the i-th element of one_list gets assigned to the j-th 
    alignment of other_list, at a cost (distance) of c(i,j) """ 
    code, dist = vq (one_list, other_list) 
    assignment_with_costs = [(one_index, other_index, dist[one_index]) 
        for (one_index, other_index) in enumerate(code)] 
    return assignment_with_costs 
    
def assigned_pairs_distances(one_list, other_list): 
    """ assigned_pairs_distances is a cost function that simply finds an 
    assignment between the two given lists using vq, and returns the sum of the 
    costs of said assignment """ 
    assignment_with_costs = vq_match(one_list.tolist(), other_list.tolist()) 
    costs = [assignment[2] for assignment in assignment_with_costs] 
    sum_of_costs = np.sum(costs) 
    return sum_of_costs 
    
def general_function_to_optimize(transform, one_list, other_list, cost_func): 
    """ general_function_to_optimize describes the function we want to 
    optimize, namely transforming the second list and returning the cost 
    between the first list and the transformed second list """ 
    transformed_other_list = transform_array(other_list, transform) 
    return cost_func(one_list, transformed_other_list) 
    
def find_transformed_list_with_optimize(original_points, other_points): 
    """ takes two lists of points and tries to transform the second one into 
    the first, returning the transformed list and the transformation vector in 
    the form (shift_x, shift_y, scale_x, scale_y) """ 
    
    # transform into numpy arrays for easier notation of calculation 
    first_list = np.array(original_points) 
    other_list = np.array(other_points) 
 
    # rewrite the general function to optimize for as a function of the 
    # transformation vector, fixing the two lists and the cost function 
    optimize_for_transform_vector = lambda transform_vector:         general_function_to_optimize(transform_vector, first_list, other_list, 
                assigned_pairs_distances) 
 
    # initial guess: the identity transformation 
    initial_transform_vector = np.array([0,0,1,1]) 
 
    # run the optimization also, check for error and throw an error if the 
    # optimization failed 
    optimization_result = optimize.minimize(optimize_for_transform_vector, 
            initial_transform_vector, method="Powell") 
    if not optimization_result.success: 
        raise Exception(optimization_result.message) 
    found_transform_vector = optimization_result.x 
 
    # transform other list, return transformed list and transform vector 
    transformed_other_list = transform_array(other_list, found_transform_vector) 
    return transformed_other_list, found_transform_vector 

def transform_dict(first_dict, other_dict, tiles):
    transformed_other_dict = {}
    transformation_dic = {}
    for tile in tiles:
        first_list = first_dict[tile]
        other_list = other_dict[tile]
        if len(first_list) == 0 or len(other_list) == 0:
            transformed_other_dict[tile] = other_list
            continue
        transformed_other_list, found_transform_vector = find_transformed_list_with_optimize(first_list, other_list)
        transformed_other_dict[tile] = transformed_other_list
        transformation_dic[tile] = found_transform_vector
    return transformed_other_dict, transformation_dic

##########################################################################################################################
#assignment for all tumor cells(aiforia-imagej)

#'a_coord_ranged' and 'a_positive_ranged' return dictionary, key is the tile number, values are coordinates
a_coord_ranged, a_positive_ranged = scale_coordinates_dic(Aiforia_dic, ImageJ_dic, tile_nrs)
unique_optimal_matches_dic = best_match(a_coord_ranged, ImageJ_dic, tile_nrs)
visualizing_assignment(a_coord_ranged, ImageJ_dic, unique_optimal_matches_dic, MATCH_THRESHOLD, label_a ='Aiforia', label_b="ImageJ", TARGET_FOLDER = "vq-assignment-visualization/aiforia-imagej")

aiforia_true_cell_detection_vq = []
aiforia_false_cell_detection_vq  = []
aiforia_missing_cell_vq = []
aiforia_negative_cells_vq = []
aiforia_score_vq = []
aiforia_total_cells_vq = []

imageJ_negative_cells_vq = []
imageJ_score_vq = []
imageJ_total_cells_vq = []

MATCH_THRESHOLD = 30
cell_count_vq = pd.DataFrame()
d_vq =[]

for tile in sorted(tile_nrs):
    unique_optimal_matches = unique_optimal_matches_dic[tile]
    true_detection = len([i for i in unique_optimal_matches if i[2] < MATCH_THRESHOLD])
    aiforia_true_cell_detection_vq.append(true_detection)
    aiforia_total_cells = len(Aiforia_dic[tile])
    aiforia_total_cells_vq.append(aiforia_total_cells)
    false_detection = len(Aiforia_dic[tile]) - true_detection
    aiforia_false_cell_detection_vq.append(false_detection)
    missing_cell = len(ImageJ_dic[tile]) - true_detection
    aiforia_missing_cell_vq.append(missing_cell)
    aiforia_positive_cells = Aiforia_positive_dic[tile]
    aiforia_negative_cells = len(Aiforia_dic[tile]) - Aiforia_positive_dic[tile]
    aiforia_negative_cells_vq.append(aiforia_negative_cells)
    aiforia_score = (aiforia_positive_cells / aiforia_total_cells)*100
    aiforia_score_vq.append(aiforia_score)
    
    imageJ_total_cells = len(ImageJ_dic[tile])
    imageJ_total_cells_vq.append(imageJ_total_cells)
    imageJ_positive_cells = len(ImageJ_positive_dic[tile])
    imageJ_negative_cells = len(ImageJ_dic[tile]) - len(ImageJ_positive_dic[tile])
    imageJ_negative_cells_vq.append(imageJ_negative_cells)
    imageJ_score = (len(ImageJ_positive_dic[tile]) / imageJ_total_cells)*100
    imageJ_score_vq.append(imageJ_score)
      
  
        
    d_vq.append((tile,
    true_detection, 
    false_detection, 
    missing_cell, 
    aiforia_total_cells, 
    aiforia_positive_cells,
    aiforia_negative_cells,
    aiforia_score,
    imageJ_total_cells,
    imageJ_positive_cells,
    imageJ_negative_cells,
    imageJ_score))

    
cell_count_vq = pd.DataFrame(d_vq, columns=(
    'Slide',
    'Aiforia true detection', # true positive
    'Aiforia false detection', #false positive
    'Aiforia missing cell', #false negative
    'Aiforia total cells', 
    'Aiforia positive cells',
    'Aiforia negative cells',
    'Aiforia score',
    'ImageJ total cells',
    'ImageJ positive cells',
    'ImageJ negative cells',
    'ImageJ score'))


last_row = cell_count_vq[['Aiforia true detection', 
    'Aiforia false detection', 
    'Aiforia missing cell', 
    'Aiforia total cells', 
    'Aiforia positive cells',
    'Aiforia negative cells',
    'ImageJ total cells',
    'ImageJ positive cells',
    'ImageJ negative cells']].sum()

last_row['Discovery rate'] = last_row['Aiforia false detection'] / (last_row['Aiforia false detection'] + last_row['Aiforia true detection'])
last_row['Sensitivity(Recall)'] = last_row['Aiforia true detection'] / (last_row['Aiforia true detection'] + last_row['Aiforia missing cell'])
last_row['Precision'] = last_row['Aiforia true detection'] / (last_row['Aiforia true detection'] + last_row['Aiforia false detection'])
cell_count_vq = cell_count_vq.append(last_row, ignore_index=True)

cell_count_vq.to_excel('vq-assignment-visualization/Aiforia_ImageJ_TotalCell_vq.xlsx', index = False)


###################################################################################################################
#assignment for positive tumor cells (aiforia-imagej)

unique_optimal_matches_positive_dic = best_match(a_positive_ranged, ImageJ_positive_dic, tile_nrs)
visualizing_assignment(a_positive_ranged, ImageJ_positive_dic, unique_optimal_matches_positive_dic, MATCH_THRESHOLD, label_a="Aiforia", label_b="ImageJ", TARGET_FOLDER = "vq-assignment-visualization/PositiveCell/aiforia-imagej")

#printing the word file for positive tumor cells (aiforia-imagej)

TARGET_FOLDER = "vq-assignment-visualization/PositiveCell"
aiforia_true_positive_detection_vq = []
aiforia_false_positive_detection_vq = []
aiforia_missing_positive_vq = []
aiforia_negative_cells_vq = []
aiforia_score_vq = []
aiforia_total_cells_vq = []
dic_a_coord_ranged_celltype = {}

imageJ_positive_vq = []
imageJ_negative_vq = []
imageJ_score_vq = []
imageJ_total_cells_vq = []

positive_cell_count_vq = pd.DataFrame()
positive_d_vq =[]

for tile in sorted(tile_nrs):
    unique_optimal_matches_positive = unique_optimal_matches_positive_dic[tile]
    aiforia_true_positive_detection = len([i for i in unique_optimal_matches_positive if i[2] < MATCH_THRESHOLD])
    aiforia_true_positive_detection_vq.append(aiforia_true_positive_detection)
    aiforia_false_positive_detection = len(a_positive_ranged[tile]) - aiforia_true_positive_detection
    aiforia_false_positive_detection_vq.append(aiforia_false_positive_detection)
    aiforia_missing_positive = len(ImageJ_positive_dic[tile]) - aiforia_true_positive_detection
    aiforia_missing_positive_vq.append(aiforia_missing_positive)
    aiforia_total_cells = len(Aiforia_dic[tile])
    aiforia_total_cells_vq.append(aiforia_total_cells)
    aiforia_negative_cells = len(Aiforia_dic[tile]) - len(a_positive_ranged[tile])
    aiforia_negative_cells_vq.append(aiforia_negative_cells)
    aiforia_score = (aiforia_positive_cells / aiforia_total_cells)*100
    aiforia_score_vq.append(aiforia_score)
    aiforia_total_positive = len(a_positive_ranged[tile])
 
    imageJ_total_cells = len(ImageJ_dic[tile])
    imageJ_total_cells_vq.append(imageJ_total_cells)
    imageJ_positive = len(ImageJ_positive_dic[tile])
    imageJ_positive_vq.append(imageJ_positive)
    imageJ_negative = imageJ_total_cells - imageJ_positive
    imageJ_negative_vq.append(imageJ_negative)
    imageJ_score = (imageJ_positive / imageJ_total_cells)*100
    imageJ_score_vq.append(imageJ_score)
    
    
    positive_d_vq.append((tile,
    aiforia_true_positive_detection, 
    aiforia_false_positive_detection, 
    aiforia_missing_positive,  
    aiforia_total_cells, 
    aiforia_total_positive,
    aiforia_negative_cells,
    aiforia_score,
    imageJ_total_cells,
    imageJ_positive,
    imageJ_negative,
    imageJ_score))

     
positive_cell_count_vq = pd.DataFrame(positive_d_vq, columns=('Slide',
    'Aiforia true positive detection', 
    'Aiforia false positive detection', 
    'Aiforia missing positive cells',  
    'Aiforia total cells', 
    'Aiforia total positive cells',
    'Aiforia total negative cells',
    'Aiforia score',
    'ImageJ total cells',
    'ImageJ positive cells',
    'ImageJ negative cells',
    'ImageJ score'))

last_row_positive = positive_cell_count_vq [[
    'Aiforia true positive detection', 
    'Aiforia false positive detection', 
    'Aiforia missing positive cells',  
    'Aiforia total cells', 
    'Aiforia total positive cells',
    'Aiforia total negative cells',
    'ImageJ total cells',
    'ImageJ positive cells',
    'ImageJ negative cells']].sum()


last_row_positive['Discovery rate'] = last_row_positive['Aiforia false positive detection'] / (last_row_positive['Aiforia false positive detection'] + last_row_positive['Aiforia true positive detection'])
last_row_positive['Sensitivity(Recall)'] = last_row_positive['Aiforia true positive detection'] / (last_row_positive['Aiforia true positive detection'] + last_row_positive['Aiforia missing positive cells'])
last_row_positive['Precision'] = last_row_positive['Aiforia true positive detection'] / (last_row_positive['Aiforia true positive detection'] + last_row_positive['Aiforia false positive detection'])
positive_cell_count_vq = positive_cell_count_vq.append(last_row_positive, ignore_index=True)

positive_cell_count_vq.to_excel('vq-assignment-visualization/Aiforia_ImageJ_PositiveCell_vq.xlsx', index = False)

###########################################################################################################################

#assignment all tumor cells(ImageScope-ImageJ)
#Manual scaling by image size and resolution

from PIL import Image
import matplotlib.pyplot as plt
import cv2
import re
from operator import getitem
from functools import partial

MATCH_THRESHOLD = 30
scope_coord_ranged = {}
imagescope_blownup_dic = {}

# simple scaling by size and resolution
# should read the image (by their tile number) and return their size and resolution,

imagej_hw = []
imagescope_hw = []
tile_imagej_num = []
tile_scope_num = []

for file in (glob.glob("./80tiles/ImageScope/*_*.tif")):
        img2 = cv2.imread(file)
        im2 = Image.open(file)
        r2 = im2.info['dpi'][0]
        tile = re.findall("\d+", file)[1:]
        tile_tile = ('-'.join(tile))
        tile_scope_num.append(tile_tile)
        h2,w2,_ = img2.shape
        tuples = h2, w2
        imagescope_hw.append(tuples)
        dict_scope = {key: value for key, value in zip(tile_scope_num, imagescope_hw)}
        table_scope = pd.DataFrame(dict_scope.items(), columns=['tilenumber', 'h2, w2'])
        table_scope['r2'] = im2.info['dpi'][0]
        table_scope['h2, w2'] = table_scope['h2, w2'].tolist()
        table_scope[['h2', 'w2']] = pd.DataFrame(table_scope['h2, w2'].tolist(), index=table_scope.index)


for file in (glob.glob("./80tiles/*_*.tif")):
        img1 = cv2.imread(file)
        im1 = Image.open(file)
        tile = re.findall("\d+", file)[1:]
        tile_tile = ('-'.join(tile))
        tile_imagej_num.append(tile_tile)
        h1,w1,_ = img1.shape
        tuplep = h1,w1
        imagej_hw.append(tuplep)
        dict_imagej = {key: value for key, value in zip(tile_imagej_num, imagej_hw)}
        table_imagej = pd.DataFrame(dict_imagej.items(), columns=['tilenumber', 'h1, w1'])
        table_imagej['h1, w1'] = table_imagej['h1, w1'].tolist()
        table_imagej[['h1', 'w1']] = pd.DataFrame(table_imagej['h1, w1'].tolist(), index=table_imagej.index)


imagej_imagescope_df = pd.merge(table_imagej, table_scope, on='tilenumber')
imagej_imagescope_df = imagej_imagescope_df.assign(coordinate = imagej_imagescope_df.tilenumber.map(ImageScope_dic))


scaled_coordinates = []
coordinate = []
scaled_dict = {}
x_y_scaled = []
x_9_coord_scaled = []
y_9_coord_scaled = []
x_y_9_scaled = []


for (k,v), (k2,v2) in zip (dict_imagej.items(), ImageScope_dic.items()):  
    x_coord = [i[0] for i in ImageScope_dic[k]]
    y_coord = [i[1] for i in ImageScope_dic[k]]

    if k[0] == '9':
        x_coord_scaled = [j / 0.0083 for j in x_coord]
        y_coord_scaled = [j / 0.0083 for j in y_coord]
        x_y_scaled = list(zip(x_coord_scaled, y_coord_scaled))
    else:
        x_coord_scaled = [j * r2 / w2 * w1 for j in x_coord]
        y_coord_scaled = [j * r2 / h2 * h1 for j in y_coord]
        x_y_scaled = list(zip(x_coord_scaled, y_coord_scaled))
    scaled_coordinates.append(x_y_scaled)
    scaled_dict[k] = x_y_scaled

scope_unique_optimal_matches_dictionary = best_match(scaled_dict, ImageJ_dic, tile_nrs)
visualizing_assignment(scaled_dict, ImageJ_dic, scope_unique_optimal_matches_dictionary, MATCH_THRESHOLD, label_a="ImageScope", label_b="ImageJ", TARGET_FOLDER = "vq-assignment-visualization/imagescope-imagej-manual")

imagescope_true_positive_detection_vq = []
imagescope_false_positive_detection_vq = []
imagescope_missing_positive_vq = []
imagescope_negative_cells_vq = []
imagescope_score_vq = []
imagescope_total_cells_vq = []


imageJ_positive_vq = []
imageJ_negative_vq = []
imageJ_score_vq = []
imageJ_total_cells_vq = []

ImageScope_cell_count_vq = pd.DataFrame()
positive_d_vq =[]

for tile in sorted(tile_nrs):
    scope_unique_optimal_matches = scope_unique_optimal_matches_dictionary[tile]
    imagescope_true_positive_detection = len([i for i in scope_unique_optimal_matches if i[2] < MATCH_THRESHOLD])
    imagescope_true_positive_detection_vq.append(imagescope_true_positive_detection)
    imagescope_false_positive_detection = len(ImageScope_dic[tile]) - imagescope_true_positive_detection
    imagescope_false_positive_detection_vq.append(imagescope_false_positive_detection)
    imagescope_missing_positive = len(ImageJ_dic[tile]) - imagescope_true_positive_detection
    imagescope_missing_positive_vq.append(imagescope_missing_positive)
    imagescope_total_cells = len(ImageScope_dic[tile])
    imagescope_total_cells_vq.append(imagescope_total_cells)
    imagescope_negative_cells = len(ImageScope_dic[tile]) - len(scaled_positive_dict[tile])
    imagescope_negative_cells_vq.append(imagescope_negative_cells)
    imagescope_positive_cells = ImageScope_positive_dic[tile]
    imagescope_total_positive = len(scaled_positive_dict[tile])
 
    imageJ_total_cells = len(ImageJ_dic[tile])
    imageJ_total_cells_vq.append(imageJ_total_cells)
    imageJ_positive = len(ImageJ_positive_dic[tile])
    imageJ_positive_vq.append(imageJ_positive)
    imageJ_negative = imageJ_total_cells - imageJ_positive
    imageJ_negative_vq.append(imageJ_negative)
    imageJ_score = (imageJ_positive / imageJ_total_cells)*100
    imageJ_score_vq.append(imageJ_score)
    
    
    positive_d_vq.append((tile,
    imagescope_true_positive_detection, 
    imagescope_false_positive_detection, 
    imagescope_missing_positive,  
    imagescope_total_cells, 
    imagescope_total_positive,
    imagescope_negative_cells,
    imageJ_total_cells,
    imageJ_positive,
    imageJ_negative,
    imageJ_score))

     
ImageScope_cell_count_vq = pd.DataFrame(positive_d_vq, columns=('Slide',
    'imagescope true positive detection', 
    'imagescope false positive detection', 
    'imagescope missing positive cells',  
    'imagescope total cells', 
    'imagescope total positive cells',
    'imagescope total negative cells',
    'ImageJ total cells',
    'ImageJ positive cells',
    'ImageJ negative cells',
    'ImageJ score'))


last_row = ImageScope_cell_count_vq [[
    'imagescope true positive detection', 
    'imagescope false positive detection', 
    'imagescope missing positive cells',  
    'imagescope total cells', 
    'imagescope total positive cells',
    'imagescope total negative cells',
    'ImageJ total cells',
    'ImageJ positive cells',
    'ImageJ negative cells']].sum()

ImageScope_cell_count_vq.to_excel('vq-assignment-visualization/ImageScope_ImageJ_TotalCell_vq.xlsx', index = False)

#########################################################################################################

#positive tumor cells(ImageScope-ImageJ)
#Manual scaling with image size adn resolution

imagescope_positive_blownup_dic = {}
imagescope_positive_blownup_ranged_dic ={}
scaled_positive_coordinates = []
scaled_positive_dict ={}

for (k,v), (k2,v2) in zip (dict_imagej.items(), ImageScope_positive_dic.items()):  
    x_positive_coord = [i[0] for i in ImageScope_positive_dic[k]]
    y_positive_coord = [i[1] for i in ImageScope_positive_dic[k]]

    if k[0] == '9':
        x_positive_coord_scaled = [j / 0.0083 for j in x_positive_coord]
        y_positive_coord_scaled = [j / 0.0083 for j in y_positive_coord]
        x_y_positive_scaled = list(zip(x_positive_coord_scaled, y_positive_coord_scaled))
    else:
        x_positive_coord_scaled = [j * r2 / w2 * w1 for j in x_positive_coord]
        y_positive_coord_scaled = [j * r2 / h2 * h1 for j in y_positive_coord]
        x_y_positive_scaled = list(zip(x_positive_coord_scaled, y_positive_coord_scaled))
    scaled_positive_coordinates.append(x_y_positive_scaled)
    scaled_positive_dict[k] = x_y_positive_scaled
    
unique_optimal_matches_positive_scope_J_dictionary = best_match(scaled_positive_dict, ImageJ_positive_dic, tile_nrs)
visualizing_assignment(scaled_positive_dict, ImageJ_positive_dic, unique_optimal_matches_positive_scope_J_dictionary, MATCH_THRESHOLD, label_a="ImageScope", label_b="imagej", TARGET_FOLDER = "vq-assignment-visualization/PositiveCell/manual-imagescope-imagej")

imagescope_true_positive_detection_vq = []
imagescope_false_positive_detection_vq = []
imagescope_missing_positive_vq = []
imagescope_negative_cells_vq = []
imagescope_score_vq = []
imagescope_total_cells_vq = []
dic_a_coord_ranged_celltype = {}

imageJ_positive_vq = []
imageJ_negative_vq = []
imageJ_score_vq = []
imageJ_total_cells_vq = []

ImageScope_positive_cell_count_vq = pd.DataFrame()
positive_scope_d_vq =[]

for tile in sorted(tile_nrs):
    unique_optimal_matches_scope_positive = unique_optimal_matches_positive_scope_J_dictionary[tile]
    imagescope_true_positive_detection = len([i for i in unique_optimal_matches_scope_positive if i[2] < MATCH_THRESHOLD])
    imagescope_true_positive_detection_vq.append(imagescope_true_positive_detection)
    imagescope_false_positive_detection = len(scaled_positive_dict[tile]) - imagescope_true_positive_detection
    imagescope_false_positive_detection_vq.append(imagescope_false_positive_detection)
    imagescope_missing_positive = len(ImageJ_positive_dic[tile]) - imagescope_true_positive_detection
    imagescope_missing_positive_vq.append(imagescope_missing_positive)
    imagescope_total_cells = len(ImageScope_dic[tile])
    imagescope_total_cells_vq.append(imagescope_total_cells)
    imagescope_negative_cells = len(ImageScope_dic[tile]) - len(scaled_positive_dict[tile])
    imagescope_negative_cells_vq.append(imagescope_negative_cells)
    imagescope_total_positive = len(scaled_positive_dict[tile])
 
    imageJ_total_cells = len(ImageJ_dic[tile])
    imageJ_total_cells_vq.append(imageJ_total_cells)
    imageJ_positive = len(ImageJ_positive_dic[tile])
    imageJ_positive_vq.append(imageJ_positive)
    imageJ_negative = imageJ_total_cells - imageJ_positive
    imageJ_negative_vq.append(imageJ_negative)
    imageJ_score = (imageJ_positive / imageJ_total_cells)*100
    imageJ_score_vq.append(imageJ_score)

    
    positive_scope_d_vq.append((tile,
    imagescope_true_positive_detection, 
    imagescope_false_positive_detection, 
    imagescope_missing_positive,  
    imagescope_total_cells, 
    imagescope_total_positive,
    imagescope_negative_cells,
    imageJ_total_cells,
    imageJ_positive,
    imageJ_negative,
    imageJ_score))
     
ImageScope_positive_cell_count_vq = pd.DataFrame(positive_scope_d_vq, columns=('Slide',
    'imagescope true positive detection', 
    'imagescope false positive detection', 
    'imagescope missing positive cells',  
    'imagescope total cells', 
    'imagescope total positive cells',
    'imagescope total negative cells',
    'ImageJ total cells',
    'ImageJ positive cells',
    'ImageJ negative cells',
    'ImageJ score'))

ImageScope_positive_cell_count_vq.to_excel('vq-assignment-visualization/ImageScope_ImageJ_PositiveCell_vq.xlsx', index = False)

################################################################################################################