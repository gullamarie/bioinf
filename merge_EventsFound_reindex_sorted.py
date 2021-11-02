"""
This script merges multiple output files from EventPointer
and creates one results file of called Events 

Event merging implementation based on theoretical basis per EP paper described in 
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6156849/



Author: Marie Gulla 

"""

#Import packages 
import os
import pandas as pd
import numpy as np
import multiprocessing
import itertools
#import multiprocessing_import_worker
#from multiprocessing import Process, Manager
from itertools import chain, combinations
from collections import Counter

# Set working dir
os.chdir('eventpointer\\sets')

#OUTPUT FORMAT 
#all_events = pd.DataFrame(columns=['event_ID', 'gene', 'event number', 'event type', 'genome position', 'P1', 'P2', 'Ref'])

######################################################################################
######################################################################################

# HELPER FUNCTIONS

######################################################################################
######################################################################################


###########################################
def all_subsets(ss):
    return chain(*map(lambda x: combinations(ss, x), range(0, len(ss)+1)))
###########################################

###########################################
def is_overlapping(list1, list2):
    """
    Test whether two path lists have overlapping sets (path numbers)
    Format:
        list1/list2 = [[1, 1], [1, 2], [2, 2], [5, 5], [11, 11]]
        Nested lists 
    """

    list1_unique = np.unique(np.array(list1).flatten())
    list1_set = set(list1_unique)
    list2_unique = np.unique(np.array(list2).flatten())
    list2_set = set(list2_unique)
    
    if ( len(list1_set.intersection(list2_set)) > 0 ) or ( len(list2_set.intersection(list1_set)) > 0 ):
        return True
    
    return False

###########################################

###########################################
def is_overlapping_path(a_path, b_path):
    patha_list = a_path.split(',')
    pathb_list = b_path.split(',')

    common_elements = 0
    len_a = len(patha_list)
    len_b = len(pathb_list)

    longest_idx = [len_a, len_b].index(max(len_a, len_b))
    
    if longest_idx == 1:
        for i in range(0,len(patha_list)):
            if patha_list[i] in pathb_list:
                common_elements += 1


    else: 
        for i in range(0,len(pathb_list)):
            if pathb_list[i] in patha_list:
                common_elements += 1

    if common_elements > 0: 
        return True 
    
    return False

###########################################


def compatible_events(event1, event2):
    """
    event1: df format with col names 
    event2: lsit format 
    """

    condition1 = False
    condition2 = False
    condition3 = False

    #all three conditions must be true to return true 

    #test1
    
    # go through all paths in events, check if is the same event
    A1 = event1['Path 1']
    A2 = event1['Path 2']
    Aref = event1['Path Reference']

    B1 = event2[-3]
    B2 = event2[-2]
    Bref = event2[-1]

    #condition 1
    # A1 contained in A2
    pathA = collapse_path([A1, B1])
    pathB = collapse_path([A2, B2])
    pathRef = collapse_path([Aref, Bref])
    

    #test 1:
    if is_contained(pathA[0], pathA[1]) or is_contained(pathA[1], pathA[0]):
        condition1 = True
    if is_contained(pathB[0], pathB[1]) or is_contained(pathB[1], pathB[0]):
        condition2 = True
    if is_overlapping(pathRef[0], pathRef[1]) or is_overlapping(pathRef[1], pathRef[0]):
        condition3 = True

    if condition1 and condition2 and condition3:
        return True

    #test 2: 
    if is_contained(pathA[0], pathB[1]) or is_contained(pathB[1], pathA[1]):
        condition1 = True
    if is_contained(pathB[0], pathA[1]) or is_contained(pathA[1], pathB[0]):
        condition2 = True
    if is_overlapping(pathRef[0], pathRef[1]) or is_overlapping(pathRef[1], pathRef[0]):
        condition3 = True

    if condition1 and condition2 and condition3:
        return True
    

    return False

###########################################
def is_contained(a_path, b_path):
    """
    a_path format: [[3, 4], [4, 4], [5, 6]]
    b_path format: [[6, 6], [6, 7]]
    """

    #store bool return arg 
    is_contained = False

    #test for wiggle -1, 0, +1
    for i in range(0,3):
        i = i - 1 #print(i-1)

        a_path_wiggle = a_path
        for j in range(0,len(a_path_wiggle)):
            a_path_wiggle[j] = [numb + i for numb in a_path_wiggle[j] ]

        #is a_path_wiggle subset of b_path? 
        if checkSubset(b_path, a_path_wiggle):
            is_contained = True
         

    return is_contained

###########################################

###########################################
def checkSubset(list1, list2):
    l1, l2 = list1[0], list2[0]
    exist = True
    for i in list2:
        if i not in list1:
            exist = False
    return exist
###########################################


###########################################
# Function to insert row in the dataframe
def Insert_row(row_number, df, row_value):
    # Starting value of upper half
    df1 = df.copy()
    start_upper = 0
  
    # End value of upper half
    end_upper = row_number
  
    # Start value of lower half
    start_lower = row_number
  
    # End value of lower half
    end_lower = df1.shape[0]
  
    # Create a list of upper_half index
    upper_half = [*range(start_upper, end_upper, 1)]
  
    # Create a list of lower_half index
    lower_half = [*range(start_lower, end_lower, 1)]
  
    # Increment the value of lower half by 1
    lower_half = [x.__add__(1) for x in lower_half]
  
    # Combine the two lists
    index_ = upper_half + lower_half
  
    # Update the index of the dataframe
    df1.index = index_
  
    # Insert a row at the end
    df1.loc[row_number] = row_value
   
    # Sort the index labels
    df1 = df1.sort_index()
  
    # return the dataframe
    return df1

###########################################

###########################################
#FUNCTION CONTAINED IN 

def contained_in(cand_path, draft_path): 
    #checks if cand_path is subst of draft path 
    cand_path_np = np.array(cand_path).flatten()
    cand_path_np = np.unique(cand_path_np)
    cand_len = len(cand_path_np)
    draft_path_np = np.array(draft_path).flatten()
    draft_path_np = np.unique(draft_path_np)
    draft_len = len(draft_path_np)

    if sum([i in draft_path_np for i in cand_path_np]) == len(cand_path_np): 
        #all elements of cand path in draft 
        return True

  
    return False
    
###########################################
# FUNCTION COLLAPSE PATH 

def collapse_path(a_path):
    """
    a_path format: ['11-11,11-12', '10-11,11-11,11-12']
    collapse list of path strings into nested list 
    
    """

    for j in range(0,len(a_path)): 
        el = a_path[j]
        el = el.split(',') #split by comma
        el = [e.split('-') for e in el] #split by hyphen 
        
        #convert to ints 
        for k in range(0,len(el)):
            numbers = el[k]
            #print('numbers', numbers)
            for l in range(0,len(numbers)):
                numbers[l] = int(numbers[l])
            #print('numbers', numbers)
            el[k] = numbers
        
        #print(el)
        a_path[j] = el

    return a_path

###########################################

def collapse_nested_lists(a_path):
    out_list = []
    for i in range(0,len(a_path)):
        element = a_path[i]
        if element not in out_list:
            out_list.append(element)

    return out_list

###########################################
# FUNCTION MERGE PATH 

def merge_path(a_path):
    '''
    a_path format: nested list 
    '''

    max_len = 0
    a_path_new = []

    for i in range(0, len(a_path)):
        #print(a_path[i])
        el = a_path[i]
        el_len = len(el)
        
        #store largest value of inner paths per path to max_len 
        if el_len > max_len: 
            max_len = el_len
        if i < (len(a_path)-1):
            if not el == a_path[i+1]:
                a_path_new.append(a_path[i])
        else: 
            a_path_new.append(a_path[i])

    #begin with longest path as draft
    start_index = max(enumerate(a_path_new), key = lambda tup: len(tup[1]))[0]
    draft = a_path_new[start_index]

    #quick_merge? Do if : if all sublists in a_path are equal: use one 
    same_as_draft = 0
    for i in range(0,len(a_path)):
        sublist = a_path[i]
        if sublist == draft:
            same_as_draft += 1
    if same_as_draft > 1: 
        return draft
    
    #loop over all other paths to merge with draft 
    #for i in range(1,len(a_path_new)-1):
    for i in range(0,len(a_path_new)):
        #print('i=', i)
        candidate_path = a_path_new[i]
        if draft == candidate_path:
            continue 
        if draft != candidate_path:
            
            #if draft and cand are nested lists 
            
            draft_flat = [item for items in draft for item in items]
            candidate_flat = [item for items in candidate_path for item in items]
            
            #combine or merge? 
            
            #combine paths ?
            #find overlap 
            #common_vals = number of common sublists in nested lists to merge 
            commons = 0
            for val in candidate_flat:
                if val in draft_flat: 
                    commons += 1

            #if number of common values are more than 50%, combine events instead of force merge 
            if commons/len(candidate_flat) > 0.5:
                #new draft:
                new_draft = draft.copy() 
                for sublist in candidate_path:
                    if sublist not in draft:
                        #add sublist 
                        new_draft.append(sublist)

                draft = sort_nested_list(new_draft)
        

            #merge paths
            else: #number of common values are less than 50%, force merge events 
                    
                #find which end to get delta adjuster 
                back_delta = draft_flat[-1] - candidate_flat[-1]
                back_delta_abs = abs(back_delta)
                front_delta = draft_flat[0] - candidate_flat[0]
                front_delta_abs = abs(front_delta)

                #adjust to the side where delta is smallest 
                #adjust to first value 
                if front_delta_abs < back_delta_abs:
                    delta = front_delta 
                    other_delta = back_delta
                #adjust to last value 
                else: 
                    delta = back_delta
                    other_delta = front_delta

                #safety check: no elements can become negative 
                #if first element of candidate - delta becomes negative, use the other delta 

                candidate_path_delta_new = [e + delta for e in candidate_flat]
                if np.any((np.array(candidate_path_delta_new) < 0)):
                    #test if other delta works 
                    test_other_delta = [e + other_delta for e in candidate_flat]
                    if not np.any((np.array(test_other_delta) < 0)):
                        #use other delta 
                        delta = other_delta
                    else: 
                        print('delta adjustment makes candidate path elements negative')
                        print(candidate_path)
                        print('a_path:', a_path)

                #adjust candidate path with delta 

                cand_delta = candidate_path.copy()
                #loop over subpaths in paths 
                for j in range(0,len(cand_delta)):
                    #print('j', j)
                    #add delta to all elements of draft 
                    values = cand_delta[j]
                    #print(values)
                    for k in range(0,len(values)):
                        number = values[k] + delta
                        #print(number)
                        values[k] = number
                    cand_delta[j] = values

                if contained_in(cand_delta, draft):
                    #cand_delta is part of draft, keep draft  
                    pass 
                    #break? 
                else: 
                    
                    #find elements in cand_delta not in draft 
                    unique_elms = []
                    for i in range(0, len(cand_delta)):
                        element = cand_delta[i]
                        #print(element)
                        if not element in draft:
                            unique_elms.append(element)
                    
                    #add these values to draft
                    draft.extend(unique_elms)
                    draft.sort()

    return draft 
            
###########################################

###########################################

def sort_nested_list(a_list):
    """
    Sorts a nested list numerically 
    format: a_list = [[3,3], [4,4], [1,1], [1,2]]
    out : a_list = [[1,1], [1,2],[3,3], [4,4]]

    """
    sorted_list = sorted(a_list, key=lambda x: x[0])


    return sorted_list

###########################################
def force_merge_path(in_event):
    '''
    format in_event: splice from pd.DataFrame 

    '''
    a_event_draft = in_event.copy()
    a_event_draft = a_event_draft.iloc[0,]

    #path1 
    path1_list = list(in_event["Path 1"].values)
    path1 = collapse_path(path1_list)
    #collapse does not work, change into collapsable format 
    path1.sort()
    path1 = collapse_nested_lists(path1)
    path1 = merge_path(path1) 
    path1 = collapse_nested_lists(path1)
    path1 = nested_list_to_string(path1)

    #path2 
    path2_list = list(in_event["Path 2"].values)
    path2 = collapse_path(path2_list)
    #collapse does not work, change into collapsable format 
    path2.sort()
    path2 = collapse_nested_lists(path2)
    path2 = merge_path(path2) 
    path2 = collapse_nested_lists(path2)
    path2 = nested_list_to_string(path2)

    #pathref 
    pathref_list = list(in_event["Path Reference"].values)
    pathref = collapse_path(pathref_list)
    #collapse does not work, change into collapsable format 
    pathref.sort()
    pathref = collapse_nested_lists(pathref)
    pathref = merge_path(pathref) 
    pathref = collapse_nested_lists(pathref)
    pathref = nested_list_to_string(pathref)

    a_event = a_event_draft.copy()

    a_event.iloc[5,] = path1
    a_event.iloc[6,] = path2
    a_event.iloc[7,] = pathref

    return a_event


###########################################
#converts nested list of ints to list of strings 

def nested_list_to_string(in_path):
    out_path = []
    for j in range(0, len(in_path)):
        numbers = in_path[j]
        numb_str_list = [str(e) for e in numbers] 
        numb_str = '-'.join(numb_str_list)
        out_path.append(numb_str)
    #print(out_path)

    #join list of paths to comma-sep string 
    out_path = ','.join(out_path)

    return out_path
###########################################

def nested_list_to_string2(in_path):
    out_path = []
    for j in range(0, len(in_path)):
        numbers = in_path[j]
        subpath = []
        for k in range(0, len(numbers)):
            values = numbers[k]
            if type(values) == list:
                values = [str(e) for e in values] 
            path = '-'.join(values)
            subpath.append(path)
        subpath = ','.join(subpath)
        out_path.append(subpath)

    return out_path
###########################################


def flatten_nested_list(a_list):
    flat_list = []
    for sublist in a_list:
        for item in sublist:
            flat_list.append(item)
    return flat_list

######################################################################################
######################################################################################
#end helper functions 


######################################################################################
######################################################################################

# MAIN CODE 

######################################################################################
######################################################################################


#list all samples 
samples = ['EventsFound_RNASeq_set1.txt', 'EventsFound_RNASeq_set2.txt',
    'EventsFound_RNASeq_set3.txt','EventsFound_RNASeq_set4.txt',
    'EventsFound_RNASeq_set5.txt', 'EventsFound_RNASeq_set6.txt', 'EventsFound_RNASeq_set7.txt']

#run loop from here: 

#number of data set results to merge 
samples_n = len(samples)
event_dfs = []

#store all files to event_dfs 
for i in range(samples_n): #sample in samples:
    print(i)
    event_df = pd.read_csv(samples[i], sep= '\t')
	#file = "%s-.txt" % (conn[sample])
    
    event_dfs.append(event_df)

# Make dfs with all events in dfs 

#start with first df 	
all_events = event_dfs[0]
# number of events in first df from set1-file: 69988

#store all events to be matched in matching dict 
#gen_pos should be keys 
matching_dict={}
all_gene_positions = all_events["Genomic Position"].tolist()

#store all present genes 
all_genes = all_events["Gene"].tolist()

#loop throgh all dfs and add new events from them
for i in range(1,samples_n):


    #loop through all dfs
    df_in = event_dfs[i]
    nrows = len(df_in)
    print("Processing df nr:", i)

    #loop through all rows/events of df
    for j in range(0,nrows):

        new_event = df_in.iloc[j,]
        new_event_list = new_event.tolist()
        gene = df_in.iloc[j,]["Gene"]
        gen_pos = df_in.iloc[j,]["Genomic Position"]
        pos = gen_pos.split('chr')[1]
        

        #exact genomic pos already in df, match events 
        if gen_pos in all_gene_positions:
            
            #get event already in dict 
            try: 
                event_in_all_events = all_events[all_events["Genomic Position"]== gen_pos].values.tolist()[0]
            except IndexError:
                print(all_events[all_events["Genomic Position"]== gen_pos])
                print(j)

            if not compatible_events(new_event, event_in_all_events):

                insert_idx = all_events[all_events["Genomic Position"] == gen_pos].index[0]

                all_events = Insert_row(insert_idx, all_events, new_event)
                all_genes.extend(gene)
                

            #  compatible
            if compatible_events(new_event, event_in_all_events):
                new_event_to_dict = [new_event_list[1]] + new_event_list[3:]
                
                if gen_pos in matching_dict.keys():
                    
                    matching_dict[gen_pos].append(new_event_to_dict)
                    #matching_dict[gen_pos] = new_values

                #genpos not in matching dict: make new dict item 
                else: 
                    
                    #if old event has gene name 
                    if 'ENSG' in event_in_all_events[0]: 
                    
                        #add info as nested list 
                        new_list_to_add = [0,0]
                        new_list_to_add[0] = event_in_all_events
                        new_list_to_add[1] = new_event_to_dict
                        #add both to matching dict 
                        matching_dict[gen_pos] = new_list_to_add

                    else: 

                        if 'ENSG' in new_event_to_dict[0]:
                            new_event_to_dict[0] = new_event_to_dict[0]
                            insert_idx = all_events[all_events["Genomic Position"] == gen_pos].index[0]

                            #update info in dict, give gene name 
                            all_events.iloc[insert_idx,1] = new_event_to_dict[0]
                            all_genes.extend(gene)
                            
                        
                        #add info as nested list 
                        new_list_to_add = [0,0]
                        new_list_to_add[0] = event_in_all_events
                        new_list_to_add[1] = new_event_to_dict
                        matching_dict[gen_pos] = new_list_to_add
            
            #compatible condtion end 

                
        #gen_pos not in df already 
        else: #new event to be added 
            if gene in all_genes:
                indeces = all_events[all_events["Gene"] == gene].index

                if len(indeces) == 0:
                    indeces = all_events[~all_events["Gene"].str.contains('ENSG')].index
                
                last_index = indeces[-1] + 1 
                #inserts new event into all_events AFTER the last item of the same gene 

                all_events = Insert_row(last_index, all_events, new_event)
                
                #update lists of gene names and gene positions 
                all_genes.extend(gene)
                all_gene_positions.extend(gen_pos)

            else: 
                #insert at first row without gene name
                first_gene_idx = ['ENSG' in gene for gene in all_events["Gene"].tolist()].index(True)
                #all_events = Insert_row(first_gene_idx, all_events, df_in.iloc[j,])
                all_events = Insert_row(first_gene_idx, all_events, new_event)
                
                #update lists of gene names and gene positions 
                all_genes.extend(gene)
                all_gene_positions.extend(gen_pos)

        #if j is n*1'000 print progress j 
        if (j % 1000 == 0):
            print("At line nr.:", j)

        
    #end j loop 
#end i loop 
#####################################

#takes 108300 s to run ~ 30h CPUh 

#make copies to work with: 
matching_dict2 = matching_dict.copy()
all_events2 = all_events.copy()
#len all_events 263546
#len matching_dict 15958


#####################################

# Merge events in matching dict to events in df 

# test version: 
key_list = list(matching_dict2.keys())#[0:30]
for i in range(0,len(key_list)):
    #print(key_list[i])
    key = key_list[i]
    events = matching_dict2[key]
    n = len(events)

    events_flat = flatten_nested_list(events)
    events_flat = [str(e) for e in events_flat]
    if any('ENSG' in substring for substring in events_flat):
        gene_name = ''
        for i in range(0,len(events_flat)):
            if 'ENSG' in events_flat[i]:
                gene_name = events_flat[i]
                break
       
    #PROCESS MERGE ALL EVENTS IN MATCHING DICT 
    pathA = []
    pathB = []
    pathRef = []
    for event in events:
        pathA.append(event[-3])
        pathB.append(event[-2])
        pathRef.append(event[-1])
    
    #collapse paths to lists of ints 
    pathA = collapse_path(pathA)
    pathB = collapse_path(pathB)
    pathRef = collapse_path(pathRef)

    #merge paths into one 
    pathA = merge_path(pathA)
    pathB = merge_path(pathB)
    pathRef = merge_path(pathRef)
    match_path = [pathA, pathB, pathRef]

    # FIND OLD PATH IN DF TO MATCH WITH 
    ##################### fix here: if empty 
    old_event = all_events2[all_events2["Genomic Position"] == key]
    old_event_idx = all_events2[all_events2["Genomic Position"] == key].index

    #print(len(old_event))
    
    if len(old_event) < 2 : 
        #if old event is not more than one element
        old_event_list = old_event.values.tolist()
        if len(old_event_list) == 1:
            old_event_list = old_event.values.tolist()[0]

        #print(len(old_event_list))
        old_event_list_short = [old_event_list[1]] + old_event_list[3:]
        old_event_idx = old_event_idx[0]
    
    # if old event is more than one element : match with gene as well 
    elif len(old_event) > 1:
        old_event = all_events2[(all_events2["Genomic Position"] == key) & (all_events2["Gene"] == gene_name)]
        #old_event[old_event["Gene"] == gene_name]
        old_event_idx = old_event.index #[0]?

        #is empty here 
        if len(old_event) == 0:
            #find INTERNAL common gene name 
            old_event = all_events2[all_events2["Genomic Position"] == key]
            old_event_idx = all_events2[all_events2["Genomic Position"] == key].index

            #find non-unique gene name in more than one event 
            genes = old_event["Gene"].to_list()
            counts = Counter(genes)
            count = counts.most_common()[0][1]
            if count > 1:
                most_common_gene = counts.most_common()[0][0]
                #one gene is most frequent, choose this 
                old_event = all_events2[(all_events2["Genomic Position"] == key) & (all_events2["Gene"] == most_common_gene)]
                old_event_idx = old_event.index
            
        
        elif len(old_event) > 1: 
            any_compatible = False
            #get all combinations of two events 
            index_pairs = []
            for subset in all_subsets(old_event.index):
                if len(subset) == 2:
                    index_pairs.append(subset)

            #loop through all combinations of index-pairs to find a compatible pair 
            for j in range(0,len(index_pairs)):
                indeces = index_pairs[j]
                
                event1 = old_event.loc[indeces[0]]
                event2 = old_event.loc[indeces[1]]

                if not compatible_events(event1, event2):
                    continue
                    ############

                elif compatible_events(event1, event2):
                    any_compatible = True
                    #print(event1, event2)

                    #merge compatible events : 
                    #store paths 
                    path_one = [event1["Path 1"], event2["Path 1"]]
                    path_two = [event1["Path 2"], event2["Path 2"]]
                    path_ref = [event1["Path Reference"], event2["Path Reference"]]
                    
                    #collapse paths 
                    path_one = collapse_path(path_one)
                    path_two = collapse_path(path_two)
                    path_ref = collapse_path(path_ref)
                    
                    #merge paths 
                    path_one = merge_path(path_one)
                    path_two = merge_path(path_two)
                    path_ref = merge_path(path_ref)

                    #convert nested list to string again 
                    path_one = nested_list_to_string(path_one)
                    path_two = nested_list_to_string(path_two)
                    path_ref = nested_list_to_string(path_ref)

                    #make merged event 
                    old_event = event1.copy()
                    old_event["Path 1"] = path_one
                    old_event["Path 2"] = path_two
                    old_event["Path Reference"] = path_ref

                    #old_event is updated
                    old_event_idx = indeces[0]
                    old_event_list = old_event_list = old_event.values.tolist()

                    #remove event2
                    #safety test: 
                    if (all_events2.iloc[indeces[1], 4] == key) & (all_events2.iloc[indeces[1], 1] == gene_name):
                        #all_events2.iloc[indeces[1], "Gene"] == gene_name and genpos = gen_pos
                        all_events2.drop(indeces[1],  axis = 0, inplace = True)

                    break #found compatible pair, stop looking 

            #no combinations of events in old_events were compatible 
            if not any_compatible:  
                # treat events as separate events 
                pathA_string = nested_list_to_string(match_path[0])
                pathB_string = nested_list_to_string(match_path[1])
                pathRef_string = nested_list_to_string(match_path[2])

                match_path_inverse = [pathRef_string, pathB_string, pathA_string]

                #count matches per event and store in event_matches
                event_matches = np.zeros((len(old_event)))
                for k in range(0,len(old_event)):
                    ev = old_event.iloc[k]
                    #count number of overlapping paths for this event 
                    overlaps = 0
                    for l in range(-3,0):
                        old_ev_path = ev.iloc[l]
                        match_ev_path = match_path_inverse[-l-1]
                        if is_overlapping_path(old_ev_path, match_ev_path):
                            overlaps += 1
                    
                    #store number of overlapping paths for this event 
                    event_matches[k] = overlaps
                        
                idx = np.argmax(event_matches)

                old_event_new = old_event.iloc[idx] #.copy()
                old_event_idx = old_event_idx[idx]
                old_event = old_event_new.copy()
                old_event_list = old_event.values.tolist()
            
        #  event was more than one element, is now one after gene matching 
        else: 
            old_event_list = old_event.values.tolist()[0]
            
        #save old_list_short version 
        if len(old_event_list) < 4: 
            print(key)
            print(old_event)
        old_event_list_short = [old_event_list[1]] + old_event_list[3:]

    

    # MERGING OLD AND NEW PATH  
    old_path = collapse_path(old_event_list[5:8])
    new_path = []
    for i in range(0,3):
        path_to_merge = [match_path[i], old_path[i]]
        path_merged = merge_path(path_to_merge)
        new_path.append(path_merged)
    
    #convert ints to string in new path 
    new_path_string = nested_list_to_string2(new_path)
    

    all_events2.at[old_event_idx,"Path 1"] = new_path_string[0]
    all_events2.at[old_event_idx,"Path 2"] = new_path_string[1]
    all_events2.at[old_event_idx,"Path Reference"] = new_path_string[2]

#end main for loop 
###########################################################################

#paths in row of df now consists of merged overlapping paths from multiple sets collapsed to one, 
# and merged with the existing paths for that genomic position in df 


#write out as text file 
#make sure path columns are not quoted 


all_events2.to_csv('EvendFound-all_non_quot.csv', encoding='utf-8', index=True, quoting = csv.QUOTE_NONE)

np.savetxt('EvendFound-all_non_np.txt', all_events2.values, fmt = "%s") 


###########################################################################
#Sort df by chr pos 
###########################################################################

###########################################################################
# Make safety copy of all_events2

###########################################################################
all_events4 = all_events2.copy()
###########################################################################
#all_events4.sort_values(['Genomic Position'])
#is sorted by gen pos 


