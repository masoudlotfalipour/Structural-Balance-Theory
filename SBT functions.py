# RESTING STATE

def get_functional_connectivity(subject, data_dir):

  # Load timeseries
  timeseries = load_timeseries(subject=subject,
                               name= task,
                              dir=data_dir,
                               runs=1)

  # Calculate correlation matrix
  fc = np.corrcoef(timeseries)

  return fc

########################################################
import numpy as np

def count_fc_links(fc_matrix):

  # Set diagonal to 0
  np.fill_diagonal(fc_matrix, 0)

  # Count number of positive links
  num_pos = np.sum(fc_matrix > 0)

  # Count number of negative links
  num_neg = np.sum(fc_matrix < 0)

  return num_pos, num_neg

########################################################
def count_triads(conn):
    n_nodes = conn.shape[0]
    n_balanced = 0
    n_imbalanced = 0

    for i in range(n_nodes):
        for j in range(i+1, n_nodes):
            for k in range(j+1, n_nodes):
                triad = np.sign(conn[i,j]) * np.sign(conn[j,k]) * np.sign(conn[k,i])
                if triad > 0:
                    n_balanced += 1
                elif triad < 0:
                    n_imbalanced += 1

    return n_balanced, n_imbalanced

#######################################################
def count_triad_links(conn):

  n_0neg = 0
  n_1neg = 0
  n_2neg = 0
  n_3neg = 0
  n_nodes = fc.shape[0]

  for i in range(n_nodes):
    for j in range(i+1, n_nodes):
      for k in range(j+1, n_nodes):

        if i == j or i == k or j == k:
          continue

        n_neg = 0
        if conn[i,j] < 0:
          n_neg += 1
        if conn[j,k] < 0:
          n_neg += 1
        if conn[k,i] < 0:
          n_neg += 1

        if n_neg == 0:
          n_0neg += 1
        elif n_neg == 1:
          n_1neg += 1
        elif n_neg == 2:
          n_2neg += 1
        elif n_neg == 3:
          n_3neg += 1

  return n_0neg, n_1neg, n_2neg, n_3neg

######################################################
def calculate_TMH(conn_matrix):
  # Calculate degrees (weighted degrees)
  degree = np.sum(conn_matrix, axis=1)
  pos_degree = np.sum(np.maximum(conn_matrix, 0), axis=1)  # Positive weighted degree
  neg_degree = np.sum(np.minimum(conn_matrix, 0), axis=1)  # Negative weighted degree

  # Calculate TMH
  TMH = np.sum(degree**2) / np.sum(degree)

  # Calculate positive and negative TMH
  posTMH = np.sum(pos_degree**2) / np.sum(np.maximum(degree, 1))  # Ensure denominator is not zero
  negTMH = np.sum(neg_degree**2) / np.sum(np.maximum(degree, 1))  # Ensure denominator is not zero

  output = np.round([TMH, posTMH, negTMH], 3)
  output_names = ['TMH', 'posTMH', 'negTMH']

  return dict(zip(output_names, output))

  ####################################################
#Task

def calculate_functional_connectivitytask(subject, task_name, condition_name, data_dir):
    # Load timeseries data and event information
    timeseries = load_timeseries(subject=subject, name=task_name, dir=data_dir, concat=False)
    evs = load_evs(subject=subject, name=task_name, condition=condition_name, dir=data_dir)
    frames = condition_frames(evs)

    task_data = []
    for i, (run_ts, run_frames) in enumerate(zip(timeseries, frames)):
      # Trim run_frames to only include valid indices
      run_frames = run_frames[run_frames < run_ts.shape[1]]
      task_data.append(run_ts[:, run_frames])

    # Concatenate the time series data across runs
    task_data = np.concatenate(task_data, axis=-1)



    # Calculate the functional connectivity matrix
    fc = np.corrcoef(task_data)


    return fc

########################################################
import numpy as np


def count_fc_links(fc_matrix):

  # Set diagonal to 0
  np.fill_diagonal(fc_matrix, 0)

  # Count number of positive links
  num_pos = np.sum(fc_matrix > 0)

  # Count number of negative links
  num_neg = np.sum(fc_matrix < 0)

  return num_pos, num_negd
########################################################
def count_triads(functional_connectivity_matrix):

  n_nodes = functional_connectivity_matrix.shape[0]
  n_balanced = 0
  n_imbalanced = 0

  for i in range(n_nodes):
    for j in range(i+1, n_nodes):
      for k in range(j+1, n_nodes):
        triad = np.sign(functional_connectivity_matrix[i,j]) * np.sign(functional_connectivity_matrix[j,k]) * np.sign(functional_connectivity_matrix[k,i])
        if triad > 0:
          n_balanced += 1
        elif triad < 0:
          n_imbalanced += 1

  return n_balanced, n_imbalanced

#######################################################
def count_triad_links(functional_connectivity_matrix):

  n_0neg = 0
  n_1neg = 0
  n_2neg = 0
  n_3neg = 0

  n_nodes = functional_connectivity_matrix.shape[0]

  for i in range(n_nodes):
    for j in range(i+1, n_nodes):
      for k in range(j+1, n_nodes):

        if i == j or i == k or j == k:
          continue

        n_neg = 0
        if functional_connectivity_matrix[i,j] < 0:
          n_neg += 1
        if functional_connectivity_matrix[j,k] < 0:
          n_neg += 1
        if functional_connectivity_matrix[k,i] < 0:
          n_neg += 1

        if n_neg == 0:
          n_0neg += 1
        elif n_neg == 1:
          n_1neg += 1
        elif n_neg == 2:
          n_2neg += 1
        elif n_neg == 3:
          n_3neg += 1

  return n_0neg, n_1neg, n_2neg, n_3neg

######################################################
def calculate_TMH(functional_connectivity_matrix):

  # Calculate degrees
  degree = np.sum(functional_connectivity_matrix, axis=1)
  pos_degree = np.sum(np.maximum(functional_connectivity_matrix, 0), axis=1)
  neg_degree = np.sum(np.minimum(functional_connectivity_matrix, 0), axis=1)

  # Calculate TMH
  TMH = np.sum(degree**2) / np.sum(degree)

  # Calculate positive and negative TMH
  posTMH = np.sum(pos_degree**2) / np.sum(np.maximum(degree, 1))
  negTMH = np.sum(neg_degree**2) / np.sum(np.maximum(degree, 1))

  output = np.round([TMH, posTMH, negTMH], 3)
  output_names = ['TMH', 'posTMH', 'negTMH']

  return dict(zip(output_names, output))

  ####################################################
