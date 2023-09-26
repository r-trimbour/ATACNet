import pandas as pd
import numpy as np
import anndata as ad
import scipy as sp
import matplotlib.pyplot as plt
import tqdm

# Create fake single-cell atac-seq data
counts = pd.DataFrame(np.random.randint(0, 100, size=(50, 200)),
                      index=['Cell_'+i for i in map(str, range(50))],
                      columns=['chr1_'+i+'_'+i for i in map(str, range(200))])

atac = ad.AnnData(counts)


def add_region_infos(AnnData,
                     sep=('_', '_'),
                     inplace=True):
    """
    Get region informations from the var_names of AnnData object.
    e.g. chr1_12345_12346 -> 'chromosome' : chr1, 'start' : 12345, 'end' : 12346
    These info will be added to var of AnnData object.
        adata.var['chromosome'] : chromosome
        adata.var['start'] : start position
        adata.var['end'] : end position

    Parameters
    ----------
    AnnData : AnnData object
        AnnData object with var_names as region names.
    sep : tuple, optional
        Separator of region names. The default is ('_', '_').
    
    Returns
    -------
    AnnData : AnnData object
        AnnData object with region informations in var.
    """
    # Check if user wants to modify AnnData inplace or return a copy
    if inplace:
        pass
    else:
        AnnData = AnnData.copy()
    regions_list = AnnData.var_names

    # Replace sep[1] with sep[0] to make it easier to split
    regions_list = regions_list.str.replace(sep[1], sep[0])

    # Split region names
    regions_list = regions_list.str.split(sep[0]).tolist()

    # Check if all regions have the same number of elements
    if set([len(i) for i in regions_list]) != set([3]):
        raise ValueError("""Not all regions have the same number of elements.
                         Check if sep is correct, it should be ({}, {}),
                         with only one occurence each in region names.""".format(sep[0], sep[1]))

    # Extract region informations from var_names
    region_infos = pd.DataFrame(regions_list,
                                index=AnnData.var_names,
                                columns=['chromosome', 'start', 'end'])

    # Convert start and end to int
    region_infos['start'] = region_infos['start'].astype(int)
    region_infos['end'] = region_infos['end'].astype(int)

    # Add region informations to var
    AnnData.var['chromosome'] = region_infos['chromosome']
    AnnData.var['start'] = region_infos['start']
    AnnData.var['end'] = region_infos['end']

    sort_regions(AnnData)
    # Return AnnData if inplace is False
    if inplace:
        pass
    else:
        return AnnData


def sort_regions(AnnData):
    """
    Sort regions by chromosome and start position.
    """
    AnnData.var.sort_values(['chromosome', 'start'], inplace=True)
    return AnnData


def get_distance_regions(AnnData, chromosomes=None):
    """
    Get distance between regions.
    """
    # Check if chromosomes is None
    if chromosomes is None:
        # Get chromosome list
        chromosomes = AnnData.var['chromosome'].unique().tolist()
    else:
        if not np.array([i in AnnData.var['chromosome'].unique() for i in chromosomes]).all():
            raise ValueError("""Chromosomes should be in AnnData.var['chromosome'].
                                Check if chromosomes is correct.""")

    # A dictionary to store distance between regions for each chromosome
    distances = {}

    # Get distance between regions for each chromosome
    for chromosome in chromosomes:
        chr_mask = AnnData.var['chromosome']==chromosome
        # Store start and end positions in two arrays
        m, n = np.meshgrid(AnnData.var['start'].values[chr_mask],
                           AnnData.var['end'].values[chr_mask])

        # Get distance between start of region m and end of region n
        distance = np.abs(m-n)
        # Substract length of the region to get distance between
        # end of region m and start of region n
        # a.k.a. distance between closest bases of two regions
        distance = (distance.T-(AnnData.var['end'].values[chr_mask]
                                - AnnData.var['start'].values[chr_mask])).T

        # Remove diagonal (distance between a region and itself)
        distance -= np.diag(distance)

        # Keep upper triangle of the distance matrix
        # (we don't want to calculate the same connection twice)
        distance = np.triu(distance, k=1)

        # Test if distance is negative
        if np.any(distance < 0):
            raise ValueError("""Distance between regions should be positive.
                            You might have overlapping regions.""")

        # Store distance in a dictionary
        distances[chromosome] = distance

    # Return distance
    return distances


def potential_connections(AnnData, threshold, chromosomes=None):
    """
    Get potential connections between regions based on distance.
    """
    # Check if chromosomes is None
    if chromosomes is None:
        # Get chromosome list
        chromosomes = AnnData.var['chromosome'].unique().tolist()
    else:
        if not np.array([i in AnnData.var['chromosome'].unique() for i in chromosomes]).all():
            raise ValueError("""Chromosomes should be in AnnData.var['chromosome'].
                                Check if chromosomes is correct.""")
    # Get distance between regions
    distances = get_distance_regions(AnnData, chromosomes=chromosomes)

    potential_connections = {}
    # Get potential connections
    for chromosome in chromosomes:
        print("Getting potential connections for chromosome {}...".format(chromosome))
        # Get potential connections for each chromosome
        distance = distances[chromosome]
        potential_chr_co = np.where((distance <= threshold)
                                         & (distance > 0))

        # Store potential connections in a dictionary
        potential_connections[chromosome] = potential_chr_co

    # Return potential connections
    return potential_connections


def corrcoef_connections(AnnData, potential_connections, as_sparse=True):
    """
    Get correlation coefficient between regions.
    """
    # Transform potential_connections into a sparse matrix
    potential_connections = global_sparse(AnnData, potential_connections)

    # Get correlation coefficient between regions
    corr_coefs = [np.corrcoef(AnnData.X[:, potential_connections.row[i]],
                              AnnData.X[:, potential_connections.col[i]])[0, 1]
                  for i in tqdm.tqdm(range(len(potential_connections.row)))]

    # Convert to sparse matrix if as_sparse is True
    if as_sparse:
        corr_coefs = sp.sparse.coo_matrix((corr_coefs,
                                           (potential_connections.row,
                                            potential_connections.col)),
                                          shape=(AnnData.shape[1],
                                                 AnnData.shape[1]))
    else:
        corr_coefs = np.array(corr_coefs)
    
    # Return correlation coefficients
    return corr_coefs




    # A dictionary storing the informations to create a sparse matrix of potential connections
    corr_coefs = {}
    corr_coefs['values'] = np.array([])
    corr_coefs['idx'] = np.array([])
    corr_coefs['idy'] = np.array([])

    # Get correlation coefficient between regions for each chromosome
    for chromosome in potential_connections.keys():
        chr_mask = AnnData.var['chromosome']==chromosome
        print("Getting correlation coefficient for chromosome {}...".format(chromosome))

        corr_coefs_chr = [np.corrcoef(AnnData.X[:,chr_mask][:, potential_connections[chromosome][0][i]],
                                      AnnData.X[:,chr_mask][:, potential_connections[chromosome][1][i]])[0, 1]
                          for i in range(len(potential_connections[chromosome][0]))]

        # Store correlation coefficient in a dictionary
        corr_coefs[chromosome] = corr_coefs_chr
    return corr_coefs


def global_sparse(AnnData, chr_idx, values=1):
    """
    Create a sparse matrix from a dictionary of np.where output on 'regions*regions' matrices.of different chromosomes.
    using global indices of an AnnData object.

    e.g.:
    AnnData.var_names = ['chr1_12345_12346', 'chr1_12347_12348', 'chr2_12345_12346', 'chr2_12347_12348']
    And we know values for  (chr1_12345_12346, chr1_12347_12348),
                            (chr1_12347_12348, chr1_12345_12346),
                            (chr2_12345_12346, chr2_12345_12346),
                            (chr2_12347_12348, chr2_12347_12348).
    
    We can create a sparse matrix with 'chr_dic' defined as:
    chr_dic = {'chr1' : np.array([[1, 0], [0, 1]]),
               'chr2' : np.array([[0, 0], [1, 1]])}

    sparse_mtx = global_sparse(AnnData, chr_dic)
    sparse_mtx
    'OUTPUT' :               chr_1_12345_12346  chr_1_12347_12348  chr_2_12345_12346  chr_2_12347_12348
    chr_1_12345_12346                 *                  1                  *                  *
    chr_1_12347_12348                 1                  *                  *                  *
    chr_2_12345_12346                 *                  *                  1                  *
    chr_2_12347_12348                 *                  *                  *                  1

    Parameters
    ----------
    AnnData : AnnData object
        AnnData object with var_names as region names.
    chr_dic : dictionary
        Dictionary of matrices (regions*regions) of different chromosomes.
    """

    data = {}
    data['values'] = np.array([])
    data['idx'] = np.array([])
    data['idy'] = np.array([])

    for chromosome in chr_idx.keys():
        chr_mask = AnnData.var['chromosome'] == chromosome
        # Get region names (needed to get global indices)
        indices_names = AnnData.var_names[chr_mask][chr_idx[chromosome][0]]
        columns_names = AnnData.var_names[chr_mask][chr_idx[chromosome][1]]
        map_indices = {AnnData.var_names[i]: i
                       for i in range(len(AnnData.var_names))}

        # Add global indices of potential connections
        data['idx'] = np.concatenate([data['idx'],
                                      indices_names.map(map_indices).values])
        data['idy'] = np.concatenate([data['idy'],
                                      columns_names.map(map_indices).values])
        # Add values of potential connections
        if type(values) == int:
            val = np.repeat(values, len(chr_idx[chromosome][0]))
        elif type(values) == dict:
            val = values[chromosome]
        else:
            raise ValueError("""values should be an int,
                             or a dict (of numpy arrays).""")
        data['values'] = np.concatenate([data['values'],
                                         val])
    # Create sparse matrix
    sparse_data = sp.sparse.coo_matrix((data['values'],
                                        (data['idx'],
                                         data['idy'])),
                                       shape=(AnnData.shape[1],
                                              AnnData.shape[1]))

    return sparse_data


    # A dictionary storing the informations to create a sparse matrix of potential connections
    potential_connections = {}
    potential_connections['corrcoef'] = np.array([])
    potential_connections['idx'] = np.array([])
    potential_connections['idy'] = np.array([])

    # Get potential connections
    for chromosome in chromosomes:
        print("Getting potential connections for chromosome {}...".format(chromosome))
 
        # Get distances
        distance = distances[chromosome]
        # Get region names (needed to get global indices)
        indices_names = AnnData.var_names[chr_mask][potential_chr_co[0]]
        columns_names = AnnData.var_names[chr_mask][potential_chr_co[1]]
        map_indices = {AnnData.var_names[i]:i for i in range(len(AnnData.var_names))}
        # Add global indices of potential connections
        potential_connections['idx'] = np.concatenate([potential_connections['idx'],
                                                       np.where(AnnData.var_names.isin(indices_names))[0]])
        potential_connections['idy'] = np.concatenate([potential_connections['idy'],
                                                       np.where(AnnData.var_names.isin(columns_names))[0]])

        # Get potential connections
        potential_chr_co = np.where((distance <= threshold)
                                    & (distance > 0))
        # Add values of potential connections
        potential_connections['corrcoef'] = np.concatenate([potential_connections['corrcoef'],
                                                            potential_chr_co])

    # Return potential connections
    return potential_connections









    # Get correlation coefficient between regions
    for i in range(len(potential_connections[0])):
        corr_coef = np.corrcoef(AnnData.X[:, potential_connections[0][i]],
                                AnnData.X[:, potential_connections[1][i]])[0, 1]
        corr_coefs.append(corr_coef)

    # Convert to sparse matrix if as_sparse is True
    if as_sparse:
        corr_coefs = sp.sparse.coo_matrix((corr_coefs,
                                           (potential_connections[0],
                                            potential_connections[1])),
                                          shape=(AnnData.shape[1],
                                                 AnnData.shape[1]))
    else:
        pass
    # Return correlation coefficients
    return corr_coefs


def sliding_graphical_lasso(AnnData,
                            scores,
                            penalties,
                            window_size,
                            ):
    """
    Extract sliding submatrix from a sparse correlation matrix.
    """
    start_slidings = [0, window_size/2]

    results = {}
    results['scores'] = np.array([])
    results['idx'] = np.array([])
    results['idy'] = np.array([])

    regions_list = AnnData.var_names
    # Get global indices of regions
    map_indices = {regions_list[i]: i for i in range(len(regions_list))}

    for k in start_slidings:
        for chromosome in AnnData.var['chromsome'].unique():
            # Get start positions of windows
            window_starts = [i for i in range(k,
                                              AnnData.var['end']
                                              [AnnData.var['chromosome']
                                                  == chromosome].max(),
                                              window_size)]

            for start in window_starts:
                end = start + window_size
                # Get global indices of regions in the window
                idx = np.where((AnnData.var['chromosome'] == chromosome)
                            & (AnnData.var['start'] >= start)
                            & (AnnData.var['start'] <= end))[0]

                # already global ?
                ## Get global indices of regions in the window
                #idx = [map_indices[i] for i in regions_list[idx]]

                # Get submatrix
                window_scores = sp.sparse.csr_matrix(scores)[idx, :][:, idx]
                window_scores = window_scores.todense()+window_scores.T.toarray()

                window_penalties = sp.sparse.csr_matrix(penalties)[idx, :][:, idx]
                window_penalties = window_penalties.todense()+window_penalties.T.toarray()

                # Initiating graphical lasso
                graph_lasso_model = quic_graph_lasso.QuicGraphicalLasso(init_method='precomputed', lam=window_penalties)

                # Fit graphical lasso
                graph_lasso_model.fit(window_scores)

                # Names of regions in the window
                window_region_names = AnnData.var_names[idx]

                # convert to sparse matrix the results
                corrected_scores = sp.sparse.coo_matrix(graph_lasso_model.covariance_)

                # Convert corrected_scores column and row indices to global indices
                idx = [map_indices[name]
                    for name in window_region_names[corrected_scores.row]]
                idy = [map_indices[name]
                    for name in window_region_names[corrected_scores.col]]

                # Add the "sub" resuls to the global sparse matrix
                results['scores'] = np.concatenate([results['scores'],
                                                            corrected_scores.data])
                results['idx'] = np.concatenate([results['idx'],
                                                idx])
                results['idy'] = np.concatenate([results['idy'],
                                                idy])

            # Create sparse matrix
            results['window_' + str(k)] = sp.sparse.coo_matrix((results['scores'],
                                            (results['idx'],
                                            results['idy'])),
                                            shape=(len(window_region_names),
                                                   len(window_region_names)))

    sliding_keys = ['window_' + str(k) for k in start_slidings]

    positive_coords = []
    negative_coords = []
    for k in sliding_keys:
        positive_coords.append({(x, y) for x, y, d in zip(results[k].row,
                                                          results[k].col,
                                                          results[k].data)
                                if d >= 0})
        negative_coords.append({(x, y) for x, y, d in zip(results[k].row,
                                                          results[k].col,
                                                          results[k].data)
                                if d <= 0})

    # Get coordinates intersection
    positive_coords = set.intersection(*positive_coords)
    negative_coords = set.intersection(*negative_coords)

    coords = pd.DataFrame(set.union((negative_coords, positive_coords)),
                          columns=['row', 'col'])

    # Get average of positive and negative coordinates
    average = sp.sparse.csr_matrix(sliding_keys[0])[coords['row'],
                                                    coords['col']]
    for k in sliding_keys[1:]:
        average += sp.sparse.csr_matrix(sliding_keys[0])[coords['row'],
                                                         coords['col']]
    average /= len(sliding_keys)

    return average
