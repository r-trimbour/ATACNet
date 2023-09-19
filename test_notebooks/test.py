import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt

# Create fake single-cell atac-seq data
counts = pd.DataFrame(np.random.randint(0,100,size=(50, 200)),
                    index=['Cell_'+i for i in map(str, range(50))],
                    columns=['chr1_'+i+'_'+i for i in map(str, range(200))])

atac = ad.AnnData(counts)

def add_region_infos(AnnData,
                     sep = ('_', '_'),
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

def group_potential_connections(AnnData, distant_constraint, chromosomes_sizes):
    """
    Group potential connections based on distance constraint.
    """

    # Create a copy of AnnData object
    AnnData = AnnData.copy()

    # Sort regions by chromosome and start position
    AnnData = sort_regions(AnnData)

    # Get chromosome list
    chromosomes = AnnData.var['chromosome'].unique().tolist()

    # Get group of regions contained in the same chromosome
    for chromosome in chromosomes:
        # Get regions in the same chromosome
        adata_chr[chromosome] = [atac[:, atac.var['chromosome']==chromosome] for chromosome in chromosomes]

        # Define windows based on distant_constraint for each chromosome
        for i in range(0, chromosomes_sizes[chromosome], distant_constraint):
            borders = [i, i+distant_constraint]
            if borders[1] > chromosomes_sizes[chromosome]:
                borders[1] = chromosomes_sizes[chromosome]
            # List of regions in the same window
            # (we accept regions that starts or end in the window)
            mask = (adata_chr[chromosome].var['start'] >= borders[0]) & (adata_chr[chromosome].var['start'] <= borders[1])
            mask += (adata_chr[chromosome].var['end'] >= borders[0]) & (adata_chr[chromosome].var['end'] <= borders[1])
            regions_in_window = adata_chr[chromosome] [mask]

            #  





    # Group potential connections
    potential_connections = pd.DataFrame(potential_connections, columns=['start', 'end'])
    potential_connections['group'] = np.nan
    group = 0
    for i in range(len(potential_connections)):
        if np.isnan(potential_connections['group'][i]):
            group += 1
            potential_connections['group'][i] = group
            for j in range(i+1, len(potential_connections)):
                if potential_connections['start'][j] >= potential_connections['start'][i] and potential_connections['end'][j] <= potential_connections['end'][i]:
                    potential_connections['group'][j] = group
                else:
                    pass
        else:
            pass

    # Add group information to AnnData object
    AnnData.var['group'] = potential_connections['group']

    return AnnData
    