# ------------------------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------------------------

def get_runs(sample,experiment):
    """
    Returns the run accessions for a given sample (biosample).
    """
    df = pep.subsample_table[pep.subsample_table.sample_name == sample]
    return list(set(df[df.Experiment == experiment].Run))

def get_experiments(sample):
    """
    Returns the experiment accessions for a given sample (biosample).
    """
    return list(set(pep.subsample_table[pep.subsample_table.sample_name == sample].Experiment))

def is_paired_end(sample):
    """
    Returns the libary layout for a given sample. Assumes that all libraries for a given sample
    use the same layout, otherwise returns a string.
    """
    layouts_represented = set(pep.subsample_table[pep.subsample_table.sample_name == sample].LibraryLayout)
    if len(layouts_represented) > 1:
        return "ERROR! make sure each sample has only paired end or single end libraries."
    return all([x == "PAIRED" for x in layouts_represented])

def get_background_sample(sample):
    """
    Returns the 'input' or WCE sample for IP type experiments or others
    where a background sample is required for analysis. Operates on the level of samples, not
    subsamples. Fails/returns NA if the column 'input' doesn't exist in the sample_table.
    """
    x=pep.sample_table[pep.sample_table.sample_name == sample].input[0]
    return x

flatten = lambda t: [item for sublist in t for item in sublist]

