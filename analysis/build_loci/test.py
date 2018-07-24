def test():
    """
    
    """
    import os
    from metaproteomics.analysis.DBInfo import DBInfo    
    from metaproteomics import utils
    from metaproteomics.analysis import build_loci    
    from metaproteomics.analysis.build_loci import Sample
    from itertools import chain
    
    DIR_NAME = "/home/gstupp/metaproteomics/analysis/build_loci"
    metadata = build_loci.read_metadata(os.path.join(DIR_NAME,"sample_metadata.csv"))        
    db_info = DBInfo("compil_mgm")
    
    samples = dict()
    for sample_name,sample_metadata in metadata.iteritems():
        samples[sample_name] = Sample(sample_name, sample_metadata.path, db_info, sample_metadata)
        
    protein_clusters = {name:build_loci.build_protein_clusters(sample) for name,sample in samples.items()}

    grouped_loci = build_loci.group_across_samples(list(chain(*protein_clusters.values())))
    
    

if __name__ == "__main__":
    test()
