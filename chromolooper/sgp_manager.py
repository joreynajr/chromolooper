class sgpair_manager(): 

    def __init__(self, snp_data=None, gene_data=None, loop_data=None):
        '''
            Constructor for an sgpair manager.

            Params:
                snp_data : str or dataframe
                    If a str, then this must be a path to a BED file with 4 columns:
                    chr, start, end, id
                gene_data : str or dataframe
                    If a str, then this must be a path to a BED file with 4 columns:
                    chr, start, end, id
                loop_data : str or dataframe
                    If a str, then this must be a path to a BEDPE file with 7 columns:
                    chrA, startA, endA, chrB, startB, endB, id. 
        '''
        self.snp_data = snp_data
        self.gene_data = gene_data
        self.loop_data = loop_data
        pass

    def sgpairs_from_closest_gene(self):
        '''
            Using SNP coordinates as input, looks for the closests gene +/- d distance
            and creates SG-Basis.
        '''
        pass
    def sgpairs_from_closest_snp(self):
        '''
            Using gene coordinates as input, looks for the closests SNP +/- d distance
            and creates a SG-Basis. 
        '''
        pass

    def sgpairs_from_loops(self):
        '''
            Intersecting loops with SNP and gene coordinates to create the SG-Basis.
        '''
        pass