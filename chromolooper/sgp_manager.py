from collections import SortedDict
from .sgls import *

class sgpair_manager(): 

    def __init__(self, snp_coords=None, gene_coords=None, loop_coords=None):
        '''
            Constructor for an sgpair manager.

            Params:
                snp_coords : str or dataframe
                    If a str, then this must be a path to a BED file with 4 columns:
                    chr, start, end, id
                gene_coords : str or dataframe
                    If a str, then this must be a path to a BED file with 4 columns:
                    chr, start, end, id
                loop_coords : str or dataframe
                    If a str, then this must be a path to a BEDPE file with 7 columns:
                    chrA, startA, endA, chrB, startB, endB, id. 
        '''
        self.snp_coords = snp_coords
        self.gene_coords = gene_coords
        self.loop_coords = loop_coords
        self.sgpair_coords = None
        self.snp_intersects = SortedDict()
        self.gene_intersects = SortedDict()
        self.loop_intersect = SortedDict()
        self.sgpair_intersect = SortedDict()
        pass

    def add_snp_coords(self, snp_coords=None):
        '''
            Add the main SNP coordinates data.
        '''
        pass

    def add_gene_coords(self, gene_coords=None):
        '''
            Add the main gene coordinates data.
        '''
        pass

    def add_loop_coords(self, loop_coords=None):
        '''
            Add the main loop coordinates data.
        '''
        pass

    def create_sgls(self):
        '''
            Intersect the SNP, gene and loop data for form SNP-Gene pairs with loops.
        '''
        pass

    def intersect_with_snps(self, other_coords=None):
        '''
            Intersect SNP coordinates with 1D or 2D dataset.
        '''
        pass

    def intersect_with_genes(self, other_coords=None):
        '''
            Intersect gene coordinates with 1D or 2D dataset.
        '''
        pass

    def intersect_with_loops(self, other_coords=None):
        '''
            Intersect loop coordinates with 1D or 2D dataset.
        '''
        pass

    def intersect_with_sgpairs(self, other_coords=None):
        '''
            Intersect SG coordinates with 1D or 2D dataset.
        '''
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

    def make_washu_coords_hub(self, fn):
        '''
            Make the corresponding WashU tracks and hub.
        '''

        # add SNP coordinates and SNP intersects first
        snp_coord_dict = make_washu_1d_dict(self.snp_coords)
        snp_intersects_dicts = []
        for snp_intrs in self.snp_intersects:
            d = make_washu_1d_dict(snp_intrs)
            snp_intersects_dicts.append(d)

        # add gene coordinates and gene intersects second
        gene_coord_dict = make_washu_1d_dict(self.gene_coords)
        gene_intersects_dicts = []
        for gene_intrs in self.gene_intersects:
            d = make_washu_1d_dict(gene_intrs)
            gene_intersects_dicts.append(d)

        # add loop coordinates and loop intersects third
        loop_coord_dict = make_washu_longrange_dict(self.loop_coords)
        loop_intersects_dicts = []
        for loop_intrs in self.loop_intersects:
            d = make_washu_longrange_dict(loop_intrs)
            loop_intersects_dicts.append(d)

        # add sgpair coordinates and sgpair intersects last
        sgpair_coord_dict = make_washu_longrange_dict(self.sgpair_coords)
        sgpair_intersects_dicts = []
        for sgpair_intrs in self.sgpair_intersects:
            d = make_washu_longrange_dict(sgpair_intrs)
            sgpair_intersects_dicts.append(d)

        # collect the hubs and save to datahub file
        hub_entries = [snp_coord_dict] + snp_intersects_dicts + \
                        [gene_coord_dict] + gene_intersects_dicts + \
                        [loop_coord_dict] + loop_intersects_dicts + \
                        [sgpair_coord_dict] + sgpair_intersects_dicts
        create_washu_hub(hub_entries, fn)