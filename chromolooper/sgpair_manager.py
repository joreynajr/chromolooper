from .sgls import *
import pandas as pd
import numpy as np
import pybedtools as pbt

# Should try using metadata/group for merging
# So you create a group tag for each resource 
# you can include disease, organ, tissue, celltype
# etc and at the end we will be able to specify which 
# meta data we should merge on. 
# OR 
# Each spg_manager instance is used with a different
# context/metadata


#############################################
######### Module functions ##################
#############################################

def get_variant_hgvs_id(): # Unfinished
    '''
        Get the variant hgvs ID.
    '''
    pass

def get_variant_rsid(): # Unfinished
    '''
        Get the variant hgvs ID.
    '''
    pass

def create_basic_snp_id(chrom, bp): # DONE
    '''
        Create a basic SNP ID that is composed of <chr>:<pos>. Sometimes called a sid.
    '''
    return('{}:{}'.format(chrom, bp))

def unversion_geneid(): # Unfinished
    '''
        Remove the version number from an ID.
    '''
    pass

def create_signal_name(): # Unfinished
    '''
        Analyze several SNPs and create a signal region that is easily
        searchable with WashU.
    '''
    pass

def compliment_values(full, subset): # Unfinished
    '''
        Find the complimentary values between full and subset. 
    '''
    full_set = set(full)
    subset_set = set(subset)
    compliment = sorted(full_set.difference(subset_set))
    return(compliment)

def rename_col_with_ints(df, idxs, renames): # DONE 
    """
    Rename the columns of a pandas dataframe using a list of integer indexes.
    
    Params:
    -------
    df: dataframe
    idxs: list of integer idxs 
    renames: list of strings to replace columns
    """
    
    d = {df.columns[i]: renames[i] for i in idxs}
    new_df = df.rename(columns=d)
    return(new_df)

def move_col(df, col, loc): # DONE
    """
    Move a column to a new index loc.
    
    Params:
    -------
    df: dataframe
    col: str/int, name or integer index of the column you want to move
    loc: int, integer index of where you want to move the column
    
    Returns:
    --------
    Dataframe is modified in place and returned (technically the same memory location).  
    """
    
    colname = col
    if type(col) == int:
        colname = df.columns[col]
    
    values = df.pop(colname)
    df.insert(loc, colname, values)
    
    return(df)

def snp_bp_to_bed(df, chrom=0, pos=1, drop=False): # DONE
    """
    Convert from a SNP dataframe to a BED dataframe for bed intersections.
    
    Params:
    -------
    df: dataframe, with chromosome and pos
    chrom: str/int, with name/index of the chromosome column
    pos: str/int, with name/index of the position column
    drop: bool, whether to drop or not the index
    
    Returns:
    --------
    
    """
    
    chrom_idx = chrom
    if type(chrom) == str: 
        chrom_idx = df.columns.get_loc(chrom)
        
    pos_idx = pos
    if type(pos) == str: 
        pos_idx = df.columns.get_loc(pos)

    # move the chrom and position if necessary  
    tdf = move_col(df, chrom_idx, 0) # move to the first col
    tdf = move_col(tdf, pos_idx, 1) # move to the second col
        
    # add start 
    start = tdf.iloc[:, 1] - 1
    tdf.insert(1, 'start', start)
    
    # rename the columns 
    tdf = rename_col_with_ints(df, [0,1,2], BED_COLS)
    
    return(tdf)


def create_loop_id(chrA, startA, endA, chrB, startB, endB, *args):
    """
        Params:
        -------
        chrA: str, chromosome of the first interaction
        startA: str, start of the first interaction
        endA: str, end of the first interaction
        chrB: str, chromosome of the second interaction
        startB: str, start of the second interaction
        endB: str, end of the second interaction
        extras: list, additional information that should can be added
    """
    
    l = [chrA, startA, endA, chrB, startB, endB] + list(args)
    l = [str(x) for x in l]
    s = '.'.join(l)
    return(s)

def create_loop_id_col(df, chrA_col=0, startA_col=1, endA_col=2,
                            chrB_col=3, startB_col=4, endB_col=5, extras=[]):
    """
    Pass a dataframe with loop data and create a loop ID columns. You can specify the index or header
    names using ALL ints or ALL strs using the *_col praams. An Errors will be thrown if done
    incorrectly. 
    Params:
    -------
    df: dataframe
    chrA_col: str/int
    startA_col: str/int
    endA_col: str/int
    chrB_col: str/int
    startB_col: str/int
    endB_col: str/int
    """
    
    l = [chrA_col, startA_col, endA_col, chrB_col, startB_col, endB_col] + list(extras)
    IDs = []
    mat = df.values
    for i in range(df.shape[0]):
        row = mat[i, l]
        ID = create_loop_id(*row)
        IDs.append(ID)

    return(IDs)











#############################################
######### sgpair_manager class ##############
#############################################
class sgpair_manager(): 

    def __init__(self, snp_data=None, gene_data=None, loop_data=None,
                    snp_coord_cols=np.arange(0,2), snp_id_col=2,
                    gene_coord_cols=np.arange(0,3), gene_id_col=3,
                    loop_coord_cols=np.arange(0,6), loop_id_col=7,
                    snp_format='bp', ref_genome='hg19'):
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
                snp_format: str
                    Choices are bp or bed which refer to using a single position for bed-like
                    positions for the SNP.
        '''

        # coverting arrays to lists
        snp_coord_cols = list(snp_coord_cols)
        gene_coord_cols = list(gene_coord_cols)
        loop_coord_cols = list(loop_coord_cols)

        # raising exceptions for incorrect snp_format specification
        if snp_format not in ['bp', 'bed']:
            raise Exception('snp_format may only be bp or bed.')

        # parse snp data to initialize snp values
        if snp_data is not None:

            # allowing -1 for convenience
            # from time to time, the snp id may be add after loading the input data and
            # using a basic id with create_basic_snp_id
            if snp_id_col == -1:
                snp_id_col = snp_data.shape[1] - 1

            # convert coordinate data to bed-like format as necessary
            self.snp_coords = snp_data.iloc[:, snp_coord_cols + [snp_id_col]]
            if snp_format == 'bp':
                self.snp_coords = snp_bp_to_bed(self.snp_coords)

            # convert the coordinate data into it's final sgp form using
            # pybedtools
            self.snp_coords = pbt.BedTool.from_dataframe(self.snp_coords)

            # rearrange the meta data with snp id at the beginning
            meta_cols = [snp_id_col] + compliment_values(np.arange(snp_data.shape[1]), snp_coord_cols + [snp_id_col])
            self.snp_meta = snp_data.iloc[:, meta_cols]


        # parse gene data to initialize gene values
        if gene_data is not None:

            # allowing -1 for convenience, same as snp data
            if gene_id_col == -1:
                gene_id_col = gene_data.shape[1] - 1

            # convert the coordinate data into it's final sgp form using
            # pybedtools
            self.gene_coords = gene_data.iloc[:, gene_coord_cols  + [gene_id_col]]
            self.gene_coords = pbt.BedTool.from_dataframe(self.gene_coords)

            # rearrange the meta data with snp id at the beginning
            meta_cols = compliment_values(np.arange(gene_data.shape[1]), gene_coord_cols)
            meta_cols = [gene_id_col] + compliment_values(np.arange(gene_data.shape[1]), gene_coord_cols + [gene_id_col])
            self.gene_meta = gene_data.iloc[:, meta_cols]

        # parse loop data to initialize loop values
        if loop_data is not None: 

            # allowing -1 for convenience, same as snp data
            if loop_id_col == -1:
                loop_id_col = loop_data.shape[1] - 1

            # convert the coordinate data into it's final sgp form using
            # pybedtools
            self.loop_coords = loop_data.iloc[:, loop_coord_cols + [loop_id_col]]
            self.loop_coords = pbt.BedTool.from_dataframe(self.loop_coords)

            # rearrange the meta data with snp id at the beginning
            meta_cols = compliment_values(np.arange(loop_data.shape[1]), loop_coord_cols)
            meta_cols = [loop_id_col] + compliment_values(np.arange(loop_data.shape[1]), loop_coord_cols + [loop_id_col])
            self.loop_meta = loop_data.iloc[:, meta_cols]

        # initialize intersection data
        self.sgpair_coords = None
        self.snp_intersects = {}
        self.gene_intersects = {}
        self.loop_intersect = {}
        self.sgpair_intersect = {}

        # initialize coloring schemes
        self.snp_color = 'blue'
        self.gene_color = 'blue'
        self.loop_color = 'blue'
        self.snp_intersect_color = 'red'
        self.gene_intersect_color = 'red'
        self.loop_intersect_color = 'red'

        # initialize meta data
        self.ref_genome = 'hg19'

    def add_snp_coords(self, snp_data, coord_cols=np.arange(0,2), format='bp'): # DONE
        '''
            Add the main SNP coordinates data.
        '''
        self.snp_coords = snp_data.iloc[:, coord_cols]
        meta_cols = compliment_values(np.arange(snp_data.shape[1]), coord_cols)
        self.snp_meta = snp_data.iloc[:, meta_cols]

    def add_gene_coords(self, gene_data, coord_cols=np.arange(0,3)): # DONE
        '''
            Add the main gene coordinates data.
        '''
        self.gene_coords = gene_data.iloc[:, coord_cols]
        meta_cols = compliment_values(np.arange(gene_data.shape[1]), coord_cols)
        self.gene_meta = gene_data.iloc[:, meta_cols]

    def add_loop_coords(self, loop_data, coord_cols=np.arange(0,7)): # DONE
        '''
            Add the main loop coordinates data.
        '''
        self.loop_coords = loop_data.iloc[:, coord_cols]
        meta_cols = compliment_values(np.arange(loop_data.shape[1]), coord_cols)
        self.loop_meta = loop_data.iloc[:, meta_cols]

    def create_sgls(self):
        '''
            Intersect the SNP, gene and loop data for form SNP-Gene pairs with loops.
        '''
        pass
















    def intersect_with_snps(self, other_coords=None):
        '''
            Intersect SNP coordinates with 1D or 2D dataset.

            If snp with snp where rsid's or hgvs id's are available, then
            the merge can be done using these values instead.
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

    def liftover(self):
        '''
            Liftover all coordinates to a new reference genome.
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

    def aggregate_with_snps(self, func, score_col=4):
        '''
            Aggregate the data on SNP id's using func and the score_col.
        '''
        pass

    def aggregate_with_genes(self, func, score_col=4):
        '''
            Aggregate the data on gene id's using func and the score_col.
        '''
        pass

    def aggregate_with_loops(self, func, score_col=8):
        '''
            Aggregate the data on loop id's using func and the score_col.
        '''
        pass

    def aggregate_with_sgpairs(self, func):
        '''
            Aggregate the data on sgpair id's using func.
        '''
        pass