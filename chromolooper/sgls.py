import os
import re
import numpy as np
import pybedtools as pbt
import subprocess as sp
import io
import pandas as pd
import json
import itertools as it
import requests

###########################################################
## set global functions
###########################################################

# set global values for chromosome numbers and strings
CHROM_NUMS = list(range(1, 23)) + ['X']
CHROM_STRS = ['chr{}'.format(x) for x in CHROM_NUMS]

# set global variables for softwares using the full path
BGZIP = '/mnt/BioApps/tabix/tabix-0.2.6/bgzip'
TABIX = '/mnt/BioApps/tabix/tabix-0.2.6/tabix'
COOLER_CLI = '/mnt/bioadhoc-temp/Groups/vd-ay/jreyna/software/mamba/envs/hichip-db/bin/cooler'
PAIRLIFTOVER = "/mnt/bioadhoc-temp/Groups/vd-ay/jreyna/software/mamba/envs/pairliftover/bin/pairLiftOver"
COOLER_CLI = '/mnt/bioadhoc-temp/Groups/vd-ay/jreyna/software/mamba/envs/hichip-db/bin/cooler'

# set global variables for softwares using the path up to the directory
# this allows for adding these paths to the PATH varibles using os.ENVIRONS
COOLER_CLI_DIR = '/mnt/bioadhoc-temp/Groups/vd-ay/jreyna/software/mamba/envs/hichip-db/bin/'
BEDTOOLS_DIR = '/mnt/BioApps/bedtools/bin/'
dir_names = [COOLER_CLI_DIR, BEDTOOLS_DIR]
for dy in dir_names:
    os.environ['PATH'] += ':{}'.format(dy)

# set global values for file columns
BED_COLS = ['chr', 'start', 'end']
BEDPE_COLS = ['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB', 'score']
BEDPE_ID_COLS = ['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB', 'score']

# LJI Config
lji_config = '/mnt/bioadhoc-temp/Groups/vd-ay/jreyna/projects/dchallenge/config/softwares/chromolooper.config.txt'

###########################################################
## utility functions
###########################################################

def rename_col_with_ints(df, idxs, renames):
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

def move_col(df, col, loc):
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

def print_enumerate_colnames(df, kind='jupyter'):
    """
    Print the columns of a dataframe with enumerated values.
    
    Params:
    -------
    df: dataframe
    kind: str, either list,  print, or jupyter format
    """
    
    l = list(enumerate(df.columns))
    if kind == 'list':
        print(l)
    elif kind == 'jupyter':
        display(l)
    else:
        for x in l:
            print(x[0], x[1])

def add_prefix_to_names(names, prefix, check=True):
    """
    Given list of names, add a prefix.

    Params:
    -------
    names: list, list of strings
    prefix: str, string that should be used as a prefix
    check: whether to check each name for the current prefix
    """
    
    l = []
    if check == True:
        for x in names:
            # stringify x
            s = str(x)

            # check if prefix is present already
            if not s.startswith(prefix):
                s = prefix + str(s)
            l.append(s)
    else:
        for x in l:
            # stringify x
            s = prefix + str(x)
            l.append(s)
    return(l)

def add_suffix_to_names(names, suffix, check=True):
    """
    Given list of names, add a suffix.
    """

    l = []
    if check == True:
        for x in names:
            # stringify x
            s = str(x)

            # check if prefix is present already
            if not s.endswith(suffix):
                s = str(s) + suffix
            l.append(s)
    else:
        for x in l:
            # stringify x
            s = str(x) + suffix
            l.append(s)
    return(l)
    
def df_to_tsv_file(df, fn, header=True):
    """
    Write a pandas dataframe to a tsv file format.

    Params:
    -------
    df: dataframe
    fn: str, output path
    header: bool
    """
    df.to_csv(fn, sep='\t', header=header, index=False)
    return(True)

def df_to_bed_file(df, fn, chr_col=0, start_col=1, end_col=2, header=False, sort=False, other_cols=[]):
    """
    Write a pandas dataframe to bed file format. Currently only uses int 
    column values.

    Params:
    -------
    df: dataframe
    fn: str, output path
    chr_col: int 
    start_col: int
    end_col: int
    other_cols: list, list of integers
    """

    cols = [chr_col, start_col, end_col] + other_cols
    if sort == True: 
        sort_cols = df.columns[[0,1,2]].tolist()
        df = df.sort_values(sort_cols)
    df.iloc[:, cols].to_csv(fn, sep='\t', header=header, index=False)
    return(True)
        
def print_list(l, sep='\n'): 
    s = [str(x) for x in l]
    s = sep.join(s)
    print(s)

def range_list(start, end):
    r = list(np.arange(start, end))
    return(r)

def optimal_samples_from_groupby(df, grp_col, score_col, func=max):
    data = []
    for grp, grp_df in df.groupby([grp_col]):
        max_val = grp_df[score_col].agg(func)
        idx = grp_df[score_col].tolist().index(max_val)
        point = grp_df.iloc[idx,:]
        data.append(point)
    data = pd.DataFrame(data)
    return(data)

def itertools_list(l, func=it.product):
    p = list(func(*l))

def get_prefix_cols(df, prefix, compliment=False):
    """
    Find columns which use the prefix.

    Params:
    -------
    df: dataframe
    prefix: str
    compliment: bool, if true returns the complimentary columns
    """
    col_matches = df.columns.str.match('^{}'.format(prefix))

    if compliment == False: 
        cols = df.columns[col_matches]
    else:
        cols = df.columns[~col_matches]
    return(cols)

def get_bedpe_prefix_cols(df, prefix):
    """
    Get the columns but set the first 6 columns to BEDPE format columns.
    """
    
    # make a list of the bedpe cols
    bedpe_cols = ['{}{}'.format(prefix, x) for x in BEDPE_COLS]
    
    # get all thee columns with the prefix
    col_matches = df.columns.str.match('^{}'.format(prefix))
    all_cols = set(df.columns[col_matches].tolist())
    
    # extract extra columns with the prefix (non-default BEDPE columns)
    extra_cols = all_cols.difference(bedpe_cols)
    extra_cols = sorted(extra_cols)
    reorder = bedpe_cols + extra_cols
    
    return(reorder)

def check_url(url):
    response = requests.get(url)
    if response.status_code == 200:
        return(1)
    else:
        return(0) 

def make_lji_url(fn, lji_base_url = 'https://informaticsdata.liai.org/'):
    """
    LJI specific URL checker.
    
    Params:
    -------
    fn: str, path to file (either relative or absolute)
    lji_base_url: str, url base path 
    
    Returns:
    url: str, url link to file
    """
    abs_path = os.path.abspath(fn)
    url = abs_path.replace('/mnt/', '')
    url = os.path.join(lji_base_url, url)
    
    if not check_url(url):
        msg = 'Absolute path is not within a web accessible directory' + \
                '(/mnt/bioadhoc-temp/ or /mnt/BioAdHoc/) or incorrect.'
        raise Exception(msg)
    else:
        return(url)


def match_multiple_regexes(s, regexes):
    if type(regexes) != list:
        regexes = [regexes]
        
    for regex in regexes:
        if regex.match(s):
            return(True)
    return(False)      

def extract_cols_with_prefixes(cols, prefixes=None):
    
    if type(prefixes) != list:
        prefixes = [prefixes]
    
    regexes = []
    for prefix in prefixes: 
        regex = re.compile('^{}'.format(prefix))
        regexes.append(regex)

    rcols = []
    for col in cols: 
        if match_multiple_regexes(col, regexes):
            rcols.append(col)
    return(rcols)

def find_bed_like_columns(df):
    """
    Return all columns that have BED interval columns.
    For this package the convention is to use chr, start, end
    """
    
    bed_like = []
    regexes = [re.compile('.*{}.*'.format(x)) for x in BED_COLS[0:3]]
    
    for col in df.columns.tolist():
        if match_multiple_regexes(col, regexes):
            bed_like.append(col)
    return(bed_like)


def find_bedpe_like_columns(df):
    """
    Return all columns which have BEDPE interval columns.
    For this package the convention is to use chrA, startA, endA, chrB, startB, endB
    """
    
    bedpe_like = []
    regexes = [re.compile('.*{}.*'.format(x)) for x in BEDPE_COLS[0:6]]
    
    for col in df.columns.tolist():
        if match_multiple_regexes(col, regexes):
            bedpe_like.append(col)
    return(bedpe_like)




###########################################################
## loop data and metadata functions
###########################################################

def create_loop_id(chrA, startA, endA, chrB, startB, endB, *args):
    r"""
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

def read_multiple_tables_to_df(fns, sep='\t', header=None, *args, **kwargs):
    
    """
    Given a list of files (fns) with the same structure, return a dataframe of everything concatenated.
    An additional file_name column will be added to indicate the data source. 
    
    Params:
    -------
    fns: list, list of file paths
    sep: str, string that should be used for splitting the rows
    """
    
    data = []
    for fn in fns:
        tdf = pd.read_table(fn, sep=sep, header=header, *args, **kwargs)
        tdf['file'] = fn
        data.append(tdf)
    df = pd.concat(data)
    return(df)

###########################################################
## SGL analysis functions
###########################################################

def snp2bed_df(df, chrom=0, pos=1, drop=False):
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

# adding an additional name to snp2bed_df which suggests this
# functionlity is general
pos2bed_df = snp2bed_df

def get_promoter_pos(start, end, strand):
    """
    Get the start position of the promoter.

    Params:
    -------
    start: int, position where the gene starts
    end: int, position where the gene ends
    strand: str, directionality of the gene

    Returns:
    --------
    promoter_pos: int, single base-pair position where the promoter starts
    """
    if strand == '+':
        promoter_pos = start
    elif strand == '-':
        promoter_pos = end
    return(promoter_pos)

def create_promoter_position_col(df, start_col=1, end_col=2, strand_col=3):
    """
    Get the start position of the promoter for a whole column.

    Params:
    -------
    df: dataframe
    start: int, position where the gene starts
    end: int, position where the gene ends
    strand: str, directionality of the gene

    Returns:
    --------
    promoter_positions: list, list of single base-pair position where
        the promoter starts
    """

    mat = df.values
    positions = []
    for i in range(df.shape[0]):
        start = mat[i, start_col]
        end = mat[i, end_col]
        strand = mat[i, strand_col]
        promoter_pos =  get_promoter_pos(start, end, strand)
        positions.append(promoter_pos)
    return(positions) 

def interval_overlap(ref_start, ref_end, query_start, query_end, rtype='bool'):
    """
    interval the overlap between ref and query. Return message specifies 
    the type of overlap.
    
    Params:
    -------
    ref_start: int
    ref_end: int
    query_start: int
    query_end: int
    rtype: string, short for return type and allows for verbose or bool
    
    Returns:
    --------
    If rtype = verbose then returns a descriptive string, otherwise returns True/False
    
    """
    if ref_end < query_start: 
        overlap = 'Outside left'
    elif ref_start > query_end:
        overlap = 'Outside right'
    elif (query_start >= ref_start) & (query_end <= ref_end):
        overlap = 'Overlap surround'
    elif (ref_start < query_start) & (ref_end >= query_start & ref_end <= query_end):
        overlap = 'Overlap left-overhang'
    elif (ref_start >= query_start) & (ref_end <= query_end):
        overlap = 'Overlap within'
    elif (ref_start >= query_start & ref_start <= query_end) & (ref_end > query_end):
        overlap = 'Overlap right-overhang'
        
    if rtype =='verbose':
        return(overlap)
    else:
        if 'Outside' in overlap:
            return(False)
        else:
            return(True)
        
    #     # quick and dirty interval 
    #     query_start, query_end = [5,10]

    #     # no-overlap - on left (False)
    #     ref_start, ref_end = [0,1]
    #     t = interval_overlap(ref_start, ref_end, query_start, query_end)
    #     print(t)

    #     # no-overlap - on right (False)
    #     ref_start, ref_end = [11,12]
    #     t = interval_overlap(ref_start, ref_end, query_start, query_end)
    #     print(t)

    #     # overlap - with left-overhang (True)
    #     ref_start, ref_end = [4, 7]
    #     t = interval_overlap(ref_start, ref_end, query_start, query_end)
    #     print(t)

    #     # overlap - within (True)
    #     ref_start, ref_end = [6, 7]
    #     t = interval_overlap(ref_start, ref_end, query_start, query_end)
    #     print(t)

    #     # overlap - with right-overhang (True)
    #     ref_start, ref_end = [6, 11]
    #     t = interval_overlap(ref_start, ref_end, query_start, query_end)
    #     print(t)

    #     # overlap - surround (True)
    #     ref_start, ref_end = [4, 11]
    #     t = interval_overlap(ref_start, ref_end, query_start, query_end)
    #     print(t)
    
    return(overlap)

def find_anchor_overlap(anchor_chrA, anchor_startA, anchor_endA,
                    anchor_chrB, anchor_startB, anchor_endB,
                    oneD_chr, oneD_start, oneD_end):
    
    """
    Find which anchor a given feature overlaps with A being the left anchor
    and B being the right.
    
    Params:
    -------
    anchor_chrA: str
    anchor_startA: int
    anchor_endA: int
    anchor_chrB: str
    anchor_startB: int
    anchor_endB: int
    oneD_chrB: str
    oneD_startB: int
    oneD_endB: int
    
    Returns:
    --------
    String with either empty (no overlap), 'A' (overlap with A only), 'B' (overlap with B only),
    'AB', overlap with both.
    """
    
    # check the A anchor 
    overlapA = ''
    if (oneD_chr == anchor_chrA) and (oneD_chr == anchor_chrA):
        check_overlap = interval_overlap(anchor_startA, anchor_endA, oneD_start, oneD_end)        
        if check_overlap == True: 
            overlapA = 'A'

    # check the B anchor 
    overlapB = ''
    if (oneD_chr == anchor_chrB) and (oneD_chr == anchor_chrB):
        check_overlap = interval_overlap(anchor_startB, anchor_endB, oneD_start, oneD_end)        
        if check_overlap == True: 
            overlapB = 'B'
        
    overlap_type = overlapA + overlapB
    return(overlap_type)

def create_anchor_overlap_col(df, anchor_chrA, anchor_startA, anchor_endA,
                    anchor_chrB, anchor_startB, anchor_endB,
                    oneD_chr, oneD_start, oneD_end):
    
    """
    Find which anchor a given feature overlaps with A being the left anchor
    and B being the right. This is build upon find_anchor_overlap() and meant
    to be used with dataframes.
    
    Params:
    -------
    df: dataframe
    anchor_chrA: int
    anchor_startA: int
    anchor_endA: int
    anchor_chrB: int
    anchor_startB: int
    anchor_endB: int
    oneD_chrB: int
    oneD_startB: int
    oneD_endB: int
    
    Returns:
    --------
    String with either empty (no overlap), 'A' (overlap with A only), 'B' (overlap with B only),
    'AB', overlap with both.
    """
    
    overlap_types = []
    mat = df.values
    for i in range(df.shape[0]):
        row = mat[i,:]
        overlap_type = find_anchor_overlap(row[anchor_chrA], row[anchor_startA], row[anchor_endA],
                            row[anchor_chrB], row[anchor_startB], row[anchor_endB],
                            row[oneD_chr], row[oneD_start], row[oneD_end])
        overlap_types.append(overlap_type)
    return(overlap_types)

def find_sgls(loops, snps, genes):
    r"""
        Params:
        -------
        loops: dataframe with BEDPE columns chrA, startA, endA, chrB, startB, endB
        snps: dataframe with snp columns chr, pos
        genes: dataframe with BED columns chr, start, end

        Params:
        -------
        sgls: dataframe of SNPs with columns chrA, startA, endA, chrB, startB, endB, snp-start, snp-end, gene-start, gene-end
    """
    pass

def loop_to_1d(loops, snps, genes):
    r"""
        Params:
        -------
        loops: dataframe with BEDPE columns chrA, startA, endA, chrB, startB, endB
        snps: dataframe with snp columns chr, pos
        genes: dataframe with BED columns chr, start, end

        Params:
        -------
        sgls: dataframe of SNPs with columns chrA, startA, endA, chrB, startB, endB, snp-start, snp-end, gene-start, gene-end
    """
    pass

def loop_to_many_1d(loops, snps, genes):
    r"""
        Params:
        -------
        loops: dataframe with BEDPE columns chrA, startA, endA, chrB, startB, endB
        snps: dataframe with snp columns chr, pos
        genes: dataframe with BED columns chr, start, end

        Params:
        -------
        sgls: dataframe of SNPs with columns chrA, startA, endA, chrB, startB, endB, snp-start, snp-end, gene-start, gene-end
    """
    pass


def get_col_integer_indexes(df, cols):
    # case i) all cols are integers, do nothing 
    if all([type(x) in [int, np.int64] for x in cols]):
        col_indexes = cols
        
    # case ii) all cols are listed are strings 
    elif all([type(x) == str for x in cols]):
        df_cols = df.columns.tolist()
        col_indexes = [df_cols.index(x) for x in cols]
    else:
        raise Exception('Please check all column are either strings or integers')
    
    return(col_indexes)


def pairtopair_dataframe(bedpe1, bedpe2, cols1=np.arange(0,7), cols2=np.arange(0,7), merge=True):
    """
    Using bedtools pairtopair, intersect the two dataframes. 
    
    Params:
    -------
    bedpe1: dataframe
    bedpe2: dataframe 
    cols1: lists, list of columns to perform the arithmetic, either ints or strs
    cols2: lists, list of columns to perform the arithmetic, either ints or strs
    merge: bool, whether to merge back data from the original dataframe
    
    Returns:
    --------
    dataframe after performing pairtopair and merge information 
    
    Usage:
    ------    
    # obtain an intersect-p2p dataframe with all previous metadata (merge=True)
    pairtopair_df(df1, df2, cols1, cols2)
    
    # obtain an intersect-p2p dataframe with just loop ids
    pairtopair_df(df1, df2, cols1, cols2, merge=False)
    """
    
    # check to make sure the pybedtools oly have 7 columns each
    if len(cols1) != 7 or len(cols2) != 7:
        raise Exception('cols1 and cols2 should only have 7 columns.')
        
    # get column indexes and names
    col_indexes1 = get_col_integer_indexes(bedpe1, cols1)
    col_indexes2 = get_col_integer_indexes(bedpe2, cols2)
    col_names1 = bedpe1.columns[col_indexes1].tolist()
    col_names2 = bedpe2.columns[col_indexes2].tolist()

    # create the BedTools objects   
    pbt1 = pbt.BedTool.from_dataframe(bedpe1.iloc[:, col_indexes1])
    pbt2 = pbt.BedTool.from_dataframe(bedpe2.iloc[:, col_indexes2])
        
    # perform the pair to pair 
    p2p = pbt1.pairtopair(pbt2)    
    p2p = p2p.to_dataframe(disable_auto_names=False, header=None).iloc[:, 0:14]
    p2p.columns = col_names1 + col_names2
    
    # merge original dataframes with pairtopair
    if merge == True: 
        merge_df = bedpe1.merge(p2p, how='left').merge(bedpe2, how='left')
        return(merge_df)
    else:
        return(p2p)



def get_grp_max(df, max_col, grp_cols=[]):
    res = df.sort_values(grp_cols, ascending=False).drop_duplicates(subset=max_col)
    return(res)

def get_grp_min(df, min_col, grp_cols=[]):
    res = df.sort_values(grp_cols, ascending=False).drop_duplicates(subset=min_col)
    return(res)






















###########################################################
## WashU Browser functions 
###########################################################

def create_washu_second_anchor_col(df, chrB_col=4, startB_col=5, endB_col=6, score_col=7):
    """ 
    Pass a dataframe with loop data and create the second anchor column. You can specify the index or header
    names using ALL ints or ALL strs using the *_col params. An Errors will be thrown if done
    incorrectly. 
    Params:
    -------
    df: dataframe
    chrB_col: str/int
    startB_col: str/int
    endB_col: str/int
    score_col: str/int
    """
    
    col_check = [chrB_col, startB_col, endB_col, score_col]
    int_check = all([type(x) == int for x in col_check])
    str_check = all([type(x) == str for x in col_check])
    
    if int_check == False and str_check == False: 
        raise Exception('*_col params are not all int or all str')

    elif int_check == True:
        col =  df.iloc[:, chrB_col] + ':' + \
                df.iloc[:, startB_col].astype(str) + '-' + \
                df.iloc[:, endB_col].astype(str) + ',' +  \
                df.iloc[:, score_col].astype(str)
    else:
        col =  df.loc[:, chrB_col] + ':' + \
                df.loc[:, startB_col].astype(str) + '-' + \
                df.loc[:, endB_col].astype(str) + ',' +  \
                df.loc[:, score_col].astype(str)
    
    return(col)


def make_washu_1d_dict(track_type, name, url, show_on_hub_load=True,
                        ensemble_style = True, height = 200,
                        color='red', display_mode='arc'):
    """
    Create dicctionary for 1D Washu track files. 
    
    Params:
    -------
    track_type: str
    name: str
    url: str
    show_on_hub_load: bool
    ensemble_style: bool
    height: int
    color: str
    display_mode: str
    """
    d = dict()
    d['type'] = track_type
    d['name'] = name
    d['url'] = url
    d['showOnHubLoad'] = show_on_hub_load
    d['options'] = {'ensembleStyle': ensemble_style,
                    'height': height,
                    'color': color,
                    'displayMode': display_mode}
    return(d)

def make_washu_longrange_dict(name, url, show_on_hub_load=True,
                              ensembleStyle = True, height = 200,
                              color='red', displayMode='arc'):    
    """
    Create dicctionary for longrang Washu track files. 
    
    Params:
    -------
    track_type: str
    name: str
    url: str
    show_on_hub_load: bool
    ensemble_style: bool
    height: int
    color: str
    display_mode: str
    """
    
    d = dict()
    d['type'] = 'longrange'
    d['name'] = name
    d['url'] = url
    d['showOnHubLoad'] = show_on_hub_load
    d['options'] = {'ensembleStyle': ensembleStyle,
                    'height': height,
                    'color': color,
                    'displayMode': displayMode}
    
    return(d)


def create_washu_hub(hub_entries, out_fn):
    """
    Create a Washu Hub from a list of tracks
    
    Params:
    -------
    hub_entries: list, list of dictionaries where each dictionary contain information
        necessary for a Washu Track
    out_fn: str, output path
    """
    
    with open(out_fn, 'w') as fw:
        s = json.dumps(hub_entries,  indent=True)
        fw.write(s) 

###########################################################
## BGZIP functions 
###########################################################

def tabix(fn, force=True, verbose=False):
    """
    Index a file that has been compressed with bgzip. Will return the path to the indexed file if 
    sucessfully completed, otherwise the error message will be printed out
    for debugging. 
    
    Params:
    -------
    fn: str, path to the bed file
    force: bool, whether to force bgzip to run or not 
    verbose: bool, will print out the subprocess message
    """
    
    # set the force parameter
    force_param = ''
    if force == True:
        force_param = '-f'
    
    # make an index file
    idx_fn = fn + '.tbi'
    cmd = '{} {} {}'.format(TABIX, force_param, fn)
    msg = sp.check_output(cmd, shell=True)
    if verbose == True:
        print(msg.decode())
    return(idx_fn)


def bgzip(fn, keep_bed=True, index=True, force=True, verbose=False):
    """
    Make a gz file with bgzip. Will return the path to the bed.gz file if 
    sucessfully completed, otherwise the error message will be printed out
    for debugging. 
    
    Params:
    -------
    fn: str, path to the bed file
    keep_bed: bool, keeps the original bed file
    index: bool, index the file with tabix as well as compress
    force: bool, whether to force bgzip to run or not 
    verbose: bool, will print out the subprocess message
    
    Returns:
    --------
    list of files with the first one being the compressed file and the
    second being the index file if index=True, otherwise None.
    """
    
    # set the force parameter
    force_param = ''
    if force == True:
        force_param = '-f'
    
    # compress the file  
    cmp_fn = fn + '.gz'
    if keep_bed == True:
        cmd = '{} -c {} > {}'.format(BGZIP, fn, cmp_fn)
    else:
        cmd = '{} {} {}'.format(BGZIP, force_param, fn)
    msg = sp.check_output(cmd, shell=True)
    
    if verbose == True:
        print(msg.decode())
        
    # index the file if index=True
    idx_fn = None
    if index == True:
        idx_fn = tabix(cmp_fn)
        
    return([cmp_fn, idx_fn])

###########################################################
## Cooler functions 
###########################################################

def cooler_py_fetch_chroms(c, chroms='*', drop_bal=True):
    """
    Fetch data from whole chromosomes from a cooler object.
    
    Params:
    -------
    c: cooler
    chroms: str/list, if a string it must either be * which represents all chromosomes,
        a single chromosome str (e.g. chr1). If a list then elements must be chromosome
        strs (e.g. ['chr1', 'chr2'])
    drop_bal: bool, drop the balance column if not necessary
        
    Returns:
    --------
    df: dataframe with interactions across all chromosomes
    """
    
    matrix_selector = c.matrix(join=True, as_pixels=True)

    # get the correct chromosomes
    select_chroms = CHROM_STRS
    if chroms == '*':
        select_chroms = CHROM_STRS
        
    # get the data
    data = []
    for chrom in select_chroms:
        data.append(matrix_selector.fetch(chrom))
    
    df = pd.concat(data)
    df.drop('balanced', axis=1, inplace=True)
    
    return(df)

def cooler_cli_mcool2bedpe(mcool, res, to_df=True, outfn=None):
    """
    Dump the cooler file using the command line cooler tool. Much faster
    than the Python API if you are trying to extract interactions across
    several chromosomes.
    
    Params:
    -------
    mcool: str, path to the mcool file 
    res: int, resolution of the dataset
    to_df: bool, whether to return a dataframe
    out_fn: str, file path for output if not dumped to dataframe

    Returns:
    --------
    dataframe_or_true: returns a dataframe when to_df=True or True when the mcool is dumped 
        to a file
    """
    
     # dump the mcool file to a dataframe
    if to_df == True:
        mcool_to_dump = '{} dump --join {}::/resolutions/{}'.format(COOLER_CLI, mcool, res)        
        dump = sp.check_output(mcool_to_dump, shell=True).decode()
        df = pd.read_table(io.StringIO(dump), header=None, names=BEDPE_COLS)
        return(df)
    
    # dump the mcool file to a file
    else:
        if outfn != None: 
            mcool_to_dump = '{} dump --join {}::/resolutions/{} > {}'.format(COOLER_CLI, mcool, res, outfn)  
            return(True)
        else:
            raise Exception('Must specify outfn when not dumping to dataframe.')

# providing an additional name for cooler_cli_mcool_bedpe
cooler_cli_mcool_dump = cooler_cli_mcool2bedpe

###########################################################
## Utility functions 
###########################################################

def update_global_config_vars(config=lji_config):
    """
    Update the global config variables using a file for dictionary.
    
    Params:
    -------
    config: dict/str, if dictionary then keys are global variables (check print_config_vars() to 
        view the current values) and if str then the config file is read and parsed accordingly. 
    """
    with open(config) as f:
        for line in f:
            if '=' in line:

                # update the SGL global vars
                var_str = line.strip()
                sgl_update_cmd = 'sgls.{}'.format(var_str)
                print(sgl_update_cmd)
                exec(var_str)

                # if value is a directory then also add to PATH variable
                if 'DIR' in line:
                    dir_name = line.split('=')[1].strip().replace("'", "")
                    os.environ['PATH'] += ':{}'.format(dir_name)
                    
def print_config_vars():
    """
    Print the global configuration variables.
    """
    global_vars = { 
            'CHROM_NUMS': CHROM_NUMS,
            'CHROM_STRS': CHROM_STRS,
            'BGZIP': BGZIP,
            'TABIX': TABIX,
            'COOLER_CLI': COOLER_CLI,
            'PAIRLIFTOVER': PAIRLIFTOVER,
            'COOLER_CLI_DIR': COOLER_CLI_DIR,
            'BEDTOOLS_DIR': BEDTOOLS_DIR,
            'BEDPE_COLS': BEDPE_COLS}
    for k, v in global_vars.items():
        print('{}: {}'.format(k,v))
