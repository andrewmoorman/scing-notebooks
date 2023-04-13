import re, subprocess, boto3, json, shlex, mysql, os, urllib, logging
import numpy as np
import pandas as pd
from pathlib import Path
import importlib
from mysql.connector import connect, Error

########## SCRIdb queries ##########

def execute_query(query, creds):
    
    user = creds['user']
    password = creds['password']

    with connect(
        host="peer-lab-db.cggxmlwgzzpw.us-east-1.rds.amazonaws.com",
        database="peer_lab_db",
        user=user,
        password=password,
    ) as connection:
        with connection.cursor(buffered=True) as cursor:
            cursor.execute(query)
            result = cursor.fetchall()
    return result

def sample_scridb_info(querys, query_col, creds):
    
    samples = get_sample(querys, query_col, creds)
    
    species = []
    sc_tech = []
    proj_id = []
    reference = []
    
    for query in querys:
        species += get_species(query, query_col, creds)
        sc_tech += get_sc_tech(query, query_col, creds)
        proj_id += get_project_id(query, query_col, creds)
        reference += get_reference(query, query_col, creds)
    
    samples['species'] = species
    samples['sc_tech'] = sc_tech
    samples['project_id'] = proj_id
    samples['reference'] = reference
    
    return samples


def get_sample(querys, query_col, creds):
    user = creds['user']
    password = creds['password']
    
    try:
        table_sample_data = "peer_lab_db.sample_data"
        query = f"""
            SELECT Sample, AWS_storage, id
            FROM {table_sample_data} 
            """
        
        if len(querys) != 1:
            query += f"WHERE {table_sample_data}.{query_col} IN {tuple(querys)}"
        else:
            query += f'WHERE {table_sample_data}.{query_col} = "{querys[0]}"'

        samples = execute_query(query, creds)
        samples = pd.DataFrame(samples)
        samples.columns = ['Sample', 'AWS_storage', 'id']
        samples = samples.set_index('Sample')
        
        samples['AWS_storage'] = samples['AWS_storage'].str.strip('/')
        return samples
    
    except Error as e:
        print(f"Error: {e}")

        
def get_species(query, query_col, creds):
    
    user = creds['user']
    password = creds['password']
    
    try:
        table_sample_data = "peer_lab_db.sample_data"
        table_species = "peer_lab_db.species"
        table_genome_idx = "peer_lab_db.genome_index"
        query = f"""
        SELECT {table_species}.Species
        FROM {table_species}
        LEFT JOIN {table_genome_idx}
        ON {table_species}.id = {table_genome_idx}.species_id
        LEFT JOIN {table_sample_data}
        ON {table_genome_idx}.id = {table_sample_data}.genomeIndex_id
        WHERE {table_sample_data}.{query_col} = "{query}"
        """
        
        species = []
        results = execute_query(query, creds)
        for result in results:
            species.append(result[0].lower())
        return species
    except Error as e:
        print(f"Error: {e}")
                       

def get_sc_tech(query, query_col, creds):
    
    user = creds['user']
    password = creds['password']

    try:
        table_sample_data = "peer_lab_db.sample_data"
        table_sc_tech = "peer_lab_db.sc_tech"
        table_genome_idx = "peer_lab_db.genome_index"
        query = f"""
        SELECT {table_sc_tech}.sc_Tech
        FROM {table_sc_tech}
        LEFT JOIN {table_genome_idx}
        ON {table_sc_tech}.id = {table_genome_idx}.scTech_id
        LEFT JOIN {table_sample_data}
        ON {table_genome_idx}.id = {table_sample_data}.genomeIndex_id
        WHERE {table_sample_data}.{query_col} = "{query}"
        """
        sc_tech = []
        results = execute_query(query, creds)
        for result in results:
            sc_tech.append(result[0])
        return sc_tech
    
    except Error as e:
        print(f"Error: {e}")
        
def get_project_id(query, query_col, creds):
    try:
        table_sample_data = "peer_lab_db.sample_data"
        table_project_data = "peer_lab_db.project_data"
        query = f"""
        SELECT {table_project_data}.projectName
        FROM {table_project_data}
        LEFT JOIN {table_sample_data}
        ON {table_project_data}.id = {table_sample_data}.projectData_id
        WHERE {table_sample_data}.{query_col} = "{query}"
        """
        
        proj_id = []
        results = execute_query(query, creds)
        for result in results:
            proj_id.append(result[0])
        return proj_id
    
    except Error as e:
        print(f"Error: {e}")


def get_reference(query, query_col, creds):
    
    user = creds['user']
    password = creds['password']

    try:
        table_sample_data = "peer_lab_db.sample_data"
        table_genome_idx = "peer_lab_db.genome_index"
        query = f"""
        SELECT {table_genome_idx}.gIndex
        FROM {table_genome_idx}
        LEFT JOIN {table_sample_data}
        ON {table_genome_idx}.id = {table_sample_data}.genomeIndex_id
        WHERE {table_sample_data}.{query_col} = "{query}"
        """
        
        reference = []
        results = execute_query(query, creds)
        for result in results:
            reference.append(result[0])
        return reference
    
    except Error as e:
        print(f"Error: {e}")
        
        
########## FASTQ map ##########

fastq_map = {
    'CellRangerVdj': ['I1','R1','R2'],
    'Hashtag': ['R1','R2'],
    'CiteSeq': ['R1','R2'],
    'AsapSeq': ['R1','R2','R3'],
    'CellRangerATAC': ['I1','R1','R2','R3'],
    'CellRangerArc': ['I1','R1','R2','R3'],
    'CellRangerGex': ['I1','R1','R2'],
    'MitoTracing': ['R1', 'R2'],    
}


########## AWS S3 functions ##########

def get_s3_objects(bucket, key, pattern, full_uri=False):
    s3r = boto3.resource("s3")
    bucket_s3 = s3r.Bucket(bucket)
    objects = []
    for obj in bucket_s3.objects.filter(Prefix=key):
        hit = pattern.search(obj.key)
        if hit:
            objects.append(obj.key)
    if full_uri:
        objects = [f"s3://{bucket}/{o}" for o in objects]
    return objects

# Get fastq file paths on S3 for each file id
# Returns dictionary from id to s3 path
# Throws exception if FASTQs don't exist for any id
def get_fastqs(
    path: str, # path to directory containing FASTQ files
    fastq_file_ids: list, # FASTQ file ids needed for this run type (e.g. I1, R1, R2, etc.)
    folder: str = "",
):
    fastq_map = dict()
    _, bucket, key, _, _ = urllib.parse.urlsplit(f"{path}/{folder}")
    for fid in fastq_file_ids:
        files = get_s3_objects(
            bucket, key.lstrip("/"),
            re.compile(f"_{fid}_\d{{3}}.fastq.gz$")
        )
        try:
            assert files, f"AssertionError: Missing `{fid}` archives!"
            fastq_map[fid] = [os.path.join("s3://", bucket, str(f)) for f in files]
        except AssertionError as err:
            logging.warning("%s\n\t %s", err, path)
            return
    return fastq_map

# Get every FASTQ in a folder
def get_all_fastqs(
    path: str, # path to directory containing FASTQ files
    folder: str = "",
):
    _, bucket, key, _, _ = urllib.parse.urlsplit(f"{path}/{folder}")
    files = get_s3_objects(
        bucket, key.lstrip("/"),
        re.compile(f".fastq.gz$")
    )
        
    try:
        fastqs = [os.path.join("s3://", bucket, str(f)) for f in files]
    except AssertionError as err:
        logging.warning("%s\n\t %s", err, path)
        return
    return fastqs


########## Reference map ##########
reference_map = {}

reference_map['CellRangerArc'] = {
    "human": "https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A.tar.gz",
    "mouse": "https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz"
}

reference_map['CellRangerAtac'] = {
    "human":"https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz",
    "mouse":"https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz"
}

reference_map['CellRangerGex'] = {
    "human":"https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz",
    "mouse":"https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz"
}

reference_map['CellRangerCellPlex'] = {
    "human":"https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz",
    "mouse":"https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz"
}

reference_map['CellRangerVdj'] = {
    "human":"GRCh38",
    "mouse":"GRCm38"
}

def update_ref(samples, prefix):
    for sample, row in samples.iterrows():
        species = row['species']
        
        if not row['reference']:
            if species in ['human', 'mouse']:
                samples.loc[sample, 'reference'] = reference_map[prefix][row['species']]

        elif prefix.startswith('CellRanger'):
            if not row['reference'].startswith('https'):
                
                if not species in ['human', 'mouse']:
                    print(f'{species} reference not in database. Manually change "reference" field')
                    samples.loc[sample, 'reference'] = np.nan
                    
                else:
                    samples.loc[sample, 'reference'] = reference_map[prefix][row['species']]
                                        
    return samples



########## SHARP functions ##########

# Priority of GEX data if multiple outputs are found in db
sharp_wl_priority_map = {
        m: ["SEQC", "CR_GEX"] for m in ["Hashtag", "CiteSeq"]
    }
# File patterns to search for in S3 for each accompanying pipeline
sharp_wl_pattern_map = {
    "SEQC": "_dense.csv$",
    "CR_GEX": "/filtered_feature_bc_matrix/barcodes.tsv.gz$",
    "CR_ATAC": "/filtered_peak_bc_matrix/barcodes.tsv"
}
sharp_wl_method_map = {
    "SEQC": "SeqcDenseCountsMatrixCsv",
    "CR_GEX": "10x",
    "CR_ATAC": "10x",
}
# Names of FASTQ inputs in WDL; order is same as fastq_file_ids
# TODO: Ask to change all inputs to "fastq{file_id}" or "uriFastq{file_id}"
sharp_fastq_inputs_map = {
    m: ["uriFastqR1", "uriFastqR2"] for m in ["Hashtag", "CiteSeq"]
}

# Get s3 path of existing GEX analysis files
from mysql.connector import connect, Error

def get_wl_dir(sample_id, creds):
    
    user = creds['user']
    password = creds['password']

    try:
        table_sample_data = "peer_lab_db.sample_data"
        table_stats_data = "peer_lab_db.stats_data"
        table_stats_data = "peer_lab_db.stats_data"
        table_hashtag_lib = "peer_lab_db.hashtag_lib"
        table_genome_index = "peer_lab_db.genome_index"
        table_sc_tech = "peer_lab_db.sc_tech"
        query = f"""
        SELECT {table_stats_data}.analysis_storage
        FROM {table_sample_data}
        LEFT JOIN {table_stats_data} 
        ON {table_stats_data}.sampleData_id = {table_sample_data}.id
        LEFT JOIN {table_hashtag_lib}
        ON {table_hashtag_lib}.sampleData_id = {table_sample_data}.id
        LEFT JOIN {table_genome_index}
        ON {table_genome_index}.id = {table_hashtag_lib}.genomeIndex_id
        LEFT JOIN {table_sc_tech}
        ON {table_sc_tech}.id = {table_genome_index}.scTech_id
        WHERE {table_sample_data}.id = {sample_id}
        """
        result = execute_query(query, creds)[0][0]
        if result: 
            return result
        # As backup, get AWS storage location directly from sample_data
        else:
            query = f"""
            SELECT AWS_storage
            FROM {table_sample_data}
            WHERE {table_sample_data}.id = {sample_id}
            """
            result = execute_query(query, creds)[0][0]
            return result
    except Error as e:
        print(f"Error: {e}")

# Get white list method and associated file
# Throws exception if no white list exists
def get_wl_params(
    sample_id: str,
    creds,
    prefix,
    wl_dir
):
    
    user = creds['user']
    password = creds['password']

    wl_params = dict()

    # wl_dir = get_wl_dir(sample_id, creds)
    wl_patterns = [sharp_wl_pattern_map[p] for p in sharp_wl_priority_map[prefix]]

    try:
        # Check white list file exists before loading info from database
        assert wl_dir, f"Empty analysis storage for sample id {sample_id}"
        _, bucket, key, _, _ = urllib.parse.urlsplit(wl_dir)
        # White list file and method is first entry found on S3 
        wl = pd.DataFrame(
            [get_s3_objects(bucket, key.strip("/"), re.compile(p)) for p in wl_patterns],
            index = sharp_wl_priority_map[prefix],
        ).dropna(how="all")
        try:
            wl_key = wl.iloc[0,0] # if empty, missing white list file
            wl_params["uri"] = os.path.join("s3://", bucket, wl_key)
            wl_params["method"] = sharp_wl_method_map[wl.index[0]]
        except IndexError:
            logging.error(
                "Path to barcodes or counts matrix of GEX data is missing!"
            )
            return

    except AssertionError:
        logging.warning(f"Path to GEX output results is missing for {sample_id}!")
        return

    return wl_params


def get_bc_params(
    sample_id,
    creds,
):
    user = creds['user']
    password = creds['password']

    bc_params = dict()

    # JSON of bc and UMI positions are stored in database
    # First check dense matrix exists before loading JSON from database
    bc_json = get_bc_json(sample_id, creds)
    bc_pos = json.loads(bc_json)
    bc_params["cb"] = bc_pos["cellbarcode"]
    bc_params["umi"] = bc_params["cb"] + bc_pos["UMIs"]

    # Get bc sequence data from database
    bcs = get_bcs(sample_id, creds)
    if not bcs:
        logging.warning(f"Barcodes data Empty:\n\t {db_connect.cur.statement}")
        return
    for bc in bcs:
        try:
            assert bc[0], "AssertionError: Missing sequence barcodes!"
            assert bc[1], "AssertionError: Missing barcode IDs"
        except AssertionError as err:
            logging.warning(f"{err}:\n\t {db_connect.cur.statement}")
            return

    barcodes = pd.DataFrame(bcs, columns=["sequence", "code", "label", "bp_shift"])
    conjugation = barcodes["code"].str.get(0)
    if conjugation.nunique() != 1:
        logging.warning(
            f"Sample has multiple hashtag barcode categories and will not be processed!"
        )
        return
    else:
        bc_params["conjugation"] = conjugation.values[0]

    if barcodes["bp_shift"].nunique() != 1:
        logging.warning(
            f"Sample {sample_id} has hashtag barcode categories, with bp-shift length/s "
            f"{barcodes['bp_shift'].unique()}, and will not be processed!"
        )
        return
    else: 
        bc_params["bp_shift"] = int(barcodes["bp_shift"][0])
        bc_params["seq_length"] = bc_params["bp_shift"] + barcodes["sequence"].apply(len).max()

    return bc_params


# Get bc sequence data from database
def get_bcs(sample_id, creds):
    user = creds['user']
    password = creds['password']
    
    try:
        table_sample_data = "peer_lab_db.sample_data"
        table_hashtag_barcodes = "peer_lab_db.hashtag_barcodes"
        table_hashtags = "peer_lab_db.hashtags"
        query = f"""
        SELECT barcode_sequence, concat(substring(category, -1), barcode), 
        demultiplex_label, bp_shift FROM {table_hashtags} 
        LEFT JOIN {table_hashtag_barcodes} 
        ON {table_hashtag_barcodes}.id = {table_hashtags}.hashtagBarcodes_id 
        WHERE {table_hashtags}.sampleData_id = {sample_id}
        """
        result = execute_query(query, creds)
        return result
    except Error as e:
        print(f"Error: {e}")

# Get bc and UMI positions from database stored in JSON format
def get_bc_json(sample_id, creds):
    
    user = creds['user']
    password = creds['password']

    try:
        table_sample_data = "peer_lab_db.sample_data"
        table_stats_data = "peer_lab_db.stats_data"
        table_stats_data = "peer_lab_db.stats_data"
        table_hashtag_lib = "peer_lab_db.hashtag_lib"
        table_genome_index = "peer_lab_db.genome_index"
        table_sc_tech = "peer_lab_db.sc_tech"
        query = f"""
        SELECT barcodes
        FROM {table_sample_data}
        LEFT JOIN {table_stats_data} 
        ON {table_stats_data}.sampleData_id = {table_sample_data}.id
        LEFT JOIN {table_hashtag_lib}
        ON {table_hashtag_lib}.sampleData_id = {table_sample_data}.id
        LEFT JOIN {table_genome_index}
        ON {table_genome_index}.id = {table_hashtag_lib}.genomeIndex_id
        LEFT JOIN {table_sc_tech}
        ON {table_sc_tech}.id = {table_genome_index}.scTech_id
        WHERE {table_sample_data}.id = {sample_id}
        """
        result = execute_query(query, creds)[0][0]
        return result
    except Error as e:
        print(f"Error: {e}")
        
# Get fastq file paths on S3 for each file id
# Returns dictionary from id to s3 path
# Throws exception if FASTQs don't exist for any id
def get_denseCountMatrix(
    path: str, # path to directory containing FASTQ files
):
    _, bucket, key, _, _ = urllib.parse.urlsplit(path)
    results = get_s3_objects(
            bucket, key.lstrip("/"),
            re.compile(f"_dense.csv$")
        )
    whitelist = []
    for result in results:
        whitelist.append(os.path.join("s3://", bucket, result))
    whitelist.sort()
    return whitelist


# Function to reformat barcode labels for Sharp
def reformat_bc_label(label):
    label = label.encode('ascii', 'namereplace').decode()
    label = label.replace("\\N", "").replace(" ", "_")
    return label


########## Misc functions ##########

# Extract FASTQ sample name from list of files
# Note: FASTQ name is file name up to lane id (e.g. L001, L002, etc.)
def get_fastqs_name(fastqs):
    fastq_name_re = r".*/(.*)_S\d+_L\d{3}_[A-Za-z]\d_\d{3}.fastq.gz$"
    fastq_names = [re.match(fastq_name_re, x)[1] for x in fastqs]
    assert len(set(fastq_names)) == 1 # make sure all names are same
    return fastq_names[0]


########## Run workflow ##########

# Numpy encoder for JSON from pandas series
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)
        
def run(
    workflow_path: str,
    execp: str,
    secrets: str,
    inputs: str,
    labels: str,
    options: str,
):
    # change working directory to the pipeline package
    oldwd = os.getcwd()
    os.chdir(workflow_path)

    # execute the pipeline command
    cmd = f"{workflow_path}/{execp} -k {secrets} -i {inputs} -l {labels} -o {options}"
    var = subprocess.run(shlex.split(cmd), universal_newlines=True, capture_output=True)
    out = var.__dict__

    # change working directory back
    os.chdir(oldwd)

    return out


# Get fastq file paths on S3 for each file id
# Returns dictionary from id to s3 path
# Throws exception if FASTQs don't exist for any id
def get_mito_whitelist(
    path: str, # path to directory containing FASTQ files
):
    _, bucket, key, _, _ = urllib.parse.urlsplit(path)
    results = get_s3_objects(
            bucket, key.lstrip("/"),
            re.compile(f".txt$")
        )
    whitelist = []
    for result in results:
        whitelist.append(os.path.join("s3://", bucket, result))
    whitelist.sort()
    return whitelist


# Get fastq file paths on S3 for each file id
# Returns dictionary from id to s3 path
# Throws exception if FASTQs don't exist for any id
def get_aws_file(
    path: str, # path to directory containing FASTQ files
    file_end: str # Extension of the file
):
    _, bucket, key, _, _ = urllib.parse.urlsplit(path)
    results = get_s3_objects(
            bucket, key.lstrip("/"),
            re.compile(f".{file_end}$")
        )
    whitelist = []
    for result in results:
        whitelist.append(os.path.join("s3://", bucket, result))
    whitelist.sort()
    return whitelist
