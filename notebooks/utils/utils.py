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
    'CellRangerGex': ['I1','R1','R2'],
    'MitoTracing': ['R1', 'R2']
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
        if prefix.startswith('CellRanger'):
            if not row['reference'].startswith('https'):

                if not row['species'] in ['human', 'mouse']:
                    print(f'{sample} species unknown')
                    samples.loc[sample, 'reference'] = np.nan
                else:
                    samples.loc[sample, 'reference'] = reference_map[prefix][row['species']]
    return samples


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