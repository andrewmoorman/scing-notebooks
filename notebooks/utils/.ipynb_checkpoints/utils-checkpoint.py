import re, subprocess, boto3, json, shlex, mysql, os, urllib, logging
import numpy as np
import pandas as pd
from pathlib import Path
import importlib
from mysql.connector import connect, Error


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


# from SCRIdb
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


def execute_query(query, user, password):
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


# Get fastq file paths on S3 for each file id
# Returns dictionary from id to s3 path
# Throws exception if FASTQs don't exist for any id
def get_fastqs(
    path: str, # path to directory containing FASTQ files
    fastq_file_ids: list = None, # FASTQ file ids needed for this run type (e.g. I1, R1, R2, etc.)
    folder: str = "",
):
    fastq_map = dict()
    _, bucket, key, _, _ = urllib.parse.urlsplit(f"{path}/{folder}")
    # User may specify exactly which files are needed
    if fastq_file_ids:
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
    # Default: get all FASTQs
    else:
        files = get_s3_objects(
                bucket, key.lstrip("/"),
                re.compile(r"_\d{3}.fastq.gz$")
        )
        fastq_map["All"] = [os.path.join("s3://", bucket, str(f)) for f in files]

    return fastq_map


# Extract FASTQ sample name from list of files
# Note: FASTQ name is file name up to lane id (e.g. L001, L002, etc.)
def get_fastqs_name(fastqs):
    fastq_name_re = r".*/(.*)_S\d+_L\d{3}_[A-Za-z]\d_\d{3}.fastq.gz$"
    fastq_names = [re.match(fastq_name_re, x)[1] for x in fastqs]
    assert len(set(fastq_names)) == 1 # make sure all names are same
    return fastq_names[0]


# Get species from database for given sample
def get_species(sample_id, user, password):
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
        WHERE {table_sample_data}.id = {sample_id}
        """
        result = execute_query(query, user, password)[0][0]
        return result.lower()
    except Error as e:
        print(f"Error: {e}")


def get_sc_tech(sample_id, user, password):
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
        WHERE {table_sample_data}.id = {sample_id}
        """
        result = execute_query(query, user, password)[0][0]
        return result
    except Error as e:
        print(f"Error: {e}")


def get_sample_id(sample_name, user, password):
    try:
        table_sample_data = "peer_lab_db.sample_data"
        query = f"""
        SELECT {table_sample_data}.id
        FROM {table_sample_data}
        WHERE {table_sample_data}.Sample="{sample_name}"
        """
        result = execute_query(query, user, password)[0][0]
        return result
    except Error as e:
        print(f"Error: {e}")


def get_project_id(sample_id, user, password):
    try:
        table_sample_data = "peer_lab_db.sample_data"
        table_project_data = "peer_lab_db.project_data"
        query = f"""
        SELECT {table_project_data}.projectName
        FROM {table_project_data}
        LEFT JOIN {table_sample_data}
        ON {table_project_data}.id = {table_sample_data}.projectData_id
        WHERE {table_sample_data}.id = {sample_id}
        """
        result = execute_query(query, user, password)[0][0]
        return result
    except Error as e:
        print(f"Error: {e}")


def get_SEQC_version(loc):
    try:
        cmd = f"aws s3 cp {loc}/seqc-results/seqc_log.txt -"
        out = subprocess.run(
            shlex.split(cmd), universal_newlines=True, capture_output=True
            ).__dict__["stdout"]
        version = re.match(r".*SEQC=v(\d+\.\d+\.\d+).*", out)[1]
        return version
    except:
        return "N/A"


def get_file_prefix(loc):
    try:
        cmd = f"aws s3 ls {loc}/seqc-results/"
        out = subprocess.run(
            shlex.split(cmd), universal_newlines=True, capture_output=True
            ).__dict__["stdout"]

        # Note: I'm expecting the aligned bam file to be in loc
        bam_pattern = re.compile(r"(.*)_Aligned\.out\.bam$")
        filename = list(filter(bam_pattern.match, out.split()))[0]
        file_prefix = re.match(bam_pattern, filename)[1]
        return file_prefix
    except:
        raise ValueError(f"BAM file not found in {loc}")


def get_cr_reference(sample_id, prefix, user, password):
    # Get species from database to decide reference
    species = get_species(sample_id, user, password)
    
    # Map to reference locations
    try:
        with open("utils/assemblies-data.json") as f:
            cr_assemblies = json.load(f)['assembliesCellRanger']
        return cr_assemblies[prefix][species]
    except:
        raise ValueError(f"Unknown Species: {species}")


def get_bc_whitelist(sample_id, user, password):
    # Get version from database to decide whitelist
    sc_tech = get_sc_tech(sample_id, user, password)
    
    # Map to reference locations
    if "V3" in sc_tech:
        return "s3://seqc-public/barcodes/ten_x_v3/flat/3M-february-2018.txt"
    elif "V2" in sc_tech:
        return "s3://seqc-public/barcodes/ten_x_v2/flat/737K-august-2016.txt"
    else:
        raise ValueError(f"Unknown Technology: {sc_tech}")


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


# Get bc sequence data from database
def get_bcs(sample_id, user, password):
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
        result = execute_query(query, user, password)
        return result
    except Error as e:
        print(f"Error: {e}")


# Create csv files and upload to S3
# Note: follow CellRanger instructions for naming columns:
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi#examples
def get_cmo_files(
    samples: pd.DataFrame,
    user: str,
    password: str,
):
    cmo_files = dict()
    for name, sample in samples.iterrows():

        # Get barcodes from database
        bcs = pd.DataFrame.from_records(
            get_bcs(sample['Sample_ID'], user, password),
            columns=["sequence", "id", "sample_id", "bp_shift"],
        )

        # CMO map file
        cmo_map = bcs[["sample_id", "id"]].copy().rename(
            {"id": "cmo_ids"}, axis=1,
        )
        cmo_map["sample_id"] = cmo_map["sample_id"].str.replace(" ","_")

        # CMO reference file
        cmo_ref = bcs[["id", "sequence"]].copy()
        cmo_ref["name"] = cmo_ref["id"]
        cmo_ref["read"] = "R2"
        cmo_ref["pattern"] = "5P(BC)"
        cmo_ref["feature_type"] = "Multiplexing Capture"
        order = ["id", "name", "read", "pattern", "sequence", "feature_type"]

        cmo_files[name] = (cmo_map, cmo_ref[order])
    
    return cmo_files
