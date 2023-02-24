import os
import pathlib
import subprocess
import sys

import zlib


def run_check_errors(cmd):
    """Run a system command, and exit if an error occurred, otherwise continue"""
    if type(cmd) == str:
        cmd = cmd.split()
    output = subprocess.run(cmd, capture_output=True, text=True)
    if output.stderr != "":
        print_cmd = " ".join(map(str, cmd))
        sys.exit(
            f"The error {output.stderr} was generated when running {print_cmd}. Exiting."
        )
    return


def collate(
    batch,
    config,
    plate,
    base_directory="../..",
    column=None,
    munge=False,
    csv_dir="analysis",
    aws_remote=None,
    aggregate_only=False,
    tmp_dir="/tmp",
    overwrite=False,
    add_image_features=True,
    image_feature_categories=["Granularity", "Texture", "ImageQuality", "Threshold"],
    printtoscreen=True,
):
    """Collate the CellProfiler-created CSVs into a single SQLite file by calling cytominer-database

    Parameters
    ----------
    batch : str
        Batch name to process
    config : str
        Config file to pass to cytominer-database
    plate : str
        Plate name to process
    base_directory : str, default "../.."
        Base directory for subdirectories containing CSVs, backends, etc; in our preferred structure, this is the "workspace" directory
    column : str, optional, default None
        An existing column to be explicitly copied to a new column called Metadata_Plate if no Metadata_Plate column already explicitly exists
    munge : bool, default False
        Whether munge should be passed to cytominer-database, if True cytominer-database will expect a single all-object CSV; it will split each object into its own table
    csv_dir : str, default 'analysis'
        The directory under the base directory where the analysis CSVs will be found. If running the analysis pipeline, this should nearly always be "analysis"
    aws_remote : str, optional, default None
        A remote AWS prefix, if set CSV files will be synced down from at the beginning and to which SQLite files will be synced up at the end of the run
    aggregate_only : bool, default False
        Whether to perform only the aggregation of existent SQLite files and bypass previous collation steps
    tmp_dir: str, default '/tmp'
        The temporary directory to be used by cytominer-databases for output
    overwrite: bool, optional, default False
        Whether or not to overwrite an sqlite that exists in the temporary directory if it already exists
    add_image_features: bool, optional, default True
        Whether or not to add the image features to the profiles
    image_feature_categories: list, optional, default ['Granularity','Texture','ImageQuality','Count','Threshold']
        The list of image feature groups to be used by add_image_features during aggregation
    printtoscreen: bool, optional, default True
        Whether or not to print output to the terminal
    """

    from pycytominer.cyto_utils.cells import SingleCells

    # Set up directories (these need to be abspaths to keep from confusing makedirs later)
    input_dir = pathlib.Path(f"{base_directory}/analysis/{batch}/{plate}/{csv_dir}")
    backend_dir = pathlib.Path(f"{base_directory}/backend/{batch}/{plate}")
    cache_backend_dir = pathlib.Path(f"{tmp_dir}/backend/{batch}/{plate}")

    aggregated_file = pathlib.Path(f"{backend_dir}/{plate}.csv")
    backend_file = pathlib.Path(f"{backend_dir}/{plate}.sqlite")
    cache_backend_file = pathlib.Path(f"{cache_backend_dir}/{plate}.sqlite")

    if not aggregate_only:
        if os.path.exists(cache_backend_file):
            if not overwrite:
                sys.exit(
                    f"An SQLite file for {plate} already exists at {cache_backend_file} and overwrite is set to False. Terminating."
                )
            else:
                os.remove(cache_backend_file)

        for eachdir in [input_dir, backend_dir, cache_backend_dir]:
            if not os.path.exists(eachdir):
                os.makedirs(eachdir, exist_ok=True)

        if aws_remote:

            remote_input_dir = f"s3://cellpainting-gallery/cpg0016-jump/source_13/workspace/analysis/{batch}/{plate}/{csv_dir}"

            remote_backend_file = f"{aws_remote}/backend/{batch}/{plate}/{plate}.sqlite"

            remote_aggregated_file = f"{aws_remote}/backend/{batch}/{plate}/{plate}.csv"

            sync_cmd = f'aws s3 sync --exclude "*" --include "*/Cells.csv" --include "*/Nuclei.csv" --include "*/Cytoplasm.csv" --include "*/Image.csv" {remote_input_dir} {input_dir}  --profile jump-cp-role-jump-cellpainting'
            if printtoscreen:
                print(f"Downloading CSVs from {remote_input_dir} to {input_dir}")
            run_check_errors(sync_cmd)

        """ingest_cmd = [
            "cytominer-database",
            "ingest",
            input_dir,
            f"sqlite:///{cache_backend_file}",
            "-c",
            config,
        ]
        if not munge:
            # munge is by default True in cytominer-database
            ingest_cmd.append("--no-munge")
        if printtoscreen:
            print(f"Ingesting {input_dir}")
        run_check_errors(ingest_cmd)"""
        from sqlalchemy import create_engine
        from sqlalchemy.pool import NullPool
        import glob
        import pandas
        engine = create_engine(f"sqlite:///{cache_backend_file}", poolclass=NullPool)
        con = engine.connect()

        image_csv = glob.glob(os.path.join(input_dir,"**/Image.csv"))[0]
        compartment_csvs = [x for x in glob.glob(os.path.join(input_dir,"**/*.csv")) if "Image" not in x]

        identifier = checksum(image_csv)
        print("Ingesting Image")
        #handle image CSV
        df_image = pandas.read_csv(image_csv)
        image_cols = list(df_image.columns)
        rename_dict = {}
        for eachcol in image_cols:
            if eachcol[:6] != "Image_":
                rename_dict[eachcol]=eachcol
            else:
                rename_dict[eachcol]=eachcol[6:]
        renamed_cols = list(rename_dict.values())
        renamed_cols.sort()
        image_new_cols = ["TableNumber"] + renamed_cols
        df_image.rename(columns=rename_dict,inplace=True)
        df_image["TableNumber"] = identifier
        df_image = df_image[image_new_cols]
        df_image.to_sql("Image",con,if_exists='append',index=False)

        for eachcompartment in ["Nuclei","Cells","Cytoplasm"]:
            print(f"Ingesting {eachcompartment}")
            comp_csv = [x for x in compartment_csvs if eachcompartment in x][0]
            df_comp = pandas.read_csv(comp_csv)
            comp_cols = list(df_comp.columns)
            # We want to keep these column names the same and first, everything else should get the compartment appended
            dont_adjust = ["ImageNumber","ObjectNumber","TableNumber"]
            rename_dict = {}
            for eachcol in compcols:
                if eachcol in dont_adust:
                    pass
                elif eachcol[:len(eachcompartment)+1]==f"{eachcompartment}_":
                    rename_dict[eachcol]=eachcol
                else:
                    rename_dict[eachcol]=f"{eachcompartment}_{eachcol}"
            renamed_cols = list(rename_dict.values())
            renamed_cols.sort()
            renamed_cols = dont_adjust + renamed_cols
            for eachcol in dont_adjust:
                rename_dict[eachcol]=eachcol
            df_comp["TableNumber"] = identifier
            df_comp.rename(columns=rename_dict,inplace=True)
            df_comp = df_comp[renamed_cols]
            print(f"Writing {eachcompartment}")
            df_comp.to_sql(eachcompartment,con,if_exists='append',index=False)

        con.close()

        if column:
            if print:
                print(f"Adding a Metadata_Plate column based on column {column}")
            alter_cmd = [
                "sqlite3",
                cache_backend_file,
                "ALTER TABLE Image ADD COLUMN Metadata_Plate TEXT;",
            ]
            run_check_errors(alter_cmd)
            update_cmd = [
                "sqlite3",
                cache_backend_file,
                f"UPDATE image SET Metadata_Plate ={column};",
            ]
            run_check_errors(update_cmd)

        if printtoscreen:
            print(f"Indexing database {cache_backend_file}")
        for eachcompartment in ["Cells", "Cytoplasm", "Nuclei", "Image"]:
            ix_cmd_compartment = [
                "sqlite3",
                cache_backend_file,
                f"CREATE INDEX IF NOT EXISTS ix_{eachcompartment.lower()}_TableNumber ON {eachcompartment}(TableNumber);",
            ]
            run_check_errors(ix_cmd_compartment)
        index_cmd_img = [
            "sqlite3",
            cache_backend_file,
            "CREATE INDEX IF NOT EXISTS table_image_idx ON Image(TableNumber, ImageNumber);",
        ]
        run_check_errors(index_cmd_img)
        for eachcompartment in ["Cells", "Cytoplasm", "Nuclei"]:
            index_cmd_compartment = [
                "sqlite3",
                cache_backend_file,
                f"CREATE INDEX IF NOT EXISTS table_image_object_{eachcompartment.lower()}_idx ON {eachcompartment}(TableNumber, ImageNumber, ObjectNumber);",
            ]
            run_check_errors(index_cmd_compartment)
        index_cmd_metadata = [
            "sqlite3",
            cache_backend_file,
            "CREATE INDEX IF NOT EXISTS plate_well_image_idx ON Image(Metadata_Plate, Metadata_Well);",
        ]
        run_check_errors(index_cmd_metadata)

        if aws_remote:

            if printtoscreen:
                print(f"Uploading {cache_backend_file} to {remote_backend_file}")
            cp_cmd = ["aws", "s3", "cp", cache_backend_file, remote_backend_file,"--profile","jump-cp-role-jump-cellpainting","--acl","bucket-owner-full-control","--metadata-directive","REPLACE"]
            run_check_errors(cp_cmd)

            if printtoscreen:
                print(
                    f"Removing analysis files from {input_dir} and {cache_backend_dir}"
                )
            import shutil

            shutil.rmtree(input_dir)

        if printtoscreen:
            print(f"Renaming {cache_backend_file} to {backend_file}")
        os.rename(cache_backend_file, backend_file)

    if printtoscreen:
        print(f"Aggregating sqlite:///{backend_file}")

    if aggregate_only and aws_remote:
        remote_backend_file = f"{aws_remote}/backend/{batch}/{plate}/{plate}.sqlite"

        remote_aggregated_file = f"{aws_remote}/backend/{batch}/{plate}/{plate}.csv"

        cp_cmd = ["aws", "s3", "cp", remote_backend_file, backend_file,"--profile","jump-cp-role-jump-cellpainting","--acl","bucket-owner-full-control","--metadata-directive","REPLACE"]
        if printtoscreen:
            print(
                f"Downloading SQLite files from {remote_backend_file} to {backend_file}"
            )
        run_check_errors(cp_cmd)

    if not os.path.exists(backend_file):
        sys.exit(f"{backend_file} does not exist. Exiting.")

    if add_image_features:
        pass
    else:
        image_feature_categories = None  # defensive but not sure what will happen if we give a list but set to False

    database = SingleCells(
        f"sqlite:///{backend_file}",
        aggregation_operation="mean",
        add_image_features=add_image_features,
        image_feature_categories=image_feature_categories,
    )
    database.aggregate_profiles(output_file=aggregated_file)

    if aws_remote:
        if printtoscreen:
            print(f"Uploading {aggregated_file} to {remote_aggregated_file}")
        csv_cp_cmd = ["aws", "s3", "cp", aggregated_file, remote_aggregated_file,"--profile","jump-cp-role-jump-cellpainting","--acl","bucket-owner-full-control","--metadata-directive","REPLACE"]
        run_check_errors(csv_cp_cmd)

        if printtoscreen:
            print(f"Removing backend files from {backend_dir}")
        import shutil

        shutil.rmtree(backend_dir)

###### CODE PULLED DIRCTLY FROM CYTOMINER-DATABASE

def checksum(pathname, buffer_size=65536):
    """
    Generate a 32-bit unique identifier for a file.
    
    :param pathname: input file
    :param buffer_size: buffer size   
    """
    with open(pathname, "rb") as stream:
        result = zlib.crc32(bytes(0))

        while True:
            buffer = stream.read(buffer_size)

            if not buffer:
                break

            result = zlib.crc32(buffer, result)

    return result & 0xffffffff
