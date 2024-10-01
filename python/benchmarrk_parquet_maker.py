#!/usr/bin/env python
# K. A. Lacek - made base
# A. H. Sullivan - built in benchmark value conversions
# Sept 2024

# Functional for samplesheet, irma summary, alleles, indels
import pandas as pd
import argparse
from pathlib import Path
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from os.path import dirname, basename, isfile, getmtime
from glob import glob
import time
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument(
    "-p", "--path", help="the file path to where the IRMA folder is located"
)
parser.add_argument(
    "-f",
    "--file",
    help="input file for conversion to CDP-compatible parquet. If empty, script creates reads table, coverage table, and all alleles table with filenames inferred from runid",
)
parser.add_argument("-o", "--outputname", help="name of parquet output file")
parser.add_argument(
    "-r",
    "--runid",
    help="runid name, if empty: runid is the name of the directory in which the script is run",
)
parser.add_argument(
    "-i", "--instrument", help="sequencing instrument name. If empty, testInstrument"
)

inputarguments = parser.parse_args()

if inputarguments.path:
    wd_path = inputarguments.path
else:
    wd_path = Path.cwd()

if inputarguments.file:
    infi = inputarguments.file
else:
    infi = ""

if inputarguments.outputname:
    outfi = inputarguments.outputname
else:
    outfi = ""

if inputarguments.runid:
    run_id = inputarguments.runid
else:
    run_id = Path.cwd()
    run_id = str(run_id).split("/")[-1]

if inputarguments.instrument:
    instrument = inputarguments.instrument
else:
    instrument = "testInstrument"


def irmatable2df(irmaFiles):
    df = pd.DataFrame()
    for f in irmaFiles:
        sample = basename(dirname(dirname(f)))
        if "insertions" not in f:
            df_prime = pd.read_csv(f, sep="\t", index_col=False)
        else:
            df_prime = pd.read_csv(f, sep="\s+", index_col=False)
        df_prime.insert(loc=0, column="Sample", value=sample)
        df = pd.concat([df, df_prime])
    return df


def irma_reads_df(irma_path):
    readFiles = glob(irma_path + "/*/tables/READ_COUNTS.txt")
    df = pd.DataFrame()
    df = irmatable2df(readFiles)
    df["Stage"] = df["Record"].apply(lambda x: int(x.split("-")[0]))
    return df


def irma_coverage_df(irma_path):
    coverageFiles = glob(irma_path + "/*/tables/*a2m.txt")
    # a2msamples = [i.split('/')[-3] for i in coverageFiles]
    # otherFiles = [i for i in glob(irma_path+'/*/tables/*coverage.txt')]
    if len(coverageFiles) == 0:
        coverageFiles = glob(irma_path + "/*/tables/*coverage.txt")
    if len(coverageFiles) == 0:
        return "No coverage files found under {}/*/tables/".format(irma_path)
    df = irmatable2df(coverageFiles)

    return df


def irma_alleles_df(irma_path, full=False):
    alleleFiles = glob(irma_path + "/*/tables/*variants.txt")
    df = irmatable2df(alleleFiles)
    if not full:
        if "HMM_Position" in df.columns:
            ref_heads = [
                "Sample",
                "Reference_Name",
                "HMM_Position",
                "Position",
                "Total",
                "Consensus_Allele",
                "Minority_Allele",
                "Consensus_Count",
                "Minority_Count",
                "Minority_Frequency",
            ]
        else:
            ref_heads = [
                "Sample",
                "Reference_Name",
                "Position",
                "Total",
                "Consensus_Allele",
                "Minority_Allele",
                "Consensus_Count",
                "Minority_Count",
                "Minority_Frequency",
            ]
        df = df[ref_heads]
        df = df.rename(
            columns={
                "Reference_Name": "Reference",
                "HMM_Position": "Reference Position",
                "Position": "Sample Position",
                "Total": "Coverage",
                "Consensus_Allele": "Consensus Allele",
                "Minority_Allele": "Minority Allele",
                "Consensus_Count": "Consensus Count",
                "Minority_Count": "Minority Count",
                "Minority_Frequency": "Minority Frequency",
            }
        )
        df["Minority Frequency"] = df[["Minority Frequency"]].applymap(
            lambda x: float(f"{float(x):.{3}f}")
        )
    return df


# file I/O
def parquetify(table, outfi):
    pd.DataFrame.to_csv(table, "temp.csv", sep="\t", index=False, header=True)
    chunksize = 100_000
    # modified from https://stackoverflow.com/questions/26124417/how-to-convert-a-csv-file-to-parquet
    csv_stream = pd.read_csv(
        "temp.csv", sep="\t", chunksize=chunksize, low_memory=False
    )
    for i, chunk in enumerate(csv_stream):
        print("Chunk", i)
        if i == 0:
            # Guess the schema of the CSV file from the first chunk
            try:
                if "indel" in infi:
                    print("indel")
                    parquet_schema = pa.schema(
                        [
                            ("Sample", pa.string()),
                            ("Sample - Upstream Position", pa.int64()),
                            ("Reference", pa.string()),
                            ("Context", pa.string()),
                            ("Length", pa.int64()),
                            ("Insert", pa.string()),
                            ("Count", pa.int64()),
                            ("Upstream Base Coverage", pa.int64()),
                            ("Frequency", pa.float64()),
                            ("runid", pa.string()),
                            ("instrument", pa.string()),
                        ]
                    )
                elif "benchmark" in infi:
                    print("benchmark")
                    parquet_schema = pa.schema(
                        [
                            ("task_id", pa.string()),
                            ("hash", pa.string()),
                            ("native_id", pa.string()),
                            ("name", pa.string()),
                            ("status", pa.string()),
                            ("exit", pa.string()),
                            ("submit", pa.timestamp("s", tz="UTC")),
                            ("duration", pa.float32()),
                            ("realtime", pa.float32()),
                            ("%cpu", pa.float32()),
                            ("peak_rss", pa.float32()),
                            ("peak_vmem", pa.float32()),
                            ("rchar", pa.float32()),
                            ("wchar", pa.float32()),
                            ("runid", pa.string()),
                            ("instrument", pa.string()),
                        ]
                    )
                elif "indel" not in infi or "benchmark" not in infi:
                    parquet_schema = pa.Table.from_pandas(df=chunk).schema
                    print("wrong")
            except:
                parquet_schema = pa.Table.from_pandas(df=chunk).schema
            # Open a Parquet file for writing
            parquet_writer = pq.ParquetWriter(
                outfi, parquet_schema, compression="snappy", version="1.0"
            )
        # Write CSV chunk to the parquet file
        table = pa.Table.from_pandas(chunk, schema=parquet_schema)
        parquet_writer.write_table(table)

    parquet_writer.close()


if ".csv" in infi:
    table = pd.read_csv(infi, header=0)
elif ".xls" in infi:
    table = pd.read_excel(infi, header=0, engine="openpyxl")
elif "run_info.txt" in infi:
    table = pd.read_csv(infi, sep="\t", header=0)
    table["runid"] = run_id
    table["instrument"] = instrument
    timestamp = getmtime(infi)
    table["timestamp"] = datetime.fromtimestamp(timestamp)
    parquetify(table, outfi)
    exit()
elif "benchmark.txt" in infi:
    table = pd.read_csv(infi, sep="\t", header=0)
    table = table.replace("%", "", regex=True)
    table["task_id"] = table["task_id"].apply(lambda x: "'" + str(x) + "'")
    table["native_id"] = table["native_id"].apply(lambda x: "'" + str(x) + "'")
    table["exit"] = table["exit"].apply(lambda x: "'" + str(x) + "'")
    ## converting to timestamp format to pyarrow parquet conversion
    for row in table["submit"]:
        row = row.split(".", 1)[0]
        timestamp = datetime.strptime(row, "%Y-%m-%d %H:%M:%S").timestamp()
        table["submit"] = int(timestamp)
    ##converting peak_rss to bytes 
    table[['peak_rss_num', 'peak_rss_size']] = table['peak_rss'].str.split(' ', n=1, expand=True)
    table.loc[table['peak_rss_size'] == 'B', 'peak_rss'] = (table['peak_rss_num'].astype(float))*1
    table.loc[table['peak_rss_size'] == 'KB', 'peak_rss'] = (table['peak_rss_num'].astype(float))*1000
    table.loc[table['peak_rss_size'] == 'MB', 'peak_rss'] = (table['peak_rss_num'].astype(float))*1000000
    table.loc[table['peak_rss_size'] == 'GB', 'peak_rss'] = (table['peak_rss_num'].astype(float))*1000000000
    table.loc[table['peak_rss_size'] == 'TB', 'peak_rss'] = (table['peak_rss_num'].astype(float))*1000000000000
    ##converting peak_vmem to bytes 
    table[['peak_vmem_num', 'peak_vmem_size']] = table['peak_vmem'].str.split(' ', n=1, expand=True)
    table.loc[table['peak_vmem_size'] == 'B', 'peak_vmem'] = (table['peak_vmem_num'].astype(float))*1
    table.loc[table['peak_vmem_size'] == 'KB', 'peak_vmem'] = (table['peak_vmem_num'].astype(float))*1000
    table.loc[table['peak_vmem_size'] == 'MB', 'peak_vmem'] = (table['peak_vmem_num'].astype(float))*1000000
    table.loc[table['peak_vmem_size'] == 'GB', 'peak_vmem'] = (table['peak_vmem_num'].astype(float))*1000000000
    table.loc[table['peak_vmem_size'] == 'TB', 'peak_vmem'] = (table['peak_vmem_num'].astype(float))*1000000000000
    ##converting rchar to bytes 
    table[['rchar_num', 'rchar_size']] = table['rchar'].str.split(' ', n=1, expand=True)
    table.loc[table['rchar_size'] == 'B', 'rchar'] = (table['rchar_num'].astype(float))*1
    table.loc[table['rchar_size'] == 'KB', 'rchar'] = (table['rchar_num'].astype(float))*1000
    table.loc[table['rchar_size'] == 'MB', 'rchar'] = (table['rchar_num'].astype(float))*1000000
    table.loc[table['rchar_size'] == 'GB', 'rchar'] = (table['rchar_num'].astype(float))*1000000000
    table.loc[table['rchar_size'] == 'TB', 'rchar'] = (table['rchar_num'].astype(float))*1000000000000
    ##converting wchar to bytes 
    table[['wchar_num', 'wchar_size']] = table['wchar'].str.split(' ', n=1, expand=True)
    table.loc[table['wchar_size'] == 'B', 'wchar'] = (table['wchar_num'].astype(float))*1
    table.loc[table['wchar_size'] == 'KB', 'wchar'] = (table['wchar_num'].astype(float))*1000
    table.loc[table['wchar_size'] == 'MB', 'wchar'] = (table['wchar_num'].astype(float))*1000000
    table.loc[table['wchar_size'] == 'GB', 'wchar'] = (table['wchar_num'].astype(float))*1000000000
    table.loc[table['wchar_size'] == 'TB', 'wchar'] = (table['wchar_num'].astype(float))*1000000000000
    table['wchar'].astype(float)

    ##converting duration time
    table['duration'] = table['duration'].str.replace(r'ms', 'x', regex=True)
    table[['duration_mix', 'duration_s_2']] = table['duration'].str.split(' ', n=1, expand=True)
    table['duration_m'] = table['duration_mix'].apply(lambda x: x if 'm' in x else None)
    table['duration_s_1'] = table['duration_mix'].apply(lambda x: x if 's' in x else None)
    table['duration_ms'] = table['duration_mix'].apply(lambda x: x if 'x' in x else None)
    table['duration_s_1'] = table['duration_s_1'].str.replace(r's', '', regex=True)
    table['duration_s_2'] = table['duration_s_2'].str.replace(r's', '', regex=True)
    table['duration_ms'] = table['duration_ms'].str.replace(r'x', '', regex=True)
    table['duration_m'] = table['duration_m'].str.replace(r'm', '', regex=True)
    table.fillna(0, inplace=True)
    table['duration_s_1'] = (table['duration_s_1'].astype(float))/1
    table['duration_s_2'] = (table['duration_s_2'].astype(float))/1
    table['duration_ms'] = (table['duration_ms'].astype(float))/1000
    table['duration_m'] = (table['duration_m'].astype(float))*60
    table['duration'] = table['duration_m'].astype(float) + table['duration_s_1'].astype(float) + table['duration_s_2'].astype(float) + table['duration_ms'].astype(float)

    ##converting realtime time
    table['realtime'] = table['realtime'].str.replace(r'ms', 'x', regex=True)
    table[['realtime_mix', 'realtime_s_2']] = table['realtime'].str.split(' ', n=1, expand=True)
    table['realtime_m'] = table['realtime_mix'].apply(lambda x: x if 'm' in x else None)
    table['realtime_s_1'] = table['realtime_mix'].apply(lambda x: x if 's' in x else None)
    table['realtime_ms'] = table['realtime_mix'].apply(lambda x: x if 'x' in x else None)
    table['realtime_s_1'] = table['realtime_s_1'].str.replace(r's', '', regex=True)
    table['realtime_s_2'] = table['realtime_s_2'].str.replace(r's', '', regex=True)
    table['realtime_ms'] = table['realtime_ms'].str.replace(r'x', '', regex=True)
    table['realtime_m'] = table['realtime_m'].str.replace(r'm', '', regex=True)
    table.fillna(0, inplace=True)
    table['realtime_s_1'] = (table['realtime_s_1'].astype(float))/1
    table['realtime_s_2'] = (table['realtime_s_2'].astype(float))/1
    table['realtime_ms'] = (table['realtime_ms'].astype(float))/1000
    table['realtime_m'] = (table['realtime_m'].astype(float))*60
    table['realtime'] = table['realtime_m'].astype(float) + table['realtime_s_1'].astype(float) + table['realtime_s_2'].astype(float) + table['realtime_ms'].astype(float)



elif ".txt" in infi or ".tsv" in infi:
    table = pd.read_csv(infi, sep="\t", header=0)


elif ".fasta" in infi:
    seq_dict = {}
    i = 0
    with open(infi, "r") as infi:
        for line in infi:
            if line[0] == ">":
                i += 1
                key = line.strip(">").split("|")[0]
                segment = line.strip().split("|")[1]
                try:
                    qc_failure = line.strip().split("|")[2]
                except:
                    qc_failure = "pass"
                seq_dict[i] = [key, segment, qc_failure]
            else:
                value = line.strip()
                seq_dict[i].append(value)
    table = pd.DataFrame.from_dict(
        seq_dict,
        orient="index",
        columns=["sample_id", "reference", "qc_decision", "sequence"],
    )

elif infi == "":
    readstable = irma_reads_df(wd_path + "/IRMA")
    covtable = irma_coverage_df(wd_path + "/IRMA")
    allelestable = irma_alleles_df(wd_path + "/IRMA")
    for t in (
        [readstable, f"{run_id}_reads.parq"],
        [covtable, f"{run_id}_coverage.parq"],
        [allelestable, f"{run_id}_alleles.parq"],
    ):
        t[0]["runid"] = run_id
        t[0]["instrument"] = instrument
        parquetify(t[0], t[1])
    exit()


table["runid"] = run_id

table["instrument"] = instrument

parquetify(table, outfi)
