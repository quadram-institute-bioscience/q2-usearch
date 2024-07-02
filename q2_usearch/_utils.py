# ----------------------------------------------------------------------------
# Copyright (c) 2024, Thanh Le Viet, Core Bioinformatics, Quadram Institute Bioscience
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import subprocess
import sqlite3

import skbio


class USearchError(Exception):
    """Exception raised for errors in the USEARCH process.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message="USEARCH command failed"):
        self.message = message
        super().__init__(self.message)


def validate_params(params: list):
    if any(param < 0 for param in params):
        raise ValueError("The parameter must be greater than or equal to 0.0.")


def run_command(cmd, verbose=True):
    """
    Execute a command line command from within Python.

    This function runs a specified command using the subprocess module. It optionally prints
    verbose messages about the execution status, including the command itself and whether
    it executed successfully or failed.

    Parameters:
    - cmd (list): The command to be executed, specified as a list of strings where the first
                  element is the command and the subsequent elements are the arguments.
    - verbose (bool, optional): If True, prints messages about the command execution, including
                                the command itself and its success or failure status. Defaults to True.

    Raises:
    - subprocess.CalledProcessError: If the command execution fails, this exception is raised
                                     and includes information about the failure.

    Note:
    - The command is executed without a shell, and therefore should be passed as a list of strings.
    - This function depends on the subprocess module for execution.
    - If verbose is True, messages are printed to stdout. In case of command execution failure,
      an error message is printed, and a subprocess.CalledProcessError is raised.
    """
    if verbose:
        print(
            "Running external command line application. This may print "
            "messages to stdout and/or stderr.\n"
            "The command being run is below. This command cannot "
            "be manually re-run as it will depend on temporary files that "
            "no longer exist.\n"
        )
        print(f"Command: {' '.join(cmd)}\n")

    try:
        subprocess.run(cmd, check=True)
        if verbose:
            print("Command executed successfully.")
    except subprocess.CalledProcessError as e:
        if verbose:
            print(f"Command execution failed with error: {e}")
        raise


def _fasta_with_sizes(input_fasta_fp, output_fasta_fp, table):
    """
    Add size annotations to sequences in a fasta file based on a given table.

    This function reads sequences from an input fasta file, retrieves their sizes from a provided table,
    and writes a new fasta file with size annotations included in the sequence headers.

    Parameters:
    - input_fasta_fp (str): File path to the input fasta file containing sequences.
    - output_fasta_fp (str): File path where the output fasta file with size annotations will be written.
    - table (object): A table object that contains size information for each sequence. The table must have
                      an 'ids' method that returns the sequence IDs and a 'sum' method that returns the size
                      for each sequence ID.

    Raises:
    - ValueError: If a sequence ID is found in the input fasta file that does not exist in the table, indicating
                  that the set of features in the sequences must be identical to the set of features in the table.

    Note:
    - The function assumes that the table object supports 'ids' and 'sum' methods for retrieving sequence IDs
      and their corresponding sizes.
    - It is expected that the sequence IDs in the fasta file and the table match exactly.
    """
    table_ids = table.ids(axis="observation")
    sizes = {id_: size for id_, size in zip(table_ids, table.sum(axis="observation"))}
    output_fasta_f = open(output_fasta_fp, "w")
    sequence_ids = set()
    for e in skbio.io.read(input_fasta_fp, constructor=skbio.DNA, format="fasta"):
        feature_id = e.metadata["id"]
        feature_seq = str(e)
        sequence_ids.add(feature_id)
        try:
            feature_size = sizes[feature_id]
        except KeyError:
            raise ValueError(
                "Feature %s is present in sequences, but not "
                "in table. The set of features in sequences must "
                "be identical to the set of features in table." % feature_id
            )
        output_fasta_f.write(
            ">%s;size=%d\n%s\n" % (feature_id, feature_size, feature_seq)
        )
    output_fasta_f.close()

    _error_on_nonoverlapping_ids(
        set(table_ids),
        sequence_ids,
        check_extra_table_ids=True,
        check_extra_sequence_ids=False,
    )


def _error_on_nonoverlapping_ids(
    table_ids, sequence_ids, check_extra_table_ids=True, check_extra_sequence_ids=True
):
    """
    Check for non-overlapping IDs between two sets and raise an error if found.

    This function compares two sets of IDs, typically representing features in a table and sequences,
    to ensure they are identical. If discrepancies are found (i.e., IDs present in one set but not the other),
    it raises a ValueError with a detailed message about the mismatch.

    Parameters:
    - table_ids (set): A set of IDs from the table.
    - sequence_ids (set): A set of IDs from the sequences.
    - check_extra_table_ids (bool, optional): If True, checks for IDs that are in the table but not in the sequences.
                                              Defaults to True.
    - check_extra_sequence_ids (bool, optional): If True, checks for IDs that are in the sequences but not in the table.
                                                 Defaults to True.

    Raises:
    - ValueError: If there are IDs present in one set but not the other, depending on the flags
                  `check_extra_table_ids` and `check_extra_sequence_ids`.

    Note:
    - This function is typically used to ensure that the features represented in a table match exactly
      with those in a set of sequences, which is a common requirement in bioinformatics analyses.
    """
    if check_extra_table_ids:
        extra_table_ids = table_ids - sequence_ids
        if len(extra_table_ids):
            raise ValueError(
                "Some feature ids are present in table, but not "
                "in sequences. The set of features in sequences "
                "must be identical to the set of features in "
                "table. Feature ids present in table but not "
                f"sequences are: {', '.join(extra_table_ids)}"
            )

    if check_extra_sequence_ids:
        extra_sequence_ids = sequence_ids - table_ids
        if len(extra_sequence_ids):
            raise ValueError(
                "Some feature ids are present in sequences, but "
                "not in table. The set of features in sequences "
                "must be identical to the set of features in "
                "table. Feature ids present in sequences but not "
                f"table are: {', '.join(extra_sequence_ids)}"
            )


def _uc_to_sqlite(uc):
    """Parse uc-style file into a SQLite in-memory database.

    This populates an in-memory database with the following schema (displayed
    below with dummy data):

        feature_id | cluster_id | count
        -----------|------------|-------
        feature1   | r1         | 204
        feature2   | r2         | 4
        feature3   | r1         | 15
        feature4   | r2         | 24
        feature5   | r2         | 16
    """
    conn = sqlite3.connect(":memory:")
    c = conn.cursor()
    # The PK constraint ensures that there are no duplicate Feature IDs
    c.execute(
        "CREATE TABLE feature_cluster_map (feature_id TEXT PRIMARY KEY,"
        "cluster_id TEXT NOT NULL, count INTEGER);"
    )
    c.execute("CREATE INDEX idx1 ON " "feature_cluster_map(feature_id, cluster_id);")
    conn.commit()
    insert_stmt = "INSERT INTO feature_cluster_map VALUES (?, ?, ?);"

    for line in uc:
        line = line.strip()
        if len(line) == 0 or line.startswith(b"#"):
            continue
        else:
            fields = line.split(b"\t")
            if fields[0] == b"S":
                sequence_id = fields[8].decode("utf-8").split(";")[0]
                c.execute(insert_stmt, (sequence_id, sequence_id, None))
            elif fields[0] == b"H":
                centroid_id = fields[9].decode("utf-8").split(";")[0]
                sequence_id = fields[8].decode("utf-8").split(";size=")
                if len(sequence_id) == 2:
                    sequence_id, count = sequence_id[0], sequence_id[1]
                else:
                    sequence_id, count = sequence_id[0], "1"
                c.execute(insert_stmt, (sequence_id, centroid_id, count))
            else:
                pass
    conn.commit()
    return conn


def _collapse_f_from_sqlite(conn):
    c = conn.cursor()
    # This query produces the following results (displayed below with dummy
    # data):
    # feature_id | cluster_id
    # -----------|------------
    # feature1   | r1
    # feature2   | r2
    # feature3   | r1
    # feature4   | r2
    # feature4   | r2
    c.execute("SELECT feature_id, cluster_id FROM feature_cluster_map;")
    id_to_centroid = dict(c.fetchall())

    if len(id_to_centroid) == 0:
        raise ValueError("No sequence matches were identified by vsearch.")

    def collapse_f(id_, x):
        return id_to_centroid[id_]

    return collapse_f


def _fasta_from_sqlite(conn, input_fasta_fp, output_fasta_fp):
    input_seqs = skbio.read(input_fasta_fp, format="fasta", constructor=skbio.DNA)
    c = conn.cursor()
    # Create a second in-memory table with the following schema (displayed
    # below with dummy data):
    # feature_id | sequence_string
    # -----------|------------------
    # feature1   | ACGTACGTACGTACGT
    # feature2   | GGGGAAAACCCCTTTT
    # feature3   | TCAGAAAATTTTTCAG
    # feature4   | AAAAAAAAAAAAAAAA
    # feature5   | GGGGGGGGGGGGGGGG
    c.execute(
        "CREATE TABLE rep_seqs (feature_id TEXT PRIMARY KEY, "
        "sequence_string TEXT NOT NULL);"
    )
    c.executemany(
        "INSERT INTO rep_seqs VALUES (?, ?);",
        [(seq.metadata["id"], str(seq)) for seq in input_seqs],
    )
    conn.commit()
    # Preemptively sort the table to deal with tie-breaking, later.
    # This is a table, not a view, because we want/need sqlite's rowid.
    c.execute(
        "CREATE TABLE sorted_feature_cluster_map AS "
        "SELECT * FROM feature_cluster_map ORDER BY cluster_id ASC,"
        "feature_id ASC;"
    )
    c.execute("CREATE INDEX idx2 ON " "sorted_feature_cluster_map(cluster_id, count);")
    conn.commit()
    # The results from this query should look like the following (displayed
    # below with dummy data):
    # cluster_id | sequence_string
    # -----------|------------------
    # r1         | ACGTACGTACGTACGT
    # r2         | AAAAAAAAAAAAAAAA
    c.execute(
        """SELECT fcm.cluster_id, rs.sequence_string, MAX(fcm.count)
                   FROM sorted_feature_cluster_map fcm
             INNER JOIN rep_seqs rs ON rs.feature_id = fcm.feature_id
               GROUP BY fcm.cluster_id
               ORDER BY fcm.cluster_id ASC;
    """
    )
    with open(output_fasta_fp, "w") as output_seqs:
        while True:
            partial_results = c.fetchmany(size=100)
            if partial_results:
                output_seqs.writelines(
                    [f">{i}\n{s}\n" for (i, s, _) in partial_results]
                )
            else:
                break
