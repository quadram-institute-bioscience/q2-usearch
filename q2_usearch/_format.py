# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
# Copyright (c) 2024, Thanh Le Viet, Core Bioinformatics, Quadram Institute Bioscience
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import csv

# import numpy as np
import qiime2.plugin.model as model
from qiime2.plugin import ValidationError


class USEARCHTextFile(model.TextFileFormat):
    def validate(self, level):
        pass  # Basic validation could be added here if needed


class UchimeStatsFmt(model.TextFileFormat):

    def _check_n_records(self, n):
        with open(str(self)) as fh:
            csv_reader = csv.reader(fh, delimiter="\t")
            for i, row in enumerate(csv_reader):
                if i == n:
                    break
                else:
                    if len(row) != 18:
                        raise ValidationError(
                            "Incorrect number of fields detected on line %d."
                            " Should be exactly 18." % (i + 1)
                        )

    def _validate_(self, level):
        record_count_map = {"min": 5, "max": float("inf")}
        self._check_n_records(record_count_map[level])


UchimeStatsDirFmt = model.SingleFileDirectoryFormat(
    "UchimeStatsDirFmt", "stats.tsv", UchimeStatsFmt
)

USEARCHDirFmt = model.SingleFileDirectoryFormat(
    "USEARCHDirFmt", "data.txt", USEARCHTextFile
)
