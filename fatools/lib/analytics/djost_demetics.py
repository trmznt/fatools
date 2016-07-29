from fatools.lib.analytics.export import export_demetics
from fatools.lib.utils import cerr, cout, random_string
from subprocess import call
import numpy as np


def run_demetics(analytical_sets, dbh, tmp_dir, mode='d.jost'):

    # demetics output cannot be managed as it directly writes output file to current working directory
    # as such, we need to run demetics under a script that will change the current directory

    file_id = random_string(3)

    script_file = '%s/demetics-%s.r' % (tmp_dir, file_id)
    data_file = '%s/data-%s.txt' % (tmp_dir, file_id)

    with open(data_file, 'w') as dataout:
        export_demetics(analytical_sets, dbh, dataout)

    with open(script_file, 'w') as scriptout:
        scriptout.write(
            'setwd("%s")\n'
            'library(DEMEtics)\n'
            'dat <- read.table("%s", header=T)\n'
            'D.Jost("dat", bias="correct", object=TRUE,format.table=FALSE,pm="pairwise", statistics="CI", bt=100)\n'
            % (tmp_dir, data_file)
        )

    ok = call( [ 'Rscript', script_file] )

    if ok != 0:
        raise RuntimeError('Rscript for DEMEtics run unsuccessfully. Please inspect the error log.')

    raise RuntimeError('please inspect here: %s' % tmp_dir)