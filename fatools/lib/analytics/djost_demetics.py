from fatools.lib.analytics.export import export_demetics
from fatools.lib.utils import cerr, cout, random_string
from subprocess import call
import numpy as np


def run_demetics(analytical_sets, dbh, tmp_dir, mode='d.jost'):

    file_id = random_string(3)

    script_file = '%s/demetics-%s.r' % (tmp_dir, file_id)
    data_file = '%s/data-%s.txt' % (tmp_dir, file_id)

    with open(data_file, 'w') as dataout:
        export_demetics(analytical_sets, dbh, dataout)

    with open(script_file, 'w') as scriptout:
        scriptout.write(
            'library(DEMEtics)\n'
            'dat <- read.'
            'D.Jost("')

    ok = call( [ 'Rscript', script_file] )

    if ok != 0:
        raise RuntimeError('Rscript for DEMEtics run unsuccessfully. Please inspect the error log.')

    raise RuntimeError('please inspect here: %s' % tmp_dir)