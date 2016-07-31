from fatools.lib.analytics.export import export_demetics
from fatools.lib.utils import cerr, cout, random_string
from subprocess import call
from collections import defaultdict
import numpy as np
import datetime


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

    today = datetime.date.today().strftime('%Y-%m-%d')

    # TODO: prepare stdout & stderr
    ok = call( [ 'Rscript', script_file] )

    if ok != 0:
        raise RuntimeError('Rscript for DEMEtics run unsuccessfully. Please inspect the error log.')

    # demetics save its output to files with names AND dates...
    #mean_file = "%s/dat.pairwise.Dest.mean.%s.txt" % (tmp_dir, today)
    ci_file = "%s/dat.pairwise.Dest.mean.ci.%s.txt" % (tmp_dir, today)

    d = defaultdict(dict)
    with open(ci_file) as infile:
        in_data = False
        for r in infile:
            r = r.strip()
            if not in_data:
                if r == 'Dest.mean Population1 Population2 Lower.0.95.CI Upper.0.95.CI':
                    in_data = True
                continue

            cols = r.split()
            d[cols[1]][cols[2]] = float(cols[0])
            d[cols[2]][cols[1]] = ( float(cols[3]), float(cols[4]) )

    raise RuntimeError('please inspect here: %s' % tmp_dir)

