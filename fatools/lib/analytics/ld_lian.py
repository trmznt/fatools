
from fatools.lib.utils import cout, cerr
from fatools.lib.analytics.export import export_flat
from subprocess import Popen, PIPE

class LianResult(object):

    def __init__( self, output, error, n ):
        self.output = output.decode('ASCII')
        self.error = error
        self.ld = ''
        self.pval = ''
        self.n = n
        self.parse()

    def __len__(self):
        return self.n

    def get_LD(self):
        return self.ld

    def get_pvalue(self):
        return self.pval

    def get_output(self):
        return self.output

    def parse(self):
        for line in self.output.split('\n'):
            if line.startswith('St. IA'):
                self.ld = line[6:].strip()
            elif line.startswith('P'):
                self.pval = line[2:].strip()


def run_lian( analytical_sets, dbh ):

    results = []

    for analytical_set in analytical_sets:

        data_set = analytical_set.allele_df.mlgt

        if len(data_set) <= 2:
            r = LianResult( output = b'', error = b'', n = len(data_set) )
            r.ld = '-'
            r.pval = '-'
            results.append( (analytical_set.label, r) )
            continue

        p = Popen(["lian"],
                stdin=PIPE, stdout=PIPE, stderr=PIPE,
                close_fds=True)
        export_flat( analytical_set, dbh, p.stdin )
    #p.stdin.write( export_flat( analytical_set ).read() )
        p.stdin.close()

        result = LianResult( output = p.stdout.read(), error = p.stderr.read(), n = len(data_set) )
        results.append( (analytical_set.label, result) )

    return results
