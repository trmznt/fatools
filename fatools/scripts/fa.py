
def init_argparser(parser=None):

    from fatools.lib.fautil import cmds

    return cmds.init_argparser(parser)



def main(args):

    from fatools.lib.fautil import cmds
    return cmds.main(args)