
#
#
# Multiplicity of Infection (MoI) calculation
# ===========================================
#
# histogram of sample with MoI:
#       MoI     Freq
#       1       10
#       2       15
#       3       8
#
# histogram of marker with MoI:
#       No of marker with allele > 1   Freq
#       1                               30
#       2                               5
#
#

def summary_moi(analytical_sets):

    moi = {}

    for analytical_set in analytical_sets:
        moi[analytical_set.label] = calculate_moi(analytical_set.allele_df)




def calculate_moi(allele_df):
    # given allele_df, return (moi_table)
    pass