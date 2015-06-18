
from collections import defaultdict

def summarize_alleles( allele_df ):
    """ return a tuple of (dict, dict):
        1dict: alleles: [ (allele, freq, count, mean_height, min_size, max_size, delta), ...]
        2dict: marker: ( [ size, ...], [ height, ....] )
    """

    allele_list = defaultdict(list)
    allele_plot = defaultdict(lambda x = None: ([], []))
    grouped = allele_df.groupby( ['marker_id', 'value'] )

    for (marker_id, allele), df in grouped:

        allele_list[marker_id].append(
            (allele, len(df), np.mean( df['height'] ), min(df['size']), max(df['size']),
                list(df['sample_id']), np.mean( df['size'] ))
        )

        code = Marker.get(marker_id).code
        allele_plot[code][0].extend( df['size'] )
        allele_plot[code][1].extend( df['height'] )
            

    # calculate other stuff

    results = {}

    for marker_id in allele_list:
        alleles = allele_list[marker_id]
        total_allele = sum( x[1] for x in alleles )
        allele_params = [
            (allele, count/total_allele, count, mean_height, min_size, max_size,
                max_size - min_size, sample_ids, mean_size )
            for (allele, count, mean_height, min_size, max_size, sample_ids, mean_size )
            in alleles
        ]

        delta_status = check_delta( allele_params)

        results[marker_id] = dict(
            code = Marker.get(marker_id).code,
            unique_allele = len(allele_params),
            total_allele = total_allele,
            alleles = allele_params,
            delta_status = delta_status )

    return (results, allele_plot)

