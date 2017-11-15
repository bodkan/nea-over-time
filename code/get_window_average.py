import sys

# sys.stderr.write("EXAMPLE COMMAND (working on a regions file): bedmap --ec --delim '\\t' --echo --count --echo-map-score --echo-overlap-size windows.10000.bed genetic_map_GRCh37_all.txt.bed | python get_window_average.py 3\n")
# sys.stderr.write('EXAMPLE COMMAND (get recomb for regions + gaps between regions): awk \'BEGIN {OFS="\\t"} c == $1 {print $1, s, $2, "gap"; print $0, "orig"} c != $1 {c=$1; s=$3; print $0, "orig"}\' <regions-file> | bedmap --ec --delim \'\\t\' --echo --count --echo-map-score --echo-overlap-size windows.10000.bed genetic_map_GRCh37_all.txt.bed | python get_window_average.py 4\n')

ncols = int(sys.argv[1])

# chr1    0       10000   1       2.981822        10000
# chr1    10000   20000   1       2.981822        10000
# chr1    20000   30000   1       2.981822        10000
# chr1    30000   40000   1       2.981822        10000
# chr1    40000   50000   1       2.981822        10000
# chr1    50000   60000   2       2.981822;2.082414       5550;4450
# chr1    60000   70000   1       2.082414        10000
# chr1    70000   80000   1       2.082414        10000
# chr1    80000   90000   3       2.082414;2.081358;3.354927      2571;5598;1831
# chr1    90000   100000  1       3.354927        10000


for line in sys.stdin:

    line = line.rstrip('\n').split('\t')

    # print( line)
    items = int(line[ncols])
    caution_flag = '.'

    if items == 0:
        tot_len = 0
        weighted_mean = 0
        
    else:
        rates = [float(r) for r in line[ncols+1].split(';')]
        lengths = [int(l) for l in line[ncols+2].split(';')]
        tot_len = sum(lengths)
        weighted_mean = sum([rates[i] / tot_len * lengths[i] for i in range(items)])
        if sum(r == 0 for r in rates) > 0: caution_flag = 'POSSIBLE_LIFTOVER_ERROR'
        pass

    print('\t'.join(str(s) for s in line[:ncols] + [items, tot_len, weighted_mean, caution_flag]))


    
