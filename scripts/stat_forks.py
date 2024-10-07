#!/usr/bin/env python3

import sys
import argparse
from signal import signal, SIGPIPE, SIG_DFL
import statistics
import pandas

signal(SIGPIPE, SIG_DFL) 

parser = argparse.ArgumentParser(description='Add BrdU statistics of bed file of features. Input file must have header')
parser.add_argument('--bedgraph', '-b', help='Input bedgraph file with BrdU probabilities (or any other suitable metric) [%(default)s]', default='-')
parser.add_argument('--fork-bed', '-f', help='Input bed file to annotate [%(default)s]', default='-')
parser.add_argument('--brdu-cutoff', '-c', help='Probability cutoff to classify BrdU [%(default)s]', default=0.5, type=float)
parser.add_argument('--version', '-v', action='version', version='%(prog)s 0.1.0')

class Feature:
    def __init__(self, feature_dict):
        self.position = None
        self.prob_brdu = [] # List of BrdU probs for this feature
        self.stats = {}
        self.feature_dict = feature_dict # Row from bed file converted to dict

    def summarise_brdu(self, brdu_cutoff):
        # self.stats['median_brdu'] = statistics.median(self.prob_brdu)
        self.stats['n_brdu'] = sum([x > brdu_cutoff for x in self.prob_brdu])
        self.stats['n_sites'] = len(self.prob_brdu)
        self.stats['pct_brdu_sites'] = self.stats['n_brdu'] / self.stats['n_sites']
    
    def print_line(self, header):
        # header: Output fields. Typically the same columns you found in input
        # plus the BrdU statistics
        outline = []
        for x in header:
            if x in self.feature_dict and x in self.stats:
                raise Exception('Duplicated column %s' % x)
            elif x in self.feature_dict:
                outline.append(self.feature_dict[x])
            elif x in self.stats:
                v = self.stats[x]
                if x == 'pct_brdu_sites':
                    v = '%.4f' % v
                outline.append(v)
            else:
                raise Exception('Column %s not found' % x)
        return '\t'.join([str(x) for x in outline])


def prepare_bed(bedfile, required_columns=['#chrom', 'start', 'end', 'rname']):
    # Read the bed file of features
    # Output a dict with the original header and a dictionary of reads.
    # 'reads' has key the read name and value a list of Features found on this
    # read
    if bedfile == '-':
        bed = pandas.read_csv(sys.stdin, sep='\t')
    else:
        bed = pandas.read_csv(bedfile, sep='\t')
    for x in required_columns:
        if x not in bed.columns:
            sys.stderr.write('Column "%s" not in bed file\n' % x)
            sys.exit(1)
    if len(bed) > 0:
        assert bed.start.dtype == 'int64'
        assert bed.end.dtype == 'int64'

    header = bed.columns
    reads = {}
    for idx,row in bed.iterrows():
        if row.rname not in reads:
            reads[row.rname] = []
        feature = Feature(row.to_dict())
        feature.position = idx
        reads[row.rname].append(feature)
    return {'header': header, 'reads': reads}

def intersect_bedgraph(feature_reads, bedgraph):
    # Populate Feature objects in feature_reads by intersecting them with
    # bedgraph of P(BrdU)
    
    #if len(feature_reads) == 0:
    #    return
        # If there are no reads in feature_reads, you could return immediately
        # end exit the program. However, if you read the bedgraph from stdin,
        # you get a pipefail and snakemake complains. So we just read the
        # bedgraph until the end. 

    if bedgraph == '-':
        bdg = sys.stdin
    else:
        bdg = open(bedgraph)

    rname_idx = None
    for line in bdg:
        line = line.strip().split('\t')
        if rname_idx is None:
            # Read bedgraph header and get the position of the required columns
            chrom_idx = line.index('#chrom')
            start_idx = line.index('start')
            prob_brdu_idx = line.index('prob_brdu')
            rname_idx = line.index('rname')
            continue
        
        rname = line[rname_idx]
        if rname in feature_reads:
            # If this read contains a feature, i.e. it is present in the bed
            # file of features, process it to extract P(BrdU) provided the
            # position of BrdU intersect a feature
            brdu_pos = int(line[start_idx])
            features = feature_reads[rname] # Features found on this read
            for feature in features:
                # If any feature intersects this position update the Feature object
                if feature.feature_dict['#chrom'] == line[chrom_idx] and \
                        brdu_pos >= feature.feature_dict['start'] and \
                        brdu_pos < feature.feature_dict['end']:
                    prob_brdu = float(line[prob_brdu_idx])
                    feature.prob_brdu.append(prob_brdu)


if __name__ == '__main__':
    args = parser.parse_args()
    
    if args.brdu_cutoff > 1 or args.brdu_cutoff < 0:
        sys.stderr.write('Invalid value for --brdu-cutoff/-c "%s" it must be between 0 and 1\n' % args.brdu_cutoff)
        sys.exit(1)

    feature_reads = prepare_bed(args.fork_bed)
    header = list(feature_reads['header'])
    feature_reads = feature_reads['reads']
    
    intersect_bedgraph(feature_reads, args.bedgraph)

    # Reshape the feature_reads dictionary from {rname:[list of features on
    # this read]} to {position:feature}
    position_dict = {} # Key is the index of the original bed file, value is the feature
    for rname in feature_reads:
        for f in feature_reads[rname]:
            f.summarise_brdu(args.brdu_cutoff)
            position_dict[f.position] = f
    # Output in the same order you saw in input
    output_order = sorted(position_dict.keys())

    header.extend(['n_sites', 'pct_brdu_sites'])
    print('\t'.join(header))

    for i in output_order:
        line = position_dict[i].print_line(header)
        print(line)
