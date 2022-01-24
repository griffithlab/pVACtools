import os
import sys
import pandas as pd
from pybiomart import Dataset
import argparse

class FilterRegtoolsResults():
    def __init__(self, input_file:str, sample_name:str, output_dir:str, score:int=10, distance:int=100, tsl:bool=True):
        self.input_file = input_file
        self.sample_name = sample_name
        self.output_dir = output_dir
        self.score = score
        self.distance = distance
        self.tsl = tsl

    def filter_regtools_results(self):
        # load in ensembl dataset
        dataset = Dataset(name='hsapiens_gene_ensembl', host='http://useast.ensembl.org')
        
        # read in regtools junctions.tsv
        junctions = pd.read_csv(self.input_file, sep='\t')

        # filter junctions by score, strand, anchor, remove NAs, reset index
        filter_junctions = junctions[(junctions['score'] > self.score) & (junctions['strand'] != '?') & (junctions['anchor'].isin(['D', 'A', 'NDA']))].dropna().reset_index()
        if len(filter_junctions) == 0:
            raise ValueError('No junctions passed score and anchor filters. Exiting...')
            
        # expand variant_info field into chrom, start, end
        filter_junctions[['variant_chrom', 'variant_start', 'variant_end']] = filter_junctions['variant_info'].str.split(':|-|,', expand=True)[[0, 1, 2]]

        # filter by chrom
        filter_junctions[filter_junctions['chrom'] == filter_junctions['variant_chrom']]

        # filter by distance: variant_start > start-distance and variant_start < end+distance 
        filter_junctions = filter_junctions[
            (filter_junctions['variant_start'].astype(int) > filter_junctions['start'].astype(int) - self.distance) & 
            (filter_junctions['variant_start'].astype(int) < filter_junctions['end'].astype(int) + self.distance)
            ] 

        # get transcript info
        tscript_dict = {i:{'gene_ids': x, 'transcripts': y} for i,(x,y) in enumerate(zip(filter_junctions['gene_ids'], filter_junctions['transcripts']))}

        # filter by transcript biotype
        pc_junctions = pd.DataFrame()
        for k,v in tscript_dict.items():
            protein_coding = dataset.query(attributes=[
                'ensembl_transcript_id', 
                'ensembl_gene_id', 
                'external_gene_name',
                'transcript_tsl',
                ], 
                filters={
                    'link_ensembl_transcript_stable_id': v['transcripts'].split(','),
                    'transcript_biotype': 'protein_coding',
                    'transcript_tsl': True,
                })
            if self.tsl:
                protein_coding = protein_coding[protein_coding['Transcript support level (TSL)'].str.contains(f'tsl1')]
            pc_junctions = pc_junctions.append(protein_coding)
            tscript_dict[k]['pc_transcripts'] = ','.join(protein_coding['Transcript stable ID'])
            tscript_dict[k]['pc_gene_ids'] = ','.join(protein_coding['Gene stable ID'].unique())
            tscript_dict[k]['pc_gene_names'] = ','.join(protein_coding['Gene name'].unique())

        # merge new values with filter_junctions
        tscript_df = pd.DataFrame(tscript_dict).transpose().replace('', float('NaN')).dropna().drop_duplicates()
        merged_junctions = filter_junctions.merge(tscript_df, on=['transcripts', 'gene_ids'])

        # filter merged_junctions to include fields for junction_to_fasta.py
        merged_junctions['strand'] = merged_junctions['strand'].replace(['+','-'], [1,-1])
        final_merged = merged_junctions[['name', 'pc_transcripts','chrom','start','end','anchor','strand','pc_gene_names']]
        
        # expand tscripts to one per line
        pd.options.mode.chained_assignment = None
        final_merged['pc_transcripts'] = final_merged.pc_transcripts.apply(lambda x: x.split(','))
        
        final_merged = final_merged.explode('pc_transcripts')

        final_merged.to_csv(os.path.join(self.output_dir, f'{self.sample_name}_filtered_junctions.tsv'), index=False, sep='\t')

        return final_merged


# write argparse
def define_parser():
        parser = argparse.ArgumentParser(
            description="Filter Regtools junctions",
        )
        parser.add_argument('input_file', help='input file path (Regtools output tsv)', type=str)
        parser.add_argument('output_dir', help='output folder path', type=str)
        parser.add_argument('-s', '--score', nargs='?', const=1, default=10, help='minimum junction coverage. Default is 10.', type=int)
        parser.add_argument('-d', '--distance', nargs='?', const=1, default=100, help='variant distance (in bps) from junction start/stop. Default is 100.', type=int)
        parser.add_argument('-t', '--tsl', nargs='?', const=1, help='option to filter transcripts via Transcript Support Level (TSL).', type=int, choices=[1,2,3,4,5])
        
        return parser


def main(args_input = sys.argv[1:]):
        parser = define_parser()
        args = parser.parse_args(args_input)

        filterRegtoolsResults(args.input_file, args.output_dir, args.score, args.distance, args.tsl).filter_regtools_results()

if __name__ == '__main__':
    main()
