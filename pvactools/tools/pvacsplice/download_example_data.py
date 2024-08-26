import sys
import os
from urllib.request import urlretrieve
import tarfile
import gzip
import shutil

from pvactools.lib.download_example_data import DownloadExampleData

def define_parser():
    return DownloadExampleData.parser('pvacsplice')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    print("Downloading regtools results")
    DownloadExampleData(args.destination_directory , 'pvacsplice').execute()

    print("Downloading reference fasta - this step may take a few minutes")
    ref_fasta_url = "http://genomedata.org/pmbio-workshop/references/genome/all/ref_genome.tar"
    ref_fasta_dest = os.path.join(args.destination_directory, 'pvacsplice_example_data', 'ref_genome.tar')
    urlretrieve(ref_fasta_url, ref_fasta_dest)
    tar = tarfile.open(ref_fasta_dest)
    tar.extractall(os.path.join(args.destination_directory, 'pvacsplice_example_data', 'ref_genome'))
    tar.close()
    os.unlink(ref_fasta_dest)
    zipped_fasta_file = os.path.join(args.destination_directory, 'pvacsplice_example_data', 'ref_genome', 'ref_genome.fa.gz')
    unzipped_fasta_file = os.path.join(args.destination_directory, 'pvacsplice_example_data', 'ref_genome.fa')
    with gzip.open(zipped_fasta_file, 'rb') as f_in, open(unzipped_fasta_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    source_fai_file = os.path.join(args.destination_directory, 'pvacsplice_example_data', 'ref_genome', 'ref_genome.fa.fai')
    dest_fai_file = os.path.join(args.destination_directory, 'pvacsplice_example_data', 'ref_genome.fa.fai')
    shutil.copy(source_fai_file, dest_fai_file)
    shutil.rmtree(os.path.join(args.destination_directory, 'pvacsplice_example_data', 'ref_genome'))

    print("Downloading reference gtf - this step may take a few minutes")
    gtf_url = "https://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.chr.gtf.gz"
    gtf_dest = os.path.join(args.destination_directory, 'pvacsplice_example_data', 'Homo_sapiens.GRCh38.105.chr.gtf.gz')
    urlretrieve(gtf_url, gtf_dest)

if __name__ == '__main__':
    main()
