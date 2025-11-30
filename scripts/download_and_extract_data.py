import tarfile
import gzip
import os
import urllib.request
import io
import argparse

def download_and_extract_data(output_dir: str) -> None:
    '''
    Download and extract the GSE166992 dataset from GEO.

    Parameters
    ----------
    output_dir : str
        Directory to store extracted data
    '''
    GEO_URL = r'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166992/suppl/GSE166992_RAW.tar'

    with urllib.request.urlopen(GEO_URL) as response:
        file_like_object = io.BytesIO(response.read())
        with tarfile.open(fileobj=file_like_object, mode='r') as tar:
            tar.extractall(path=output_dir)

    for file in os.listdir(output_dir):
        path = os.path.join(output_dir, file)
        out_path = path[:-3]  # remove '.gz'
        with gzip.open(path, 'rb') as read:
            with open(out_path, 'wb') as f_out:
                f_out.writelines(read)

        os.remove(path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', type=str, default='data', help='Directory to store extracted data')
    args = parser.parse_args()

    download_and_extract_data(args.output_dir)