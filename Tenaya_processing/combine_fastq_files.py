import os
import argparse
import subprocess
from subprocess import check_call, CalledProcessError
from pathlib import Path


def parse_arguments(parser=None):
    if not parser:
        parser = argparse.ArgumentParser()

    parser.add_argument("--input_dir", help="name of the file to tweet")
    parser.add_argument("--output_dir", help="Where to put the combined files")
    args = parser.parse_args()

    return args


DIR_TEMPLATE = "Ivey-MM-2338-{sample}_L00{lane}"
FILE_TEMPLATE = "Ivey-MM-2338-{sample}_S{sample_int}_L00{lane}_R1_001.fastq.gz"
OUTPUT_TEMPLATE = "{output_dir}output/S{sample}/{sample}.fastq.gz"
LANES = [str(x+1) for x in range(4)]
SAMPLES = [('{:02}'.format(x+1)) for x in range(18)]


class DefaultClass(object):
    def __init__(self, args):
        self.input_dir = args.input_dir
        self.output_dir = args.output_dir
        if not self.input_dir.endswith(os.path.sep):
            self.input_dir += os.path.sep
        if not self.output_dir.endswith(os.path.sep):
            self.output_dir += os.path.sep
        self.subdirectories = {}

        self.p = Path(self.input_dir)

    def get_subdirectories(self):
        directories = [x for x in self.p.iterdir() if x.is_dir()]
        for directory in directories:
            self.subdirectories[directory.parts[-1][:20]] = directory.parts[-1]

    def merge_samples(self):
        for sample in SAMPLES:
            cat_set = ['cat']
            for lane in LANES:
                directory = self.subdirectories[DIR_TEMPLATE.format(sample=sample, lane=lane)]
                filename = FILE_TEMPLATE.format(sample=sample, sample_int=int(sample), lane=lane)
                location = f"{self.input_dir}{directory}{os.path.sep}{filename}"
                cat_set.append(location)

            print(" ".join(cat_set))
            self.p = Path(f"{self.output_dir}{os.path.sep}output{os.path.sep}S{sample}")
            self.p.mkdir(parents=True, exist_ok=True)
            output = OUTPUT_TEMPLATE.format(output_dir=self.output_dir, sample=sample)
            try:
                with open(output, 'w') as f:
                    subprocess.check_call(cat_set, stdout=f)
            except CalledProcessError as cpe:
                print("Error while concattenating gzip files.", cpe)
                return

    def run(self):
        self.get_subdirectories()
        self.merge_samples()

    def close(self):
        pass


def main():
    args = parse_arguments()
    default_class = DefaultClass(args)
    default_class.run()


if __name__ == '__main__':
    main()
