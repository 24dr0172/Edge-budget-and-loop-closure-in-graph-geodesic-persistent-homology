"""Copy a supporting note between explicitly supplied external paths."""
from pathlib import Path
import argparse,shutil
ROOT=Path(__file__).resolve().parents[1]
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();source=args.source.resolve();output=args.output.resolve()
    if any(p==ROOT or ROOT in p.parents for p in (source,output)):
        parser.error('Keep source notes and their copies outside this numerical release.')
    if output.exists():parser.error('Choose a new output path.')
    output.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,output)
if __name__=='__main__':main()
