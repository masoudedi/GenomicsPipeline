import argparse
from pathlib import Path
from pack.cardinal import *
from pack.gath import GatherMetrics

parser = argparse.ArgumentParser(prog='ILYOME', description='Ilyome pipeline for whole-exome sequencing.')
parser.add_argument("--version", action="version", version="Ilyome v.3.3")
parser.add_argument('-i', metavar="--input", required=True, type=Path, help="Path to the input files to be analyzed")
parser.add_argument("-o", metavar="--output", required=True, type=Path, help="Path to the output files to be saved")
parser.add_argument('-r', metavar="--reference", required=True,  type=Path, help="Path to the reference file")
parser.add_argument('-l', metavar="--interval", required=True,  type=Path, help="Path to the interval file")
parser.add_argument('-f', metavar="--forks", default=20, type=int, help="Number of forks")
parser.add_argument('-padding', default=300, help="padding for the interval")
args = parser.parse_args()

Mapper(args.i, args.o, args.r, forks=args.f)
Sorter(args.i, args.o, forks=args.f)
GvcfCaller(args.i, args.o, args.r, args.l, args.padding, forks=args.f)
GetVCF(args.i, args.o, args.r, forks=args.f)
VCFmetrics(args.i, args.o, args.dbsnp, forks=args.f)
Depth(args.i, args.o, args.r, args.l, forks=args.f)
Coverage(args.i, args.o, args.r, args.l, args.l, forks=args.f)
GatherMetrics(args.i)
