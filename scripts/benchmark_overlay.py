# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Attach frozen control summaries to a complete primary benchmark grid."""
import argparse
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"backend"))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("primary",help="Completed primary experiment ID")
    parser.add_argument("controls",nargs="+",help="Completed or stopped control experiment IDs")
    args=parser.parse_args()
    from app.services.benchmark import overlay_controls
    result=overlay_controls(args.primary,args.controls)
    print(f"Saved comparison {result['id']}: {result['completed_runs']} original runs; primary grid complete.")


if __name__=="__main__":main()
