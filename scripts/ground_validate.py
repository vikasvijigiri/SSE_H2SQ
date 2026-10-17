# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Start or inspect the local random-parameter ground-state validation study."""
import argparse
import json
from pathlib import Path
import urllib.request


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--request',type=Path,help='Optional JSON request with lattice, seed, budgets and custom parameters')
    parser.add_argument('--id',help='Read an existing study instead of launching one')
    parser.add_argument('--api',default='http://127.0.0.1:8000/api/ground-validation')
    args=parser.parse_args()
    request=urllib.request.Request(f'{args.api}/{args.id}') if args.id else urllib.request.Request(args.api,data=args.request.read_bytes() if args.request else b'{}',headers={'Content-Type':'application/json'},method='POST')
    with urllib.request.urlopen(request,timeout=30) as response:job=json.load(response)
    print(json.dumps({k:job[k] for k in ['id','status','completed_runs','total_runs','verdict']},indent=2))


if __name__=='__main__':main()
