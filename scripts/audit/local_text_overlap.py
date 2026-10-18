# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""A narrow, reproducible exact-shingle screen against supplied local drafts.

This is not a plagiarism detector. It does not search publisher databases,
books, theses, or the full web and cannot assess semantic paraphrase.
"""
import argparse
from pathlib import Path
import re


def prose(path, section_start=None, section_end=None):
    text=Path(path).read_text(encoding="utf-8",errors="ignore")
    if section_start:
        start=text.index(section_start)
        end=text.index(section_end,start) if section_end else len(text)
        text=text[start:end]
    text=re.sub(r"(?<!\\)%.*"," ",text)
    text=re.sub(r"\\begin\{(?:align\*?|equation\*?|eqnarray\*?)\}.*?\\end\{(?:align\*?|equation\*?|eqnarray\*?)\}"," ",text,flags=re.S)
    text=re.sub(r"\$.*?\$"," ",text,flags=re.S)
    text=re.sub(r"\\[A-Za-z]+\*?(?:\[[^]]*\])?"," ",text)
    text=re.sub(r"[{}\\]"," ",text)
    return re.findall(r"[a-z0-9]+",text.lower())


def shingles(tokens,n):
    return {tuple(tokens[i:i+n]) for i in range(max(0,len(tokens)-n+1))}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("candidate",type=Path);p.add_argument("baseline",type=Path)
    p.add_argument("--start",default=None);p.add_argument("--end",default=None)
    p.add_argument("--words",type=int,default=8)
    args=p.parse_args()
    candidate=prose(args.candidate,args.start,args.end)
    baseline=prose(args.baseline)
    source=shingles(candidate,args.words); reference=shingles(baseline,args.words)
    matches=source & reference
    percent=100*len(matches)/len(source) if source else 0.
    print(f"candidate_tokens={len(candidate)}")
    print(f"unique_{args.words}-word_shingles={len(source)}")
    print(f"exact_local_matches={len(matches)}")
    print(f"matched_shingle_fraction={percent:.4f}%")
    print(f"baseline={args.baseline}")
    print("scope=exact normalized word sequences in this local baseline only; no semantic or external-corpus search")
    for phrase in sorted(matches): print("MATCH:"," ".join(phrase))


if __name__=="__main__": main()
