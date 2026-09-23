#!/usr/bin/env python3
"""Byte-lock the slice: MANIFEST.sha256 over data+results+code+paper+tools+docs (excludes .git, work/, itself)."""
import hashlib, os
skip_dirs = {".git", "work"}
man = []
for root, dirs, files in os.walk("."):
    dirs[:] = [d for d in dirs if d not in skip_dirs]
    for fn in sorted(files):
        p = os.path.join(root, fn).lstrip("./")
        if p in ("MANIFEST.sha256",): continue
        h = hashlib.sha256(open(os.path.join(root, fn), "rb").read()).hexdigest()
        man.append(f"{h}  {p}")
man.sort(key=lambda l: l.split("  ", 1)[1])
open("MANIFEST.sha256", "w").write("\n".join(man) + "\n")
print(len(man), "files byte-locked")
