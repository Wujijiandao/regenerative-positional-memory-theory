#!/usr/bin/env python3
import argparse,hashlib,json
from pathlib import Path
ap=argparse.ArgumentParser(); ap.add_argument('--data-dir',default='data'); args=ap.parse_args()
root=Path(__file__).resolve().parents[1]; manifest=json.loads((root/'data_manifest/source_data.json').read_text())
ok=True
for name,m in manifest.items():
    p=root/m['target']
    if not p.exists(): print(f'MISSING  {p}'); ok=False; continue
    h=hashlib.sha256(p.read_bytes()).hexdigest(); good=h==m['sha256']; print(('OK      ' if good else 'MISMATCH'),p); ok &= good
raise SystemExit(0 if ok else 1)
