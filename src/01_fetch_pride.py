"""Extract the PD protein table and study metadata from PXD072052.

The search-engine output on PRIDE is a single 223 GB ZIP64 archive, but the two
members we need total under 3 MB uncompressed. The PRIDE HTTPS mirror advertises
Accept-Ranges, so we read the ZIP64 central directory from the tail and inflate
only those members.
"""

import re
import struct
import sys
import xml.etree.ElementTree as ET
import zlib

import pandas as pd
import requests

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import DATA_RAW  # noqa: E402

ZIP_URL = (
    "https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD072052/"
    "Remeasured_TRY_0463_500.zip"
)
WANTED = {"_Proteins.txt": "Proteins.txt", ".pdStudy": "study.pdStudy"}


class RemoteZip:
    def __init__(self, url):
        self.url = url
        self.session = requests.Session()
        head = self.session.head(url, timeout=120, allow_redirects=True)
        head.raise_for_status()
        if head.headers.get("Accept-Ranges") != "bytes":
            raise RuntimeError("server does not support range requests")
        self.size = int(head.headers["Content-Length"])
        self.entries = self._central_directory()

    def _range(self, start, end):
        r = self.session.get(
            self.url, headers={"Range": f"bytes={start}-{end}"}, timeout=300
        )
        r.raise_for_status()
        return r.content

    def _central_directory(self):
        tail = self._range(max(0, self.size - 66000), self.size - 1)
        loc = tail.rfind(b"PK\x06\x07")
        if loc >= 0:
            eocd64_off = struct.unpack("<Q", tail[loc + 8 : loc + 16])[0]
            eocd64 = self._range(eocd64_off, eocd64_off + 55)
            cd_size, cd_off = struct.unpack("<QQ", eocd64[40:56])
        else:
            eocd = tail.rfind(b"PK\x05\x06")
            cd_size, cd_off = struct.unpack("<II", tail[eocd + 12 : eocd + 20])
        cd = self._range(cd_off, cd_off + cd_size - 1)

        entries, p = {}, 0
        while p < len(cd) and cd[p : p + 4] == b"PK\x01\x02":
            (*_, csize, usize, nlen, elen, clen, _, _, _, offset) = struct.unpack(
                "<HHHHHHIIIHHHHHII", cd[p + 4 : p + 46]
            )
            name = cd[p + 46 : p + 46 + nlen].decode("utf-8", "replace")
            extra = cd[p + 46 + nlen : p + 46 + nlen + elen]
            q = 0
            while q + 4 <= len(extra):
                hid, hsz = struct.unpack("<HH", extra[q : q + 4])
                if hid == 1:
                    body, k = extra[q + 4 : q + 4 + hsz], 0
                    if usize == 0xFFFFFFFF:
                        usize = struct.unpack("<Q", body[k : k + 8])[0]
                        k += 8
                    if csize == 0xFFFFFFFF:
                        csize = struct.unpack("<Q", body[k : k + 8])[0]
                        k += 8
                    if offset == 0xFFFFFFFF:
                        offset = struct.unpack("<Q", body[k : k + 8])[0]
                q += 4 + hsz
            entries[name] = {"offset": offset, "csize": csize, "usize": usize}
            p += 46 + nlen + elen + clen
        return entries

    def read(self, name):
        e = self.entries[name]
        header = self._range(e["offset"], e["offset"] + 29)
        nlen, elen = struct.unpack("<HH", header[26:30])
        start = e["offset"] + 30 + nlen + elen
        blob = self._range(start, start + e["csize"] - 1)
        return zlib.decompressobj(-15).decompress(blob)


def parse_study(xml_bytes):
    """Sample sheet from the pdStudy XML, ordered by SampleNumber -> F-number."""
    root = ET.fromstring(xml_bytes.decode("utf-8-sig"))
    rows = []
    for sample in root.iter("Sample"):
        num = int(sample.get("SampleNumber"))
        name = sample.get("Name")
        batch = re.search(r"TRY_(\d+)", name).group(1)
        prep = re.search(r"_(P\d+)", name).group(1)
        fraction = "EV" if "_EV_" in name else "SOL"
        rows.append(
            {
                "f_number": f"F{num}",
                "sample_number": num,
                "run_name": name,
                "batch": f"TRY_{batch}",
                "prep": prep,
                "fraction": fraction,
            }
        )
    return pd.DataFrame(rows).sort_values("sample_number").reset_index(drop=True)


def main():
    targets = {out: DATA_RAW / out for out in WANTED.values()}
    if all(p.exists() for p in targets.values()):
        print("cached, skipping remote fetch")
    else:
        rz = RemoteZip(ZIP_URL)
        print(f"archive {rz.size:,} bytes, {len(rz.entries)} members")
        for suffix, out in WANTED.items():
            name = next(n for n in rz.entries if n.endswith(suffix))
            data = rz.read(name)
            (DATA_RAW / out).write_bytes(data)
            print(f"  {out:16} <- {name.split('/')[-1]}  ({len(data):,} B)")

    samples = parse_study((DATA_RAW / "study.pdStudy").read_bytes())
    samples.to_csv(DATA_RAW / "samples.tsv", sep="\t", index=False)
    print(f"\n{len(samples)} runs: "
          f"{(samples.fraction == 'EV').sum()} EV / {(samples.fraction == 'SOL').sum()} SOL, "
          f"{samples.prep.nunique()} preps, {samples.batch.nunique()} batches")
    print(samples.to_string(index=False))


if __name__ == "__main__":
    main()
