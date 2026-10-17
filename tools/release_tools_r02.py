#!/usr/bin/env python3
"""zip | rebuild | manifest.  zip: deterministic source-only zip (no build/, dist/, pdf, status, manifest).
rebuild: author rebuild from that zip in a temp dir, compare certificates (seconds stripped) to the reference.
manifest: written last."""
import hashlib, json, shutil, subprocess, sys, tempfile, zipfile  # noqa: E401
from pathlib import Path
R = Path(__file__).resolve().parents[1]
DIST = R / "dist"
EXCL_DIRS = {"build", "dist", "archive", "__pycache__", "reproducibility", ".pytest_cache"}


def files():
    out = []
    for p in sorted(R.rglob("*")):
        rel = p.relative_to(R)
        if p.is_file() and not (set(rel.parts) & EXCL_DIRS) and rel.parts[0] != "certificates" and \
                p.suffix not in (".pdf", ".aux", ".log", ".out") and rel.name not in ("RELEASE_STATUS.json", "RELEASE_MANIFEST_EXTERNAL.txt"):
            out.append(rel.as_posix())
    return out


def mkzip(path):
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for f in files():
            i = zipfile.ZipInfo("XVI-R/" + f, (1980, 1, 1, 0, 0, 0))
            i.external_attr = 0o644 << 16
            i.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(i, (R / f).read_bytes())


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if k != "seconds"}
    if isinstance(o, list):
        return [strip(v) for v in o]
    return o


def main(cmd):
    DIST.mkdir(exist_ok=True)
    if cmd == "zip":
        mkzip(DIST / "XVI-R-R02-source-only.zip")
        print(sha(DIST / "XVI-R-R02-source-only.zip"))
    elif cmd == "rebuild":
        z = DIST / "XVI-R-R02-source-only.zip"
        t = Path(tempfile.mkdtemp(prefix="xvir_rb_"))
        zipfile.ZipFile(z).extractall(t)
        r = subprocess.run([sys.executable, "tools/build_r02.py"], cwd=t / "XVI-R", capture_output=True, text=True)
        (R / "reproducibility").mkdir(exist_ok=True)
        (R / "reproducibility/author_rebuild_r02.log").write_text(r.stdout + r.stderr)
        ref = sorted((R / "certificates").glob("*.json"))
        match = [p.name for p in ref if (t / "XVI-R/certificates" / p.name).is_file() and
                 strip(json.loads(p.read_text())) == strip(json.loads((t / "XVI-R/certificates" / p.name).read_text()))]
        rep = {"zip_sha256": sha(z), "build_ok": "BUILD_OK" in r.stdout, "matched": len(match), "total": len(ref),
               "mismatched": sorted({p.name for p in ref} - set(match)),
               "verdict": "PASS" if "BUILD_OK" in r.stdout and len(match) == len(ref) else "FAIL"}
        (R / "reproducibility/author_rebuild_r02.json").write_text(json.dumps(rep, indent=1) + "\n")
        print(json.dumps(rep))
        shutil.rmtree(t, ignore_errors=True)
    elif cmd == "manifest":
        z = DIST / "XVI-R-R02.zip"
        with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED) as zz:
            for f in files() + sorted(p.relative_to(R).as_posix() for p in (R / "certificates").glob("*.json")) + \
                    ["manuscript/R02.pdf", "RELEASE_STATUS.json", "reproducibility/author_rebuild_r02.json"]:
                i = zipfile.ZipInfo("XVI-R/" + f, (1980, 1, 1, 0, 0, 0))
                i.external_attr = 0o644 << 16
                zz.writestr(i, (R / f).read_bytes(), zipfile.ZIP_DEFLATED)
        pages = subprocess.run(["pdfinfo", str(R / "manuscript/R02.pdf")], capture_output=True, text=True).stdout
        pages = [ln.split()[-1] for ln in pages.splitlines() if ln.startswith("Pages")][0]
        rb = json.loads((R / "reproducibility/author_rebuild_r02.json").read_text())
        txt = (f"XVI-R R02 external release manifest (written last)\n"
               f"XVI-R-R02-source-only.zip {sha(DIST / 'XVI-R-R02-source-only.zip')}\nXVI-R-R02.zip {sha(z)}\n"
               f"R02.pdf pages {pages} sha {sha(R / 'manuscript/R02.pdf')}\n"
               f"author rebuild {rb['verdict']} {rb['matched']}/{rb['total']}\n"
               f"imports XVI-Q Q11 source-only 1888f3277b793ee8cbf3954b3819bec03e41d8c4dc59757744dfda8b147a35ad\n"
               f"no /workspace sync performed\n")
        (R / "RELEASE_MANIFEST_EXTERNAL.txt").write_text(txt)
        print(txt)


if __name__ == "__main__":
    main(sys.argv[1])
