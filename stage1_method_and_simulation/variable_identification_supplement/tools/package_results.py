"""Assemble the variable identification publication files and detailed runtime archives."""

import csv
from contextlib import contextmanager
import hashlib
import os
from pathlib import Path
import shutil
import sys
import uuid
import zipfile


SUPPLEMENT = "variable_identification_supplement"
RESULTS_DIRECTORY = Path("github", "stage1_method_and_simulation", SUPPLEMENT, "results")
ASSET_DIRECTORY = Path("github_release_assets", "stage1_runtime", SUPPLEMENT)
ASSET_LIMIT = 1024 ** 3


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 ** 2), b""):
            digest.update(block)
    return digest.hexdigest()


def files_under(root):
    return sorted(path for path in root.rglob("*") if path.is_file())


def filesystem_path(path):
    path = path.resolve()
    if os.name != "nt" or str(path).startswith("\\\\?\\"):
        return path
    name = str(path)
    prefix = "\\\\?\\UNC\\" if name.startswith("\\\\") else "\\\\?\\"
    return Path(prefix + (name[2:] if name.startswith("\\\\") else name))


def file_row(path, relative):
    return {"file": relative, "bytes": path.stat().st_size, "sha256": sha256(path)}


def write_csv(path, rows, fields=("file", "bytes", "sha256")):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def batches(files):
    batch = []
    total = 0
    for path in files:
        size = path.stat().st_size
        if size > ASSET_LIMIT:
            raise ValueError("Runtime record exceeds archive allocation: " + path.name)
        if batch and total + size > ASSET_LIMIT:
            yield batch
            batch, total = [], 0
        batch.append(path)
        total += size
    if batch:
        yield batch


@contextmanager
def staging_directory(work):
    temporary = work / ("package_" + uuid.uuid4().hex)
    temporary.mkdir()
    try:
        yield temporary
    finally:
        previous_unrestored = (temporary / "previous_delivery").exists() and not (work / "delivery").exists()
        if temporary.resolve().parent == work.resolve() and temporary.name.startswith("package_") and not previous_unrestored:
            shutil.rmtree(str(temporary))


def package(work):
    work = filesystem_path(work)
    analysis, runtime = work / "analysis", work / "runtime"
    if not (analysis / "index.html").is_file() or not (runtime / "fit_and_test_status.csv").is_file():
        raise ValueError("Completed analysis files are required")
    delivery = work / "delivery"
    with staging_directory(work) as temporary:
        staging = temporary / "delivery"
        reports = staging / RESULTS_DIRECTORY
        assets = staging / ASSET_DIRECTORY
        reports.parent.mkdir(parents=True)
        assets.mkdir(parents=True)
        shutil.copytree(str(analysis), str(reports))
        write_csv(reports / "results_manifest.csv",
                  [file_row(path, path.relative_to(reports).as_posix()) for path in files_under(reports)])
        inputs = [path for path in files_under(runtime)
                  if path.parent == runtime or (path.parent == runtime / "checkpoints" and path.suffix == ".rds")]
        asset_rows, runtime_rows = [], []
        for index, batch in enumerate(batches(inputs), 1):
            asset = assets / ("stage1_variable_identification_{:03d}.zip".format(index))
            with zipfile.ZipFile(str(asset), "w", compression=zipfile.ZIP_STORED, allowZip64=True) as archive:
                for path in batch:
                    name = (Path(SUPPLEMENT, "runtime") / path.relative_to(runtime)).as_posix()
                    archive.write(str(path), name)
                    row = file_row(path, name)
                    row["asset"] = asset.name
                    runtime_rows.append(row)
            with zipfile.ZipFile(str(asset)) as archive:
                if archive.testzip() is not None:
                    raise IOError("Archive integrity check failed: " + asset.name)
            asset_rows.append(file_row(asset, asset.name))
        write_csv(assets / "runtime_manifest.csv", runtime_rows, ("file", "bytes", "sha256", "asset"))
        asset_rows.append(file_row(assets / "runtime_manifest.csv", "runtime_manifest.csv"))
        write_csv(assets / "release_assets.csv", asset_rows)

        previous = temporary / "previous_delivery"
        if delivery.exists():
            delivery.rename(previous)
        try:
            staging.rename(delivery)
        except OSError:
            if previous.exists():
                previous.rename(delivery)
            raise


if __name__ == "__main__":
    package(Path(sys.argv[1]).resolve())
