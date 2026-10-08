#!/usr/bin/env python3
"""MorphPot analysis results archive (Python 3.9+, standard library only).

Run from the MorphPot repository:
  python3 move_morphpot_test_results.py
  python3 move_morphpot_test_results.py --apply

Default: preview only. --apply copies, SHA256-verifies, then removes originals.
Destination: /Volumes/KIOXIA/MorphPot/test_results/<sample>/<kind>/<folder>
Unknown directories are left alone. Input-only pose/mesh directories stay.
--include NAME explicitly selects an otherwise unrecognized top-level folder.
--sample ID groups all selected folders under ID (use only for one individual).
Historical absolute paths in JSON/CSV are preserved, not rewritten.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import uuid
from datetime import datetime, timezone

PROTECTED = {'morphpot', 'tests', 'test', 'docs', 'scripts', 'src', 'build',
             'dist', 'venv', 'env', '__pycache__', 'input', 'inputs', 'data',
             'samples', 'deliverables', 'node_modules'}
MARKERS = {'rim_qa.json': 'rim_sections', 'whole_model.json': 'whole_model',
           'phase_experiment.json': 'phase_validation',
           'dimension_validation.json': 'dimension_validation',
           'thickness_diagnostics.json': 'thickness_diagnostics',
           'validation.json': 'validation'}


def safe_name(value):
    return re.sub(r'[^\w.\-]+', '_', value).strip('._') or 'unassigned'


def sample_from_name(name):
    if name.lower().startswith(('phase_', 'whole_', 'curvature_', 'validation')):
        return None
    name = re.sub(r'^Sample[-_]', '', name, flags=re.I)
    match = re.search(r'_(?:RadialSections|whole_model|rim_standardization|'
                      r'modefied|modified|phase_validation|validation)(?:_|$)', name, re.I)
    return safe_name(name[:match.start()]) if match else None


def classify(folder):
    markers = []
    input_samples = set()
    # Only inspect four levels; no recursive traversal of arbitrarily deep trees.
    for root, dirs, files in os.walk(folder, followlinks=False):
        relative = Path(root).relative_to(folder)
        dirs[:] = [d for d in dirs if not d.startswith('.') and
                   not (Path(root)/d).is_symlink()]
        if len(relative.parts) >= 3:
            dirs[:] = []
        for filename in files:
            if filename not in MARKERS and filename != 'metadata.json':
                continue
            path = Path(root)/filename
            if path.is_symlink():
                continue
            try:
                if path.stat().st_size > 5_000_000:
                    continue
                data = json.loads(path.read_text(encoding='utf-8'))
            except (OSError, ValueError):
                continue
            if not isinstance(data, dict):
                continue
            if filename in MARKERS:
                markers.append(MARKERS[filename])
            elif 'RadialSection' in str(data.get('program', '')):
                markers.append('rim_sections')
            for key in ('input', 'input_mesh', 'source_mesh'):
                value = data.get(key)
                if isinstance(value, str) and value.lower().endswith(('.ply', '.obj', '.stl')):
                    input_samples.add(safe_name(Path(value).stem))
    if not markers:
        return None, None
    # Outer run containers take priority over their nested QA/validation outputs.
    for kind in ('phase_validation', 'whole_model', 'rim_sections',
                 'dimension_validation', 'thickness_diagnostics', 'validation'):
        if kind in markers:
            break
    sample = sample_from_name(folder.name)
    if not sample and len(input_samples) == 1:
        sample = next(iter(input_samples))
    return kind, sample or ('multiple_samples' if len(input_samples) > 1 else 'unassigned')


def snapshot(folder):
    """Digest all regular files and retain directory topology; reject symlinks."""
    records = {}
    directories = []
    for root, dirs, files in os.walk(folder, followlinks=False):
        for name in dirs + files:
            path = Path(root)/name
            if path.is_symlink():
                raise ValueError(f'Symlink is not supported: {path}')
        directories.extend(str((Path(root)/d).relative_to(folder)) for d in dirs)
        for name in files:
            path = Path(root)/name
            if not path.is_file():
                raise ValueError(f'Not a regular file: {path}')
            digest = hashlib.sha256()
            with path.open('rb') as handle:
                for chunk in iter(lambda: handle.read(1024*1024), b''):
                    digest.update(chunk)
            records[str(path.relative_to(folder))] = [path.stat().st_size, digest.hexdigest()]
    return {'files': records, 'directories': sorted(directories)}


def tracked_paths(source):
    result = subprocess.run(['git', '-C', str(source), 'ls-files', '-z'],
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    return set(os.fsdecode(p).split('/')[0] for p in result.stdout.split(b'\0') if p) if result.returncode == 0 else set()


def save_log(path, log):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(log, ensure_ascii=False, indent=2), encoding='utf-8')
    temporary.replace(path)


def copy_contents(source, destination):
    """Copy bytes and topology only, without filesystem-specific attributes."""
    destination.mkdir()
    for root, dirs, files in os.walk(source, followlinks=False):
        target = destination/Path(root).relative_to(source)
        for name in dirs:
            if (Path(root)/name).is_symlink():
                raise ValueError('Source gained a symlink during copy')
            (target/name).mkdir()
        for name in files:
            original = Path(root)/name
            if original.is_symlink() or not original.is_file():
                raise ValueError(f'Unsupported file during copy: {original}')
            with original.open('rb') as incoming, (target/name).open('xb') as outgoing:
                shutil.copyfileobj(incoming, outgoing, length=1024*1024)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--source', type=Path, default=Path.cwd())
    parser.add_argument('--destination', type=Path, default=Path('/Volumes/KIOXIA/MorphPot'))
    parser.add_argument('--include', action='append', default=[], metavar='FOLDER', help='Explicit top-level output folder; repeatable')
    parser.add_argument('--sample', help='Override sample ID for every selected output folder')
    parser.add_argument('--allow-tracked-results', action='store_true', help='Allow tracked output folders (Git will see deletions)')
    parser.add_argument('--apply', action='store_true', help='Execute the previewed moves')
    args = parser.parse_args()
    source = args.source.expanduser().resolve()
    destination = args.destination.expanduser().resolve()
    if not source.is_dir():
        raise ValueError(f'Source directory does not exist: {source}')
    if destination == source or source in destination.parents or destination in source.parents:
        raise ValueError('Source and destination must be separate directory trees')
    includes = set(args.include)
    for name in includes:
        if Path(name).name != name or name.startswith('.') or name.lower() in PROTECTED:
            raise ValueError(f'Invalid/protected --include folder: {name}')
        if not (source/name).is_dir():
            raise ValueError(f'--include folder is missing: {name}')
    try:
        tracked = tracked_paths(source)
    except FileNotFoundError:
        tracked = set()
    plan, skipped = [], []
    for folder in sorted(source.iterdir()):
        if not folder.is_dir() or folder.name.startswith('.'):
            continue
        if folder.is_symlink() or folder.name.lower() in PROTECTED:
            skipped.append((folder.name, 'protected / symlink'))
            continue
        kind, sample = classify(folder)
        if folder.name in includes and not kind:
            kind, sample = 'other_results', sample_from_name(folder.name) or 'unassigned'
        if not kind:
            skipped.append((folder.name, 'no recognized output marker'))
            continue
        if folder.name in tracked and not args.allow_tracked_results:
            skipped.append((folder.name, 'contains Git-tracked files; explicit opt-in required'))
            continue
        target = destination/'test_results'/safe_name(args.sample or sample)/kind/folder.name
        if target.exists():
            raise ValueError(f'Destination already exists (no overwrite): {target}')
        plan.append({'source': str(folder), 'destination': str(target), 'status': 'planned'})
    print('MOVE PLAN' if args.apply else 'PREVIEW ONLY (add --apply to execute)')
    for item in plan:
        print(f"\n{item['source']}\n  -> {item['destination']}")
    print(f'\nSelected: {len(plan)} folders')
    for name, reason in skipped:
        print(f'KEEP: {name} ({reason})')
    if not args.apply or not plan:
        return
    # Avoid creating a fake external-drive directory on the system disk.
    if len(destination.parts) > 2 and destination.parts[1] == 'Volumes':
        volume = Path('/Volumes')/destination.parts[2]
        if not os.path.ismount(volume):
            raise ValueError(f'External volume is not mounted: {volume}')
    ancestor = destination
    while not ancestor.exists():
        ancestor = ancestor.parent
    total = sum(p.stat().st_size for item in plan for p in Path(item['source']).rglob('*') if p.is_file())
    if shutil.disk_usage(ancestor).free < total + 64*1024*1024:
        raise ValueError('Insufficient destination free space')
    destination.mkdir(parents=True, exist_ok=True)
    logdir = destination/'move_logs'
    logdir.mkdir(exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    logpath = logdir/f'move_{stamp}_{uuid.uuid4().hex[:8]}.json'
    log = {'created_utc': stamp, 'source': str(source), 'items': plan,
           'note': 'Existing absolute paths in result metadata have not been rewritten.'}
    save_log(logpath, log)
    for item in plan:
        original, target = Path(item['source']), Path(item['destination'])
        staging = target.parent/f'.{target.name}.copying-{uuid.uuid4().hex[:8]}'
        try:
            before = snapshot(original)
            target.parent.mkdir(parents=True, exist_ok=True)
            item['staging'] = str(staging)
            item['status'] = 'copying'
            save_log(logpath, log)
            copy_contents(original, staging)
            if snapshot(staging) != before or snapshot(original) != before:
                raise ValueError('SHA256 verification failed or source changed during copying')
            if target.exists():
                raise ValueError(f'Destination appeared during copying: {target}')
            staging.rename(target)
            item['status'] = 'copied_verified_source_retained'
            save_log(logpath, log)
            shutil.rmtree(original)
            item['status'] = 'moved'
            save_log(logpath, log)
            print(f"MOVED: {original.name}")
        except Exception as error:
            item['error'] = str(error)
            save_log(logpath, log)
            print(f'STOPPED. Check log: {logpath}', file=sys.stderr)
            raise
    print(f'Completed. Move log: {logpath}')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
