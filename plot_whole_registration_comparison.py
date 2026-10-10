#!/usr/bin/env python3
"""Explain whole-shape registration using saved exploration outputs, no mesh rerun.

Requires pottery_whole_validation.py, NumPy, SciPy and Matplotlib. Raw and
registered distances use IDENTICAL raw-arc stations and the SAME saved model.
Registration is replayed from recorded accepted warp fields; rejected directions
remain unchanged. This is descriptive registration fit, not accuracy validation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from pottery_whole_validation import csv_write, nearest_segments, resample

VERSION = '0.1.0'
METHOD = 'neck_segmented_warp'
FRAMES = ('original', 'registered')
COLORS = {'original': '#385aab', 'registered': '#ce4d36'}


def replay(points, record):
    """Replay the exact recorded shared field without changing sampling stations."""
    points = np.asarray(points, dtype=float)
    if record['status'] == 'rejected':
        return points.copy()
    if record['status'] != 'applied':
        raise ValueError('Unsupported warp status: ' + record['status'])
    knots = np.asarray(record['knots_z_mm'], dtype=float)
    delta = np.asarray(record['delta_rz_mm'], dtype=float)
    if knots.shape != (3,) or delta.shape != (3, 2) or np.any(np.diff(knots) <= 0):
        raise ValueError('Invalid recorded warp field')
    result = points.copy()
    active = points[:, 1] > knots[0]
    z = np.clip(points[active, 1], knots[0], knots[-1])
    result[active] += np.column_stack([
        PchipInterpolator(knots, delta[:, j], extrapolate=False)(z) for j in (0, 1)
    ])
    return result


def describe(values):
    values = np.asarray(values, dtype=float)
    if not len(values) or not np.isfinite(values).all():
        raise ValueError('Distance sample is empty or non-finite')
    return dict(n_stations=len(values), mean_mm=float(values.mean()),
                std_mm=float(values.std(ddof=0)), min_mm=float(values.min()),
                median_mm=float(np.median(values)), rms_mm=float(np.sqrt(np.mean(values**2))),
                p95_mm=float(np.percentile(values, 95)), max_mm=float(values.max()))


def load_group(root, metadata, group):
    selected = metadata[group + '_angles_deg']
    records = {float(r['angle_deg']): r for r in metadata['warp_' + group]}
    landmark_file = 'landmarks.json' if group == 'training' else 'holdout_landmarks.json'
    landmarks = {float(r['angle_deg']): r for r in json.loads(
        (root / 'analysis' / landmark_file).read_text(encoding='utf-8'))}
    profiles = []
    with np.load(root / 'analysis' / f'raw_{group}_profiles.npz', allow_pickle=False) as data:
        angles = np.asarray(data['angles_deg'], dtype=float)
        for angle in selected:
            matches = np.flatnonzero(np.isclose(angles, angle, atol=1e-8, rtol=0))
            if len(matches) != 1 or float(angle) not in records or float(angle) not in landmarks:
                raise ValueError(f'Missing or ambiguous profile/warp/landmark: {group} {angle}')
            k = int(matches[0])
            profiles.append(dict(angle=float(angle), record=records[float(angle)],
                                 landmarks=landmarks[float(angle)],
                                 **{s: np.array(data[f'p{k:03d}_{s}_mm']) for s in ('outer', 'inner')}))
    if not profiles:
        raise ValueError('No selected profiles: ' + group)
    return profiles


def measure(profiles, models, interval_mm):
    rows = []
    for profile in profiles:
        lm = profile['landmarks']
        guard = max(lm['base_outer_mm'][1], lm['base_inner_mm'][1]) + .1 * (
            lm['lip_mm'][1] - lm['base_outer_mm'][1])
        for surface in ('outer', 'inner'):
            raw = profile[surface]
            s, points = resample(raw, interval_mm)
            length = np.linalg.norm(np.diff(raw, axis=0), axis=1).sum()
            base_arc = s if surface == 'outer' else length - s
            regions = np.where(points[:, 1] <= guard, 'base',
                               np.where(base_arc >= lm[surface]['arc_mm'], 'neck_to_lip', 'body'))
            registered = replay(points, profile['record'])
            before, _ = nearest_segments(points, models[surface])
            after, _ = nearest_segments(registered, models[surface])
            for k in range(len(points)):
                rows.append(dict(angle_deg=profile['angle'], surface=surface,
                    warp_status=profile['record']['status'], raw_arc_mm=float(s[k]),
                    region=str(regions[k]), raw_r_mm=float(points[k, 0]), raw_z_mm=float(points[k, 1]),
                    registered_r_mm=float(registered[k, 0]), registered_z_mm=float(registered[k, 1]),
                    original_distance_mm=float(before[k]), registered_distance_mm=float(after[k])))
    return rows


def summarise(rows, group):
    summary = []
    for frame in FRAMES:
        for surface in ('outer', 'inner', 'combined'):
            for region in ('all', 'base', 'body', 'neck_to_lip'):
                subset = [r for r in rows if (surface == 'combined' or r['surface'] == surface)
                          and (region == 'all' or r['region'] == region)]
                if not subset:
                    continue
                summary.append(dict(dataset=group, frame=frame, surface=surface, region=region,
                    n_profiles=len({r['angle_deg'] for r in subset}),
                    **describe([r[frame + '_distance_mm'] for r in subset])))
    return summary


def overlay(profiles, models, target, path, title, detail=False):
    fig, axes = plt.subplots(1, 2, figsize=(12, 8), sharex=True, sharey=True)
    all_points = []
    for ax, frame in zip(axes, FRAMES):
        for p in profiles:
            color = '#397bb8' if p['record']['status'] == 'applied' else '#b98836'
            for side in ('outer', 'inner'):
                q = p[side] if frame == 'original' else replay(p[side], p['record'])
                all_points.append(q)
                ax.plot(*q.T, color=color, lw=.45, alpha=.25)
        for q in models.values():
            ax.plot(*q.T, color='#ae1838', lw=1.6)
        for label in ('neck', 'lip'):
            r, z = target[label]
            ax.scatter(r, z, marker='*', color='black', s=40, zorder=5)
        ax.plot(0, 0, marker='+', color='black', ms=8)
        ax.set(aspect='equal', xlabel='Radius r [mm]', ylabel='Input height z [mm]')
        ax.grid(alpha=.15)
    applied = sum(p['record']['status'] == 'applied' for p in profiles)
    axes[0].set_title('Before: original r,z; no landmark alignment')
    axes[1].set_title(f'After: constrained shared warp\nApplied {applied}/{len(profiles)}; '
                      f'rejected {len(profiles)-applied} kept unchanged')
    q = np.vstack(all_points + list(models.values()))
    if detail:
        low = float(target['neck'][1]) - 10
        high = max(float(target['lip'][1]) + 5, q[:, 1].max() + 2)
        upper = q[q[:, 1] >= low]
        xmin, xmax = upper[:, 0].min() - 3, upper[:, 0].max() + 3
    else:
        xmin, xmax = min(0., q[:, 0].min()) - 3, q[:, 0].max() + 3
        low, high = min(0., q[:, 1].min()) - 3, q[:, 1].max() + 3
    axes[0].set_xlim(xmin, xmax)
    axes[0].set_ylim(low, high)
    legend = [Line2D([0], [0], color='#397bb8', label='Source: warp accepted'),
              Line2D([0], [0], color='#b98836', label='Source: warp rejected'),
              Line2D([0], [0], color='#ae1838', label='Same representative model in both panels'),
              Line2D([0], [0], marker='*', color='black', ls='', label='Target neck / lip')]
    fig.legend(handles=legend, loc='lower center', ncol=2, fontsize=8)
    fig.suptitle(title + ' | same profiles, axes and model')
    fig.tight_layout(rect=(0, .075, 1, .95))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def statistical_plot(rows, summary, group, bin_mm, path, title):
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    metrics = ['mean_mm', 'std_mm', 'rms_mm', 'p95_mm', 'max_mm']
    x = np.arange(len(metrics))
    for j, frame in enumerate(FRAMES):
        row = next(r for r in summary if r['frame'] == frame and r['surface'] == 'combined'
                   and r['region'] == 'all')
        axes[0, 0].bar(x + (j-.5)*.36, [row[k] for k in metrics], width=.36,
                       color=COLORS[frame], label=frame)
    axes[0, 0].set(xticks=x, xticklabels=['Mean', 'SD', 'RMS', 'P95', 'Max'],
                   ylabel='Distance [mm]', title='Pooled station statistics (both walls)')
    angle_rows = []
    for angle in sorted({r['angle_deg'] for r in rows}):
        subset = [r for r in rows if r['angle_deg'] == angle]
        for frame in FRAMES:
            angle_rows.append(dict(dataset=group, angle_deg=angle, frame=frame,
                warp_status=subset[0]['warp_status'],
                **describe([r[frame + '_distance_mm'] for r in subset])))
    for frame in FRAMES:
        selected = [r for r in angle_rows if r['frame'] == frame]
        axes[0, 1].plot([r['angle_deg'] for r in selected], [r['rms_mm'] for r in selected],
                        '.-', color=COLORS[frame], lw=.9, ms=3, label=frame)
    axes[0, 1].set(xlabel='Half-section azimuth [deg]', ylabel='RMS [mm]', title='Same sections by azimuth')
    bins = np.floor(np.array([r['raw_z_mm'] for r in rows]) / bin_mm).astype(int)
    height_rows = []
    for b in np.unique(bins):
        indices = np.flatnonzero(bins == b)
        for frame in FRAMES:
            height_rows.append(dict(dataset=group, frame=frame, raw_z_low_mm=float(b*bin_mm),
                raw_z_high_mm=float((b+1)*bin_mm),
                **describe([rows[k][frame + '_distance_mm'] for k in indices])))
    for frame in FRAMES:
        selected = [r for r in height_rows if r['frame'] == frame]
        z = np.array([(r['raw_z_low_mm'] + r['raw_z_high_mm'])/2 for r in selected])
        mean = np.array([r['mean_mm'] for r in selected])
        sd = np.array([r['std_mm'] for r in selected])
        axes[1, 0].plot(z, mean, color=COLORS[frame], label=frame)
        axes[1, 0].fill_between(z, np.maximum(0, mean-sd), mean+sd, color=COLORS[frame], alpha=.12)
        d = np.sort([r[frame + '_distance_mm'] for r in rows])
        axes[1, 1].plot(d, np.arange(1, len(d)+1)/len(d), color=COLORS[frame], label=frame)
    axes[1, 0].set(xlabel='RAW input height z [mm]', ylabel='Mean distance +/- SD [mm]',
                   title=f'Fixed raw-height bins ({bin_mm:g} mm); SD is not a CI')
    axes[1, 1].set(xlabel='Distance [mm]', ylabel='Cumulative proportion', title='Distance distribution')
    for ax in axes.flat:
        ax.grid(axis='y', alpha=.2)
        ax.legend(fontsize=8)
    fig.suptitle(title + ' | same raw-arc stations; descriptive fit, not accuracy')
    fig.tight_layout(rect=(0, 0, 1, .95))
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return angle_rows, height_rows


def run(root, output=None, interval_mm=1., height_bin_mm=5.):
    root = Path(root)
    output = Path(output) if output is not None else root / 'validation' / 'registration_comparison'
    if not np.isfinite(interval_mm) or interval_mm <= 0 or not np.isfinite(height_bin_mm) or height_bin_mm <= 0:
        raise ValueError('Intervals must be finite and positive')
    metadata = json.loads((root / 'exploration.json').read_text(encoding='utf-8'))
    if metadata['preflight']['status'] != 'proceed':
        raise ValueError('Exploration did not pass preflight')
    with np.load(root / 'analysis' / METHOD / 'profile_distribution.npz', allow_pickle=False) as data:
        models = {s: np.array(data[s + '_median_mm']) for s in ('outer', 'inner')}
    output.mkdir(parents=True, exist_ok=True)
    all_summary = []
    applied_counts = {}
    for group in ('training', 'holdout'):
        profiles = load_group(root, metadata, group)
        applied_counts[group] = dict(total=len(profiles),
            applied=sum(p['record']['status'] == 'applied' for p in profiles))
        rows = measure(profiles, models, interval_mm)
        summary = summarise(rows, group)
        all_summary.extend(summary)
        csv_write(output / f'{group}_paired_station_distances.csv', rows)
        for detail in (False, True):
            name = 'registration_comparison_detail' if detail else 'registration_comparison'
            overlay(profiles, models, metadata['target_landmarks_mm'],
                    output / f'{group}_{name}.png', root.name + ' / ' + group, detail)
        angles, heights = statistical_plot(rows, summary, group, height_bin_mm,
            output / f'{group}_registration_statistics.png', root.name + ' / ' + group)
        csv_write(output / f'{group}_angular_statistics.csv', angles)
        csv_write(output / f'{group}_height_statistics.csv', heights)
    csv_write(output / 'registration_distance_summary.csv', all_summary)
    note = dict(version=VERSION, method=METHOD, source_input_sha256=metadata['input_sha256'],
                interval_mm=interval_mm, height_bin_mm=height_bin_mm, counts=applied_counts,
                distance='unsigned source point to same-surface representative line segment',
                sampling='identical raw-arc stations in both frames; transformed stations are not resampled',
                height_bins='assigned from RAW z in both frames', std_ddof=0,
                units='mm', transform='recorded constrained shared warp only; no pose transform',
                rejected='identity, retained in both panels and all statistics')
    (output / 'registration_comparison.json').write_text(json.dumps(note, indent=2), encoding='utf-8')
    lines = ['# 全体断面の原位置・位置合わせ後の比較', '', f'対象：{root.name}。単位：mm。', '',
        '左は特徴点による座標変形を行っていないr,z断面、右は保存済みの頸部・口唇の限定変形を再現した断面です。'
        '断面平面を共通の半径–高さ平面へ投影しており、原入力の3D頂点座標を直接描く図ではありません。'
        '軸と入力高さを保持し、追加の姿勢正規化・平行移動・相似変換は行っていません。', '',
        '左右は同じ断面群・同じ軸範囲・同じ代表モデルです。青は変形を適用した方向、黄褐色は拒否して原形状を保持した方向、'
        '赤は共通の代表輪郭、黒星は目標頸部・口唇です。底部保護域は固定です。', '']
    for group in ('training', 'holdout'):
        count = applied_counts[group]
        lines += [f'## {group}', '', f"断面数 {count['total']}、変形適用 {count['applied']}、拒否 {count['total']-count['applied']}。", '',
            f'![全体比較]({group}_registration_comparison.png)', '',
            f'![頸部・口縁拡大]({group}_registration_comparison_detail.png)', '',
            f'![距離統計]({group}_registration_statistics.png)', '',
            '| 範囲 | 状態 | 平均 | SD | RMS | p95 | 最大 | 測点数 |',
            '| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
        for r in all_summary:
            if r['dataset'] == group and r['surface'] == 'combined':
                lines.append(f"| {r['region']} | {r['frame']} | {r['mean_mm']:.3f} | {r['std_mm']:.3f} | "
                             f"{r['rms_mm']:.3f} | {r['p95_mm']:.3f} | {r['max_mm']:.3f} | {r['n_stations']} |")
        lines.append('')
    lines += ['## 統計と読み方', '',
        '距離は各表面の測点から、同じ表面の代表輪郭線分までの符号なし最近傍距離です。器厚・C2C距離ではありません。'
        '内外面を合算した棒グラフ、方位別RMS、原位置の高さ区間別平均±SD、累積分布を示します。'
        '高さ区間と部位は変形前の断面から固定し、左右で同じ測点を使用します。SDは母標準偏差で、信頼区間ではありません。', '',
        '全体統計は測点をプールした値です。弧長が長い断面・部位ほど測点数が多くなります。'
        '元の検証CSVは変形後の輪郭を再標本化しているため、この同一測点による比較とは小幅な差が生じ得ます。', '',
        '右図で散らばりが小さくなることは、特徴点を揃える処理の効果を示します。'
        '真の器形への精度改善、歪みの完全除去、器厚保存を証明するものではありません。'
        '原位置で上部の偏差が大きいかどうかも実測結果から読み取り、その傾向を強制していません。', '',
        'trainingはモデル生成に使用した方向、holdoutは同じメッシュの未使用方向です。'
        'holdoutも共通の目標特徴点に合わせるため、独立した実物計測による検証ではありません。', '']
    (output / 'RegistrationComparison.md').write_text('\n'.join(lines), encoding='utf-8')
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('roots', nargs='+', type=Path)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--interval-mm', type=float, default=1.)
    parser.add_argument('--height-bin-mm', type=float, default=5.)
    args = parser.parse_args()
    if args.output_dir and len({root.name for root in args.roots}) != len(args.roots):
        parser.error('Input names must be unique with --output-dir')
    for root in args.roots:
        destination = args.output_dir / root.name if args.output_dir else None
        print(run(root, destination, args.interval_mm, args.height_bin_mm))


if __name__ == '__main__':
    main()
