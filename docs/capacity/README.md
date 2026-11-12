# 土器の容量計算

ArchMorphの土器解析プロジェクトMorphPotが、PotteryVolumeCalculatorの容量計算・断面・表面痕跡・破片境界ツールの保守先です。[English guide](README_EN.md) · [開発報告](../../DevelopmentReport.md) · [移管元とライセンス](../../THIRD_PARTY_NOTICES.md)

| 目的 | 入口 | 計算結果 |
|---|---|---|
| 3Dメッシュの最大保持容量 | `MorphPot.py volume`（voxel v1.4.0） | 最初に外へ溢れる直前の液体体積、溢水高、品質診断 |
| 断面による容量推定 | [`MorphPot.py sections`](../../README.md#radial-sections-guide) | 内面形状から single / optimized / angular / ellipse の4方式で推定 |
| 2D実測図による容量推定 | `pottery_drawing_capacity.py`（v0.6.0、Tk GUI） | 既知距離で縮尺を校正し、左右の内面線をトレースして Drawing-Single 容量を推定 |

器壁そのものの体積と、液体の保持容量は別です。断面法とボクセル法では内腔・口縁の扱いも異なるため、方式と条件を結果に併記します。

## セットアップと単一モデルの計算

MorphPotのリポジトリ直下で実行します。Python仮想環境を推奨します。macOS/Linuxは `source .venv/bin/activate`、Windows PowerShellは `.\.venv\Scripts\Activate.ps1` で有効化できます。

```bash
python -m venv .venv
# 仮想環境を有効化してから実行
python -m pip install -r requirements.txt
python MorphPot.py volume --help
python MorphPot.py volume "/path/to/pot_rev.ply" --unit auto --pitch 1 --output-dir pot_capacity
python MorphPot.py volume "/path/to/pot_rev.ply" --unit auto --validate --output-dir pot_pitch_validation
```

通常は姿勢正規化済みの三角形メッシュPLYを使います。+Zが上で、内面と底があり、口から外へ液体が出られる資料が対象です。欠損・穴・内面のないメッシュでは、実物の保持容量を推定できない場合があります。入力座標の単位は `--unit auto`（既定）で `<stem>.asset.json` または同じフォルダの `transform.json` から取得します。メタデータがなければ、実際の単位を `--unit mm|cm|m` で指定してください。姿勢Transformは再適用しません。ピッチは入力座標の単位に関係なくmmです。

従来の `python vessel_voxel_volume.py model.obj --unit m --pitch 1` も使えます。直接起動時、単位メタデータのない入力で `--unit` を省略すると従来どおりmmです。PyMeshLab照合が必要な場合は `python -m pip install -r requirements-qa.txt` と `--diagnose-env` を使い、照合失敗をエラーにするには `--require-pymeshlab` を付けます。

`--validate` は2・1・0.5 mmの独立実行で、`--pitch` と同時には指定できません。単一ピッチの `result.json`、複数ピッチの `validation_summary.csv` / `.json` で容量・溢水高・失敗理由を確認します。細分化で収束するかを確認し、外挿値は診断に用います。0.5 mmは計算資源を多く使うため、まず1–2 mmで実行することを勧めます。詳細なアルゴリズムと出力構成は[英語ガイド](README_EN.md)、日本語の詳細手順・FAQは[移管前の固定スナップショット](https://github.com/kotdijian/PotteryVolumeCalculator/blob/93fbdbad24f22c38423d9e1b4780622b2ef8ccab/README.md)を参照してください。旧手順のインストール先はMorphPotに読み替えてください。

<a id="aom-pose-comparison"></a>

## AOM出力を使った姿勢別容量比較

[ArtefactsOrthoMaker](https://github.com/kotdijian/ArtefactsOrthoMaker#pottery-pose)の「資料を書き出す」で、同一資料・同一入力メッシュをRim・Base・Sliceの各姿勢で出力します。初回は資料フォルダ、それ以降は `revisions/` 以下の各出力フォルダを指定します。各フォルダに `item_manifest.json`、`transform.json`、正規化済みPLYが必要です。

```bash
python MorphPot.py volume --compare-poses "AOM/output/cup" "AOM/output/cup/revisions/BASE_REV" "AOM/output/cup/revisions/SLICE_REV" --pitch 1 --pose-baseline rim --output-dir cup_pose_comparison
```

- 同じ資料ID・元メッシュ・座標単位の、異なる姿勢のフォルダを2つ以上指定します。出力先は未作成のディレクトリを指定します。
- 単位はAOM出力から取得します。PLYとTransformのファイルサイズ・SHA-256を照合し、行列は姿勢間角度の計算だけに読み取ります。
- 基準を省略すると、計算に成功したRim、Base、採用Sliceの順です。参考Sliceは `--include-reference-slice` 指定時のみ計算し、基準にはできません。保留・未評価のSliceは計算しません。
- `pose_comparison.csv` / `.json` に角度差・容量差・溢水高・底面からの液面高・状態を、姿勢別サブフォルダに個別QAを保存します。成功が2姿勢未満、正式な比較対象がない、または計算失敗がある場合は、集計を残してエラー終了します。

姿勢によって重力方向に対する口縁の最低位置が変わるため、保持容量も変わり得ます。この比較は姿勢の影響を見るためのもので、単独で正しい姿勢や真の容量を決めるものではありません。

## 検証記録と旧版資料

- [ボクセル溢水法 v1.3.1 開発・検証報告](validation/PotteryVolumeCalculator_development_validation_report_v1.3.1.md) — 0015 Jinmenのピッチ収束と溢水高。v1.4.0の姿勢比較の実資料検証とは区別します。
- [表面痕跡解析 v0.1.0の自己検証](validation/surface-trace-v0.1.0.md) · [MorphPotの口縁・形状検証](../../ValidationReport.md)
- 旧版手順：[放射断面 v0.4.0](archive/radial-sections-v0.4.0.md) · [図面容量 v0.5.0](archive/drawing-capacity-v0.5.0.md) · [断面・図面ワークフロー v0.6.0（英語）](archive/section-workflows-v0.6.0.md)
- [0015 Jinmenの入力メッシュ](https://github.com/kotdijian/PotteryVolumeCalculator/blob/93fbdbad24f22c38423d9e1b4780622b2ef8ccab/0015Jinmen_small.ply)と[旧版の放射断面出力例](https://github.com/kotdijian/PotteryVolumeCalculator/tree/93fbdbad24f22c38423d9e1b4780622b2ef8ccab/examples/0015Jinmen_small_RadialSections_30deg)は、再現用の固定スナップショットを参照してください。
