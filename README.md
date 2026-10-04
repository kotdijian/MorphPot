# MorphPot

**姿勢正規化済みの土器3Dモデルを対象とする、容量・器軸・断面・表面形態解析プロジェクトです。**

基本入力は[ArtefactsOrthoMaker](https://github.com/kotdijian/ArtefactsOrthoMaker)が出力する三角形メッシュPLYです。既に適用された姿勢Transformを再適用しません。

現在は容量計算、水平・放射断面、断面容量、実測図容量、3Dガイド、破片再構成の既存モジュールを継承した開発版です。土器表面痕跡と接合境界候補の解析も収録しています。各処理はCLIまたは継承元GUIから独立して実行でき、MorphPot全体の統合GUIは今後開発します。

Surface Enhancement Labは共通手法の継承用としてexperimental/に収録しています。土器の投影方向・単位に適合させた統合は今後の予定です。

石器の平面形態・連続断面は[LithMorph](https://github.com/kotdijian/LithMorph)で開発します。

実装範囲、モジュール構成、開発計画、検証結果は[DevelopmentReport](DevelopmentReport.md)を参照してください。移管元は[SOURCE_PROVENANCE.json](SOURCE_PROVENANCE.json)に記録しています。

本リポジトリは開発版であり、全資料での精度検証を完了した正式安定版ではありません。

## PotteryRadialSections 操作ガイド

`pottery_radial_sections.py`（v0.6.0）は、水平断面の楕円中心から器軸のXY位置を推定し、その軸を通る放射状の縦断面を抽出します。断面からの容量計算も実装しています。現在はCLIで操作します。

### 1. 入力モデルと環境を準備する

OrthoMakerで姿勢正規化した、内外面を含む土器の三角形メッシュPLYを用意してください。このモジュールでは**Zを高さ、器軸方向をZと平行**に扱います。これは現行モジュールの処理規約です。任意軸の指定は今後対応します。

モデルを再回転・再正規化する必要はありません。器軸の傾きは診断値として出力し、自動補正しません。内面のないモデルや断面が閉じない破片では、内面による器軸推定・容量計算が成立しない場合があります。

リポジトリのフォルダで、次を実行します（macOS／Linux）。

```bash
python3 -m venv venv
source venv/bin/activate
python -m pip install -r requirements.txt
```

Windowsでは環境の有効化を `venv\Scripts\Activate.ps1` に読み替えてください。この断面解析にはGUI用の追加依存は不要です。

### 2. 基本操作：器軸・断面・容量をまとめて計算する

```bash
python MorphPot.py sections "/path/to/pot001_rev.ply" --unit auto
```

`/path/to/pot001_rev.ply` は実際の入力ファイルのパスに置き換えてください。空白を含むパスは引用符で囲みます。

`auto` は同名の `pot001_rev.asset.json`、または隣接する従来の `transform.json` から単位を取得します。metadataの読込契約は[DevelopmentReport](DevelopmentReport.md#4-入出力単位の現段階)を参照してください。単位が不明なら停止します。metadataがない場合は、**座標値の実際の単位**を指定します。

```bash
python MorphPot.py sections "/path/to/pot001_rev.ply" --unit m
```

正規化済みPLYでも単位がmmとは限りません。`--unit` は入力単位の宣言であり、姿勢Transformを再適用するオプションではありません。`--z-step-mm` など末尾が `-mm` の間隔は、入力単位にかかわらずmmで指定します。

既定設定では、水平断面20箇所から内面の楕円中心を求め、外れ値除去後の平均XYを器軸位置とします。縦断面の角度間隔は30°で、全断面6方向・半断面12方向を出力します。容量は4方式を計算し、参照用PNGも生成します。

### 3. 目的に合わせて設定する

**容量を計算せず、器軸と断面だけを取得：**

```bash
python MorphPot.py sections "/path/to/pot001_rev.ply" --unit auto --volume-mode none
```

**水平断面を5mm間隔、縦断面を15°間隔にする：**

```bash
python MorphPot.py sections "/path/to/pot001_rev.ply" --unit auto --z-step-mm 5 --angle-step 15 --output-dir "results/pot001_sections15"
```

**内外面を比較して器軸を確認：**

```bash
python MorphPot.py sections "/path/to/pot001_rev.ply" --unit auto --axis-surface both --volume-mode none
```

`both` は内外面の中心差を確認する設定で、最終器軸には内面中心を使います。外面のみで器軸位置を求める場合は `--axis-surface outer` を指定します。外面モードでも器厚から内面を自動復元する機能はありません。

**容量方式を選択：**

```bash
python MorphPot.py sections "/path/to/pot001_rev.ply" --unit auto --volume-mode single angular --single-angle 45 --volume-z-step-mm 0.5
```

| 容量方式 | 内容 |
| --- | --- |
| `single` | 指定した1枚の縦断面の対向する内面profileから計算。`--single-angle` で方向を指定 |
| `optimized` | 複数方向の内面半径から外れ値を除き、代表的な回転体として計算 |
| `angular` | 方位ごとの内面半径を角度方向に積分して計算 |
| `ellipse` | 水平断面の内面楕円面積を高さ方向に積分して計算 |
| `all` / `none` | 4方式すべてを計算／容量計算を省略 |

これらは断面profileによる容量推定です。voxel法は別の `MorphPot.py volume` で実行します。断面容量には内面・底・口縁の取得状態が影響するため、方式間の差と結果の `status` を確認してください。

### 4. 主なオプション

| オプション | 既定値 | 用途 |
| --- | --- | --- |
| `--axis-surface` | `inner` | 器軸推定に使う面：`inner` / `outer` / `both` |
| `--z-sections` | `20` | 高さを等分し、各区間中央に水平断面を配置 |
| `--z-step-mm` | 未指定 | 水平断面を絶対間隔で配置。`--z-sections` とは併用不可 |
| `--angle-step` | `30` | 縦断面の角度間隔。180°を割り切れる正の値を指定（例：30、15、10、5） |
| `--start-angle` | `0` | 放射断面の開始方位。0°は+X、90°は+Y |
| `--center-outlier-mad` | `3.0` | 水平断面中心の外れ値除去に使うMAD倍率 |
| `--contour-spacing-mm` | `0.5` | 楕円fit用の輪郭再標本化間隔 |
| `--sample-spacing-mm` | `0.5` | 出力断面点群の標本化間隔 |
| `--volume-mode` | `all` | 容量方式。複数指定可 |
| `--single-angle` | `0` | `single` 用縦断面の方位。定期抽出角度にない方向も計算可 |
| `--volume-z-step-mm` | `0.5` | 断面容量積分用の高さ刻み。器軸推定用の水平断面間隔とは別 |
| `--min-angular-valid-fraction` | `0.75` | `optimized` / `angular` の各高さで必要な有効方向の割合 |
| `--output-dir` | 自動設定 | 出力先を明示 |
| `--no-visualization` | 無効 | 参照PNGの生成を省略 |

角度間隔を細かくすると断面数が増えます。点群間隔や容量の高さ刻みを小さくすると、処理時間・出力サイズが増えます。入力メッシュ以上の細部が復元されるわけではありません。

### 5. 出力と確認の順序

既定の出力先は、入力PLYと同じ場所の `<入力名>_RadialSections_<角度間隔>deg/` です。例：`pot001_rev_RadialSections_30deg/`。同じ入力・角度間隔で再実行すると同じ出力先を使うため、条件比較には別々の `--output-dir` を指定してください。

| 出力 | 内容 |
| --- | --- |
| `metadata.json` | 入力単位、器軸位置、処理条件、容量結果 |
| `sections_summary.csv` | 各断面の方向、交差線分数・標本点数 |
| `axis_estimation/axis_summary.csv` | 採用断面数、中心のばらつき、診断用の傾き |
| `axis_estimation/horizontal_sections.csv` | 水平断面ごとの内外面楕円fit・中心・QA |
| `axis_estimation/horizontal_sections/` | 水平交差点とfit楕円のPLY |
| `axis_estimation/rotation_axis_edges.ply` | 推定器軸の線分PLY |
| `full_sections/` | 器軸を通る全縦断面のCSV、線分PLY、点群PLY |
| `radial_half_sections/` | 各方位の半断面のCSV、線分PLY、点群PLY |
| `visualization/` | 断面斜視図、XZ／YZ器軸検証図、中心分布、集計表PNG |
| `volume/` | 容量一覧CSV／JSON、方式別profile、内面点群・参照PNG（計算方式による） |

1. `visualization/axis_validation_xz.png`、`axis_validation_yz.png`、`axis_centers_xy.png` で器軸と中心の散らばりを確認します。
2. `axis_summary.csv` の採用断面数・RMS中心偏差・傾きを確認し、内外面の中心差は `horizontal_sections.csv` で確認します。
3. 全断面・半断面のPLYを入力メッシュと重ねて確認します。**出力PLYは入力と同じXYZ座標系・単位**です。点群PLYと線分PLYは三角形メッシュではありません。
4. `volume/volume_summary.csv` または `.json` で容量（L）、方式、`status`、失敗理由を確認します。計算処理の終了だけで全方式の成功を判断しないでください。

`--no-visualization` 指定時はPNGを生成しないため、CSV・PLYで確認します。古い出力を残したフォルダでは、過去のPNGなどが残る場合があります。

### 6. エラーが出た場合・詳細ヘルプ

| 状況 | 確認すること |
| --- | --- |
| 単位不明で停止 | metadataの配置、または実際の単位を `--unit mm` / `cm` / `m` で指定 |
| assetのハッシュ／単位不一致 | 対応するPLYとmetadataを揃える。明示単位指定でもassetの不一致は無視しない |
| 器軸推定に失敗 | Z方向の姿勢、内外面の有無、水平断面の閉合状態を確認。外面の位置推定が目的なら `outer` を検討 |
| `angle step must divide ...` | `--angle-step` を180°を割り切れる値へ変更 |
| 容量の一部が `failed` | 方式別の `note` を確認。内面profile、対向面、水平楕円の有効数などを確認 |

モジュール単体の詳細ヘルプは次のコマンドで表示します。

```bash
python pottery_radial_sections.py --help
```

単体実行も可能ですが、`--unit` は必須で `auto` は使えません。metadataによる自動解決には上記の `MorphPot.py sections` を使ってください。

```bash
python pottery_radial_sections.py "/path/to/pot001_rev.ply" --unit m --angle-step 30
```
