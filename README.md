# MorphPot

**姿勢正規化済みの土器3Dモデルを対象とする、容量・器軸・断面・表面形態解析プロジェクトです。**

基本入力は[ArtefactsOrthoMaker](https://github.com/kotdijian/ArtefactsOrthoMaker)が出力する三角形メッシュPLYです。既に適用された姿勢Transformを再適用しません。

現在は容量計算、水平・放射断面、断面容量、実測図容量、3Dガイド、破片再構成の既存モジュールを継承した開発版です。土器表面痕跡と接合境界候補の解析も収録しています。各処理はCLIまたは継承元GUIから独立して実行でき、MorphPot全体の統合GUIは今後開発します。

Surface Enhancement Labは共通手法の継承用としてexperimental/に収録しています。土器の投影方向・単位に適合させた統合は今後の予定です。

石器の平面形態・連続断面は[LithMorph](https://github.com/kotdijian/LithMorph)で開発します。

実装範囲、モジュール構成、開発計画は[DevelopmentReport](DevelopmentReport.md)、口縁部の検証方法・結果・残る課題は[検証報告](ValidationReport.md)を参照してください。移管元は[SOURCE_PROVENANCE.json](SOURCE_PROVENANCE.json)に記録しています。

本リポジトリは開発版であり、全資料での精度検証を完了した正式安定版ではありません。

## 目次

- [器形条件に基づく全体モデル生成法の探索](#whole-shape-exploration)

- [口縁付近の器厚交点・対応の診断](#whole-thickness-diagnostics)

- [全体モデルの寸法・器厚検証](#whole-vessel-dimension-validation)

- [全体代表モデルの独立検証スクリプト](#whole-vessel-validation)
- [個体の全体代表モデル（独立CLI）](#whole-vessel-model)
- [PotteryRadialSections 操作ガイド](#radial-sections-guide)
- [1. 入力モデルと環境を準備する](#setup)
- [2. 基本操作：器軸・断面・容量をまとめて計算する](#basic-usage)
- [3. 目的に合わせて設定する](#purpose-settings)
- [4. 主なオプション](#options)
- [5. 出力と確認の順序](#outputs)
- [6. エラーが出た場合・詳細ヘルプ](#troubleshooting)
- [7. 口縁部中央線を右側へ揃え、個体内標準モデルを作る（v0.10.0）](#rim-midline)
- [8. 抽出元の断面と中間点を重ねて検証する（v0.8.2）](#source-section-overlay)
- [9. 内外面を含む標準断面モデル（v0.9.0）](#standard-section-model)
- [10. 変換別の出力構成・終端モード・形状候補（v0.10.0）](#standardization-modes)
- [Similarityとaffineの違い](#similarity-vs-affine)
- [フォルダとファイル群](#standardization-files)
- [CLIモードと検出条件](#end-detection)
- [現時点の未実装と次段階](#next-development)
- [11. 標準断面の距離測定の再検証（v0.10.1）](#thickness-measurement)
- [12. 口縁部の検証版（現行 RimValidation 0.2.2-dev）](#rim-validation)
- [13. 全断面サマリーと1°出力による開始角度・間隔検証（RimValidation 0.2.0-dev）](#phase-validation)
- [14. 間隔依存性の直接比較画像（RimValidation 0.2.1-dev）](#interval-dependence)
- [15. 曲率・旋回角の比較（RimValidation 0.2.2-dev）](#curvature-turning)
- [16. 検証プログラムの選択と現在の検証結果（2026-10-06 JST）](#validation-programs)

検証プログラムの外部参照には、以下の固定URLを使用できます。バージョン表記を含む見出しが変わっても同じアンカーを使用します。

- [検証プログラム一覧・現在の結果](https://github.com/kotdijian/MorphPot/blob/main/README.md#validation-programs)
- [器厚・元断面との比較（validate／compare／sweep）](https://github.com/kotdijian/MorphPot/blob/main/README.md#rim-validation)
- [開始位置・角度間隔の検証（phases）](https://github.com/kotdijian/MorphPot/blob/main/README.md#phase-validation)
- [間隔依存性の画像追加](https://github.com/kotdijian/MorphPot/blob/main/README.md#interval-dependence)
- [曲率・旋回角の追加検証](https://github.com/kotdijian/MorphPot/blob/main/README.md#curvature-turning)

詳細な検証方法・結果は [ValidationReport.md](ValidationReport.md) を参照してください。

<a id="radial-sections-guide"></a>

## PotteryRadialSections 操作ガイド

`pottery_radial_sections.py`（v0.10.1）は、水平断面の楕円中心から器軸のXY位置を推定し、その軸を通る放射状の縦断面を抽出します。断面からの容量計算と、右側基準の口縁部中央線の抽出・個体内標準化を実装しています。現在はCLIで操作します。

<a id="setup"></a>

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

<a id="basic-usage"></a>

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

<a id="purpose-settings"></a>

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

<a id="options"></a>

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
| `--no-section-overlay` | 無効 | XY重ね合わせPLYと中央値ポリラインの自動出力を省略 |

角度間隔を細かくすると断面数が増えます。点群間隔や容量の高さ刻みを小さくすると、処理時間・出力サイズが増えます。入力メッシュ以上の細部が復元されるわけではありません。

<a id="outputs"></a>

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
| `section_overlay/all_full_sections_xy_points.ply` | 全縦断面を同一XY平面に重ねた方向別カラー点群 |
| `section_overlay/all_full_sections_xy_edges.ply` | 全縦断面を同一XY平面に重ねた線分 |
| `section_overlay/median_full_section_xy.ply` | 左右の内外面を含む中央値の閉じたポリライン（頂点＋edge） |
| `section_overlay/median_full_section_xy.csv` | 頂点順序、外面／内面の区分、入力単位座標・mm座標 |
| `section_overlay/overlay_qa.json` | 投影規約、口唇位置、中央値の採用／除外断面と理由 |

**v0.7.0のXY重ね合わせ・中央値輪郭：** 既定で自動出力します。通常の断面PLYとは異なり、Xは器軸からの符号付き距離、Yは元モデルのZ（高さ）、Zは0です。器軸はX=0に揃え、高さの原点・入力単位は保持します。左右を平均して対称形にする処理は行いません。`--no-visualization` や `--volume-mode none` でも生成します。

中央値は各全断面の正側／負側の口唇で輪郭を区切り、「正側口唇→外面・外底→負側口唇」と「負側口唇→内面・内底→正側口唇」を、それぞれ正規化弧長で再標本化した対応点の座標中央値です。同じ高さで半径を中央値にする方式とは異なります。頂点0は正側口唇、最後のedgeが頂点0へ戻ります。PLYを読み込むソフトがedge表示に対応している必要があります。

口唇は各側の幾何学的最高点として自動検出し、平坦な口唇の線分では中央を使います。考古学的な口唇の同定を保証するものではありません。波状口縁、装飾、欠損がある場合は口唇位置を確認してください。中央値には、連続した閉輪郭で内外底を区別できる断面のみを採用します。開曲線・複数成分・分岐は補間して接続せず、除外理由をQAへ記録します。有効断面が2方向未満なら中央値PLY／CSVは生成しません。重ね合わせPLYには除外された断面も含めます。

1. `visualization/axis_validation_xz.png`、`axis_validation_yz.png`、`axis_centers_xy.png` で器軸と中心の散らばりを確認します。
2. `axis_summary.csv` の採用断面数・RMS中心偏差・傾きを確認し、内外面の中心差は `horizontal_sections.csv` で確認します。
3. 全断面・半断面のPLYを入力メッシュと重ねて確認します。**出力PLYは入力と同じXYZ座標系・単位**です。点群PLYと線分PLYは三角形メッシュではありません。
4. `volume/volume_summary.csv` または `.json` で容量（L）、方式、`status`、失敗理由を確認します。計算処理の終了だけで全方式の成功を判断しないでください。
5. `section_overlay/overlay_qa.json` の有効方向数・口唇位置を確認し、中央値を重ね合わせ点群と比較します。少数の有効方向からの中央値が全体を代表するとは限りません。

`--no-visualization` 指定時はPNGを生成しないため、CSV・PLYで確認します。古い出力を残したフォルダでは、過去のPNGなどが残る場合があります。

<a id="troubleshooting"></a>

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

<a id="rim-midline"></a>

### 7. 口縁部中央線を右側へ揃え、個体内標準モデルを作る（v0.10.0）

目的は、一個体の口縁部周辺の半断面を比較し、位置・向き・倍率・限定的な歪みを除いた標準モデルを作ることです。左右の半断面を右側へ揃え、口唇から指定した幾何的な移行位置＋バッファーまでに限定します。既定autoは半径方向の停滞・反転と高さ方向の水平化を比較します。胴部全体は位置合わせに使いません。

v0.7.0の `median_full_section_xy.ply` は内外面を含む全断面の中央値輪郭です。v0.8.1の**中央線**は、口縁部周辺の内外面を輪郭順序を保った局所DTWで対応付け、対応点の中間を取った線です。高さの単調性を要求しません。対応は正規化弧長位置の25%帯に制限し、各壁を最大512点で標本化します。厳密なmedial axisや法線器厚の中間線ではありません。口唇位置、対応、比較終端をQAとPLYで確認してください。

通常の `sections` 実行で、次の処理を自動で行います。

1. 全縦断面の内外面輪郭から左右の半断面を取り出す。
2. 左半断面を器軸に対して反転し、右半断面と同じ向きへ揃える。Xは右向きの半径、Yは元モデルのZ高さ、投影Zは0。
3. 元内外面で終端候補を検出して方位間の支持を集計し、個体全体に共通の終端方針を決める。その局所区間で中央線を作り、同じ種類の移行を中央線でも確認する。半径方向モードは上下の増減を条件にしない。水平移行モードは高さ変化の絶対値を使う。
4. 変化点の先に局所的な**対応内外面間距離の2倍**の弧長バッファーを加え、そこで切る。固定mm指定も可能。
5. 各比較区間を正規化弧長で同じ129点へ再標本化する。
6. 採用された右側中央線群の座標中央値を固定参照とし、相似変換と制約付きアフィン変換でそれぞれ位置合わせする。
7. 位置合わせ後の座標中央値を標準中央線として出力し、点ごとの残差と変換量を保存する。

v0.9.0では、口唇の最高点を中央線の始点にしません。最高点は内外面輪郭を分けるためにだけ使います。初期対応線の口唇近傍の巻き込みを避け、安定した局所中央線の方向を求め、その外向き延長と元輪郭の最初の交点を「端部・始点」（頂点0）とします。交点から安定区間までは直線補間し、以後はDTW中間点を保持します。交点を求められない場合は最高点に戻さず、除外理由をQAへ記録します。全断面中央値輪郭（section_overlay）の従来の最高点登録は別処理として保持しています。

始点推定の支持位置は概ね対応壁間距離の1.5倍以上（平滑化長・標本間隔による下限あり）、方向推定区間は壁間距離または平滑化長以上です。始点、旧最高点、支持位置、方向、交差した元線分をrim_qa.jsonへ保存します。`raw_tip_extensions_xy.ply` の黒い線は新始点から支持位置までの延長確認用です。口唇形状が強く曲がる場合は支持位置と交点を視認確認してください。

半径方向モードは上下方向の増減を停止条件にしません。水平移行モードは傾きの減少を判定し、上→下だけを絶対的な追跡条件にはしません。停止判定用の平滑化は中央線本体へ適用しません。`--rim-end-mm` はこの新始点からの弧長となるため、旧版と同じ値でも比較区間が変わります。

**バッファーを3mmに指定：**

```bash
python MorphPot.py sections "/path/to/pot001_rev.ply" --unit auto --volume-mode none --rim-buffer-mm 3
```

**変化点がない直立口縁などで、口唇から15mmの弧長を手動指定：**

```bash
python MorphPot.py sections "/path/to/pot001_rev.ply" --unit auto --volume-mode none --rim-end-mm 15
```

変化点が見つからない場合や、バッファー分の線が足りない場合は、胴部へ自動延長せずその半断面を除外します。`--rim-end-mm` は口唇からの弧長であり、垂直深さではありません。設定すると自動の変化点・バッファー指定に代わります。

| オプション | 既定値 | 内容 |
| --- | --- | --- |
| `--no-rim-standardization` | 無効 | 口縁部抽出・標準化を省略。全断面overlayの省略とは独立 |
| `--rim-end-mode` | `auto` | 共通自動推定／`radial_turn`／`horizontal_transition`／`manual` |
| `--rim-horizontal-angle-deg` | `10` | 水平化の傾き許容角。先行する斜め区間はこの角度+10°以上 |
| `--rim-buffer-mm` | 自動 | 変化点の先に追加する弧長mm。0も指定可 |
| `--rim-buffer-thickness-ratio` | `2` | 自動バッファーの対応内外面間距離倍率 |
| `--rim-end-mm` | 未指定 | 口唇からの比較終端弧長を手動指定 |
| `--rim-smooth-mm` | `2` | 方向判定の平滑化・持続長。mm |
| `--rim-x-progress-tol` | `0.1` | 水平進行の許容幅。dx/ds≥−指定値が持続すれば停滞・反転と判定 |
| `--rim-turn-angle-deg` | 未指定 | 旧CLI互換。sin(角度)を水平進行許容幅に使う。上下方向は判定しない |
| `--rim-points` | `129` | 比較区間の対応点数。個体間比較では揃える |
| `--rim-endpoint-weight` | `5` | 始点・終点それぞれの重み（内点1点を1とする）。完全一致の拘束ではない |
| `--rim-affine-anisotropy` | `0.1` | 方向別伸縮係数の制約。伸縮項はexp(±k)、abs(k)≤log(1+指定値) |
| `--rim-affine-shear` | `0.1` | せん断係数の絶対上限 |
| `--rim-affine-penalty` | `1` | 方向別伸縮・せん断を相似変換へ近づける正則化重み |
| `--rim-outlier-mad` | `3.5` | 標準断面の追加inliersモデルに使う断面単位の外れ値閾値 |
| `--no-rim-bimodal` | 無効 | A/B二峰性スクリーニングを省略 |
| `--rim-bimodal-min-profiles` | `20` | 二峰性スクリーニングの最少半断面数（10以上） |
| `--rim-bimodal-bic-delta` | `10` | 二成分モデルに必要なBIC改善量 |

アフィン制約の初期値は実資料による校正前の暫定値です。相似変換に対して反転を許さず、アフィン変換でも正の行列式を保持します。ほぼ直線の中央線ではアフィン変形を十分に推定できないため、相似変換へ戻しQAに記録します。参照は固定された右側中央線群の中央値であり、反復的な一般化Procrustes解析ではありません。

出力先は `rim_standardization/` です。v0.10.0から未変換・候補・共通QAは直下、相似変換関連は `standard_similarity/`、アフィン変換関連は `standard_affine/` へ保存します。以下の表はファイル名を示し、フォルダの詳細は第10節にまとめています。

| 出力 | 内容 |
| --- | --- |
| `raw_midlines_xy.ply` / `raw_midlines.csv` | 採用された比較区間の右向き中央線。左右・方位・対応点順序をCSVに保存 |
| `similarity_midlines_xy.ply` / `affine_midlines_xy.ply` | 相似／制約付きアフィン位置合わせ後の全中央線 |
| `standard_similarity_midline_xy.ply` / `standard_affine_midline_xy.ply` | 個体内標準中央線。始点から終点までの開ポリライン |
| `standard_*_midline.csv` | 標準中央線の入力単位／mm座標、点ごとのRMS偏差 |
| `transforms.json` | 各半断面の変換行列、mm平行移動、主伸縮倍率、補正前後RMS、制約到達・fallback |
| `rim_qa.json` | 検出条件、比較範囲、始終点、採用・除外理由、標準化条件 |
| `standard_*_unit_shape.csv` | 標準中央線を重心で中心化し、centroid size=1にした無次元座標 |

PLY座標の単位は入力と同じです。`transforms.json` の行列はmmの投影座標に作用し、入力PLYへの姿勢Transformではありません。実測器厚・容量・元メッシュはこの標準化では更新しません。

v0.8.1では、高さの逆行による中央線生成失敗を修正しました。内外面それぞれの水平停滞位置から局所対応の範囲を決め、検出用の余裕（原則2×平滑化長）を加えた局所線で中央線を生成し、最終的には中央線の停止位置＋バッファーで切ります。胴部全体を対応付け・位置合わせには使いません。

**次段階の個体間解析：** 無次元座標はサイズを除いたProcrustes比較の入力候補です。複数個体の回転位置合わせ、一般化Procrustes解析、一致度・変異の集計はまだ実装していません。対応点数だけでなく、口唇・変化点・比較終端の相同性を確認してから比較します。完成個体から成形・乾燥・焼成の原因を分離したり、焼成前の形を復元したりする処理ではありません。


<a id="source-section-overlay"></a>

### 8. 抽出元の断面と中間点を重ねて検証する（v0.8.2）

口縁部解析時に `rim_standardization/` へ自動出力します。追加オプションは不要です。

| ファイル | 内容・表示色 |
| --- | --- |
| `raw_source_outer_xy.ply` | 元断面の局所外面輪郭・赤 |
| `raw_source_inner_xy.ply` | 元断面の局所内面輪郭・青 |
| `raw_midlines_xy.ply` | 比較区間の中央線・緑 |
| `raw_pair_connectors_xy.ply` | 中央線の各点に対応する内外面点を結ぶ線・灰色 |
| `raw_paired_points.csv` | 方位、左右、点番号、対応内外面点と中間点のmm座標 |

同じ接頭辞のPLYを一緒に読み込むと重ねて確認できます。`similarity_` と `affine_` にも同じ4種類を出力し、それぞれの中央線と同じ変換を元輪郭・対応点にも適用します。PLYは頂点とedge要素の線データです。edgeと頂点色の表示に対応したビューアを使用してください。

元輪郭は原断面の頂点を保持した局所区間で、停止判定のための余裕も含みます。灰色の線は最終比較区間だけに出力し、その中点が緑の中央線上にあることを確認します。左右は右側向きに統一済みで、X=右向き半径、Y=元のZ高さ、Z=0、PLY単位は入力と同じです。CSVはmmです。輪郭ごとの並びは `raw_midlines.csv` のprofile_id順です。

中央線と対応点は共通の中央線弧長で補間します。v0.8.1の切り出し後の再標本化から変更したため、途中の標本位置には微小な差が出る場合があります。標準中央値線は複数輪郭の集約結果なので、単一の元断面の中点としては解釈しません。


<a id="standard-section-model"></a>

### 9. 内外面を含む標準断面モデル（v0.9.0）

相似・アフィンの各位置合わせ後に、標準中央線に内外面までの距離を集約した断面を自動出力します。中央線だけを出力する既存ファイルも保持します。

| 出力（`*` は `similarity` または `affine`） | 内容 |
| --- | --- |
| `standard_*_section_all_xy.ply` | 幾何的に有効な取得半断面をすべて使う標準断面 |
| `standard_*_section_inliers_xy.ply` | 断面単位の外れ値を除いた標準断面 |
| `standard_*_section_A_xy.ply` / `standard_*_section_B_xy.ply` | 二峰性が支持された場合の二群の標準断面 |
| `standard_*_section_<区分>_outer_xy.ply` / `_inner_xy.ply` / `_midline_xy.ply` | 同じモデルの外面・内面・中央線。赤／青／緑 |
| `standard_*_section_<区分>.csv` | 3線の入力単位・mm座標、内外面距離の中央値と分位点、法線、断面数、元輪郭法線交点の取得数 |
| `*_section_models.json` | 使用profile_id、外れ値スコア、二峰性判定と群、交点未取得時の代替件数 |
| `*_section_distribution.csv` | 各半断面・対応位置の内外面距離、実交点／代替推定の区別、採用と群。全モデルの法線による標本 |
| `*_bimodal_histogram.csv` | 判定を実行できた場合の第一主成分スコアの頻度表 |
| `raw_tip_extensions_xy.ply` / `similarity_tip_extensions_xy.ply` / `affine_tip_extensions_xy.ply` | 始点と局所支持位置の延長確認線・黒 |

**標準断面の生成**：各群の変換済み中央線から対応位置の座標中央値を求め、その接線に直交する法線を設定します。v0.10.1では距離の測定には各断面自身の中央線法線を使い、変換済みの元内外輪郭との両方向の最初の交点までの距離を取得します。各位置の距離を集約した後、中央値距離を標準中央線の法線に沿って両側へ配置して標準外面・内面とします。測定法線と配置法線を分離します。口唇近傍の元輪郭も使うため、中央線の始点延長に伴う直線補間だけから器厚を作りません。

法線交点が取得できない箇所は、変換済みの内外面対応点とその中点との差を、その断面自身の中央線法線へ投影した距離で補います。代替件数と箇所をJSON／CSVへ記録します。距離分布の2.5／25／50／75／97.5分位点は取得断面の経験分布であり、独立標本による信頼区間ではありません。新始点では内外面が一致し距離0です。標準断面PLYは表示用の閉輪郭とし、比較終端の内外面を結ぶ線は人工的な切断線です。外面・内面別PLYにはその線を含めません。

**全利用と外れ値除外**：オリジナル・復元・接合部の識別はまだありません。ラベルによる除外は行わず、輪郭・対応・始点・比較区間を求められた取得断面を等重みで使います。幾何的に処理できない断面は従来どおりQAに理由を記録し、利用数を明示します。allにも中央値集約を使うため、外れ値が少なければallとinliersの形は一致することがあります。

inliersは、中央線座標と両面距離の対応位置ごとの中央値からの偏差を1.4826×MADで標準化し、断面内95パーセンタイルのスコアが既定3.5を超える半断面を除外します。MADの下限は0.05mmまたは中央器厚の5%の大きい方です。6半断面未満では外れ値判定をせず、inliersも全利用とし理由を記録します。残数2未満の場合も全利用へ戻し、状態を記録します。基準座標・倍率は全断面の右側固定参照のままとし、inliersのために再フィットはしません。

**二峰性とA/B**：全断面の中央線座標・内外面距離を、列ごとに5〜95パーセンタイルへ制限してから第一主成分へ投影します。一成分と二成分の1次元Gaussian mixtureを決定的な複数初期値のEMで比較します。既定では半断面20本以上、BIC改善10以上、両群それぞれ全体の20%以上かつ5本以上、平均差／プール標準偏差2以上、混合密度に二つのピーク、EM収束の全条件でA/Bを出力します。Aは主成分スコアの低い群、Bは高い群です。二峰性を支持した場合、外れ値スコアは各群内で計算し、少数群をまとめて外れ値扱いすることを避けます。A/Bモデル自体には各群の全メンバーを使います。

これは相関する一個体の断面群に対する探索的な判定で、全ての局所二峰性を検出する検定ではありません。第一主成分以外の変異、三群以上、片側だけにある変異は別途検討対象です。判定条件は今後の実資料による校正が必要です。A/Bをオリジナル／復元に対応付けません。Gaussian mixtureとBICのモデル比較の参考：[scikit-learn公式例](https://scikit-learn.org/stable/auto_examples/mixture/plot_gmm_selection.html)。追加のscikit-learn依存は導入していません。

```bash
python MorphPot.py sections "/path/to/pot001_rev.ply" --unit auto --angle-step 5 --volume-mode none
```

この段階の標準断面は口縁部の2Dモデルです。個体間Procrustes解析、標準3D回転メッシュ、表面セグメンテーション・アノテーションは未実装です。次段階ではfragment boundaryの拡張と表面ラベルから、取得断面をオリジナル・復元等で分けて再集約します。


<a id="standardization-modes"></a>

### 10. 変換別の出力構成・終端モード・形状候補（v0.10.0）

<a id="similarity-vs-affine"></a>

#### Similarityとaffineの違い

| 項目 | similarity（相似変換） | affine（制約付きアフィン変換） |
| --- | --- | --- |
| 平行移動・回転 | 許可 | 許可 |
| 一様な倍率補正 | 許可 | 許可 |
| 方向別伸縮・せん断 | 許可しない | 制約・正則化の範囲で許可 |
| 角度・縦横比 | 保持 | 補正により変わり得る |
| 最適化 | 正の回転・倍率で右側固定参照へ合わせる | 相似解から開始し、方向別伸縮・せん断を追加調整 |
| 推定不能時 | 通常の相似fit | ほぼ直線の中央線／最適化失敗時は相似へ戻りstatusに記録 |

両者は同じ未変換中央線と、採用右側中央線群の座標中央値を固定参照として計算します。端点の重みは内点の5倍です。相似変換の出力へさらにアフィン変換を適用する二段階処理ではありません。アフィンの伸縮項はexp(±k)、abs(k)≤log(1.1)、せん断係数は±0.1が既定です。一様倍率は別に推定され、せん断を含む最終主伸縮倍率はtransforms.jsonで確認します。

各方式の変換を元内外輪郭にも適用して中央線・距離分布を集約するため、外れ値やA/B判定の結果も方式ごとに異なる場合があります。いずれも個体内標準化であり、個体間のサイズ除去Procrustes解析はまだ実装していません。

<a id="standardization-files"></a>

#### フォルダとファイル群

| 保存先 | ファイル群と役割 |
| --- | --- |
| `rim_standardization/` 直下 | `raw_midlines_xy.ply`／`raw_midlines.csv`：右向き・未変換中央線。`raw_source_outer_xy.ply`／`raw_source_inner_xy.ply`：元局所輪郭。`raw_pair_connectors_xy.ply`／`raw_paired_points.csv`：内外対応と中間点。`raw_tip_extensions_xy.ply`：端部延長確認線 |
| 同直下 | `rim_qa.json`：取得・除外、始終点、設定、方式別QA。`shape_candidates.json`：中立的な形状候補、支持方位、位置のばらつき、推奨・実際の終端方針。`endpoint_candidates.csv`／`endpoint_candidates_xy.ply`：候補・始点・最終終端の座標と色付き点。`transforms.json`：両方式の変換一覧 |
| `standard_similarity/` | `similarity_midlines_xy.ply`：変換後全中央線。`similarity_source_outer_xy.ply`／`similarity_source_inner_xy.ply`：同変換後の元輪郭。`similarity_pair_connectors_xy.ply`／`similarity_tip_extensions_xy.ply`：同変換後の検証線。`transforms.json`：相似の変換だけ |
| `standard_affine/` | 上記の`similarity`を`affine`に置き換えた変換後の輪郭・検証線と、アフィンだけの`transforms.json` |
| 各方式サブフォルダ | `standard_<方式>_midline_xy.ply`／`_midline.csv`：全採用中央線の標準中央値と偏差。`standard_<方式>_unit_shape.csv`：重心中心化・centroid size=1の形状 |
| 各方式サブフォルダ | `standard_<方式>_section_all_xy.ply`／`_inliers_xy.ply`：全採用／外れ値除外の標準断面。支持された場合のみ`_A_xy.ply`／`_B_xy.ply`。各区分の`_outer_xy.ply`／`_inner_xy.ply`／`_midline_xy.ply`と`.csv`：構成3線、距離分位点・元輪郭交点取得数 |
| 各方式サブフォルダ | `<方式>_section_models.json`：断面の採用・除外スコア、二峰性判定、A/B、代替距離の件数。`<方式>_section_distribution.csv`：全断面・対応位置の距離標本。`<方式>_bimodal_histogram.csv`：判定実行時のPC1頻度表 |

正式な綴りは`affine`です。旧版の直下にあった方式別の生成ファイルは再実行時に削除し、サブフォルダで再生成します。無効化・検出失敗時にも所有する旧生成ファイルを消去し、利用者の別名ファイルは保持します。v0.9.0までの直下パスを読む外部スクリプトは変更してください。

候補点PLYの色：橙=半径方向の停滞・反転、紫=水平移行、黒=新始点、緑=バッファーを加えた最終比較終端。PLYは入力単位の共通未変換XY平面です。CSVの`stage=wall_screening`は内外面候補座標・弧長の平均で、始点基準は元輪郭の最高点分割です。`stage=final_midline`は端部延長で決めた新始点からの中央線候補・弧長です。二段階の候補は厳密な同一点ではありません。profile_indexはrim_qa.jsonの行順で、採用中央線のIDは同JSON内のaccepted_profile_idです。

<a id="end-detection"></a>

#### CLIモードと検出条件

```bash
# メッシュの断面群から中立的な方針を推定
python MorphPot.py sections "/path/to/pot001_rev.ply" --unit auto --rim-end-mode auto --volume-mode none
# 半径方向の停滞・反転を終端とする
python MorphPot.py sections "/path/to/pot001_rev.ply" --unit auto --rim-end-mode radial_turn --volume-mode none
# 斜めの器壁から水平部への移行を終端とする
python MorphPot.py sections "/path/to/pot001_rev.ply" --unit auto --rim-end-mode horizontal_transition --volume-mode none
# 新始点からの比較弧長で指定
python MorphPot.py sections "/path/to/pot001_rev.ply" --unit auto --rim-end-mode manual --rim-end-mm 15 --volume-mode none
```

GUIは不要です。候補PLYとQAを既存ビューアで確認します。`--rim-end-mm`を指定するとautoや他のモードより優先してmanualとなります。manualモードはこの弧長指定が必須です。

- **radial_turn**：左方向進行が持続した後、正規化接線の半径成分が−rim-x-progress-tol以上になる状態が持続した位置。上下増減・90°以上のベクトル変化は要求しません。
- **horizontal_transition**：左進行を保った斜めの区間（水平に対する傾きが既定20°以上）が持続した後、高さ方向変化が小さくなり、傾き既定10°以内が持続する位置。最初から水平な口唇だけでは停止しません。上下の符号は制限しません。
- 両モードともrim-smooth-mmの長さで方向を平滑化・持続判定し、口唇直後の最低guard区間を設けます。移行位置から先のバッファーは従来どおり対応壁間距離2倍または固定mmです。中央線本体は停止判定用に平滑化しません。

内外面それぞれで候補を求め、同じ種類の候補が両面にある場合を半断面の支持とします。全方位の有効な元輪郭から支持割合を集計し、autoでは支持が最も多い終端方針を個体全体に適用します。同率なら候補弧長中央値が小さい方を選びます。支持差が0.15以下なら推奨が曖昧であることも記録します。方位ごとに異なる意味の終端を混ぜて標準モデルを作ることは避け、選ばれた方針に対応できない方位は除外理由を残します。設定したモードと自動推奨は両方を記録します。

形状候補は`radial_turn`（半径方向の停滞・反転）、`horizontal_transition`（高さ方向変化の減少・水平化）、`near_vertical`（直立傾向）、`compound`（複合）、`no_clear_transition`（明瞭な移行なし）、`undetermined`（判定保留）。直立傾向は口唇直後の局所区間でabs(半径接線成分)≤0.15の割合が75%以上の場合です。複数カテゴリを併記し、器種名・用途名へ強制分類しません。

支持割合は両面を処理できた半断面数に対する割合で、除外前の幾何候補も含みます。候補弧長の中央値・IQR、高さ位置のIQR、支持方位を記録します。校正された確率や独立標本の信頼区間ではありません。自動推奨は幾何的な処理方針であり、個体間比較区間の考古学的な相同性を保証しません。

<a id="next-development"></a>

### 現時点の未実装と次段階

丸底などで水平化が短い場合、片面だけに候補がある場合、複数輪郭や不連続輪郭は判定保留・除外となり得ます。既定閾値は実資料で校正中です。CLIで同じ終端方針を明示し、比較範囲を確認してください。GUIによる候補選択・終端編集、オリジナル／復元表面のセグメンテーションとアノテーション、個体間Procrustes解析、標準断面の3D回転モデルは今後の実装です。


<a id="thickness-measurement"></a>

### 11. 標準断面の距離測定の再検証（v0.10.1）

v0.10.0までは標準中央線の法線を各個別中央点に置いて距離を測定していました。個別の器壁方向と標準方向に差があると斜めの横断距離となり、平行壁では器厚が1/cos(方向差)倍になる問題を修正しました。

現在は、(1)変換済みの各個別中央線の法線を計算、(2)その法線と同じ個別元輪郭の交点距離を取得、(3)対応位置の距離分位点を集約、(4)標準中央線の法線に沿って中央値距離を配置、の順です。all/inliers/A/Bの選択にかかわらず各断面の測定方向を固定します。交点未取得時の対応点投影も個別法線を使います。これはDTW対応点から得た中央線を基準とする法線横断距離であり、厳密な表面法線器厚やmedial-axis厚ではありません。

`<方式>_section_distribution.csv`にmeasurement_normal_x／yを追加し、各区分の断面CSVにはmeasurement_vs_model_angle_p50_deg／max_deg、モデルJSONには全位置の角度差p95／maxと測定・配置フレームを記録します。抽出中央線とallの標準中央線の計算法は保持します。距離が変わるため外れ値スコアやA/B判定、対応する群の標準線は変わる場合があります。

検証例：一定厚4mm、方向差±30°の個別平行壁で、旧方式4.6188mm→新方式4.0000mm。提供0°–180°断面は左右2本のみで、修正前後の標準器厚差は最大約0.068mm。画像だけで今回の全方位の膨らみをこの問題に帰属しません。

`*_pair_connectors_xy.ply`は内外対応点を結ぶ線であり、元輪郭そのものではありません。口唇の新始点から支持位置までは対応点の直線補間を含むため、その外周と、元輪郭法線交点で再構成した丸い端部は一致しないことがあります。標準断面の肥厚を確認する際は、同じ方式サブフォルダの`*_source_outer_xy.ply`／`*_source_inner_xy.ply`も重ねて確認してください。比較終端を結ぶ人工切断線は外面・内面別PLYで除いて確認できます。多方位での局所方向、対応弧長、交点の枝選択の妥当性は継続検証対象です。


<a id="rim-validation"></a>

### 12. 口縁部の検証版（現行 RimValidation 0.2.2-dev）

公開版に向けた検証専用CLI `pottery_rim_validation.py` を追加しました。RadialSections本体はv0.10.1のままです。入力は既に生成した断面出力フォルダ、または分割数検証用のメッシュです。元のPLYへ姿勢Transformを再適用しません。

```bash
# 既存出力を始点から1mm間隔で検証（similarity/affineの両方）
python pottery_rim_validation.py validate sample_RadialSections_5deg \
  --output-dir validation --interval-mm 1

# 角度間隔だけを変えて再抽出し、最後の5度を比較基準にする
python pottery_rim_validation.py sweep sample.ply \
  --output-dir division_validation --steps 30 15 10 5 --interval-mm 1 \
  -- --unit m --z-step-mm 5 --rim-end-mode radial_turn

# 既存の複数条件を比較（最後のフォルダが基準）
python pottery_rim_validation.py compare run_30deg run_15deg run_5deg \
  --output-dir comparison --interval-mm 1 --mode similarity
```

`--mode similarity|affine|both`、`--selection all|inliers|A|B`（既定all）、`--no-plots` に対応します。A/Bは生成済みの場合のみ指定できます。sweepでは容量計算・従来の概観画像を省き、検証用の断面重ね合わせ画像は生成します。必要な単位・終端などの抽出オプションは `--` の後に渡します。

**位置の定義**：標準中央線の始点から指定したmm間隔で測点を置きます。始点そのものは単一の口唇接点なので器厚評価に含めません。`absolute_mm` は各変換済み中央線でも同じ弧長、`corresponding_u` は標準測点の対応インデックス比率uを各断面に移します。断面の長さが違うため両者は一致するとは限りません。どちらも各断面自身の局所法線と、元のメッシュ交線から保持した内外面ポリラインとの両方向交点から器厚を再測定します。交点なし・対象長さ外は欠測で、pair projectionによる代替や外挿は行いません。入力mm換算の元器厚も逆対応位置で測定し、相似・アフィンのサイズ補正と区別します。

| 出力 | 内容 |
| --- | --- |
| `validation_XX/standard_MODE_SELECTION/station_measurements.csv` | 断面別・測点別・位置定義別のモデル厚、変換後実測厚、元実測厚、誤差、中央線位置差、欠測理由 |
| `station_statistics.csv` | 各位置の有効数、平均・最小・最大・母標準偏差、モデルに対するbias/RMSE |
| `profile_geometry.csv` | 中央線RMS、始点・終端位置差、元と変換後の弧長、変換の主伸縮率 |
| `profile_NNN_overlay_xy.ply` / `.png` | 灰色の元メッシュ断面内外面と赤色の標準モデル内外面の重ね合わせ。コネクタや人工終端閉鎖線は使わない |
| `profile_NNN_model_xy.ply` | 同じ座標・入力単位の標準モデル単独PLY |
| `validation.json` | 比較定義、条件、断面数、制約 |
| `MODE_SELECTION_division_comparison.csv` | 分割条件間の共通弧長測点における中央線位置差、始点を揃えた位置差、厚み差、接線角度差 |
| `MODE_SELECTION_division_summary.json` | 分割条件間のRMS・最大差、弧長差、終端差、両条件の設定 |

測定ポリラインには終端検出の余裕範囲も残っています。画像の黒点は標準中央線上の測点です。PLYは入力単位、表・画像はmmです。72方向を設定していても抽出無効方向があれば検証数は減ります。分割間隔の比較では角度の開始位置、器軸推定条件、終端方針、中央線点数などを揃え、抽出側の `rim_qa.json` と検証有効数を併読してください。中央線の補間点数（`--rim-points`）は角度分割数とは別の条件です。必要なら点数だけ変更した出力をcompareで検証できます。

この比較はモデル構築に利用した同じメッシュ交線に対する記述的検証です。独立した実物計測に対する精度評価ではありません。高密度分割を正解とはみなさず、変化量を出力します。「変わらない」の許容幅は測定分解能・研究目的から別途指定する必要があります。autoで終端方針が変わった場合は分割数だけの比較にならないため、設定を確認して同じ方針で再実行します。


<a id="phase-validation"></a>

### 13. 全断面サマリーと1°出力による開始角度・間隔検証（RimValidation 0.2.0-dev）

`validate` は各方式・選択群のフォルダに `profile_all_overlay.ply` と `profile_all_overlay.png` を追加します。選択された全半断面の元内外面（灰）と標準モデルの内外面（赤）を同じ座標に重ねます。個別の `profile_NNN_overlay_xy.ply` / `.png` も補足資料用に引き続き生成します。`--no-plots` を指定した場合だけPNGを省きます。

```bash
# 180平面（0〜179°など）を試行した1°間隔の解析出力を指定
python pottery_rim_validation.py phases 0015Jinmen_small_RadialSections_1deg \
  --output-dir phase_validation --steps 5 10 15 --interval-mm 5

# どの間隔でも開始角度を0〜4°の5セットだけに限定する場合
python pottery_rim_validation.py phases 0015Jinmen_small_RadialSections_1deg \
  --output-dir phase_validation_5starts --steps 5 10 15 --phases 0 1 2 3 4 \
  --interval-mm 5
```

既定では5°の5位相、10°の10位相、15°の15位相、計30セットを生成します。欠落のない場合、各セットはそれぞれ72・36・24半断面です。左右は同じ断面平面のセットに属します。開始角度は元の1°グリッドの最小角からの相対オフセットです。明示的な `--phases` で10°・15°を5セットに絞る場合、その間隔の全位相を網羅する検証ではありません。

**再計算範囲**：1°出力の元内外輪郭、中央線、内外対応点を抽出します。器軸・口唇始点・各断面の終端範囲は1°抽出時のものに固定し、各セットの右側中央値基準、similarity/affine変換、外れ値判定、標準中央線と標準断面モデルを再計算します。既存の変換済みモデルを間引く処理ではありません。これは標準モデル構築の開始角度・間隔依存性を切り分ける実験で、セットごとに器軸推定や終端検出から再実行する総合的な感度検証とは区別します。

| 出力 | 内容 |
| --- | --- |
| `interval_05deg/phase_00deg/rim_standardization/` 等 | セットごとに再構築した元輪郭、変換、モデル、分布、QA。profile_idはセット内ID、source_profile_idで1°元IDへ対応 |
| 各セットの `validation/standard_MODE_SELECTION/` | 測点別比較、個別重ね合わせ、`profile_all_overlay.ply` / `.png` |
| `dense_reference_validation/` | 1°全断面モデルと全元断面の検証・サマリー重ね合わせ |
| `comparison_05deg/` 等 | 各位相モデルと1°全断面基準との位置・器厚差、RMS・最大差 |
| `summary_05deg_MODE_SELECTION/phase_station_statistics.csv` 等 | 全位相間の器厚平均・標準偏差・最小・最大、中央線の最大位相間距離 |
| 同フォルダの `phase_models_overlay.ply` / `.png` | 色分けした位相別モデルと、黒の1°基準モデルの重ね合わせ |
| `phase_experiment.json` | 間隔、開始角度、予定数・有効数、使用角度・元ID、固定／再計算条件 |

180平面の試行を `rim_qa.json` で確認できる1°出力が必要です。無効断面は補完せず、予定数と有効数を併記します。PLYは元入力単位、表とPNGはmmです。各セットの比較基準は同じ1°全断面モデルです。位相間ばらつきと間隔ごとの1°基準との差を併読します。位相セットや左右断面を独立標本とみなした信頼区間・有意差検定は行いません。


<a id="interval-dependence"></a>

### 14. 間隔依存性の直接比較画像（RimValidation 0.2.1-dev）

`phases` の終了時に、2種類以上の間隔を比較し、PNGが有効な場合は `interval_dependence/interval_dependence_similarity_all.png` と `interval_dependence_affine_all.png` を自動出力します。選択群がinliers等なら末尾の名前が変わります。

各画像は4パネルです。横軸は比較基準の始点からの弧長mm、間隔5/10/15°などを色分けします。器厚差（間隔モデル−共通1°基準モデル、符号あり）、中央線位置差、始点を揃えた中央線位置差、接線角度差を表示します。線は各位置での位相間中央値、半透明帯は位相間の最小〜最大です。帯は信頼区間ではありません。始点合わせ前後の位置差を比較して基準位置のずれと残る形状差を見分けます。

**画像だけ再出力**する場合、モデル再計算は不要です。既存のphase_validation内の比較CSVを使います。

```bash
python pottery_rim_interval_plot.py phase_validation --steps 5 10 15

# 出力先を指定、similarityだけを生成
python pottery_rim_interval_plot.py phase_validation \
  --output-dir interval_images --mode similarity --steps 5 10 15
```

`pottery_rim_interval_plot.py` は独立スクリプトです。この1ファイルとNumPy・Matplotlibだけで動作し、morphpotパッケージや入力メッシュは不要です。必要な入力は `comparison_05deg/similarity_all_division_comparison.csv` 等。`--steps` を省くと比較フォルダを検出します。既定は両方式、`--selection` はall/inliers/A/B、`--dpi` の既定は180。生成物はPNGのみです。既存CSV・モデルには書き込みません。

異なる基準モデルや基準断面数のCSVを混ぜると停止します。全間隔で共通する記録済み測点のみ描画し、欠測はゼロに置き換えません。短い位相モデルがあると測点ごとの有効位相数は減る場合があるため、詳細は元の比較CSVとphase_station_statistics.csvで確認します。共通測点以外の補間・外挿は行いません。図は間隔を変えたときの変化量を表示し、「安定」の許容差や統計的有意差を自動判定しません。


<a id="curvature-turning"></a>

### 15. 曲率・旋回角の比較（RimValidation 0.2.2-dev）

標準中央線の曲がり方そのものを比較するPNG・測点CSV・集計CSV・条件JSONを追加しました。phasesの画像生成時に `curvature_turning/` へ自動出力します。既存結果に追加分だけを生成する場合は、独立スクリプトを実行します。モデルの再計算・メッシュ抽出は行いません。

```bash
python pottery_rim_curvature_plot.py phase_validation --steps 5 10 15

# 全間隔を開始角度0〜4°の共通5セットに揃える場合
python pottery_rim_curvature_plot.py phase_validation \
  --steps 5 10 15 --phases 0 1 2 3 4 --output-dir curvature_common5

# 1°基準の保存場所が変わった場合は明示的に指定
python pottery_rim_curvature_plot.py phase_validation \
  --reference-dir ../PotteryVolumeCalculator/0015Jinmen_small_RadialSections_1deg
```

必要なのは既存phase_validationの `phase_experiment.json`、各セットの標準断面CSV（mid_x_mm・mid_y_mm）、元の1°出力の標準断面CSVです。参照先はmanifestから取得し、移動・別環境では `--reference-dir` で指定します。`pottery_rim_curvature_plot.py` はPython・NumPy・SciPy・Matplotlibだけで動作する単独ファイルです。

各方式・選択群について次を出力します。既定は両方式・allで、`--mode` / `--selection` に対応します。

| ファイル | 内容 |
| --- | --- |
| `curvature_turning_similarity_all.png`（affineも同様） | 上段：始点からの接線旋回角、符号付き曲率、累積絶対旋回角。下段：それぞれの1°基準との差 |
| `…_stations.csv` | 各間隔・開始角度・測点の曲率、符号付き／絶対旋回角、基準との差、端部フラグ |
| `…_summary.csv` | 共通区間の総旋回角、最大絶対曲率と位置、基準に対するピーク比、内部区間の曲率RMS差 |
| `curvature_turning_MODE_SELECTION.json` | 曲率・旋回角の定義、平滑化条件、共通終端、制約 |

黒線は1°基準、色線は位相間中央値、帯は最小〜最大で信頼区間ではありません。全間隔・全位相で共通する弧長区間を使います。

**測定方法**：各中央線を同じ弧長グリッド（既定0.25mm）へ補間し、同じ物理幅（既定2mm）の3次Savitzky–Golay局所多項式から1・2階微分を求めます。`--grid-mm` / `--smooth-mm` で変更でき、使用した実効幅はJSONに記録します。曲率は `(x'y''−y'x'')/(x'^2+y'^2)^(3/2)`、単位は1/mm。正負は始点から胴側へ進む向きの反時計回り／時計回りです。符号付き旋回角はunwrapした接線角から始点角を引き、絶対旋回角は各グリッド間の接線角変化の絶対値を累積します。逆向きの湾曲が相殺されないよう、両方を報告します。

端部は片側の多項式当てはめになるため、窓幅の半分を灰色表示・フラグ付けし、曲率ピークと曲率RMSの評価から除外します。総旋回角は共通区間の始点・終点の推定も含むため端部条件に影響されます。曲率は平滑化幅に敏感で、補間は元中央線の情報量を増やしません。2mm・4mmなどの幅でも比較してください。最大曲率比は区間内の最大値同士の比で、同じ解剖学的位置のピーク同士を自動対応付けした値ではありません。モデルの湾曲消失を自動判定しません。


<a id="validation-programs"></a>

### 16. 検証プログラムの選択と現在の検証結果（2026-10-06 JST）

検証CLIは `pottery_rim_validation.py`（0.2.2-dev）、画像の追加生成は `pottery_rim_interval_plot.py` と `pottery_rim_curvature_plot.py`（各0.1.0）です。RadialSections本体は0.10.1のままです。詳しい実験条件・結果・制約は[ValidationReport.md](ValidationReport.md)にまとめています。

| 目的 | プログラム／コマンド | 主な確認出力 |
| --- | --- | --- |
| モデルと元メッシュ交線の器厚・位置比較 | `pottery_rim_validation.py validate` | `station_measurements.csv`、`station_statistics.csv`、`profile_all_overlay.png`／`.ply`、個別重ね合わせ |
| 保存済みモデル同士を比較 | 同 `compare` | `MODE_SELECTION_division_comparison.csv`、`…_summary.json`。最後の入力を基準とする |
| メッシュから角度間隔を変えて再抽出 | 同 `sweep` | 各抽出・検証フォルダと条件間比較 |
| 1°出力から開始位置・間隔依存性を検証 | 同 `phases` | 位相別検証、`phase_models_overlay.png`、比較CSV、`phase_experiment.json` |
| 間隔依存性の画像だけ追加 | `pottery_rim_interval_plot.py` | 位置・器厚・接線差の4パネルPNG |
| 曲率・旋回角だけ追加 | `pottery_rim_curvature_plot.py` | 6パネルPNG、測点・集計CSV、条件JSON |

既存出力の場所は、現在の作業ディレクトリからの相対パスまたは絶対パスで指定します。入力は `rim_standardization/rim_qa.json` を持つ抽出出力フォルダ、またはその `rim_standardization` 自体です。フォルダ名を推測せず、実在する場所を指定してください。

```bash
# MorphPotから、隣のPotteryVolumeCalculatorにある保存済み出力を検証
python pottery_rim_validation.py validate "../PotteryVolumeCalculator/0015Jinmen_small_RadialSections_5deg" --output-dir validation --interval-mm 5

# 各間隔で開始位置0〜4°の共通5セットを作成
python pottery_rim_validation.py phases "../PotteryVolumeCalculator/0015Jinmen_small_RadialSections_1deg" --output-dir phase_validation --steps 5 10 15 --phases 0 1 2 3 4 --interval-mm 5

# 既存モデルから2mm／4mm窓の曲率・旋回角を追加生成
python pottery_rim_curvature_plot.py phase_validation --steps 5 10 15 --phases 0 1 2 3 4 --smooth-mm 2 --output-dir curvature_common5_smooth2
python pottery_rim_curvature_plot.py phase_validation --steps 5 10 15 --phases 0 1 2 3 4 --smooth-mm 4 --output-dir curvature_common5_smooth4
```

**ここまでの結果**：既知厚4mmの合成壁では、旧共通法線測定の4.6188mmを個別法線測定で4.0000mmへ修正できました。提供された実データCSVのモデル厚4mmに最も近い対応点（モデル4.0268mm）では、64半断面の平均4.0300mm、最小2.5136mm、最大4.9239mm、母標準偏差0.4072mmでした。72半断面の集計ではなく、同じ対応点における変換後の法線器厚です。モデル厚は内外距離の中央値を別々に集約するため、器厚の平均とは一致する必要がありません。

実土器 `0015Jinmen_small` の提示図（similarity／all、各間隔共通5位相）では、約23mmの主要曲率ピークと全体の約90°の符号付き旋回が5°・10°・15°で保持されています。累積絶対旋回角の1°基準との差は、2mm窓では概ね5°で20〜25°、10°で40°、15°で30°、4mm窓では5°で4〜7°、10°で10°、15°で15°へ縮小しました。これらは図からの概算で、精密な値は各 `…_summary.csv` で確認します。4〜7°は終端付近での推移を含む範囲です。

現段階で主要な湾曲の消失を示す明瞭な所見はありません。細かな方向変化には角度間隔・開始位置・平滑化幅の依存性があり、この個体の比較条件では5°が比較的1°基準に近い傾向です。粗い間隔ほど常に悪化するとは限らず、5°を全個体共通の既定値とする根拠はまだありません。4mm窓を真値とせず、2mmと4mmを異なる形状スケールとして併記します。

1°基準は真値ではなく、位相帯は最小〜最大で信頼区間ではありません。実土器の位相比較は利用者が実行し提示した画像の評価で、開発環境で元1°データから再計算した結果とは区別します。独立実測、未使用断面での照合、別個体・異なる器形、器軸／終端推定を含めた感度検証を次段階とします。


<a id="whole-vessel-model"></a>

## 個体の全体代表モデル（WholeModel 0.1.0-dev）

`pottery_whole_model.py` を追加しました。口縁部標準化とは別のCLI・別の出力フォルダで、個体全体の内外面と底部を持つ軸対称の代表メッシュを生成します。口縁部の相似／アフィン変換を延長せず、部位モデルの接合も行いません。姿勢Transformを再適用せず、物理寸法を保持します。

### 全体モデルの生成方法

1. metadataまたは明示指定から入力単位を解決。現段階の軸方向は入力Zで、XY中心だけを水平内面楕円から推定します。傾き・中心移動は診断し、自動補正しません。
2. 軸を通る全体縦断面を取得し、閉じた輪郭を左右半断面へ分離します。開曲線・分岐・複数輪郭・底部の軸端点が不明な断面は除外します。
3. 各半断面の外側底部軸端点→口唇最高点と、口唇最高点→内側底部軸端点を、それぞれ正規化弧長で対応付けます。元の輪郭を伸縮・回転する処理ではありません。
4. 対応点の半径・高さの座標中央値で内外面を集約し、底部の軸区間で接続します。生成輪郭の交差、回転メッシュの閉鎖性・面向き・正の材料体積を確認します。
5. 代表輪郭を分析軸まわりに回転し、元の座標位置・入力単位でPLYを出力します。全断面重ね合わせと方位別の最近傍輪郭距離も保存します。

これは全体代表モデルの第一段階です。弧長対応は部位の形態学的相同性を保証せず、口唇最高点は現在の口縁中央線の延長始点とは別定義です。相似・アフィン補正、部位抽出、中央線・法線器厚、個体間分類はこのCLIでは行いません。全体と部位別の比較情報を後から関連付ける基準として使用します。

### 実行例

```bash
# 従来のrequirementsだけで実行（新しい依存パッケージは不要）
python pottery_whole_model.py "0015Jinmen_small.ply" --unit m --output-dir whole_model --angle-step 5 --profile-spacing-mm 0.5

# 保存済みの解析軸中心を明示する場合（数値の単位は常にmm）
python pottery_whole_model.py "0015Jinmen_small.ply" --unit m --output-dir whole_model_fixed_axis --axis-xy-mm 138.148217 -107.444779

# 対応点の間隔を細かくした比較
python pottery_whole_model.py "0015Jinmen_small.ply" --unit m --output-dir whole_model_fine --profile-spacing-mm 0.25
```

`--unit auto|mm|cm|m`（既定auto）、`--angle-step`（既定5°、180°を整数分割）、`--start-angle`、`--axis-surface inner|outer`、`--axis-sections`（既定40）、`--axis-xy-mm X Y`、`--revolution-sections`（既定180）、`--max-memory-mb`（既定1024）、`--require-watertight`、`--no-plots` に対応します。プロフィール標本間隔は各輪郭の中央値弧長から対応点数を決める目安で、各断面の実際の点間隔は異なります。voxelピッチではありません。出力先は空フォルダを指定します。

非watertight入力を自動修復せず、各断面の閉鎖性を検証します。入力トポロジーはJSONに記録し、`--require-watertight` を指定すれば入力段階で厳格に停止します。予定半断面の80%未満しか有効でない場合も停止します。左右は同じ断面平面として一緒に除外します。受理された方位は等重みで、欠落角度を補完しません。

### 出力と実データでの確認

| ファイル | 内容 |
| --- | --- |
| `whole_representative.ply` | 内外面・口縁・底部を持つ全体代表メッシュ。元座標・入力単位 |
| `representative_profile.csv`／`representative_profile_xy.ply` | 全体代表断面。CSVはr,zのmm、PLYはX=半径・Y=入力高さ・Z=0、入力単位 |
| `profiles_all_overlay.png`／`profiles_all_overlay_xy.ply` | 元全半断面（灰）と代表輪郭（赤） |
| `angular_deviations.csv`／`.png` | 方位別の最近傍輪郭距離の平均・RMS・p95・最大。器厚誤差ではない |
| `profile_distribution.npz` | 方位別の対応済み内外輪郭、中央値、使用方位。mm |
| `section_qa.csv` | 予定角度ごとの採否と理由、受理輪郭の弧長 |
| `whole_model.json` | 入力SHA-256、単位、解析軸と診断、条件、有効数、トポロジー、寸法・制約 |

提供された `0015Jinmen_small(1).ply`（350,070頂点・700,000面、m）で生成しました。入力は重複頂点統合後に開放境界0本・非多様体エッジ1本。5°の72半断面中66本を採用し、40°／220°、85°／265°、160°／340°は複数輪郭のため除外しました。0.5mm設定の代表器高278.650mm、最大径227.546mm、メッシュ249,482頂点・498,960面。出力はwatertightかつ面向き整合です。0.25mm設定では最大径227.660mm、器壁材料体積の変化約0.0077%。これはこの個体での補間密度比較で、角度間隔依存性や器厚精度の保証ではありません。JSONの材料体積は土器の保持容量ではありません。

軸対称化で消える方位変異は元断面と偏差に残します。復元面とオリジナル面の識別は未導入です。把手、台脚、複数材料輪郭、底部軸上の穴などは現方式の対象外または除外対象です。次段階は部位の対応、全体構成特徴、個体間比較と、未使用断面・独立実測による評価です。


<a id="whole-vessel-validation"></a>

## 全体代表モデルの独立検証（WholeValidation 0.1.0-dev）

`pottery_whole_validation.py` はこの1ファイルだけで使える検証CLIです。MorphPotパッケージをimportせず、標準ライブラリ・NumPy・SciPyで保存済み出力を検証します。未使用方位のメッシュ抽出にはtrimesh・networkx、PNG生成にはMatplotlibが必要です。既存requirementsで実行できます。全体モデルの生成CLI、口縁部検証CLIとは別立てです。

```bash
# 保存済み全体モデルだけで中央値・平均・トリム平均を比較
python pottery_whole_validation.py whole_model --output-dir whole_validation --interval-mm 1

# 同じ元PLYから未使用方位を抽出（5°モデルなら既定は2.5°の開始オフセット）
python pottery_whole_validation.py whole_model --output-dir whole_validation_holdout --input-mesh "0015Jinmen_small.ply" --holdout-offset-deg 2.5 --interval-mm 1

# 集約方式・トリム率・検証密度を指定
python pottery_whole_validation.py whole_model --output-dir whole_validation_trim --methods median trimmed_mean --trim-fraction 0.1 --interval-mm 0.5 --no-plots
```

入力は `whole_model.json` と `profile_distribution.npz` を持つWholeModel出力フォルダです。元PLYを指定した場合、SHA-256が生成時と一致することを確認してmetadataの単位・解析軸を使います。未使用方位は元の試行グリッドと重複すると拒否します。既定の方位間隔は生成時と同じ、開始オフセットはその半分です。`--holdout-angle-step` で変更できます。閉じた単一の半断面のみ採用し、有効率80%未満は停止します。元メッシュや元モデルを変更せず、出力は別の空フォルダを指定します。

中央値・算術平均・トリム平均は同じ内外輪郭の正規化弧長対応点で再計算します。`--trim-fraction 0.1` は各座標分布の下側・上側それぞれ10%を除く指定で、断面単位の外れ値除外ではありません。座標ごとに除く断面は異なる場合があります。平均／トリム平均の曲線は比較記述子であり、交差・3D閉鎖性を検証した代替メッシュとしては出力しません。

距離は半径–高さ平面における点→輪郭線分の正確な最近傍距離で、内面・外面を分けて測ります。比較点を一定mm弧長で置き、平均・RMS・中央値・p95・最大を報告します。半径×弧長区間で回転面の面積を近似した重み付き平均・RMSも併記します。方位ごとの表にはモデル→個別輪郭の逆方向RMS・p95、双方向最大距離も記録します。CloudCompareの符号付きC2Mや器厚の測定ではありません。距離の向きと表面の指定を揃えて評価します。

| ファイル／フォルダ | 内容 |
| --- | --- |
| `distance_summary.csv` | 学習利用断面／未使用方位 × 集約方式 × 内面・外面・合算の距離統計 |
| `aggregation_comparison.csv` | 各方式の対応点座標、中央値との差r,zと距離 |
| `comparison_models.npz` | 3方式の内外輪郭、mm |
| `aggregation_models_overlay.png` | 3方式の輪郭重ね合わせ |
| `training/median/` 等 | 利用断面との測点距離CSV、方位別統計CSV、全断面重ね合わせPNG／PLY |
| `holdout/median/` 等 | 未使用方位との同じ比較出力。元PLY指定時のみ |
| `holdout_section_qa.csv` | 未使用方位の採否・理由 |
| `validation.json` | 条件・有効数・定義・限界 |

training評価は保存済みの補間輪郭、holdout評価は元メッシュの新しい交線を用います。未使用方位も同じ元メッシュを使い、軸は固定するため、独立した実物計測による精度検証や軸推定を含む総合検証ではありません。ランドマーク対応、部位別区分、軸感度、開始位相・角度間隔の一括実験は未実装です。最小距離の方式を自動的に最良と判定せず、指標ごとの変化を比較します。


<a id="whole-vessel-dimension-validation"></a>

## 全体モデルの寸法・器厚検証（WholeDimensionValidation 0.1.0-dev）

`pottery_whole_dimension_validation.py` を追加しました。同じフォルダの `pottery_whole_validation.py` を共通処理として使用します。この2ファイルと既存requirementsで動作し、morphpotパッケージは不要です。入力はWholeModelの出力フォルダです。元モデル・メッシュは変更しません。

```bash
python pottery_whole_dimension_validation.py \
  "0015Jinmen_whole_model_v0.1.0/whole_model" \
  --output-dir whole_dimension_validation --interval-mm 1

# 未使用角度も検証（元PLYのSHA-256を確認、単位と解析軸を継承）
python pottery_whole_dimension_validation.py \
  "0015Jinmen_whole_model_v0.1.0/whole_model" \
  --output-dir whole_dimension_validation_holdout \
  --input-mesh "../PotteryVolumeCalculator/0015Jinmen_small.ply" \
  --holdout-offset-deg 2.5 --interval-mm 1
```

| 測定項目 | 定義 |
| --- | --- |
| 口唇径 `lip_diameter_mm` | 最高口唇点の半径×2。各半断面の直径換算値で、対向半断面を結んだ実径や内側開口径ではない |
| 最大径 | 外面最大半径×2 |
| 最大径の高さ | 外面最小Zから最大半径点まで。最大点が複数ならそのZ範囲の中点を使用、範囲も保存 |
| 器高 | 内外面の最高Z−外面最小Z |
| 器厚 | 外側底部軸端点→口唇の外面接線に対する内向き法線と、内面の最初の正方向交点までの距離 |

器厚は各元断面自身の法線で測り、モデルの法線を元断面へ一律に適用しません。接線の推定幅は `--tangent-window-mm`（既定2mm）、外面の測点間隔は `--interval-mm`（既定1mm）。端点から弧長 `--endpoint-margin-mm`（既定2mm）以内は除外します。複数交点は最初を使用してフラグを付け、未取得は欠測とします。内面法線との入射角のcosが0.3未満となる浅い交差も器厚から除外し、交点までの生距離・cos値を別列に保持します。欠測区間を跨ぐ補間は行いません。

モデルと元断面の器厚比較は、底部軸端点から口唇までの正規化外面弧長で対応付けます。同じ絶対弧長や形態学的ランドマークでの比較ではありません。測点ごとに有効半断面数、平均・中央値・母標準偏差・最小・最大、モデル−実測中央値を保存します。分布は受理断面の等重みで、断面単位の外れ値除去はしません。中央値・平均・座標ごとの両側10%トリム平均の3モデルを比較します。

| 出力 | 内容 |
| --- | --- |
| `model_dimensions.csv` | 集約方式別のモデル寸法 |
| `source_dimensions.csv` | 方位別の元断面寸法 |
| `dimension_summary.csv` | 寸法ごとの断面間統計とモデル−中央値 |
| `METHOD_model_thickness.csv` | モデルの器厚測点、交点、欠測・複数交点フラグ |
| `source_thickness.csv` | 各元断面の器厚と測点座標、QA |
| `thickness_comparison.csv` | 正規化弧長対応による実測器厚分布とモデル差 |
| `training_dimensions.png`／`holdout_dimensions.png` | 方位別寸法とモデルの比較 |
| `training_thickness.png`／`holdout_thickness.png` | モデル器厚、実測中央値・最小〜最大の比較 |
| `dimension_validation.json` | 条件・定義・有効数・限界 |
| `holdout_section_qa.csv` | 元PLY指定時の未使用角度の採否・理由 |

出力先は元フォルダ外の空フォルダを指定します。`--no-plots`でCSVのみ生成できます。全数の最小〜最大は信頼区間ではありません。法線が遠方の内面へ達する場合があり、複数交点や極端な器厚は座標・フラグと合わせて確認します。真の内側開口径の同定、口唇そのものの器厚、部位ランドマーク、復元面の分離、軸感度は未対応です。


<a id="whole-thickness-diagnostics"></a>

## 口縁付近の器厚交点・対応の診断（WholeThicknessDiagnostics 0.1.0-dev）

`pottery_whole_thickness_diagnostics.py` は保存済みWholeModelから追加分だけを実行するCLIです。モデルを再生成せず、3集約モデルと各実測半断面で測定線を可視化します。既存の `morphpot.rim_standardization.make_midline` を使用するため、MorphPotリポジトリ全体と既存requirementsが必要です。

```bash
python pottery_whole_thickness_diagnostics.py \
  "0015Jinmen_whole_model_v0.1.0/whole_model" \
  --output-dir whole_thickness_diagnostics \
  --rim-length-mm 40 --pair-window-mm 8 --interval-mm 1

# 未使用角度も追加
python pottery_whole_thickness_diagnostics.py \
  "0015Jinmen_whole_model_v0.1.0/whole_model" \
  --output-dir whole_thickness_diagnostics_holdout \
  --input-mesh "../PotteryVolumeCalculator/0015Jinmen_small.ply" \
  --holdout-offset-deg 2.5 --rim-length-mm 40 --pair-window-mm 8
```

| 方法・色 | 内容 |
| --- | --- |
| outer_unrestricted・橙 | 既存の外面法線測定。離れた位置への生の測定線も保持 |
| outer_constrained・青 | 最高口唇点からの内外面弧長差が指定窓内の交点だけを採用 |
| rim_midline・緑 | 既存DTW・中央線延長始点アルゴリズムで局所中央線を抽出し、各中央線の法線で両壁を測定 |
| 灰 | 変形していない元の内外面断面 |

`--rim-length-mm` は既存中央線処理のmanual比較範囲を指定（既定40mm）。全体輪郭の外面最高口唇点からも同じ長さの区間を診断します。中央線の始点は延長交点であり最高口唇点とは異なるため、比較CSVでは外面上の口唇からの弧長へ対応させます。相似・アフィン変換をかけず、元のr,zで測ります。既存の標準化済み口縁出力をそのまま取り込む方式ではありません。

外面法線の制約は `abs(outer_lip_arc−inner_lip_arc) <= pair_window_mm`。中央線法線では、DTWで対応した各壁の位置から指定窓内にある交点だけを採用します。窓の既定8mmは検証用パラメータで、4・8・12mmなど別出力先で感度を確認してください。交点なし・浅い入射cos<0.3は欠測とし、中央線始点は厚さ観測から除外します。中央線のペア投影値は別列に残しますが、欠測厚さの代用にはしません。中央線処理が失敗した断面も法線結果は残し、理由をQAへ記録します。

| 出力 | 内容 |
| --- | --- |
| `models/median/`等 | 各モデルの測定CSV、断面＋交点線分PLY、局所重ね合わせと器厚比較PNG |
| `training/angle_XXX.XXX/` | 各実測半断面の同じ出力 |
| `holdout/angle_XXX.XXX/` | 元PLY指定時の未使用角度の同じ出力 |
| `training_all_rays_overlay.png`／`…_xy.ply` | 全実測断面と測定線の重ね合わせ（models・holdoutも同様） |
| `all_measurements.csv` | 全測点の方法・厚さ・生距離・内外交点座標・対応弧長差・入射cos・QA |
| `model_method_comparison.csv` | 共通外面弧長位置でのモデル測定3方式の比較。欠測区間を補間・外挿しない |
| `measurement_summary.csv` | 方法別の有効数・最大厚・対応窓超過数 |
| `midline_qa.csv` | 局所中央線の採否・既存始点診断 |
| `thickness_diagnostics.json` | 条件・定義・限界 |

PLYはX=半径、Y=入力高さ、Z=0で入力単位。CSV・PNGはmmです。PNGの線分には交点を取得したが器厚から除外した線も含み、CSVのstatusで判別します。出力先は元モデルフォルダ外の空フォルダを指定。`--no-plots`でもCSV・PLYは生成します。異なる法線方向の距離は同じ定義の厚さではなく、方法間差をすべて誤差と判断しません。

実データの中央値モデルでは、無制約24.346mmの線が外面口唇弧長2.685mmから内面22.935mmへ到達し、対応弧長差20.250mmでした。8mm窓ではこの交点を除外し、局所区間の最大有効値は制約外面法線6.724mm、中央線法線6.499mm。これは真の器厚の確定ではなく、遠い内面との交差が急増に寄与していたことの診断です。


<a id="whole-shape-exploration"></a>

## 器形条件に基づく全体モデル生成法の探索（0.1.0-exploration）

`pottery_whole_shape_exploration.py`は、正規化済みPLYから全長弧長対応、頸部で分割する区間別対応、区間別対応＋底部固定の限定変形の3方式を生成・比較する開発検証用CLIです。現在のカテゴリ`necked_jar`は頸部のくびれを持つ器形の暫定的な幾何条件であり、考古学的な器種・型式の自動同定ではありません。

```bash
python -m pip install -r requirements.txt
python pottery_whole_shape_exploration.py "/path/to/normalized.ply" --category necked_jar --output-dir exploration/individual_01
```

入力PLYの付属metadataから単位を読み取り、姿勢Transformは再適用しません。metadataがない場合は`--unit m`等を明示してください。出力先は新規または空ディレクトリを指定します。条件に不適合・未確定なら診断を保存して非ゼロ終了し、モデル生成を中止します。

- `analysis/global_arc/`：全長弧長対応によるモデル。
- `analysis/neck_segmented/`：頸部前後を別々に対応づけたモデル。座標変形なし。
- `analysis/neck_segmented_warp/`：同じ空間変位場を内外面に適用したモデル。底部保護域を固定し、変位・せん断・高さひずみを制限。
- `validation/`：全断面とモデルの重ね合わせ、距離・寸法・器厚のCSV、3方式の比較画像。
- `exploration.json`：入力ハッシュ、単位、軸推定、ランドマーク支持率、変形の適用・拒否方向。

制限を超える変形は拒否し、その方向は原形状のまま集計に残します。変形後のモデルを全方向の歪み除去済みと扱わないでください。内外面共通の変形でも器厚保存は保証しません。

生成済み結果から部位別比較と器厚診断を追加する場合：

```bash
python summarize_whole_shape_exploration.py exploration/individual_01 exploration/individual_02 --output-dir exploration/research
python plot_whole_thickness_diagnostics.py exploration/individual_01 exploration/individual_02
```

全体モデル探索の検証出力には、考古学者向けの説明図を追加しました。`validation/registration_comparison/RegistrationComparison.md`に、特徴点で位置合わせしない原位置断面（左）と、保存済みの頸部・口唇の限定変形を適用した断面（右）を同じ縮尺で掲載します。全体図・頸部口縁拡大図・統計グラフをtrainingとholdoutに分けて出力します。

保存済みの探索結果から追加分だけを生成できます（メッシュの再解析は不要）：

```bash
python plot_whole_registration_comparison.py exploration/individual_01 exploration/individual_02

# 出力を別フォルダへまとめる場合：その下に個体名のサブフォルダを作成
python plot_whole_registration_comparison.py exploration/individual_01 exploration/individual_02 \
  --output-dir registration_figures --interval-mm 1 --height-bin-mm 5
```

| 追加出力 | 内容 |
| --- | --- |
| `training_registration_comparison.png`／`holdout_…` | 原位置と位置合わせ後の全断面を左右に比較。共通の代表モデル・軸範囲 |
| `*_registration_comparison_detail.png` | 同じ比較の頸部・口縁拡大 |
| `*_registration_statistics.png` | 平均・母SD・RMS・p95・最大、方位別RMS、原位置高さ別平均±SD、距離累積分布 |
| `*_paired_station_distances.csv` | 同じ元弧長測点の変形前後座標・距離・部位・変形状態 |
| `registration_distance_summary.csv` | 内外面・合算、全体・底部・胴部・頸部〜口唇別の距離統計 |
| `*_angular_statistics.csv`／`*_height_statistics.csv` | グラフの方位別・原位置高さ別統計 |
| `RegistrationComparison.md`／`registration_comparison.json` | 説明レポートと条件・定義 |

対象は`exploration.json`を持つ器形探索出力です。従来のWholeModel／WholeValidationだけでは特徴点変形の記録がないため、この比較CLIの入力にはなりません。原位置は各断面を共通の半径–高さ平面へ投影したr,zで、追加の特徴点変形を行っていない状態です。左右で同じモデル・同じ原弧長測点を使い、底部を固定、変形拒否方向は原形状で残します。統計は符号なしの測点→同表面モデル線分距離で、器厚ではありません。変形後は再標本化しないため、既存の距離CSVとは小幅な差が生じ得ます。変形後の散らばり減少を、真の形態への精度改善とは解釈しません。

集計先`research`には部位別RMS・寸法差のCSVと`regional_and_dimensions.png`を出力します。器厚診断は各個体の`validation/model_thickness_diagnostic.png`へ出力します。解析一次出力、個別検証、解釈に用いる派生図表を分けて整理する構成です。

0015Jinmen・MK18-Hajikiのbase姿勢PLYで3方式ずつ生成済み。今回の比較では従来方式を置き換える明確な改善は確認できませんでした。器厚測定は既存の無制約法線交点を使う診断であり、口唇付近に遠方交点が含まれます。詳細な条件・実行例・探索判断は[WholeShapeExploration.md](WholeShapeExploration.md)を参照してください。
