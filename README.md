# MorphPot

**姿勢正規化済みの土器3Dモデルを対象とする、容量・器軸・断面・表面形態解析プロジェクトです。**

基本入力は[ArtefactsOrthoMaker](https://github.com/kotdijian/ArtefactsOrthoMaker)が出力する三角形メッシュPLYです。既に適用された姿勢Transformを再適用しません。

現在は容量計算、水平・放射断面、断面容量、実測図容量、3Dガイド、破片再構成の既存モジュールを継承した開発版です。土器表面痕跡と接合境界候補の解析も収録しています。各処理はCLIまたは継承元GUIから独立して実行でき、MorphPot全体の統合GUIは今後開発します。

Surface Enhancement Labは共通手法の継承用としてexperimental/に収録しています。土器の投影方向・単位に適合させた統合は今後の予定です。

石器の平面形態・連続断面は[LithMorph](https://github.com/kotdijian/LithMorph)で開発します。

実装範囲、モジュール構成、開発計画、検証結果は[DevelopmentReport](DevelopmentReport.md)を参照してください。移管元は[SOURCE_PROVENANCE.json](SOURCE_PROVENANCE.json)に記録しています。

本リポジトリは開発版であり、全資料での精度検証を完了した正式安定版ではありません。

## PotteryRadialSections 操作ガイド

`pottery_radial_sections.py`（v0.8.1）は、水平断面の楕円中心から器軸のXY位置を推定し、その軸を通る放射状の縦断面を抽出します。断面からの容量計算と、右側基準の口縁部中央線の抽出・個体内標準化を実装しています。現在はCLIで操作します。

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
| `--no-section-overlay` | 無効 | XY重ね合わせPLYと中央値ポリラインの自動出力を省略 |

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

### 7. 口縁部中央線を右側へ揃え、個体内標準モデルを作る（v0.8.1）

目的は、一個体の口縁部周辺の半断面を比較し、位置・向き・倍率・限定的な歪みを除いた標準モデルを作ることです。左右の半断面を右側へ揃え、口唇から左への進行の停滞・反転位置＋バッファーまでに限定します。胴部全体は位置合わせに使いません。

v0.7.0の `median_full_section_xy.ply` は内外面を含む全断面の中央値輪郭です。v0.8.1の**中央線**は、口縁部周辺の内外面を輪郭順序を保った局所DTWで対応付け、対応点の中間を取った線です。高さの単調性を要求しません。対応は正規化弧長位置の25%帯に制限し、各壁を最大512点で標本化します。厳密なmedial axisや法線器厚の中間線ではありません。口唇位置、対応、比較終端をQAとPLYで確認してください。

通常の `sections` 実行で、次の処理を自動で行います。

1. 全縦断面の内外面輪郭から左右の半断面を取り出す。
2. 左半断面を器軸に対して反転し、右半断面と同じ向きへ揃える。Xは右向きの半径、Yは元モデルのZ高さ、投影Zは0。
3. 口唇を始点として中央線を追跡し、左への進行が持続した後、その進行が停滞・反転する位置を検出する。上下方向の傾き・増減は条件にしない。
4. 変化点の先に局所的な**対応内外面間距離の2倍**の弧長バッファーを加え、そこで切る。固定mm指定も可能。
5. 各比較区間を正規化弧長で同じ129点へ再標本化する。
6. 採用された右側中央線群の座標中央値を固定参照とし、相似変換と制約付きアフィン変換でそれぞれ位置合わせする。
7. 位置合わせ後の座標中央値を標準中央線として出力し、点ごとの残差と変換量を保存する。

右側の口唇を頂点0として、輪郭を下方へ追跡します。口唇検出は前述の幾何学的最高点に基づき、考古学的同定の保証はありません。水平進行判定用の平滑化だけを行い、保存する中央線の形状にはその平滑化を適用しません。

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

アフィン制約の初期値は実資料による校正前の暫定値です。相似変換に対して反転を許さず、アフィン変換でも正の行列式を保持します。ほぼ直線の中央線ではアフィン変形を十分に推定できないため、相似変換へ戻しQAに記録します。参照は固定された右側中央線群の中央値であり、反復的な一般化Procrustes解析ではありません。

出力先は `rim_standardization/` です。

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
