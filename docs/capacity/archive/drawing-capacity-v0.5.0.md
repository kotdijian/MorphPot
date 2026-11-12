# PotteryDrawingCapacity v0.5.0 と関連する断面計算

> Historical instructions; use the [current MorphPot capacity guide](../README.md) for maintained commands.

土器3Dメッシュからの回転軸・縦断面・容量計算に加えて、v0.5.0では**2D実測図・断面図から `single` 法で容量を計算する Drawing-Single GUI**を追加しました。

このパッケージには2つの入口があります。

```text
pottery_radial_sections.py   3D PLY/OBJ 用
pottery_drawing_capacity.py  2D 図面用 GUI
```

**Pythonやコンソール、Powershellなどの使用方法について[MorphPot容量ガイド](../README.md)をご参照ください**

---

# A. 2D図面からの Drawing-Single 容量計算

## A-1. 対応入力

優先する利用形態に合わせ、次をサポートします。

- ラスタ画像: `PNG / JPEG / TIFF / BMP / WEBP`
- `SVG`: Adobe Illustrator等から書き出したSVGを含む
- `PDF`: 指定ページを画像化して利用

SVG/PDFについても、v0.5.0ではベクトルpathを直接容量計算へ渡さず、**画面表示用にレンダリングした画像をdigitize**します。したがってSVG/PDFの内部単位やDPIを縮尺として信用せず、すべての形式でGUIによる実寸校正を行います。

### 前提

- 左右両側の**内面縁（inner profile）**が図化されている
- 図面の水平・垂直がすでに定まっている
- 回転軸は垂直である
- 単純な1価関数 `r(z)` として扱える器形を初期対象とする

図面線の完全自動認識はv0.5.0には入れていません。研究用の初版として、**縮尺・回転軸・左右内面を人が明示的に指定**し、その入力をCSV/JSONへ保存する方式を採用しています。

## A-2. インストール

```bash
python3 -m pip install -r requirements.txt
```

macOS/Windowsの標準的なPython.org配布版ではGUIにTkinterを利用できます。

## A-3. 起動

図面を起動時に指定する場合:

```bash
python3 pottery_drawing_capacity.py pottery_section.png
```

ファイル選択から開始する場合:

```bash
python3 pottery_drawing_capacity.py
```

PDFの場合は複数ページならGUIでページ番号を指定できます。コマンドから明示することもできます。

```bash
python3 pottery_drawing_capacity.py report.pdf --page 12
```

## A-4. GUIの操作順

1. **図面を開く**
2. **1 縮尺 (2点)** を押す
3. スケールバーなど既知距離の両端を2点クリック
4. その実寸を `mm` で入力
5. **2 回転軸X** を押し、垂直な回転軸上を1点クリック
   - 水平・垂直が既定なのでX座標のみ使用します
6. **3 左内面** を押し、内腔底部から口縁へ順にクリック
7. `Enter` または右クリックで左内面を確定
8. **4 右内面** も同様に指定
9. **計算・保存**

`Ctrl+Z`または「1点戻す」で直前のdigitize点を戻せます。`+ / -`またはマウスホイールで拡大縮小できます。

## A-5. 縮尺

2点間の画像距離を `D_px`、ユーザー入力実寸を `D_mm` とすると、

```text
mm_per_pixel = D_mm / D_px
```

を採用します。図面・SVG・PDFの内部DPIや用紙サイズは容量計算の実寸には使いません。

## A-6. Drawing-Single 容量計算

回転軸から左右内面までの半径を

```text
r_left(z)
r_right(z)
```

とし、3Dモードの `single` と同じ二側面定義で

```text
r_equivalent(z) = sqrt((r_left(z)^2 + r_right(z)^2) / 2)
A(z) = pi/2 * (r_left(z)^2 + r_right(z)^2)
```

を用います。

左右プロファイルの共通Z範囲を計算領域とし、上端は**左右のうち低い口縁高**、すなわちspill levelとします。プロファイルを既定 `0.5 mm` 間隔でZ方向に補間し、断面積を台形積分します。

```bash
--z-step-mm 0.5
```

結果には合成容量だけでなく、左側のみ・右側のみをそれぞれ回転した診断容量も保存します。

## A-7. 出力

既定では入力ファイルと同じ場所に、

```text
<入力名>_DrawingCapacity/
```

を作成します。

```text
source_render.png
  GUIで実際にdigitizeした表示画像

drawing_profile_raw.csv
  左右のクリック点、画像座標、軸からの距離、Z座標

drawing_profile_resampled.csv
  容量計算に使用したZごとの左右半径・等価半径・断面積

drawing_volume_summary.csv
drawing_volume_summary.json
  容量、左右個別容量、spill level等

drawing_metadata.json
  元ファイルSHA256、縮尺指定2点、mm/px、回転軸X、前提条件等

drawing_qc_overlay.png
  元図上に縮尺、回転軸、左右digitize線、spill levelを重ねた監査図

drawing_profile_plot.png
  左右内面を回転軸基準のr-z図として表示
```

SVG/PDFを入力した場合も、`source_render.png`を保存するため、どのレンダリング画像上で測定したかを後から確認できます。

## A-8. SELF TEST

半径50 mm、高さ100 mmの円筒を模擬した内部テスト:

```bash
python3 pottery_drawing_capacity.py --self-test
```

理論値は

```text
0.785398163 L
```

で、v0.5.0のテストでは数値誤差範囲で一致します。

---

# B. 3Dメッシュからの回転軸・縦断面・容量計算

この節はv0.5.0時点の断面計算の記録です。同梱の `pottery_radial_sections.py`（v0.4.0）の実行手順は[放射断面ガイド](radial-sections-v0.4.0.md)を参照してください。

## 1. 前提となる3Dモデル

- 完形またはほぼ完形の土器
- Z軸が土器の上下方向
- PLYまたはOBJ
- 入力座標単位は `mm` / `cm` / `m`
- 容量計算には、現段階では**内面が3Dメッシュとして取得されていること**が必要
- 回転軸推定は `inner / outer / both` を選択可能

PLYはtexture UVによる頂点再分割を避けるため、`fix_texture=False`でgeometry topologyを保持して読み込みます。

---

## 2. 回転軸推定

### 水平断面

既定では器高を20等分し、各区間中央、すなわち相対器高

```text
2.5%, 7.5%, 12.5%, ... , 97.5%
```

でXY水平断面を取得します。

```bash
--z-sections 20
```

絶対間隔を使う場合：

```bash
--z-step-mm 20
```

`--z-sections`と`--z-step-mm`は排他的です。

### inner / outer contour

水平断面のmesh-plane intersectionを閉曲線化し、面積と包含関係からouter / inner contourを識別します。輪郭は弧長方向に等間隔再サンプリングしてrobust ellipse fittingを行います。

既定：

```bash
--axis-surface inner
--contour-spacing-mm 0.5
```

`outer`では外面楕円中心を軸推定に使います。`both`では内外両方をfitしますが最終軸はinner中心から求めます。

### outerモードと器厚

`outer`モードの**回転軸推定には器厚を使用しません**。

外面しか利用できない資料について、後続の容量復元で単一器厚またはZ位置付き器厚プロファイルを入力して内面を推定する機能は、別段階として追加する想定です。

v0.5.0の容量計算は、入力メッシュ中に実際の内面が存在する場合を対象とします。

### 中心外れ値と最終軸

各水平断面のellipse centerを `(Xi, Yi)` とし、中心分布からMADで2次元外れ値を除外します。既定閾値は、中心中央値からの距離について

```text
median(distance) + 3 × 1.4826 × MAD(distance)
```

です。

残った中心の算術平均を最終軸とします。

```text
Xaxis = mean(Xi)
Yaxis = mean(Yi)
rotation axis = (Xaxis, Yaxis, Z)
```

Zに伴う中心ドリフトと推定傾斜角もQAとして保存しますが、モデルの自動傾き補正は行いません。

---

## 3. 縦断面

推定回転軸を通る縦断面を一定角度間隔で取得します。既定30°では、6枚の全断面と12方向の放射半断面を作ります。

```text
full planes:
0°/180°, 30°/210°, 60°/240°,
90°/270°, 120°/300°, 150°/330°

radial directions:
0°, 30°, 60°, ... , 330°
```

CSV、edge PLY、point-cloud PLYを保存します。PLYは**入力メッシュと同じXYZ座標系・単位**なのでCloudCompare等で元モデルへ重ねられます。

---

# 4. 容量計算モード

## 4.1 `single` — 指定した1縦断面の回転体

```bash
--volume-mode single --single-angle 0
```

`--single-angle 0`は、0°/180°を通る**1枚の全縦断面平面**を選択します。

1枚の全断面には軸の両側に2つの内面半径プロファイルがあります。v0.5.0では各高さで、

```text
r_equivalent = sqrt((r_angle² + r_opposite²) / 2)
```

を等価半径とし、

```text
A(z) = π r_equivalent(z)²
```

として回転体容量を計算します。

これは左右2側の回転体断面積を平均したものと等価です。CSVには両側を単独で回転した場合の容量も診断値として保存します。

特定方向による容量差を確認する比較用・従来的回転体法として位置づけます。

---

## 4.2 `optimized` — 複数縦断面からの最適回転体

```bash
--volume-mode optimized
```

30°間隔なら各Zで12方向の内面半径を取得します。

局所的な異常半径の影響を抑えるため、各高さでmedian/MADを用いたrobust選択を行い、残った半径の平均から

```text
r_opt(z)
```

を求めます。

```text
A(z) = π r_opt(z)²
```

として積分します。

これは**資料全体を最も代表する軸対称モデル**を作る方式です。

---

## 4.3 `angular` — 角度方向積分

```bash
--volume-mode angular
```

軸対称を仮定せず、各方向の実測内面半径 `r(z, θ)` を使います。

水平断面積を

```text
A(z) = 1/2 ∫ r(z,θ)² dθ
```

として数値積分し、さらにZ方向に積分します。

30°間隔なら12方向を使います。

この方式は、楕円化・円形化されない**器形の非対称性を容量へ直接反映**できるため、v0.5.0では主要な計算法の一つとして扱います。

一部方向の半径が欠ける高さでは、既定で75%以上の方向が有効な場合に限り、角度方向の周期的線形補間を行います。

```bash
--min-angular-valid-fraction 0.75
```

---

## 4.4 `ellipse` — 水平楕円積分

```bash
--volume-mode ellipse
```

回転軸推定時に得られたinner ellipseの長半径 `a(z)`、短半径 `b(z)`から、

```text
A(z) = π a(z) b(z)
```

を求め、Z方向に積分します。

縦断面半径とは独立した水平断面解析なので、**QA・相互検証用**として有効です。

現在は軸推定用水平断面（既定20断面）のellipseを補間して用います。したがって`angular`より粗い独立検証法という位置づけです。

---

## 5. spill level（容量上限）

容量は口縁より上まで積分しません。

各放射方向について、内面と外面の2交点が連続して得られる最長Z区間を内腔プロファイルとみなし、その上端を各方向のrim heightとして推定します。

`optimized / angular / ellipse`では、正常な放射プロファイルのうち**最も低いrim height**をfirst-spill levelの近似値として使用します。

```text
spill = minimum reliable radial rim height
```

したがって波状口縁などでは、最も低い口縁位置が容量上限になります。

この推定精度はZサンプリング間隔に依存し、概ね `± volume_z_step_mm / 2` の離散化幅があります。

既定：

```bash
--volume-z-step-mm 0.5
```

---

# 6. モードの選択

### 4方式すべて — 既定

```bash
python3 pottery_radial_sections.py pottery.ply --unit m
```

または明示して、

```bash
python3 pottery_radial_sections.py pottery.ply --unit m \
  --volume-mode all
```

### 1方式だけ

```bash
--volume-mode single
```

```bash
--volume-mode optimized
```

```bash
--volume-mode angular
```

```bash
--volume-mode ellipse
```

### 複数方式だけ

例えば、optimizedとangularだけを比較する場合：

```bash
python3 pottery_radial_sections.py pottery.ply --unit m \
  --volume-mode optimized angular
```

singleとellipseだけなら：

```bash
python3 pottery_radial_sections.py pottery.ply --unit m \
  --volume-mode single ellipse \
  --single-angle 90
```

### 容量計算をしない

回転軸推定と縦断面抽出だけを行う場合：

```bash
--volume-mode none
```

---

## 7. 主な実行例

### 標準

```bash
python3 pottery_radial_sections.py 0015Jinmen_small.ply --unit m
```

### 15°間隔でangular積分

```bash
python3 pottery_radial_sections.py 0015Jinmen_small.ply \
  --unit m \
  --angle-step 15 \
  --volume-mode angular
```

30°→15°→10°と細かくして、角度方向の収束を確認できます。

### singleで60°/240°断面を使用

```bash
python3 pottery_radial_sections.py 0015Jinmen_small.ply \
  --unit m \
  --volume-mode single \
  --single-angle 60
```

### Z積分を1 mm間隔にする

```bash
--volume-z-step-mm 1.0
```

### outer軸推定＋容量計算

```bash
python3 pottery_radial_sections.py pottery.ply \
  --unit m \
  --axis-surface outer \
  --volume-mode angular optimized
```

この場合も容量計算自体には入力メッシュの内面を使用します。器厚からの内面推定はまだ行いません。

---

## 8. 出力構成

```text
<model>_RadialSections_30deg/
├── metadata.json
├── sections_summary.csv
│
├── axis_estimation/
│   ├── horizontal_sections.csv
│   ├── axis_summary.csv
│   ├── all_horizontal_section_points.ply
│   ├── all_fitted_ellipse_points.ply
│   ├── section_centers_points.ply
│   ├── rotation_axis_edges.ply
│   └── horizontal_sections/
│
├── full_sections/
│   ├── section_000_180.csv
│   ├── section_000_180_edges.ply
│   ├── section_000_180_points.ply
│   └── ...
│
├── radial_half_sections/
│   ├── ray_000.csv
│   ├── ray_000_edges.ply
│   ├── ray_000_points.ply
│   └── ...
│
├── volume/
│   ├── volume_summary.csv
│   ├── volume_summary.json
│   ├── radial_inner_profiles.csv
│   ├── radial_inner_profile_points.ply
│   ├── single_section_000_180_profile.csv      # single選択時
│   ├── optimized_rotation_profile.csv          # optimized選択時
│   ├── angular_area_profile.csv                # angular選択時
│   ├── horizontal_ellipse_area_samples.csv     # ellipse選択時
│   ├── volume_method_comparison.png
│   └── radial_inner_profiles.png
│
└── visualization/
    ├── axis_validation_xz.png
    ├── axis_validation_yz.png
    ├── axis_centers_xy.png
    ├── horizontal_sections_oblique.png
    ├── horizontal_sections_table.png
    ├── axis_summary_table.png
    ├── full_sections_oblique.png
    └── radial_half_sections_oblique.png
```

---

## 9. `volume_summary.csv`

選択した方式ごとに1行を保存します。主な項目は、

- `method`
- `status`
- `volume_mm3`
- `volume_ml`
- `volume_l`
- `bottom_z_mm`
- `spill_z_mm`
- `z_step_mm`
- `valid_z_cells`
- `angle_deg / opposite_angle_deg`（single）
- `side_a_volume_l / side_b_volume_l`（singleの診断値）
- `radial_directions`
- `angular_step_deg`
- `min_valid_rays`
- `horizontal_ellipse_sections`
- `note`

です。

---

## 10. 3D監査用PLY

次のPLYはすべて入力メッシュと同じXYZ座標系・単位です。

- `axis_estimation/all_horizontal_section_points.ply`
- `axis_estimation/all_fitted_ellipse_points.ply`
- `axis_estimation/section_centers_points.ply`
- `axis_estimation/rotation_axis_edges.ply`
- `full_sections/*_points.ply / *_edges.ply`
- `radial_half_sections/*_points.ply / *_edges.ply`
- `volume/radial_inner_profile_points.ply`

`radial_inner_profile_points.ply`は、容量計算に実際に使われた方向別inner profileを3D空間上に示します。元メッシュへ重ねることで、どの位置を内面として取得したか監査できます。

---

## 11. 0015Jinmen_small.plyでのv0.5.0検証

既定設定：

```bash
python3 pottery_radial_sections.py 0015Jinmen_small.ply --unit m
```

回転軸：

```text
horizontal sections attempted : 20
valid inner ellipse centers    : 19
used for final axis            : 19
axis X                         : 138.092372 mm
axis Y                         : -107.491406 mm
RMS center offset              : 2.3259 mm
diagnostic centerline tilt     : 1.2695°
```

容量（30°、Z step 0.5 mm）：

```text
single    0°/180° : 6.238819 L
optimized          : 5.986216 L
angular            : 6.017601 L
ellipse            : 6.000679 L
```

`single`の0°/180°では、各側を単独で回転した診断値は約6.315 L、6.272 Lでした。

3つの複数断面方式は約5.99–6.02 Lに集中し、単一断面方式は方向依存性によってより大きな値を示しました。この差自体が、単一縦断面による回転体近似の不確実性を評価する情報になります。

---

## 12. 現段階の制約

- 回転軸はZ平行で、推定されたcenterline tiltを自動補正しない
- 容量計算は入力メッシュにinner surfaceが存在する完形・ほぼ完形資料向け
- `outer`から器厚情報でinner surfaceを推定する機能は未実装
- inner profileは、各Zで正の放射方向に2つ以上のmesh交点がある場合、最小正半径をinnerと判定する
- 複雑な内面突起、二重口縁、注口、把手等では追加分類が必要になる可能性がある
- `ellipse`は既定20水平断面に基づくため、主計算法より粗いQA法
- spill heightは縦断面の角度間隔とZ stepに依存する

研究用の正式値として使用する前に、`--angle-step`と`--volume-z-step-mm`を変えた収束確認を推奨します。
