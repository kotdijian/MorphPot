# PotteryRadialSections v0.4.0

> Historical instructions; use the [current MorphPot capacity guide](../README.md) for maintained commands.

土器3Dメッシュ（PLY / OBJ）から、**水平断面の楕円中心で回転軸を推定し、放射状の縦断面を抽出し、その内面プロファイルから容量を計算する**Python CLIプログラムです。

v0.4.0では容量計算を追加し、次の4方式を**個別または複数同時に選択**できます。

1. `single`：指定した1つの縦断面を回転体化
2. `optimized`：複数方向の縦断面をrobustに合成した最適回転体
3. `angular`：軸対称を仮定せず、方向別半径を角度方向に積分
4. `ellipse`：水平断面のinner ellipse面積をZ方向に積分する独立QA

既定の `--volume-mode all` では4方式すべてを計算します。

**Pythonやコンソール、Powershellなどの使用方法について[MorphPot容量ガイド](../README.md)をご参照ください**

---

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

v0.4.0の容量計算は、入力メッシュ中に実際の内面が存在する場合を対象とします。

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

1枚の全断面には軸の両側に2つの内面半径プロファイルがあります。v0.4.0では各高さで、

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

この方式は、楕円化・円形化されない**器形の非対称性を容量へ直接反映**できるため、v0.4.0では主要な計算法の一つとして扱います。

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

## 11. 0015Jinmen_small.plyでのv0.4.0検証

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
