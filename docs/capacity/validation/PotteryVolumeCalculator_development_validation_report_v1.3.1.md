# PotteryVolumeCalculator v1.3.1 開発・検証レポート

**3D土器モデルに対する最大保持容量のvoxel推定**\
**対象資料：0015Jinmen_small.ply**\
**作成日：2026年8月16日**

---

## 1. 本レポートの目的

PotteryVolumeCalculatorは、PLY / OBJ形式の土器3Dメッシュから、**液体が最初に外部へ流出する直前まで保持できる最大内容量**を推定するPython CLIプログラムである。

READMEは利用者向けの操作説明を中心とするのに対し、本レポートでは、次の事項を技術的に整理する。

- 容量計算アルゴリズムを構築した経緯
- 現行アルゴリズムの処理手順と判定原理
- QA・監査設計
- 実資料`0015Jinmen_small.ply`を用いた2.0 / 1.0 / 0.5 mm pitch検証
- 検証中に判明したtexture UV seam問題
- v1.3.1での`archaeological/raw_*`出力停止
- 大量資料の一括処理へ進む前の課題

v1.3.1は、v1.3.0で確立した容量計算ロジックを変更せず、**考古学的接合境界として意味づけできないpre-weld raw boundary派生出力を停止する暫定修正版**である。

---

## 2. 検証用サンプルデータ

本レポートで使用した3Dモデルの出典は次のとおりである。

> **人面墨書土器（奈良時代）**\
> 静岡県磐田市御殿・二之宮遺跡出土／磐田市埋蔵文化財センター所蔵\
> レガシズ3Dよりダウンロード\
> https://lega-shizu.com/legashizu3d/archives/data/117

レガシズ3D掲載ページでは、奈良時代の人面墨書土器であり、御殿・二之宮遺跡出土、磐田市埋蔵文化財センター所蔵とされ、モデル識別情報としてLS0015が記載されている。

今回の入力ファイルは`0015Jinmen_small.ply`で、容量計算時には**+Zを上方向**とする姿勢に配置した。

### 2.1 監査対象ファイル

| ファイル | 用途 | SHA-256 |
|---|---|---|
| `0015Jinmen_small.ply` | 入力3Dモデル | `d98e8369b293fcc002303088bc3f48103be20d4fe5d10c41311af015c4a7785d` |
| `vessel_voxel_volume(20260816-021426).py` | v1.3.0実検証コード | `71a00ded854033ffbba00c199dae0564c3d0d9d6b75cdd720d24a35cfb0d0a7a` |
| `environment.txt` | 実行環境記録 | `0b7c4adb5dd29a0581efad93fa40ed24fd28c6e34ec51e3e6af71f908715e4d5` |
| `pip_freeze.txt` | Python package記録 | `1512b99227d973233e1e703440d730d55fcfa4422d16c3e79a602bbf5221760a` |
| `0015Jinmen_small_PotteryVolume_v1.zip` | validation出力一式 | `393a24acc6540cf32e50eca26628ad978f7301df4d1c982b1e070bdd0f4801b4` |
| `vessel_voxel_volume.py` v1.3.1 | raw出力停止版 | `3ac55b3ef5161e7d1478667861ad137361e7da78ae7be21c4aeb7108a8e381c3` |

---

## 3. 入力モデルの概要

元PLYは350,070 vertices、700,000 triangular facesを持ち、vertex RGBとface-corner texture coordinateを含む。座標値はm単位であり、外形寸法はおよそ次のとおりである。

| 軸 | extent |
|---|---:|
| X | 231.424 mm |
| Y | 229.222 mm |
| Z | 283.793 mm |

Zが器高方向となる直立姿勢である。

[元モデルの図版（PDF p. 3）](https://github.com/kotdijian/PotteryVolumeCalculator/blob/93fbdbad24f22c38423d9e1b4780622b2ef8ccab/docs/validation/PotteryVolumeCalculator_development_validation_report_v1.3.1.pdf#page=3)

**図1　`0015Jinmen_small.ply`の3D形状。** PLY内のvertex RGBを用いた表示。外部texture画像は用いていない。

[正面の図版（PDF p. 4）](https://github.com/kotdijian/PotteryVolumeCalculator/blob/93fbdbad24f22c38423d9e1b4780622b2ef8ccab/docs/validation/PotteryVolumeCalculator_development_validation_report_v1.3.1.pdf#page=4)

**図2　正面方向からの正投影。**

---

## 4. 容量の定義

本プログラムでは容量を、

> **土器を+Z方向に直立させ、液体の自由表面を上昇させたとき、液体が初めて外部へ流出する直前まで保持できる最大体積**

と定義する。

この定義では、口縁全体が水平である必要はない。波状口縁や最低口縁部がある場合には、内部free-spaceと外部が最初に連結する高さがspill levelとなる。

また、底部から連続した残存資料では、本来の完形時容量ではなく、**現存形状で保持可能な上限容量**を求めることになる。

---

## 5. アルゴリズム構築の経緯

### 5.1 水平cap方式の問題

初期案では、土器口部を水平面で閉じて閉曲面を作り、その内部体積を求める方法を検討した。しかし考古資料では口縁が水平とは限らず、波状口縁や局所欠損も一般的である。この方法ではcapの高さを人為的に決める必要があり、最大保持容量の定義と一致しない。

### 5.2 whole-mesh surface voxelization

そこで、内面だけを明示的に抽出するのではなく、**土器全体の三角形surfaceをvoxel化し、surface voxelを液体が通過できない障壁として扱う**方式へ移行した。

計算対象はmeshのsolid volumeではなく、surface voxelに囲まれたfree voxelのうち、液体が連続して占有できる領域である。

### 5.3 capを廃止しspill-level探索へ

口部を人工的に閉じる代わりに、水面高さ`Z = k`以下だけを移動可能なfree-spaceとして内部seedからflood fillを行い、外部へ初めて到達する高さを探索する方式へ変更した。

これにより非水平口縁を特別扱いする必要がなくなり、「実際に最初に溢れる高さ」を計算原理に直接組み込める。

### 5.4 exact-coordinate weld

復元3Dモデルでは、同一位置にあるvertexがtopology上分離している場合がある。容量計算用meshでは、

- unreferenced vertexを除去
- **XYZが完全一致するvertexのみ統合**
- tolerance = 0
- 近接vertexは融合しない
- vertex位置を移動しない
- faceを追加・削除しない

という限定的なcleanupを採用した。

### 5.5 cavity seed誤認への対策

surface voxelizationでは、内面と外面の間の本来は粘土である領域が、surface間のfree voxelとして見える場合がある。この領域を内容空間と誤認しないため、seed候補に3条件を課した。

1. XY断面内で外周に接しないenclosed free-spaceである
2. 全高さを許せば外部へ到達できる - 口部を通じて最終的に外へ出られる
3. seed高さ以下に制限すると外部へ到達しない - seed位置ではまだ漏れていない

この条件により、器壁内部の閉じたfree-spaceをseedから除外する。

### 5.6 Z軸姿勢

液面上昇はZ方向に定義されるため、入力モデルの姿勢は計算結果に直接影響する。検証過程では長軸がZでないモデルで異常に低いspillが得られ、姿勢修正後にspillが口縁近くへ移動した。このため、**+Z = 上方向**は入力条件として明示した。

### 5.7 multi-pitch validation

surface voxelには有限の厚さがあるため、容量値はpitchに依存する。v1.3.0以降では詳細検証用に、

- 2.0 mm
- 1.0 mm
- 0.5 mm

の3解像度を独立した全工程として一括実行するvalidation modeを実装した。

---

## 6. 現行の容量計算フロー

```text
PLY / OBJ
   │
   ├─ Trimeshで読込（process=False）
   │
   ├─ loader-state topology QA
   │
   ├─ exact-coordinate weld
   │
   ├─ weld後 topology QA
   │
   ├─ 座標をmmへ統一
   │
   ├─ surface voxelization
   │
   ├─ surface connectivity QA
   │
   ├─ validated cavity seed探索
   │
   ├─ Z <= k に制限した6-connected flood fill
   │
   ├─ safe / spill境界をbinary search
   │
   ├─ safe componentのvoxel数を計数
   │
   └─ V = N × pitch³
```

### 6.1 flood fill

各候補水位`k`について、free voxelのうち`Z <= k`だけを通行可能とする。内部seedから6-neighborでbinary propagationし、連結領域がpaddingしたvoxel grid外周へ到達するかを判定する。

### 6.2 binary search

水位が高くなるほど許可されるfree-spaceは単調に増えるため、

- leakしない最高水位 = safe level
- 初めてleakする水位 = spill level

を二分探索できる。

### 6.3 容量

safe levelにおける内部連結領域のfluid voxel数を`N`、voxel pitchを`p` mmとすると、容量は`N × p³` mm³である。

---

## 7. QA設計

容量値だけでなく、次の情報を保存して計算過程を監査できるようにした。

- exact weld前後のvertex / face数
- geometry preserved exact
- boundary / non-manifold edge数
- winding consistency
- voxel grid shape
- surface voxel数
- 6 / 18 / 26 connectivityによるsurface component数
- seed位置とseed validation記録
- safe / spill level
- fluid surface PLY
- spill-level region PLY
- pitch別result.json
- validation_summary.csv / json

PyMeshLabは独立cross-checkとして位置づけ、容量計算用meshの生成には使用しない。

---

## 8. 実行環境

実資料validationを実施した環境は次のとおりである。

| 項目 | バージョン |
|---|---|
| PotteryVolumeCalculator | 1.3.0 |
| Python | 3.14.6 |
| macOS | 26.5.2 arm64 |
| NumPy | 2.5.2 |
| SciPy | 1.18.0 |
| Trimesh | 5.0.0 |
| Pillow | 12.3.0 |
| PyMeshLab | 2025.7.post1 |

PyMeshLabはimport可能であったがplugins loaded = 0で、必要filterは利用できなかった。このためPyMeshLab cross-checkは計算に寄与せず、NumPy + Trimesh基準QAで処理を継続した。

---

## 9. topology前処理の検証

v1.3.0の実行では、Trimesh読込後に358,804 vertices / 700,000 facesとなり、exact-coordinate weld後は350,067 vertices / 700,000 facesとなった。

- removed vertices: 8,737
- faces changed: 0
- geometry preserved exact: True
- weld後 boundary edges: 0
- weld後 non-manifold edges: 1

ただし、このpre-weld topologyについては後述するtexture UV seam問題が判明したため、**現在は考古学的なraw boundary情報として使用しない**。

---

## 10. texture UV seam問題とv1.3.1の対応

元PLY header上のvertex数は350,070である。ところがTrimeshの通常読込では358,804 verticesへ増加し、17,290 boundary edgesが検出された。

監査の結果、元PLYはface-corner texture coordinateを持っており、Trimeshがtexture UVの不連続を表現するために同一位置vertexを再indexすることが原因と判明した。

元PLY topologyを保持して読むと、

| topology | vertex | boundary edge | connected component | non-manifold edge |
|---|---:|---:|---:|---:|
| Trimesh既定texture処理後 | 358,804 | 17,290 | 159 | 0 |
| 元PLY topology保持 | 350,070 | **0** | **39** | **1** |

となり、後者はMeshLabでの確認値と整合する。

### 10.1 v1.3.1の暫定方針

このためv1.3.1では、接合境界と誤認される可能性がある次の出力を**生成しない**。

```text
archaeological/raw_boundary_*
archaeological/raw_components.*
boundary_before_after_comparison.json
spill_vs_raw_fragment_boundaries.ply
spill_boundary_proximity.json
```

`archaeological/`には当面、exact-coordinate weld後のQA用情報だけを残す。

```text
archaeological/
├── after_exact_weld_boundary_stats.json
├── after_exact_weld_components.csv
├── after_exact_weld_components.json
└── after_exact_weld_components_colored.ply
```

これらも実破片を自動同定した結果ではなく、**weld後meshのtopology QA**として扱う。

raw boundary関連処理は容量計算には使用していなかったため、この停止によってspill探索・fluid voxel計数・pitch validationのアルゴリズムは変化しない。

### 10.2 v1.3.1 smoke test

小型synthetic vesselを0.5 mm pitchで実行し、容量計算が正常終了するとともに、`archaeological/raw_*`およびraw-spill比較ファイルが生成されないことを確認した。

---

## 11. 中央縦断面とspill level

[断面の図版（PDF p. 11）](https://github.com/kotdijian/PotteryVolumeCalculator/blob/93fbdbad24f22c38423d9e1b4780622b2ef8ccab/docs/validation/PotteryVolumeCalculator_development_validation_report_v1.3.1.pdf#page=11)

**図3　中央縦断面とspill level。** 0.5 mm pitchではspill upper level = 241 mm、1.0 / 2.0 mm pitchでは242 mmである。

---

## 12. multi-pitch validation結果

| pitch | safe level | spill upper level | fluid voxels | capacity | total time |
|---:|---:|---:|---:|---:|---:|
| 2.0 mm | 240.0 mm | 242.0 mm | 725,340 | **5.802720 L** | 14.84 s |
| 1.0 mm | 241.0 mm | 242.0 mm | 5,917,481 | **5.917481 L** | 33.25 s |
| 0.5 mm | 240.5 mm | 241.0 mm | 47,583,957 | **5.947995 L** | 132.08 s |

pitchを細かくすると容量は、

```text
5.802720 → 5.917481 → 5.947995 L
```

と単調に増加した。

- 2.0 → 1.0 mm: +114.761 mL
- 1.0 → 0.5 mm: +30.514 mL
- 1.0 mmと0.5 mmの差: **0.513%（0.5 mm値に対して）**
- spill upper levelの全変動幅: **1.0 mm**

[容量収束の図版（PDF p. 12）](https://github.com/kotdijian/PotteryVolumeCalculator/blob/93fbdbad24f22c38423d9e1b4780622b2ef8ccab/docs/validation/PotteryVolumeCalculator_development_validation_report_v1.3.1.pdf#page=12)

**図4　voxel pitchと最大保持容量。** pitch細分化に伴い容量値の変化量が縮小している。

[spill収束の図版（PDF p. 13）](https://github.com/kotdijian/PotteryVolumeCalculator/blob/93fbdbad24f22c38423d9e1b4780622b2ef8ccab/docs/validation/PotteryVolumeCalculator_development_validation_report_v1.3.1.pdf#page=13)

**図5　voxel pitchとspill upper level。** 3解像度の差は1 mm以内である。

### 12.1 経験的収束診断

3点から`V(p)=V∞+ap^q`を仮定した診断では、

- empirical order `q` = **1.911**
- extrapolated `V∞` = **5.959046 L**
- 0.5 mm実計算値との差 = **11.1 mL**

となる。

この外挿値はモデル仮定に基づく**補助診断値**であり、測定値として扱わない。本資料の詳細計算値としては、実際に計算した、

> **0.5 mm pitch: 5.947995 L ≒ 5.95 L**

を用いるのが適切である。

---

## 13. fluid surfaceの形状確認

[fluidの図版（PDF p. 14）](https://github.com/kotdijian/PotteryVolumeCalculator/blob/93fbdbad24f22c38423d9e1b4780622b2ef8ccab/docs/validation/PotteryVolumeCalculator_development_validation_report_v1.3.1.pdf#page=14)

**図6　0.5 mm pitchで得たfluid surfaceと元モデルの重ね合わせ。** fluid regionは底部から器内へ連続し、上面は口縁直下のspill levelで水平となる。

---

## 14. voxel connectivity

| pitch | grid shape | surface voxels | components 6 / 18 / 26 |
|---:|---|---:|---|
| 2.0 mm | 117 × 116 × 142 | 113,536 | 1 / 1 / 1 |
| 1.0 mm | 232 × 231 × 285 | 455,599 | 1 / 1 / 1 |
| 0.5 mm | 464 × 459 × 569 | 1,822,698 | 6 / 3 / 2 |

0.5 mmではsurface voxelが複数componentに分離するが、validated seed、spill search、fluid regionはいずれも正常に求まり、容量値も1 mmから連続的に収束した。

今後はcomponent数だけでなく、各componentのvoxel数とlargest-component fractionを保存することで、小規模な孤立componentと重大なsurface breakを区別しやすくなる。

---

## 15. 実用上の評価

今回の資料では、

- spillが口縁近傍で安定
- fluid surfaceが器内形状と整合
- 2 → 1 → 0.5 mmで容量が単調収束
- 1 mmと0.5 mmの差が約0.51%

であった。

したがって、本資料については容量計算の中核である、

> **surface voxelization + validated cavity seed + height-limited flood fill + first-spill search**

が良好に機能したと評価できる。

運用上は、

- **詳細検証:** 0.5 mm
- **大量処理の候補:** 1.0 mm

という使い分けが妥当である。ただし大量一括処理へ移る前に、異なる器形・寸法の代表資料でもvalidationを行う必要がある。

---

## 16. 出力監査上の注意

同じoutput directoryを再利用すると、以前のrunで生成された`error.json`やdebug用PLYが残る場合がある。正式な検証では、

- 新しいoutput directoryを使用する
- または既存outputを退避・削除してから実行する

ことを推奨する。

v1.3.1 READMEにもこの注意を明記した。

---

## 17. 大量一括処理へ進む前の課題

### 17.1 PLY topology読込の整理

次版では、容量・topology QA用途について元PLY face topologyを保持する読込方法を明示し、texture処理とgeometry topologyを分離することが望ましい。

### 17.2 接合境界抽出の再設計

接合境界は単純なboundary edgeではなく、

- 元PLY topology
- face-corner UV discontinuity
- surface geometry
- curvature / dihedral angle
- vertex RGB / texture

を別々の情報層として解析し、複合的に判定する必要がある。これは容量計算から独立した機能として設計する。

### 17.3 run manifest

一括処理では、各資料について、

- input file hash
- program version
- unit
- pitch mode
- start/end time
- status
- result file

を1行ずつ記録するmanifestが必要になる。

---

## 18. 結論

`0015Jinmen_small.ply`に対するmulti-pitch validationでは、最大保持容量は、

- 2.0 mm: **5.802720 L**
- 1.0 mm: **5.917481 L**
- 0.5 mm: **5.947995 L**

となり、spill upper levelも241–242 mmに安定した。したがって本資料の詳細計算値は、**約5.95 L**と評価できる。

一方、監査により、pre-weldで検出された17,290 boundary edgesは、考古学的接合境界ではなくtexture UV seamによる頂点再indexの影響を強く受けていることが判明した。このためv1.3.1では`archaeological/raw_*`出力を停止した。

容量計算と接合境界抽出を明確に分離したことにより、次段階では、まず接合境界抽出を独立アルゴリズムとして検討し、その後にsingle-pitch modeを用いたフォルダ内一括処理へ進む構成が適切である。

---

## 付録A　主要監査ファイル

```text
0015Jinmen_small.ply
vessel_voxel_volume(20260816-021426).py
environment.txt
pip_freeze.txt
0015Jinmen_small_PotteryVolume_v1.zip
```

## 付録B　v1.3.1の主要出力構成

```text
<stem>_PotteryVolume_v1/
├── archaeological/
│   ├── after_exact_weld_boundary_stats.json
│   ├── after_exact_weld_components.csv
│   ├── after_exact_weld_components.json
│   └── after_exact_weld_components_colored.ply
├── processed/
├── qa/
├── pitch_2mm/
├── pitch_1mm/
├── pitch_0p5mm/
├── validation_summary.csv
└── validation_summary.json
```

**`archaeological/raw_*`はv1.3.1では生成しない。**
