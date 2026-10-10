# 全体代表モデル生成法の探索：0.1.0-exploration

## 目的と対象

器形カテゴリに共通の解析基盤を検討するため、全長弧長対応・頸部区間別対応・区間別対応＋限定変形を同じ断面集合で比較する。現在の対象は頸部のくびれを持つ器形で、カテゴリ名は`necked_jar`。考古学的な器種同定、カテゴリごとの標準方式の確定、AOM連携は今回の実装範囲外。

既存の口縁中央線モデルとは別立ての全体モデルである。本実装の口唇節点は全体断面分割に用いる最高点であり、口縁中央線延長による始点とは定義が異なる。

## 実行

リポジトリのルートで実行する。

```bash
python -m pip install -r requirements.txt
python pottery_whole_shape_exploration.py "/path/to/Sample-0015Jinmen_modefied_base/Sample-0015Jinmen_modefied_rev.ply" --category necked_jar --output-dir exploration/0015Jinmen
python pottery_whole_shape_exploration.py "/path/to/MK18-Hajiki_2M_rev.ply" --category necked_jar --output-dir exploration/MK18_Hajiki
python summarize_whole_shape_exploration.py exploration/0015Jinmen exploration/MK18_Hajiki --output-dir exploration/research
python plot_whole_thickness_diagnostics.py exploration/0015Jinmen exploration/MK18_Hajiki
```

出力先は新規または空フォルダ。入力は既に正規化されたPLYとし、付属`transform.json`等から単位を解決する。metadataなしなら`--unit m`等を明示する。Transformは適用しない。復元部分は現時点では全利用する。

5°間隔で訓練断面を取得し、2.5°だけずらして同一メッシュから検証断面を新規取得する。`--angle-step`を変更すると検証開始角度もその半分に変わる。軸XYは水平断面40段の内面から推定し、傾きは診断のみ。目標ランドマーク・中央値は訓練断面から決定するが、軸推定には入力全体を使うため、独立個体での汎化検証とは異なる。

## 器形条件と変形

頸部候補は半径列を0.5 mm間隔で標本化し、約5 mm窓で平滑化して検出する。平滑化は検出専用。候補は弧長比0.45〜0.985、prominence 1 mm以上、胴部との半径差3 mm以上、口唇への半径増加1 mm以上。内外面の候補高さ差8 mm超、拮抗候補などは未確定とする。取得成功率80%以上・取得断面の支持率80%以上が生成を続ける条件。未確定・不適合なら診断を残して非ゼロ終了する。

限定変形は`r'=r+Δr(z)`、`z'=z+Δz(z)`。内外面に共通の空間変位場を適用する。底部保護域と頸部・口唇を節点としてPCHIP補間し、底部保護域以下を固定する。底部保護域は内外の軸上点の高い方から器高の10%を加えた高さ。境界で位置は連続だが微分の連続は保証しない。

変位5 mm、せん断勾配0.15、高さひずみ0.15を暫定上限とし、500高さ点で評価する。節点の高さ順序が不適切な場合も拒否する。拒否方向は除外せず原形状のまま集計する。全域アフィン変換・相似変換とは異なり、器厚保存を保証しない。

## 出力

| パス | 内容 |
|---|---|
| `analysis/raw_*_profiles.npz` | 元の訓練・検証断面（mm） |
| `analysis/landmarks.json`、`holdout_landmarks.json` | 各方向の支持状態とランドマーク |
| `analysis/global_arc/` | 全長対応モデル |
| `analysis/neck_segmented/` | 座標を動かさない区間別対応モデル |
| `analysis/neck_segmented_warp/` | 区間別対応＋限定変形モデル |
| 各方式の`whole_representative.ply` | 入力単位・入力姿勢での全体回転メッシュ |
| 各方式の`profile_distribution.npz` | 対応断面群と座標中央値（mm） |
| `validation/<方式>/*_original_all_overlay.*` | 全実測断面とモデルの重ね合わせ |
| 変形方式の`*_registered_all_overlay.*` | 変形後断面との重ね合わせ |
| `validation/distance_summary.csv` | 元座標・変形後を分けた片方向距離 |
| `validation/dimensions.csv` | 元断面・変形後断面・モデルの寸法 |
| `validation/thickness_stations.csv` | 法線交点測定と状態・交点位置 |
| `validation/method_comparison.png` | 3方式の全体・口縁拡大比較 |
| `validation/model_thickness_diagnostic.png` | 器厚診断スクリプトで追加する法線到達先 |
| `research/` | 集計スクリプトで追加する部位別RMS・寸法差CSVと図 |
| `exploration.json` | ハッシュ、単位、軸、採用方向、変形の適用・拒否 |

全体PLYは元入力との重ね合わせ用。`*_xy.ply`はX=半径、Y=入力Z、Z=0の投影断面で、全体メッシュとは座標配置が異なる。PLYは入力単位、CSV・NPZはmm。

`analysis`はsupplementary候補となる解析一次出力、`validation`は個別の方法検証、`research`は研究論文・レポートの解釈に用いる派生比較資料として整理する。現段階の図表は方法探索であり、考古学的な型式分類を実証したものではない。

## 2サンプルでの実行結果（2026-10-08）

両入力はbase姿勢・m単位。最大72半断面のうち、訓練の共通採用数は0015で65、MK18で66。検証は66、70。限定変形の適用数は訓練で13/65、65/66、検証で12/66、69/70だった。

元の検証断面からモデル線分への1 mm間隔の片方向距離について、内外面を合わせて集計したRMS：

| 個体 | 全長 mm | 区間別 mm | 区間別＋変形 mm |
|---|---:|---:|---:|
| 0015Jinmen | 2.446 | 2.451 | 2.424 |
| MK18-Hajiki | 1.319 | 1.323 | 1.324 |

頸部〜口唇に限定したRMSは0015で4.710 / 4.768 / 4.665 mm、MK18で1.951 / 1.953 / 1.941 mm。小幅な差であり、今回の結果から従来方式を置き換える明確な優位性は確認できない。0015は拒否方向が多く、全方向補正済みとは扱えない。MK18は変形が成立しやすいが最大径・最大径高さも変わるため、カテゴリ内の追加個体と特徴量再現性を検証して採否を判断する。

器厚診断は既存の無制約な外面法線と内面交点を利用する。0015の口唇付近では法線が離れた頸部側へ到達し、約31〜36 mmを返す。`status=ok`でも正しい局所器厚とは限らず、この値をモデルの実際の肥厚と解釈してはならない。対応範囲制約・口縁中央線測定との統合は未実装。

## コードの検証

```bash
python -m pip install -r requirements-dev.txt
python -m pytest tests/test_whole_shape_exploration.py tests/test_whole_vessel.py tests/test_whole_validation.py tests/test_whole_dimensions.py -q
```

合成断面で頸部検出・無頸形不支持・区間接続・内外共通変位場・底部固定・節点位置・過大変形拒否を検証する。実サンプル6メッシュでは閉鎖性・法線整合・正体積と画像を確認した。これらは器厚対応・局所形態・カテゴリの妥当性を全面的に保証しない。


## 原位置と特徴点位置合わせ後の説明図（2026-10-10）

`plot_whole_registration_comparison.py`を追加した。新規の探索実行では、holdout断面が存在する場合に`validation/registration_comparison/`へ自動出力する。既存出力からはこのスクリプトだけを実行できる。

```bash
python plot_whole_registration_comparison.py /path/to/0015Jinmen /path/to/MK18_Hajiki
```

`RegistrationComparison.md`は全体図、頸部・口縁拡大、統計グラフ、部位別統計表をtrainingとholdoutごとに掲載する。左図は原位置r,z、右図は保存済みの共通変位場を再現したr,z。これは頸部で対応点を分割するだけの`neck_segmented`とは異なる座標変形である。左右は同じ採用断面群・同じ軸範囲・同じ`neck_segmented_warp`代表輪郭を表示し、拒否方向を黄褐色で残す。底部保護域を固定する。

距離評価は元弧長1 mm間隔の同じ測点を使用し、右側で変形後の再標本化をしない。どちらも同じ代表輪郭線分への符号なし最近傍距離で、内外面の対応は保持する。平均、母SD、RMS、p95、最大、方位別RMS、原位置高さ5 mm区間別平均±SD、累積分布を図示する。部位区分も元断面から固定する。測点をプールするため長い断面・部位ほど重みが大きい。既存距離表の再標本化とは定義が少し異なる。

この出力は解析手法に馴染みのない読者へ位置合わせの意味を説明するもの。散らばりの減少を精度改善・完全な歪み除去・器厚保存の証明とは扱わず、改善のない方向もそのまま示す。全体モデルや既存の一次検証CSVは変更しない。
