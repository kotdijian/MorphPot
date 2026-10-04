# MorphPot DevelopmentReport

更新日：2026-10-05 JST  
状態：0.1.0-dev1 — 既存計算モジュールの移管・初期CLI統合

## 1. 開発方針
土器と石器を独立したアプリとして開発する。MorphPotは土器容量、器軸位置、連続断面、表面特徴を担当する。姿勢・原点の正規化はOrthoMakerが上流で行い、MorphPotはそのPLYを基本入力とする。

座標軸は全ツール共通の固定規約にせず、処理ごとに定義する。現行の継承アルゴリズムはZを器高／器軸としている。任意軸や投影方向を指定するadapterの共通化は今後行う。実測モデル、復元器形、ガイド、fit後配置は区別する。

## 2. 今回の移管とマージ
- RadialSections v0.7.0：全縦断面のXY重ね合わせ（点群・線分PLY）を既定で自動出力。器軸基準の符号付き距離と元Z高さへ投影し、入力単位を保持。
- 同版：左右を保持した内外面の中央値輪郭を、口唇で分割した外面・内面の正規化弧長登録による座標中央値として出力。正側口唇を頂点0とする閉ポリラインPLY／CSVとQAを追加。開曲線・複数成分・分岐は除外し、採用2方向未満なら中央値を生成しない。
- PotteryVolumeCalculatorのcommit d756a5667f860705fa5ffaef8e4fb81dcd93e059からvoxel容量v1.3.2、破片境界v0.3.0、多解像度特徴v0.1.0、表面痕跡v0.1.0を移管。
- 同commitのRadialSections v0.4.0を配布v0.6.0系へ統合。基本断面取得部分は同系統で、容量積分はpottery_volume_core.pyへ分離された版を採用。
- DrawingCapacityも配布v0.6.0へ揃え、3D断面／2D図面で同じ数値コアを使う。
- PotteryReconstruction3Dに収録されたSectionGuide v0.1.2をコピー。
- 配布済みReconstruction v0.2.5とGuideFit v0.1.0を、後段の独立モジュールとして収録。
- Surface Enhancement Lab v0.7.3をexperimental/へ収録。継承用コードであり、MorphPotの統合機能にはまだ接続していない。
- 移管元ファイルは残し、上流repoとOrthoMaker v1.0.0 tagは変更していない。
- 新CLI MorphPot.pyは必要なアルゴリズムのみを読み込む。GUIは必須依存にしない。

## 3. 実装済み範囲
| モジュール | 実装・役割 |
|---|---|
| MorphPot.py | volume / sections / surface-trace / fragment-boundaryの入口 |
| morphpot/asset_metadata.py | 単位・メタデータ解決、対応PLYハッシュ検証。Transform行列の適用は行わない |
| vessel_voxel_volume.py | 全体voxel化、内腔seed検証、高さ制限flood fill、first-spill直前の最大保持容量、複数pitch検証、QA |
| pottery_radial_sections.py | 水平断面楕円fit、内外輪郭、MAD中心選択、Z平行器軸位置、放射全断面／半断面、容量4方式 |
| morphpot/section_overlay.py | 全断面のXY重ね合わせ、左右の口唇登録、内外面中央値の閉ポリライン、採用／除外QA |
| pottery_volume_core.py | single / optimized / angular / ellipseの共通数値積分 |
| pottery_drawing_capacity.py | 2D実測図の縮尺校正・内面デジタイズ・Drawing-Single容量、既存Tk GUI |
| pottery_section_guide.py | 実測図から独立内外左右profile、水平ring、OBJ／PLYガイド、既存Tk GUI |
| pottery_reconstruction_3d.py | 現配置保持の破片解析、器形、ガイド登録・QA。既存Qt／PyVista GUI |
| pottery_guidefit_reconstruction.py | Guide-derived器形と破片内外同時fit、放射移動・小tilt、観測／fit結果分離 |
| pottery_multiscale_features.py | 内外面分類、同面graph diffusion、物理mm尺度の多解像度特徴 |
| pottery_surface_trace.py | 表面特徴応答と予備的形態状態、CSV／PLY／JSON出力 |
| pottery_fragment_boundary.py | geometryのみの異常スコア、反対面支持・厚み異常、境界候補とQA |
| experimental/surface_enhancement_lab.py | Base／Top／線状性・連続性・合成・vertex color bakeの継承用Lab |

表面痕跡・接合境界はgeometryの記述・候補スコアであり、考古学的技法や実際の接合境界を自動確定する機能ではない。粗いメッシュでは要求スケールと実現スケールに差が生じるため、QA出力を参照する。

## 4. 入出力・単位の現段階
MorphPot.pyはnormalized PLYを受け取り、次の順序で単位を解決する。
1. 同名の<stem>.asset.jsonがあればschema_version=1.0、triangle_mesh、SHA-256、coordinates.unit／unit_to_mmを検証する。
2. 明示--unit mm|cm|mがあればそれを使う。ただしasset JSONとの矛盾やハッシュ不一致を黙って無視しない。
3. autoでは隣接transform.jsonのsource.input_unit／unit_to_mmを読む。coordinate_values_rescaledがtrueの場合は自動推定しない。
4. 単位不明なら停止し、明示指定を求める。bboxから単位を決めない。

旧transform.jsonのsource.sha256はrawモデルのハッシュであり、normalized PLYとの比較には使わない。読込時にpose matrixを再適用しない。

asset JSONは現時点では単位・ファイル同一性の読込機能のみ。完全なschema、フレーム／変換履歴の出力、metadata生成、BagIt保存は未実装。legacy metadataはハッシュで出力PLYへ結び付けられない制限がある。

従来の各.pyを直接実行した場合は従来CLIの単位指定・既定値が残る。auto解決を使う場合はMorphPot.pyを使う。個別アルゴリズムの出力名・座標単位・QA構造は保持している。

## 5. 開発者用実行例
READMEは概要のみとし、現段階の実行例はこの文書に集約する。

初回環境：
~~~bash
python -m venv venv
source venv/bin/activate
python -m pip install -r requirements.txt
~~~
Windowsでは有効化をvenv\\Scripts\\Activate.ps1へ読み替える。

~~~bash
python MorphPot.py volume /path/pot001_rev.ply --unit auto --pitch 1
python MorphPot.py sections /path/pot001_rev.ply --unit auto --angle-step 30
python MorphPot.py surface-trace /path/pot001_rev.ply --unit auto
python MorphPot.py fragment-boundary /path/pot001_rev.ply --unit auto
~~~
隣接metadataがなければ--unit mm/cm/mを指定する。各計算のオプション詳細は対応.pyの--help。

図面・再構成GUIは独立入口を維持：
~~~bash
python pottery_drawing_capacity.py
python pottery_section_guide.py
python pottery_reconstruction_3d.py
python pottery_guidefit_reconstruction.py
~~~

requirements-gui.txtは再構成Qt GUI、PDF／図面対応の追加依存。TkinterはPython本体側のTk対応が必要。
requirements-surface.txtは実験Labの追加依存。
requirements-qa.txtは任意PyMeshLab照合。PyMeshLabはvoxel計算の必須コアではない。

Surface Labを直接起動する場合は、現行Labのmm前提・投影規約に適合した入力が必要。m入力のauto換算や土器投影adapterはまだ接続していない。

## 6. 検証結果
今回の環境：Python3.12、NumPy2.5.3、SciPy1.17.0、Trimesh5.1.1、Pillow12.3.0、NetworkX3.7。正式な対応版固定は今後行う。

実施済み：
- v0.7.0追加検証：既存metadata3＋overlay7＝pytest 10件passed。断面の方位・線分順序、mm/cm/m同値性、外れ方向に対する中央値、連続edgeと始点、開曲線・複数成分・余剰開線の除外、失敗／無効化時の古い出力削除を確認。
- 人工開口容器のCLIで、容量・PNGを省略した場合も6方向の重ね合わせと472頂点の中央値輪郭が出力されることを確認。実資料の波状口縁・欠損資料による精度検証は未実施。
- 全Pythonファイルのcompileall。
- 新metadata契約テスト3件passed（旧rawハッシュ／行列の非再適用、asset一致と衝突検出、単位欠落・不整合）。
- PotteryVolumeCore self-test：人工円柱の理論容量との一致。
- DrawingCapacity self-test：同円柱の容量一致。
- SectionGuide self-test：軸・ring・単位換算。
- Reconstruction self-test：器軸、厚み、ガイド登録。
- GuideFit self-test：前後RMSEと固定並進／周方向条件。
- 512face人工開口容器でMorphPot CLIのvoxel容量、断面4容量方式、surface-trace、fragment-boundaryを最後まで実行し出力確認。

voxelの粗いpitchによる値は近似であり、人工容器smokeはパイプライン接続確認。考古資料での容量精度検証やGUI実機確認を今回やり直したものではない。Surface Lab GUIはcompileのみでレンダリング／bakeの実機回帰を実施していない。

## 7. 今後の予定
優先順：
1. asset metadataのschema・出力生成、現在フレームと変換履歴、処理条件・入力ハッシュの統一。
2. 実資料でvoxel／profile容量の収束・差分検証。入力単位m/cm/mmの同値性確認。
3. 土器水平／縦断面の形態指標を拡張（現行はLithMorphと同じ全断面指標群を実装済みではない）。
   口唇の現行検出は側ごとの最高点（平坦な口唇は中央）。波状口縁・装飾・欠損を含む資料での検証、必要に応じ手動の口唇指定や対応付けの改善を行う。
4. 各処理の軸・観察方向adapter。Surface手法を共有し、土器用投影と石器用投影を分ける。
5. 共通GUIを整え、独立Tk／QtモジュールのUIを計算コアから分離。
6. 破片器形／ガイド器形／実測器形の結果区分、外面＋器厚から内面推定。
7. QAの統一、CI、必要に応じBagItで成果保存。

瓦はOrthoMaker瓦モード導入後の別サブプロジェクト。OrthoMakerの大規模レンダリング・展開最適化は保留。

## 8. 出典・ライセンス
新repoのMIT LICENSEを保持。過去CC0の出典や版履歴はSOURCE_PROVENANCE.jsonおよびTHIRD_PARTY_NOTICES.mdに記録する。移管元全体のGUIアプリや大型サンプルを丸ごと複製せず、必要コードのみ収録した。
