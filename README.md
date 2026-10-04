# MorphPot

**姿勢正規化済みの土器3Dモデルを対象とする、容量・器軸・断面・表面形態解析プロジェクトです。**

基本入力は[ArtefactsOrthoMaker](https://github.com/kotdijian/ArtefactsOrthoMaker)が出力する三角形メッシュPLYです。既に適用された姿勢Transformを再適用しません。

現在は容量計算、水平・放射断面、断面容量、実測図容量、3Dガイド、破片再構成の既存モジュールを継承した開発版です。土器表面痕跡と接合境界候補の解析も収録しています。各処理はCLIまたは継承元GUIから独立して実行でき、MorphPot全体の統合GUIは今後開発します。

Surface Enhancement Labは共通手法の継承用としてexperimental/に収録しています。土器の投影方向・単位に適合させた統合は今後の予定です。

石器の平面形態・連続断面は[LithMorph](https://github.com/kotdijian/LithMorph)で開発します。

実装範囲、モジュール構成、開発計画、検証結果は[DevelopmentReport](DevelopmentReport.md)を参照してください。移管元は[SOURCE_PROVENANCE.json](SOURCE_PROVENANCE.json)に記録しています。

本リポジトリは開発版であり、全資料での精度検証を完了した正式安定版ではありません。

