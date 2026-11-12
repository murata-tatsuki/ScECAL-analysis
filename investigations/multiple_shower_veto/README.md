# 前方の複数成分による veto の簡易調査

5、40、100 GeV の現行 `Calib_Hit` を読み、分離した信号成分が前方の複数層で反復する事象を調べる。
正式な解析マクロ、元データ、較正、既存の resolution 処理は変更しない。
結果は `investigations/result/multiple_shower_veto/` に保存する。

## 判定

- 前方 layer 0–7。偶数層は短辺方向 Y、奇数層は短辺方向 X の投影を独立に使用する。
- 同じ短辺座標の信号を合算する。memory hit と長辺方向の異なる segment も投影上は合算する。
- 投影 bin の信号が 0.5 MIP 以上の点を、隣接 pitch 5.3 mm で連結する。1 MIP は 0.305 MeV。
- 二つの成分のそれぞれが 1 MIP 以上、かつその層の投影総信号の 15% 以上であることを要求する。
- 二成分の重心間隔が 15.9 mm 以上で、同じ投影の4層中3層以上に二つの位置が各 ±5.3 mm 以内で反復すると nominal 候補とする。
- bad-channel mask の穴による分断を抑えるため、mask 適用前後の両方で同じ二成分の反復を要求する。
- X または Y で成立した場合に veto。両投影で成立した件数も保存する。
- 必要層数2、間隔21.2 mm、成分比10%/20%の条件も事前固定して比較する。resolution の改善で条件を最適化しない。

これは複数成分の**形状 tag**である。複数粒子入射、上流で分岐した単電子 shower、異時刻 hit の混在を識別したものではない。
重なった粒子、近接した成分、早い層で信号が十分にない事象への検出効率は、この調査だけでは測れない。

## 入力・既存の選別

入力は既存 HG–LG tail 調査の manifest に列挙された2023年電子 data と `e_ssa/threshold` 単電子 MC。
manifest の旧 prefix は現サーバーの `data_beamtest` に明示的に読み替え、実在、entry 数、採用 entry の Run/Event/Time を検査する。
各 sample/energy から2万事象を、全ファイルに entry 数比例で割り当て、各ファイル内の層化ランダム連続区間から採用する。
入力 path・size・mtime・区間と event の識別情報を保存する。

信号和は現行 resolution と同じ非負 good hit と mask を使用する。追加の総 energy 窓、track veto、leakage cut は使用しない。
CoG は layer 9、10 の energy 加重 X/Y それぞれの矩形条件で、200/20/10/5 mm。
200 mm は既存表示規則に合わせて図中 `nocut` とするが、現行 resolution と同様、両層に正の信号がある条件は残る。

data は保存済み HG–LG tail-rich 判定を ID 照合して除外する。
MC の raw/calib は旧 tail 調査後に内容が更新されていたため、旧 MC event flag は再利用しない。
`audit_mc_tail.py` で、保存済み common-channel/residual-center/noise モデルを固定して、**現行 decode と calib**から判定し直す。
これは現行 MC でモデルを学習し直した完全な tail 再解析ではない。入力・モデル・pedestal の由来と検査結果は audit JSON に保存する。

## 幅の評価

- 全体の RMS/mean と、`(q84-q16)/(2*median)` を併記する。
- Gaussian core は主ピークを初期化して2段階で求めた簡易 fit。before の最終 fit 範囲を after に固定する。
- 主ピーク fit の探索範囲は event cut ではない。正式な resolution マクロと同じ fit ではないため、正式値の置換には使わない。
- nominal の RMS/mean と68%幅の before/after 差には、同じ event を再抽出する400回の paired bootstrap 区間を保存する。
  区間は採用 sample 内の event 統計だけで、run/spill 内相関、sampling、較正、モデルの系統誤差を含まない。
- MC の tag 率は単電子応答に対する参考値であり、data の候補純度や真の contamination 率ではない。0件でも偽判定率が厳密に0とはいえない。

## 再現

ROOT 入力を読む uproot、awkward、NumPy、SciPy、Matplotlib を使用する。
プロジェクトのルートから実行する。既存の同名結果を全面再生成する場合は別の `--out` を指定する。

```bash
python3 investigations/multiple_shower_veto/check.py --energy 5
python3 investigations/multiple_shower_veto/check.py --energy 40
python3 investigations/multiple_shower_veto/check.py --energy 100
# 保存したモデルを固定して現行MCの採用entryを再判定する。
python3 investigations/multiple_shower_veto/audit_mc_tail.py --energy 5 40 100 --target 20000
# audit CSVをID照合して適用し、保存済み特徴から最終集計・図を再生成する。
python3 investigations/multiple_shower_veto/check.py --energy 5 --summarize-only
python3 investigations/multiple_shower_veto/check.py --energy 40 --summarize-only
python3 investigations/multiple_shower_veto/check.py --energy 100 --summarize-only
python3 investigations/multiple_shower_veto/check.py --combine
```

`comparison/{metrics,tag_rates,paired_changes,run_rates}.csv` に全 energy・CoG・比較条件をまとめる。
`comparison/summary20mm.csv` は CoG 20 mm の代表値、`comparison/width_summary.png` は4条件の68%幅の比較。
`data/diagnostics/<energy>GeV/` と `mc/diagnostics/<energy>GeV/` に event 特徴、energy 分布、候補図を保存する。
条件別の数値と histogram は `<sample>/<condition>/<CoG>mm/<energy>GeV/`、入力と由来は `inputs/<energy>GeV/` に保存する。
`complete.json` の `current_mc_tail_audited=true` を満たした3energyだけを `--combine` で統合できる。
