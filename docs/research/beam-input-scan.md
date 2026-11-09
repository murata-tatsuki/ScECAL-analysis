# ビーム入射分布探索：ディレクトリ構造と出力仕様

記録日：2026-10-06。会話で合意した保存規則をまとめた仕様書。
ここに示す探索用コード・出力階層は今後の実装対象であり、MC探索を実行済みという意味ではない。

## 基本規則

- イベントROOT・集計ROOTは、`<condition>/<energy>GeV/*.root` の形で保存する。
- `100GeV/` などのエネルギーディレクトリが並ぶ親には、エネルギーディレクトリ以外のファイル・ディレクトリを置かない。
- ROOT用の各 `<energy>GeV/` の中にはROOTファイルだけを置く。summary・図・設定・ログを混在させない。
- 探索と独立検証を区別する `scan` / `validation` は、条件ディレクトリより上に配置する。
- 図、summary、入力設定・再現記録、Slurmログは、それぞれ専用の保存先で管理する。
- 既存のMCイベントROOTの工程・粒子・条件階層とファイル対応を継承する。
- エネルギー表記は既存と同じ `100GeV`、`0.5GeV` などとする。

## 条件の識別

| 名前 | 意味 | 例 |
|---|---|---|
| `label` | 固定した生成・digitization・校正・SSA・評価条件での探索一式 | `ssa_fixed` |
| `phase` | 候補探索または独立seedによる検証 | `scan`, `validation` |
| `condition` | 入射位置分布の条件 | `point`, `sx0p50_sy0p50` |
| `energy` | 電子のビームエネルギー（GeV） | `100` |
| ファイル番号 | 各条件内の生成単位。工程間で同じ番号を維持 | `100GeV_0.root` |

`sx0p50_sy0p50` は、元の入射位置分布の中心からのX・Yの距離をそれぞれ0.50倍する条件。
`point` は入射位置を固定する条件。中心位置、倍率、seed、生成イベント数の正確な値は設定とmanifestに保存する。

## MCイベントROOT

保存基点は既存の `simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/`。
`generate` は既存と同じく `threshold` 階層を持たず、`decode` / `calib` / `ssaCalib` は既存の `threshold` 階層を継承する。

```text
Result_MC/
  generate/e-/e_ssa/beam_scan/<label>/<phase>/<condition>/
    100GeV/
      100GeV_0.root
      100GeV_1.root
    200GeV/
      200GeV_0.root

  decode/e-/e_ssa/threshold/beam_scan/<label>/<phase>/<condition>/
    100GeV/
      100GeV_0.root
      100GeV_1.root

  calib/e-/e_ssa/threshold/beam_scan/<label>/<phase>/<condition>/
    100GeV/
      100GeV_0.root
      100GeV_1.root

  ssaCalib/e-/e_ssa/threshold/beam_scan/<label>/<phase>/<condition>/
    100GeV/
      100GeV_0.root
      100GeV_1.root
```

SSA後の具体例：

```text
Result_MC/ssaCalib/e-/e_ssa/threshold/beam_scan/ssa_fixed/scan/sx0p50_sy0p50/
  100GeV/
    100GeV_0.root
    100GeV_1.root
```

生成からSSAまで同じファイル名を維持し、対応する入力・イベントをmanifestで追跡する。
イベントROOTのTTree名・branch形式は既存処理を継承する。確認した既存100 GeVファイルでは、生成は `MC_Truth`、校正後は `Calib_Hit`、SSA後は `SSA_Hit`。

## 解析コード・解析出力

探索コードの予定保存先：`analysis/ssa_shower_shape/beam_input_scan/`。
解析出力の基点：`analysis/result/ssa_shower_shape/beam_input_scan/<label>/`。

```text
analysis/result/ssa_shower_shape/beam_input_scan/<label>/
  root/
    data/<phase>/reference/<energy>GeV/
      histograms.root
      diagnostics.root
    mc/<phase>/<condition>/<energy>GeV/
      histograms.root
      diagnostics.root
    comparison/<phase>/<condition>/<energy>GeV/
      comparison.root

  plots/
    <phase>/<condition>/<energy>GeV/
      beam_closure.png
      beam_closure_log.png
      layer0_xy_comparison.png
      rms_overview_count.png
      rms_overview_positive_energy.png

  summary/
    <phase>/<energy>GeV/
      ranking.tsv
      layer_statistics.tsv
      best_candidates.json

  inputs/
    <energy>GeV/
      configuration.json
      conditions.tsv
      files.tsv
      run.log
      sources/
```

dataの `reference` は比較基準。探索・独立検証に使うイベント範囲を設定に記録する。
`stage` / `frame` / `weight` を条件ディレクトリとエネルギーの間へ追加しない。
複数の集計条件は `histograms.root` 内部の `stage/frame/weight/` ディレクトリで区別する。
定義は既存解析の `pre` / `post_split`、`lab` / `aligned`、`count` / `positive_energy` などを継承する。
入射分布のROOTを保存する場合も、ROOT用の `<condition>/<energy>GeV/` 以下に `beam_input.root` として置き、manifestから参照する。点入射は固定座標を設定に記録する。

## 出力形式

| ファイル | 内容 |
|---|---|
| `histograms.root` | 集計ヒストグラム、規格化前の重み、採用イベント数 |
| `diagnostics.root` | イベント・信号・再構成の診断情報 |
| `comparison.root` | 再編集できる比較canvas |
| `*.png` | 条件比較・分布・層別幅の図 |
| `ranking.tsv` | 条件、X/Y倍率、採用イベント数、一致度、bootstrapによる不確かさ、順位 |
| `layer_statistics.tsv` | data/MC、条件、stage、frame、weight、layer、平均、標準偏差、包含幅、裾の割合 |
| `best_candidates.json` | 採用候補、パラメータ、評価結果、設定・ROOTへの参照 |
| `configuration.json` | 入力分布、生成・校正・SSA条件、選別、bin、規格化、評価方法、計算量上限 |
| `conditions.tsv` | 実行条件、倍率、中心位置、探索・検証の区分 |
| `files.tsv` | 条件、phase、ファイル番号、seed、イベント数、各工程のROOTパス、処理状態 |
| `run.log` | 実行コマンド、主要条件、実行状況 |
| `sources/` | 使用したコード・設定の控えと版・checksumの記録 |

数値表はヘッダー付きのタブ区切りTSVとする。単位を列名に含め、位置・幅には `mean_x_mm`、`rms_y_mm` などを使う。

図の見出しにはエネルギー、入射条件、処理段階、layer、探索・検証の区分を表示する。
凡例は図上部の専用欄に配置し、ヒストグラム・グラフと重ならないようにする。
同じ条件のlinear/log表示で、比較対象・規格化・binを一致させる。

## Slurmログと実行管理

Slurmログの基点：`analysis/log/ssa_shower_shape/beam_input_scan/<label>/`。

```text
analysis/log/ssa_shower_shape/beam_input_scan/<label>/
  submissions_<投入ID>.tsv
  <phase>/<condition>/<energy>GeV/
    generate_<jobID>.out
    generate_<jobID>.err
    decode_<jobID>.out
    decode_<jobID>.err
    calib_<jobID>.out
    calib_<jobID>.err
    ssaCalib_<jobID>.out
    ssaCalib_<jobID>.err
    analysis_<jobID>.out
    analysis_<jobID>.err
```

- 異なる生成・校正・SSA・評価条件の探索は、意味の分かる別 `label` で管理する。
- 探索候補を追加する場合は同じ探索内に `condition` を追加し、設定・manifestを更新する。
- 同じ保存先への同時書込みを避ける。
- 同じ条件の再描画は既存の解析ROOTを再利用して、対応する図・比較canvasを更新する。
- MCイベントROOT・入力dataと、実行結果の由来を示す設定・ソース記録は保持する。

## 根拠・参照元

- ユーザー指定：`条件/*GeV/*root` とし、エネルギーディレクトリが並ぶ親にはそれ以外を置かない。summaryなどは別ディレクトリで管理する。
- [プロジェクトの配置規則](../../AGENTS.md)。特に調査資料、解析結果、Slurmログの保存先。
- [既存SSA shower shape解析の定義・出力](../../analysis/ssa_shower_shape/README.md)。
- 既存 `Result_MC/{generate,decode,calib,ssaCalib}/e-/e_ssa/` の工程別配置と100 GeVイベントROOTを確認して、形式を継承した。
