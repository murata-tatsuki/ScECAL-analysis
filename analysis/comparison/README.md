# CoG と channel 座標による data / simulation 比較

layer 9・10 の重心でイベントを選び、そのイベントの指定範囲内 channel のみを集計する。
どちらの範囲も x・y 方向の半幅 [mm]。境界は含まない。

- `cog_range > 0`: layer 9・10 の両方で `abs(CoG_X) < cog_range` かつ
  `abs(CoG_Y) < cog_range` を要求する。重心を定義できない層（エネルギー和が0）が
  あるイベントは除外する。`cog_range=0` は CoG カットなし。
- CoG は **channel カット前**の全有効 hit のエネルギー重み付き重心。
  従来どおり負の hit energy と、設定に応じた bad channel を除く。
- 集計する hit は `abs(Hit_X) < channel_range` かつ
  `abs(Hit_Y) < channel_range` を満たす layer 0〜29 の channel。
  `channel_range` は正の値を指定する。200 mm は検出器全体を含む。

CoG の図も channel カット前の重心を表示する。
`beforeEventCut` / `afterEventCut` は **どちらも CoG と channel カット適用後**で、
従来の energy–Nhit イベントカットの前後を表す。
energy–Nhit カットは範囲内 channel の energy と Nhit で計算する。
comparison は従来どおり `beforeEventCut` を比較する。
同じ Nhit での energy は各 ROOT の `beforeEventCut/hit_vs_e_<energy>GeV`
（x: energy、y: Nhit）から確認できる。

## 実行

ROOT と g++ を使える環境で、このディレクトリから実行する。

```bash
make
# コマンドの確認のみ（ジョブは投入しない）
DRY_RUN=1 bash analysis_paralell_allCondition.sh
# エネルギーごとの ROOT を生成する SLURM ジョブを投入
bash analysis_paralell_allCondition.sh
# 上記の全ジョブが完了してから data / simulation を比較
bash compare_ds_allCondition.sh
```

**両方の shell で同じ配列を設定する。** 初期設定は次の8通り。

```bash
cog_ranges=(0 20 10 5)
channel_ranges=(68.2 22.5)
```

CoG なし/20/10/5 mm × channel 68.2/22.5 mm の全組合せを回す。
凡例も配列の順に表示し、上から下へ CoG カットが厳しくなるようにする。
例えば `cog_ranges=(0 20 10 5 2.5)` と `channel_ranges=(22.5)` なら
22.5 mm の channel 集計で CoG の条件だけを比較する。小数も指定できる。
両 shell では小数の表記も一致させる（`5` と `5.0` は別の出力パス）。
200 mm の CoG カットは、CoG を定義できないイベントも除外するため、0 と異なる。
calibration の入力先と simulation の threshold 条件は既存設定を引き継いでいる。

## 出力

出力先は `../result/comparison/`。
`condition=cog<cog_range>mm_channel<channel_range>mm` として保存する。

- `data/<condition>/<energy>GeV_<condition>.root`
- `simulation/<threshold>/<condition>/<energy>GeV_<condition>.root`
- `figures/data/<condition>/`、`figures/simulation/<threshold>/<condition>/`
- 比較 ROOT: `comparison_dataThre_cog_channel.root`
- 比較図: `figures/comparison/dataThre_cog_channel/`

例: `data/cog5mm_channel22.5mm/100GeV_cog5mm_channel22.5mm.root`。
各単一エネルギー ROOT の最上位に `cog_range_mm` と `channel_range_mm` も保存する。
以前の channel のみの結果、`../result/resolution/` の結果とは別に保存する。
SLURM ログは `comparison/job/test-<job ID>.out` / `.err`。

比較図は CoG/channel の組合せごとに色・形を割り当て、同条件の data は
塗りつぶし＋実線、simulation は白抜き＋破線で表示する。
凡例の値は `5 / 22.5` のように **CoG / channel** の順とし、単位と順序は
凡例の見出しにまとめる。CoG の0は `no cut / 22.5` のように `no cut` と表示する。
number of hits の表示範囲は全条件の分布に合わせる。

resolution 図の凡例は右上パネルだけに表示する。fit の式は共通の見出しにまとめ、
各行は条件と係数 `(a, c)` または `(a, b, c)` を表示する。
条件が多い場合は2列にし、行数に応じて文字サイズを調整する。

- `resolution.png`: noise 項なし、linear–linear
- `resolution_loglog.png`: noise 項なし、log–log
- `resolution_with_noise.png`: noise 項あり、linear–linear
- `resolution_with_noise_loglog.png`: noise 項あり、log–log

log–log は同じデータ点・fit 結果を使い、再fit は行わない。
energy と resolution がともに正かつ有限な点を表示する。
y 軸の上限は2に固定し、下限は正の値・誤差から設定する。
fit の除外範囲（25 < E < 35、75 < E < 85 GeV）でも曲線は同じ係数で連続して描画する。
非正値の点は log 図だけで非表示とし、元データ・linear 図は維持する。
linear–linear は従来の表示範囲を維持する。
比較 ROOT にも同名（拡張子なし）の4つの canvas を保存する。

さらに、同じ CoG / channel 条件の data と simulation だけを重ねた図を追加する。
条件ごとに noise 項なし・ありを1枚の PNG にまとめる。data は青、simulation は橙。

| | 左: linear–linear | 右: log–log |
| --- | --- | --- |
| 上段 | noise 項なし | noise 項なし |
| 下段 | noise 項あり | noise 項あり |

各パネルの右上に data / simulation の2本の fit 関数を表示する。
係数名 `(a, b, c)` の一覧ではなく、`σ(E)/E = 0.221/√E ⊕ 0.279/E ⊕ 0.010`
のように数値を代入した式とする（`⊕` は二乗和の平方根）。
小さな係数は指数表記で表示する。軸範囲は全条件図と揃え
（linear は x: 0〜130 GeV、y: 0〜0.3、log の y 上限は2）、点・fit 係数も同じ値を使う。
全条件を重ねた既存4枚も引き続き出力する。

ファイル名は `resolution_pair_<condition>_<threshold>.png`。
例: `resolution_pair_cog0mm_channel22.5mm_threshold.png`。
現在の8条件では合計8枚を、既存図と同じ比較図ディレクトリに保存する。
simulation の threshold が複数ある場合は、それぞれ data と1対1で比較する。
比較 ROOT にも PNG と同名（拡張子なし）の4パネル canvas を追加する。
`bash compare_ds_allCondition.sh` で既存図と一緒に生成する。


## 単独実行・従来形式

```text
SingleEnergyAnalysis output.root Nenergy E1...En E1Nfiles...EnNfiles input1.root ... cog_range_mm channel_range_mm b_only_best figure_path
```

末尾の範囲引数は CoG、channel の順。
以前の channel のみの呼出し（`cog_range_mm` 省略）も CoG カットなしとして利用可能。
`MultiEnergyAnalysis` は旧 `100GeV_22.5mm.root` 形式も読み込める。

Gaussian fit の方法・下限15 MeVは従来どおり。
狭い範囲・低エネルギー・厳しい CoG 条件では、統計と fit の成立を確認すること。

## 検証

```bash
python3.12 test_channel_cut.py
```

人工の calibration ROOT を実際の `SingleEnergyAnalysis` で処理する。
CoG なしの従来形式と新形式、CoG 2.5/5 mm × channel 22.5 mm を検証する。
channel 範囲外の hit が CoG に効くこと、エネルギー重み、矩形条件、両層の条件、
境界、重心を定義できない層、bad channel、負のエネルギーを含める。
イベントカット前後の総エネルギー・Nhit・30層の分布・SiPM別分布を確認する。
テストの ROOT・図・ログは `../result/comparison/test_channel_cut/` に保存する。
