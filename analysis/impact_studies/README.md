# Impact studies

Slurmログと投入記録の既定保存先は `analysis/log/impact_studies/` です。明示した `--log-dir` は優先されます。

data/MCの不一致について、条件を一つ変えたときの影響を確認する解析を置くディレクトリ。
C++/ROOTとshellのみを使用し、図はPNGのみ出力する。

## ディレクトリと命名規則

コードは検証項目ごとにC++ファイルを分ける。今回の `HLInterceptStudy.cc` に対応する
study名は `hl_intercept`。今後も検証項目を追加し、study名で出力を分離する。

```
analysis/impact_studies/
  HLInterceptStudy.cc    切片項だけを除いたMCと既存MC/dataの比較
  StudyCommon.hh        入出力、統計量、既存校正定数の読み込み
  EventHistograms.hh    イベント分布からrange/bin幅を決定、概要図を描画
  ReplotEventEnergy.cc  完了済み結果のイベント分布を再描画
  SummarizeImpact.cc    完了した各energyの結果を集約
  run.sh                1 energyの実行
  run_all.sh            全energyを順次実行し、終了後に集約
  energies.sh           PS/SPSの共通energy一覧
  submit.sh             energyごとにジョブ投入
  summarize.sh          energy間の比較図を作成
  replot.sh             完了済み結果の概要ヒストグラム・PNGを更新

analysis/result/impact_studies/
  <study>/<mc-tag>/<label>/<energy>GeV/
    inputs/{data,mc}.tsv
    inputs/configuration.txt
    study.root
    summary.tsv
    paired_effect.tsv
    channels.tsv
    event_histogram_binning.tsv
    figures/
      event_energy.png
      event_energy_full_range.png
      event_energy_change.png
      event_energy_change_full_range.png
      channel_energy_change.png
      channels/layer<L>/chip<C>/channel<CH>_cell<ID>_gain_consistency.png
  <study>/<mc-tag>/<label>/energy_scan/
    energy_scan.tsv
    mean_energy_vs_beam.png
    mc_over_data_vs_beam.png
    rms_over_mean_vs_beam.png
    paired_energy_change_vs_beam.png
```

- `study`: 検証項目。今回は `hl_intercept`。現在実装されているのはこの項目のみ。
- `mc-tag`: 入力MCの条件ディレクトリ名。既定値 `threshold`。
- `label`: 実行条件を区別する名前。既定値 `default`。例 `production_v1`、`validation`。
- `energy`: 入射beam energy。イベントの再構成energyとは区別する。
- 図名は観測量を表し、検証項目・条件・beam energyは親ディレクトリで区別する。
- 比較曲線の名前は `data`、`mc_nominal`、`mc_intercept_zero`。

## 実行

ROOTの `root-config` とg++が使用できる環境で実行する。必要な実行ファイルは自動ビルド。

```bash
cd /home/murata_t/ScECAL_BeamTest/analysis/impact_studies

# 入力確認のみ（出力やジョブを作らない）
bash submit.sh --dry-run

# 100 GeV、全イベント・全channelをジョブ投入
bash submit.sh --energy 100 --label production_v1

# PS 0.5,1,2,3,4,5 GeV + SPS 10–250 GeVの計19ジョブと自動集約
bash submit.sh --energy all --summarize --label production_v1

# ジョブ終了後、完了したenergyを集約
bash summarize.sh --label production_v1

# ジョブを使わず100 GeVを実行
bash run.sh --energy 100 --label direct_v1

# ジョブを使わず全19 energyを順次実行し、比較図も作成
bash run.sh --energy all --label direct_all_v1
# 同じ処理: bash run_all.sh --label direct_all_v1

# 動作確認。詳細PNGだけ2channelに制限する（全channelの集計は維持）
bash run.sh --energy 100 --max-events 1000 \
 --plot-channels 240008,920012 --label validation
```

描画設定の変更を既存結果に反映する場合は、次を実行する。

```bash
# production_v1の完了済み全energyを再描画（100だけなら --energy 100）
bash replot.sh --energy all --label production_v1
```

再描画は概要PNG 4枚と `event_histogram_binning.tsv`、`summary.tsv` の範囲外件数を更新する。
平均・RMSなどの統計量は変えず、読み直したイベント数・平均・RMSが元の結果と一致することを検証する。
既存の `study.root` は変更せず、現在の図に対応するカウント分布を `event_histograms.root` に保存する。
元の `study.root` 内のヒストグラムは以前のbin設定なので、再描画後は後者を使用する。
MCは保存済み `paired_events` を利用する。旧結果にdataイベント総energyがない場合に限り、
元のmanifestに従ってdataの `CellID` と `Hit_Energy` を読み直す。
以降の再描画では保存したdata総energyも再利用する。rawの照合・channel profile集計・channel PNGは再実行しない。
実行中または未完了の結果は更新しない。`all` は指定labelの完了済みenergyを対象とする。

`--max-events` は各sampleについてファイルを跨いだ上限。既定値0は全イベント。
`--plot-channels` は詳細図の選択だけを変更する。既定値allで観測された全channelを作図する。
総energy、全channelのROOT profile、channel一覧TSVは常に全channelを使用する。
大量のPNGが不要なら `--no-plots` を指定できる。

PSだけは `bash submit.sh --beam ps`、SPSだけは `--beam sps`。
複数energyの指定は `--energies 5,30,100`。`--partition`、`--time`、`--mem` も指定可能。
`--energy all` はenergy指定なしと同じ一覧を使用し、`--beam ps/sps` を付ければそのbeam内の全energyとなる。
`--summarize` は指定した全解析ジョブに `afterok` 依存する集約ジョブを追加する。
解析に失敗した場合は集約を実行しない。集約は今回指定したenergyだけを使用し、
別の実行で残っている他energyを混ぜない。入力確認は全energyについて投入前に行う。
`run.sh --energy all` / `run_all.sh` はジョブ投入を行わず、実行中のマシンで順次処理する。
この環境の `bjobs` はSlurm互換コマンドで、実際の投入は `sbatch` を使用する。
`submit.sh` は全入力の確認とビルド後に投入する。ジョブの自動投入はしない。

入力既定値は既存 `channel_response` と同じ。低energy MCも `sps` 以下にある。
入力・校正定数の変更は `bash run.sh --help` を参照。
異なるbeam、定数、イベント数等で再実行する際は別の `--label` を使う。
完了済み `study.root` は上書きせず、同じ出力への同時実行も拒否する。
失敗した実行は `study.root.tmp` のままで、energy集約には入らない。
失敗時のPNG/TSVは未完了として扱う。条件や描画対象を変える再試行は別labelを推奨。

## 今回の計算

元のCalib_Hitを入力として、校正の加算切片項だけを取り除く。
`interFactorIntercept = -Intercept / Slope` とすると、channelごとに

```
offset_MeV = 0.305 * interFactorIntercept / mipFactor
LG_zero    = Hit_LG_Energy - offset_MeV
HG_zero    = Hit_HG_Energy
E_zero     = raw_HG < 2600 ? HG_zero : LG_zero
```

これは既存の校正処理で切片を0にした場合と、保持されるhitについて代数的に等価。
温度補正、傾き、MIP定数、pedestal、HG/LG切り替え、hit/event選択は既存値を使用する。
元のROOTファイルや本番calibrationのソースは変更しない。別の完全なCalib_Hitファイルを
再生成する処理ではなく、比較用ROOTに集計値とイベントごとの比較treeを保存する。

raw/calibは同じbasenameで対応させ、全処理entryのRun_Num・Event_Time・TriggerID/Event_Numを
照合する。hitはメモリcellを含む完全なCellIDで結合し、図は物理channel単位にまとめる。
entry順が違う、hit対応が曖昧、配列長が異なる場合は停止する。

MCでは、指定した定数とraw ADC・保存済み温度からHG/LG energyを再計算し、
**元の保存energyを再現できることを全処理hitで検証する**。
許容差は `1e-5 MeV + 1e-6 * abs(expected)`。Hit_EnergyのHG/LG選択もdata/MCで検証する。
定数・校正処理が異なれば停止し、成功した結果として保存しない。

定数の品質カットと代替値は現在の `ECAL_Analysis_LCIO/src/Calibration.cxx` の
`*_simulation` 関数を再現する。特に既存実装ではMIPの負のsentinel値も比のカットを
通る場合があること、および切片配列の `{1}` 初期化は先頭だけ1で他は0であることを
そのまま再現する。これらの扱いを独自に修正すると、切片以外の条件まで変わるため。
pedestalの読み込み時のfloat精度への丸めも既存calibrationに合わせる。
`channels.tsv` に使用offsetと代替値使用フラグを記録する。

## 出力の読み方

- `event_energy.png`: data、元MC、切片0 MCの校正後hit energyの総和（入射energyへの換算なし）。
  各sampleの全イベント数で正規化する。各sampleの0.1–99.9百分位と0を含む共通rangeに余白を付ける。
  元MCの `(P84-P16)/2` に対し約10 binを目安に、実用上約800 binを上限としてbin幅を決定する。
  ビームenergyをそのまま横軸のスケールには使わない。範囲外の数を図と `summary.tsv` に記録する。
- `event_energy_full_range.png`: 負の値・外れ値を含む全イベントのrangeを対数縦軸で表示する。
  主図のrange外の裾も確認できる。全範囲図は約400 binを目安とする。
- `event_energy_change.png`: 同一MCイベントにおける `E_zero - E_nominal` とその相対差。
  相対差は元の総energyが0のイベントだけ除外する。各分布に合わせてrange/bin幅を決定する。
- `event_energy_change_full_range.png`: 変化量・相対差の全範囲を対数縦軸で表示する。
- `event_histogram_binning.tsv`: 各histogramのbin数・range・bin幅・全entry数・範囲外件数。
  主要部分の図では範囲外のイベントも正規化の分母に含め、選別による分布変更は行わない。
- `channel_energy_change.png`: channelごとの総energy変化をMC全イベント数で割ったマップ。
- `*_gain_consistency.png`: 左はraw HG ADCに対する平均 `E_LG-E_HG`、右は切り替え付近の
  平均selected hit energy。縦線はraw HG=2600、エラーバーは平均値の標準誤差。
  低信号のLGノイズやHG飽和の領域まで一致することを要求する図ではない。
- `channels.tsv`: offset、MC hit数、LG選択hit数、総変化量、data/元MC/切片0 MCの
  残差平均と標準誤差。残差の集計範囲は共通のraw HG `[800,2200)` ADC。
  対象hit数0の行の平均・標準誤差0は「測定値なし」を意味する。
- `summary.tsv`: イベント数、平均、平均の標準誤差、RMS、RMS/平均、表示範囲外の数。
  統計量にはヒストグラム範囲外のイベントも含む。
- `paired_effect.tsv`: 同一イベント差の平均と標準誤差、RMS、LG hit数、元校正値の検証結果。
- `study.root`: カウント分布、全観測channelのprofile、`paired_events` tree。
  treeには入力ファイル番号・entry・イベントID・hit数・LG hit数・補正前後の総energy・差を保存。
  入力ファイル番号は `inputs/mc.tsv` の0始まりの行番号。
  新規解析ではdata総energyも `data_events/energy_GeV` に保存する。
- `energy_scan/`: 完了済みenergyだけを集約する。途中実行でも部分的な一覧を作れる。
  元MCと切片0 MCの差の誤差は同一イベント差から計算する。
  手動の `summarize.sh --energies 0.5,1,5,100` は指定した結果がすべて揃っている場合のみ集約する。

対象はlayers 0–29の接続channel（chip5はchannel0–29、他は0–35）。
追加のCoG・beam位置・シャワー選択・bad-channel除去は行わない。空イベントも含める。
ここでの **RMS/平均は分布全体の指標で、Gaussian fit等で求めた既存のenergy resolutionとは異なる**。
まず切片だけの影響をこの条件で確認し、既存resolutionとの比較には同じ選択・fitを適用する必要がある。

## 検証

```bash
make test
```

人工ROOT入力でHG=2599/2600の境界、正負の切片、定数未取得時の挙動、温度25℃、
hit並べ替え、空イベント、補正前後のイベント差、描画・集約を確認する。
元校正値の不整合を検出して停止すること、完了結果の上書きを拒否することも確認する。
`bash test_all.sh` は全19 energyの選択、集約ジョブの依存関係、途中失敗時の停止、
指定energy不足時の集約拒否を確認する。実際のジョブ投入やROOT解析は行わない。
