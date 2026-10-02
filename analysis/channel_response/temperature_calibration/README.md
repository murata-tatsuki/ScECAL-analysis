# PSの温度・HG/LG校正の検証

同じphysical channelのHG-only MIP MPVを温度別・run別に測り、温度補正の検証後にHG–LG接続とpedestalを確認する解析です。既存の`mip_closure/CollectMuon.cc`の飛跡選択、`mip_closure/mip_tools.py`のLandau–Gaussian fit、`calibration_residual/CalibrationAudit.hh`の定数抽出・fallbackを参照・再利用しています。

## 実行

ROOT / PyROOT、NumPy、C++17以上が必要です。この環境では`root-config`と`python3.12`を使用します。追加のPythonパッケージは不要です。

```bash
cd ~/ScECAL_BeamTest/analysis/channel_response/temperature_calibration
make
make check

# PS muon 10 GeVをMIP検証、electron 0.5–5 GeVをHG/LG検証に使用
mkdir -p ../../log/channel_response/temperature_calibration/inputs
python3.12 make_manifest.py --output ../../log/channel_response/temperature_calibration/inputs/ps_temperature_inputs.tsv

params=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share
bash run.sh \
  --manifest ../../log/channel_response/temperature_calibration/inputs/ps_temperature_inputs.tsv \
  --pedestal "$params/pedestal2023_SPS.root" \
  --mip "$params/all_auto_muon_v4_trackfit.root" \
  --hl "$params/all_hl_electron2023_hlratio_v3.root" \
  --threshold ../../result/threshold/threshold.root \
  --output ../../result/channel_response/temperature_calibration/ps_default
```

Slurmではmanifestを上記の共有領域に保存してください。`/tmp`はノードごとに異なるため、投入元で存在していても計算ノードで読み取れない場合があります。既存manifestを再利用する場合は生成コマンドを省略してください。

上記は既存解析と同じ**候補定数**を渡す例です。SPSと名付けられたpedestalをPSの正しい定数と認定するものではありません。実際に使用した定数とproduction記録が分かる場合はそれを指定してください。校正期間・定数が異なるデータは別manifest・別outputに分けます。

`bash submit.sh`に同じ引数を渡すと共通の`../execute.sl`経由でSlurmに投入します。`--dry-run`は入力パス・選択設定の確認のみで、ROOTデータの完全検査やジョブ投入はしません。`--events-per-file 500 --channels 2300008,2500008,2700008`で短い動作確認ができます。これは**raw event**の全期間等間隔サンプルであり、500個のmuon trackを保証しません。測定用は0（全event、既定）を使います。

`--output`は未作成directoryを指定してください。成功時だけ`COMPLETE`を作ります。失敗時の詳細は`file_N/collector.log`に残ります。

PNGは既定で共通の`../adc_energy/FastPng.hh`をPyROOTから呼んで高速保存します。円形マーカーの描画範囲を小さくし、可逆PNG圧縮をlevel 1に設定しています。通常の`SaveAs()`よりファイルサイズが大きくなる場合があります。追加オプションは不要です。channelごとの描画・エンコード時間とファイルサイズを`png_timings.tsv`へ保存し、100枚ごとに進捗、最後に合計時間をログに出します。`--no-plots`なら描画・PNG保存を省略します。

## 入力manifest

header付きTAB区切りで`kind, raw, calib, track`を指定します。相対パスはmanifestのdirectory基準です。

| kind | 入力 | 用途 |
|---|---|---|
| `mip` | Raw_Hit、Calib_Hit、T_Eventの3ファイル | 飛跡選択したHGのMIP MPV |
| `response` | Raw_Hit、Calib_Hit、trackは`-` | 同一hitのHG–LG残差 |
| `pedestal` | Raw_Hit、残りは`-` | 専用pedestal runなどのHitTag=0 |

すべてのkindでrawの`HitTag=0`をpedestal診断に使います。同じrawを複数行に書くと二重計数を防ぐためエラーになります。MIP定数の取得時温度を直接監査するには、その定数を作成したmuon runも入力し、production記録と照合してください。現在のPSデータの温度だけから定数の測定温度は推定しません。

raw/calibは全eventの`Run_Num, Event_Time, TriggerID/Event_Num`を照合し、順番が違う場合は一意なkeyで結合します。繰返しkeyは全列の順序一致を確認できた場合だけ許容します。trackはRun_Numが無いため、1 rawファイル内で`Event_Time, TriggerID`の一意性を要求します。hitはmemory cellを含むfull CellIDで結合した後にphysical channelへまとめます。破損からrecoverされたROOT、曖昧な対応、vector長不一致はエラーです。

## 温度と式

現行の`ECAL_Analysis_LCIO/src/Calibration.cxx`と同じ候補式は次の通りです。

```text
F(T) = 1 - alpha * (T - Tref)
HG = (rawHG - pedHG) * F(T) * 0.305 / MIP
LG = ((rawLG - pedLG) * F(T) * HLgain + HLoffset) * 0.305 / MIP
Tref = 20 C
alpha = 1.6/135 (layer 4–27), 3.5/230 (layer 0–3,28–29)
```

29℃では係数はそれぞれ約0.89333と0.86304です。符号を正しいと仮定して合否を決めるのではなく、まずこの式と保存branchの一致を検証します。`--reference-temperature`、`--coefficient-10`、`--coefficient-15`で仮説を変更できます。`--compare-models`は`2-F`（符号反転）、`1/F`（逆数）でも同じhit集合をfitします。これらは診断用の代替式で、推奨校正式ではありません。

`--mip-measurement-temperature`はMIP取得時の実測温度、`--mip-constant-temperature`は取得時補正なども含めてその定数が表す温度です。両者を区別し、`--mip-temperature-evidence`でproduction記録を保存できます。未指定は`unconfirmed`、指定値の一致／不一致も`declared_match/mismatch`とし、自動的なprovenance証明とはしません。

信号の温度軸には保存された`NewTemperature`を使います。ただしrawにそのlayerの16センサーの有効値（0<T<60℃）が無いhitは、温度依存fitから除きます。校正時の欠測代入20℃を実測値として扱わないためです。独立に再構成した温度との差も保存します。pedestalは保存calibにhitが無いのでrawセンサーとstrip位置から逆距離重みで再構成します。`TemperatureGeometry.hh`は現行`EBUdecode.cxx`の配置のsnapshotです。元コードの`localTemperature`未初期化は再現せず0で初期化します。これは読取り・診断側だけの処置で、productionコードの修正ではありません。

## 選択と統計

既存の定数抽出はchi2/NDFのcutのみで、MIP値`-10`などの非正値が採用される場合があります。本解析は定数の読取り方と`calibration_channels.tsv`を維持し、そのようなchannelをMIP／HG–LG／closure評価から除外して他のchannelの解析を続けます。除外理由・file/run/温度bin・channel・hit数・採用MIP値を`excluded_channels.tsv`に記録します。定数の自動置換は行いません。MIPが不正でもpedestal診断は実行します。pedestal定数自体が無い場合はrawの平均・SEMを残し、基準との差はNaN、statusは`missing_pedestal`にします。

除外数はログと`summary.json`にも保存します。`COMPLETE`は処理の完了を表し、全channelの校正が有効であるという意味ではありません。全信号が除外された場合は`signal_status=no_valid_signal_constants`になります。途中終了した出力は削除せず、再実行には新しい`--output`（例：`ps_fast_png_v2`）を指定してください。

- MIP: 各投影7layer以上、投影角<0.05 rad、保存objective/NDFの平方根<3 mm、strip内を全通過するtrackを要求。`hitCellnew`の負値は欠測layer。HGのみ、raw HG<2200（`--hg-max`、最大2600）、斜入射は`1/sqrt(1+sx²+sy²)`で補正します。LG・補正後energyに基づく選択はしません。
- MPV: channel × file × run ×温度binで、raw／corrected／保存HGを同じhit集合、共通fit区間のLandau–Gaussian畳込みでfit。値はLandau成分のMPVで、畳込み分布の最大位置とは異なります。`threshold + 3 sigma`と校正時の10 ADC cutより上、HG上限より下の共通区間を使います。thresholdはpedestal差引後ADCの`expErfThre`（既存mip_closureと同じ定義）が必要です。chip平均のthreshold fallbackは既存helperを踏襲しています。
- HG–LG: raw HGの線形域800–2200 ADCを100 ADC幅に分け、温度・run・channelを維持して`LG-HG`を集計。保存energy、候補式、補正なし、切片も温度補正する順序変更を並べます。平均・SEM・中央値・16/84%点を出力し、異なるADC binを同じ傾きfitに混ぜません。
- pedestal: `HitTag=0`のHG/LG平均とSEM、候補pedestalとの差をrun・温度・channel・**memory cell別**に集計。hit選別のない専用ランの無条件pedestal fitとは異なり、非hit条件によるバイアスや信号混入があり得ます。記録の無いchannelのpedestalは測れません。
- 温度傾き: 重み付き直線で`within_run`、`pooled_runs`、run別切片を持つ`run_adjusted`を出します。3点未満、温度幅<0.5℃、fit不成立は未判定です。runごとにほぼ一定温度しかなければ`run_temperature_confounded`となり、run差と温度差を分離できません。run内温度幅が無いデータを増やすだけでは解決しません。

誤差は独立pointを仮定した統計誤差です。reduced chi2>1なら傾き誤差を拡大します。fit間・event内相関、温度測定誤差、定数の系統誤差は含みません。raw/correctedは同一hitなので両者の傾き差の有意度を独立誤差から計算しないでください。飛跡の位置・角度分布がrunで違う場合は、狭い選択条件だけで完全に除去できません。合否は自動判定しません。

## 出力と読む順序

| ファイル | 確認事項 |
|---|---|
| `reference_audit.tsv` | 定数取得温度・定数が表す温度・基準温度・29℃の係数 |
| `calibration_channels.tsv` | 採用MIP/HL/pedestal、品質cut、fallback |
| `excluded_channels.tsv` | 不正なMIP／欠測pedestalにより信号評価から除外したchannel・hit数・理由 |
| `closure.tsv` | 現行20℃式と保存HG/LGの差、仮説式との差、保存温度とセンサー再構成差 |
| `mip_fits.tsv`, `mip_spectra.root` | 同じchannelの補正前後MPV・fit誤差・失敗理由・histogram/model |
| `temperature_slopes.tsv` | MPV、固定ADC binの残差、memory別pedestalの温度傾き |
| `residual_bins.tsv` | HG–LGの温度/run/ADC別残差と補正順序比較 |
| `pedestal.tsv` | HG/LGのrun・温度・memory別pedestal |
| `trends.root`, `figures/channel_*.png` | MPV・HG–LG残差・HG/LG pedestalそれぞれの温度軸／run軸 |
| `png_timings.tsv` | channel別のPNG描画時間・エンコード時間・ファイルサイズ |
| `counters.tsv`, `file_N/*` | event数、選択track数、欠測温度、結合方式、flat ROOT |
| `inputs/*`, `summary.json`, `COMPLETE` | 入力・設定・コード/定数hash・処理完了記録 |

まず保存branchとのclosureと温度の出所を確認し、HG-only MIPの傾きが補正後に残るか見ます。残るなら基準温度・係数・符号を調べます。HGが安定してもHG–LG残差が残る場合はLG pedestal、HL ratio・切片、補正順序を確認します。補正を外してdata/MCが近づくことだけでは補正の誤りと判断しません。本解析はdata内部の検証であり、MCとの一致を判定基準に使用しません。

`make check`は温度依存・run差・pedestal変動を持つ人工ROOTで、補正後MPVの回復、残差/ pedestal傾き、温度欠測、event並替え、重複hitの拒否、既存出力の保護を検証します。
