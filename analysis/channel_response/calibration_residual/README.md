# HG/LG校正残差と切替点の診断

保存済みの校正energyから、channelごとの `E_LG - E_HG` とpedestal差引き後HG ADCの関係を調べます。
run・温度依存、切替点の接続、使用候補定数とfallbackを追跡する解析です。
解析本体・描画はC++/ROOT、実行制御はshellです。実行にはPythonを使いません。

```bash
cd ~/ScECAL_BeamTest/analysis/channel_response/calibration_residual

# 入力を確認。書き込み・build・投入なし
bash run.sh --energy 100 --dry-run

# 最初の動作確認：全data runから各1000 eventを均等に抽出
bash run.sh --energy 100 --sample data --channels 240008,420008,920012 \
  --events-per-file 1000 --min-hits 10 --label smoke

# data全channel・全eventの数値とROOTを保存
bash submit.sh --energy 100 --sample data --plots 0 --label data_full

# 特定channelの図を全eventから作成（dataと現在のMC）
bash submit.sh --energy 100 --channels 240008,420008,920012 --label selected_channels

# 全19 energy。各energyを独立したSlurm jobで実行
bash submit.sh --energy all --sample data --plots 0 --label all_data
```

`run.sh`は端末で実行、`submit.sh`はSlurmに投入します。`submit.sh --dry-run`は入力検証と投入コマンド表示だけです。
既定は100 GeV、全channel、全event、dataとMC両方、PNGありです。
全energyはPS: 0.5,1,2,3,4,5、SPS: 10,20,30,40,50,60,70,80,100,120,150,200,250 GeV。
`--energies 30,50,100`も指定できます。ログは`../../log/channel_response/calibration_residual/`です。
`--mem`・`--time`・`--partition`・`--log-dir`は`submit.sh`で指定できます。

既存の解析・本番decode/calibration・MCの定数は変更しません。
校正energyを再選択したり、MCのHL切片を追加で差し引いたりしません。

## 入力と対応付け

- data: `/megraid01/users/data_beamtest/ECAL_data/analysed/2023/<ps|sps>/{decode,calib}/e-/<energy>GeV/`
- MC: `/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/{decode,calib}/e-/e_ssa_noHLintercept/threshold/<energy>GeV/`
- 既定定数: 同MCディレクトリの`Analysis_edit/share/`にある`pedestal2023_SPS.root`、`all_auto_muon_v4_trackfit.root`、`all_hl_electron2023_hlratio_v3.root`。

data/MCそれぞれにdecode・calibのファイルまたはディレクトリを指定できます。
ディレクトリでは同じbasenameを対応させ、欠落はエラーにします。
`--data-pedestal`・`--data-mip`・`--data-hl`、MCについても同様の引数で定数を個別に指定できます。
他のMC系列を調べる場合は`--mc-decode`と`--mc-calib`を両方指定します。

`(Run_Num, Event_Time, TriggerID/Event_Num)`を照合後、memory cellを含むfull CellIDでhitを対応させます。
イベント順が違う場合は両tree全体で一意性・同じキー集合を確認します。
重複キーは、entry数と**全entryのキー列**が一致した場合だけentry順で対応します。
calibrationにCycleIDがないため、同じキー内部だけを並べ替えた外部生成ファイルは識別できません。
full CellIDの重複、対応するraw `HitTag=1`の欠落、vector長の不一致はエラーです。
`NewTemperature`は必須で、非有限温度は除外数を記録します。
集計時にのみmemory cellを除いてphysical CellIDにまとめます。

対象は接続されたlayer 0–29、chip 0–5（chip5はchannel0–29）。
既定ではbad-channel・CoG・位置・正energyのカットを追加しません。
calibration側ですでに落とされたhitは解析に含まれません。
`--channels`はphysical CellIDのカンマ区切り、`--layers`は`4,9-10`等です。

## run・温度・統計

集計単位は **入力ファイル × 実際のRun_Num × 温度bin × physical CellID** です。
温度は保存されたhitごとの`NewTemperature`を用い、既定は0.5℃幅の`[k*0.5,(k+1)*0.5)`。
異なるrunや温度を一つの平均に混ぜません。同じrunが複数ファイルに分かれていても、ファイルを分けたまま出力します。
MCもファイル別に保存します。比較時は`files.tsv`のファイル対応と温度binを確認してください。
温度20℃が実測・再構成された値か、校正処理による欠測時の代替値かは、保存枝だけでは識別できません。

各ADC binに以下を保存します。

- hit数、平均HG ADC、残差の平均とSEM。
- 残差の16%点、中央値、84%点、中心68%幅 `q84-q16`。
- 統計不足の表示。既定は30 hit未満を`low_statistics`とし、図と中央値binのfitから外します。TSVには残します。

分位点は元の残差を全て使って昇順に並べ、位置`(N-1)*p`で線形補間するexact/type-7定義です。
分布を打ち切らず、負energy・tailも含めます。q16–q84は**分布の幅**であり、中央値の誤差ではありません。
SEMはhitが独立という仮定に基づきます。event内相関や系統誤差は含みません。
x軸は温度補正前のpedestal差引きHG ADC、y軸は保存された全校正後energyです。
HGの表示範囲外hitは`adc_outside`として記録します。これらはADC別図には含まず、切替fitやoverlap fitは各指定範囲で独立に集計します。

既定のHG ADC範囲は-100–4200、50 ADC/bin。
実際の切替位置はchannelごとに`2600 - pedestal_HG`となり、残差図に縦線を引きます。
raw HG 800–2200の範囲では残差のhit単位線形回帰の切片・傾きと、中央値binの傾きを保存します。
切片はpedestal差引きHG=0への外挿です。この範囲の外挿や傾きだけで原因を断定しないでください。
曲がりは残差図のADC依存から確認します。非線形性を一つの数値で自動判定する処理はありません。

## raw HG 2600 ADC付近の接続

既定はraw HG `[2400,2800)`を25 ADC/binで調べます。
保存されたHG・LG・選択energyそれぞれの平均/SEMと中央値/q16–q84を同じ図に描きます。
`Hit_Energy`がraw HG `<2600`でHG、`>=2600`でLGと一致するかも、`selection_mismatch`で検査します。

選択energyの段差は、切替左右を**別々に直線fitし、同じraw HG=2600へ外挿したenergyの差（右−左）**です。
左右の信号量の違いを段差として数えないため、左右の単純平均差は使いません。

- `selected_jump_OLS_MeV`: 全hitの左右OLSから求めた差。
- `selected_jump_OLS_SE_MeV`: 独立hit・線形近似を仮定した形式的標準誤差。
- `selected_jump_median_bins_MeV`: 各bin中央値を等重みでfitした差。tailへの感度比較用で、統計誤差の推定はしません。
- `paired_delta_at_switch_OLS_MeV`: 切替窓全体の同一hit `E_LG-E_HG`をfitし、2600で評価した値。

左右それぞれに`min-hits`を満たすbinが3個以上、最近接の有効binが切替から1 bin幅以内に必要です。
不足時は`switch_status=insufficient_support`とし、上記の段差値は`nan`です。0で代用しません。
左右のhit数、切替点との最近接距離も記録します。
fit窓の線形性は仮定です。`--switch-window 100`と`200`などで安定性を比較してください。
HG飽和が始まる領域の曲がりや選択biasは、単純な段差とは区別して図で確認します。

## 定数・適用期間・fallbackの監査

`inputs/calibration_provenance.tsv`にdata/MC別の定数の絶対パス、SHA-256、申告された適用期間、根拠ファイルを保存します。
既定の適用期間は`UNCONFIRMED`です。ファイル名の年や対象runの時刻を、校正定数の有効期間とみなしません。
`--data-period 'Run105-107; ...'`等で既知の適用条件を記録し、`--data-provenance PATH`で生成時ログ/設定を添付できます。
証拠ファイルの有無だけで自動的に「生成時定数確認済み」には変更しません。
実際に処理したrunとEvent_Timeの最小・最大は`files.tsv`に別途保存します。Event_Timeは保存値そのままです。

採用係数を再現した`calibration_channels.tsv`には全6300 channelについて次を保存します。

- pedestalと有無。現在の校正と同じfloat丸め。観測channelのpedestal欠落時は停止します。
- 元MIP値、χ²/NDF、採用MIP、品質棄却/欠落/sentinel=1によるfallback、その採用値が非正か。
- 元HL slope/intercept、採用する逆ratioとHG換算offset、品質棄却/欠落、ratioと切片それぞれのfallback。
- HL棄却時に切片がfallbackされず初期値0のまま残るchannel。

現行`Calibration.cxx`の`*_simulation`関数はdataにも使用されています。
この抽出方針を再現し、MIP品質cutは10µm: χ²/NDF≤2、15µm: ≤1.7、HL slopeは0.02–0.05。
MIP代替平均には現行どおり、品質cutを通った非正MPVも含みます。ここで定数を修正することはしません。
HL切片の初期値`array={1}`と`Init()`で切片を再初期化しない挙動も記録対象です。
pedestal未収録channelは監査表に`missing_production_zero`と記録しますが、この解析では観測された場合にゼロで補って進みません。

現行校正ソース・header・global_configが存在する場合は`inputs/current_*`に参考コピーを保存します。
これは生成当時の履歴ではありません。ソース・実行バイナリ・定数のhash、入力一覧、コマンド、ROOT versionも保存します。
指定定数と保存温度から再計算したHG/LGと保存energyの差を`channels.tsv`に記録します。
このclosureがずれる場合、定数または適用式の対応を確認してからfallbackと残差の因果関係を解釈してください。
異なる期間に異なる定数を使う場合は、対象ファイル・定数を明示して別`--label`で実行します。

## 出力と実行上の注意

```text
analysis/result/channel_response/calibration_residual/<mc-tag>/<label>/<energy>GeV/
  inputs/                         入力一覧・設定・定数/ソースの由来
  data/ または mc/
    residuals.root                各群の平均/分位点graph、切替結果・定義
    residual_bins.tsv             ADC別残差統計
    switch_bins.tsv               切替付近のHG/LG/選択energy統計
    channels.tsv                  群ごとの傾き・段差・closure・除外数
    calibration_channels.tsv      採用候補定数とfallback（全channel）
    files.tsv                     ファイル、run、時刻範囲、照合方式
    figures/file_N/run_R/temperature_bin_K/channel_ID.png
  data.log / mc.log
  COMPLETE                        指定sampleの処理が全て成功した印
```

PNGは群ごとに6 panel：残差の全範囲・overlap拡大・切替付近拡大、切替energy平均、切替energy中央値、選択energyの左右fit。
全パネルの描画点が0個になる群は、canvas生成とPNG保存を省略します。
判定は`--min-hits`適用後の点数で行い、残差または切替energyのどれかに点があればPNGを保存します。
省略する群もROOT・TSVの統計は保存します。既存の空PNGは削除しません。
ROOT graphも統計条件を満たすbinのみ含みます。低統計binを含む全数値はTSVにあります。
既存出力ディレクトリの上書きは拒否します。再実行は別labelを使ってください。
PNGの別プロセス再描画は実装していません。`--plots 0`は数値集計だけにし、図が必要なら別labelでchannelを絞って再実行します。

PNGは既定で`../adc_energy/FastPng.hh`による高速出力を使用します。追加オプションは不要です。
点・bin・誤差・画像寸法を維持し、丸いmarkerの描画を小領域で処理します。
PNGの可逆圧縮は速度を優先するため、従来の`SaveAs()`よりファイル容量が増えます。
100 GeV dataの1000 event・10 channelを指定した小標本では、生成された9枚すべてのpixelと統計表が一致し、
PNG描画・保存時間は1.29秒から0.82秒でした。全量解析の短縮率を保証する値ではありません。
ログは1000枚保存ごとの進捗と、入力ファイルごとの`PNG: saved=... skipped_empty=... raster_seconds=... encode_seconds=...`を出します。
描画とPNG圧縮・保存の時間を分けて記録し、1枚ごとの保存ログは出しません。
更新した実行ファイルは新しく起動するプロセスから使用されます。起動済みプロセスは旧版で処理を続けます。

分位点を正確に求めるため、各入力ファイルのhit残差を一時的にメモリに保持します。
残差値だけで対象hitあたり8 byte、切替窓内はHG/LG/選択energyの24 byteが追加され、map/vector等の管理領域も必要です。
ファイルの処理完了ごとに解放します。大きいdataは`--layers`や`--channels`で分割し、`--mem`を調整できます。
全channel・全MCファイルでPNGを作ると大量になるため、全体集計は`--plots 0`、図はchannelを指定する実行例を推奨します。
`--events-per-file`はファイル全体を均等に覆う中点sampleです。0なら全event。
イベント照合のためのキー列は、sample数によらず全entryを読みます。
`--max-files`は辞書順で先頭Nファイルに制限し、runを落とす可能性があるため、run比較時は原則0にしてください。

## 検証

```bash
make
make check   # テスト時だけpython3.12 + PyROOTが必要
```

人工ROOTで既知のoffset/slope・切替段差を注入し、run×温度分離、exact分位点、中央値と平均の段差、
memory cellを含むhit照合、並べ替えイベント、重複キーの許容/拒否、対応hit欠落の検出、
fallback、片側しかない場合の統計不足判定、ROOTとPNGの生成、全点なしのPNG省略と残差のみ/切替のみのPNG保存を検査します。
テスト用ファイルは`/tmp`の一時領域に作り、終了時に削除します。

残差がほぼ一定ならpedestal/HL切片、ADCに比例するならHL ratio、高ADCで曲がるなら非線形性/HG飽和、
温度やrunで変わるなら温度依存/適用条件が調査候補です。原因の確定や校正定数の自動更新は行いません。
