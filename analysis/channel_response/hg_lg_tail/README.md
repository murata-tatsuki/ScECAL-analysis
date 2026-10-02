# HG–LG低側tailとイベント内の同時発生

HG–LGの低側tailがどのchannelにあり、どのeventに集中し、負のLG energyや総energyとどう関連するかをdata/MCで調べるC++/ROOT解析です。hit集計とevent集計に**同じchannel集合・同じ外れ値判定**を用います。全channelを対象とし、既定の実行energyは100 GeVです。

## 実行

Slurm実行wrapperは共通の`../execute.sl`、C++・ROOTビルド設定は`../root.mk`を使用します。

```bash
cd ~/ScECAL_BeamTest/analysis/channel_response/hg_lg_tail

# 入力と保存先の確認のみ。ファイル作成・ジョブ投入なし
bash submit.sh --energy 100 --dry-run

# 100 GeVをジョブ投入：全channel、全event、全入力ファイル
bash submit.sh --energy 100

# 手元で実行する場合
bash run.sh --energy 100

# 必要になったときだけ全energyを投入
bash submit.sh --energy all

# energyを複数指定
bash submit.sh --energies 5,100 --label comparison_v1
```

Slurmの`sbatch`を使用します。`bjobs`または`squeue -u "$USER"`で確認できます。1 energyにつき1ジョブ・1 CPUです。各ジョブが集計とPNG保存まで実行するので、異なるenergyはスケジューラが許す範囲で並列に処理されます。同じenergy内は順次処理します。ログは`analysis/log/channel_response/hg_lg_tail/`、投入したjob IDはその中の`submissions_*.tsv`に保存します。`--partition`、`--mem`、`--time`、`--log-dir`も指定できます。指定しないメモリ・時間制限はクラスタの既定値です。

全energyはPSの0.5, 1, 2, 3, 4, 5 GeVと、SPSの10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200, 250 GeVです。dataのPS/SPSをenergyから自動選択します。MCは既存の配置どおり、PS相当energyも`Result_MC/{decode,calib}/e-/sps/<mc-tag>/`から読みます。

少量の診断例：

```bash
bash run.sh --energy 100 --max-files 3 --events-per-file 2000 --label sample_3files
```

`--events-per-file N`は各ファイルを10区間に分けて決定的に抽出します。`--max-events N`はdata・MCそれぞれについて、全ファイル合計の上限です。これらと`--max-files`の既定値0は制限なしです。少量サンプルでは共通channel数が減るので、必要なら`--min-hits`を明示的に調整してください。どのファイルから何eventを読んだかは保存されます。

## 判定と集計

pedestal差し引き後を`H=HG_raw−pedHG`、`L=LG_raw−pedLG`とします。校正のHG換算係数を`A=1/Slope`、`B=−Intercept/Slope`として、残差は

```
r = L − H/A + B/A
z = (r − channel・sampleごとのrの中央値) / σ
σ = sqrt(pedSigmaLG² + (pedSigmaHG/A)² + 1/12)
```

です。dataとMCそれぞれのchannel中央値を引くため、旧MCのdigitizationに含まれていない切片など、**一定のoffsetを除いたtail**を比較します。中央値自体は`channels.tsv`に残します。σはpedestalと量子化からの基準幅で、tailを含むRMSで割り直しません。「5σ」はこの基準に対する診断指標です。信号依存の傾き差・非線形性はこの処理では取り除きません。

- 対象は既存comparisonのbad-channel maskを除いた、MIP係数が正の接続channel。
- HGの比較範囲は**raw HG ADC**の`800 <= HG_raw < 2200`。
- data・MCの両方に比較範囲のhitが100以上あるchannelだけを共通集合に採用。
- `z < −5`を低側tail、`z > +5`を高側tailとして保存。
- 1 eventで低側tailが**異なる物理channelに3個以上**あれば`tail_rich`。同じchannelの異なるmemory cellはhit数には加算しますが、channel数は1と数えます。

これらは`--hg-min`、`--hg-max`、`--min-hits`、`--nsigma`、`--min-tail-channels`で変更できます。共通channel集合と中央値は、event選別前の読み込んだ全eventから決めます。channel集計も選別前です。event集計は選別前後の両方を保存します。従来の調査ではhit側とevent側の判定が異なっていたため、以前の「112 event」と同じ件数を再現する仕様ではありません。

イベント選別は`--selection`で指定します。

- `all`：event選別なし。
- `cog20`：layer 9と10の両方で、正のgood energyによるCoGがx・yとも±20 mm以内。
- `shower-cog20`（既定）：上記に加えて、event energyがMC参照値の0.5〜1.5倍。

MC参照値は読み込んだMCのgood event energy中央値です。`--reference-energy GEV`で固定値も指定できます。参照値はROOTに保存します。

総energyはgood channelの保存済み`Hit_Energy >= 0`を合計し、GeVで表示します。負のLG採用energyは総energyに加算せず、個数と符号付き合計を別に記録します。HG/LG採用は既存校正と同じraw HG 2600 ADCで切り替えます。LG負energyは比較範囲外も含む全good channelが対象なので、同一eventにある低側tailとの関連を確認できます。tail判定のHG範囲はHG採用領域にあるため、tail hit自体がLG採用されるという意味ではありません。

`tail_rich`を除いたenergy分布・RMS/meanは原因を調べるための診断です。正式なevent除去条件や補正を導入する処理ではありません。

## 保存先と図

既存解析と同じ`analysis/result/`以下に保存します。

```
analysis/result/channel_response/hg_lg_tail/<mc-tag>/<label>/<energy>GeV/
  inputs/                     入力manifest、設定、校正ファイル等のchecksum、maskの写し
  run.log
  tail_study.root              event・tail hitのTTree、histogram、参照値
  channels.tsv                channelごとの中央値・幅・tail率計算用の件数
  events.tsv                  event ID、energy、tail多重度、負のLG hit数等
  summary.tsv                 選別前後・tail-rich/poorの件数、mean、RMS、RMS/mean
  tail_rates.tsv              共通channel数、tail件数、channel構成をそろえたtail率
  validation.tsv              ADC→energy照合とhit/event集計の一致確認
  histogram_ranges.tsv        energy図の範囲・bin数・underflow/overflow
  figures/
    hg_lg_tail_overview.png
    hg_lg_tail_channel_map.png
    event_tail_correlations.png
    event_energy.png
    event_energy_full_range.png
```

既定は`mc-tag=threshold`、`label=default`です。別条件での実行には`--label`を変更します。完了済みの`tail_study.root`がある保存先には上書きしません。同一保存先の同時実行もロックで防止します。`--result-dir PATH`で保存先の基点を変更できます。

| PNG | 内容 |
|---|---|
| `hg_lg_tail_overview` | 左：中央値を引いたzのdata/MC比較。channelごとの重みを`min(Ndata,Nmc)`にそろえる。右：選別後eventの低側tail channel数。ともに規格化・対数縦軸。 |
| `hg_lg_tail_channel_map` | data/MCのlayer・channel別の低側tail率。両者の色範囲は共通。共通集合外も表示値0になるので、`channels.tsv`の`common`列で識別。 |
| `event_tail_correlations` | 左列data、右列MC。上段：低側tail channel数と総energy。下段：低側tail channel数と負energyのLG hit数。選別後event、色はevent件数、data/MCで軸・色範囲共通。 |
| `event_energy` | 選別後のdata・MC、およびdataのtail-rich event除去後。各分布を規格化。表示範囲は分位点から決定。 |
| `event_energy_full_range` | 上と同じ比較を全energy範囲・対数縦軸で表示。 |

選別後eventがdata/MCとも0の場合、event energyと相関の図は生成しません。統計量はhistogramの範囲・bin幅に依存しない方法で計算し、表示範囲外も含めます。ROOTのhistogramは規格化前の内容です。PDFは生成しません。

ROOTの`events`は全読み込みeventを保存し、`selected`と`tail_rich`、full CellIDのtail配列を持ちます。`tail_hits`は共通channelの全低側・高側tailを保存します。`sample=0`がdata、`1`がMC。`file_index`・`entry`を`inputs/{data,mc}.tsv`に対応させて元eventに戻れます。Run/Trigger/BCID/GainTag、raw HG/LG、校正後HG/LG energy、zも記録します。`event_index`は`events`のentry番号です。

## 入力、整合性確認、処理量

入力・校正パラメータは既存解析と同じ既定値を使用します。変更する場合は`bash run.sh --help`にある`--data-decode`、`--data-calib`、`--mc-decode`、`--mc-calib`、`--pedestal`、`--hl`、`--mip`を指定します。decode/calibにディレクトリを指定すると同名ROOTを対応づけます。個々のROOTファイルを指定することもできます。maskを使わない診断は`--bad-channel-source none`です。

イベント数・Run/Event_Time/Triggerの一致、full CellIDの対応、重複、保存済みHG/LG energyとraw ADCから再計算したenergyの一致を検査します。使用する校正式・係数のfallbackは`../../impact_studies/StudyCommon.hh`を共用しますが、impact studyの実行結果は不要です。校正条件や入力の組合せが合わない場合は停止します。decode→calibrationのhit残存率は解析しません。

元のROOTは必要branchだけを各sampleにつき1回読みます。中央値の確定後にhitとeventを同じ定義で分類するため、比較範囲内のhitを一時ROOTに保存し、そこを再読込します。中央値は近似せず、残差をsampleごとにメモリに保持して算出します。比較範囲内のhitがN個なら残差の数値だけで約8N byte、加えてvectorの余剰容量等が必要です。全event実行ではI/O・メモリ・一時領域が増えます。重いchannelごとのPNG一括描画は行いません。

一時ROOTは各energyの結果ディレクトリ内の`.hg_lg_tail.XXXXXX.tmp/events.root`に保存します。
既定の保存先は`analysis/result/channel_response/hg_lg_tail/<mc-tag>/<label>/<energy>GeV/`です。
正常終了・エラー終了・HUP/INT/TERMによる終了時に、この一時ディレクトリごと削除します。
完成したROOT・PNG・集計表は結果ディレクトリに残します。`--scratch-root PATH`で一時ディレクトリの親を変更できます。
SIGKILLやノード停止の場合は終了処理が動かないため、一時ファイルが残ることがあります。

## 実装時の検証

既知のtail・切片差を入れた合成ROOTで、10 tail hitと2 tail-rich eventの検出、物理channel重複排除、負LGとの対応、event不一致の拒否、一時cacheの削除、PNGのみの保存を確認しました。100 GeVのdata・MC各200 eventでも、保存済み校正energyとの一致、共通740 channelのhit/event集計一致、5枚のPNG出力を確認しました。

PSを含む全19 energyは入力確認と投入コマンドのdry-runのみ実施しました。実装時の検証用コード・ROOT・図・ログは削除済みです。全energyの解析やジョブ投入は実施していません。

## strip実座標のlayer別tail率マップ

既存の `hg_lg_tail_channel_map.png`（layer×channel番号）はそのまま保存し、通常実行時に次の4枚も追加保存します。合計9枚のPNGになります。

```text
figures/hg_lg_tail_layer_map_data_layers00_15.png
figures/hg_lg_tail_layer_map_data_layers16_29.png
figures/hg_lg_tail_layer_map_mc_layers00_15.png
figures/hg_lg_tail_layer_map_mc_layers16_29.png
hg_lg_tail_layer_maps.root
```

各canvasは4×4で、左上からlayer番号順です。2枚目はlayer 16–29の14枚と凡例です。横軸x・縦軸yはmm、色は選別前hitの低側tail率（0.01 = 1%）。同じenergyのdata/MC・全30 layerで色の最小値0と最大値を共通にします。白は対象channelのtail率0、灰色はbad channelや統計不足などで共通集合に入らないstripです。

座標変換は `analysis/calibration_parameters/EBUdecode.cxx` を共用し、既存の校正マップと同じ5.3 mm×45.4 mmのpitch binを使います。偶数layerは5×42、奇数layerは42×5 binです。各binは1本のstripに対応し、strip間の隙間もpitchに含みます。ROOTにはdata/MCのlayer別TH2D（計60個）、共通channelのmask（30個）、入射energy・共通色上限・geometryの説明を保存します。元の `tail_study.root` の集計は変更しません。

すでに完了した結果には、event再解析なしで追加図だけを保存できます。

```bash
cd ~/ScECAL_BeamTest/analysis/channel_response/hg_lg_tail
bash plot_layers.sh --energy 50
# 別の保存条件の場合
bash plot_layers.sh --energy 50 --mc-tag threshold --label comparison_v1
```

`plot_layers.sh`は `tail_study.root` と `channels.tsv` を読み、追加4枚と追加ROOTだけを生成・更新します。元の5枚は変更しません。保存先の基点は `--result-dir PATH` で指定できます。未完了の結果には実行できません。

layer別マップの表示範囲はx・yとも±120 mmです。各padのframeがピクセル上でも正方形になるよう余白を調整し、x・yの1 mmを同じ長さで表示します。stripのbin境界・tail率は変更しません。
