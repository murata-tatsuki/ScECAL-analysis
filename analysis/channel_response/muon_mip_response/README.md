# MuonのMIP応答再現

100 GeV muonについて、MC truth → pedestal差引きHG ADC → 通過長・温度補正 → 保存校正値の対応をchannelごとに調べるC++/shell解析です。dataとMCで同じtrack推定・選択・fitを使います。図はPNGのみ、分布・fit関数はROOT、数値と未判定理由はTSVに保存します。

## 実行

```bash
cd ~/ScECAL_BeamTest/analysis/channel_response/muon_mip_response

# 入力パスと設定の確認。出力・ビルド・ジョブ投入は行わない。
bash run.sh --dry-run

# 全muon run、全event、全channelをSlurmへ投入する（既定で高速版・4 process・4 CPU）。
bash submit.sh

# 同じ解析をローカルで実行する。
bash run.sh

# 保存済みの分布・fitから全種類の既存PNGだけ更新する。
bash submit.sh --plots-only

# MC内部のtruth→ADC再現を先に調べる。
bash submit.sh --sample mc --label mc_only

# 別の条件を保存する例。
bash submit.sh --runs 44 --label run44
```

既定は全run・全event・全channelです。`--max-events N` は各ファイルから等間隔にN eventを抽出し、0は全eventです。`--channels CELLID,...` はfit・channel図だけを限定し、収集・histogram保存は全channelを保持します。例えば`400000`はlayer4/chip0/channel0の物理CellIDです。

`--dependency-fits 0` は位置・温度・run別MPVの追加fitを省略します。`--collect-only` は入力監査とhit収集まで行います。全オプションは `bash run.sh --help`、ジョブ資源は `bash submit.sh --help` を参照してください。C++コンパイラ、ROOT、ASImage、flock、ジョブ投入時はSlurmが必要です。

既存の`flat.root`または`COMPLETE`があるlabelには書き込みません。異なる条件を残すときは`--label`を指定してください。ジョブは入力の集計後にchannelごとのfitと描画を複数processへ分配し、最後にROOT・TSVを統合します。全energyのelectron解析や本番digitization・定数更新は実行しません。

## 処理時間と並列実行

既定は `--workers 4`、Slurmの割当ては4 CPU・16 GBです。workerはchannelごとに担当を分け、同じchannelのdata・MC、全段階、依存性のfitとPNGを同じprocessで処理します。共有ROOTへ同時書き込みはせず、終了後に統合します。収集・geometry matching・histogram作成・最後の統合は親processで行うため、全体がCPU数に比例して速くなるとは限りません。

`--fit-cache 1` が既定です。binの値・誤差・axis・entriesとfit範囲・規格化・有効統計下限・品質条件が完全に一致する場合だけ、同じprocess内でfit結果を再利用します。名前が違うだけの分布でも利用できます。近い分布や異なるfit条件は再利用しません。cacheはworkerごとに32 MiBを上限とし、`--fit-cache-mb N` で変更、`--fit-cache 0` で無効化できます。

Crossingsはgeometry集計用と応答集計用で必要なbranchだけ読み、依存性のhistogramも応答集計時に作成して追加の全hit走査を省きます。TruthHits・Events・校正定数も必要なbranchだけ読みます。依存性histogramを先に保持する分、メモリ使用量は増えます。

高速化した処理が通常の既定動作です。`bash submit.sh` または `bash run.sh` だけで並列化・fit cache・必要branchのみの読み込みが有効になり、`analysis/result/channel_response/muon_mip_response/default/` に保存します。高速化のための別labelや追加directoryは不要です。

`--reuse-flat` は収集完了と保存された設定を確認します。run・event制限・track条件などの収集条件は元の設定を引き継ぎ、異なる条件を明示するとエラーにします。元のflat.rootはread-onlyで使い、出力側からシンボリックリンクで参照します。既に実行中のジョブは変更前のコードで継続し、この変更は次回実行から適用されます。

`--plots-only` は同じlabelの既存 `spectra.root`・`fits.tsv`・`dependency_fits.tsv` をread-onlyで使い、既存のresponse・MPV profiles・mean profiles・channel mapのPNGを更新します。収集・再fitは行わず、数値のROOT・TSVは更新しません。`--channels CELLID,...` で更新するchannelを絞れます。channel mapは全channel指定時だけ更新します。描画コード・コマンド・入力の前後checksumは既存の `inputs/plot_*` に保存し、run.logへ追記します。fit条件の変更optionはこのmodeでは受け付けません。

workerの一時出力はプロジェクトの `.scratch/muon_mip_response_runs/` にまとめ、統合成功後に今回の分だけ削除します。失敗時は保存先をrun.logへ表示して保持します。`performance.tsv` にworkerごとのfit要求・cache hit・実際のROOT fit呼び出し回数を保存します。`dependency_prepared_bins` は親で用意した総bin数を各行へ記載した値で、worker間で合算しません。

## 入力

既定のdataは `/home/murata_t/data_beamtest/ECAL_data/analysed/2023/sps/{decode,calib}/mu-/100GeV/`。MCは `/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/` 以下の

- `{decode,calib}/mu-/mu_track/threshold/100GeV/ECAL_Run<run>.root`
- `generate/mu-/mu_track/100GeV/ECAL_Run<run>.root`

を使います。`--data-base`、`--mc-base`、`--mc-tag`で変更できます。保存済みtrackFit3Dは不要です。

任意入力は `--input-manifest PATH` で指定できます。ヘッダーなし、絶対パス、タブ区切りの6列です。

```text
sample(0=data/1=MC)  run  decode.root  calib.root  -  truth.root(dataは-)
```

manifestにも`--sample`と`--runs`を適用します。MCの比較用runは生成元runを示すファイル名/manifestの値です。MC ROOT内の`Run_Num`が0でも、rawとcalibの保存番号は照合し、`source_run`として残します。dataはmanifestと保存run番号の一致を要求します。

MIP・pedestal・HG/LG定数は現在の`Analysis_edit/share/`を既定に使い、data用・MC校正用・MC生成用MIPを別々に指定できます。2026-10-10更新の実装と同じ有限値・正MPV・NDF・χ²/NDFの採用条件、接続channelの平均補完を再現します。thresholdは `analysis/result/threshold/threshold.root` の `expErfThre`、欠損値は現行MCと同じchip内補完を記録します。

入力ROOTのtree・entry数・復旧状態・更新時刻、定数hash、現在のproducerソースを保存します。現在のソースの控えは、入力ROOT生成時の実行ソースを証明するものではありません。MIP定数の基準温度・取得時の選択・通過長補正の由来も、入力監査だけでは確定できません。

欠損・読めないtree・復旧ROOTがあるrunは、data/MCの対ごと除外して監査表に残します。`--allow-recovered`で復旧ROOTを含められますが診断扱いです。

## Trackと検出割合

Raw_HitのHitTag=1のstrip中心から、奇数層のXと偶数層のYを推定します。信号量で重みを付けません。各対象layer全体を除外してからtrackを求め、既定では各投影7層以上、角度0.05 rad以内、残差RMS 3 mm以内を要求します。

45×5×2 mmのstrip内を全厚通過する予測を使い、長手1 mm・幅0.25 mmの余白を設けます。対象layerにhitがない通過予測も`Crossings`へ保存します。重複memory cellや対応の曖昧さは応答fitから除外し、別に数えます。

`efficiency.tsv`は **他層でtrackが求まる条件下の、予測stripでの一意なhit取得割合** です。track位置の誤差・隣接stripへの割当ても含み、hardware単独の取得効率ではありません。marginやtrack条件を変えて安定性を確認してください。

## 応答比較とfit

同一channelで次の6段階を保存します。

| stage | 量 |
|---|---|
| `pedestal_adc` | HG − pedestal |
| `path_adc` | pedestal差引きHG × 2 mm / track通過長 |
| `temperature_adc` | path ADC × 現行20℃補正係数 |
| `saved_hg_mip` | 保存HG energy / 0.305 MeV × 2 mm / 通過長 |
| `detected_truth_mev` | 同じ検出・track選別済みMC hitのtruth edep × 2 mm / 通過長 |
| `linear_prediction_adc` | 同じMC hitのtruth edep / 0.305 MeV × 生成に採用するMIP × 2 mm / 通過長 |

温度補正は現行式 `C=1−(T−20)×α`、αは10 μmで1.6/135、15 μmで3.5/230です。HG<2600を応答比較に使用し、HG=4000の強制LG切替は入れません。保存HG・採用energyが現在の定数で再現するかは、応答選別とは別に照合します。

`unmatched`は通過位置・角度をそろえる前、`matched`はrun/channel/長手位置/幅位置/sx/syの共通binで最小件数を目標に重みを付けた分布です。bin軸にADC・energyは使いません。既定binは5 mm、1.5 mm、0.01、0.01、両試料3 hit以上です。取得・HG選別済みhitについてのmatchingであり、取得前の母集団が一致したことを保証しません。

data/MCで同じADC bin、同じfit範囲のLandau×Gaussianを使います。MPVはLandauのMPVパラメータで、畳み込み後の最大位置とは区別します。下端はthreshold+3σを通過長・温度で変換した保守的な共通下限、上端はHG切替に届かない共通範囲です。

fit区間の有効統計100未満、fit失敗、共分散不良、χ²/NDF>3、パラメータ境界、MPVが下限付近の場合は未判定です。温度補正後ADCについて2σ・設定値・4σの下端変更も調べます。「3%以内で安定」は診断目安で、確定校正値の採用条件ではありません。

別の`TruthHits`には、再構成trackを要求しない生成muonのtruthを保存します。通過長はprimaryの直線とstrip形状から得る近似です。散乱やsecondaryの実際のGeant4 step通過長ではなく、primaryが幾何的に通過しないdepositは未判定として数えます。

## 出力先と見方

```text
analysis/result/channel_response/muon_mip_response/<label>/
  inputs/                         manifest、設定、定数監査、hash、現在のソース
  run.log                         実行記録
  flat.root                       Crossings / TruthHits / Events
  spectra.root                    全channel分布、flow、誤差、fit関数
  validation.tsv                  対応付け・欠落・重複・track収集件数
  comparison.tsv                  channel別MIP再現、data/MC比、fit判定
  fits.tsv                        各段階のMPV・幅・誤差・fit範囲・不可理由
  plateau_stability.tsv            threshold下限変更への依存
  dependency_fits.tsv             位置・温度・run bin別のMPVと判定
  geometry_bins.tsv               共通binと重み
  geometry_coverage.tsv           matchingで残った割合
  efficiency.tsv                  track条件付きのhit取得割合
  calibration_closure.tsv         保存校正値の再現誤差
  group_summary.tsv               SiPM種類・個別/平均補完別の集計
  common_group_summary.tsv        同じ有効channel集合のdata/MC比較
  fit_status_counts.tsv           fit成功・未判定理由別の件数
  event_summary.tsv               保存muon event energy・hit数の平均
  analysis_status.tsv             共通fitの有無と処理件数
  performance.tsv                 worker別cache使用・実際のfit呼び出し数
  comparison/<matched|unmatched>/100GeV/layer<layer>/chip<chip>/
    cell<CellID>_response.png      6枚の分布比較
    cell<CellID>_profiles.png      位置・温度・run別MPV
    cell<CellID>_mean_profiles.png 同じ変数に対する平均HG応答
  comparison/muon_mip_channel_map.png
                                  data/MC×matching前後のMPV map

analysis/log/channel_response/muon_mip_response/
  slurm-<jobID>.out / .err         Slurmログ
  submissions.tsv                 投入記録
```

全種類のPNGはcanvasを2250×1500 pixel（従来の1500×1000から縦横それぞれ1.5倍）で描画します。

channel PNGはlayer・chipごとのdirectoryに保存し、channelはファイル名のCellIDで識別します。CellIDは `layer × 100000 + chip × 10000 + channel` です。例えばlayer0・chip2・channel13は `comparison/matched/100GeV/layer0/chip2/cell20013_response.png` になります。`unmatched` と3種類のchannel PNGも同じ階層です。

channel図の6枚は左上から、pedestal差引きADC、通過長補正ADC、温度補正ADC、保存校正後MIP、生成truthと検出条件付きtruth、MCの線形予測と最終ADCです。黒=data、赤=MCが基本ですが、下段中央・右は凡例に示すMC内部の比較です。実線は分布（dataは黒、MCは赤）、破線は採用fitの曲線（dataは灰色、MCは淡赤）です。fitのMPV±誤差、Landau幅 `w_L`、Gaussianの `sigma_G`、χ²/ndfをlegendに表示します。`w_L`はLandauのscale parameterで、標準偏差ではありません。`ok`は表示せず、不採用の場合だけlegendの別行へ `insufficient_...` などの理由を表示します。分布はflowを含む総件数で規格化するため、表示範囲内の面積が1とは限りません。

全channelを処理しますが、既定のchannel PNGは正常な主fitがあるchannelに限ります。未判定もROOT・TSVには残ります。特定channelを`--channels`で指定すると、未判定の図も保存します。位置・温度・run別のMPV fitは主fitが有効なchannel、または明示されたchannelに実施し、bin統計不足は未判定として残します。MPV profilesの長手位置・幅方向位置・温度には、fitに使ったbinの両端まで横棒を表示します。点はbin内の重み付き平均位置なので、横棒は左右非対称です。既定幅は長手5 mm・幅方向1.5 mm・温度1℃で、strip pitchの5.3 mmとは別です。これは測定誤差ではなくfit bin範囲で、runの横棒は0です。縦棒はMPVのfit誤差です。平均profileはMPVではなく、今回のfit bin横棒の対象は `_profiles.png` です。

`analysis_status.tsv`が`no_common_valid_matched_fit`なら、処理が完了してもdata/MCのMIP一致を確認できていません。MC-onlyでは`unmatched`のfitとMC内部比較を見ます。MPVが1からずれても、基準条件・選別・fit安定性を確認してから原因を判断し、crosstalk補正やMIP更新を自動適用しません。

`analysis_status.tsv`は保存HG・採用energyの再現不一致も集計し、不一致がある場合は`saved_calibration_mismatch`を優先します。この場合は定数の指定・生成時の条件を確認してください。`comparison.tsv`のtruth・保存MIPの値は、それぞれのfitが有効な場合だけ出力し、無効ならNaNです。channel mapの0/空欄は有効fitなしを表し、応答が0という測定ではありません。
