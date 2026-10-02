# MuonでのMIP規格化・fallback校正の検証

Slurmログの保存先は `analysis/log/channel_response/mip_closure/` です。

2023 SPSの100 GeV muonを対象に、保存されたMC edep → digitization後HG ADC → 校正後MIPをchannel別に調べます。data・MCに同じ解析・校正方針を適用します。既存ROOT、定数、digitizer、校正ソースは変更しません。

## 実行

ROOTとPyROOTが使用できる環境の `python3.12`、NumPy、C++コンパイラが必要です。

```bash
cd ~/ScECAL_BeamTest/analysis/channel_response/mip_closure
make check

# 入力の存在、tree、復旧状態を点検。出力・ビルド・ジョブ投入なし。
bash run.sh --dry-run

# 警告のないRun44で処理を確認。channel別fitは統計不足になり得る。
bash run.sh --runs 44 --max-tracks 3000 \
  --output ../../result/channel_response/mip_closure/run44

# 復旧ROOTも読み取る限定的な診断。正式な全run比較とは区別する。
bash run.sh --runs 29,40,57 --allow-recovered --max-tracks 3000 \
  --output ../../result/channel_response/mip_closure/diagnostic

# 新しい出力先で全run・全保存trackを処理。復旧ROOTは既定で除外。
bash submit.sh --output ../../result/channel_response/mip_closure/full

# trackFit3Dの状態と独立に、全生成muonの方向・位置を監査。
python3.12 audit_beam.py \
  --output ../../result/channel_response/mip_closure/beam_audit
```

出力先が存在すると停止します。`--max-tracks N` は各track treeから等間隔にN件まで抽出し、0なら全件です。抽出前のイベントが保存されていない場合は復元できません。`--no-fit` は収集・定数・energy・beam監査まで実行します。全オプションは `bash run.sh --help` を参照してください。

## 入力と実際のbeam設定

- data: `ECAL_data/analysed/2023/sps/{decode,calib,trackFit3D}/mu-/100GeV/`
- MC: `Result_MC/{decode,calib,trackFit3D}/mu-/mu_track/threshold/100GeV/`
- MC truth: `Result_MC/generate/mu-/mu_track/100GeV/`
- MIP: `Analysis_edit/share/all_auto_muon_v4_trackfit.root`
- pedestal: `Analysis_edit/share/pedestal2023_SPS.root`
- plateau範囲用のthreshold: `analysis/result/threshold/threshold.root` の `expErfThre`。data・MC別に指定可能。旧 `threshold_rev.root` は異なる形式なので使用しません。

`--data-base`、`--mc-base`、`--{data,mc}-{mip,pedestal,threshold}`、`--digi-mip` で変更できます。生成時のMIPと校正時のMIPは別入力です。

2026-09-30に確認した生成マクロは `/gps/direction 0 0 1` でした。`PrimaryGeneratorAction.cc` はrun別の `h_hitmap_smooth` から(x,y)をサンプリングし、z=-60 mmに置きます。このmapを作る `beamsize_track.cc` はdataのtrackの**z=0切片**にpitch幅の一様乱数を加えています。個々の最初のhitを対応させる方法でも、位置と角度の同時分布を再現する方法でもありません。垂直入射ならsourceからz=0までx,yは変わりません。

`audit_beam.py` は保存truthの全primary方向・位置と、現在の生成マクロ・位置mapの平均と幅を出力します。現在のソースやmapが保存ROOTの生成時と完全に同一だったことまでは保証しません。

同日にMC `trackFit3D` を確認すると28run中27runでROOT復旧が必要でした。Run44には復旧警告がありませんでした。既定では問題のあるrunをdata・MCの対ごと除外して `inputs/excluded_runs.tsv` に残します。`--allow-recovered` は読めた保存entryによる診断を許可し、`STATUS.json` に明示します。全runの確定評価には完全なtrackファイルの用意が必要です。警告がないことだけで全生成イベントの処理完了を保証するものでもありません。

## 共通の選択・位置・通過長

1. 各run内の `(Event_Time, TriggerID/Event_Num)` でraw、校正後、trackを対応付けます。entry順では対応付けません。MC truthだけはdigitizerが `TriggerID=i_entry` とする規約を使い、ID範囲・時刻を検査します。
2. 最終 `trackFitPars` が有限、NDF>0、各向き7層以上、投影角 `|atan(sx)|, |atan(sy)| < 0.05 rad`、`sqrt(objective/NDF)<3 mm` を既定とします。後者は既存3D fitの無規格化距離目的関数であり、確率的なχ²ではありません。角度と距離は変更できます。
3. 45×5×2 mmのactive strip内を厚さ方向に通過する条件を両試料に適用します。長手側1 mm、幅側0.25 mmの余白を確保し、厚さによる横方向移動も含めます。2 mm厚に対する通過長は最終傾きから `2*sqrt(1+sx²+sy²)` と再計算します。保存 `realLength` は診断用に比較するだけです。
4. MIP用hitはtrackの採用channel、校正hit中の物理channel重複なし、raw HG<2600、raw hit tagありを要求します。rawとのhit対応はmemory cellを含む完全CellIDで行います。SSAで複製された信号を重ねて数えません。
5. MIP fitはrun・channel・strip長手位置・幅方向位置・sx・syの共通binで行います。既定のbin幅は5 mm・1.5 mm・0.01・0.01です。両試料3hit以上のbinで、共通目標数 `min(Ndata,Nmc)` を各試料の件数で割った重みを使います。ADC・energyを重み作成に使いません。非共通領域は両方から除外します。

同一cutだけでは角度分布は一致しないため、bin単位の対応と重みを保存します。bin内部には分布差が残り得るので、`matching_diagnostics.tsv` とbin幅変更で確認してください。MCに存在しない真の角度分布を重みで作ることはできません。既存SSA・track前段の選択による偏りも残ります。

## MPVの定義とfit

data・MCともLandau×Gaussianの**Landau MPVパラメータ**を同じ関数・bin・範囲決定法でfitします。畳み込み後の最大位置は別列 `convolution_mode` に保存します。raw pedestal差引きADC、2 mm通過長補正後ADC、現行温度補正式を併用したADC、保存HG MIP、MC truth edepを対応付けます。

実際のfit対象 `adc` は `(HG-pedestal)*2/length`、`adc20` はこれに現行校正の `1-(T-20)*coefficient` を掛けたものです。係数は10 µmで1.6/135、15 µmで3.5/230です。補正式を独自に変更せず、保存HG枝の再現性を `saved_calibration_check.tsv` で検査します。入力MIP定数が20℃基準として妥当かは別途確認が必要です。

hardware threshold近傍の未知の効率をfitで吸収しないよう、両試料共通の下端を `threshold+3*sigma` と校正のHG>10 ADC条件から決めます。通過長・温度補正後の変動を含め、全採用hitでその下端以上となる範囲を使います。threshold欠損は同chipの有限な正幅の値から平均し、平均元もなければfitを不可とします。上端は10 µmで800 ADC、15 µmで2000 ADCを基本とし、HG切替を越えません。

fit範囲内の有効統計100件未満、失敗status、共分散不良、MPV≤0、NDF≤0、非有限値、パラメータ境界、χ²/NDF>3、MPVがfit範囲外なら正常fitに数えません。thresholdがMPVに迫るchannelは**未判定**になります。`--plateau-sigmas` とfit条件依存を確認し、足りない統計を低い品質条件で補った結果を確定値として扱わないでください。

校正後MIP MPVは `adc20 MPV / 採用MIP定数` で求めます。定数で割るだけの線形変換なので再fitによる違いを導入しません。保存HG枝は独立にfitします。data/MC比の集計には両方で正常fitになった同じchannelを使います。MC truth MPVは同じ**検出・選択済みhit**の条件付き分布です。

`mc_adc_mpv_over_input` が主診断です。さらに `mc_truth_mpv_over_0p305` と、その比でADC比を割った値を保存し、edep規格化の影響を区別します。入力MPVの不確かさ・温度・位置・fitモデルの系統誤差はfit統計誤差には含めません。

`primary_pe_scale_trial_only = 1/(MC ADC MPV / 入力MIP)` は**再digitizationの初期試行値**です。crosstalkを有効にしたまま一次光電子数へこの係数を試し、再生成・同じ解析でMPVが入力に戻るか反復確認するための値です。MPVと一次光電子数の関係は厳密な線形ではないため、この列を確定補正として自動適用しません。本コードはdigitizationを変更・再実行しません。

## fallbackと同一hitでのenergy

- `legacy`: 現行校正のχ²/NDF cutのみと `MIP==1` のfallback sentinelを再現します。負MPVや0/0の比較動作も監査対象として保持します。
- `digi`: 現行Extractの負MPV除外、χ²cutなし、fallback sentinelに加え、digitizer内のMIP≤5 → 90/400 ADCの上書きも再現します。
- `clean`: MPV>0、NDF>0、MPV・幅・Gaussian幅・χ²の有限性、χ²≥0と従来のχ²/NDF cutを要求します。接続channelの正常値だけを10/15 µm別に平均し、異常・欠損channelに代入します。data・MCの両方で同じルールです。

総energy比較では選択muonイベントの**保存Calib_Hitにあるlayer0–29の全接続hit**を固定します。MIP fit用の位置・HG選択や重みを総energyに持ち込みません。元のhit選択・HG/LG切替を保持します。HGはrawから旧/新MIPで再計算し、LGは保存HGによって旧校正の再現を確認した上で保存energyに `旧MIP/新MIP` を掛けます。HG/LG係数・threshold・crosstalkは変更しません。

旧MPVが負でも有限なら負energyをbaselineに保持します。負energyの正常化とfallback平均変更は逆方向に働き得るため、`channel_energy.tsv` の `n_negative_legacy` と差分を確認してください。旧MPVが0/非有限、rawが欠損/重複、保存HGの再現に失敗したhitがあれば、そのイベント全体の比較を無効にします。比較可能hitだけの小計は別列として残し、総energyと混同しません。追加のbad-channel maskは導入せず、保存校正前段の選択はそのままです。

## 主な出力

| ファイル | 内容 |
|---|---|
| `inputs/`、`STATUS.json`、`validation.tsv` | 入力・除外run・復旧状態・source/定数hash・実行設定・join件数 |
| `*_constants.tsv`、`*_fallback_means.tsv`、`*_changed_channels.tsv` | 採用値・採用理由・fallback平均・変更channel |
| `flat.root` | 選択イベントの全校正hit、raw、track、truthの対応と全primary |
| `matched_hits.root`、`matching_*.tsv` | MIP用hitと共通重み、共通領域のcoverage、残る分布差 |
| `spectra.root`、`fits.tsv` | channel別分布・fit曲線・誤差・status・不可理由 |
| `channel_closure.tsv` | ADC MPV/入力MIP、旧/新校正MIP、data/MC、truth規格化 |
| `group_summary.tsv`、`calibrated_group_summary.tsv` | 10/15 µm・個別/fallback別の正常channel集計 |
| `event_energy.tsv`、`energy_summary.tsv`、`channel_energy.tsv` | 同一hitでのenergy変化と負energyの寄与 |
| `saved_calibration_check.tsv` | 現行候補定数で保存HGを再現できるか |
| `beam_truth.tsv`、`beam_reconstructed.tsv` | 入射方向と再構成傾き・通過長の診断 |
| `closure_summary.png/pdf`、`report.txt` | 比較図と結果の読み方 |
| `channel_axis.tsv` | 図の連番と元のCellID・layer・chip・channelの対応 |

図はPNG・PDF・SVGの3形式で保存します。PDF・SVGは拡大しても文字や線が潰れません。

- `closure_summary.*`: 10/15 µm別のMC ADC MPV / 入力MIPの分布。個別校正/fallbackを色と線種で区別し、channel数と中央値を表示します。1 channelを1票とし、全有効値を横軸範囲に含めます。
- `mc_channel_map.*`: 横軸をlayer内の連番 `36*chip + channel`（0〜209）、縦軸をlayerとした応答マップ。灰色は有効なMC ADC fitなし、青/白/赤は1未満/1付近/1超です。色は0.5〜1.5で飽和し、範囲外channel数を図に明記します。個別/fallbackを合わせて表示します。
- `data_mc_pairs.*`: 両方でfitが成立したchannelを1行ずつ表示するdata/MC比。layer/chip/channelを明記し、dataとMCを独立としたfit統計誤差の伝播を表示します。系統誤差は含みません。
- `visualization.json`: 描画ソースと入力TSVのhash、描画channel数。再描画はfitや校正を変更しません。

chip 0〜4は各36 channel、chip 5は30 channelです。`channel_axis.tsv`には全体連番
`210*layer + 36*chip + channel`（30 layerで0〜6299）とCellIDの対応も残します。

保存済み結果から図だけを再生成できます（再fit不要）。

```bash
python3.12 visual_report.py ../../result/channel_response/mip_closure/run40_channel_index_20260930
```

`STATUS.json` の完了は処理完了を示します。`valid_paired_channels=0` や小さい共通coverageは、MIP一致を確認できたことを意味しません。

### 全muon runの統合

`--runs all --max-tracks 0`で、見つかった全runの全保存trackを処理します。
run/channel/通過位置/傾きごとにdataとMCの共通領域を求めた後、重み付きhitを全runから集め、
各channelに対して一つの分布をfitします。run別MPVの平均ではありません。
energyの旧/新校正比較は、全runの選択イベントの同じhit集合で行います。

```bash
bash run.sh --runs all --max-tracks 0 --allow-recovered \
  --output ../../result/channel_response/mip_closure/all_muon_runs_20260930
```

`--allow-recovered`では復旧されたROOT入力も診断として含めます。
採用runとentry数は`inputs/files.tsv`、除外は`inputs/excluded_runs.tsv`に記録されます。
全run統合はメモリ使用量が大きいため、十分なメモリを確保した計算ノードで実行してください。

data/MC比較は20 channelごとに分割し、`data_mc_pairs.pdf`に全ページを保存します。
複数ページの場合、PNG/SVGは`data_mc_pairs_page_001.*`等に分け、
`data_mc_pairs.png/svg`には1ページ目を保存します。
