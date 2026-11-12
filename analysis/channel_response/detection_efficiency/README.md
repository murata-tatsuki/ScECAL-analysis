# Muonによるdata/MC hit取得割合

他layerのmuon trackから予測した物理stripに、decodeの`HitTag=1`が記録される割合を測るC++17/ROOT解析です。electron occupancyから効率を逆算せず、MCの本番間引きや校正・digitization・入力ROOTの変更は行いません。data/MCは同じ収集・統計処理を使います。

## 実行

```bash
cd /home/murata_t/ScECAL_BeamTest/analysis/channel_response/detection_efficiency
bash run.sh --build-only
bash run.sh --label pilot --runs 40,44 --energy 100 --max-events 500
bash run.sh --label pilot --plots-only
bash submit.sh --label production --runs 40,44 --dry-run
# 入力を明示する場合（muon_mip_responseと同じ6列形式）
bash run.sh --label explicit --input-manifest /absolute/manifest.tsv --max-events 500
# 収集と集計を別実行。再利用元の定義・設定・checksumを確認する。
bash run.sh --label collected --runs 44 --max-events 500 --collect-only
bash run.sh --label analyzed --reuse-crossings ../../result/channel_response/detection_efficiency/collected
```

既定は100 GeV・全muon run・全event。**今回の小規模検証では`--runs`と`--max-events`を明示しています。** `--max-events`は各入力fileから等間隔抽出。5/20 GeVなども`--energy`で指定できますが、実在するmuon入力がなければ停止します。2026-10-13確認時はdataは100 GeVのみ、MCは100/10 GeVでした。electron energyの存在からmuon入力を推測しません。

`--help`に全optionを表示します。`--dry-run`はパス・指定を確認して終了し、生成・ビルド・投入しません。`submit.sh --dry-run`は実際の`sbatch`引数も表示します。現在のSlurm `all` partition、既存`../execute.sl`を読み取りで再利用し、既定2 CPU/8 GB/4時間。`--partition --mem --time --workers`で変更可能です。収集・統計は各1 process、描画は既定2 processで分割。割当CPUを超えるworkerを拒否し、同一ROOTへの並列書込みはありません。

`validation_100gev`と`validation_reuse`は実装時に作成済みの検証結果です。完了報告の実行記録を同じlabelで再実行すると、上書き防止のため停止します。既存の数値を確認するだけなら再実行は不要です。図だけ更新する場合は`bash run.sh --label validation_reuse --plots-only`を使います。

新しく同じ少量条件をSlurmで実行する場合は、未使用labelで次のように確認・投入します。`--dry-run`だけではジョブは投入されません。

```bash
bash submit.sh --label pilot_slurm --runs 40,44 --energy 100 --max-events 500 --workers 2 --dry-run
bash submit.sh --label pilot_slurm --runs 40,44 --energy 100 --max-events 500 --workers 2
```

新規解析は既存labelを拒否します。`--plots-only`は既存`statistics.root`と保存設定だけを読み、`redraw_XXXXXXXX/`以下にPNGを作り、元の数値・図を上書きしません。元decode/calibration/truthやcrossings.rootは読みません。描画前後の数値checksumを保存・比較します。

## 入力と由来

既定data: `/home/murata_t/data_beamtest/ECAL_data/analysed/2023/sps/{decode,calib}/mu-/100GeV/`。
既定MC: `/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/`下の

- `{decode,calib}/mu-/mu_track/threshold/100GeV/ECAL_Run<run>.root`
- `generate/mu-/mu_track/100GeV/ECAL_Run<run>.root`

`--data-base --mc-base --mc-tag --input-manifest`で変更できます。electron `e-/e_ssa`入力は使用しません。manifestはヘッダー任意、絶対パスの6列TSV:

```text
sample(0=data/1=MC)  run  decode.root  calib.rootまたは-  -  truth.root(dataは-)
```

sample/runは一意、canonical decode pathの重複を拒否します。複数fileに分割した同じrunを黙って合算しません。MC truthの全entryについて粒子種`abs(PDG)==13`と要求energy（MeVへの換算）を確認します。

- rawの`Run_Num,Event_Time,CycleID,TriggerID`をscalar passで索引化し、重複identityは全出現を除外します。data runはmanifestとの一致を要求します。MC内のRun_Num=0はsource_runとして保持します。
- calibにCycleIDがないため、両側で一意な`Run_Num,Event_Time,TriggerID/Event_Num`のみを照合します。並べ替え・欠落・曖昧なIDを監査し、entry順によるjoinはしません。**calibrationとの照合可否で取得効率の分母を減らしません。hit残存率は計算しません。** 校正hit/energyを読みません。
- truthには保存event IDがなく、効率へのevent joinはしません。現在のMCDigi sourceは参考として保存しますが、過去の生成設定やCT補正適用の証明にはしません。ファイルの同じproduction・CTの歴史的設定は生成receiptが未確認なら`unverified/unknown`です。入力ファイル一覧、ROOT schema、stat、現行producer source、解析source/hashを保存します。
- 既存muon `flat.root/Crossings`とMIP `crossings.root`は分母の考え方の参考です。隣接strip、群除外の再track、全identity、今回の割当て不確かさの情報を満たさないため、直接再利用しません。今回の`--reuse-crossings`は自解析のschema・収集完了・ROOT/config/manifest/source checksumが一致する場合だけ可能です。集計条件も引き継ぎ、異なる指定は拒否します。

## Track、分母、hit定義

共有`../muon_mip_response/TrackGeometry.hh`と実座標`../../calibration_parameters/EBUdecode.cxx`を変更せず再利用します。全30layer、chip0–4は36strip、chip5は30strip、6300接続channel。chip5のchannel30–35の180未接続channelも一覧に状態を出します。

HitTag=1のstrip中心を非加重でtrackに使用。奇数層でX、偶数層でY。**対象layer全体をtrack fitの前に除外**します。対象stripのADC・hit有無・校正energyは選択に使いません。maskしたstripはtrackからも除きます。保存済みtrackFitは不要です。

既定は各投影7layer以上、各投影angle≤0.05 rad、各投影RMS≤3 mm。45×5×2 mmのactive stripを全厚通過し、長手1 mm/幅0.25 mmのmarginを要求します。幾何gapや境界はhardware欠損と扱いません。位置・通過長はtrackから予測した値で、独立な外部trackerによる入射位置ではありません。

track位置の不確かさを、他層のOLS leverageと`max(5.3²/12, RMS²*n/(n−2))`の位置分散から近似します。既定1σの予測区間がactive境界を越す割当てを曖昧とし、明確な効率分母から除きます（`--ambiguity-sigma`）。これはstrip量子化を含む診断近似で、alignment・散乱systematicを測定したものではありません。2層以上で第1/第2候補の残差差が0.25 mm未満ならtrack候補の曖昧さも記録・除外します。検出できない複数trackの曖昧さは残り得ます。

`geometry_predicted`はmargin込み全厚通過、`predicted`はそこから曖昧な割当てとmaskを除いた分母。`Crossings`にはtrack不成立・領域外も含め全対象attemptを保存します。`cutflow.tsv`はsample/run/variant/group別の段階的な分母です。

- `observed_any`: 対象stripに1件以上HitTag=1。
- `observed_unique`: 対象stripに記録がちょうど1件。
- `duplicate_strip_crossings`: 2件以上の記録がある通過。
- `extra_memory_hits`: 同じ物理stripの異なるfull CellIDの2件目以降。`repeated_full_ID_records`は同一full ID重複の別計数。
- `neighbor_only`: 対象にhitがなく、gridで辺を共有する上下左右の接続・unmasked stripにhitがある通過。chip番号順の隣接ではありません。
- `ambiguous`・maskは別計数で、原因カテゴリは互いに排他的とは限りません。段階比較にはcutflowを使用します。

`Events`の`raw_records`は全decode vector長、`hit_records`はHitTag=1記録数、`physical_strips`は接続channelにまとめた物理strip数。SSA仮想hit数は入力がなく未測定（-1）です。memory重複はtrackでは同じ座標にまとめ、取得分子ではany/uniqueを分けます。未接続hit・mask hitも別に記録します。

muonにelectron低LG tail/energy選別を適用せず、sample別設定と除外0件を保存します。CoGなし、温度・通過長による信号補正なし、負/ゼロenergyによる除外なし。ADC/energy branch自体を読まないため、未検出hitのADCを測定したことにはなりません。「検出hitのADC」を無条件効率の横軸として扱う図は作りません。

## 比較と統計

`unadjusted`は全予測通過。`common`はvariant/run/channel/長手位置/幅位置/sx/sy/RMSの共通binにdata/MCそれぞれ既定3通過以上ある領域での比較です。bin幅は10 mm/1 mm/0.02/0.02/1 mm、floor bin、原点0。ADCや取得結果をmatchingに使いません。

各共通binでmin(data分母, MC分母)件を、event identityの固定hash順で各試料から選びます。**解析上の通過条件を比較するための選択であり、MCの補正・本番間引きではありません。** 同じchannelとtrack条件binの母数をそろえますが、bin内の残存差・track biasは残ります。共通化前後の絶対分母・分子とcoverageを`geometry_bins.tsv`と集計表に残します。muon MCのrunは生成元runラベルで、dataと同時取得したrunとは主張しません。

channel/chip/layerは1event当たり最大1予測なので、any/uniqueのWilson区間（既定95%）を保存。data−MCは独立2標本のNewcombe区間。run/全体/依存性では1event内の複数layerを独立試行と見なさず、event clusterの比率誤差と、有効event数`N²/Σ n_event²`のWilson区間の広い側を使用します。これは保守的な近似区間で厳密なcoverage保証ではなく、全hit取得/欠損時にも不確かさをゼロとしません。推定値・分母は低統計でも保存しますが、既定20通過/20event未満をinsufficient、分母0はNaN/unmeasuredとします。mask・未接続は別状態です。fitは不要で、fit失敗値のゼロ置換もありません。

stabilityの既定4条件: nominal、margin 0/0 mm、margin 2/0.75 mm、各投影9層/angle0.035/RMS2 mm。`--stability 0`でnominalのみ。元のnominal指定より緩いtight条件にはしません。各条件で分母を再構築します。

## 同時欠損とchipの限界

mode0は各対象layerだけ除外する診断です。他の対象layerのhitがtrackに使われるため、複数欠損には選択biasがあります。mode1–3は既定`0–3`, `10–13`, `26–29`の**群全体をtrack前に除外**し、同じtrackから群内を予測します。群は`--groups`で変更可能です。groupごとのcutflowを記録し、残るlayerのtrack条件に対する条件付き測定である点は共通です。

event当たり予測数と欠損数を絶対件数で保存。独立欠損期待はsample/run/variant/mode/layer別の欠損率を、当該eventを除いて推定し、そのeventの予測layer集合にPoisson-binomialを適用します。各layerの残り分母がmin-count未満なら期待値は未測定。その条件を満たす同じevent集合の観測と期待を別保存します。期待のSDは固定推定率下の条件付き分散で、欠損率推定誤差・geometry/track systematicを含みません。significanceの確定値ではありません。位置・channel組成の違いでも独立layer期待からずれ得ます。layer pairの両欠損と各周辺率の積もTSVに保存します。

今回の入力schemaには独立なchip readout状態の解釈や、同じchipの複数独立通過を作る情報がありません。HitTag/BCID等はhit由来です。単一trackはlayer当たり1stripしか予測せず、他channelにhitがないことからchip全欠損とは判定しません。測定結果はchip別の条件付き取得割合までです。

## 出力と図

`analysis/result/channel_response/detection_efficiency/<label>/`。CoG階層を形式だけ追加しません。Slurmログは`analysis/log/channel_response/detection_efficiency/slurm-%j.{out,err}`。

| ファイル/図 | 内容 |
|---|---|
| crossings.root | 全対象attempt・identity・track・geometry・取得/memory/隣接状態、Events |
| statistics.root | 集計Summary、Dependencies、予測数/欠損数/独立期待histogram |
| efficiency.tsv / comparison.tsv | channel/chip/layer/run/全体の分母・分子・any/unique効率・CI・差・状態 |
| cutflow.tsv / input_audit.tsv | 選別段階、保存event ID照合・truth muon/energy確認 |
| unmeasured_channels.tsv | 未接続・mask・分母0・統計不足を含む一覧 |
| geometry_bins.tsv / dependencies.tsv / dependency_comparison.tsv | 共通binの元/採用母数、予測位置・角度・track品質の分母付き依存性 |
| simultaneous.tsv / simultaneous_comparison.tsv / event_distributions.tsv / layer_pairs.tsv | 同時欠損・event分布・独立期待、対応母数と不確かさ |
| chip_readout_audit.tsv | chip一括欠損を測定できない理由 |
| variants.tsv / groups.tsv / sample_selection.tsv | 実際に適用した設定・sample別選別 |
| inputs/ / run.log | manifest、設定、source archive/hash、入力stat、由来の限界、コマンド |
| comparison/{unadjusted,common}_nominal/<energy>GeV/layer_efficiency.png | layer効率/差/絶対予測数の3panel |
| 同chip_efficiency.png / run_efficiency.png | chip集計/生成元run別。同じ3panel |
| 同layerN/chipN/channel_efficiency.png | channel別の効率・差・絶対分母。通過のあるchipのみPNG |
| 同layerN/{efficiency_map,denominator_map}.png | 実45×5 mm strip位置。data/MC/差。全layer同じx/y縮尺、data/MC共通色軸、未測定効率は灰色 |
| 同longitudinal_mm等のPNG | 予測通過位置・角度・RMS・使用層数に対する効率・差・分母 |
| 同stability.png | 4条件の効率と分母変化。相関する選別間の差を独立有意差とは解釈しない |
| {data,mc}/unadjusted_nominal/<energy>GeV/missing_modeN.png | 予測layer数、欠損数の観測と独立期待（絶対event数） |

PNGのみ。ROOT/TSV・再現情報を保持します。高速化は必要branchのみ、scalar identity監査＋vector収集1回、trackの不要な再計算削減、統計は保存Crossingsの2pass、描画は保存集計から並列実行。入力をchannelごとに再走査しません。

## 今回の検証（2026-10-13）

実入力は100 GeV muonのrun40/44、各file500event、data/MC計2000event。identity重複検出とtruthの粒子種/energyのscalar監査は各入力の全identityを確認しますが、track解析は指定した少量だけです。収集は`validation_100gev/`、保存Crossingsを再利用した最終数値は`validation_reuse/`にあります。取得効率TSVは再利用前後で一致しました。

確認済みPNGは`validation_reuse/redraw_GJaVF2dc/`以下です。初期描画にはmapのpalette/縮尺表示に不足があり、数値を変更せず修正しました。既存出力を上書きせず、修正済み描画を別の再描画directoryに保持しています。今後の`run.sh`/`--plots-only`では修正済み描画が既定です。

- ビルド、実入力収集・統計・PNG、Slurm dry-run、plots-onlyを確認。plots-onlyでstatistics ROOT・TSV・設定の前後checksumは一致。
- 336,000対象attemptの適格性・memory内訳・集計一致・共通channel母数一致・分母0のNaN・未接続状態を独立の一時C++検証で照合。
- 人工decodeで、未検出通過を含む分母、群内hit変更から独立なtrack、隣接のみ、memory重複、同一full ID重複、mask、重複event排除を確認。calibrationの順序反転・event欠落でもID照合でき、分母は不変。人工入力にはADC/energy branch自体を含めず、対象信号を選択に読んでいないことも確認。
- 異なるmin-layersを指定した再利用、plots-onlyでのmax-events変更は意図どおり拒否。
- 既存geometryの未使用温度補間関数にコンパイラ警告があります。本解析は温度関数を呼ばず、共通ソースを変更していません。

標準条件の無調整ではdata 6,949/9,364 = 0.7421、MC 22,144/25,094 = 0.8824。ただしtrack条件を満たすevent数・通過分布が異なります。**標準の共通binに残ったのは各3通過だけで統計不足です。hardware効率差やelectron hit過剰の原因を確定した結果ではありません。** 多数のchannelも低統計/未測定です。

`verification.tsv`、`synthetic_validation.txt`、`validation.tsv`、`inputs/slurm_dry_run.txt`、再描画先の`numerical_{before,after}.sha256`が検証記録です。一時検証コード・人工入力は完了後に削除します。全量/全energy解析、実Slurm投入、歴史的MC production/CT設定の確証、chip一括欠損の実測は未実施/未確認です。
