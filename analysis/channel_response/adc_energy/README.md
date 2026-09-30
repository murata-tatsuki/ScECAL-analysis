# ChannelごとのHG–LG・ADC–energy比較

decodeの`Raw_Hit`とcalibrationの`Calib_Hit`を対応させ、data / simulationを
物理channelごとに比較する。保存されているenergyを読み、再校正やfitは行わない。
hit残存率の解析は含まない。

## ジョブ投入

ROOT、g++、Bash、GNU coreutils / awk、flock、rsyncが必要。Python / PyROOTは不要。
イベント処理・描画はC++、入力の組合せと実行制御はshellで行う。

```bash
cd /home/murata_t/ScECAL_BeamTest/analysis/channel_response/adc_energy

# 入力の存在と投入コマンドを確認。ジョブや出力は作らない。
bash submit.sh --energy 100 --dry-run

# 100 GeVを集計し、完了後にlayer別の描画タスクを実行
bash submit.sh --energy 100 --label fast_v1

# PS + SPS全19 energyの集計→全570 layerタスク（最大32同時）
bash submit.sh --energy all --plot-concurrency 32 --label fast_v1

# 既に完成した集計ROOTから描画だけを投入（イベント再読込なし）
bash submit.sh --energy 100 --plots-only

# PSの0.5, 1, 2, 3, 4, 5 GeVだけを投入
bash submit.sh --beam ps

# SPSの10–250 GeVだけを投入
bash submit.sh --beam sps

# energyを指定してまとめて投入
bash submit.sh --energies 0.5,1,5,30,50,100,150,250

# 状態確認
bjobs
# または
squeue -u "$USER"
```

標準のenergyは以下の19点。

- PS: 0.5, 1, 2, 3, 4, 5 GeV。
- SPS: 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200, 250 GeV。

`--beam auto`（初期値）ではenergyからdataのPS / SPSを選ぶ。
`--beam ps` / `--beam sps`で片方だけをまとめて投入できる。
`--energy`または`--energies`を指定すると、そのenergyの一覧を優先する。
各energyの集計ジョブでdataとMCのヒストグラムROOTを先に確定する。
指定した全energyの集計ジョブが正常終了すると、`afterok`依存の描画ジョブ配列が開始する。
1タスクは1 energy・1 layerで、同じ投入全体の同時描画数を`--plot-concurrency`（既定32）に制限する。
19 energy × 30 layerなら570タスク。channel指定時はそのchannelが属するlayerだけを作る。
実際の並列数は空き資源に従う。共有資源への負荷を下げる場合は4や8などに設定する。
集計に失敗した場合は描画を開始しない。完成した一部energyを描画するには
そのenergyを指定して`--plots-only`を使用する。
channel / event数の制限は初期設定では行わない（対象layerは0–29）。

この環境の`bsub` / `bjobs`はSlurmの互換コマンド。
`submit.sh`は`sbatch`と共通の`../execute.sl`を使い、各集計・描画タスクに1 CPUを割り当てる。
全energyについて`run.sh --dry-run`で入力パス・ファイル対応・引数を確認し、
成功した場合に投入前の`make`を1回実行する。ジョブには投入元の環境変数と作業ディレクトリを引き継ぐ。
ROOTやg++の環境を設定済みのログイン端末から実行する。

queueに相当するpartitionは`all`。`--partition`で変更できる。
実行時間とメモリはclusterの設定を使い、必要に応じて`--time` / `--mem`を指定できる。
`--label`、`--mc-tag`、`--no-plots`など、`run.sh`の解析オプションもそのまま渡せる。

ログは集計が`analysis/channel_response/job/adc_energy/<energy>GeV-fill-<jobID>.out/.err`、
描画が`plot-<arrayJobID>_<taskID>.out/.err`。
工程・energy・ジョブIDは`submissions_<日時>_<PID>.tsv`、
配列taskIDとenergy/layerの対応は`plot_tasks_<日時>_<PID>.tsv`（0始まりの行番号）に記録する。
`--log-dir`でログ保存先を変更できる。解析結果は`analysis/result/channel_response/adc_energy/`に保存する。
同じ出力名で再実行すると、その解析のROOT・PNGを更新する。
異なるrun・選択・試験条件は`--label Run105`などで出力を分ける。

## 端末で直接実行する場合

```bash
# 100 GeVの全イベント・全有効配線channel、layer 0–29
bash run.sh --energy 100

# 0.5 GeVではPS dataを自動選択
bash run.sh --energy 0.5

# 図の確認用に小さく試す場合のみ、channelとevent数を制限する
bash run.sh --energy 100 --channels 420008,920012,1420010 \
  --max-events 1000 --label smoke
```

`run.sh`は必要なら`make`を実行する。`submit.sh`と異なり、端末上で直接処理する。
直接実行では集計後に描画を順次実行するため、layer並列化には`submit.sh`を使用する。

## ローカル描画と出力負荷

`plot.sh`はdata/MCの集計ROOTをノードの`/tmp/channel_response_cache_<UID>/`にコピーする。
パス・サイズ・更新日時・inodeからキャッシュ名を決め、flockで同時コピーを防止する。
同じノードの後続layerタスクはコピーを再利用する。入力更新時は別キャッシュになる。
コピー中に入力が変われば停止する。既存の集計ROOTを書き換える処理と同時には実行しない。

PNGは`/tmp/channel_response_plot.XXXXXX/`で完成させ、layer単位でrsync転送する。
作業用ディレクトリは終了時に削除する。入力キャッシュは再利用のためノード上に残る。
`--scratch-root PATH`で別の高速なローカル領域を指定できる（NFSは避ける）。

出力先の`.plot_<hash>.lock`は同じlayer指定の描画を重複実行しないために使う。
ロックを取得したジョブが正常終了・エラー終了・HUP/INT/TERMで終了した際に削除する。
ロック取得に失敗したジョブは、実行中の別ジョブのロックを削除しない。
SIGKILLやノード停止時は終了処理が動かないため、ファイルが残ることがある。

キャッシュ容量はそのノードで使用した入力ROOTの合計サイズと世代数に依存する。
同じUIDの描画ジョブが使っていないことを確認したうえで古いキャッシュを削除できる。

描画ログに入力準備・描画・共有領域への転送の所要秒数を出す。
さらに`PNG timing`で、画像の描画とPNGエンコードを分けて計測する。

PNG生成には`FastPng.hh`を使う。ROOTの丸いmarkerの描画は、点ごとに画像全体を
処理するため、大きな複数panelの図ではprofile描画が重くなる。marker周囲の小領域を
背景ごと取り出し、同じROOTの円描画処理を適用して戻すことで、この負荷を減らす。
点・bin・誤差を間引かず、PNGのpixel寸法・panel配置・軸・色も維持する。
通常のROOT描画とのpixel一致を検証している。集計ROOTと任意のcomparison ROOTの
内容は変更しない。channel名の正規表現もファイル内の各keyごとに作り直さない。

PNGは可逆圧縮の速度を優先し、ROOTのcompression値10（PNG level 1）を使用する。
解像度・画質は変わらないが、圧縮後のファイル容量は従来より増える。
PNG枚数は同じで、全energyの総所要時間はノード負荷やI/Oにも依存する。
実装時の検証用ROOT・PNG・ログ・比較コードは検証後に削除する。

100 GeVの既存集計ROOTから12 channel（24 PNG）をローカルで描画した比較では、
旧版18.04秒→新版3.33秒、再計測で17.09秒→4.02秒だった（約4〜5倍）。
起動時間を含み、入力の共有領域からのコピーと最終転送は含まない。
24枚すべての全pixelが旧版と一致し、PNG容量は合計約2.5 MB→4.0 MBになった。
これは少量の描画ベンチマークであり、全energy解析の時間を測ったものではない。

既存の集計結果から新しい描画処理だけを実行するには、同じ`--label`等を指定して
`bash submit.sh --energy 100 --plots-only`を使う。元eventの再読込・再集計は不要。
同じ保存先に旧描画ジョブが実行中の場合は、その完了後に再描画する。
実行ファイルは原子的に差し替えるので、起動済みのプロセスは旧版を継続し、
差し替え後に起動するプロセスから新版が使われる。

従来はPNGと同じ図をcomparison ROOTにも保存していたが、現在はPNGのみを既定とする。
data/MCの**集計ROOT（全TH2D、profile、誤差、flow bin）は常に保存する**。
キャンバスも必要な場合だけ`--save-canvases`を付ける。

## 初期設定と入力

既存の比較で使っている2023 PS / SPSのdataと、MCの`threshold`条件を初期値とする。
`--beam auto`では0.5, 1, 2, 3, 4, 5 GeVのdataは`ps`、その他は`sps`を使う。

- data decode: `ECAL_data/analysed/2023/<ps|sps>/decode/e-/<energy>GeV/`
- data calib: `ECAL_data/analysed/2023/<ps|sps>/calib/e-/<energy>GeV/`
- MC decode: `CEPCScECAL_SML_Portable_update_new/Result_MC/decode/e-/sps/threshold/<energy>GeV/`
- MC calib: `CEPCScECAL_SML_Portable_update_new/Result_MC/calib/e-/sps/threshold/<energy>GeV/`
- pedestal: 現在の校正設定と同じ`Analysis_edit/share/pedestal2023_SPS.root`

完全なパスは`bash run.sh --help`と`--dry-run`で確認できる。
**低energyのMCも保存先は`sps`階層**なので、dataがPSでもMCのパスは上記を使う。
pedestalの初期値は現在の`ECAL_Analysis_LCIO/run/global_config`と
`run_simulation/global_config`の指定に合わせている。
別の校正条件を比較するときは、入力とpedestalを明示する。

`--data-decode`、`--data-calib`、`--mc-decode`、`--mc-calib`で
ディレクトリ、または1個のROOTファイルを指定できる。
ディレクトリでは同じbasenameのdecode / calibを組にする。対応ファイルの欠落はエラー。
dataとMCのファイル名・イベントIDが同じである必要はない。

pedestalは`--data-pedestal`と`--mc-pedestal`で別々に指定できる。
`ChnLevel/CellID, PedHighMean, PedLowMean`を読み、double / floatのvectorに対応する。
対象channelのpedestalがなければエラーにし、平均値などで補わない。
入力ファイルの生成時と同じpedestalであることを確認して指定する。

## 対応付けと選択

- 各decode / calibファイルの組の中で、`(Run_Num, Event_Time, TriggerID/Event_Num)`を照合。
- entry順序が検証できたファイルでは直接entry対応を使い、イベントキーのmap検索を省く。
- 続いてmemory cellの桁を含む**完全なCellID**でhitを対応させる。
- hit対応・重複検査には再利用する配列とイベント世代番号を使う。
  channel選択も配列化し、イベントごとのmap/set生成を省く。照合・重複検査は省略しない。
  ROOTの読み込みキャッシュは有効なbranchだけに設定する。
- 集計時のみmemory cellの桁を除き、`layer*100000 + chip*10000 + channel`でまとめる。
- 一意なevent keyならROOT内のイベント順が違ってもよい。hit順も問わない。
- 実dataには同じevent keyの繰り返しがあり、calib側にはCycleIDが保存されていない。
  この場合は、現在のcalibrationがdecodeのentry数・順序を保存する性質を使う。
  **両ファイルのentry数と全entryのevent key列が一致することを検査してから**同じentryを対応付ける。
  重複key内だけを並べ替えた外部生成ファイルは対象外とする。
- 対応したcalibration hitがdecodeのHitTag=1に存在することも確認する。
  重複keyがありentry列が一致しない場合、full CellIDの重複、対応イベントの欠落はエラー。
- 信号はdata / MCとも`HitTag == 1`。HG–LGはdecodeの信号を使用する。
- ADC–energyは対応したhitの保存`Hit_HG_Energy`・`Hit_LG_Energy`・`Hit_Energy`を使用する。
- 初期設定はlayer 0–29、chip 0–5（chip 5はchannel 0–29）。
- CoG・energy・channel位置のカット、bad-channel maskは適用しない。
  問題channelも見られるようにするためで、既存のshower選択とは異なる。
- `--channels 420008,920012`、`--layers 4,9-10`で対象を指定できる。
- `--max-events`は**sample全体**の上限。各ファイルの上限ではない。
  例えば1000なら最初のrunだけで上限に達する場合がある。run別の比較にはファイルを明示する。

## 図とROOTの内容

各図は左列data、中央列simulation、右列はADC binごとの平均と標準誤差の重ね描き。
data / MCで軸・binを共通にする。2次元分布はhit数のlog色表示で、
同じ行では色の上限も共通。分布の面積やイベント数では規格化していない。

### `hg_lg`

1. pedestal差し引き後のHG ADC vs LG ADC、全範囲
2. 低ADC側の拡大：HG -50–3000 ADC、LG -5–130 ADC

### `adc_energy`

1. HG ADC vs `Hit_HG_Energy`
2. LG ADC vs `Hit_LG_Energy`
3. HG ADC vs 最終採用`Hit_Energy`
4. LG ADC vs 最終採用`Hit_Energy`

3・4段目は**HG/LG選択にかかわらず全対応hit**を含む。
HGが飽和して同じADCに異なるLG由来energyが並ぶ場合も、そのまま表示する。
`GainTag`から採用枝を推定したり、raw HG 2600などを仮定して枝を選び直したりしない。

x軸はpedestalのみを引いたADCで、温度補正前。y軸は既存の全校正を適用済みのenergy。
温度による線の幅も残る。負のenergyも除去しない。

表示初期範囲はHG -100–4200、LG -50–3200 ADC、HG energy -5–40 MeV、
LG / 最終energy -5–400 MeV。`--hg-max`、`--lg-max`、`--hg-energy-max`、
`--energy-max`で上限を変更できる。bin数は`--bins`（初期値128）。
flow binもROOTに保存し、図に`outside`として範囲外hit数を表示する。
Profileはbin中心による近似ではなく元の値から計算し、y方向の表示範囲外も含む。
x方向の範囲外はProfileのflow binに保存され、図では描かれない。
空のProfile binは0の測定として描かない。

histogramの蓄積は疎なbinカウントで行い、各channelの出力時だけTH2Dに展開する。
イベントや点の間引きは行わない。全6300 channelに統計があればPNGは2枚ずつとなるため、
先にROOTだけを作り、図は対象channelに絞って作ることもできる。

```bash
bash run.sh --energy 100 --no-plots
bash run.sh --energy 100 --plots-only --channels 420008,920012
```

`--plots-only`ではbin・表示範囲は保存ROOTの設定を使用する。
異なる範囲のヒストグラムを作るにはfillから再実行する。

## 出力先

既存解析と同じ`analysis/result/`配下に保存する。
`tag=<energy>GeV`、label指定時は`tag=<energy>GeV_<label>`。

```text
analysis/result/channel_response/adc_energy/
  data/<tag>.root
  simulation/<mc-tag>/<tag>.root
  comparison_<mc-tag>_<tag>.root
  comparisons/<mc-tag>/<tag>/layer<L>.root
  figures/comparison/<mc-tag>/<tag>/layer<L>/chip<C>/
    channel<ch>_cell<CellID>_hg_lg.png
    channel<ch>_cell<CellID>_adc_energy.png
  inputs/<tag>_<mc-tag>_data.tsv
  inputs/<tag>_<mc-tag>_mc.tsv
  inputs/<tag>_<mc-tag>_configuration.txt
```

図はPNGのみ保存する。sample ROOTは`channel_<CellID>/`に
上記6種類のTH2D・対応する`*_profile`・`*_outside`、pedestal値を保存する。
comparison ROOTは`--save-canvases`を指定した場合だけ保存する。
複数layerを一つのプロセスで描く直接実行は従来の`comparison_<mc-tag>_<tag>.root`、
layer分割ジョブは`comparisons/<mc-tag>/<tag>/layer<L>.root`となる。
過去に作成したcomparison ROOTは自動削除しないため、PNGのみで再描画しても更新されない。
入力一覧・選択・照合方法・処理イベント数をsample ROOTにも記録する。
設定と実行コマンドは`configuration.txt`にshellの引用形式で保存する。
dataのPS / SPSの選択も`beam`として記録する。
標準のPS / SPSはenergyが重複しないため、出力名は従来の`<energy>GeV`のまま。

## コードの構成

- `submit.sh`: 集計ジョブと依存関係付き描画配列の投入・並列数制限。
- `plot_array.sh`: taskIDをenergy/layerに変換するwrapper。
- `plot.sh`: 入力ROOTのノード共有キャッシュ、ローカル描画、PNG転送。
- `FastPng.hh`: 同じROOTのmarker描画を小領域で実行する高速PNG出力。
- `run.sh`: 引数・入力ファイルの確認、TSV入力一覧の作成、buildと実行。
- `ChannelResponse.cc`: decode / calibの対応付け、pedestal差引き、TH2D / TProfileの保存。
- `PlotChannelResponse.cc`: ROOTからdata / MCのcanvasとPNGを作成。
- `../execute.sl`: 両解析共通のSlurm実行用wrapper。
- `../root.mk`: 両解析共通のC++・ROOTビルド設定。

## 既存結果の再描画

すでに実行中の旧ジョブは新コードを配置しても旧実行ファイルで動き続ける。
新旧のジョブが同じPNGへ書き込まないよう、新しい実行は別`--label`で分離する。
既存ROOTから`--plots-only`で再描画する場合は、その出力への旧描画ジョブが終了してから行う。

## HG–LG tailとイベント内同時発生の解析

[`../hg_lg_tail/`](../hg_lg_tail/README.md)で、同じ外れ値定義を使ってchannelごとの低LG tail、event内の同時発生、負のLG energy、総energyへの影響をまとめて比較できる。
C++とshellで実行し、結果は`analysis/result/channel_response/hg_lg_tail/`以下に保存する。既定は100 GeVのみで、全energyの投入にも対応する。
