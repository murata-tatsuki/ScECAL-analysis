# 保存済みADC分布とfitのchip別描画

C++/ROOTの描画専用コード。元の`analysis/sipm_gain`の収集・fit出力を読み、chipごとに`TCanvas::Divide(6,6)`でchannel 0–35を左上から横方向に並べます。chip 5の未接続channel 30–35は`Not connected`、その条件にADCがないchannelは`No ADC in this condition`と表示します。

各padは、黒のraw HG ADCヒストグラムに青のmulti-Gaussian・緑のenvelope-combを`Draw("L SAME")`で重ねます。candidateは実線、reviewは破線です。未採用のfitも保存曲線があれば描き、C=candidate、R=review、I=insufficient_dataを併記します。最終採用gainやその帯の図ではありません。

FFTは、保存済みFFT gainの基本周波数 `f0 = 1/gain` 付近を通す逆FFTの再構成曲線を赤で重ねます。candidateは実線、reviewは破線です。これは表示用の再構成であり、確率分布のfitやgainの再推定ではありません。自己相関はgainの数値のみを紫で記載します。`--fft-reconstruction 0`で赤の曲線を非表示、`--spectral-labels 0`でFFT/自己相関の数値だけを非表示にできます。

### FFT再構成の定義

- FFT手法自身の保存済みROIヒストグラムと`env`を読み、元と同じ `(ADC counts - env) / sqrt(max(env,1))`、Hann窓、8倍以上のゼロ埋めで複素FFTを計算します。保存済みpowerとの一致を検証します。powerのみから位相を推測せず、元ADCから位相・振幅を取り戻します。
- 通す基本周波数帯は `f0 ± min(2/(N*ΔADC), 0.45*f0)`。中心から帯域半幅の半分までは重み1、その外側はcosineで0まで落とします。高調波や別のgain候補は加えません。選択されたgain付近の周期を見やすくする設定で、鋭いピークの全形状を再現するものではありません。
- 元の複素係数にこの重みを掛け、実数への逆FFTを行いFFT点数で割ります。Hann窓と白色化を戻し、保存済み`env`を加えてADC countsに戻します。
- Hann窓が0.25以上の中央部分だけを描画し、両端での誤差増幅を避けます。ピーク高さのfit・積分値の再規格化・負値の切り捨てはしません。負値があれば縦軸にも表示します。
- gainが得られていないchannelには再構成を作らず、理由をmetadataに残します。元ADC・既存fit・FFT gain・採用判定は変更しません。

4手法の入力は同じADC分布なので、rawヒストグラムは1つを表示します。各手法のdiagnostics中のADC binがそのraw分布と一致することを確認してから重ねます。再fit、bootstrap、再binning、pedestal subtraction、校正定数の更新は行いません。既存fitのソース照合を保つため、独立した隣のディレクトリに置いています。

## 入力

- `analysis/result/sipm_gain/<label>/data/<kind>/<sample>/collection/spectra.root`の`Spectra`と`s<sid>`：4096 binのraw HG ADC。
- 同じsampleの各手法の`part_*/diagnostics.root`：`s<sid>/multi_gaussian`、`envelope_comb`、`fft`、`autocorrelation`。
- 各partの`gain_estimates.tsv`と`FIT_COMPLETE`：gain・status・元のcollectionとの対応。

PS/SPS・particle・LED日時/DACまたはbeam run・温度区間・time bin・HitTag/GainTag・memory cellを混ぜません。既定はHitTag=1、GainTag=1、memory cell=0。`--gain-tag 0`も選べますが、収集済みの枝は`HG_Charge`のままでLG分布への切り替えにはなりません。

既定`--verify-inputs metadata`は完了記録・file size/mtime・ADC binの一致を確認し、元のSHA256と検証範囲を出力に記録します。`--verify-inputs sha256`は消費するraw/diagnostics/数値TSVの全ファイルhashも照合するため、共有ディスクの読込量が増えます。入力には共有lock、出力にはchip単位の排他lockを使います。

## ビルドと代表channelの確認

ROOTとC++17の環境で実行します。

```bash
cd /home/murata_t/ScECAL_BeamTest
make -C analysis/sipm_gain_adc_plots -j2

# channel 330029を含むlayer 3 / chip 3。25.75℃が属する25.5–26.0℃区間。
./analysis/sipm_gain_adc_plots/Plot --label full_grouped \
  --sample ps_led_20230516_1400_Calib_dac2450 --layer 3 --chip 3 --temperature 25.75

# DAC 2500も別canvas。
./analysis/sipm_gain_adc_plots/Plot --label full_grouped \
  --sample ps_led_20230516_1400_Calib_dac2500 --layer 3 --chip 3 --temperature 25.75
```

`--chip`を省くとそのlayerの全chip、`--temperature`を省くと全温度区間をそれぞれ別canvasに出します。`--layer`は必須です。

表示範囲は既定`--adc-range data`（全ての非ゼロADC binとfit曲線の範囲）。`--adc-range fit`は保存fit範囲、`--adc-range full`は0–4095の全域です。どの指定でもROOTには4096 binの内容・誤差を保存し、変えるのは表示範囲だけです。

## 出力

```text
analysis/result/sipm_gain_adc_plots/<output-label>/
  inputs/source/, rendering_sources.json
  data/<kind>/<sample>/temperature_<lower>-<upper>C/
    hit_1_gain_1_memory_0/layer_03/
      chip_03.png
      chip_03.root
      chip_03.json
```

`--output-label`の既定は入力と同じlabel。温度不明は`temperature_unavailable`、time binが0以外ならlayer配下に`time_bin_<N>/`を追加します。

ROOTの`adc_fits`は36分割canvas。ROOTで開く際に`.rootlogon`等が色・文字サイズを強制変更する環境では、`root -n -l chip_03.root`で起動すると保存時の表示を確認できます。`channel_0`～`channel_35`には`raw_adc`、利用可能な2つのfit曲線とFFT/自己相関のグラフ、`metadata_json`を保存します。さらに`fft_reconstruction`（赤線のADC再構成）、`fft_envelope`（加え戻した背景）、`fft_filter`（周波数ごとの重み）、`fft_filtered_windowed_residual`（逆FFT直後、Hann/白色化補正前）を保存します。metadataにはFFT gain、周波数帯、窓、表示範囲、元powerとの一致を記録します。`raw_adc`はbin全域を保持しますが、初期表示範囲はPNGと同じです。必要なら`GetXaxis()->UnZoom()`で解除できます。ROOT全体と隣のJSONにも入力パス・条件・ソースhash・実行コマンド・各channelの状態を記録します。

同じ入力・コード・表示設定で完了済みなら出力hashを確認してスキップします。表示設定や描画コードを変えて再実行すると、そのchipのPNG/ROOTを更新します。元ADC、元fit、採用gainの出力は更新しません。

## Slurmで全channelを描画

ビルドと代表例の確認を済ませ、ROOT環境が有効な端末から都度投入します。これは描画だけを行う独立工程です。

```bash
# 1 sample × 30 layer、各layerの全chip・全温度区間。まず準備のみ。
bash analysis/sipm_gain_adc_plots/submit.sh \
  --label full_grouped --sample ps_led_20230516_1400_Calib_dac2450 --max-concurrent 8

# 上記と同じ範囲を実際に投入。
bash analysis/sipm_gain_adc_plots/submit.sh \
  --label full_grouped --sample ps_led_20230516_1400_Calib_dac2450 --max-concurrent 8 --execute

# 全LED sample。必要なら--max-concurrent 32に変更。
bash analysis/sipm_gain_adc_plots/submit.sh \
  --label full_grouped --kind led --sample all --max-concurrent 8 --execute
```

`--kind electron`、`--kind muon`、`--kind all`も選べます。全体の同時実行数が指定を超えないよう、1000task以下の配列に分け、次の配列は前の正常終了後に開始します。ログ・task一覧・投入コマンドは`analysis/log/sipm_gain_adc_plots/<output-label>/submission.*/`。全規模のジョブは実装検証では投入しません。

## 実データ確認

2023 PS LEDの20230516_1400、DAC 2450・2500、layer 3/chip 3とchip 5の一部条件で9 canvasを作成しました。36分割・全ADC bin内容/誤差・全曲線点・channel配置・未接続表示・完了済み出力の再利用を確認済みです。結果は`analysis/result/sipm_gain_adc_plots/full_grouped/validation.json`、ログは対応する`analysis/log/sipm_gain_adc_plots/full_grouped/validation/`。全LEDの3840task一覧も確認しましたが、本番ジョブは投入していません。


FFT再構成追加後は、LED・muon・electron各1 canvasで再検証しました。9 channelのFFT再構成、303252 binのraw内容/誤差、9783点の既存グラフの不変を確認済みです。人工分布の直接DFT/逆DFTとも比較し、位相・振幅・帯域外抑制、gain不明、power不一致、平坦分布を確認しています。検証用出力は `analysis/result/sipm_gain_adc_plots/full_grouped_fft_preview/`、結果は `fft_validation.json`、ログは対応する `analysis/log/sipm_gain_adc_plots/full_grouped_fft_preview/validation/` にあります。全規模の描画は同じコマンドを再実行してください。再fitは不要です。


Slurm workerは`sbatch --chdir`で指定されたproject作業ディレクトリを使います。Slurmがworkerを`/var/spool/slurmd`へ複製するため、worker自身の`BASH_SOURCE`からproject位置を求めてはいけません。`DependencyNeverSatisfied`は先行配列の正常終了条件が満たせない状態で、待機しても開始しません。失敗原因を修正した後、旧待機ジョブを取り消して再投入します。投入済みworkerはSlurm側に複製済みなので、ファイルを修正するだけでは旧待機ジョブは直りません。
